#!/usr/bin/env python3

"""
Disk-backed read-name registry mapping long read names <-> compact integer IDs.

Primary backend: LMDB (if available). Fallback: SQLite3 (bundled with Python).

API:
  - ReadNameStore(db_path): open or create store at path (without extension).
  - get_or_add(name: str) -> int: returns stable integer ID for given name, creating if new.
  - get_id(name: str) -> Optional[int]: returns ID if exists else None.
  - get_name(id: int) -> Optional[str]: reverse lookup.
  - close(): close any handles.

Thread/process safety: Backends use their own locking. This class is intended for
single-process write (append) with read access from same process.
"""

import os
import sqlite3
from typing import Optional

try:
    import lmdb  # type: ignore
    _HAS_LMDB = True
except Exception:
    lmdb = None  # type: ignore
    _HAS_LMDB = False


class ReadNameStore:
    def __init__(self, db_path_base: str):
        """
        db_path_base: base path for the store; for LMDB a directory, for SQLite a .sqlite file is created.
        """
        # Allow forcing backend via env var for portability/debugging
        force_backend = os.environ.get("LRAA_READSTORE_BACKEND", "").lower()
        self._use_lmdb = _HAS_LMDB and (force_backend not in ("sqlite", "sql", "s"))
        self._closed = False
        if self._use_lmdb:
            # LMDB expects a directory; create if not exists
            self._path = db_path_base if os.path.isdir(db_path_base) else db_path_base + ".lmdb"
            os.makedirs(self._path, exist_ok=True)
            # map_size: 1GB default, grows by re-open if needed; big enough for many names
            self._env = lmdb.open(self._path, map_size=1 << 30, max_dbs=4, subdir=True, lock=True)
            self._db_id2name = self._env.open_db(b"id2name")
            self._db_name2id = self._env.open_db(b"name2id")
            # init counter
            with self._env.begin(write=True) as txn:
                if txn.get(b"__next_id__") is None:
                    txn.put(b"__next_id__", b"1")
        else:
            # SQLite fallback
            self._path = db_path_base if db_path_base.endswith(".sqlite") else db_path_base + ".sqlite"
            os.makedirs(os.path.dirname(self._path) or ".", exist_ok=True)
            self._conn = sqlite3.connect(self._path, timeout=60.0)
            self._conn.execute(
                "CREATE TABLE IF NOT EXISTS read_names (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT UNIQUE)"
            )
            self._conn.execute("PRAGMA journal_mode=WAL;")
            self._conn.execute("PRAGMA synchronous=NORMAL;")
            self._conn.commit()

    def get_or_add(self, name: str) -> int:
        if self._use_lmdb:
            key_name = name.encode("utf-8", errors="ignore")
            # Use a single write transaction to avoid nested writer deadlocks
            with self._env.begin(write=True) as txn:
                val = txn.get(key_name, db=self._db_name2id)
                if val is not None:
                    return int(val.decode("ascii"))
                next_id_bytes = txn.get(b"__next_id__")
                next_id = int(next_id_bytes.decode("ascii")) if next_id_bytes else 1
                new_id = next_id
                txn.put(b"__next_id__", str(next_id + 1).encode("ascii"))
                # write both directions under same txn
                txn.put(key_name, str(new_id).encode("ascii"), db=self._db_name2id)
                txn.put(str(new_id).encode("ascii"), key_name, db=self._db_id2name)
                return new_id
        else:
            cur = self._conn.cursor()
            try:
                cur.execute("INSERT OR IGNORE INTO read_names(name) VALUES (?)", (name,))
                self._conn.commit()
            except sqlite3.Error:
                pass
            # fetch id
            cur.execute("SELECT id FROM read_names WHERE name=?", (name,))
            row = cur.fetchone()
            if row is None:
                raise RuntimeError("Failed to insert or fetch read name")
            return int(row[0])

    def get_id(self, name: str) -> Optional[int]:
        if self._use_lmdb:
            key_name = name.encode("utf-8", errors="ignore")
            with self._env.begin(db=self._db_name2id) as txn:
                val = txn.get(key_name)
                return int(val.decode("ascii")) if val is not None else None
        else:
            cur = self._conn.cursor()
            cur.execute("SELECT id FROM read_names WHERE name=?", (name,))
            row = cur.fetchone()
            return int(row[0]) if row else None

    def get_name(self, rid: int) -> Optional[str]:
        if self._use_lmdb:
            key_id = str(rid).encode("ascii")
            with self._env.begin(db=self._db_id2name) as txn:
                val = txn.get(key_id)
                return val.decode("utf-8") if val is not None else None
        else:
            cur = self._conn.cursor()
            cur.execute("SELECT name FROM read_names WHERE id=?", (int(rid),))
            row = cur.fetchone()
            return str(row[0]) if row else None

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            if self._use_lmdb:
                self._env.close()
            else:
                self._conn.close()
        except Exception:
            pass
