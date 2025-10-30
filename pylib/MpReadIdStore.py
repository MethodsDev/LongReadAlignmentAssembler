#!/usr/bin/env python3

"""
Disk-backed mapping from MultiPath IDs -> stream of read IDs (ints).

Primary backend: LMDB (if available). Fallback: SQLite3.

API:
  - MpReadIdStore(db_path_base): open/create store.
  - append(mp_id: str, read_id: int): append one association.
  - iter_read_ids(mp_id: str): iterator yielding read IDs for the multipath.
  - count(mp_id: str) -> int: number of associated read IDs.
  - close(): close handles.

Note: For LMDB, we encode entries as keys f"{mp_id}:{seq}" -> read_id to avoid dupsort complexity.
      A per-mp sequence counter is maintained under key "__seq__:{mp_id}".
"""

import os
import sqlite3
from typing import Iterator

try:
    import lmdb  # type: ignore
    _HAS_LMDB = True
except Exception:
    lmdb = None  # type: ignore
    _HAS_LMDB = False


class MpReadIdStore:
    def __init__(self, db_path_base: str):
        force_backend = os.environ.get("LRAA_READSTORE_BACKEND", "").lower()
        self._use_lmdb = _HAS_LMDB and (force_backend not in ("sqlite", "sql", "s"))
        self._closed = False
        if self._use_lmdb:
            self._path = db_path_base if os.path.isdir(db_path_base) else db_path_base + ".lmdb"
            os.makedirs(self._path, exist_ok=True)
            relaxed = os.environ.get("LRAA_LMDB_RELAXED", "1").lower() not in ("0", "false", "no")
            self._env = lmdb.open(
                self._path,
                map_size=1 << 30,
                max_dbs=4,
                subdir=True,
                lock=True,
                sync=False if relaxed else True,
                map_async=True if relaxed else False,
                metasync=False if relaxed else True,
            )
            self._db_pairs = self._env.open_db(b"mp_to_reads")
            self._db_seq = self._env.open_db(b"mp_seq")
            # batching state
            try:
                self._batch_n = int(os.environ.get("LRAA_STORE_BATCH_N", "1000"))
            except Exception:
                self._batch_n = 1000
            self._batch_n = max(1, self._batch_n)
            self._ops_since_commit = 0
            self._wtxn = None
        else:
            self._path = db_path_base if db_path_base.endswith(".sqlite") else db_path_base + ".sqlite"
            os.makedirs(os.path.dirname(self._path) or ".", exist_ok=True)
            self._conn = sqlite3.connect(self._path, timeout=60.0)
            self._conn.execute(
                "CREATE TABLE IF NOT EXISTS mp_reads (mp_id TEXT, read_id INTEGER)"
            )
            self._conn.execute("CREATE INDEX IF NOT EXISTS idx_mp_id ON mp_reads(mp_id)")
            self._conn.execute("PRAGMA journal_mode=WAL;")
            self._conn.execute("PRAGMA synchronous=NORMAL;")
            self._conn.commit()

    def _begin_wtxn(self):
        if self._wtxn is None:
            self._wtxn = self._env.begin(write=True)

    def _commit_wtxn_if_needed(self, force: bool = False):
        if self._wtxn is None:
            return
        if force or (self._ops_since_commit >= self._batch_n):
            try:
                self._wtxn.commit()
            except Exception:
                try:
                    self._wtxn.abort()
                except Exception:
                    pass
            finally:
                self._wtxn = None
                self._ops_since_commit = 0

    def append(self, mp_id: str, read_id: int):
        if self._use_lmdb:
            mp_key = mp_id.encode("utf-8")
            self._begin_wtxn()
            txn = self._wtxn
            assert txn is not None
            seq_key = b"__seq__:" + mp_key
            seq_val = txn.get(seq_key, db=self._db_seq)
            next_seq = int(seq_val.decode("ascii")) + 1 if seq_val else 1
            txn.put(seq_key, str(next_seq).encode("ascii"), db=self._db_seq)
            key = mp_key + b":" + str(next_seq).encode("ascii")
            txn.put(key, str(int(read_id)).encode("ascii"), db=self._db_pairs)
            self._ops_since_commit += 1
            self._commit_wtxn_if_needed(False)
        else:
            self._conn.execute(
                "INSERT INTO mp_reads(mp_id, read_id) VALUES (?, ?)", (mp_id, int(read_id))
            )
            self._conn.commit()

    def iter_read_ids(self, mp_id: str) -> Iterator[int]:
        if self._use_lmdb:
            mp_key_prefix = mp_id.encode("utf-8") + b":"
            with self._env.begin(db=self._db_pairs) as txn:
                with txn.cursor() as cur:
                    if cur.set_range(mp_key_prefix):
                        for k, v in cur:
                            if not k.startswith(mp_key_prefix):
                                break
                            try:
                                yield int(v.decode("ascii"))
                            except Exception:
                                continue
                    else:
                        return
        else:
            cur = self._conn.cursor()
            cur.execute("SELECT read_id FROM mp_reads WHERE mp_id=?", (mp_id,))
            for row in cur:
                yield int(row[0])

    def count(self, mp_id: str) -> int:
        if self._use_lmdb:
            # Count by iterating; sufficient for reporting
            return sum(1 for _ in self.iter_read_ids(mp_id))
        else:
            cur = self._conn.cursor()
            cur.execute("SELECT COUNT(1) FROM mp_reads WHERE mp_id=?", (mp_id,))
            row = cur.fetchone()
            return int(row[0]) if row else 0

    def close(self):
        if self._closed:
            return
        self._closed = True
        try:
            if self._use_lmdb:
                try:
                    self._commit_wtxn_if_needed(True)
                except Exception:
                    pass
                self._env.close()
            else:
                self._conn.close()
        except Exception:
            pass

    def flush(self):
        try:
            if self._use_lmdb:
                self._commit_wtxn_if_needed(True)
            else:
                # noop for sqlite
                pass
        except Exception:
            pass
