import logging
import re
from typing import Dict, Tuple, Optional

logger = logging.getLogger(__name__)

# Priority ranking for gffcompare class codes (lower value means higher specificity)
_DEFAULT_CLASS_CODE_PRIORITY: Dict[str, int] = {
    "=": 0,
    "c": 1,
    "k": 2,
    "m": 3,
    "n": 4,
    "j": 5,
    "e": 6,
    "o": 7,
    "s": 8,
    "x": 9,
    "i": 10,
    "y": 11,
    "p": 12,
    "r": 13,
    "u": 14,
}


def get_default_class_code_priority() -> Dict[str, int]:
    """Return a copy of the default class code priorities."""

    return dict(_DEFAULT_CLASS_CODE_PRIORITY)


def parse_gffcompare_mappings(
    gffcompare_tracking_filename: str,
    class_code_priority: Optional[Dict[str, int]] = None,
) -> Dict[str, Tuple[str, str]]:
    """Parse gffcompare tracking output into target-id → (ref_gene_id, ref_trans_id).

    When multiple assignments exist for a target id, favor the one with the most
    specific (lowest-ranked) class code according to the provided priority map.
    """

    if class_code_priority is None:
        class_code_priority = _DEFAULT_CLASS_CODE_PRIORITY

    logger.info("-parsing gffcompare output: %s", gffcompare_tracking_filename)

    mappings: Dict[str, Tuple[str, str, int]] = {}

    def get_priority(code: str) -> int:
        return class_code_priority.get(code, 100)

    def consider_assignment(target_id: Optional[str], ref_gene_id: str, ref_trans_id: str, priority: int) -> None:
        if not target_id or target_id == "-":
            return

        existing = mappings.get(target_id)
        if existing is None or priority < existing[2]:
            mappings[target_id] = (ref_gene_id, ref_trans_id, priority)

    with open(gffcompare_tracking_filename, "rt") as fh:
        for raw_line in fh:
            line = raw_line.rstrip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                continue

            _, _, ref_info, compare_code, target_info = parts[:5]

            if ref_info == "-":
                continue

            try:
                ref_gene_id, ref_trans_id = ref_info.split("|", 1)
            except ValueError:
                logger.warning("Skipping malformed ref info: %s", ref_info)
                continue

            priority = get_priority(compare_code)

            # Trackings can include multiple q#: entries separated by ',' or ';'.
            for entry in re.split(r"[,;]", target_info):
                candidate = entry.strip()
                if not candidate:
                    continue

                # Remove leading queue token (e.g. q1:...)
                token_split = candidate.split(":", 1)
                if len(token_split) == 2:
                    candidate = token_split[1]

                fields = candidate.split("|")
                if len(fields) < 2:
                    continue

                target_gene_id = fields[0]
                target_trans_id = fields[1]

                consider_assignment(target_gene_id, ref_gene_id, ref_trans_id, priority)
                consider_assignment(target_trans_id, ref_gene_id, ref_trans_id, priority)

    # Strip priority metadata before returning
    return {target_id: (ref_gene_id, ref_trans_id) for target_id, (ref_gene_id, ref_trans_id, _) in mappings.items()}


def get_ref_gene_names(ref_gtf: str) -> Dict[str, str]:
    """Extract dictionaries mapping reference ids to gene symbols."""

    logger.info("-extracting gene_names and identifiers from reference gtf: %s", ref_gtf)

    ref_id_to_gene_name: Dict[str, str] = {}

    with open(ref_gtf, "rt") as fh:
        for line in fh:
            vals = line.split("\t")
            if len(vals) < 9:
                continue

            if vals[2] != "transcript":
                continue

            info = vals[8]

            m = re.search('gene_id "([^\"]+)";.*transcript_id "([^\"]+)";', info)
            if not m:
                continue

            gene_id = m.group(1)
            transcript_id = m.group(2)
            gene_name = gene_id

            m2 = re.search(' gene_name "([^\"]+)";', info)
            if m2:
                gene_name = m2.group(1)

            ref_id_to_gene_name[transcript_id] = gene_name
            ref_id_to_gene_name[gene_id] = gene_name

    return ref_id_to_gene_name


def resolve_gene_symbol(
    target_id: Optional[str],
    ref_id_to_gene_name: Dict[str, str],
    gffcompare_mapping: Optional[Dict[str, Tuple[str, str]]] = None,
    prefer_transcript_first: bool = True,
) -> Optional[str]:
    """Resolve the gene symbol for a target id using direct or mapped lookups."""

    if not target_id:
        return None

    if target_id in ref_id_to_gene_name:
        return ref_id_to_gene_name[target_id]

    if not gffcompare_mapping or target_id not in gffcompare_mapping:
        return None

    ref_gene_id, ref_trans_id = gffcompare_mapping[target_id]

    lookup_order = []
    if prefer_transcript_first:
        lookup_order = [ref_trans_id, ref_gene_id]
    else:
        lookup_order = [ref_gene_id, ref_trans_id]

    for ref_id in lookup_order:
        if ref_id and ref_id in ref_id_to_gene_name:
            return ref_id_to_gene_name[ref_id]

    return None
