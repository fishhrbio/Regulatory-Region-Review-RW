"""Exact base-pair arithmetic; all intervals are zero-based, half-open.

Published states are inputs, not predictions trained by this project.
Copyright (c) 2026 Robert Waite. See LICENSE and data/SOURCES.md.
"""
from dataclasses import dataclass
import gzip
import re


@dataclass(frozen=True, order=True)
class Interval:
    start: int
    end: int
    state: int


@dataclass(frozen=True)
class Atom:
    start: int
    end: int
    gm: int
    k: int

    @property
    def bp(self):
        return self.end - self.start


def validate_intervals(rows):
    """Reject invalid, unsorted or overlapping intervals, including duplicates."""
    for i, row in enumerate(rows):
        if not (0 <= row.start < row.end and 1 <= row.state <= 15):
            raise ValueError(f"Invalid interval: {row}")
        if i and rows[i - 1].end > row.start:
            raise ValueError(f"Unsorted/overlapping intervals: {rows[i-1]}, {row}")


def read_track(path, chromosome):
    """Stream whole BED for structural QC; retain only requested chromosome.

    No missing state is imputed. Malformed records fail the run. Selection is
    sorted explicitly; source order and selected-coordinate gaps are reported.
    """
    rows = []
    total = 0
    with gzip.open(path, "rt", encoding="utf-8") as source:
        for number, line in enumerate(source, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or not fields[0] or not re.fullmatch(r"(?:[1-9]|1[0-5])_.+", fields[3]):
                raise ValueError(f"Malformed BED record {path}:{number}")
            start, end = int(fields[1]), int(fields[2])
            state = int(fields[3].split("_", 1)[0])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid coordinates {path}:{number}")
            total += 1
            if fields[0] == chromosome:
                rows.append(Interval(start, end, state))
    if not rows:
        raise ValueError(f"No records for {chromosome}")
    was_sorted = rows == sorted(rows)
    duplicates = len(rows) - len(set(rows))
    rows.sort()
    validate_intervals(rows)
    qc = {
        "whole_file_rows": total, "selected_rows": len(rows),
        "other_chromosome_rows": total - len(rows),
        "selected_bp": sum(x.end-x.start for x in rows),
        "selected_span_start": rows[0].start, "selected_span_end": rows[-1].end,
        "internal_gap_bp": sum(b.start-a.end for a,b in zip(rows,rows[1:])),
        "internal_gap_count": sum(b.start>a.end for a,b in zip(rows,rows[1:])),
        "input_selected_sorted": was_sorted, "duplicate_selected_rows": duplicates,
        "overlapping_selected_rows": 0, "malformed_or_missing_records": 0,
        "observed_states": sorted({x.state for x in rows})
    }
    return rows, qc


def intersect_tracks(gm, k):
    """Linear two-pointer sweep. Return atoms only where BOTH tracks cover.

    Every overlap is split at a boundary from either track. O(n+m) time,
    O(n+m) output; no per-base genome table or rounding to 200-bp bins.
    """
    validate_intervals(gm)
    validate_intervals(k)
    result = []
    i = j = 0
    while i < len(gm) and j < len(k):
        a, b = gm[i], k[j]
        start, end = max(a.start,b.start), min(a.end,b.end)
        if start < end:
            result.append(Atom(start, end, a.state, b.state))
        if a.end <= b.end:
            i += 1
        if b.end <= a.end:
            j += 1
    return result


def eligible_atoms(atoms, excluded_states):
    excluded = set(excluded_states)
    return [a for a in atoms if a.gm not in excluded and a.k not in excluded]


def class_overlap(atoms, states):
    """Within common eligible territory, union partitions into 3 exclusive sets.

    Jaccard and directional retention are undefined (None) for a zero
    denominator. Absence of any annotated bases is never reported as overlap 0.
    """
    selected = set(states)
    gm_only = shared = k_only = 0
    for a in atoms:
        g, k = a.gm in selected, a.k in selected
        if g and k:
            shared += a.bp
        elif g:
            gm_only += a.bp
        elif k:
            k_only += a.bp
    union = gm_only + shared + k_only
    g_total, k_total = gm_only+shared, k_only+shared
    return {"gm_only_bp": gm_only, "shared_bp": shared, "k_only_bp": k_only,
            "union_bp": union, "gm_bp": g_total, "k_bp": k_total,
            "jaccard": shared/union if union else None,
            "gm_retained": shared/g_total if g_total else None,
            "k_retained": shared/k_total if k_total else None}


def state_matrix(atoms):
    matrix = [[0]*15 for _ in range(15)]
    for a in atoms:
        matrix[a.gm-1][a.k-1] += a.bp
    return matrix


def summarize(atoms, classes, excluded_states):
    retained = eligible_atoms(atoms, excluded_states)
    territory_bp = sum(a.bp for a in retained)
    if not territory_bp:
        raise ValueError("Eligible territory is empty")
    matrix = state_matrix(retained)
    assert sum(map(sum,matrix)) == territory_bp
    coverage = []
    for state in range(1,16):
        gm_bp = sum(matrix[state-1])
        k_bp = sum(row[state-1] for row in matrix)
        coverage.append({"state": state, "gm_bp": gm_bp, "k_bp": k_bp,
                         "gm_fraction": gm_bp/territory_bp, "k_fraction": k_bp/territory_bp})
    return {"eligible_bp": territory_bp, "atom_count": len(retained),
            "excluded_bp": sum(a.bp for a in atoms)-territory_bp,
            "excluded_states": list(excluded_states), "coverage": coverage,
            "overlap": {name:class_overlap(retained,states) for name,states in classes.items()},
            "state_matrix_bp": matrix,
            "same_state_fraction": sum(matrix[i][i] for i in range(15))/territory_bp}
