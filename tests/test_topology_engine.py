from __future__ import annotations

from tapis.topology.events import detect_splice_events
from tapis.topology.graph import build_annotation_graph


def test_detects_exon_skipping_and_intron_retention() -> None:
    annotation = {
        "tx_ref": [(100, 150), (200, 250), (300, 350)],
    }
    observed = {
        "read_skip": [(100, 150), (300, 350)],
        "read_ir": [(100, 250), (300, 350)],
    }
    graph = build_annotation_graph(annotation, chromosome="chr1", strand="+")
    events = detect_splice_events(graph, observed)

    kinds = {(event.kind, event.detail) for event in events}
    assert ("exon_skipping", "200-250") in kinds
    assert ("intron_retention", "151-199") in kinds


def test_detects_alternative_splice_site() -> None:
    annotation = {
        "tx_ref": [(100, 150), (200, 250), (300, 350)],
    }
    observed = {
        "read_alt": [(100, 150), (210, 250), (300, 350)],
    }
    graph = build_annotation_graph(annotation, chromosome="chr1", strand="+")
    events = detect_splice_events(graph, observed)

    kinds = {(event.kind, event.detail) for event in events}
    assert ("alternative_splice_site", "151-210") in kinds
