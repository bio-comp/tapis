from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from typing import Mapping, Sequence

import networkx as nx

from tapis.topology.enums import SpliceEventKind
from tapis.topology.graph import Exon, normalize_exons


@dataclass(frozen=True)
class SpliceEvent:
    kind: SpliceEventKind
    read_id: str
    chromosome: str
    strand: str
    detail: str
    anchor_left: Exon | None = None
    anchor_right: Exon | None = None


def detect_splice_events(
    graph: nx.DiGraph,
    observed_transcripts: Mapping[str, Sequence[Exon]],
) -> list[SpliceEvent]:
    chromosome = str(graph.graph.get("chromosome", "unknown"))
    strand = str(graph.graph.get("strand", "+"))
    event_index: dict[
        tuple[SpliceEventKind, str, str, Exon | None, Exon | None], SpliceEvent
    ] = {}
    annotation_junctions = [(left, right) for left, right in graph.edges()]

    for read_id, observed in sorted(observed_transcripts.items()):
        exon_chain = normalize_exons(observed)
        _collect_exon_skipping_events(
            graph,
            read_id,
            exon_chain,
            chromosome,
            strand,
            event_index,
        )
        _collect_alternative_splice_site_events(
            read_id,
            exon_chain,
            annotation_junctions,
            chromosome,
            strand,
            event_index,
        )
        _collect_intron_retention_events(
            graph,
            read_id,
            exon_chain,
            chromosome,
            strand,
            event_index,
        )

    return sorted(
        event_index.values(),
        key=lambda event: (event.kind.value, event.read_id, event.detail),
    )


def _collect_exon_skipping_events(
    graph: nx.DiGraph,
    read_id: str,
    exons: list[Exon],
    chromosome: str,
    strand: str,
    event_index: dict[
        tuple[SpliceEventKind, str, str, Exon | None, Exon | None], SpliceEvent
    ],
) -> None:
    for left, right in pairwise(exons):
        if graph.has_edge(left, right):
            continue
        if not graph.has_node(left) or not graph.has_node(right):
            continue
        try:
            paths = list(nx.all_shortest_paths(graph, left, right))
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            continue
        if not paths:
            continue
        path = sorted(paths)[0]
        if len(path) <= 2:
            continue
        skipped = path[1:-1]
        detail = ",".join(f"{start}-{end}" for start, end in skipped)
        _add_event(
            event_index=event_index,
            event=SpliceEvent(
                kind=SpliceEventKind.EXON_SKIPPING,
                read_id=read_id,
                chromosome=chromosome,
                strand=strand,
                detail=detail,
                anchor_left=left,
                anchor_right=right,
            ),
        )


def _collect_alternative_splice_site_events(
    read_id: str,
    exons: list[Exon],
    annotation_junctions: list[tuple[Exon, Exon]],
    chromosome: str,
    strand: str,
    event_index: dict[
        tuple[SpliceEventKind, str, str, Exon | None, Exon | None], SpliceEvent
    ],
) -> None:
    for left, right in pairwise(exons):
        donor = left[1]
        acceptor = right[0]
        has_same_donor = any(
            a_left[1] == donor and a_right[0] != acceptor
            for a_left, a_right in annotation_junctions
        )
        has_same_acceptor = any(
            a_right[0] == acceptor and a_left[1] != donor
            for a_left, a_right in annotation_junctions
        )
        if not has_same_donor and not has_same_acceptor:
            continue

        detail = f"{donor + 1}-{acceptor}"
        _add_event(
            event_index=event_index,
            event=SpliceEvent(
                kind=SpliceEventKind.ALTERNATIVE_SPLICE_SITE,
                read_id=read_id,
                chromosome=chromosome,
                strand=strand,
                detail=detail,
                anchor_left=left,
                anchor_right=right,
            ),
        )


def _collect_intron_retention_events(
    graph: nx.DiGraph,
    read_id: str,
    exons: list[Exon],
    chromosome: str,
    strand: str,
    event_index: dict[
        tuple[SpliceEventKind, str, str, Exon | None, Exon | None], SpliceEvent
    ],
) -> None:
    for exon in exons:
        for left, right in graph.edges():
            intron_start = left[1] + 1
            intron_end = right[0] - 1
            if intron_start > intron_end:
                continue
            if exon[0] <= intron_start and exon[1] >= intron_end:
                detail = f"{intron_start}-{intron_end}"
                _add_event(
                    event_index=event_index,
                    event=SpliceEvent(
                        kind=SpliceEventKind.INTRON_RETENTION,
                        read_id=read_id,
                        chromosome=chromosome,
                        strand=strand,
                        detail=detail,
                        anchor_left=left,
                        anchor_right=right,
                    ),
                )


def _add_event(
    event_index: dict[tuple[SpliceEventKind, str, str, Exon | None, Exon | None], SpliceEvent],
    event: SpliceEvent,
) -> None:
    key = (event.kind, event.read_id, event.detail, event.anchor_left, event.anchor_right)
    event_index[key] = event
