from __future__ import annotations

from itertools import pairwise
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import networkx as nx

from tapis.topology.enums import GraphEdgeKind, GraphNodeKind

Exon = tuple[int, int]
TranscriptMap = dict[str, list[Exon]]


def build_annotation_graph(
    transcripts: Mapping[str, Sequence[Exon]],
    chromosome: str,
    strand: str,
) -> nx.DiGraph:
    graph = nx.DiGraph(chromosome=chromosome, strand=strand)
    for transcript_id, exon_chain in sorted(transcripts.items()):
        exons = normalize_exons(exon_chain)
        for exon in exons:
            graph.add_node(
                exon,
                kind=GraphNodeKind.EXON.value,
                start=exon[0],
                end=exon[1],
            )
        for left, right in pairwise(exons):
            if graph.has_edge(left, right):
                graph[left][right]["transcripts"].add(transcript_id)
                continue
            graph.add_edge(
                left,
                right,
                kind=GraphEdgeKind.SPLICE.value,
                donor=left[1],
                acceptor=right[0],
                transcripts={transcript_id},
            )
    return graph


def build_annotation_graph_from_gtf(
    gtf_path: str | Path,
    chromosome: str = "unknown",
    strand: str = "+",
) -> nx.DiGraph:
    transcripts = parse_gtf_exons(gtf_path)
    return build_annotation_graph(transcripts, chromosome=chromosome, strand=strand)


def parse_gtf_exons(gtf_path: str | Path) -> TranscriptMap:
    transcript_exons: TranscriptMap = {}
    for raw_line in Path(gtf_path).read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9 or fields[2].lower() != "exon":
            continue
        attributes = parse_gtf_attributes(fields[8])
        transcript_id = attributes.get("transcript_id")
        if transcript_id is None:
            continue
        start = int(fields[3])
        end = int(fields[4])
        transcript_exons.setdefault(transcript_id, []).append((start, end))

    for transcript_id, exons in transcript_exons.items():
        transcript_exons[transcript_id] = normalize_exons(exons)
    return transcript_exons


def extract_read_exons_from_bam(
    bam_path: str | Path,
    min_mapq: int = 0,
) -> TranscriptMap:
    import pysam

    transcripts: TranscriptMap = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam_stream:
        for read in bam_stream.fetch(until_eof=True):
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue
            blocks = read.get_blocks()
            if not blocks:
                continue
            if read.query_name is None:
                continue
            transcripts[read.query_name] = normalize_exons(blocks)
    return transcripts


def normalize_exons(exons: Iterable[Exon]) -> list[Exon]:
    normalized = sorted((int(start), int(end)) for start, end in exons)
    return normalized


def parse_gtf_attributes(attributes: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for chunk in attributes.split(";"):
        item = chunk.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        parsed[key] = value.strip().strip('"')
    return parsed
