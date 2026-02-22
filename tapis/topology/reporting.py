from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

from tapis.topology.enums import SpliceEventKind
from tapis.topology.events import SpliceEvent, detect_splice_events
from tapis.topology.graph import (
    Exon,
    build_annotation_graph,
    extract_read_exons_from_bam,
    parse_gtf_exons,
)

EVENT_ORDER = (
    SpliceEventKind.INTRON_RETENTION,
    SpliceEventKind.EXON_SKIPPING,
    SpliceEventKind.ALTERNATIVE_SPLICE_SITE,
)


def compute_events_from_transcripts(
    annotation_transcripts: Mapping[str, Sequence[Exon]],
    observed_transcripts: Mapping[str, Sequence[Exon]],
    chromosome: str,
    strand: str,
) -> list[SpliceEvent]:
    graph = build_annotation_graph(annotation_transcripts, chromosome=chromosome, strand=strand)
    return detect_splice_events(graph, observed_transcripts)


def compute_as_statistics_from_gtf(gtf_path: str | Path) -> Counter[SpliceEventKind]:
    transcripts = parse_gtf_exons(gtf_path)
    counts: Counter[SpliceEventKind] = Counter()
    if len(transcripts) <= 1:
        return counts

    canonical_transcript_id = _select_canonical_transcript(transcripts)
    canonical = {canonical_transcript_id: transcripts[canonical_transcript_id]}
    observed = {
        transcript_id: exons
        for transcript_id, exons in transcripts.items()
        if transcript_id != canonical_transcript_id
    }
    events = compute_events_from_transcripts(
        annotation_transcripts=canonical,
        observed_transcripts=observed,
        chromosome="gtf",
        strand="+",
    )
    for event in events:
        counts[event.kind] += 1
    return counts


def write_as_statistics_from_gtf(gtf_path: str | Path, output_path: str | Path) -> None:
    counts = compute_as_statistics_from_gtf(gtf_path)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("event\tcount\n")
        for event_kind in EVENT_ORDER:
            handle.write(f"{event_kind.value}\t{counts.get(event_kind, 0)}\n")


def write_event_table(events: Sequence[SpliceEvent], output_path: str | Path) -> None:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    ordered_events = sorted(
        events,
        key=lambda event: (event.kind.value, event.read_id, event.detail),
    )
    with output.open("w", encoding="utf-8") as handle:
        handle.write("kind\tread_id\tchromosome\tstrand\tdetail\n")
        for event in ordered_events:
            handle.write(
                f"{event.kind.value}\t{event.read_id}\t"
                f"{event.chromosome}\t{event.strand}\t{event.detail}\n"
            )


def write_event_table_from_gtf_and_bam(
    gtf_path: str | Path,
    bam_path: str | Path,
    output_path: str | Path,
    min_mapq: int = 0,
) -> None:
    transcripts = parse_gtf_exons(gtf_path)
    observed = extract_read_exons_from_bam(bam_path, min_mapq=min_mapq)
    chromosome, strand = _infer_reference_context_from_gtf(gtf_path)
    events = compute_events_from_transcripts(
        annotation_transcripts=transcripts,
        observed_transcripts=observed,
        chromosome=chromosome,
        strand=strand,
    )
    write_event_table(events, output_path)


def _select_canonical_transcript(transcripts: dict[str, list[tuple[int, int]]]) -> str:
    def _key(item: tuple[str, list[tuple[int, int]]]) -> tuple[int, int, str]:
        transcript_id, exons = item
        span = sum(end - start + 1 for start, end in exons)
        return (len(exons), span, transcript_id)

    return max(transcripts.items(), key=_key)[0]


def _infer_reference_context_from_gtf(gtf_path: str | Path) -> tuple[str, str]:
    for raw_line in Path(gtf_path).read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        chromosome = fields[0]
        strand = fields[6] if fields[6] in {"+", "-"} else "+"
        return chromosome, strand
    return "unknown", "+"
