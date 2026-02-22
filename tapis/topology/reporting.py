from __future__ import annotations

from collections import Counter
from pathlib import Path

from tapis.topology.enums import SpliceEventKind
from tapis.topology.events import detect_splice_events
from tapis.topology.graph import build_annotation_graph, parse_gtf_exons

EVENT_ORDER = (
    SpliceEventKind.INTRON_RETENTION.value,
    SpliceEventKind.EXON_SKIPPING.value,
    SpliceEventKind.ALTERNATIVE_SPLICE_SITE.value,
)


def compute_as_statistics_from_gtf(gtf_path: str | Path) -> Counter[str]:
    transcripts = parse_gtf_exons(gtf_path)
    counts: Counter[str] = Counter()
    if len(transcripts) <= 1:
        return counts

    canonical_transcript_id = _select_canonical_transcript(transcripts)
    canonical = {canonical_transcript_id: transcripts[canonical_transcript_id]}
    observed = {
        transcript_id: exons
        for transcript_id, exons in transcripts.items()
        if transcript_id != canonical_transcript_id
    }
    graph = build_annotation_graph(canonical, chromosome="gtf", strand="+")
    events = detect_splice_events(graph, observed)
    for event in events:
        counts[event.kind] += 1
    return counts


def write_as_statistics_from_gtf(gtf_path: str | Path, output_path: str | Path) -> None:
    counts = compute_as_statistics_from_gtf(gtf_path)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("event\tcount\n")
        for event_name in EVENT_ORDER:
            handle.write(f"{event_name}\t{counts.get(event_name, 0)}\n")


def _select_canonical_transcript(transcripts: dict[str, list[tuple[int, int]]]) -> str:
    def _key(item: tuple[str, list[tuple[int, int]]]) -> tuple[int, int, str]:
        transcript_id, exons = item
        span = sum(end - start + 1 for start, end in exons)
        return (len(exons), span, transcript_id)

    return max(transcripts.items(), key=_key)[0]
