from __future__ import annotations

from pathlib import Path

import pysam

from tapis.topology.enums import SpliceEventKind
from tapis.topology.reporting import write_event_table_from_gtf_and_bam


def _write_read(
    out: pysam.AlignmentFile,
    name: str,
    reference_start: int,
    cigar: list[tuple[int, int]],
    query_length: int,
) -> None:
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = "A" * query_length
    read.flag = 0
    read.reference_id = 0
    read.reference_start = reference_start
    read.mapping_quality = 60
    read.cigar = cigar
    read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    out.write(read)


def test_write_event_table_from_gtf_and_bam(tmp_path: Path) -> None:
    gtf_path = tmp_path / "annotations.gtf"
    gtf_path.write_text(
        "\n".join(
            [
                (
                    "chr1\ttapis\texon\t100\t149\t.\t+\t.\t"
                    'gene_id "g1"; transcript_id "tx_ref"; exon_number "1";'
                ),
                (
                    "chr1\ttapis\texon\t200\t249\t.\t+\t.\t"
                    'gene_id "g1"; transcript_id "tx_ref"; exon_number "2";'
                ),
                (
                    "chr1\ttapis\texon\t300\t349\t.\t+\t.\t"
                    'gene_id "g1"; transcript_id "tx_ref"; exon_number "3";'
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    bam_path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # Exon-skipping read: exon1 -> exon3
        _write_read(
            out=out,
            name="read_skip",
            reference_start=99,
            cigar=[(0, 50), (3, 150), (0, 50)],
            query_length=100,
        )
        # Intron-retention read covering exon1+intron+exon2
        _write_read(
            out=out,
            name="read_ir",
            reference_start=99,
            cigar=[(0, 150)],
            query_length=150,
        )

    out_path = tmp_path / "AS_events.tsv"
    write_event_table_from_gtf_and_bam(gtf_path=gtf_path, bam_path=bam_path, output_path=out_path)

    lines = out_path.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "kind\tread_id\tchromosome\tstrand\tdetail"
    body = "\n".join(lines[1:])
    assert SpliceEventKind.EXON_SKIPPING.value in body
    assert SpliceEventKind.INTRON_RETENTION.value in body
