from __future__ import annotations

from pathlib import Path

from tapis.topology.enums import SpliceEventKind
from tapis.topology.reporting import write_as_statistics_from_gtf


def test_writes_as_statistics_from_gtf(tmp_path: Path) -> None:
    gtf_path = tmp_path / "assembled.gtf"
    gtf_path.write_text(
        "\n".join(
            [
                (
                    "chr1\ttapis\texon\t100\t150\t.\t+\t.\t"
                    'gene_id "g1"; transcript_id "g1_t1"; exon_number "1";'
                ),
                (
                    "chr1\ttapis\texon\t200\t250\t.\t+\t.\t"
                    'gene_id "g1"; transcript_id "g1_t1"; exon_number "2";'
                ),
                (
                    "chr1\ttapis\texon\t300\t350\t.\t+\t.\t"
                    'gene_id "g1"; transcript_id "g1_t1"; exon_number "3";'
                ),
                (
                    "chr1\ttapis\texon\t100\t150\t.\t+\t.\t"
                    'gene_id "g1"; transcript_id "g1_t2"; exon_number "1";'
                ),
                (
                    "chr1\ttapis\texon\t300\t350\t.\t+\t.\t"
                    'gene_id "g1"; transcript_id "g1_t2"; exon_number "2";'
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    out_path = tmp_path / "AS_statistics.txt"

    write_as_statistics_from_gtf(gtf_path, out_path)

    lines = out_path.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0] == "event\tcount"
    assert any(line.startswith(f"{SpliceEventKind.EXON_SKIPPING.value}\t") for line in lines[1:])
