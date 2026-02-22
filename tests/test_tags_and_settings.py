from __future__ import annotations

from tapis.enums import AnchorCigarOp, SamCigarOp
from tapis.settings import TapisSettings
from tapis.tags import ReadEvidence, parse_read_tags


def test_cigar_enums_match_sam_codes() -> None:
    assert SamCigarOp.MATCH.value == 0
    assert SamCigarOp.DIFF.value == 8
    assert AnchorCigarOp.AMATCHL.value == 101
    assert AnchorCigarOp.ADELETE.value == 12


def test_parse_legacy_pacbio_tags() -> None:
    bundle = parse_read_tags({"XS": "+", "XL": 2, "XR": 3, "XD": "GT-AG"})
    assert bundle.pacbio_schema_version == "legacy_v1_v2"
    assert bundle.read_evidence == ReadEvidence.THREE_PRIME_ANCHORED
    assert bundle.single_cell is None


def test_parse_modern_single_cell_tags() -> None:
    bundle = parse_read_tags(
        {
            "CB": "AACCTT",
            "CR": "AACCTT",
            "UB": "GGTTA",
            "UR": "GGTTA",
            "gp": 1,
            "nb": 0,
            "rc": 1,
            "XM": "GGTTA",
            "rq": 0.99,
            "iz": 12,
        }
    )
    assert bundle.pacbio_schema_version == "modern_v3_v4"
    assert bundle.single_cell is not None
    assert bundle.single_cell.CB == "AACCTT"


def test_settings_defaults_keep_single_cell_off() -> None:
    settings = TapisSettings()
    assert settings.single_cell.enabled is False
    assert settings.hybrid_correction.enabled is False
