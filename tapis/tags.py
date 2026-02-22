from __future__ import annotations

from enum import Enum
from typing import Any, Mapping

from pydantic import ValidationError

from tapis.enums import ReadEvidence
from tapis.models.read_tags import ReadTagBundle, SingleCellTagBundle


class PacBioTag(str, Enum):
    XS = "XS"
    XL = "XL"
    XR = "XR"
    XD = "XD"
    XF = "XF"
    CR = "CR"
    CB = "CB"
    UR = "UR"
    UB = "UB"
    XM = "XM"
    XC = "XC"
    XA = "XA"
    nc = "nc"
    oc = "oc"
    gp = "gp"
    nb = "nb"
    rc = "rc"
    ic = "ic"
    is_ = "is"
    XO = "XO"
    XG = "XG"
    rq = "rq"
    iz = "iz"
    it = "it"
    im = "im"


LEGACY_TAGS = {"XS", "XL", "XR", "XD", "XF"}
MODERN_TAGS = {
    "CR",
    "CB",
    "UR",
    "UB",
    "XM",
    "XC",
    "XA",
    "nc",
    "oc",
    "gp",
    "nb",
    "rc",
    "ic",
    "is",
    "XO",
    "XG",
    "rq",
    "iz",
    "it",
    "im",
}
SINGLE_CELL_TAGS = {
    "CR",
    "CB",
    "UR",
    "UB",
    "XM",
    "XC",
    "XA",
    "nc",
    "oc",
    "gp",
    "nb",
    "rc",
    "ic",
    "is",
    "XO",
    "XG",
    "rq",
    "iz",
    "it",
    "im",
}


def parse_read_tags(tags: Mapping[str, Any]) -> ReadTagBundle:
    tag_map: dict[str, Any] = {str(key): value for key, value in tags.items()}

    single_cell = _parse_single_cell_bundle(tag_map)
    return ReadTagBundle(
        XS=_as_str(tag_map.get("XS")),
        XL=_as_int(tag_map.get("XL")),
        XR=_as_int(tag_map.get("XR")),
        XD=_as_str(tag_map.get("XD")),
        XF=_as_int(tag_map.get("XF")),
        pacbio_schema_version=_infer_schema_version(tag_map),
        read_evidence=_infer_read_evidence(tag_map),
        single_cell=single_cell,
    )


def _parse_single_cell_bundle(tag_map: Mapping[str, Any]) -> SingleCellTagBundle | None:
    subset = {tag: tag_map[tag] for tag in SINGLE_CELL_TAGS if tag in tag_map}
    if not subset:
        return None
    if "is" in subset:
        subset["is_"] = subset.pop("is")
    try:
        return SingleCellTagBundle(**subset)
    except ValidationError:
        return None


def _infer_schema_version(tag_map: Mapping[str, Any]) -> str:
    if MODERN_TAGS.intersection(tag_map):
        return "modern_v3_v4"
    if LEGACY_TAGS.intersection(tag_map):
        return "legacy_v1_v2"
    return "unknown"


def _infer_read_evidence(tag_map: Mapping[str, Any]) -> ReadEvidence:
    if _as_bool(tag_map.get("XF")):
        return ReadEvidence.FULL_LENGTH

    xs = _as_str(tag_map.get("XS"))
    xl = _as_int(tag_map.get("XL"))
    xr = _as_int(tag_map.get("XR"))

    if xs == "+":
        if xr is not None and xr <= 5:
            return ReadEvidence.THREE_PRIME_ANCHORED
        if xl is not None and xl <= 5:
            return ReadEvidence.FIVE_PRIME_ANCHORED
    elif xs == "-":
        if xl is not None and xl <= 5:
            return ReadEvidence.THREE_PRIME_ANCHORED
        if xr is not None and xr <= 5:
            return ReadEvidence.FIVE_PRIME_ANCHORED

    return ReadEvidence.UNANCHORED


def _as_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _as_str(value: Any) -> str | None:
    if value is None:
        return None
    return str(value)


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y"}
    return False
