from __future__ import annotations

from pydantic import BaseModel, ConfigDict

from tapis.enums import ReadEvidence


class SingleCellTagBundle(BaseModel):
    CR: str | None = None
    CB: str | None = None
    UR: str | None = None
    UB: str | None = None
    XM: str | None = None
    XC: str | None = None
    XA: str | None = None
    nc: int | None = None
    oc: str | None = None
    gp: int | None = None
    nb: int | None = None
    rc: int | None = None
    ic: int | None = None
    is_: int | None = None
    XO: str | None = None
    XG: str | None = None
    rq: float | None = None
    iz: int | None = None
    it: str | None = None
    im: str | None = None

    model_config = ConfigDict(extra="ignore", populate_by_name=True)


class ReadTagBundle(BaseModel):
    XS: str | None = None
    XL: int | None = None
    XR: int | None = None
    XD: str | None = None
    XF: int | None = None
    pacbio_schema_version: str = "unknown"
    read_evidence: ReadEvidence = ReadEvidence.UNANCHORED
    single_cell: SingleCellTagBundle | None = None

    model_config = ConfigDict(extra="allow")
