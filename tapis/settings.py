from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, Field


class SingleCellConfig(BaseModel):
    enabled: bool = False
    require_cb: bool = False
    min_gp: int = Field(default=0, ge=0, le=1)
    barcode_whitelist: Path | None = None
    umi_mode: str = "auto"


class HybridCorrectionConfig(BaseModel):
    enabled: bool = False
    method: str = "none"
    short_read_bam: Path | None = None
    min_support: int = Field(default=2, ge=1)


class TapisSettings(BaseModel):
    single_cell: SingleCellConfig = Field(default_factory=SingleCellConfig)
    hybrid_correction: HybridCorrectionConfig = Field(default_factory=HybridCorrectionConfig)
