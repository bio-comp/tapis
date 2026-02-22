from __future__ import annotations

from enum import Enum


class GraphNodeKind(str, Enum):
    EXON = "exon"


class GraphEdgeKind(str, Enum):
    SPLICE = "splice"


class SpliceEventKind(str, Enum):
    INTRON_RETENTION = "intron_retention"
    EXON_SKIPPING = "exon_skipping"
    ALTERNATIVE_SPLICE_SITE = "alternative_splice_site"
