from __future__ import annotations

from enum import Enum, IntEnum


class SamCigarOp(IntEnum):
    MATCH = 0
    INSERT = 1
    DELETE = 2
    GAP = 3
    SOFT_CLIP = 4
    HARD_CLIP = 5
    PAD = 6
    EQUAL = 7
    DIFF = 8


class AnchorCigarOp(IntEnum):
    AINSERT = 11
    ADELETE = 12
    AMATCHL = 101
    AMATCHR = 102


class ReadEvidence(str, Enum):
    FULL_LENGTH = "full_length"
    THREE_PRIME_ANCHORED = "three_prime_anchored"
    FIVE_PRIME_ANCHORED = "five_prime_anchored"
    UNANCHORED = "unanchored"
