import dataclasses
from typing import NamedTuple


class Configuration(NamedTuple):
    knockdown_factor: float
    sigma_ul: float
    ul_factor: float


class DimensionsStringer(NamedTuple):
    DIM1: float
    DIM2: float
    DIM3: float
    DIM4: float


class DimensionsPanel(NamedTuple):
    width: float
    length: float
    height: float


class DimensionsPly(NamedTuple):
    thickness: list[float]
    angles: list[float]


class MaterialProperties(NamedTuple):
    E1: float
    E2: float
    G12: float
    nu12: float
    R1t: float
    R1c: float
    R2t: float
    R2c: float
    R21: float


class ReserveFactors(NamedTuple):
    RF_FF: list[float]
    RF_IFF: list[float]
    RF_strength: list[float]

class BiaxialReserveFactors(NamedTuple):
    RF_FF: list[float]


class Stringer(NamedTuple):
    e_id: int
    strain: float
    reserve_factor: ReserveFactors
    mode: str


class Skin(NamedTuple):
    e_id: int
    sigma_1: float
    sigma_2: float
    tau: float
    reserve_factor: ReserveFactors
    mode: str


class LoadCase(NamedTuple):
    id: int
    Panels: list[Skin]
    Stringers: list[Stringer]