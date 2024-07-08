from dataclasses import dataclass
from typing import NamedTuple


class Updatable(object):
    def update(self, new):
        for key, value in new.items():
            if hasattr(self, key):
                setattr(self, key, value)


class Configuration(NamedTuple):
    knockdown_factor: float = 0.9
    sigma_ul: float = 650
    ul_factor: float = 1.5


class DimensionsStringer(NamedTuple):
    DIM1: float
    DIM2: float
    DIM3: float
    DIM4: float


class DimensionsPanel(NamedTuple):
    a: float
    b: float
    t: float


class DimensionsPly(Updatable):
    thickness: [float] = [1.104, 1.104, 1.104, 1.104, 1.104, 1.104, 1.104, 1.104]
    angles: [float] = [45, -45, 0, 90, 90, 0, -45, 45]


class MaterialProperties(Updatable):
    E_hom: float
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
    RF_FF: float
    RF_IFF: float
    RF_strength: float


@dataclass
class Stringer(Updatable):
    e_id: int
    strain: float
    stress: float
    reserve_factor: ReserveFactors = ReserveFactors(0, 0, 0)
    mode: str = ""
    sig_combined: float = 0.0
    sig_crip: float = 0.0
    RF_combined_buckling: float = 0.0

@dataclass
class Panel(Updatable):
    e_id: int
    sigma_1: float
    sigma_2: float
    tau: float
    mode: str = ""
    sig_xx_avg: float = 0.0
    sig_yy_avg: float = 0.0
    sig_xy_avg: float = 0.0
    sig_crit_shear: float = 0.0
    sig_crit_biax: float = 0.0
    RF_panel_buckling: float = 0.0

@dataclass
class PanelLayers(NamedTuple):
    e_id: int
    layer: int
    sig_xx: float
    sig_yy: float
    sig_xy: float
    reserve_factor: ReserveFactors = 0
    mode: str = ""


@dataclass
class LoadCase(Updatable):
    id: int
    PanelsLayers: [PanelLayers]
    Panels: [Panel]
    Stringers: [Stringer]


class IO(NamedTuple):
    output_file: str = './input/ASE_Project2024_task2_1_Template_3771075.xlsx'
    sheet_name_output: str = 'ASE_Project2024_task2_1_Templat'
    delimiter: str = ''
