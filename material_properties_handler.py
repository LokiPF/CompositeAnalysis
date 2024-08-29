import numpy as np

from configuration import MaterialProperties, DimensionsStringer, DimensionsPanel, Configuration


def calc_homogenized_E_stringers(mat: MaterialProperties, dim: DimensionsStringer, A, t):
    """
    Function to get the homogenized E for the stringers
    :param mat:
    :param dim:
    :return: homogenized E
    """
    return A[0,0] / t
    A1 = dim.DIM1 * dim.DIM3  # Area flange
    A2 = (dim.DIM2 - dim.DIM3) * dim.DIM4  # Area web
    return (mat.E1 * A1 + mat.E2 * A2) / (A1 + A2)


def create_material_properties():
    material_properties = MaterialProperties()
    material_properties.E1 = 133382.5
    material_properties.E2 = 10260.19
    material_properties.G12 = 5130.1
    material_properties.nu12 = 0.33
    material_properties.R1t = 3050
    material_properties.R1c = 1500
    material_properties.R2t = 300
    material_properties.R2c = 50
    material_properties.R21 = 100
    return material_properties


def create_dimensions_stringer():
    return DimensionsStringer(70, 42.944, 2.944, 2.944)


def create_dimensions_panel():
    return DimensionsPanel(600, 400, 8 * 1.104)


def create_configuration():
    return Configuration()

def print_Engineering_Constants(E_x_flange, E_x_web, G_xy_flange, G_xy_web, dim):
    A1 = dim.DIM1 # Area Web
    A2 = (dim.DIM2 - dim.DIM4) # Area flange
    print(f"E: {(A2 * E_x_flange + A1 * E_x_web) / (A1 + A2)}")
    print(f"G: {(A2 * G_xy_flange + A1 * G_xy_web) / (A1 + A2)}")

