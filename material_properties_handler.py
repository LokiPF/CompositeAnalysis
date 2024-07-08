from configuration import MaterialProperties, DimensionsStringer, DimensionsPanel, Configuration


def calc_homogenized_E_stringers(mat: MaterialProperties, dim: DimensionsStringer):
    """
    Function to get the homogenized E for the stringers
    :param mat:
    :param dim:
    :return: homogenized E
    """
    A1 = dim.DIM1 * dim.DIM3  # Area flange
    A2 = (dim.DIM2 - dim.DIM3) * dim.DIM4  # Area web
    return (mat.E1 * A1 + mat.E2 * A2) / (A1 + A2)


def create_material_properties():
    material_properties = MaterialProperties()
    material_properties.E1 = 130656.82
    material_properties.E2 = 10050.52
    material_properties.G12 = 5025.26
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
    return DimensionsPanel(600, 400, 8*1.104)

def create_configuration():
    return Configuration()


print(calc_homogenized_E_stringers(create_material_properties(), create_dimensions_stringer()))
