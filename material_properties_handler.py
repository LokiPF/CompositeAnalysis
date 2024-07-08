from configuration import MaterialProperties, DimensionsStringer


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
