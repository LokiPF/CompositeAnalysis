from copy import deepcopy

import numpy as np

from configuration import MaterialProperties, DimensionsPly, Stringer, Panel, ReserveFactors, DimensionsPanel, \
    DimensionsStringer, LoadCase, Configuration
from excel_handler import read_excel_input, write_excel_template, parse_ADB_matrix, parse_stringer_strength, \
    parse_constants
from material_properties_handler import create_material_properties, create_dimensions_stringer, create_dimensions_panel, \
    create_configuration, print_Engineering_Constants

E_x_flange = 0
E_x_web = 0
G_xy_flange = 0
G_xy_web = 0


def calculate_Q_matrix(mat: MaterialProperties, nu21):
    Q11 = mat.E1 / (1 - mat.nu12 * nu21)
    Q22 = mat.E2 / (1 - mat.nu12 * nu21)
    Q12 = mat.nu12 * mat.E2 / (1 - mat.nu12 * nu21)
    Q66 = mat.G12
    return np.array([[Q11, Q12, 0], [Q12, Q22, 0], [0, 0, Q66]])


def calculate_Q_bar_matrix(Q, theta):
    cos_theta = np.cos(np.radians(theta))
    sin_theta = np.sin(np.radians(theta))
    cos2_theta = cos_theta ** 2
    sin2_theta = sin_theta ** 2
    cos4_theta = cos_theta ** 4
    sin4_theta = sin_theta ** 4
    sin2cos2_theta = sin2_theta * cos2_theta

    Q11, Q12, _ = Q[0]
    _, Q22, _ = Q[1]
    Q66 = Q[2, 2]

    Q11_bar = Q11 * cos4_theta + 2 * (Q12 + 2 * Q66) * sin2cos2_theta + Q22 * sin4_theta
    Q12_bar = (Q11 + Q22 - 4 * Q66) * sin2cos2_theta + Q12 * (sin4_theta + cos4_theta)
    Q22_bar = Q11 * sin4_theta + 2 * (Q12 + 2 * Q66) * sin2cos2_theta + Q22 * cos4_theta
    Q16_bar = (Q11 - Q12 - 2 * Q66) * sin_theta * (cos_theta ** 3) + (Q12 - Q22 + 2 * Q66) * (sin_theta ** 3) * cos_theta
    Q26_bar = (Q11 - Q12 - 2 * Q66) * cos_theta * (sin_theta ** 3) + (Q12 - Q22 + 2 * Q66) * (cos_theta ** 3) * sin_theta
    Q66_bar = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * sin2cos2_theta + Q66 * (sin4_theta + cos4_theta)

    Q_bar = np.array([
        [Q11_bar, Q12_bar, Q16_bar],
        [Q12_bar, Q22_bar, Q26_bar],
        [Q16_bar, Q26_bar, Q66_bar]
    ])

    return Q_bar


def calculate_ABD_matrices(Q_bars, z_values):
    n = len(Q_bars)
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))

    for k in range(n):
        z_k = z_values[k]
        z_k1 = z_values[k + 1]
        delta_z = z_k1 - z_k
        delta_z2 = z_k1 ** 2 - z_k ** 2
        delta_z3 = z_k1 ** 3 - z_k ** 3

        Q_bar = Q_bars[k]

        A += Q_bar * delta_z
        B += Q_bar * delta_z2 / 2
        D += Q_bar * delta_z3 / 3

    return deepcopy(np.round(A, 5)), deepcopy(np.round(B, 5)), (np.round(D, 5))


def calculate_z_values(ply_thicknesses):
    total_thickness = sum(ply_thicknesses)
    z_values = [-total_thickness / 2]
    for thickness in ply_thicknesses:
        z_values.append(z_values[-1] + thickness)
    return np.round(z_values, 3)


def transform_strains(theta, epsilon_x, epsilon_y, gamma_xy):
    theta_rad = np.deg2rad(theta)
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)

    T = np.array([
        [c ** 2, s ** 2, s * c],
        [s ** 2, c ** 2, -s * c],
        [-2 * s * c, 2 * s * c, c ** 2 - s ** 2]
    ])

    epsilon_xy = np.array([epsilon_x, epsilon_y, gamma_xy])
    epsilon_12 = T @ epsilon_xy

    return epsilon_12


def calculate_stresses(Q, epsilon_12):
    sigma_12 = Q @ epsilon_12
    return sigma_12


def puck_failure_criteria(sigma2, tau21, p, material_props):
    """
    Calculate the Puck failure indices for inter-fiber fracture (IFF) in a composite material.

    Args:
    sigma2 (float): Stress component sigma2.
    tau21 (float): Shear stress component tau21.
    p_pos (float): Inclination parameter for positive sigma2.
    p_neg (float): Inclination parameter for negative sigma2.
    material_props (MaterialProperties): Material properties.

    Returns:
    float: Failure index f_E^{IFF}.
    """
    # Material properties
    R_perp = material_props.R2t
    R_parallel = material_props.R21
    R_perp_neg = material_props.R2c

    # Calculate R_perp^A
    R_perp_A = R_perp_neg / (2.0 * (1.0 + p))

    # Calculate sigma_2C and tau_21C
    sigma2C = -R_perp_A
    tau21C = R_parallel * np.sqrt(1.0 + 2.0 * p)

    mode = ""
    if sigma2 >= 0:
        # Mode A
        f_E_IFF = np.sqrt(
            (tau21 / R_parallel) ** 2 +
            (1 - p * R_perp / R_parallel) ** 2 * (sigma2 / R_perp) ** 2) + p * (sigma2 / R_parallel)
        mode = "A"
    elif abs(tau21 / sigma2) >= abs(tau21C) / R_perp_A:
        # Mode B
        f_E_IFF = (1.0 / R_parallel) * (
                np.sqrt(tau21 ** 2 + (p * sigma2) ** 2) +
                p * sigma2
        )
        mode = "B"
    else:
        # Mode C
        f_E_IFF = (
                          (tau21 / (2.0 * (1.0 + p) * R_parallel)) ** 2 +
                          (sigma2 / R_perp_neg) ** 2
                  ) * (R_perp_neg / -sigma2)
        mode = "C"

    return float(f_E_IFF), mode


def calculate_RF_FF(sigma1, material_props):
    """
    Calculate the reserve factor RF_FF for fiber fracture (FF) in a composite material.

    Args:
    sigma1 (float): Stress component sigma1.
    material_props (MaterialProperties): Material properties.

    Returns:
    float: Reserve factor RF_FF.
    """
    if sigma1 > 0:
        R_parallel = material_props.R1t
    else:
        R_parallel = -material_props.R1c
    if sigma1 == 0:
        RF_FF = 0
    else:
        RF_FF = R_parallel / sigma1
    return float(RF_FF)


def task_c(mat: MaterialProperties, dim_ply: DimensionsPly, load_cases: [LoadCase]):
    column = 2
    for load_case in load_cases:
        nu21 = mat.nu12 * mat.E2 / mat.E1
        Q_matrix = calculate_Q_matrix(mat, nu21)
        # --- Calculate reserve factors for Panel ---
        for panel in load_case.PanelsLayers:
            if panel.e_id == 1:
                RF_ff = calculate_RF_FF(1.5 * panel.sig_xx, mat)
                f_iff, mode = puck_failure_criteria(1.5 * panel.sig_yy, 1.5 * panel.sig_xy, 0.25, mat)
                RF_iff = 1 / f_iff
                RF_strength = RF_ff
                if RF_iff < RF_strength:
                    RF_strength = RF_iff
                panel.reserve_factor = ReserveFactors(RF_ff, RF_iff, RF_strength)
                panel.mode = mode
                # TODO: add error message when certain constraints missed
        # --- Calculate reserve factors for stringers ---
        for stringer in load_case.Stringers:
            if stringer.e_id == 40:
                stringer_list = []
                for angle in dim_ply.angles:
                    epsilon_12 = transform_strains(angle, 1.5 * stringer.strain, 0, 0)
                    sigma_12 = calculate_stresses(Q_matrix, epsilon_12)
                    RF_ff = calculate_RF_FF(sigma_12[0], mat)
                    f_iff, mode = puck_failure_criteria(sigma_12[1], sigma_12[2], 0.25, mat)
                    if f_iff == 0:
                        RF_iff = 0
                    else:
                        RF_iff = 1 / f_iff
                    RF_strength = RF_ff
                    if RF_iff < RF_strength:
                        RF_strength = RF_iff
                    stringer.reserve_factor = ReserveFactors(RF_ff, RF_iff, RF_strength)
                    stringer.mode = mode
                    stringer_list.append(deepcopy(stringer))
                parse_stringer_strength(stringer_list, column)
                column += 6
                # TODO: add error message when certain constraints missed


def rotate_stress_90(sigma_x, sigma_y, tau_xy):
    """
    Rotate the stress tensor by 90 degrees if the highest compressive stress is sigma_y.

    Args:
    sigma_x (float): Stress component sigma_x.
    sigma_y (float): Stress component sigma_y.
    tau_xy (float): Shear stress component tau_xy.

    Returns:
    tuple: Transformed stress components (sigma_x', sigma_y', tau_xy').
    """
    if sigma_y < sigma_x:
        # 90 degree rotation transformation
        sigma_x_prime = sigma_y
        sigma_y_prime = sigma_x
        tau_xy_prime = -tau_xy
    else:
        # No transformation needed
        sigma_x_prime = sigma_x
        sigma_y_prime = sigma_y
        tau_xy_prime = tau_xy

    return sigma_x_prime, sigma_y_prime, tau_xy_prime


def avg_panel(panels: [Panel]):
    numerator_xx, numerator_yy, numerator_xy = 0, 0, 0
    for panel in panels:
        numerator_xx += panel.sigma_1
        numerator_yy += panel.sigma_2
        numerator_xy += panel.tau
    return numerator_xx / len(panels), numerator_yy / len(panels), numerator_xy / len(panels)


def avg_stringer(stringer: list[Stringer]):
    numerator = 0
    for s in stringer:
        numerator += s.stress
    return numerator / len(stringer)


def biaxial_loading_stress(a, b, t, D, sigma_x, sigma_y):
    # Extract stiffness coefficients from the D matrix
    D11, D12, D22, D66 = D[0, 0], D[0, 1], D[1, 1], D[2, 2]

    # Initialize the minimum sigma_cr with a large number
    sigma_cr = np.inf

    # Calculate alpha
    alpha = a / b
    # Calculate beta
    beta = sigma_y / sigma_x

    # Iterate over m and n ranges to find the minimum sigma_cr
    for m in range(1, 11):
        for n in range(1, 11):
            # Calculate the critical biaxial stress
            term1 = (np.pi ** 2) / (b ** 2 * t)
            term2 = 1 / ((m / alpha) ** 2 + beta * n ** 2)
            term3 = D11 * (m / alpha) ** 4
            term4 = 2 * (D12 + D66) * (m * n / alpha) ** 2
            term5 = D22 * n ** 4

            sigma_x_cr_bi = term1 * term2 * (term3 + term4 + term5)
            if sigma_x_cr_bi < sigma_cr and sigma_x_cr_bi > 0:
                sigma_cr = sigma_x_cr_bi
    return sigma_cr


def stiffness_ratio(D):
    # Extract stiffness coefficients from the D matrix
    D11, D12, D22, D66 = D[0, 0], D[0, 1], D[1, 1], D[2, 2]

    # Calculate delta
    delta = np.sqrt(D11 * D22) / (D12 + 2 * D66)

    return delta


def shear_loading_stress(dim_ply, dim_pan, mat, a, b, t, D):
    # Extract stiffness coefficients from the D matrix
    A, B, D = calc_ABD_matrix_panel(mat, dim_pan, dim_ply)
    D11, D12, D22, D66 = D[0, 0], D[0, 1], D[1, 1], D[2, 2]

    # Calculate delta
    delta = stiffness_ratio(D)

    if delta >= 1:
        term1 = 4 / (t * b ** 2)
        term2 = np.power(D11 * D22 ** 3, 0.25)
        term3 = 8.12 + 5.05 / delta
        tau_cr = term1 * term2 * term3
    else:
        term1 = 4 / (t * b ** 2)
        term2 = np.sqrt(D22 * (D12 + 2 * D66))
        term3 = 11.7 + 0.532 * delta + 0.938 * delta ** 2
        tau_cr = term1 * term2 * term3

    return tau_cr


def task_e(config: Configuration, load_cases: list[LoadCase], dim_stringers: DimensionsStringer,
           dim_panels: DimensionsPanel, mat: MaterialProperties, dim_ply: DimensionsPly):
    # --- Calculating the volumes ---
    A1 = dim_stringers.DIM1 * dim_stringers.DIM3  # Area flange
    A2 = (dim_stringers.DIM2 - dim_stringers.DIM3) * dim_stringers.DIM4  # Area web
    volume_stringers = (A1 + A2) * dim_panels.a
    volume_panels = dim_panels.b * dim_panels.a * dim_panels.t
    A, B, D = calc_ABD_matrix(mat, dim_ply, True)
    skin_A, skin_B, skin_D = calc_ABD_matrix_panel(mat, dim_panels, dim_ply)
    # --- Calculating weighted stress average
    for load_case in load_cases:
        for i in range(4):
            stringer = load_case.Stringers[i * 3: i * 3 + 3]
            panel = load_case.Panels[i * 6 + 3: i * 6 + 9]
            avg_stringer_xx = avg_stringer(stringer) * volume_stringers
            avg_sigma_panel = avg_panel(panel)[0] * volume_panels
            sigma_combined = (avg_sigma_panel + avg_stringer_xx) / (volume_panels + volume_stringers)
            sigma_cr, sigma_crippling, E = euler_johnson(config, dim_panels, dim_stringers, skin_A, skin_D, A, D)
            RF = np.abs(sigma_cr / (sigma_combined * 1.5))
            load_case.Stringers[i].RF_combined_buckling = RF
            load_case.Stringers[i].sig_combined = sigma_combined
            load_case.Stringers[i].sig_crip = sigma_crippling


def sigma_crip(config: Configuration, dim: DimensionsStringer):
    """One Edge Free"""
    # sigma_crip_flange = config.sigma_ul * 1.63 / (np.power((dim.DIM1 * 0.5 / dim.DIM3), 0.717))
    sigma_crip_web = config.sigma_ul * 1.63 / (np.power(((dim.DIM2 - dim.DIM3) / dim.DIM4), 0.717))
    # sigma_crip_tot = (sigma_crip_flange * dim.DIM1 * dim.DIM3 + sigma_crip_web * (dim.DIM2 - dim.DIM3) * dim.DIM4) / (
    #        dim.DIM1 * dim.DIM3 + (dim.DIM2 - dim.DIM3) * dim.DIM4)
    return float(sigma_crip_web)


def E_hom_avg_x(t, A, free_lateral_deformation: bool) -> float:
    if free_lateral_deformation:
        # noinspection PyTypeChecker
        return 1 / np.linalg.inv(A)[0, 0] / t
    else:
        # noinspection PyTypeChecker
        return (A[0, 0]) / t


def euler_johnson(config: Configuration, dim_panel: DimensionsPanel, dim_stringer: DimensionsStringer, A_panel, D_panel,
                  A_stringer, D_stringer):
    A_skin = dim_panel.t * dim_panel.b
    A_stringer_flange = dim_stringer.DIM1 * dim_stringer.DIM3
    A_stringer_web = (dim_stringer.DIM2 - dim_stringer.DIM3) * dim_stringer.DIM4
    A_tot = A_skin + A_stringer_flange + A_stringer_web
    z_bar_denom = A_skin * E_hom_avg_x(dim_panel.t, A_panel / 0.9, False) + A_stringer_web * E_hom_avg_x(
        dim_stringer.DIM3, A_stringer, True) + A_stringer_flange * E_hom_avg_x(
        dim_stringer.DIM3, A_stringer, False)
    z_bar_numerator = -dim_panel.t / 2 * A_skin * E_hom_avg_x(dim_panel.t, A_panel / 0.9,
                                                              False) + dim_stringer.DIM3 / 2 * A_stringer_flange * E_hom_avg_x(
        dim_stringer.DIM3, A_stringer, False) + (
                              dim_stringer.DIM3 + (
                              dim_stringer.DIM2 - dim_stringer.DIM3) / 2) * A_stringer_web * E_hom_avg_x(
        dim_stringer.DIM3, A_stringer, True)
    z_bar = z_bar_numerator / z_bar_denom
    Az2_stringer_web = A_stringer_web * np.square(
        dim_stringer.DIM3 + (dim_stringer.DIM2 - dim_stringer.DIM3) / 2 - z_bar)
    Az2_stringer_flange = A_stringer_flange * np.square(
        dim_stringer.DIM3 / 2 - z_bar)
    Az2_panel = A_skin * np.square(-dim_panel.t / 2 - z_bar)
    I_panel = np.power(dim_panel.t, 3) * dim_panel.b / 12 + Az2_panel
    I_stringer_flange = np.power(dim_stringer.DIM3, 3) * dim_stringer.DIM1 / 12 + Az2_stringer_flange
    I_stringer_web = np.power(dim_stringer.DIM2 - dim_stringer.DIM3,
                              3) * dim_stringer.DIM4 / 12 + Az2_stringer_web
    I = I_panel + I_stringer_flange + I_stringer_web
    r_gyr = np.sqrt(I / A_tot)
    sigma_crippling = sigma_crip(config, dim_stringer)
    lamda = dim_panel.a / r_gyr
    E, E_flange, E_web, E_panel, EI, parse_E_hom = calc_E(A_panel, D_panel, A_stringer, D_stringer, dim_stringer, dim_panel,
               I_stringer_flange - Az2_stringer_flange, I_stringer_web - Az2_stringer_web,
               I_panel - Az2_panel, Az2_stringer_web,
               Az2_stringer_flange, Az2_panel)
    lamda_crit = calc_lamda_crit(E, sigma_crippling)
    # print(sigma_crippling, E, lamda_crit)
    sigma_cr = calc_sigma_cr(lamda, lamda_crit, E, sigma_crippling)
    parse_constants(parse_E_hom, E_web, E_panel, EI, z_bar, r_gyr, lamda, lamda_crit)
    return sigma_cr, sigma_crippling, E


def calc_E(A_panel, D_panel, A_stringer, D_stringer, dim_stringers: DimensionsStringer, dim_panels: DimensionsPanel,
           I_stringer_flange, I_stringer_web,
           I_panel, Az2_stringer_web, Az2_stringer_flange, Az2_panel):
    global E_x_web
    global E_x_flange
    global G_xy_flange
    global G_xy_web
    #A_stringer = np.divide(A_stringer, 0.9)
    #D_panel = np.divide(D_stringer, 0.9)

    A_inv = np.linalg.inv(A_stringer)
    D_inv = np.linalg.inv(D_stringer)

    # Calculate G_xy_flange and G_xy_web
    G_xy_flange = A_stringer[2, 2] / dim_stringers.DIM3
    G_xy_web = 1 / (A_inv[2, 2] * dim_stringers.DIM3)

    E_x_b_flange = 0.9 * A_stringer[0, 0] / dim_stringers.DIM3  # 12 / (dim_stringers.DIM3 ** 3) * D[0, 0]
    E_x_b_panel = A_panel[0, 0] / dim_panels.t  # 12 / (dim_panels.t ** 3) * A[0, 0]
    E_x_b_web = 0.9 * 1 / (A_inv[0, 0] *
                           dim_stringers.DIM3)  # 12 / (D_inv[0, 0] * (dim_stringers.DIM4 ** 3))
    E_y_b_flange = 0.9 * D_stringer[0, 0] * 12 / (
            dim_stringers.DIM3 ** 3)  # 12 / (dim_stringers.DIM3 ** 3) * A_stringer[1, 1]
    E_y_b_panel = (D_panel[0, 0]) * 12 / (dim_panels.t ** 3)  # 12 / (dim_panels.t ** 3) * A_stringer[1, 1]
    E_y_b_web = E_x_b_web  # 12 / (D_inv[1, 1] * (dim_stringers.DIM4 ** 3))
    # numerator = (I_stringer_web + Az2_stringer_web * E_x_b_web + E_y_b_flange * I_stringer_flange +
    #             Az2_stringer_flange * E_x_b_flange + E_y_b_panel * I_panel + Az2_panel * E_x_b_panel)
    denominator = I_stringer_web + Az2_stringer_web + I_stringer_flange + Az2_stringer_flange + I_panel + Az2_panel
    numerator = E_y_b_flange * I_stringer_flange + E_x_b_flange * Az2_stringer_flange + E_y_b_panel  * I_panel + E_x_b_panel * Az2_panel + E_y_b_web * I_stringer_web + E_x_b_web * Az2_stringer_web
    E_y_b = numerator / denominator
    E_x_web = E_x_b_web / 0.9
    E_x_flange = E_x_b_flange / 0.9
    return E_y_b, E_x_b_flange, E_x_b_web, 12 / (dim_panels.t ** 3) * (D_panel[0, 0]), numerator, 0.9 * 12 / (dim_stringers.DIM3 ** 3) * D_stringer[0, 0]


def calc_lamda_crit(E, sigma_crip):
    lamda_crit = np.pi * np.sqrt((2 * E) / sigma_crip)
    return lamda_crit


def calc_sigma_cr(lamda, lamda_crit, E, sigma_crip):
    if lamda < lamda_crit:
        return sigma_crip - 1 / E * np.square(sigma_crip / (2 * np.pi) * lamda)
    return (np.square(np.pi) * E) / np.square(lamda)


def calc_ABD_matrix(mat: MaterialProperties, dim_ply: DimensionsPly, write_excel_flag: bool = False):
    # --- Calculate A, B, D ---
    nu21 = mat.nu12 * mat.E2 / mat.E1
    Q_matrix = calculate_Q_matrix(mat, nu21)
    z_values = calculate_z_values(dim_ply.thickness)
    Q_bars = []
    for k in range(len(dim_ply.angles)):
        Q_bars.append(calculate_Q_bar_matrix(Q_matrix, dim_ply.angles[k]))
    A, B, D = calculate_ABD_matrices(Q_bars, z_values)
    if write_excel_flag:
        parse_ADB_matrix(A, B, D)
    return A, B, D


def calc_ABD_matrix_panel(mat: MaterialProperties, dim_panel: DimensionsPanel, dim_ply: DimensionsPly,
                          write_excel_flag: bool = False):
    # --- Calculate A, B, D ---
    nu21 = mat.nu12 * mat.E2 / mat.E1
    Q_matrix = calculate_Q_matrix(mat, nu21)
    z_values = calculate_z_values(
        [dim_panel.t / 8, dim_panel.t / 8, dim_panel.t / 8, dim_panel.t / 8, dim_panel.t / 8, dim_panel.t / 8,
         dim_panel.t / 8, dim_panel.t / 8])
    Q_bars = []
    for k in range(len(dim_ply.angles)):
        Q_bars.append(calculate_Q_bar_matrix(Q_matrix, dim_ply.angles[k]))
    A, B, D = calculate_ABD_matrices(Q_bars, z_values)
    if write_excel_flag:
        parse_ADB_matrix(A, B, D)
    return 0.9 * A, 0.9 * B, 0.9 * D


def task_d(load_cases: [LoadCase], dim_stringers: DimensionsStringer,
           dim_panels: DimensionsPanel, mat: MaterialProperties, dim_ply: DimensionsPly):
    # --- Calculating the volumes ---
    A1 = dim_stringers.DIM1 * dim_stringers.DIM3  # Area flange
    A2 = (dim_stringers.DIM2 - dim_stringers.DIM3) * dim_stringers.DIM4  # Area web
    volume_stringers = (A1 + A2) * dim_panels.b
    volume_panels = dim_panels.b * dim_panels.a * dim_panels.t
    A, B, D = calc_ABD_matrix(mat, dim_ply)
    for load_case in load_cases:
        for i in range(5):
            panel = load_case.Panels[i * 6: i * 6 + 6]
            avg_sigma_panel = avg_panel(panel)
            A, B, D = calc_ABD_matrix_panel(mat, dim_panels, dim_ply)
            sigma_cr = biaxial_loading_stress(dim_panels.a, dim_panels.b, dim_panels.t, D, avg_sigma_panel[0],
                                              avg_sigma_panel[1])
            RF_biax = np.abs((sigma_cr) / (avg_sigma_panel[0] * 1.5))
            tau_cr = shear_loading_stress(dim_ply, dim_panels, mat, dim_panels.a, dim_panels.b, dim_panels.t, D)
            RF_shear = np.abs((tau_cr) / (avg_sigma_panel[2] * 1.5))
            RF_combined = np.abs(1 / (1 / RF_biax + np.square(1 / RF_shear)))
            load_case.Panels[i].RF_panel_buckling = RF_combined
            load_case.Panels[i].sig_xx_avg = avg_sigma_panel[0]
            load_case.Panels[i].sig_yy_avg = avg_sigma_panel[1]
            load_case.Panels[i].sig_xy_avg = avg_sigma_panel[2]
            load_case.Panels[i].sig_crit_shear = tau_cr
            load_case.Panels[i].sig_crit_biax = sigma_cr


if __name__ == '__main__':
    load_cases = read_excel_input()
    material_properties = create_material_properties()
    dimensions_stringer = create_dimensions_stringer()
    dimensions_pannel = create_dimensions_panel()
    dimensions_ply = DimensionsPly()
    A, B, D = calc_ABD_matrix(material_properties, dimensions_ply, dimensions_stringer)
    configuration = create_configuration()
    task_c(material_properties, dimensions_ply, load_cases)
    task_d(load_cases, dimensions_stringer, dimensions_pannel, material_properties, dimensions_ply)
    task_e(configuration, load_cases, dimensions_stringer, dimensions_pannel, material_properties, dimensions_ply)
    write_excel_template(load_cases)
    print_Engineering_Constants(E_x_web, E_x_flange, G_xy_flange, G_xy_web, dimensions_stringer)
    pass
"""# task_c  # TODO
    # task_d  # TODO
    # task_e  # TODO"""
