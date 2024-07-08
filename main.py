import numpy as np

from configuration import MaterialProperties, DimensionsPly, Stringer, Skin, ReserveFactors, DimensionsPanel, \
    DimensionsStringer, LoadCase, Configuration


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
    Q16_bar = (Q11 - Q12 - 2 * Q66) * sin_theta * cos_theta ** 3 + (Q12 - Q22 + 2 * Q66) * sin_theta ** 3 * cos_theta
    Q26_bar = (Q11 - Q12 - 2 * Q66) * cos_theta * sin_theta ** 3 + (Q12 - Q22 + 2 * Q66) * cos_theta ** 3 * sin_theta
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

    return A, B, D


def calculate_z_values(ply_thicknesses):
    total_thickness = sum(ply_thicknesses)
    z_values = [-total_thickness / 2]
    for thickness in ply_thicknesses:
        z_values.append(z_values[-1] + thickness)
    return z_values


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
    tau21C = R_parallel * np.sqrt(1.0 + 2.0 * p * (R_perp_A / R_parallel))

    mode = ""
    if sigma2 >= 0:
        # Mode A
        f_E_IFF = np.sqrt(
            (tau21 / R_parallel) ** 2 +
            (1 - p) * (R_perp / R_parallel) ** 2 * (sigma2 / R_perp) ** 2 +
            p * (sigma2 / R_parallel)
        )
        mode = "A"
    elif sigma2 < 0 and abs(sigma2 / tau21) <= abs(R_perp_A / tau21C):
        # Mode B
        f_E_IFF = (1.0 / R_parallel) * (
                np.sqrt(tau21 ** 2 + (p * sigma2) ** 2) +
                p * sigma2
        )
        mode = "B"
    else:
        # Mode C
        f_E_IFF = np.sqrt(
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
        R_parallel = material_props.R1c

    RF_FF = R_parallel / sigma1
    return float(RF_FF)


def task_c(mat: MaterialProperties, dim_ply: DimensionsPly, stringers: list[Stringer], skins: list[Skin]):
    # --- Calculate A, B, D ---
    nu21 = mat.nu12 * mat.E2 / mat.E1
    Q_matrix = calculate_Q_matrix(mat, nu21)
    z_values = calculate_z_values(dim_ply.thickness)
    Q_bars = []
    for k in range(len(dim_ply.angles)):
        Q_bars.append(calculate_Q_bar_matrix(Q_matrix, dim_ply.angles[k]))
    A, B, D = calculate_ABD_matrices(Q_bars, z_values)
    # --- Calculate reserve factors for stringers ---
    print("Calculating reserve factors for stringers...", end='')
    for stringer in stringers:
        if stringer.e_id == 40:
            for angle in dim_ply.angles:
                epsilon_12 = transform_strains(angle, stringer.strain, 0, 0)
                sigma_12 = calculate_stresses(Q_matrix, epsilon_12)
                RF_ff = calculate_RF_FF(sigma_12[0], mat)
                f_iff, mode = puck_failure_criteria(sigma_12[1], sigma_12[2], 0.25, mat)
                RF_iff = 1 / f_iff
                stringer.reserve_factor = ReserveFactors(RF_ff, RF_iff, [])
                stringer.mode = mode
                # TODO: add strengths
                # TODO: add error message when certain constraints missed
    # --- Calculate reserve factors for skin ---
    print("Calculating reserve factors for skin...", end='')
    for skin in skins:
        if skin.e_id == 1:
            for angle in dim_ply.angles:
                RF_ff = calculate_RF_FF(skin.sigma_1, mat)
                f_iff, mode = puck_failure_criteria(skin.sigma_2, skin.tau, 0.25, mat)
                RF_iff = 1 / f_iff
                skin.reserve_factor = ReserveFactors(RF_ff, RF_iff, [])
                skin.mode = mode
                # TODO: add strengths
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


def avg_panel(skins: list[Skin]):
    numerator_xx, numerator_yy, numerator_xy = 0, 0, 0
    for skin in skins:
        numerator_xx += skin.xx
        numerator_yy += skin.yy
        numerator_xy += skin.xy
    return numerator_xx / len(skins), numerator_yy / len(skins), numerator_xy / len(skins)


def avg_stringer(stringer: list[Stringer]):
    numerator = 0
    for s in stringer:
        numerator += s.stress
    return numerator / len(stringer)


def biaxial_loading_stress(a, b, t, D11, D12, D22, D66, sigma_x, sigma_y):
    sigma_cr = 9999999999999999999
    # Calculate alpha
    alpha = a / b
    # Calculate beta
    beta = sigma_y / sigma_x
    for m in range(10):
        for n in range(10):
            # Calculate the critical biaxial stress
            term1 = (np.pi ** 2) / (b ** 2 * t)
            term2 = 1 / ((m / alpha) ** 2 + beta * n ** 2)
            term3 = D11 * (m / alpha) ** 4
            term4 = 2 * (D12 + D66) * (m * n / alpha) ** 2
            term5 = D22 * n ** 4

            sigma_x_cr_bi = term1 * term2 * (term3 + term4 + term5)
            if sigma_x_cr_bi < sigma_cr:
                sigma_cr = sigma_x_cr_bi

    return sigma_cr


def stiffness_ratio(D11, D12, D22, D66):
    # Calculate delta
    delta = np.sqrt(D11 * D22) / (D12 + 2 * D66)
    return delta


def shear_loading_stress(a, b, t, D11, D12, D22, D66):
    # Calculate delta
    delta = stiffness_ratio(D11, D12, D22, D66)

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


def averaged_stress(load_cases: list[LoadCase], dim_stringers: DimensionsStringer,
                    dim_panels: DimensionsPanel):
    # --- Calculating the volumes ---
    A1 = dim_stringers.DIM1 * dim_stringers.DIM3  # Area flange
    A2 = (dim_stringers.DIM2 - dim_stringers.DIM3) * dim_stringers.DIM4  # Area web
    volume_stringers = (A1 + A2) * dim_panels.length
    volume_panels = dim_panels.length * dim_panels.width * dim_panels.height
    # --- Calculating weighted stress average
    for load_case in load_cases:
        for i in range(4):
            stringer = load_case.Stringers[i * 3: i * 3 + 3]
            panel = load_case.Panels[i * 6 + 3: i * 6 + 9]
            avg_stringer_xx = avg_stringer(stringer) * volume_stringers
            avg_sigma_panel = avg_panel(panel)[0] * volume_panels
            sigma_combined = (avg_sigma_panel + avg_stringer_xx) / (volume_panels + volume_stringers)
            sigma_cr = euler_johnson(dim_panels, dim_stringers)
            RF = (sigma_cr * 0.9) / (sigma_combined * 1.5) # TODO


def sigma_crip_flange(config: Configuration, dim: DimensionsStringer):
    """One Edge Free"""
    sigma_crip_flange_left = config.sigma_ul * 1.63 / (np.power((dim.DIM1 / dim.DIM3), 0.717))
    sigma_crip_web = config.sigma_ul * 1.63 / (np.power(((dim.DIM2 - dim.DIM3) / dim.DIM4), 0.717))
    sigma_crip_tot = (sigma_crip_flange * dim.DIM1 * dim.DIM3 + sigma_crip_web * (dim.DIM2 - dim.DIM3) * dim.DIM4) / (
            dim.DIM1 * dim.DIM3 + (dim.DIM2 - dim.DIM3) * dim.DIM4)
    return sigma_crip_tot


def euler_johnson(dim_panel: DimensionsPanel, dim_stringer: DimensionsStringer):
    A_skin = dim_panel.height * dim_panel.width
    A_stringer_flange = dim_stringer.DIM1 * dim_stringer.DIM3
    A_stringer_web = (dim_stringer.DIM2 - dim_stringer.DIM3) * dim_stringer.DIM4
    A_tot = A_skin + A_stringer_flange + A_stringer_web
    z_bar_numerator = -dim_panel.height / 2 * A_skin + dim_stringer.DIM3 / 2 * A_stringer_flange + (
            dim_stringer.DIM3 + (dim_stringer.DIM2 - dim_stringer.DIM3) / 2) * A_stringer_web
    z_bar = z_bar_numerator / A_tot
    I_skin = np.power(dim_panel.height, 3) * dim_panel.width / 12 + A_skin * np.square(-dim_panel.height / 2 - z_bar)
    I_stringer_flange = np.power(dim_stringer.DIM3, 3) * dim_stringer.DIM1 / 12 + A_stringer_flange * np.square(
        dim_stringer.DIM1 / 2 - z_bar)
    I_stringer_web = np.power(dim_stringer.DIM2 - dim_stringer.DIM3,
                              3) * dim_stringer.DIM4 / 12 + A_stringer_web * np.square(
        dim_stringer.DIM3 + (dim_stringer.DIM2 - dim_stringer.DIM3) / 2 - z_bar)
    I = I_skin + I_stringer_flange + I_stringer_web
    r_gyr = np.sqrt(I / A_tot)
    lamda = dim_panel.length / r_gyr
    sigma_cr = calc_sigma_cr() # TODO
    return sigma_cr


def calc_lamda_crit(E, sigma_crip):
    lamda_crit = np.sqrt((2 * np.pi ** 2 * E) / sigma_crip)
    return lamda_crit


def calc_sigma_cr(lamda, lamda_crit, E, sigma_crip):
    if lamda < lamda_crit:
        return sigma_crip - 1 / E * np.square(sigma_crip / (2 * np.pi) * lamda)
    return (np.square(np.pi) * E) / np.square(lamda)


def task_d(load_cases: list[LoadCase], dim_stringers: DimensionsStringer,
           dim_panels: DimensionsPanel):
    # --- Calculating the volumes ---
    A1 = dim_stringers.DIM1 * dim_stringers.DIM3  # Area flange
    A2 = (dim_stringers.DIM2 - dim_stringers.DIM3) * dim_stringers.DIM4  # Area web
    volume_stringers = (A1 + A2) * dim_panels.length
    volume_panels = dim_panels.length * dim_panels.width * dim_panels.height
    for load_case in load_cases:
        for i in range(5):
            panel = load_case.Panels[i * 6 + 3: i * 6 + 9]
            avg_sigma_panel = avg_panel(panel)
            sigma_cr = biaxial_loading_stress()  # TODO
            RF_biax = (sigma_cr * 0.9) / (avg_sigma_panel[0] * 1.5)
            tau_cr = shear_loading_stress()  # TODO
            RF_shear = (tau_cr * 0.9) / (avg_sigma_panel[1] * 1.5)
            RF_combined = np.abs(1 / (1 / RF_biax + np.square(1 / RF_shear)))
