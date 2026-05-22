import numpy as np
import numpy.typing as npt
import pandas as pd
import math, cmath
import scipy.special as special
import scipy.constants as conts
from typing import Any, Union, overload
import mpmath

# distutils: language=Py3

a0 = 0.529177210903e-10  # in m
el_charge = conts.e  # in C
h = conts.h / el_charge  # in eV*s
Ryd_ener = 13.605693122990  # in eV
alpha = conts.fine_structure  # Sommerfeld fine structure constant
alpha_sq = pow(alpha, 2)
el_mass = conts.electron_mass  # in kg
prefactor = pow(el_charge, 2) / (2 * pow(math.pi, 2) * el_mass)
speed_of_light = conts.c  # m/s
r_0 = pow(el_charge, 2) / el_charge / pow(speed_of_light, 2)
r_0_for_nu = pow(el_charge, 2) / el_charge
angstrom2eV = 1.23984193 * 10000  # eV*µm * µm/Angstrom
barn2bohr = 2.80028520539078e7
constant_factor = 4 * math.pi * math.pi * el_mass / h / h  # p_0 according to Hönl's Waller paper (1933) Page 646
ipi = math.pi * 1j
point_5_pi = mpmath.mpc(0.0, 0.5) * mpmath.pi
half_plus_i = mpmath.mpc(0.5, 1.0)
two_thirds = mpmath.mpf(2.0)/mpmath.mpf(3.0)
ten_thirds = mpmath.mpf(10.0)/mpmath.mpf(3.0)

foa = Union[float, npt.NDArray]
ioa = Union[int, npt.NDArray]

Zeff_data = pd.read_csv("zeff.csv")
l_dict = {0: "s", 1: "p", 2: "d", 3: "f"}


def dummy_func(*args, **kwargs) -> Any:
    pass


@overload
def sugiura_exps(z: float, n_0: int) -> float: ...


@overload
def sugiura_exps(z: npt.NDArray, n_0: int) -> npt.NDArray: ...


def sugiura_exps(z: Union[float, npt.NDArray], n_0: int) -> Union[float, npt.NDArray]:
    z_1 = z - 1.0
    if z <= 1.0:
        temp = np.arctan(np.sqrt(1-z)) / (np.sqrt(1-z)) - 2.0 / 3.0 * (z-1) * special.hyp2f1(0.75, 1, 1.75, (z-1)**2)
        part2 = np.exp(-4 * n_0 * temp)
        part3 = 1.0
    else:
        sqrt_var = np.sqrt(z_1)
        part2 = np.exp(-4 * n_0 / sqrt_var * np.arctan(sqrt_var))
        part3 = 1 - np.exp(-2 * n_0 * math.pi / sqrt_var)
    return part2 / part3
    # The one above is with the correction for the part when z is <1 where we continue the function using the hyp2f1
    return np.exp(-4 * n_0 / np.sqrt(z - 1) * np.arctan(np.sqrt(z - 1))) / (1 - np.exp(-2 * n_0 * math.pi / np.sqrt(z - 1)))


@overload
def f_s_1_hoenl(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_1_hoenl(z: float) -> float: ...


def f_s_1_hoenl(z: foa) -> foa:
    part1 = 64 / (3 * np.power(z, 3))
    part2 = sugiura_exps(z, 1)
    return part1 * part2


@overload
def f_s_2_1_hoenl(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_2_1_hoenl(z: float) -> float: ...


def f_s_2_1_hoenl(z: foa) -> foa:
    part1 = 256 * (z - 2) / (15 * z**5)
    part2 = sugiura_exps(z, 1)
    return part1 * part2


@overload
def f_s_2_2_hoenl(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_2_2_hoenl(z: float) -> float: ...


def f_s_2_2_hoenl(z: foa) -> foa:
    part1 = 256 * (4 * z - 3) / (15 * z**5)
    part2 = sugiura_exps(z, 1)
    return part1 * part2


@overload
def f_s_1_EM(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_1_EM(z: float) -> float: ...


def f_s_1_EM(z: foa) -> foa:
    part1 = 512 * (z + 3) / (3 * z**4)
    part2 = sugiura_exps(z, 2)
    return part1 * part2


@overload
def f_s_2_1_EM(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_2_1_EM(z: float) -> float: ...


def f_s_2_1_EM(z: foa) -> foa:
    part1 = 512 * (z - 2) * (z + 3) / (15 * z**6)
    part2 = sugiura_exps(z, 2)
    return part1 * part2


@overload
def f_s_2_2_EM(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_2_2_EM(z: float) -> float: ...


def f_s_2_2_EM(z: foa) -> foa:
    part1 = 2048 * pow(z - 1, 2) * (z + 3) / (15 * z**7)
    part2 = sugiura_exps(z, 2)
    return part1 * part2


@overload
def f_p_1_EM(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_p_1_EM(z: float) -> float: ...


def f_p_1_EM(z: foa) -> foa:
    part1 = 512 * (z + 8.0 / 3.0) / (3 * z**5)
    part2 = sugiura_exps(z, 2)
    return part1 * part2


@overload
def f_p_2_1_EM(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_p_2_1_EM(z: float) -> float: ...


def f_p_2_1_EM(z: foa) -> foa:
    return (z - 2) / z**2 / 5 * f_p_1_EM(z)


@overload
def f_p_2_2_EM(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_p_2_2_EM(z: float) -> float: ...


def f_p_2_2_EM(z: foa) -> foa:
    return 2 * (11 * z - 6) * (z + 3) / (z + 8.0 / 3.0) / z**2 / 15 * f_p_1_EM(z)


@overload
def f_s_1_WA(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_1_WA(z: float) -> float: ...


def f_s_1_WA(z: foa) -> foa:
    part1 = 64 * (3 * z + 4) ** 2 * (z + 8) / z**6
    part2 = sugiura_exps(z, 3)
    return part1 * part2


@overload
def f_s_2_1_WA(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_2_1_WA(z: float) -> float: ...


def f_s_2_1_WA(z: foa) -> foa:
    part1 = 256 * pow(3 * z + 4, 2) * (z + 8) * (z - 2) / (45 * z**8)
    part2 = sugiura_exps(z, 3)
    return part1 * part2


@overload
def f_s_2_2_WA(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_s_2_2_WA(z: float) -> float: ...


def f_s_2_2_WA(z: foa) -> foa:
    part1 = 256 * (3 * z - 4) * (z + 8) * (4 * z + 5) * (3 * z - 4) / (45 * z**8)
    part2 = sugiura_exps(z, 3)
    return part1 * part2


@overload
def f_p_1_WA(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_p_1_WA(z: float) -> float: ...


def f_p_1_WA(z: foa) -> foa:
    part1 = 512 * (3 * z**2 + 26 * z + 28) / z**6
    part2 = sugiura_exps(z, 3)
    return part1 * part2


@overload
def f_d_1_WA(z: npt.NDArray) -> npt.NDArray: ...


@overload
def f_d_1_WA(z: float) -> float: ...


def f_d_1_WA(z: foa) -> foa:
    part1 = 1024 / 5.0 * (5 * z**2 + 46 * z + 48) / (z**7)
    part2 = sugiura_exps(z, 3)
    return part1 * part2


#def Rnl(r: np.ndarray, b: float, l: int, n: int) -> np.ndarray:
#    part1 = np.sqrt((2 * b) ** 3 * (special.factorial(n - l - 1)) / (2 * n * (special.factorial(n + l))))
#    part2 = np.exp(-b * r)
#    part3 = np.power(2 * b * r, l)
#    part4 = special.genlaguerre(n - l - 1, 2 * l + 1)(2 * b * r)
#    return part1 * part2 * part3 * part4


# 2pi m/s^2
@overload
def q(nu: float) -> float: ...


@overload
def q(nu: npt.NDArray) -> npt.NDArray: ...


def q(nu: foa) -> foa:
    return 2.0 * math.pi * nu / speed_of_light


def b(n_0: int, l_0: int, Z: int, j: str) -> float:
    Z_eff = get_Zeff(Z, n_0, l_0, j)
    return Z_eff / (n_0 * a0)


def get_Zeff(Z: int, n_0: int, l_0: int, j: str):
    if l_0 != 0:
        search_string = f"{n_0}{l_dict[l_0]}{j}"
    else:
        search_string = f"{n_0}{l_dict[l_0]}"
    
    # Check if Z exists in the dataframe
    if Z not in Zeff_data["Z"].values:
        raise ValueError(f"Z={Z} not found in Zeff_data")
    
    # Check if the column exists
    if search_string not in Zeff_data.columns:
        raise ValueError(f"Column '{search_string}' not found in Zeff_data")
    
    Z_eff = Zeff_data.loc[Zeff_data["Z"] == Z, search_string].values
    # Check if any value was found
    if len(Z_eff) == 0:
        raise ValueError(f"No value found for Z={Z} and column '{search_string}'")
    # Get the value
    return float(Z_eff[0])
  
def Energie(Z: int, n0: int, l: int, j: float) -> float:
    # this is equation (86) in the paper for E_{n_0,l,j}(Z)
    if j == 0.5:
      Z_eff = get_Zeff(Z, n0, l, "1/2")
    elif j == 1.5:
      Z_eff = get_Zeff(Z, n0, l, "3/2")
    elif j == 2.5:
      Z_eff = get_Zeff(Z, n0, l, "5/2")
    res = (
        -2
        / alpha_sq
        * (
            1
            - np.power(
                1
                + np.power(
                    alpha * Z_eff / (n0 - j - 0.5 + np.power((j + 0.5) ** 2 - alpha_sq * Z_eff ** 2, 0.5)),
                    2,
                ),
                -0.5,
            )
        )
    )
    return res * Ryd_ener
    

@overload
def xn(nu_in: float, n_0: int, el_nr: int, l_0: int, n: int) -> float: ...


@overload
def xn(nu_in: npt.NDArray, n_0: int, el_nr: int, l_0: int, n: int) -> npt.NDArray: ...


def xn(nu_in: foa, n_0: int, el_nr: int, l_0: int, n: int) -> foa:
    return pow(n_0 * q(nu_in) / b(n_0, l_0, el_nr), n)


@overload
def x2(nu_in: float, n_0: int, el_nr: int, l_0: int) -> float: ...


@overload
def x2(nu_in: npt.NDArray, n_0: int, el_nr: int, l_0: int) -> npt.NDArray: ...


def x2(nu_in: foa, n_0: int, el_nr: int, l_0: int) -> foa:
    return xn(nu_in, n_0, el_nr, l_0, 2)


def z_EE(E: float, E2: float) -> float:
    return (E + abs(E2)) / abs(E2)


def z_kb(k: float, b: float) -> float:
    return pow(k, 2) / pow(b, 2) + 1


def z_nprime(n_prime: float, n_0: int) -> float:
    return n_0 * n_0 / pow(n_prime, 2) + 1


def n_prime_from_z(z: float, n_0: int) -> complex:
    return n_0 / cmath.sqrt(z - 1)


def z_nunu(nu_j: float, nu_2: float) -> float:
    return nu_j / nu_2


def k(E: float) -> float:
    return 2 * math.pi / h * math.sqrt(2 * el_mass * E)


def n_prime(E: float, Z: int, n: int, l: int) -> float:
    return 2 * b(n, l, Z) / k(E)


# introduces m^3 / pi
def N0_square(b_: float) -> float:
    return pow(b_, 3) / math.pi


def N0(b_: float) -> float:
    return math.sqrt(N0_square(b_))


def product_n_prime_from_z(n_0: int, z: float, l: int) -> float:
    n_p = n_prime_from_z(z, n_0)
    fact = 1.0
    for nu in range(1, l + 1):
        fact *= (n_p * n_p).real + nu * nu
    denom = 1 - math.exp(-2 * math.pi * abs(n_p))
    return fact / denom


# Introduces factors 2pi m_e /h^2 and m
def N_square_from_z(l: int, m: int, b_: float, n_0: int, z: float) -> float:
    if m > l:
        return 0
    result = (2 * l + 1) * math.factorial(l - m) / math.factorial(l + m) * constant_factor / math.pi * n_0 * b_ * product_n_prime_from_z(n_0, z, l)
    if m >= 1:
        result *= 2
    return result


def N(l: int, m: int, b_: float, n_0: int, z: float) -> float:
    return math.sqrt(N_square_from_z(l, m, b_, n_0, z))


def N_square(l: int, m: int, b_: float, n_0: int, z: float) -> float:
    return N_square_from_z(l, m, b_, n_0, z)


def N_lm_from_z(l: int, m: int, z: float, b_: float, n_0: int) -> float:
    return N(l, m, b_, n_0, z)


def N_lm_square_from_z(l: int, m: int, z: float, b_: float, n_0: int) -> float:
    return N_square(l, m, b_, n_0, z)


kpcor = [
    0.0,
    0.0,
    0.0,
    1e-3,
    0.0,
    1e-3,
    1e-3,
    2e-3,
    3e-3,
    4e-3,
    4e-3,
    6e-3,
    8e-3,
    8e-3,
    11e-3,
    12e-3,
    14e-3,
    17e-3,
    20e-3,
    22e-3,
    25e-3,
    28e-3,
    31e-3,
    35e-3,
    39e-3,
    42e-3,
    48e-3,
    52e-3,
    57e-3,
    61e-3,
    67e-3,
    73e-3,
    79e-3,
    85e-3,
    92e-3,
    99e-3,
    106e-3,
    114e-3,
    122e-3,
    130e-3,
    138e-3,
    147e-3,
    156e-3,
    166e-3,
    175e-3,
    186e-3,
    196e-3,
    207e-3,
    219e-3,
    230e-3,
    242e-3,
    255e-3,
    267e-3,
    281e-1,
    294e-3,
    308e-3,
    323e-3,
    338e-3,
    354e-3,
    369e-3,
    386e-3,
    402e-3,
    419e-3,
    437e-3,
    455e-3,
    474e-3,
    493e-3,
    512e-3,
    532e-3,
    553e-3,
    574e-3,
    596e-3,
    617e-3,
    640e-3,
    663e-3,
    687e-3,
    711e-3,
    736e-3,
    762e-3,
    799e-3,
    814e-3,
    842e-3,
    870e-3,
    899e-3,
    928e-3,
    957e-3,
    988e-3,
    1018e-3,
    1050e-3,
    1083e-3,
    1115e-3,
    1149e-3,
    1184e-3,
]
relcor = [
    0.0,
    0.0,
    0.0,
    1.000000047e-03,
    1.000000047e-03,
    2.000000095e-03,
    3.000000026e-03,
    4.999999888e-03,
    7.000000216e-03,
    8.999999613e-03,
    1.099999994e-02,
    1.400000043e-02,
    1.799999923e-02,
    2.099999972e-02,
    2.600000054e-02,
    2.999999933e-02,
    3.500000015e-02,
    4.100000113e-02,
    4.699999839e-02,
    5.299999937e-02,
    5.999999866e-02,
    6.800000370e-02,
    7.500000298e-02,
    8.399999887e-02,
    9.300000221e-02,
    0.101999998,
    0.112999998,
    0.123000003,
    0.135000005,
    0.145999998,
    0.158999994,
    0.172000006,
    0.186000004,
    0.200000003,
    0.215000004,
    0.231000006,
    0.246999994,
    0.263999999,
    0.282000005,
    0.300000012,
    0.319000006,
    0.338000000,
    0.358999997,
    0.379999995,
    0.400999993,
    0.423999995,
    0.446999997,
    0.470999986,
    0.495999992,
    0.521000028,
    0.546999991,
    0.574999988,
    0.601999998,
    0.630999982,
    0.660000026,
    0.689999998,
    0.721000016,
    0.753000021,
    0.786000013,
    0.819000006,
    0.853999972,
    0.888999999,
    0.925000012,
    0.962000012,
    1.00000000,
    1.03900003,
    1.07900000,
    1.11899996,
    1.16100001,
    1.20400000,
    1.24800003,
    1.29299998,
    1.33800006,
    1.38499999,
    1.43299997,
    1.48199999,
    1.53199995,
    1.58299994,
    1.63600004,
    1.68900001,
    1.74300003,
    1.79900002,
    1.85599995,
    1.91400003,
    1.97300005,
    2.03299999,
    2.09500003,
    2.15700006,
    2.22099996,
    2.28699994,
    2.35299993,
    2.42100000,
    2.49000001,
]


def exp_from_z(z: float, n_0: int) -> float:
    z_1 = z - 1
    temp = 0.0
    sqrt_var = np.sqrt(abs(z_1))
    if z <= 1:
        temp -= 2.0 / 3.0 * z_1 * special.hyp2f1(0.75, 1, 1.75, z_1 * z_1)
    temp += np.arctan(sqrt_var) / sqrt_var
    return math.exp(-2 * n_0 * temp)
    return math.exp(-2 * n_0 / math.sqrt(z - 1) * math.atan(math.sqrt(z - 1)))


def exp_squared_from_z(z: float, n_0: int) -> float:
    return exp_from_z(z, n_0) ** 2
    return math.exp(-4 * n_0 / math.sqrt(z - 1) * math.atan(math.sqrt(z - 1)))


def kronecke_delta(a: Union[int, float], b: Union[int, float]) -> int:
    if a == b:
        return 1
    else:
        return 0


def one_minus_delta_edge(Z: int, n0: int, l: int, j: float) -> float:
    if n0 == 1:
        if l != 0:
            raise ValueError("unknown orbital type!")
        else:
            Z_s_sq = pow(get_Zeff_1s(Z), 2)
            e_ion = get_ionization_energy_1s(Z)
            delta_K = alpha_sq * Z_s_sq / 4 - e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
            return -delta_K
    elif n0 == 2:
        if l == 0:
            Z_s_sq = pow(get_Zeff_2s(Z), 2)
            e_ion = get_ionization_energy_2s(Z)
            delta_l1 = alpha_sq * Z_s_sq * 0.3125 - 4 * e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
            return -delta_l1
        elif l == 1:
            if j == 0.5:
                Z_s_sq = pow(get_Zeff_2p_1_2(Z), 2)
                e_ion = get_ionization_energy_2p1_2(Z)
                delta_l2 = alpha_sq * Z_s_sq * 0.3125 - 4 * e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
                return -delta_l2
            elif j == 1.5:
                Z_s_sq = pow(get_Zeff_2p_3_2(Z), 2)
                e_ion = get_ionization_energy_2p3_2(Z)
                delta_l2 = alpha_sq * Z_s_sq * 0.0625 - 4 * e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
                return -delta_l2
            else:
                raise ValueError("unknown orbital type!")
        else:
            raise ValueError("unknown orbital type!")
    elif n0 == 3:
        if l == 0:
            Z_s_sq = pow(get_Zeff_3s(Z), 2)
            e_ion = get_ionization_energy_3s(Z)
            delta_m1 = alpha_sq * Z_s_sq * 0.25 - 9 * e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
            return -delta_m1
        elif l == 1:
            if j == 0.5:
                Z_s_sq = pow(get_Zeff_3p_1_2(Z), 2)
                e_ion = get_ionization_energy_3p_1_2(Z)
                delta_m2 = alpha_sq * Z_s_sq * 0.25 - 9 * e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
                return -delta_m2
            elif j == 1.5:
                Z_s_sq = pow(get_Zeff_3p_3_2(Z), 2)
                e_ion = get_ionization_energy_3p_3_2(Z)
                delta_m3 = alpha_sq * Z_s_sq / 12 - 9 * e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
                return -delta_m3
            else:
                raise ValueError("unknown orbital type!")
        elif l == 2:
            if j == 1.5:
                Z_s_sq = pow(get_Zeff_3d_3_2(Z), 2)
                e_ion = get_ionization_energy_3d_3_2(Z)
                delta_m4 = alpha_sq * Z_s_sq / 12 - 9 * e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
                return -delta_m4
            elif j == 2.5:
                Z_s_sq = pow(get_Zeff_3d_5_2(Z), 2)
                e_ion = get_ionization_energy_3d_5_2(Z)
                delta_m5 = alpha_sq * Z_s_sq / 36 - 9 * e_ion / (Ryd_ener * Z_s_sq)  # 21a in Hoenl
                return -delta_m5
            else:
                raise ValueError("unknown orbital type!")
        else:
            raise ValueError("unknown orbital type!")
    else:
        raise ValueError("unknown orbital type!")


elements = [
    "DUMMY",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
]


def get_ionization_energy_1s(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        13.6,
        24.6,
        54.7,
        111.5,
        188.0,
        284.2,
        409.9,
        543.1,
        696.7,
        870.2,
        1070.8,
        1303.0,
        1559.6,
        1839.0,
        2145.5,
        2472.0,
        2822.4,
        3205.9,
        3608.4,
        4038.5,
        4492,
        4966,
        5465,
        5989,
        6539,
        7112,
        7709,
        8333,
        8979,
        9659,
        10367,
        11103,
        11867,
        12658,
        13474,
        14326,
        15200,
        16105,
        17038,
        17998,
        18986,
        20000,
        21044,
        22117,
        23220,
        24350,
        25514,
        26711,
        27940,
        29200,
        30491,
        31814,
        33169,
        34561,
        35985,
        37441,
        38925,
        40443,
        41991,
        43569,
        45184,
        46834,
        48519,
        50239,
        51996,
        53789,
        55618,
        57486,
        59390,
        61332,
        63314,
        65351,
        67416,
        69525,
        71676,
        73871,
        76111,
        78395,
        80725,
        83102,
        85530,
        88005,
        90524,
        93105,
        95730,
        98404,
        101137,
        103922,
        106755,
        109651,
        112601,
        115606,
    ]
    if Z <= 92:
        return ionization_energies[Z - 1]
    else:
        # extrapoalte from fourth order regression of tabulated data
        a1 = 4.00729846e-04
        a2 = -2.27683735e-02
        a3 = 1.30249774e01
        a4 = -6.43729353e01
        c = 1.99592659e02
        return float(a1 * pow(Z, 4) + a2 * pow(Z, 3) + a3 * pow(Z, 2) + a4 * Z + c)
    import np as np
    from scipy.optimize import curve_fit

    def func(x, a1, a2, a3, a4, c):
        return a1 * x**4 + a2 * x**3 + a3 * x**2 + a4 * x + c

    x = range(2, 93)
    result, resultcov = curve_fit(func, x, ionization_energies[1:])
    print(result)
    print(resultcov)
    coef = np.polyfit(x, ionization_energies[1:], 3)
    funct = np.polyval(coef, Z)
    print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, ionization_energies[1:], marker="*")
    plt.plot(x, func(x, *result), "r-")
    plt.show()


def get_ionization_energy_2s(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.00,
        37.3,
        41.6,
        0.000,
        48.5,
        63.5,
        88.7,
        117.8,
        149.7,
        189.0,
        230.9,
        270.0,
        326.3,
        378.6,
        438.4,
        498.0,
        560.9,
        626.7,
        696.0,
        769.1,
        844.6,
        925.1,
        1008.6,
        1096.7,
        1196.2,
        1299.0,
        1414.6,
        1527.0,
        1652.0,
        1782,
        1921,
        2065,
        2216,
        2373,
        2532,
        2698,
        2866,
        3043,
        3224,
        3412,
        3604,
        3806,
        4018,
        4238,
        4465,
        4698,
        4939,
        5188,
        5453,
        5714,
        5989,
        6266,
        6549,
        6835,
        7126,
        7428,
        7737,
        8052,
        8376,
        8708,
        9046,
        9394,
        9751,
        10116,
        10486,
        10870,
        11271,
        11682,
        12100,
        12527,
        12968,
        13419,
        13880,
        14353,
        14839,
        15347,
        15861,
        16388,
        16939,
        17493,
        18049,
        18639,
        19237,
        19840,
        20472,
        21105,
        21757,
    ]
    # Coefficients from fit to data
    a1 = 1.44070614e01
    a2 = 1.49753416e-01
    a3 = 1.15499698e-02
    a4 = 1.61081136e00
    a5 = -2.27553738e01
    b1 = 4.99999128e01
    c = 1.27894569e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 3:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapolate from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
        import np as np
        from scipy.optimize import curve_fit

        def func(x, a1, a2, a3, a4, a5, b1, c):
            return a1 * np.exp2(a2 * (x - b1)) + a3 * x * x * x + a4 * x * x + a5 * x + c

        x = range(10, 93)
        x2 = range(1, 103)
        y = ionization_energies[9:]
        result, resultcov = curve_fit(
            func,
            x,
            y,
            p0=[8.36040351e02, 5.75329301e-02, 2.15148012e00, 0, 0, 2.51302235e01, 0],
            bounds=([0, 0, -np.inf, -np.inf, -np.inf, -50, -500], [2000, 1000, 1000, 2000, 2000, 50, 500]),
        )
        print(result)
        print(resultcov)
        # coef = np.polyfit(x,ionization_energies[9:],3)
        # funct = np.polyval(coef,Z)
        # print(coef)
        import matplotlib.pyplot as plt

        plt.scatter(x, y, marker="*")
        plt.plot(x2, func(x2, *result), "r-")
        plt.show()


def get_ionization_energy_2p1_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.00,
        0.0,
        0.0,
        0.000,
        21.7,
        30.65,
        49.78,
        72.95,
        99.82,
        136,
        163.6,
        202,
        250.6,
        297.3,
        349.7,
        403.6,
        460.2,
        519.8,
        583.8,
        649.9,
        719.9,
        793.2,
        870.0,
        952.3,
        1044.9,
        1143.2,
        1248.1,
        1359.1,
        1474.3,
        1596,
        1730.9,
        1864,
        2007,
        2156,
        2307,
        2465,
        2625,
        2793,
        2967,
        3146,
        3330,
        3524,
        3727,
        3938,
        4156,
        4380,
        4612,
        4852,
        5107,
        5359,
        5624,
        5891,
        6164,
        6440,
        6722,
        7013,
        7312,
        7617,
        7930,
        8252,
        8581,
        8918,
        9264,
        9617,
        9978,
        10349,
        10739,
        11136,
        11544,
        11959,
        12385,
        12824,
        13273,
        13734,
        14209,
        14698,
        15200,
        15711,
        16244,
        16785,
        17337,
        17907,
        18484,
        19083,
        19693,
        20314,
        20948,
    ]
    # Coefficients from fit to data
    a1 = 1.30802855e01
    a2 = 1.51798476e-01
    a3 = 1.14195659e-02
    a4 = 1.58049679e00
    a5 = -2.74880264e01
    b1 = 4.99597550e01
    c = 1.48993379e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 5:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
        import np as np
        from scipy.optimize import curve_fit

        def func(x, n1, n2, n3, n4, n5, m1, o):
            return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

        x = range(10, 93)
        x2 = range(1, 103)
        y = ionization_energies[9:]
        result, resultcov = curve_fit(
            func, x, y, p0=[a1, a2, a3, a4, a5, b1, c], bounds=([0, 0, -np.inf, -np.inf, -np.inf, -50, -500], [2000, 1000, 1000, 2000, 2000, 50, 500])
        )
        print(result)
        print(resultcov)
        # coef = np.polyfit(x,ionization_energies[9:],3)
        # funct = np.polyval(coef,Z)
        # print(coef)
        import matplotlib.pyplot as plt

        plt.scatter(x, y, marker="*")
        plt.plot(x2, func(x2, *result), "r-")
        plt.show()


def get_ionization_energy_2p3_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.00,
        0.0,
        0.0,
        0.000,
        21.6,
        30.81,
        49.50,
        72.55,
        99.42,
        135,
        163.6,
        200,
        248.4,
        294.6,
        346.2,
        398.7,
        453.8,
        512.1,
        574.1,
        638.7,
        706.8,
        778.1,
        852.7,
        932.7,
        1021.8,
        1116.4,
        1217.0,
        1323.6,
        1433.9,
        1550,
        1678.4,
        1804,
        1940,
        2080,
        2223,
        2371,
        2520,
        2677,
        2838,
        3004,
        3173,
        3351,
        3538,
        3730,
        3929,
        4132,
        4341,
        4557,
        4786,
        5012,
        5247,
        5483,
        5723,
        5964,
        6208,
        6459,
        6716,
        6977,
        7243,
        7514,
        7790,
        8071,
        8358,
        8648,
        8944,
        9244,
        9561,
        9881,
        10207,
        10535,
        10871,
        11215,
        11564,
        11919,
        12284,
        12658,
        13035,
        13419,
        13814,
        14214,
        14619,
        15031,
        15444,
        15871,
        16300,
        16733,
        17166,
    ]
    # Coefficients from fit to data
    a1 = 6.26079141e-01
    a2 = 2.48380219e-01
    a3 = 3.01836854e-03
    a4 = 2.15937528e00
    a5 = -4.18326000e01
    b1 = 6.02739274e01
    c = 2.60817688e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 5:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
        from scipy.optimize import curve_fit

        def func(x, n1, n2, n3, n4, n5, m1, o):
            return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

        x = range(10, 93)
        x2 = range(7, 103)
        y = ionization_energies[9:]
        result, resultcov = curve_fit(
            func,
            x,
            y,
            p0=[a1, a2, a3, a4, a5, b1, c],
            bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
            maxfev=100000000,
        )
        print(result)
        print(resultcov)
        # coef = np.polyfit(x,ionization_energies[9:],3)
        # funct = np.polyval(coef,Z)
        # print(coef)
        import matplotlib.pyplot as plt

        plt.scatter(x, y, marker="*")
        plt.plot(x2, func(x2, *result), "r-")
        plt.show()


def get_ionization_energy_3s(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.00,
        0.0,
        0.0,
        0.000,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        29.3,
        34.8,
        44.3,
        51.1,
        58.7,
        66.3,
        74.1,
        82.3,
        91.3,
        101.0,
        110.8,
        122.5,
        139.8,
        159.5,
        180.1,
        204.7,
        229.6,
        257,
        292.8,
        326.7,
        358.7,
        392.0,
        430.3,
        466.6,
        506.3,
        544,
        586.1,
        628.1,
        671.6,
        719.0,
        772.0,
        827.2,
        884.7,
        946,
        1006,
        1072,
        1148.7,
        1211,
        1293,
        1362,
        1436,
        1511,
        1575,
        0.0,
        1723,
        1800,
        1881,
        1968,
        2047,
        2128,
        2207,
        2307,
        2398,
        2491,
        2601,
        2708,
        2820,
        2932,
        3049,
        3174,
        3296,
        3425,
        3562,
        3704,
        3851,
        3999,
        4149,
        4317,
        4482,
        4652,
        4822,
        5002,
        5182,
        5367,
        5548,
    ]
    # Coefficients from fit to data
    a1 = 1.43433111e02
    a2 = 5.93352595e-02
    a3 = 3.63930366e-03
    a4 = 1.96661938e-01
    a5 = -8.55996542e00
    b1 = 2.98052563e01
    c = 6.63407582e-05
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 11:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_3p_1_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.00,
        0.0,
        0.0,
        0.000,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        15.9,
        18.3,
        25.4,
        28.3,
        32.6,
        37.2,
        42.2,
        47.2,
        52.7,
        58.9,
        68.0,
        77.3,
        91.4,
        103.5,
        124.9,
        146.9,
        166.5,
        189,
        222.2,
        248.7,
        280.3,
        310.6,
        343.5,
        376.1,
        411.6,
        447.6,
        483.5,
        521.3,
        559.9,
        603.8,
        652.6,
        703.2,
        756.5,
        812.7,
        870.8,
        931,
        1002.1,
        1071,
        1137,
        1209,
        1274,
        1337,
        1403,
        1471,
        1541,
        1614,
        1688,
        1768,
        1842,
        1923,
        2006,
        2090,
        2173,
        2264,
        2365,
        2469,
        2575,
        2682,
        2792,
        2909,
        3027,
        3148,
        3279,
        3416,
        3554,
        3696,
        3854,
        4008,
        4159,
        4327,
        4490,
        4656,
        4830,
        5001,
        5182,
    ]
    # Coefficients from fit to data
    a1 = 1.70056228e02
    a2 = 5.84139911e-02
    a3 = -1.06353074e-02
    a4 = 1.50935984e00
    a5 = -6.65267026e01
    b1 = 2.65575820e00
    c = 4.99863585e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 13:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_3p_3_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.00,
        0.0,
        0.0,
        0.000,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        15.7,
        18.3,
        25.4,
        28.3,
        32.6,
        37.2,
        42.2,
        47.2,
        52.7,
        59.9,
        66.2,
        75.1,
        88.6,
        100.0,
        120.8,
        141.2,
        160.7,
        182,
        214.4,
        239.1,
        270.0,
        298.8,
        329.8,
        360.6,
        394.0,
        417.7,
        461.4,
        496.5,
        532.3,
        573.0,
        618.4,
        665.3,
        714.6,
        766.4,
        820.0,
        875,
        940.6,
        1003,
        1063,
        1128,
        1187,
        1242,
        1297,
        1357,
        1420,
        1481,
        1544,
        1611,
        1676,
        1741,
        1812,
        1885,
        1950,
        2024,
        2108,
        2194,
        2281,
        2367,
        2457,
        2551,
        2645,
        2743,
        2847,
        2957,
        3066,
        3177,
        3302,
        3426,
        3538,
        3663,
        3792,
        3909,
        4046,
        4174,
        4303,
    ]
    # Coefficients from fit to data
    a1 = 2.09817052e02
    a2 = 5.25431495e-02
    a3 = -1.34329853e-02
    a4 = 1.70347853e00
    a5 = -7.71441584e01
    b1 = -4.18496392e00
    c = 4.99999983e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 18:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_3d_3_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        10.2,
        18.7,
        29.8,
        41.7,
        55.5,
        70,
        95.0,
        113.0,
        136.0,
        157.7,
        181.1,
        205.0,
        231.1,
        257.6,
        284.2,
        311.9,
        340.5,
        374.0,
        411.9,
        451.4,
        493.2,
        537.5,
        583.4,
        630.8,
        689.0,
        740.5,
        795.7,
        853,
        902.4,
        948.3,
        1003.3,
        1052,
        1110.9,
        1158.6,
        1221.9,
        1276.9,
        1333,
        1392,
        1453,
        1515,
        1576,
        1639,
        1716,
        1793,
        1872,
        1949,
        2031,
        2116,
        2202,
        2291,
        2385,
        2485,
        2586,
        2688,
        2798,
        2909,
        3022,
        3136,
        3248,
        3370,
        3491,
        3611,
        3728,
    ]
    # Coefficients from fit to data
    a1 = 2.19640503e02
    a2 = 5.02997104e-02
    a3 = -1.30856531e-02
    a4 = 1.59217604e00
    a5 = -8.09352821e01
    b1 = -8.92936626e00
    c = 4.99999918e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 21:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_3d_5_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        10.1,
        18.7,
        29.2,
        41.7,
        54.6,
        69,
        93.8,
        112,
        134.2,
        155.8,
        178.8,
        202.3,
        227.9,
        253.9,
        280.0,
        307.2,
        335.2,
        368.3,
        405.2,
        443.9,
        484.9,
        528.2,
        573.0,
        619.3,
        676.4,
        726.6,
        780.5,
        836,
        883.8,
        928.8,
        980.4,
        1027,
        1083.4,
        1127.5,
        1189.6,
        1241.1,
        1292.6,
        1351,
        1409,
        1468,
        1528,
        1589,
        1662,
        1735,
        1809,
        1883,
        1960,
        2040,
        2122,
        2206,
        2295,
        2389,
        2484,
        2580,
        2683,
        2787,
        2892,
        3000,
        3105,
        3219,
        3332,
        3442,
        3552,
    ]
    # Coefficients from fit to data
    a1 = 2.19703182e02
    a2 = 4.99554795e-02
    a3 = -1.31912578e-02
    a4 = 1.60333067e00
    a5 = -8.08169197e01
    b1 = -8.82796315e00
    c = 4.99877841e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 21:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_4s(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        27.5,  # Kr
        30.5,
        38.9,
        43.8,
        50.6,
        56.4,
        63.2,
        69.5,
        75.0,
        81.4,
        87.1,
        97.0,  # Ag
        109.8,
        122.9,
        137.1,
        153.2,
        169.3,
        186.0,
        213.2,
        232.3,
        253.5,
        274.7,
        291.0,
        304.5,
        319.2,
        0.0,  # Pm
        347.2,
        360.0,
        378.6,
        396.0,
        414.2,
        432.4,
        449.8,
        470.9,
        480.5,  # Yb
        506.8,
        538.0,
        563.4,
        594.1,
        625.4,
        658.2,
        691.1,
        725.4,
        762.1,
        802.2,
        846.2,
        891.8,
        939.0,
        995.0,
        1042.0,
        1097.0,
        1153.0,
        1208.0,
        1269.0,
        1330.0,
        1387.0,
        1439.0,
    ]
    a1 = 2.76936897e02
    a2 = 4.81427814e-02
    a3 = -7.62990406e-03
    a4 = 5.40154613e-01
    a5 = -3.47989808e01
    b1 = -3.32677945e-01
    c = 1.11829908e-02
    # Coefficients from fit to data

    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 21:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
                if result < 0:
                    result = 0.0
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)


def get_ionization_energy_4p_1_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        14.1,  # Kr
        16.3,
        21.3,
        24.4,
        28.5,
        32.6,
        37.6,
        42.3,
        46.3,
        50.5,
        55.7,
        63.7,  # Ag
        63.9,
        73.5,
        83.6,
        95.6,
        103.3,
        123.0,
        146.7,
        172.4,
        192.0,
        205.8,
        223.2,
        236.3,
        243.3,
        242.0,  # Pm
        265.6,
        284.0,
        286.0,
        322.4,
        333.5,
        343.5,
        366.2,
        385.9,
        388.7,  # Yb
        412.4,
        438.2,
        463.4,
        490.4,
        518.7,
        549.1,
        577.8,
        609.1,
        642.7,
        680.2,
        720.5,
        761.9,
        805.2,
        851.0,
        886.0,
        929.0,
        980.0,
        1058.0,
        1080.0,
        1168.0,
        1224.0,
        1271.0,
    ]
    # Coefficients from fit to data
    a1 = 8.06479135e01
    a2 = 1.22162542e-01
    a3 = -1.07034361e-03
    a4 = 3.15458384e-01
    a5 = -1.53071387e01
    b1 = 6.70610396e01
    c = 1.94681955e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 31:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_4p_3_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        14.1,  # Kr
        15.3,
        20.1,
        23.1,
        27.1,
        30.8,
        35.5,
        39.9,
        43.2,
        47.3,
        50.9,
        58.3,  # Ag
        63.9,
        73.5,
        83.6,
        95.6,
        103.3,
        123.0,
        145.5,
        172.4,
        178.6,
        196.0,
        206.5,
        217.6,
        224.6,
        242.0,  # Pm
        247.4,
        257.0,
        271.0,
        284.1,
        293.2,
        308.2,
        320.2,
        332.6,
        339.7,  # Yb
        359.2,
        380.7,
        400.9,
        423.6,
        446.8,
        470.7,
        495.8,
        519.4,
        546.3,
        576.6,
        609.5,
        643.5,
        678.8,
        705.0,
        740.0,
        768.0,
        810.0,
        879.0,
        890.0,
        966.4,
        1007.0,
        1043.0,
    ]
    # Coefficients from fit to data
    a1 = 6.73916652e00
    a2 = 1.16287758e-01
    a3 = -2.33797834e-03
    a4 = 4.43185875e-01
    a5 = -1.95460837e01
    b1 = 3.45882667e01
    c = 2.38936562e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 31:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_5s(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,  # Kr
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,  # Ag
        0.0,
        0.0,
        0.0,  # Sn
        0.0,
        0.0,
        0.0,  # I
        23.3,  # Xe
        22.7,
        30.3,
        34.3,
        37.8,
        37.4,
        37.5,
        0.0,  # Pm
        37.4,
        32.0,
        36.0,
        45.6,
        49.9,
        49.3,
        50.6,
        54.7,
        52.0,  # Yb
        57.3,
        64.2,
        69.7,
        75.6,
        83.0,
        84.0,
        95.2,
        101.7,
        107.2,
        127,
        136.0,
        147.0,
        169.3,
        177.0,
        195.0,
        214.0,
        234.0,
        254.0,
        272.0,
        290.0,
        310.0,
        321.0,
    ]
    # Coefficients from fit to data
    a1 = 2.19703182e02
    a2 = 4.99554795e-02
    a3 = -1.31912578e-02
    a4 = 1.60333067e00
    a5 = -8.08169197e01
    b1 = -8.82796315e00
    c = 4.99877841e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 37:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_5p_1_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,  # Kr
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,  # Ag
        0.0,
        0.0,
        0.0,  # Sn
        0.0,
        0.0,
        0.0,  # I
        13.4,  # Xe
        14.2,
        17.0,
        19.3,
        19.8,
        22.3,
        21.1,
        0.0,  # Pm
        21.3,
        22.0,
        28.0,
        28.7,
        26.3,
        30.8,
        31.4,
        31.8,
        30.3,  # Yb
        33.6,
        38.0,
        42.2,
        45.3,
        45.6,
        58.0,
        63.0,
        65.3,
        74.2,
        83.1,
        94.6,
        106.4,
        119.0,
        132.0,
        148.0,
        164.0,
        182.0,
        200.0,
        215.0,
        229.0,
        232.0,
        257.0,
    ]
    # Coefficients from fit to data
    a1 = 2.19703182e02
    a2 = 4.99554795e-02
    a3 = -1.31912578e-02
    a4 = 1.60333067e00
    a5 = -8.08169197e01
    b1 = -8.82796315e00
    c = 4.99877841e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 21:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_ionization_energy_5p_3_2(Z: int) -> float:
    # table from X-ray data booklet, in eV
    ionization_energies = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,  # Kr
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,  # Ag
        0.0,
        0.0,
        0.0,  # Sn
        0.0,
        0.0,
        0.0,  # I
        12.1,  # Xe
        12.1,
        14.8,
        16.8,
        17.0,
        22.3,
        21.1,
        0.0,  # Pm
        21.3,
        22.0,
        28.0,
        22.6,
        26.3,
        24.1,
        24.7,
        25.0,
        24.1,  # Yb
        26.7,
        29.9,
        32.7,
        36.7,
        34.6,
        44.5,
        48.0,
        51.7,
        57.2,
        64.5,
        73.5,
        83.3,
        92.6,
        104.0,
        115.0,
        127.0,
        140.0,
        153.0,
        167.0,
        182.0,
        232.0,
        192.0,
    ]
    # Coefficients from fit to data
    a1 = 2.19703182e02
    a2 = 4.99554795e-02
    a3 = -1.31912578e-02
    a4 = 1.60333067e00
    a5 = -8.08169197e01
    b1 = -8.82796315e00
    c = 4.99877841e02
    if Z <= 92:
        result = ionization_energies[Z - 1]
        if result == 0.0:
            if Z >= 21:
                import numpy as np

                result = float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
        return result
    else:
        # extrapoalte from fourth order regression of tabulated data
        import numpy as np

        return float(a1 * np.exp2(a2 * (Z - b1)) + a3 * Z**3 + a4 * Z**2 + a5 * Z + c)
    from scipy.optimize import curve_fit

    def func(x, n1, n2, n3, n4, n5, m1, o):
        return n1 * np.exp2(n2 * (x - m1)) + n3 * x * x * x + n4 * x * x + n5 * x + o

    x = range(19, 93)
    x2 = range(7, 103)
    y = ionization_energies[18:]
    result, resultcov = curve_fit(
        func,
        x,
        y,
        p0=[a1, a2, a3, a4, a5, b1, c],
        bounds=([0, 0, -100, -100, -1000, -100, 0], [2000, 100, 100, 1000, 1000, 100, 500]),
        maxfev=100000000,
    )
    print(result)
    print(resultcov)
    # coef = np.polyfit(x,ionization_energies[9:],3)
    # funct = np.polyval(coef,Z)
    # print(coef)
    import matplotlib.pyplot as plt

    plt.scatter(x, y, marker="*")
    plt.plot(x2, func(x2, *result), "r-")
    plt.show()


def get_line_width_K(Z: int) -> float:
    # returns line width in eV
    # according to Kraus & Oliver J. Phys. Chem. Ref. Data, Vol.8, No.2, 1979
    Linewidths = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.24,
        0.30,
        0.36,
        0.42,
        0.48,
        0.53,
        0.59,
        0.64,
        0.68,
        0.74,
        0.81,
        0.86,
        0.94,
        1.01,
        1.08,
        1.16,
        1.25,
        1.33,
        1.44,
        1.55,
        1.67,
        1.82,
        1.96,
        2.14,
        2.33,
        2.52,
        2.75,
        2.99,
        3.25,
        3.52,
        3.84,
        4.14,
        4.52,
        4.91,
        5.33,
        5.77,
        6.24,
        6.75,
        7.28,
        7.91,
        8.49,
        9.16,
        9.89,
        10.6,
        11.4,
        12.3,
        13.2,
        14.1,
        15.1,
        16.2,
        17.3,
        18.5,
        19.7,
        21.0,
        22.3,
        23.8,
        25.2,
        26.8,
        28.4,
        30.1,
        31.9,
        33.7,
        35.7,
        37.7,
        39.9,
        42.1,
        44.4,
        46.8,
        49.3,
        52.0,
        54.6,
        57.4,
        60.4,
        63.4,
        66.6,
        69.8,
        73.3,
        76.8,
        80.4,
        84.1,
        88.0,
        91.9,
        96.1,
        100,
        105,
        109,
        114,
        119,
        124,
        129,
        135,
        140,
        145,
        150,
    ]
    return Linewidths[Z - 1]


def get_line_width_Lp_3_2(Z: int) -> float:
    # returns line width in eV
    # according to Kraus & Oliver J. Phys. Chem. Ref. Data, Vol.8, No.2, 1979
    Linewidths = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.1,
        0.20,
        0.41,
        0.73,
        1.03,
        1.26,
        1.49,
        1.58,
        1.63,
        1.92,
        2.07,
        2.21,
        2.34,
        2.41,
        2.54,
        2.62,
        2.76,
        2.79,
        2.89,
        3.06,
        3.28,
        3.38,
        3.53,
        3.79,
        3.94,
        4.11,
        4.28,
        4.44,
        4.67,
        4.71,
        4.78,
        3.94,
        4.25,
        4.36,
        4.58,
        4.73,
        4.93,
        4.88,
        4.87,
        5.00,
        2.97,
        3.13,
        3.32,
        3.46,
        3.64,
        3.78,
        3.92,
        4.06,
        4.21,
        4.34,
        4.52,
        4.67,
        4.80,
        4.91,
        5.05,
        5.19,
        5.25,
        5.33,
        5.43,
        5.47,
        5.53,
        5.54,
        5.63,
        5.58,
        5.61,
        6.18,
        7.25,
        8.30,
        9.39,
        10.5,
        11.3,
        12.0,
        12.2,
        12.4,
        12.6,
        12.8,
        13.1,
        13.3,
        13.4,
        13.6,
        13.7,
        14.3,
        14.0,
        14,
        13.5,
        13.3,
        13.6,
        13.8,
        14.0,
        14.3,
        14.4,
        14.8,
        15.1,
        15.9,
    ]
    return Linewidths[Z - 1]


def get_line_width_Lp_1_2(Z: int) -> float:
    # returns line width in eV
    # according to Kraus & Oliver J. Phys. Chem. Ref. Data, Vol.8, No.2, 1979
    Linewidths = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.001,
        0.004,
        0.015,
        0.032,
        0.054,
        0.083,
        0.126,
        0.152,
        0.17,
        0.19,
        0.24,
        0.26,
        0.29,
        0.34,
        0.37,
        0.43,
        0.52,
        0.62,
        0.72,
        0.83,
        0.95,
        1.03,
        1.13,
        1.21,
        1.31,
        1.43,
        1.54,
        1.65,
        1.78,
        1.87,
        1.97,
        2.08,
        2.23,
        2.35,
        2.43,
        2.57,
        2.62,
        2.72,
        2.84,
        3.00,
        3.12,
        3.25,
        3.40,
        3.51,
        3.57,
        3.68,
        3.80,
        3.89,
        3.97,
        4.06,
        4.15,
        4.23,
        4.32,
        4.43,
        4.55,
        4.66,
        4.73,
        4.79,
        4.82,
        4.92,
        5.02,
        5.15,
        5.33,
        5.48,
        5.59,
        5.69,
        5.86,
        6.00,
        6.17,
        6.32,
        6.48,
        6.67,
        6.83,
        7.01,
        7.20,
        7.47,
        7.68,
        7.95,
        8.18,
        8.75,
        9.32,
        9.91,
        10.5,
        10.9,
        11.4,
        11.8,
        12.2,
        12.7,
        13.1,
        13.6,
        14.0,
        14.4,
    ]
    return Linewidths[Z - 1]


def get_line_width_Ls(Z: int) -> float:
    # returns line width in eV
    # according to Kraus & Oliver J. Phys. Chem. Ref. Data, Vol.8, No.2, 1979
    Linewidths = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.00,
        0.001,
        0.004,
        0.014,
        0.033,
        0.054,
        0.087,
        0.128,
        0.156,
        0.17,
        0.19,
        0.22,
        0.24,
        0.27,
        0.32,
        0.36,
        0.43,
        0.48,
        0.56,
        0.65,
        0.76,
        0.82,
        0.94,
        1.00,
        1.08,
        1.17,
        1.27,
        1.39,
        1.50,
        1.57,
        1.66,
        1.78,
        1.91,
        2.00,
        2.13,
        2.25,
        2.40,
        2.50,
        2.65,
        2.75,
        2.87,
        2.95,
        3.08,
        3.13,
        3.25,
        3.32,
        3.41,
        3.48,
        3.60,
        3.65,
        3.75,
        3.86,
        3.91,
        4.01,
        4.12,
        4.17,
        4.26,
        4.35,
        4.48,
        4.60,
        4.68,
        4.80,
        4.88,
        4.98,
        5.04,
        5.16,
        5.25,
        5.31,
        5.41,
        5.50,
        5.65,
        5.81,
        5.98,
        6.13,
        6.29,
        6.41,
        6.65,
        6.82,
        6.98,
        7.13,
        7.33,
        7.43,
        7.59,
        7.82,
        8.04,
        8.26,
        8.55,
        8.75,
        9.04,
        9.33,
        9.61,
        9.90,
        10.1,
    ]
    return Linewidths[Z - 1]