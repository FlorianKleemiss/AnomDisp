import os
from constants_and_atomic_properties import *
import mpmath
from matrix_coefficients_angular import S_ort, S_par, P_ort, P_par, D_ort, D_par
import time
import interpolators as interp_mod
from interpolators import (
    I_integral,
    _build_adaptive_interpolator,
    _build_interpolator,
    _pool_worker_init,
    enable_precomputed_adaptive_interpolators,
    get_parts_for_typ_m,
    init_precomputed_interpolators,
    resolve_interpolators_file,
)

digits = 10
mpmath.mp.dps = digits

def norm_fact(l:int, m:int, sqrt_zone, b:float, n0: int):
    if m == 0:
        part1 = 1.0
    else:
        part1 = 2.0
    if m > l:
      return 0.0
    part1 *= (2.0 * l + 1) * mpmath.factorial(l - m)/mpmath.factorial(l + m) * b**4 * 2 * n0
    zm1 = (sqrt_zone**2).real
    negative = zm1 < 0
    if l == 0:
        pass
    else:
        product = 1.0
        for s in range(1, l+1):
            product *= (s**2 + n0**2 / zm1)
        part1 *= product
    part1 *= zm1**l
    if negative:
        return float(part1)
    else:
        denom = 1.0 / (1.0 - mpmath.exp(-2 * mpmath.pi * n0 / sqrt_zone))
        return float(part1 * denom)

def Rnl(args):
    sqrt_zone, br, l, n0 = args
    brz = sqrt_zone * br
    return -mpmath.power(2 * br, l) / mpmath.factorial(2 * l + 1) * mpmath.exp(1j * brz) * mpmath.hyp1f1(-n0*1j/sqrt_zone + l + 1, 2 * l + 2, -2j * brz)

def polynomial(n0, m, typ, br, part):
    if n0 == 1:
        return br**2
    elif n0 == 2:
        if typ == "s":
            return (2.0 - br) * br**2
        elif typ == "p":
            if m == 0:
                if part == 1:
                    return -2.0 * br**2
                elif part == 2:
                    return br**3
            elif m == 1:
                return br**3
            elif m == 2:
                return 0.5 * br**3
    elif n0 == 3:
        if typ == "s":
            return (3.0 - ten_thirds * br + two_thirds * br**2) * br**2
        elif typ == "p":
            sqr23 = mpmath.sqrt(mpmath.mpf(2.0)/mpmath.mpf(3.0))
            if m == 0:
                if part == 1:
                    return -sqr23 * (2 * (2 - br) * br**2)
                elif part == 2:
                    return sqr23 * br**3 * (3 - br)
            elif m == 1:
                return sqr23 * br**3 * (3 - br)
            elif m == 2:
                return sqr23 * 0.5 * br**3 * (3 - br)
        elif typ == "d":
            sqr23 = mpmath.sqrt(mpmath.mpf(2.0)/mpmath.mpf(3.0))
            if m == 0:
                if part == 1:
                    return -sqr23 * 2.0*br**3
                elif part == 2:
                    return sqr23 * (br**4)
            elif m == 1:
                if part == 1:
                    return -sqr23 * br**3
                elif part == 2:
                    return sqr23 * 0.25*br**4
            elif m == 2:
                return sqr23 * 0.5*br**4
            elif m == 3:
                return sqr23 * 0.25*br**4
            elif m == 4:
                if part == 1:
                    return mpmath.sqrt(mpmath.mpf(0.5)) * (2.0-br)/3.0 * br**3
                elif part == 2:
                    return mpmath.sqrt(mpmath.mpf(0.5)) * br**4
    elif n0 == 4:
        if typ == "s":
            return (4.0 - 7.0*br + 3.0*br**2 - (1.0/3.0)*br**3) * br**2
        elif typ == "p":
            sqr5 = mpmath.sqrt(mpmath.mpf(5.0))
            if m == 0:
                if part == 1:
                    return -sqr5 * 2*(1-br+br**2/5.0)*br**2
                elif part == 2:
                    return sqr5 * (2-1.4*br+0.2*br**2)*br**3
            elif m == 1:
                return sqr5 * (2-1.4*br+0.2*br**2)*br**3 # ARE WE SURE ABOUT THIS ONE?
            elif m == 2:
                return sqr5 * (1.0 - 0.7*br + 0.1*br**2) * br**3
    elif n0 == 5:
        if typ == "s":
            return (5.0 - 12.0*br + 8.0*br**2 - (28.0/15.0)*br**3 + (2.0/15.0)*br**4) * br**2

    raise ValueError("Polynomial not implemented for n0={}, typ={}, m={}, part={}".format(n0, typ, m, part))

def integrator(b_val, sqrt_zone, q_val, l=1, m=1, upper_val=1E-9, threshold=1E-6, n0=3, typ="s"):
    if upper_val == 0.0:
        return mpmath.mpf(0.0), mpmath.mpf(0.0)

    parts = get_parts_for_typ_m(typ, m)

    b_val = mpmath.mpf(b_val)
    zone = mpmath.mpc(sqrt_zone)
    q_val = mpmath.mpf(q_val)
    res_val = mpmath.mpf(0.0)
    res_err = mpmath.mpf(0.0)

    for part in parts:
        def integrand(r):
            br = b_val * r
            Rnl_val = Rnl((zone, br, l, n0))
            inner_int = I_integral((r * q_val, l, m, typ, part))
            return Rnl_val * inner_int * polynomial(n0, m, typ, br, part) * mpmath.exp(-br)
        if upper_val == -1:
            nodes = [0, mpmath.inf]
        else:
            nodes = [0, 0.25 * upper_val, 0.5 * upper_val, 0.75 * upper_val, upper_val]
        val, err = mpmath.quadts(integrand, nodes, error=True)
    
        # Optional: retry with higher precision if the reported error is large
        if err > threshold * abs(val):
            base_dps = 2 * digits
            if upper_val != -1:
                upper_val *= 1.1
                nodes = [0, 0.25 * upper_val, 0.5 * upper_val, 0.75 * upper_val, upper_val]
            else:
                print("Retrying integral to infinity with higher precision")
                nodes = [0, 0.25E-10, 0.5E-10, 0.75E-10, upper_val] # Subdivide more finely
            with mpmath.workdps(base_dps):
                val, err = mpmath.quadts(integrand, nodes, error=True)

        res_val += val
        res_err += err
    
    res_val *= 1.0 / b_val
    res_err *= 1.0 / b_val
    return res_val, res_err

def interpolgrator(b_val, sqrt_zone, q_val, l=1, m=1, upper_val=-1, threshold=1E-6, n0=3, typ="s", precomputed_data=None, interpolator_mode="adaptive"):
    if upper_val == 0.0:
        return mpmath.mpf(0.0), mpmath.mpf(0.0)

    parts = get_parts_for_typ_m(typ, m)
    if interpolator_mode == "adaptive":
        available_cache = interp_mod.PRECOMPUTED_INTERPOLATORS_ADAPTIVE if precomputed_data is None else precomputed_data
        inty = {
            part: _build_adaptive_interpolator(n0, l, m, typ, part, precomputed_data=available_cache)
            for part in parts
        }
    else:
        available_cache = interp_mod.PRECOMPUTED_INTERPOLATORS if precomputed_data is None else precomputed_data
        inty = {
            part: _build_interpolator(n0, l, m, typ, part, precomputed_data=available_cache)
            for part in parts
        }

    b_val = mpmath.mpf(b_val)
    zone = mpmath.mpc(sqrt_zone)
    q_val = mpmath.mpf(q_val)
    res_val = mpmath.mpf(0.0)
    res_err = mpmath.mpf(0.0)

    for part in parts:
        def integrand(r):
            br = b_val * r
            Rnl_val = Rnl((zone, br, l, n0))
            inner_int = inty[part](r * q_val)
            return Rnl_val * inner_int * polynomial(n0, m, typ, br, part) * mpmath.exp(-br)
        if upper_val == -1:
            upper_val = 1.5E1/b_val # This is a heuristic upper limit for the integral to infinity, can be adjusted, now aiming for the exp(-br) to become smaller than ~1E-20

        nodes = [0, 0.25 * upper_val, 0.5 * upper_val, 0.75 * upper_val, upper_val]
        val, err = mpmath.quadts(integrand, nodes, error=True)

        # Optional: retry with higher precision if the reported error is large
        if err > threshold * abs(val):
            base_dps = 2 * digits
            if upper_val != -1:
                upper_val *= 1.1
                nodes = [0, 0.25 * upper_val, 0.5 * upper_val, 0.75 * upper_val, upper_val]
            else:
                #print("Retrying integral to infinity with higher precision")
                nodes = [0, 0.25E-10, 0.5E-10, 0.75E-10, upper_val] # Subdivide more finely
            with mpmath.workdps(base_dps):
                val, err = mpmath.quadts(integrand, nodes, error=True)

        res_val += val
        res_err += err

    res_val *= 1.0 / b_val
    res_err *= 1.0 / b_val
    return res_val, res_err

def Xlm(args):
    """Computes value of the integral for a given Xlm incoming matrix coefficient
    need as args: n0, l, m, typ, b, k, q
    """
    n0, l_, m_, typ, b_, sqrt_zone, q_ = args
    result = integrator(b_, sqrt_zone, q_, l=l_, n0=n0, typ=typ, m=m_)
    #print(f"Computed l={l_val} integral up to r={r:e}, value={result[0]}, error estimate={result[1]}")
    return complex(result[0])

def Xlm_interpol(args):
    """Computes value of the integral for a given Xlm incoming matrix coefficient
    need as args: n0, l, m, typ, b, k, q
    """
    n0, l_, m_, typ, b_, sqrt_zone, q_ = args    
    result = interpolgrator(b_, sqrt_zone, q_, l=l_, n0=n0, typ=typ, m=m_, interpolator_mode="adaptive")
    #print(f"Computed l={l_val} integral up to r={r:e}, value={result[0]}, error estimate={result[1]}")
    return complex(result[0])

def matrix_product(args):
    """Compute the matrix product for given parameters.
    needs as args: n0, l, typ, b, k, q, theta_0
    """
    n0, l_, typ, b_, sqrt_zone, q_, theta_0, interpol, orth = args
    if typ == "s":
        averaging_factor = 1.0
        range_m = range(1,2)
        if orth:
            angle_matrix = S_ort(l_, theta_0)
        else:
            angle_matrix = S_par(l_, theta_0)
        
    elif typ == "p":
        averaging_factor = 1.0/3.0
        range_m = range(0,3)
        if orth:
            angle_matrix = P_ort(l_, theta_0)
        else:
            angle_matrix = P_par(l_, theta_0)
        
    elif typ == "d":
        averaging_factor = 1.0/5.0
        range_m = range(0,5)
        if orth:
            angle_matrix = D_ort(l_, theta_0)
        else:
            angle_matrix = D_par(l_, theta_0)
        
    if all(v == 0.0 for v in angle_matrix.flatten()):
        return 0.0
    
    normy = np.array([norm_fact(l_, m_, sqrt_zone, b_, n0) for m_ in range_m])  # Apply the norm factor to Xlm_val_conj
    norm_by_m = {m_: norm for m_, norm in zip(range_m, normy)}
    if interpol:
        Xlm_val = np.array([
            interpolgrator(b_, sqrt_zone, q_, l=l_, n0=n0, typ=typ, m=m_, interpolator_mode="adaptive")[0]
            if norm_by_m[m_] != 0.0 else 0.0
            for m_ in range_m
        ])
    else:
        Xlm_val = np.array([integrator(b_, sqrt_zone, q_, l=l_, n0=n0, typ=typ, m=m_)[0] for m_ in range_m])
    Xlm_val_conj = np.conj(Xlm_val).T
    Xlm_val_conj *= normy
    res = 0
    if Xlm_val.size == 1:
        res = (Xlm_val_conj[0] * angle_matrix[0] * Xlm_val[0])[0]
    elif Xlm_val.size >= 2:
        res =  Xlm_val_conj @ angle_matrix @ Xlm_val
    else:
        raise ValueError("Xlm_val has unexpected size:", Xlm_val.size)
    
    # Convert mpmath.mpc to Python float - extract real part
    try:
        return float(res.real * averaging_factor * mpmath.pi / 2)
    except (TypeError, AttributeError):
        # Fallback: convert to complex first, then extract real
        return float(complex(res).real * averaging_factor * mpmath.pi / 2)

def compute_func_for_r(args):
    """Helper function for parallel computation"""
    r, b_val, sqrt_zone, q_val, l_val, typ, m, n0 = args
    parts = get_parts_for_typ_m(typ, m)

    b_val = mpmath.mpf(b_val)
    zone = mpmath.mpc(sqrt_zone)
    q_val = mpmath.mpf(q_val)
    res_val = mpmath.mpf(0.0)

    for part in parts:
        br = b_val * r
        Rnl_val = Rnl((zone, br, l_val, n0))
        inner_int = I_integral((r * q_val, l_val, m, typ, part))
        val = Rnl_val * inner_int * polynomial(n0, m, typ, br, part) * mpmath.exp(-br)

        res_val += val
    return res_val / b_val * mpmath.pi

def compute_for_z(args):
    """Helper for parallel computation over z values; computes sum of matrix_product over all l values."""
    z_val, n0_, l_range, typ_, b_val, q_val, theta_0, orth = args
    sqrt_zone_ = mpmath.sqrt(z_val - 1.0)
    mp_args = [(n0_, l_, typ_, b_val, sqrt_zone_, q_val, theta_0, True, orth) for l_ in l_range]
    values = np.array(list(map(matrix_product, mp_args)))
    return float(values.sum())

def compute_for_z_no_interpolator(args):
    """Helper for parallel computation over z values; computes sum of matrix_product over all l values."""
    z_val, n0_, l_range, typ_, b_val, q_val, theta_0, orth = args
    sqrt_zone_ = mpmath.sqrt(z_val - 1.0)
    mp_args = [(n0_, l_, typ_, b_val, sqrt_zone_, q_val, theta_0, False, orth) for l_ in l_range]
    values = np.array(list(map(matrix_product, mp_args)))
    return float(values.sum())


def pool_map_with_progress(pool, func, items, label="pool"):
    """Run pool.imap while showing progress, elapsed time and ETA."""
    total = len(items)
    if total == 0:
        return []

    start = time.perf_counter()
    results = []
    update_step = max(1, total // 200)

    for i, result in enumerate(pool.imap(func, items), start=1):
        results.append(result)
        if i == 1 or i == total or i % update_step == 0:
            elapsed = time.perf_counter() - start
            rate = i / elapsed if elapsed > 0 else 0.0
            eta = (total - i) / rate if rate > 0 else float("inf")
            percent = 100.0 * i / total
            eta_str = "{:.1f}s".format(eta) if eta != float("inf") else "inf"
            print(
                "\r{}: {}/{} ({:.1f}%) | elapsed {:.1f}s | ETA {}".format(
                    label, i, total, percent, elapsed, eta_str
                ),
                end="",
                flush=True,
            )
    print()
    return results

def integral(args):
    n0, l_max, typ, b_val, q_val, z_0, delt, theta_0, orth = args
    l_range = range(0, l_max + 1)
    def integrant_for_z(z_):
        g = compute_for_z((z_, n0, l_range, typ, b_val, q_val, theta_0, orth))
        pre = 2.0 * z_ / (z_**2 - z_0**2)
        return pre * g

    upper_lim = mpmath.inf
    #upper_lim = 3.0
    z_start = mpmath.mpf(1.0) - mpmath.mpf(delt)

    # Pole from 1 / (z^2 - z0^2) handled here only at z = +z0,
    # assuming integration starts at positive z values.
    pole = mpmath.mpf(z_0)
    poles_in_range = [pole] if (mpmath.isfinite(pole) and pole > z_start) else []

    if not poles_in_range:
        return mpmath.quad(integrant_for_z, [z_start, upper_lim], maxdegree=4)

    # Single-pole case: split into three pieces
    # 1) left of pole, 2) symmetric window around pole (PV-style), 3) right of pole.
    pole = poles_in_range[0]
    window = min(mpmath.mpf("1e-2") * abs(pole), pole - z_start)
    if window <= 0:
        return mpmath.quad(integrant_for_z, [z_start, upper_lim], maxdegree=4)

    left_end = pole - window
    right_start = pole + window

    left_val = mpmath.mpf("0.0")
    if left_end > z_start:
        left_val = mpmath.quad(integrant_for_z, [z_start, left_end], maxdegree=4)

    # Symmetric principal-value contribution around the pole:
    # PV∫[p-w,p+w] f(z)dz = ∫[0,w] (f(p+t) + f(p-t)) dt
    def around_integrand(t):
        return integrant_for_z(pole + t) + integrant_for_z(pole - t)

    around_val = mpmath.quad(around_integrand, [0, window], maxdegree=6)
    right_val = mpmath.quad(integrant_for_z, [right_start, upper_lim], maxdegree=4)

    return left_val + around_val + right_val

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    from multiprocessing import Pool
    import scipy.constants as sc   
    from matrix_coefficients_v2 import g_poly, convert_fractions_to_float, evaluate_at_z
    import json
    from fest_verdrahtet import get_Ionization_Energy

    interpolator_file = resolve_interpolators_file()
    init_precomputed_interpolators(interpolator_file)
    adaptive_interpolator_file = "interpolator_adaptive.json"
    if os.path.exists(adaptive_interpolator_file):
        adaptive_count = enable_precomputed_adaptive_interpolators(adaptive_interpolator_file)
        print("Loaded {} adaptive interpolators from {}".format(adaptive_count, adaptive_interpolator_file))
    else:
        print("Adaptive interpolator file not found: {}".format(adaptive_interpolator_file))
    #Te
    Z=52
    
    #Tb
    Z = 65
    
    #W
    Z = 74
    typ = "s"
    l_max = 7
    m = 1
    n0 = 3
    l0 = -1
    if typ == "s":
        l0 = 0
    elif typ == "p":
        l0 = 1
    elif typ == "d":
        l0 = 2
    else:
        print("unknown type, nope!")
        exit(-1)
    #Typical q range for wl = 3.0 to 0.15 : 2.0E10 to 4.2E11
    nu_ag = sc.c / 0.560900e-10
    q_ = q(nu_ag)
    b_ = 1E11
    orth = False
    if typ == "s":
        j_string = "1/2"
        j_val = 0.5
        b_ = b(n0, 0, Z, j_string)
    elif typ == "p":
        j_string = "1/2"
        j_val = 0.5
        b_ = b(n0, 1, Z, j_string)
    elif typ == "d":
        j_string = "3/2"
        j_val = 1.5
        b_ = b(n0, 2, Z, j_string)

    E_n0l = -(get_Zeff(Z, n0, l0, j_string)**2 / n0**2 * Ryd_ener)
    IE = get_Ionization_Energy(Z, n0, l0, j_val)
    En0lj = Energie(Z, n0, l0, j_val)
    delta_z = (abs(En0lj) - IE) / abs(E_n0l)

    #for n0 = 1 we use
    #z = 0.617823
    #z = 0.889204

    #for n0 = 2 we use
    #z = 2.82399
    #z = 0.630988

    #for n0 = 3 we use
    #z = 0.394502
    #z = 8.66827

    #for n0 = 4 we use
    #z = 27.1351
    #z = 0.207952

    #for n0 = 5 we use
    #z = 222.761
    #z = 0.120932
    
    #z = 2.95835
    #z = 0.617245
    
    #z = 9.55004
    z = 1.001
    #z = 0.376221
    
    #z = 32.4982
    
    #z=11.6647
    
    #z = 2.0
    
    sqrt_zone = mpmath.sqrt(z - 1.0)
    z_n0l = h * nu_ag / abs(E_n0l)
    print("Performing calculations for b_:", b_, "q:", q_, "sigma:", Z - b_*n0*a0, "n0:", n0, "l0:", l0, "delta_z:", delta_z)
    print("1 - delta_z:", 1 - delta_z, " zn0l: ", z_n0l, "ionization energy:", IE, "En0lj:", En0lj, "E_n0l:", E_n0l)
    
    num_processes = 24
    print(f"Using up to {num_processes} CPU cores for parallel computation")
    
    mode = "Intensity"
    ran = None
    if mode == "values":
        ax1, ax2 = plt.subplots(2,1)[1]
    elif mode == "matrix_coefficient":
        ran = range(0, l_max)
        ax1, ax2 = plt.subplots(2,1)[1]
    elif mode == "matrix_product":
        ran = range(0, l_max)

    if mode == "values":
        z_range = np.linspace(0.7, 1.1, 100)
        for l_val in range(0, l_max + 1):
            print("l_val:", l_val)

            # Parallel computation
            r = 0.9 * 1E-10
            args_list = [(r, b_, mpmath.sqrt(z-1), q_, l_val, typ, m, n0) for z in z_range]
            with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
                values = np.array(pool_map_with_progress(pool, compute_func_for_r, args_list, label="values-mode pool"))
            
            #def norm_fact(l:int, m:int, sqrt_zone, b:float, n0: int):
            
            norm_vals = np.array([norm_fact(l_val, m, mpmath.sqrt(z-1), b_, n0) for z in z_range])

            ax1.plot(z_range, [v.real for v in values], label=f"l={l_val}")
            ax2.plot(z_range, [v.real for v in norm_vals], label=f"l={l_val}")
        ax1.set_xlabel("z")
        ax1.set_ylabel("Integrand value real")
        ax2.set_xlabel("z")
        ax2.set_ylabel("Norm factor value")
    
    elif mode == "normalizazion_factor":
        ax1 = plt.subplots(1,1)[1]
        l_ = 4
        m = 0
        z_range = np.linspace(0.7, 1.1, 150)
        norm_vals = np.array([norm_fact(l_, m, mpmath.sqrt(z-1), b_, n0) for z in z_range])
        from matrix_coefficients_v2 import horner_poly_dict
        
        norm_list = None
        with open("coefficients_norm.txt", "r") as file:
            norm_list = file.readlines()
            
        #has keys of the form (n0, l, m) and values of the form [norm1, norm2, norm3, ...]
        norm_dict = {}
        for line in norm_list:
            vals = convert_fractions_to_float(line).split(',')
            if len(vals) > 3:
                key = (int(vals[0]), int(vals[1]), int(vals[2]))
                norms = [float(val.strip()) for val in vals[3:]]
                norm_dict[key] = norms

        normy = norm_dict[(n0, l_, m)]
        poly_dict = {}
        for i,c in enumerate(normy):
            poly_dict[str(len(normy) - i - 1)] = c
        
        old_norm_vals = [-2.0 * horner_poly_dict(poly_dict, z) * b_**4 for z in z_range]
        ax1.plot(z_range, [v.real for v in norm_vals], label=f"Normalization factor for given l", linewidth=0, marker="*", markersize=5)
        ax1.plot(z_range, old_norm_vals, label=f"Normalization factor (old polynomial) for given l", linewidth=0, marker="o", markersize=3)
        ax1.set_xlabel("z")
        ax1.set_ylabel("Normalization factor value")
        
    elif mode == "matrix_coefficient":
        l_ = 1
        m = 0
        p_max = 10
        z_range = np.linspace(0.7, 1.1, 150)
            
        norm_list = None
        with open("coefficients_norm.txt", "r") as file:
            norm_list = file.readlines()

        #has keys of the form (n0, l, m) and values of the form [norm1, norm2, norm3, ...]
        norm_dict = {}
        for line in norm_list:
            vals = convert_fractions_to_float(line).split(',')
            if len(vals) > 3:
                key = (int(vals[0]), int(vals[1]), int(vals[2]))
                norms = [float(val.strip()) for val in vals[3:]]
                norm_dict[key] = norms

        polynomial_values = json.load(open("polynomials.json", "r"))
        from matrix_coefficients_v2 import w_schlange, horner_poly_dict
        ws = w_schlange(n0, l_, p_max, q_, b_, typ, polynomial_values)
        values_hoenl = np.zeros_like(z_range, dtype=complex)
        for i,z in enumerate(z_range):
            G_value = 0
            if z <= 1.0:
                G_value = mpmath.atan(mpmath.sqrt(1-z)) / mpmath.sqrt(1-z) - 2.0*(z-1) / 3.0 * mpmath.hyp2f1(0.75, 1.0, 1.75, (z-1)**2)
            else:
                G_value = mpmath.atan(mpmath.sqrt(z-1)) / mpmath.sqrt(z-1)
            values_hoenl[i] = -horner_poly_dict(ws[0], z)/z**(n0+1+p_max) / b_**2 * mpmath.pi *mpmath.exp(-2*n0*G_value)
        
        args_list = [(n0, l_, m, typ, b_, mpmath.sqrt(z-1), q_) for z in z_range]
        with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
            #values = np.array(pool.map(Xlm, args_list))
            values2 = np.array(pool_map_with_progress(pool, Xlm_interpol, args_list, label="integral-mode pool"))

        ax1.plot(z_range, [v.real for v in values_hoenl], label=f"hoenl", linewidth=0, marker="*", markersize=5)
        ax2.plot(z_range, [v.imag for v in values_hoenl], label=f"hoenl", linewidth=0, marker="*", markersize=5)
        ax1.plot(z_range, [v.real for v in values2], label=f"new", linewidth = 0, marker='o', markersize=3)
        ax2.plot(z_range, [v.imag for v in values2], label=f"new", linewidth = 0, marker='o', markersize=3)
        
        ax1.set_xlabel("z")
        ax2.set_xlabel("z")
        ax1.set_ylabel("Xlm value real")
        ax2.set_ylabel("Xlm value imag")

    elif mode == "matrix_product":
        ax1, ax2 = plt.subplots(2,1)[1]
        z_range = np.linspace(1-delta_z, 20.0, 197)
        g = np.zeros_like(z_range)
        g_hoenl = np.zeros_like(z_range)
        z_args = [(z_val, n0, list(ran), typ, b_, q_, 0.0, orth) for z_val in z_range]
        with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
            g = np.array(pool_map_with_progress(pool, compute_for_z, z_args, label="matrix-product pool"))
        #with Pool(processes=num_processes) as pool:
        #    g_full = np.array(pool.map(compute_for_z_no_interpolator, z_args))

        ##This was comparison to g_l from Hoenl like approach
        norm_list = None
        with open("coefficients_norm.txt", "r") as file:
            norm_list = file.readlines()

        #has keys of the form (n0, l, m) and values of the form [norm1, norm2, norm3, ...]
        norm_dict = {}
        for line in norm_list:
            vals = convert_fractions_to_float(line).split(',')
            if len(vals) > 3:
                key = (int(vals[0]), int(vals[1]), int(vals[2]))
                norms = [float(val.strip()) for val in vals[3:]]
                norm_dict[key] = norms

        polynomial_values = json.load(open("polynomials.json", "r"))
        #z, n0, l, q_, b_, p_max, t0, par_ort: bool, typ:str = "s"
        p_max = 10
        for i, z_val in enumerate(z_range):
            sqrt_zone = mpmath.sqrt(z_val - 1.0)
            poly = [n0, g_poly(n0, q_, b_, l_max=l_max-1, p_max=p_max, t0=0.0, par_ort=False, typ=typ, polynomial_values=polynomial_values, norm_dict=norm_dict)]
            g_hoenl[i] = evaluate_at_z(poly, z_val, p_max) * mpmath.pi

        #print("g values for given z:", g)
        #print("Difference between g and g_hoenl for given z:", g-g_hoenl)
        #print("ratio between g and g_hoenl for given z:", g/g_hoenl)
        ax1.scatter(z_range, g, label=f"Matrix product without interpolator for given z", linestyle='dotted')
        ax1.scatter(z_range, g_hoenl, label=f"Hoenl-like approach for given z", linestyle='dashed')
        ax2.scatter(z_range, (g-g_hoenl)/g, label=f"Relative Difference between matrix product and hoenl-like approach for given z", linestyle='dotted')
        ax1.set_xlabel("z")
        ax1.set_ylabel("Matrix product value")
        #set y max to 1.5* highest value of g
        ax1.set_ylim(0, 1.5*max(g))
        ax2.set_xlabel("z")
        ax2.set_ylabel("Difference")
        ax2.set_ylim(-1, 1)

    elif mode == "single_fp_fdp":
        ax1 = plt.subplots(1,1)[1]
        l_range = range(0, l_max + 1)

        #g_val = compute_for_z((z_n0l, n0, list(l_range), typ, b_, q_, 0.0))
        #print("f'' = " + str(float(g_val/2/mpmath.pi)))
        
        #exit(0)
        
        
        #z_range = np.linspace(1-delta_z, max(3,2*z_n0l), 100)
        #integrant_values = np.zeros_like(z_range, dtype=complex)
        #z_args = [(z_val, n0, list(l_range), typ, b_, q_, 0.0, orth) for z_val in z_range]
        #with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
        #    g_values = pool_map_with_progress(pool, compute_for_z, z_args, label="single_fp_fdp pool")
        #for i, (z_val, g) in enumerate(zip(z_range, g_values)):
        #    z_mp = mpmath.mpf(z_val)
        #    pre = 2.0 * z_mp / (z_mp**2 - z_n0l**2)
        #    val = complex(pre * g)
        #    print(f"z: {z_val:.3f}, integrand value: {val}")
        #    integrant_values[i] = val
        
        #ax1.plot(z_range, [float(1/2/mpmath.pi)*v.real for v in g_values], label=f"f''(z) real part", linestyle='dashed', marker='o', markersize=3)
        
        theta_range = np.linspace(0, np.pi, 100)
        
        reference1 = np.sin(theta_range)
        reference2 = np.sin(2*theta_range)
        reference3 = np.ones_like(theta_range)
        reference4 = np.cos(theta_range)
        reference5 = np.cos(2*theta_range)
        fp_values = np.zeros_like(theta_range)
        fdp_values = np.zeros_like(theta_range)
        
        args = [(n0, l_max, typ, b_, q(nu_ag), h * nu_ag / abs(E_n0l), delta_z, theta_0, orth) for theta_0 in theta_range]
        
        with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
            results = [float(result) for result in pool_map_with_progress(pool, integral, args, label="fp pool")]
            
        args2 = [(z_n0l, n0, list(l_range), typ, b_, q_, theta_0, orth) for theta_0 in theta_range]
        with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
            g_values = pool_map_with_progress(pool, compute_for_z, args2, label="g for fp pool")
        
        for i,r in enumerate(results):
            fp_values[i] = r
            reference1[i] *= results[0]
            reference2[i] *= results[0]
            reference3[i] *= results[0]
            reference4[i] *= results[0]
            reference5[i] *= results[0]
            
        for i,g in enumerate(g_values):
            fdp_values[i] = g
        
        ax1.plot(theta_range, fp_values, label=f"Integrand for single fp", linestyle='dashed', marker='o', markersize=3)
        ax1.plot(theta_range, fdp_values, label=f"Integrand for single fdp", linestyle='dashed', marker='s', markersize=3)
        ax1.plot(theta_range, reference1, label=f"Reference sin(theta)", linestyle='dotted')
        ax1.plot(theta_range, reference2, label=f"Reference sin(2*theta)", linestyle='dotted')
        ax1.plot(theta_range, reference3, label=f"Reference constant", linestyle='dotted')
        ax1.plot(theta_range, reference4, label=f"Reference cos(theta)", linestyle='dotted')
        ax1.plot(theta_range, reference5, label=f"Reference cos(2*theta)", linestyle='dotted')
        ax1.set_xlabel("θ (radians)")
        
        #g_val = compute_for_z((z_n0l, n0, list(l_range), typ, b_, q_, 0.0))
        #args = (n0, l_max, typ, b_, q(nu_ag), h * nu_ag / abs(E_n0l), delta_z, 0.0)
        #result = integral(args)
        #print(f"Computed anom disp correction value: {result} +  i*{g_val}")
        #exit(0)
        #ax1.plot(z_range, [v.real for v in integrant_values], label=f"Integrand for single fp", linestyle='dashed', marker='o', markersize=3)

    elif mode == "Intensity":
        #It has four contributions: (f_0 + f'_par)^2, (f_0 + f'_orth)^2, (f''_par)^2, (f''_orth)^2
        
        #W
        Z = 74
        #        a1     b1   a2    b2    a3     b3    a4   b4     c      f'    f"    mu         r   masse
        #SFAC  W 29.082 1.72 15.43 9.226 14.433 0.322 5.12 57.056 9.887 -0.668 6.765 28634.925 1.37 183.842
        #f_0 = c + sum_i a_i * exp(-b_i * cos(theta_0/2) / lambda)
        lambda_ag = 0.560900e-10
        
        def f_0_calculator(theta_0):
            f_0 = 9.887
            f_0 += 29.082 * np.exp(-1.72e-10 * np.cos(theta_0/2) / lambda_ag)
            f_0 += 15.43 * np.exp(-9.226e-10 * np.cos(theta_0/2) / lambda_ag)
            f_0 += 14.433 * np.exp(-0.322e-10 * np.cos(theta_0/2) / lambda_ag)
            f_0 += 5.12 * np.exp(-57.056e-10 * np.cos(theta_0/2) / lambda_ag)
            return f_0
        
        l_max = 7
        nu_ag = sc.c / lambda_ag
        q_ = q(nu_ag)
        
        theta_range = np.linspace(0, np.pi, 100)
        f_p_tot_par = np.zeros_like(theta_range)
        f_dp_tot_par = np.zeros_like(theta_range)
        f_p_tot_orth = np.zeros_like(theta_range)
        f_dp_tot_orth = np.zeros_like(theta_range)
        
        f_real_par = np.zeros_like(theta_range)
        f_real_orth = np.zeros_like(theta_range)
        f_imag_real = np.zeros_like(theta_range)
        f_imag_orth = np.zeros_like(theta_range)
        # n0, "typ", "l0", "j_string", j_val, electrons in this shell
        electrons = [
            (1, "s", 0, "1/2", 0.5, 2),
            (2, "s", 0, "1/2", 0.5, 2),
            (2, "p", 1, "1/2", 0.5, 2),
            (2, "p", 1, "3/2", 1.5, 4),
            (3, "s", 0, "1/2", 0.5, 2),
            (3, "p", 1, "1/2", 0.5, 2),
            (3, "p", 1, "3/2", 1.5, 4),
            (3, "d", 2, "3/2", 1.5, 4),
            (3, "d", 2, "5/2", 2.5, 6),
            (4, "s", 0, "1/2", 0.5, 2),
            (4, "p", 1, "1/2", 0.5, 2),
            (4, "p", 1, "3/2", 1.5, 4),
            (5, "s", 0, "1/2", 0.5, 2)
        ]
        
        for i,t in enumerate(theta_range):
            f_0 = f_0_calculator(t)
            f_real_par[i] += f_0
            f_real_orth[i] += f_0
        
        for e in electrons:
            n0_e, typ_e, l0_e, j_string_e, j_val_e, electrons_in_shell = e
            b_ = b(n0_e, l0_e, Z, j_string_e)
            E_n0l = -(get_Zeff(Z, n0_e, l0_e, j_string_e)**2 / n0_e**2 * Ryd_ener)
            IE = get_Ionization_Energy(Z, n0_e, l0_e, j_val_e)
            En0lj = Energie(Z, n0_e, l0_e, j_val_e)
            delta_z = (abs(En0lj) - IE) / abs(E_n0l)
            z_n0l = h * nu_ag / abs(E_n0l)
            
            orth = False
            args = [(n0, l_max, typ, b_, q(nu_ag), h * nu_ag / abs(E_n0l), delta_z, theta_0, orth) for theta_0 in theta_range]
        
            with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
                results = [float(result) for result in pool_map_with_progress(pool, integral, args, label=f"fp pool for shell n={n0_e}, typ={typ_e}, j={j_string_e}, orth={orth}")]

            args2 = [(z_n0l, n0, list(range(0, l_max + 1)), typ, b_, q_, theta_0, orth) for theta_0 in theta_range]
            with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
                g_values = pool_map_with_progress(pool, compute_for_z, args2, label=f"g for fp pool for shell n={n0_e}, typ={typ_e}, j={j_string_e}, orth={orth}")
            
            for i,r in enumerate(results):
                f_p_tot_par[i] += electrons_in_shell * r
                f_dp_tot_par[i] += electrons_in_shell * g_values[i]
                
            orth = True
            args = [(n0, l_max, typ, b_, q(nu_ag), h * nu_ag / abs(E_n0l), delta_z, theta_0, orth) for theta_0 in theta_range]
        
            with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
                results = [float(result) for result in pool_map_with_progress(pool, integral, args, label=f"fp pool for shell n={n0_e}, typ={typ_e}, j={j_string_e}, orth={orth}")]

            args2 = [(z_n0l, n0, list(range(0, l_max + 1)), typ, b_, q_, theta_0, orth) for theta_0 in theta_range]
            with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
                g_values = pool_map_with_progress(pool, compute_for_z, args2, label=f"g for fp pool for shell n={n0_e}, typ={typ_e}, j={j_string_e}, orth={orth}")
            
            for i,r in enumerate(results):
                f_p_tot_orth[i] += electrons_in_shell * r
                f_dp_tot_orth[i] += electrons_in_shell * g_values[i]
            
        for i in range(len(theta_range)):
            f_real_par[i] += f_p_tot_par[i]
            f_real_orth[i] += f_p_tot_orth[i]
            f_imag_real[i] += f_dp_tot_par[i]
            f_imag_orth[i] += f_dp_tot_orth[i]
            
        Intensity = np.zeros_like(f_real_par)
        Intensity_ohne_winkel = np.zeros_like(f_real_par)
        
        for i in range(len(theta_range)):
            Intensity[i] = 0.5 * (f_real_par[i]**2 + f_real_orth[i]**2 + f_imag_real[i]**2 + f_imag_orth[i]**2)
            Intensity_ohne_winkel[i] = 0.5 * Intensity[0] * (1+np.cos(theta_range[i])**2)

        ax1 = plt.subplots(1,1)[1]
        ax1.plot(theta_range, Intensity, label=f"Total intensity", linestyle='dashed', marker='o', markersize=3)
        ax1.plot(theta_range, Intensity_ohne_winkel, label=f"Intensity without angle", linestyle='dotted', marker='s', markersize=3)
        ax1.set_xlabel("θ (radians)")
        ax1.set_ylabel("Intensity (Depolarized)")
            
    elif mode == "fp_fdp":
        ax1 = plt.subplots(1,1)[1]
        #First we need to choose an incoming wavelength
        wls = np.linspace(0.15e-10, 2.5e-10, 30)  # Wavelengths from 0.15 Å to 2.5 Å
        nu_range = sc.c / wls  # Convert wavelengths to frequencies
        #nu_range = [nu_ag]  # For testing with a single frequency

        fp_values = np.zeros_like(nu_range)
        fdp_values = np.zeros_like(nu_range)

        args = [(n0, l_max, typ, b_, q(nu_val), h * nu_val / abs(E_n0l), delta_z, 0.0, orth) for nu_val in nu_range]

        with Pool(processes=num_processes, initializer=_pool_worker_init, initargs=(interpolator_file,)) as pool:
            results = [float(result) for result in pool_map_with_progress(pool, integral, args, label="fp/fdp pool")]
        #results = []
        #for arg in args:
        #    result = integral(arg)
        #    results.append(float(result))
        for i,r in enumerate(results):
            fp_values[i] = r

        print(wls, fp_values)
        ax1.plot(nu_range, fp_values, label=f"f'", linestyle='dashed')
        #ax1.plot(nu_range, fdp_values, label=f"f''", linestyle='dashed')
        ax1.set_xlabel("ν (Hz)")
        ax1.set_ylabel("Value")
        #exit(0)

    #min_y, max_y = ax1.get_ylim()
    plt.legend(ncol=3, fontsize='small')
    plt.tight_layout()
    plt.show()
