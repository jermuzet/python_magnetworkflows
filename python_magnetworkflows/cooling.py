"""
Cooling Models
"""

from math import exp, log, log10, sqrt

import freesteam


def rho(Tw: float, P: float) -> float:
    """
    compute water volumic mass in ???

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    # print(f'rho={Steam.rho}')
    return Steam.rho


def Cp(Tw: float, P: float) -> float:
    """
    compute water specific heat in ???

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    # print(f'Cp={Steam.cp}')
    return Steam.cp


def viscosity(Tw: float, P: float) -> float:
    """
    compute water viscosity in ??

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    return Steam.mu


def k(Tw: float, P: float) -> float:
    """
    compute water heat conductivity in ??

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    return Steam.k


def Montgomery(
    Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float
) -> float:
    """
    compute heat exchange coefficient in ??

    Tw: K
    Umean: m/s
    Dh: meter

    see: Montgomery p38 eq 3.3 Watch out for Unit change
    Montgomery formukq is given for Length==Centimeter, T in Celsius
    HMFL introduce an additional fuzzy factor
    """
    return 1426.404 * (1 + 1.5e-2 * (Tw - 273)) * exp(log(U) * 0.8) / exp(log(Dh) * 0.2)


def Dittus(Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float) -> float:
    params = (0.023, 0.8, 0.4)
    return hcorrelation(params, Tw, Pw, dPw, U, Dh, L)


def Colburn(Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float) -> float:
    params = (0.023, 0.8, 0.3)
    return hcorrelation(params, Tw, Pw, dPw, U, Dh, L)


def Silverberg(
    Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float
) -> float:
    params = (0.015, 0.85, 0.3)
    return hcorrelation(params, Tw, Pw, dPw, U, Dh, L)


def Blasius(Re: float, Dh: float, f: float, rugosity: float) -> float:
    return 0.316 / exp(log(Re) * 0.25)


def Filonenko(Re: float, Dh: float, f: float, rugosity: float) -> float:
    return 1 / (1.82 * log10(Re) - 1.64) ** 2


def Colebrook(self, Re: float, Dh: float, f: float, rugosity: float) -> float:
    val = 1 / sqrt(f)
    while True:
        nval = -2 * log10(rugosity / (3.7 * Dh) + 2.51 / Re * val)
        error = abs(1 - nval / val)
        val = nval

        if error <= 1.0e-3:
            break

    return 1 / val ** 2


def Swanee(self, Re: float, Dh: float, f: float, rugosity: float) -> float:
    val = f
    while True:
        nval = 1.3254 / log(rugosity / (3.75 * Dh) + 5.74 / exp(log(Re) * 0.9)) ** 2
        error = abs(1 - nval / val)

        val = nval

        if error <= 1.0e-3:
            break

    return val


def Uw(
    Tw: float,
    Pw: float,
    dPw: float,
    Dh: float,
    L: float,
    friction: str = "Colebrook",
    Pextra: float = 1.0,
    fguess: float = 0.055,
    uguess: float = 0,
    rugosity: float = 0,
) -> float:
    """
    compute water velocity
    """
    friction_method = {
        "Blasius": Blasius,
        "Filonenko": Filonenko,
        "Colebrook": Colebrook,
        "Swanee": Swanee,
    }

    U = uguess
    f = fguess
    while True:
        Re = Reynolds(Tw, Pw, U, Dh, L)
        nf = friction_method[friction](Re, Dh, f, rugosity)

        nU = sqrt(2 * dPw / (rho(Tw, Pw) * (Pextra + nf * L / Dh)))
        error_U = abs(1 - nU / U)
        error_f = abs(1 - nf / f)
        U = nU
        f = nf
        if error_U <= 1.0e-3 and error_f <= 1.0e-3:
            break
    return U


def Nusselt(params: tuple, Tw: float, Pw: float) -> float:
    (alpha, n, m) = params
    Re = Reynolds()
    Pr = Prandlt()
    Nu = alpha * exp(log(Re) * n) * exp(log(Pr) * m)
    return Re


def Reynolds(Tw: float, Pw: float, U: float, Dh: float, L: float) -> float:
    Re = rho(Tw, Pw) * U * Dh / viscosity(Tw, Pw)
    return Re


def Prandlt(Tw: float, Pw: float) -> float:
    Pr = viscosity(Tw, Pw) * Cp(Tw, Pw) / k(Tw, Pw)
    return Pr


def hcorrelation(
    params: tuple, Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float
) -> float:
    """
    compute heat exchange coeff in W/m²/K

    h = alpha * Re() * Pr()*m / Dh

    params (alpha, n, m): coeffs of correlations
    Tw: K
    Pw: bar
    dPw: bar # given by Pump
    Uw: meter/second
    Dh: meter
    L: meter
    """

    (alpha, n, m) = params
    nU = Uw(
        Tw, Pw, dPw, Dh, L, "Colebrook", Pextra=1, fguess=0.055, uguess=U, rugosity=0
    )

    Re = Reynolds(Tw, Pw, nU, Dh, L)
    Pr = Prandlt(Tw, Pw)

    h = alpha * exp(log(Re) * n) * exp(log(Pr) * m) / Dh
    return h
