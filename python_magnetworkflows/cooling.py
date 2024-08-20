"""
Cooling Models
"""

from typing import List
from math import exp, log, log10, sqrt

from iapws import IAPWS97


def steam(Tw: float, P: float):
    """
    return steam object

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Mpa = 0.1 * P
    # print(f"steam: {Mpa} MPa, {Tw} K", flush=True)
    # _steam = IAPWS97(P=Mpa, T=Tw)
    # print(f"_steam: rho={_steam.rho}, cp={_steam.cp* 1.0e3}", flush=True)
    return IAPWS97(P=Mpa, T=Tw)


# def rho(Tw: float, P: float) -> float:
#     """
#     compute water volumic mass in ???

#     Tw: Temperature in Kelvin
#     P: Pressure in Bar
#     """
#     Mpa = 0.1 * P
#     Steam = IAPWS97(P=Mpa, T=Tw)
#     # print(f"rho({Mpa} MPa, {Tw} K)={Steam.rho}")
#     return Steam.rho

# def Cp(Tw: float, P: float) -> float:
#     """
#     compute water specific heat in ???kJ/kg·K

#     Tw: Temperature in Kelvin
#     P: Pressure in Bar
#     """
#     Mpa = 0.1 * P
#     Steam = IAPWS97(P=Mpa, T=Tw)
#     # print(f"Cp({Mpa} MPa, {Tw} K)={Steam.cp} kJ/kg·K")
#     return Steam.cp * 1.0e3  # watch out units

# def viscosity(Tw: float, P: float) -> float:
#     """
#     compute water viscosity in ??

#     Tw: Temperature in Kelvin
#     P: Pressure in Bar
#     """
#     Mpa = 0.1 * P
#     Steam = IAPWS97(P=Mpa, T=Tw)
#     # print(f"mu({Mpa} MPa, {Tw} K)={Steam.mu}")
#     return Steam.mu

# def k(Tw: float, P: float) -> float:
#     """
#     compute water heat conductivity in ??

#     Tw: Temperature in Kelvin
#     P: Pressure in Bar
#     """
#     Mpa = 0.1 * P
#     Steam = IAPWS97(P=Mpa, T=Tw)
#     # print(f"k({Mpa} MPa, {Tw} K)={Steam.k}")
#     return Steam.k


def Montgomery(
    Tw: float,
    Pw: float,
    dPw: float,
    U: float,
    Dh: float,
    L: float,
    friction: str,
    fuzzy: float = 1.0,
) -> float:
    """
    compute heat exchange coefficient in ??

    Tw: K
    Umean: m/s
    Dh: meter

    see: Montgomery p38 eq 3.3 Watch out for Unit change
    Montgomery formula is given for Length==Centimeter, T in Celsius
    HMFL introduce an additional fuzzy factor
    """

    # fuzzy = 1.7
    h = (
        fuzzy
        * 1426.404
        * (1 + 1.5e-2 * (Tw - 273))
        * exp(log(U) * 0.8)
        / exp(log(Dh) * 0.2)
    )
    # print(f"hcorrelation(Montgomery): h={h}")
    return h


def Dittus(
    Tw: float,
    Pw: float,
    dPw: float,
    U: float,
    Dh: float,
    L: float,
    friction: str,
    fuzzy: str = 1.0,
) -> float:
    params = (0.023, 0.8, 0.4)
    h = hcorrelation(params, Tw, Pw, dPw, U, Dh, L, friction, "Dittus")
    return h


def Colburn(
    Tw: float,
    Pw: float,
    dPw: float,
    U: float,
    Dh: float,
    L: float,
    friction: str,
    fuzzy: str = 1.0,
) -> float:
    params = (0.023, 0.8, 0.3)
    h = hcorrelation(params, Tw, Pw, dPw, U, Dh, L, friction, "Colburn")
    return h


def Silverberg(
    Tw: float,
    Pw: float,
    dPw: float,
    U: float,
    Dh: float,
    L: float,
    friction: str,
    fuzzy: str = 1.0,
) -> float:
    params = (0.015, 0.85, 0.3)
    h = hcorrelation(params, Tw, Pw, dPw, U, Dh, L, friction, "Silverberg")
    return h


def Constant(Re: float, Dh: float, f: float, rugosity: float) -> float:
    cf = 0.055
    # print(f"Constant={cf}")
    return cf


"""
To be implemented

friction == "karman":
            iterate = True
            Cf  = math.pow(1.93*math.log10(Reynolds*math.sqrt(f))-0.537,-2)
        elif friction == "rough":
            iterate = True
            eps = 2.5e-2 # mm
            rstar = 1/math.sqrt(8.) * (Reynolds*math.sqrt(f))*eps/dh
            brstar = 1/(1.930*math.sqrt(f)) + math.log10(1.9/math.sqrt(8.) * eps/dh)
            ###print "brstar=%g" % brstar

            # Cf = math.pow(-1.930*math.log(1.90/(Reynolds*math.sqrt(f))*(1+0.34*rstar*math.exp(-11./rstar))),-2.)
            Cf = math.pow(-2.00*math.log10(2.51/(Reynolds*math.sqrt(f))*(1+rstar/3.3)),-2.)
            
        # Gnielinski breaks when a tends to 1
        # elif friction == "gnielinski":
        #     a = diameter_ratio
        #     Re = Reynolds * ( (1.+a**2) * math.log(a)+(1-a**2) / ( (1.-a)**2 * math.log(a) )) 
        #     Cf = math.pow(1.8*math.log10(Re)-1.5,-2)
        # # print ("%s Cf=%g" % (friction,Cf) )
"""


def Blasius(Re: float, Dh: float, f: float, rugosity: float) -> float:
    cf = 0.316 / exp(log(Re) * 0.25)
    # print(f"Blasius={cf}")
    return cf


def Filonenko(Re: float, Dh: float, f: float, rugosity: float) -> float:
    cf = 1 / (1.82 * log10(Re) - 1.64) ** 2
    # print(f"Filonenko={cf}")
    return cf


def Colebrook(Re: float, Dh: float, f: float, rugosity: float) -> float:
    val = 1 / sqrt(f)
    isOK = False
    it = 0
    max_err = 1.0e-3
    while it < 10:
        nval = -2 * log10(rugosity / (3.7 * Dh) + 2.51 / Re * val)
        error = abs(1 - nval / val)
        val = nval
        # print(f"Colebrook: val={val}, error={error}, it={it}")

        it += 1
        if error <= max_err:
            isOk = True
            break

    if isOk != True:
        raise RuntimeError(f"ColeBrook: cf failed to converge")
    cf = 1 / val**2
    # print(f"Colebrook={cf}, it={it}")
    return cf


def Swanee(Re: float, Dh: float, f: float, rugosity: float) -> float:
    val = f
    isOK = False
    it = 0
    max_err = 1.0e-3
    while it < 10:
        nval = 1.3254 / log(rugosity / (3.75 * Dh) + 5.74 / exp(log(Re) * 0.9)) ** 2
        error = abs(1 - nval / val)

        val = nval
        # print(f"Colebrook: val={val}, error={error}, it={it}")
        it += 1

        if error <= max_err:
            isOK = True
            break

    if isOK != True:
        raise RuntimeError(f"Swanee: cf failed to converge")
    cf = val
    # print(f"Swanee={cf}, it={it}")
    return cf


def Uw(
    Steam,
    dPw: float,
    Dh: float,
    L: float,
    friction: str = "Colebrook",
    Pextra: float = 1.0,
    fguess: float = 0.055,
    uguess: float = 0,
    rugosity: float = 0.012e-3,
) -> tuple:
    """
    compute water velocity
    """
    friction_method = {
        "Constant": Constant,
        "Blasius": Blasius,
        "Filonenko": Filonenko,
        "Colebrook": Colebrook,
        "Swanee": Swanee,
    }

    U = uguess
    f = fguess

    isOk = False
    it = 0
    error_U = 0
    error_f = 0
    max_err_U = 1.0e-3
    max_err_f = 1.0e-3
    while it < 10:
        Re = Reynolds(Steam, U, Dh, L)
        nf = friction_method[friction](Re, Dh, f, rugosity)

        dPw_Pascal = dPw * 1.0e5
        nU = sqrt(2 * dPw_Pascal / (Steam.rho * (Pextra + nf * L / Dh)))  # Faux!!!
        error_U = abs(1 - nU / U)
        error_f = abs(1 - nf / f)
        # print(
        #     f"Uw: U={U:.3f}, nU={nU:.3f}, f={f}, nf={nf}, Re(Tw={Tw:.3f}, Pw={Pw:.3f}, U={U:.3f}, Dh={Dh:.3e}, L={L:.3e})={Re:.3f}, dPw={dPw:.3f}, Pextra={Pextra:.3f}, error_U={error_U}, error_f={error_f}, it={it}"
        # )
        U = nU
        f = nf

        it += 1
        if error_U <= max_err_U and error_f <= max_err_f:
            isOk = True
            break

    # print(f"Uw={U:.3f}, Cf={f:.3e} ({friction}), rugosity={rugosity:.3e}, Re={Re:.3f}")
    if isOk != True:
        raise RuntimeError(f"Uw: max it reached")
    return U, f


def Nusselt(params: tuple, Re: float, Pr: float) -> float:
    """Compute Nusselt nb from Reynolds (Re) and Prandtl (Pr)"""
    (alpha, n, m) = params
    Nu = alpha * exp(log(Re) * n) * exp(log(Pr) * m)
    # print(f"Nu={Nu}, Pr={Pr}, Re={Re}")
    return Nu


def Reynolds(Steam, U: float, Dh: float, L: float) -> float:
    """Compute Reynolds as Re = rho*U*Dh/mu"""
    Re = Steam.rho * U * Dh / Steam.mu
    # print(f"Re={Re}")
    return Re


def Prandlt(
    Steam,
) -> float:
    """Compute Prandlt as Pr = mu*cp/k"""
    Pr = Steam.mu * Steam.cp * 1.0e3 / Steam.k
    # print(f"Pr={Pr}")
    return Pr


def hcorrelation(
    params: tuple,
    Tw: float,
    Pw: float,
    dPw: float,
    U: float,
    Dh: float,
    L: float,
    friction: str = "Constant",
    model: str = "Montgomery",
    rugosity: float = 0.012e-3,
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
    Steam = steam(Tw, Pw)
    nU = Uw(
        Steam,
        dPw,
        Dh,
        L,
        friction,
        Pextra=1,
        fguess=0.055,
        uguess=U,
        rugosity=rugosity,
    )

    Re = Reynolds(Steam, nU, Dh, L)
    Pr = Prandlt(Steam)

    h = alpha * exp(log(Re) * n) * exp(log(Pr) * m) / Dh
    # print(f"hcorrelation({model}): friction={friction}, h={h}, Pr={Pr}, Re={Re}")
    del Steam
    del nU
    del Re
    del Pr
    return h


def getDT(flow: float, Power: float, Tw: float, P: float) -> float:
    """
    compute dT as Power / rho *Cp * Flow(I)
    """
    Steam = steam(Tw, P)
    _DT = Power / (Steam.rho * Steam.cp * 1.0e3 * flow)
    del Steam
    return _DT


def getHeatCoeff(
    Dh: float,
    L: float,
    U: float,
    Tw: float,
    Pw: float,
    dPw: float,
    model: str = "Montgomery",
    friction: str = "Constant",
    fuzzy: float = 1.0,
):
    correlation = {
        "Montgomery": Montgomery,
        "Dittus": Dittus,
        "Colburn": Colburn,
        "Silverberg": Silverberg,
    }

    return correlation[model](Tw, Pw, dPw, U, Dh, L, friction, fuzzy)


def getTout(
    T: List[float], VolMass: List[float], SpecHeat: List[float], Q: List[float]
) -> float:
    Tout = 0
    rhoCpQ = 0
    # print(f"Sum(Qi)={sum(Q)}")
    for i, (Ti, RHOi, CPi, Qi) in enumerate(zip(T, VolMass, SpecHeat, Q)):
        # print(f"i:{i}, (Ti:{Ti}, RHOi:{RHOi}, CPi:{CPi}, Qi:{Qi})")
        Tout += Ti * RHOi * CPi * Qi
        rhoCpQ += RHOi * CPi * Qi

    Tout /= rhoCpQ
    return Tout
