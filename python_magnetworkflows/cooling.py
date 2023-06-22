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
    # print(f'rho({P} bar, {Tw} K)={Steam.rho}')
    return Steam.rho


def Cp(Tw: float, P: float) -> float:
    """
    compute water specific heat in ???

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    # print(f'Cp({P} bar, {Tw} K)={Steam.cp}')
    return Steam.cp


def viscosity(Tw: float, P: float) -> float:
    """
    compute water viscosity in ??

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    # print(f'mu({P} bar, {Tw} K)={Steam.mu}')
    return Steam.mu


def k(Tw: float, P: float) -> float:
    """
    compute water heat conductivity in ??

    Tw: Temperature in Kelvin
    P: Pressure in Bar
    """
    Steam = freesteam.steam_pT(P * 1.0e5, Tw)
    # print(f'k({P} bar, {Tw} K)={Steam.k}')
    return Steam.k


def Montgomery(
    Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float, friction: str
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
    h = 1426.404 * (1 + 1.5e-2 * (Tw - 273)) * exp(log(U) * 0.8) / exp(log(Dh) * 0.2)
    print(f"hcorrelation(Montgomery): h={h}")
    return h


def Dittus(
    Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float, friction: str
) -> float:
    params = (0.023, 0.8, 0.4)
    h = hcorrelation(params, Tw, Pw, dPw, U, Dh, L, friction, "Dittus")
    return h


def Colburn(
    Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float, friction: str
) -> float:
    params = (0.023, 0.8, 0.3)
    h = hcorrelation(params, Tw, Pw, dPw, U, Dh, L, friction, "Colburn")
    return h


def Silverberg(
    Tw: float, Pw: float, dPw: float, U: float, Dh: float, L: float, friction: str
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
    isOK= False
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
        raise RuntimeError(f'ColeBrook: cf failed to converge')    
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
        raise RuntimeError(f'Swanee: cf failed to converge')    
    cf = val
    # print(f"Swanee={cf}, it={it}")
    return cf


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
    rugosity: float = 0.012e-3,
) -> float:
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
        Re = Reynolds(Tw, Pw, U, Dh, L)
        nf = friction_method[friction](Re, Dh, f, rugosity)
        
        dPw_Pascal = dPw * 1.e+5
        nU = sqrt(2 * dPw_Pascal / (rho(Tw, Pw) * (Pextra + nf * L / Dh))) # Faux!!!
        error_U = abs(1 - nU / U)
        error_f = abs(1 - nf / f)
        U = nU
        f = nf
        # print(f'Uw: nU={nU}, nf={nf}, Re={Re}=Re(Tw={Tw}, Pw={Pw}, U={U}, Dh={Dh}, L={L}), dPw={dPw}, Pextra={Pextra}, error_U={error_U}, error_f={error_f}, it={it}')

        it += 1
        if error_U <= max_err_U and error_f <= max_err_f:
            isOk = True
            break

    print(f"Uw={U}, Cf={f} ({friction}), rugodity={rugosity}, Re={Re}")
    if isOk != True:
        raise RuntimeError(f'Uw: max it reached')
    return U


def Nusselt(params: tuple, Tw: float, Pw: float) -> float:
    (alpha, n, m) = params
    Re = Reynolds()
    Pr = Prandlt()
    Nu = alpha * exp(log(Re) * n) * exp(log(Pr) * m)
    # print(f"Nu={Nu}, Pr={Pr}, Re={Re}")
    return Re


def Reynolds(Tw: float, Pw: float, U: float, Dh: float, L: float) -> float:
    Re = rho(Tw, Pw) * U * Dh / viscosity(Tw, Pw)
    # print(f"Re={Re}")
    return Re


def Prandlt(Tw: float, Pw: float) -> float:
    Pr = viscosity(Tw, Pw) * Cp(Tw, Pw) / k(Tw, Pw)
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
    rugosity: float = 0.012e-3
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
    nU = Uw(Tw, Pw, dPw, Dh, L, friction, Pextra=1, fguess=0.055, uguess=U, rugosity=rugosity)

    Re = Reynolds(Tw, Pw, nU, Dh, L)
    Pr = Prandlt(Tw, Pw)

    h = alpha * exp(log(Re) * n) * exp(log(Pr) * m) / Dh
    print(f"hcorrelation({model}): friction={friction}, h={h}, Pr={Pr}, Re={Re}")
    return h
