import json
import pandas as pd
import math

import warnings
from pint import UnitRegistry, Unit, Quantity
# Ignore warning for pint
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])

# Pint configuration
ureg = UnitRegistry()
ureg.default_system = 'SI'
ureg.autoconvert_offset_to_baseunit = True

Vpump0 = 1000 # rpm
Vpmax = 2840 # rpm
F0_l_per_second = 0 # l/s
Fmax_l_per_second = 140 # l/s
Pmax = 22 # bar
Pmin = 4 # bar
Imax = 28000 # A

# Flow params
def flow_params(filename: str):
    with open(filename, 'r') as f:
        print(f"Load flow params from {f.name}")
        flow_params = json.loads(f.read())

    global Vpump0, Vpmax, Fmax_l_per_second, Pmax, Pmin, Imax
    
    Vpump0 = flow_params['Vp0']['value'] # rpm
    Vpmax = flow_params['Vpmax']['value'] # rpm
    F0_l_per_second = flow_params['F0']['value'] # l/s
    Fmax_l_per_second = flow_params['Fmax']['value'] # l/s
    Pmax =  flow_params['Pmax']['value'] # bar
    Pmin =  flow_params['Pmin']['value'] # bar
    Imax =  flow_params['Imax']['value'] # Amperes

    pass

def update_U(params: dict, marker: str, target: float, val: float):
    return float(params[marker]['U']) * target/val
    pass

def getCurrent(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Intensity_{marker}_integrate"].iloc[-1]

def setCurrent(marker: str, params: dict, objectif):
    return float(params[marker]['N']) * objectif

def getPower(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Power_{marker}_integrate"].iloc[-1]

def setPower():
    pass

def getMeanT(df: pd.DataFrame, marker: str):
    return df[f"Statistics_MeanT_{marker}_mean"].iloc[-1]

def setMeanT():
    pass

def getMaxT(df: pd.DataFrame, marker: str):
    return df[f"Statistics_MaxT_{marker}_max"].iloc[-1]

def setMaxT():
    pass

def getFlux(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Flux_{marker}_integrate"].iloc[-1]

# For Heat exchange

def setFlux():
    pass

def vpump(objectif: float) -> float:
    Vpump = Vpmax
    Vpump0 = 1000
    if objectif <= Imax:
        Vpump = Vpump0+(Vpmax-Vpump0)*(objectif/Imax)**2

    return Vpump

def flow(objectif: float) -> float: 
    """
    compute flow in m^3/s
    """

    units = [ ureg.liter/ureg.second, ureg.meter*ureg.meter*ureg.meter/ureg.second]
    F0 = Quantity(F0, units[0]).to(units[1]).magnitude
    Fmax = Quantity(Fmax_l_per_second, units[0]).to(units[1]).magnitude
    return F0 + Fmax * vpump(objectif)/Vpmax

def pressure(objectif: float) -> float:
    """
    compute pressure in bar ???
    """
    return (Pmin + (Pmax-Pmin) * (vpump(objectif)/Vpmax)**2)

def umean(objectif: float, section: float) -> float:
    """
    compute umean in m/s ???
    """
    # print("flow:", flow(objectif), section)
    return flow(objectif)/section

def rho(Tw: float, P: float) -> float:
    """
    compute water volumic mass in ???
    TODO link with freesteam
    """
    return 1.e+3

def Cp(Tw: float, P: float) -> float:
    """
    compute water specific heat in ???
    TODO link with freesteam
    """
    return 2840

def montgomery(Tw: float, Umean: float, Dh: float) -> float:
    """
    compute heat exchange coefficient in ??

    Tw: K
    Umean: m/s
    Dh: meter
    """
    
    return 1426*(1+1.5e-2*(Tw-273))*math.exp(math.log(Umean)*0.8)/math.exp(math.log(Dh)*0.2)

def getDT(objectif: float, Power: float, Tw: float, P: float) -> float:
    # compute dT as Power / rho *Cp * Flow(I)
    return Power/ (rho(Tw, P) * Cp(Tw, P) * flow(objectif))

def setDT():
    pass

def getHeatCoeff(Dh: float, Umean: float, Tw: float):
    # compute h as Montgomery()
    # P = pressure(objectif)
    # dTw = setDT(objectif, Power, Tw, P) 
    return montgomery(Tw, Umean, Dh)

def setHeatCoeff():
    pass


