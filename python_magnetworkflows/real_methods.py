"""
def getCurrent(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Intensity_{marker}_integrate"].iloc[-1]

def getPower(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Power_{marker}_integrate"].iloc[-1]

# Stress
def getMinHoop(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Stress_{marker}_min"].iloc[-1]

def getMeanHoop(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Stress_{marker}_mean"].iloc[-1]

def getMaxHoop(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Stress_{marker}_max"].iloc[-1]

def getMinVonMises(df: pd.DataFrame, marker: str):
    return df[f"Statistics_VonMises_{marker}_min"].iloc[-1]

def getMeanVonMises(df: pd.DataFrame, marker: str):
    return df[f"Statistics_VonMises_{marker}_mean"].iloc[-1]

def getMaxVonMises(df: pd.DataFrame, marker: str):
    return df[f"Statistics_VonMises_{marker}_max"].iloc[-1]

# Temperature
def getMeanT(df: pd.DataFrame, marker: str):
    return df[f"Statistics_MeanT_{marker}_mean"].iloc[-1]

def getMaxT(df: pd.DataFrame, marker: str):
    return df[f"Statistics_MaxT_{marker}_max"].iloc[-1]

def getMinT(df: pd.DataFrame, marker: str):
    return df[f"Statistics_MaxT_{marker}_min"].iloc[-1]

def getFlux(df: pd.DataFrame, marker: str):
    return df[f"Statistics_Flux_{marker}_integrate"].iloc[-1]
"""
import pandas as pd
from typing import List
from .waterflow import waterflow, rho, Cp

# For Heat exchange


def getDT(
    flow: float, Power: float, Tw: float, dTw: float, P: float, relax: float = 0.0
) -> float:
    # compute dT as Power / rho *Cp * Flow(I)
    DT = Power / (rho(Tw, P) * Cp(Tw, P) * flow)
    # print(
    #     f"GETDT={DT} objectif: {objectif}, flow: {flow}, Power: {Power}, Tw: {Tw}, P: {P}"
    # )
    return (1 - relax) * DT + relax * dTw


def getHeatCoeff(
    waterflow: waterflow, Dh: float, U: float, Tw: float, hw: float, relax: float = 0.0
):
    # compute h as Montgomery()
    # P = pressure(objectif)
    # dTw = setDT(objectif, Power, Tw, P)
    return (1 - relax) * waterflow.montgomery(Tw, U, Dh) + relax * hw


# Temperature
def getMeanT(df: pd.DataFrame, marker: str):
    return df[f"Statistics_MeanT_{marker}_mean"].iloc[-1]


def getMaxT(df: pd.DataFrame, marker: str):
    return df[f"Statistics_MaxT_{marker}_max"].iloc[-1]


def getMinT(df: pd.DataFrame, marker: str):
    return df[f"Statistics_MaxT_{marker}_min"].iloc[-1]


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
