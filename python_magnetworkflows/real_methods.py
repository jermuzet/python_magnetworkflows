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
from .waterflow import waterflow
from .cooling import rho, Cp
from .cooling import Montgomery, Dittus, Colburn, Silverberg

# For Heat exchange


def getDT(
    objectif: float, waterflow: waterflow, Power: float, Tw: float, P: float
) -> float:
    # compute dT as Power / rho *Cp * Flow(I)
    return Power / (rho(Tw, P) * Cp(Tw, P) * waterflow.flow(objectif))


def getHeatCoeff(
    Dh: float,
    L: float,
    U: float,
    Tw: float,
    Pw: float,
    dPw: float,
    model: str = "Montgomery",
):
    correlation = {
        "Montgomery": Montgomery,
        "Dittus": Dittus,
        "Colburn": Colburn,
        "Silverberg": Silverberg,
    }
    return correlation[model](Tw, Pw, dPw, U, Dh, L)


def getTout(
    T: List[float], VolMass: List[float], SpecHeat: List[float], Q: List[float]
) -> float:
    Tout = 0
    rhoCpQ = 0
    for i, (Ti, RHOi, CPi, Qi) in enumerate(zip(T, VolMass, SpecHeat, Q)):
        Tout += Ti * RHOi * CPi * Qi
        rhoCpQ += RHOi * CPi * Qi

    Tout /= rhoCpQ
    return Tout


# Temperature
def getMean(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Mean{field}_{marker}_mean"].iloc[-1]


def getMax(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Max{field}_{marker}_max"].iloc[-1]


def getMin(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Max{field}_{marker}_min"].iloc[-1]
