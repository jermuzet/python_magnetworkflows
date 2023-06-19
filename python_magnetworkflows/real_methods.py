
import pandas as pd
from typing import List
from .waterflow import waterflow
from .cooling import rho, Cp
from .cooling import Montgomery, Dittus, Colburn, Silverberg

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
    Dh: float,
    L: float,
    U: float,
    Tw: float,
    hw: float,
    Pw: float,
    dPw: float,
    model: str = "Montgomery",
    friction: str = "Constant",
    relax: float = 0.0,
):
    correlation = {
        "Montgomery": Montgomery,
        "Dittus": Dittus,
        "Colburn": Colburn,
        "Silverberg": Silverberg,
    }

    h = correlation[model](Tw, Pw, dPw, U, Dh, L, friction)
    return (1 - relax) * h + relax * hw


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


# Temperature
def getMean(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Mean{field}_{marker}_mean"].iloc[-1]


def getMax(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Max{field}_{marker}_max"].iloc[-1]


def getMin(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Max{field}_{marker}_min"].iloc[-1]
