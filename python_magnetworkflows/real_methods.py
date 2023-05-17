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

from .waterflow import waterflow, rho, Cp

# For Heat exchange


def getDT(
    objectif: float, waterflow: waterflow, Power: float, Tw: float, P: float
) -> float:
    # compute dT as Power / rho *Cp * Flow(I)
    return Power / (rho(Tw, P) * Cp(Tw, P) * waterflow.flow(objectif))


def getHeatCoeff(waterflow: waterflow, Dh: float, U: float, Tw: float):
    # compute h as Montgomery()
    # P = pressure(objectif)
    # dTw = setDT(objectif, Power, Tw, P)
    return waterflow.montgomery(Tw, U, Dh)
