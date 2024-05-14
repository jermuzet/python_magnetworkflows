from typing import List
import pandas as pd


# Temperature
def getMean(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Mean{field}_{marker}_mean"].iloc[-1]


def getMax(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Max{field}_{marker}_max"].iloc[-1]


def getMin(df: pd.DataFrame, marker: str, field: str):
    return df[f"Statistics_Max{field}_{marker}_min"].iloc[-1]
