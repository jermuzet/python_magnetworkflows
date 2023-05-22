"""
Define targets for Axi simulation

structure of a target:
name: 'I'
csv: name of the csv file where name data are recorded
rematch: regexp to recover name data from csv
params from parameters section,
control_params from parameters section,
value: name of the method to compute name
"""

from typing import List, Union, Optional, Type

import os

import json
import pandas as pd

# TODO create/modify targetdefs on the fly
# add target defs for Bitter: replace H by B
# add target defs for Supra?

"""
targetdefs = {
    "B0": {
        "csv": "magnetic.measures/values.csv",
        "rematch": "Points_B0_expr_Bz",
        "params": [],
        "control_params": [],
        "value": (),
        "unit": "T",
    },
    "Inductance": {
        "csv": "magnetic.measures/values.csv",
        "rematch": "Statistics_MagneticEnergy_integrate",
        "params": [],
        "control_params": [],
        "value": (),
        "unit": "T",
    },
    "MinHoop": {
        "csv": "elastic.measures/values.csv",
        "rematch": "Statistics_Stress_H\\w+_min",
        "params": [],
        "control_params": [],
        "value": (getMinHoop),
        "unit": "Pa",
    },
    "MinVonMises": {
        "csv": "elastic.measures/values.csv",
        "rematch": "Statistics_VonMises_H\\w+_min",
        "params": [],
        "control_params": [],
        "value": (getMinVonMises),
        "unit": "Pa",
    },
    "MeanHoop": {
        "csv": "elastic.measures/values.csv",
        "rematch": "Statistics_Stress_H\\w+_mean",
        "params": [],
        "control_params": [],
        "value": (getMeanHoop),
        "unit": "Pa",
    },
    "MeanVonMises": {
        "csv": "elastic.measures/values.csv",
        "rematch": "Statistics_VonMises_H\\w+_mean",
        "params": [],
        "control_params": [],
        "value": (getMeanVonMises),
        "unit": "Pa",
    },
    "MaxHoop": {
        "csv": "elastic.measures/values.csv",
        "rematch": "Statistics_Stress_H\\w+_max",
        "params": [],
        "control_params": [],
        "value": (getMaxHoop),
        "unit": "Pa",
    },
    "MaxVonMises": {
        "csv": "elastic.measures/values.csv",
        "rematch": "Statistics_VonMises_H\\w+_max",
        "params": [],
        "control_params": [],
        "value": (getMaxVonMises),
        "unit": "Pa",
    },
    "MinTH": {
        "csv": "heat.measures/values.csv",
        "rematch": "Statistics_T_H\\w+_min",
        "params": [],
        "control_params": [],
        "value": (getMinT),
        "unit": "K",
    },
    "MinT": {
        "csv": "heat.measures/values.csv",
        "rematch": "Statistics_Stat_T_min",
        "params": [],
        "control_params": [],
        "value": (getMinT),
        "unit": "K",
    },
    "MeanTH": {
        "csv": "heat.measures/values.csv",
        "rematch": "Statistics_T_H\\w+_mean",
        "params": [],
        "control_params": [],
        "value": (getMeanT),
        "unit": "K",
    },
    "MeanT": {
        "csv": "heat.measures/values.csv",
        "rematch": "Statistics_Stat_T_mean",
        "params": [],
        "control_params": [],
        "value": (getMeanT),
        "unit": "K",
    },
    "MaxTH": {
        "csv": "heat.measures/values.csv",
        "rematch": "Statistics_T_H\\w+_max",
        "params": [],
        "control_params": [],
        "value": (getMaxT),
        "unit": "K",
    },
    "MaxT": {
        "csv": "heat.measures/values.csv",
        "rematch": "Statistics_Stat_T_max",
        "params": [],
        "control_params": [],
        "value": (getMaxT),
        "unit": "K",
    },
}


def setTarget(targetdefs: dict, name: str, params: dict, objectif: float, debug: bool = False):
    # print(f"setTarget: workingdir={ os.getcwd() } name={name}")
    targets = {}
    for key in params:
        if targetdefs[name]["value"]:
            I_target = targetdefs[name]["value"][1](key, params, objectif)
            if debug:
                print(f"{name} objectif={objectif}, setvalue={I_target}")
            targets[key] = I_target

    if debug:
        print(f"targets: {targets}")
    return targets


def getTargetUnit(targetdefs: dict, name: str, debug: bool = False) -> str:
    # print(f"getTargetUnit: workingdir={ os.getcwd() } name={name}")
    defs = targetdefs[name]
    return defs["unit"]
"""


def getTarget(targetdefs: dict, name: str, e, debug: bool = False) -> pd.DataFrame:
    # print(f"getTarget: workingdir={ os.getcwd() } name={name}")

    defs = targetdefs[name]
    if debug:
        print(f"defs: {defs}")
        print(f"csv: {defs['csv']}")
        print(f"rematch: {defs['rematch']}")

    filename = defs["csv"]
    try:
        with open(filename, "r") as f:
            if debug:
                print(f"csv: {f.name}")
            filtered_df = post(f.name, defs["rematch"], debug)

        # rename columns
        dict_columns = {}
        for col in filtered_df.columns:
            cname = col.replace(f'{defs["post"]["type"]}_', "")
            cname = cname.replace(f'_{defs["post"]["math"]}', "")
            dict_columns[col] = cname
        filtered_df.rename(columns=dict_columns, inplace=True)

    except:
        return pd.DataFrame()

    if debug and e.isMasterRank():
        print(filtered_df)
        for key in filtered_df.columns.values.tolist():
            print(key)

    return filtered_df


def getparam(param: str, parameters: dict, rmatch: str, debug: bool = False) -> list:
    """ """
    if debug:
        print(f"getparam: {param} ====== Start")

    # print(f'getparams: param={param}, rmatch={rmatch}, parameters={parameters}')
    n = 0
    val = []

    import re

    regex_match = re.compile(rmatch)
    for p in parameters.keys():
        if regex_match.fullmatch(p):
            # marker = p.split(param + "_")[-1]
            if debug:
                # print(f"match {p}: {marker}")
                print(f"{p}: {parameters[p]}")
            val.append(p)

    if debug:
        print(f"val: {val}")
        print(f"getparam: {param} ====== Done")
    return val


def post(csv: str, rmatch: str, debug: bool = False):
    """
    extract data for csv result files

    eg:
    rmatch= "Intensity_\\w+_integrate"
    csv = ["cfpdes.heat.measures.csv", "cfpdes.magnetic.measures.csv"]
    """
    if debug:
        print(f"post: workingdir={ os.getcwd() }")
        print(f"post: csv={csv}")

    # Retreive current intensities
    df = pd.DataFrame()
    if debug:
        print("post: loading {csv_}")
    with open(csv, "r") as f:
        _df = pd.read_csv(f, sep=",", engine="python")
        if debug:
            for key in _df.columns.values.tolist():
                print(key)

        tmp_df = _df.filter(regex=(rmatch))
        if debug:
            print(f"tmp_df: {tmp_df}")

        df = pd.concat([df, tmp_df], axis="columns")

    return df
