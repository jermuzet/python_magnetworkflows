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

from .real_methods import update_U, getCurrent, setCurrent
from .real_methods import getPower, setPower
from .real_methods import getMinHoop, setMinHoop, getMeanHoop, setMeanHoop, getMaxHoop, setMaxHoop
from .real_methods import getMinVonMises, setMinVonMises, getMeanVonMises, setMeanVonMises, getMaxVonMises, setMaxVonMises
from .real_methods import getMinT, setMinT, getMeanT, setMeanT, getMaxT, setMaxT
from .real_methods import getFlux, setFlux
from .real_methods import getHeatCoeff, setHeatCoeff
from .real_methods import getDT, setDT

# TODO create/modify targetdefs on the fly
# add target defs for Bitter: replace H by B
# add target defs for Supra?

targetdefs = {
    "B0": {
        "csv": 'magnetic.measures/values.csv',
        "rematch": 'Points_B0_expr_Bz',
        "params": [],
        "control_params": [],
        "value": (),
        "unit": "T"
        },
    "Inductance": {
        "csv": 'magnetic.measures/values.csv',
        "rematch": 'Statistics_MagneticEnergy_integrate',
        "params": [],
        "control_params": [],
        "value": (),
        "unit": "T"
        },
    "I": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_Intensity_H\w+_integrate',
        "params": [('N','N_H\w+')],
        "control_params": [('U', 'U_H\w+', update_U)],
        "value": (getCurrent, setCurrent),
        "unit": "A"
        },
    "PowerH": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_Power_\w+_integrate',
        "params": [],
        "control_params": [],
        "value": (getPower, setPower),
        "unit": "W"
    },
    "Power": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_Power_integrate',
        "params": [],
        "control_params": [],
        "value": (getPower, setPower),
        "unit": "W"
    },
    "MinHoop": {
        "csv": 'elastic.measures/values.csv',
        "rematch": 'Statistics_Stress_H\w+_min',
        "params": [],
        "control_params": [],
        "value": (getMinHoop, setMinHoop),
        "unit": "Pa"
    },
    "MinVonMises": {
        "csv": 'elastic.measures/values.csv',
        "rematch": 'Statistics_VonMises_H\w+_min',
        "params": [],
        "control_params": [],
        "value": (getMinVonMises, setMinVonMises),
        "unit": "Pa"
    },
    "MeanHoop": {
        "csv": 'elastic.measures/values.csv',
        "rematch": 'Statistics_Stress_H\w+_mean',
        "params": [],
        "control_params": [],
        "value": (getMeanHoop, setMeanHoop),
        "unit": "Pa"
    },
    "MeanVonMises": {
        "csv": 'elastic.measures/values.csv',
        "rematch": 'Statistics_VonMises_H\w+_mean',
        "params": [],
        "control_params": [],
        "value": (getMeanVonMises, setMeanVonMises),
        "unit": "Pa"
    },
    "MaxHoop": {
        "csv": 'elastic.measures/values.csv',
        "rematch": 'Statistics_Stress_H\w+_max',
        "params": [],
        "control_params": [],
        "value": (getMaxHoop, setMaxHoop),
        "unit": "Pa"
    },
    "MaxVonMises": {
        "csv": 'elastic.measures/values.csv',
        "rematch": 'Statistics_VonMises_H\w+_max',
        "params": [],
        "control_params": [],
        "value": (getMaxVonMises, setMaxVonMises),
        "unit": "Pa"
    },
    "MinTH": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_T_H\w+_min',
        "params": [],
        "control_params": [],
        "value": (getMinT, setMinT),
        "unit": "K"
    },
    "MinT": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_Stat_T_min',
        "params": [],
        "control_params": [],
        "value": (getMinT, getMinT),
        "unit": "K"
    },
    "MeanTH": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_T_H\w+_mean',
        "params": [],
        "control_params": [],
        "value": (getMeanT, setMeanT),
        "unit": "K"
    },
    "MeanT": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_Stat_T_mean',
        "params": [],
        "control_params": [],
        "value": (getMeanT, getMeanT),
        "unit": "K"
    },
    "MaxTH": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_T_H\w+_max',
        "params": [],
        "control_params": [],
        "value": (getMaxT, setMaxT),
        "unit": "K"
    },
    "MaxT": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_Stat_T_max',
        "params": [],
        "control_params": [],
        "value": (getMaxT, getMaxT),
        "unit": "K"
    },
    "Flux": {
        "csv": 'heat.measures/values.csv',
        "rematch": 'Statistics_Flux_Channel\d+_integrate',
        "params": [],
        "control_params": [],
        "value": (getFlux, setFlux),
        "unit": "W"
    },
    "HeatCoeff": {
        "csv": '',
        "rematch": '',
        "params": [('Dh','Dh\d+'), ('Sh','Sh\d+'), ('hw','hw'),('h','h\d+')],
        "control_params": [],
        "value": (getHeatCoeff, setHeatCoeff),
        "unit": "W/m2/K"
    },
    "DT": {
        "csv": '',
        "rematch": '',
        "params": [('Tw','Tw'), ('Tw','Tw\d+'), ('dTw','dTw\d+')],
        "control_params": [],
        "value": (getDT, setDT),
        "unit": "K"
    },
}

def setTarget(name: str, params: dict, objectif: float, debug: bool = False):
    # print(f"setTarget: workingdir={ os.getcwd() } name={name}")
    targets = {}
    for key in params:
        if targetdefs[name]['value']:
            I_target = targetdefs[name]['value'][1](key, params, objectif)
            if debug:
                print(f"{name} objectif={objectif}, setvalue={I_target}")
            targets[key] = I_target

    if debug: print(f"targets: {targets}")
    return targets

def getTargetUnit(name: str, debug: bool = False) -> str:
    # print(f"getTargetUnit: workingdir={ os.getcwd() } name={name}")
    defs = targetdefs[name]
    return defs['unit']

def getTarget(name: str, e, debug: bool = False) -> pd.DataFrame:
    # print(f"getTarget: workingdir={ os.getcwd() } name={name}")

    defs = targetdefs[name]
    if debug:
        print(f"defs: {defs}")
        print(f"csv: {defs['csv']}")
        print(f"rematch: {defs['rematch']}")

    filename = defs['csv']
    try:
        with open(filename, "r") as f:
            if debug: print(f"csv: {f.name}")
            filtered_df = post(f.name, defs['rematch'], debug)
    except:
        return pd. DataFrame()

    if debug and e.isMasterRank():
        print(filtered_df)
        for key in filtered_df.columns.values.tolist():
            print(key)

    return filtered_df

def Merge(dict1: dict, dict2: dict, debug: bool = False) -> dict:

    if debug :
        print(f"dict1: {dict1}")
        print(f"dict2: {dict2}")

    if isinstance(dict2, type(None)):
        return dict1

    for key1 in dict1:
        if key1 in dict2:
            dict2[key1].update(dict1[key1])
        else:
            dict2[key1] = dict1[key1]

    if debug :
        print(f"dict1: {dict1}")
        print(f"dict2: {dict2}")

    if debug :
        print(f"res: {dict2}")
    return dict2

def getparam(param:str, parameters: dict, rmatch: str, debug: bool = False ):
    """
    """
    if debug:
        print(f"getparam: {param} ====== Start")

    n = 0
    val = {}

    import re
    regex_match = re.compile(rmatch)
    for p in parameters.keys() :
        if regex_match.fullmatch(p):
            marker = p.split(param + '_')[-1]
            if debug:
                print(f"match {p}: {marker}")
            val[marker] = { param: parameters[p]}
            if debug:
                print(f"{p}: {parameters[p]}")

    if debug:
        print(f"val: {val}")
        print(f"getparam: {param} ====== Done")
    return (val)


def post(csv: str, rmatch: str, debug: bool = False):
    """
    extract data for csv result files

    eg:
    rmatch= "Intensity_\w+_integrate"
    csv = ["cfpdes.heat.measures.csv", "cfpdes.magnetic.measures.csv"]
    """
    if debug:
        print(f"post: workingdir={ os.getcwd() }")
        print(f"post: csv={csv}")

    # Retreive current intensities
    df = pd.DataFrame()
    if debug: print("post: loading {csv_}")
    with open(csv, 'r') as f:
        _df = pd.read_csv(f, sep=",", engine='python')
        if debug:
            for key in _df.columns.values.tolist():
                print(key)

        tmp_df = _df.filter(regex=(rmatch))
        if debug: print(f"tmp_df: {tmp_df}")

        df = pd.concat([df, tmp_df], axis='columns')

    return df

def update(cwd: str, name: str, jsonmodel: str, paramsdict: dict, params: List[str], bcparams: dict, objectif: float, debug: bool=False):
    # Update tensions U
    import re

    pwd = os.getcwd()
    os.chdir(cwd)
    if debug:
        print(f"update: workingdir={ os.getcwd() }")

    with open(jsonmodel, 'r') as jsonfile:
        dict_json = json.loads(jsonfile.read())
        parameters = dict_json['Parameters']

    for key in paramsdict:
        for p in params:
            if debug:
                print(f"param: {p}")
                print(f"init {p}_{key} = {parameters[f'{p}_{key}']}")
                print(f"after {p}_{key} = {paramsdict[key][p]}")
            parameters[f'{p}_{key}'] = paramsdict[key][p]

    for key in bcparams[name]:
        if debug:
            print(f"param: {key}")
            print(f"init {key} = {parameters[f'{key}']}")
            print(f"after {key} = {bcparams[name][key]['value']}")
        parameters[key] = bcparams[name][key]['value']

    new_name_json =  jsonmodel.replace('.json', f'-{name}{str(objectif)}A.json')

    with open(new_name_json, 'w+') as jsonfile:
        jsonfile.write(json.dumps(dict_json, indent=4))

    # cfg = jsonmodel.replace(".json", ".cfg")
    # new_cfg =  jsonmodel.replace('.json', f'-I{str(objectif)}A.cfg')
    # cregexp = re.compile(r"^directory=\.$")
    # jregexp = re.compile(r".json")
    # with open(cfg, 'r') as cfgfile:
    #     content = cfgfile.read()
    # print("content:", cregexp.sub(f'-I{str(objectif)}A', content))
    # print("content:", jregexp.sub(f'-I{str(objectif)}A.json', content))

    os.chdir(pwd)
    return new_name_json
