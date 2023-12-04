"""
Run feelpp model for one config
"""
from typing import List

import os
from math import fabs

import json

import re
import pandas as pd
from tabulate import tabulate

from .params import getparam
from .solver import solve

# from ..units import load_units


def oneconfig(
    fname: str,
    e,
    f,
    feelpp_directory: str,
    jsonmodel: str,
    meshmodel: str,
    args,
    targets: dict,
    postvalues: dict,
    parameters: dict,
):
    """
     Run a simulation until currents are reached

     e: feelpp environment
     jsonmodel:
     meshmodel:
     targets: dict of target (name: objectif , csv, ...)
     postvalues:
     parameters: all jsonmodel parameters

    returns a tuple:
    table:
    dict_df:
    e:
    """

    if e.isMasterRank():
        print("\n\nONECONFIG START", flush=True)

    table_values = []
    table_headers = []

    # pwd = os.getcwd()
    basedir = os.path.dirname(args.cfgfile)
    if e.isMasterRank():
        print(f"oneconfig: jsonmodel={jsonmodel}, basedir={basedir}", flush=True)

    dict_df = {}
    for target, values in targets.items():
        if e.isMasterRank():
            print(f"{target}: {values['objectif']}")
        table_values.append(float(values["objectif"]))
        table_headers.append(f'{target}[{values["unit"]}]')
        dict_df[target] = {
            "PowerM": pd.DataFrame(),
            "PowerH": pd.DataFrame(),
            "Flux": pd.DataFrame(),
            "HeatCoeff": pd.DataFrame(),
            "DT": pd.DataFrame(),
            "statsT": {
                "MinT": pd.DataFrame(),
                "MaxT": pd.DataFrame(),
                "MeanT": pd.DataFrame(),
            },
            "statsTH": {
                "MinTH": pd.DataFrame(),
                "MaxTH": pd.DataFrame(),
                "MeanTH": pd.DataFrame(),
            },
            "target": values["objectif"],
            "flow": 0,
            "Tout": 0,
        }

    # capture actual params per target:
    params = {}
    for key, values in targets.items():
        params[key] = []
        for p in values["control_params"]:
            if args.debug and e.isMasterRank():
                print(f"extract control params for {p[0]}")
            tmp = getparam(p[0], parameters, p[1], args.debug)
            params[key] += tmp

    (table, dict_df, e) = solve(
        fname,
        e,
        f,
        feelpp_directory,
        jsonmodel,
        meshmodel,
        args,
        targets,
        postvalues,
        params,
        parameters,
        dict_df,
    )
    for target, values in dict_df.items():
        if args.debug and e.isMasterRank():
            print(f"\n\nresult for {target}:", flush=True)
        for key, df in values.items():
            if isinstance(df, pd.DataFrame):
                if args.debug and e.isMasterRank():
                    print(f"\t{key}:")
                df["I"] = f'I={dict_df[target]["target"]}A'
                df.set_index("I", inplace=True)
                if args.debug and e.isMasterRank():
                    print(df)
                # df.to_markdown(tablefmt="psql") # requires pqndqs >= 1.0.0
            if key in ["statsT", "statsTH"]:
                for keyT, dfT in df.items():
                    dfT["I"] = f'{keyT}_I={dict_df[target]["target"]}A'
                    dfT.set_index("I", inplace=True)

    return (table, dict_df, e)
