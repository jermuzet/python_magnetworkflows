"""
Run feelpp model for one config
"""

from typing import List

import os
import pandas as pd

from .params import getparam
from .solver import solve

# from ..units import load_units


def oneconfig(
    fname: str,
    e,
    f,
    fields,
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
     postvalues: dict for postprocessing quantities
     parameters: all jsonmodel parameters

    returns a tuple:
    table:
    dict_df:
    e:
    """
    comm = e.worldCommPtr()

    if e.isMasterRank():
        print("\n\nONECONFIG START", flush=True)

    # pwd = os.getcwd()
    basedir = os.path.dirname(args.cfgfile)
    if e.isMasterRank():
        print(f"oneconfig: jsonmodel={jsonmodel}, basedir={basedir}", flush=True)

    dict_df = {}
    table_values = []
    table_headers = []
    e.worldComm().barrier()
    if e.isMasterRank():
        for target, values in targets.items():
            print(f"{target}: {values['objectif']}", flush=True)
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
                "target": float(values["objectif"]),
                "flow": 0.0,
                "Tout": 0.0,
                "Uw": pd.DataFrame(),
                "L": float(values["inductance"]),
            }

            if "thmagel" in args.cfgfile:
                dict_df[target]["statsDispl"] = {
                    "MinDispl": pd.DataFrame(),
                    "MaxDispl": pd.DataFrame(),
                    "MeanDispl": pd.DataFrame(),
                }
                dict_df[target]["statsStress"] = {
                    "MinStress": pd.DataFrame(),
                    "MaxStress": pd.DataFrame(),
                    "MeanStress": pd.DataFrame(),
                }
                dict_df[target]["statsVonMises"] = {
                    "MinVonMises": pd.DataFrame(),
                    "MaxVonMises": pd.DataFrame(),
                    "MeanVonMises": pd.DataFrame(),
                }
                dict_df[target]["statsDisplH"] = {
                    "MinDisplH": pd.DataFrame(),
                    "MaxDisplH": pd.DataFrame(),
                    "MeanDisplH": pd.DataFrame(),
                }
                dict_df[target]["statsStressH"] = {
                    "MinStressH": pd.DataFrame(),
                    "MaxStressH": pd.DataFrame(),
                    "MeanStressH": pd.DataFrame(),
                }
                dict_df[target]["statsVonMisesH"] = {
                    "MinVonMisesH": pd.DataFrame(),
                    "MaxVonMisesH": pd.DataFrame(),
                    "MeanVonMisesH": pd.DataFrame(),
                }

    table_values = comm.localComm().bcast(table_values, root=0)
    table_headers = comm.localComm().bcast(table_headers, root=0)
    dict_df = comm.localComm().bcast(dict_df, root=0)

    # capture actual params per target:
    params = {}
    for key, values in targets.items():
        params[key] = []
        for p in values["control_params"]:
            if args.debug and e.isMasterRank():
                print(f"extract control params for {p[0]}")
            tmp = getparam(p[0], parameters, p[1], args.debug)
            params[key] += tmp

    e.worldComm().barrier()
    (table, dict_df, e) = solve(
        fname,
        e,
        f,
        fields,
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

    if e.isMasterRank():
        print("solve done")

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
            if key in [
                "statsT",
                "statsTH",
                "statsDispl",
                "statsStress",
                "statsDisplH",
                "statsStressH",
                "statsVonMises",
                "statsVonMisesH",
            ]:
                for keyT, dfT in df.items():
                    dfT["I"] = f'{keyT}_I={dict_df[target]["target"]}A'
                    dfT.set_index("I", inplace=True)

    if e.isMasterRank():
        print("end of oneconfig")

    e.worldComm().barrier()
    return (table, dict_df, e)
