"""
Run feelpp model
"""

import sys
import os
import argparse
import configparser

import pandas as pd
import json

from .waterflow import waterflow
from .real_methods import getDT, getHeatCoeff
from .real_methods import getMinT, getMeanT, getMaxT

from .oneconfig import oneconfig


def main():
    # print(f'sys.argv({type(sys.argv)})={sys.argv}')

    epilog = (
        "Setup for cfpdes static formulations\n"
        "Support: only magnet of insert type\n"
        "Workflow: actually fix current and compute cooling BCs using Montgomery correlation with an average water velocity\n"
        "\n"
        "Before running adapt flow_params to your magnet setup \n"
    )

    command_line = None
    parser = argparse.ArgumentParser(description="Cfpdes model", epilog=epilog)
    parser.add_argument("cfgfile", help="input cfg file (ex. HL-31.cfg)")

    # TODO: make a group: oneconfig
    # TODO make current a dict: magnet_name: value, targetkey: value
    parser.add_argument("--mdata", help="specify current data", type=json.loads)

    parser.add_argument(
        "--cooling",
        help="choose cooling type",
        type=str,
        choices=["mean", "grad", "meanH", "gradH"],
        default="mean",
    )
    parser.add_argument(
        "--eps",
        help="specify requested tolerance (default: 1.e-3)",
        type=float,
        default=1.0e-3,
    )
    parser.add_argument(
        "--itermax",
        help="specify maximum iteration (default: 10)",
        type=int,
        default=10,
    )
    parser.add_argument("--debug", help="activate debug", action="store_true")
    parser.add_argument("--verbose", help="activate verbose", action="store_true")

    args = parser.parse_args()
    if args.debug:
        print(args)

    pwd = os.getcwd()
    print(f"cwd: {pwd}")

    # Load units:
    # TODO force millimeter when args.method == "HDG"
    # units = load_units('meter')

    # Load cfg as config
    jsonmodel = ""
    meshmodel = ""
    feelpp_config = configparser.ConfigParser()
    with open(args.cfgfile, "r") as inputcfg:
        feelpp_config.read_string("[DEFAULT]\n[main]\n" + inputcfg.read())
        feelpp_directory = feelpp_config["main"]["directory"]
        print(f"feelpp_directory={feelpp_directory}")

        basedir = os.path.dirname(args.cfgfile)
        print(f"basedir={basedir}")
        if not basedir:
            basedir = "."

        jsonmodel = feelpp_config["cfpdes"]["filename"]
        jsonmodel = jsonmodel.replace(r"$cfgdir/", f"{basedir}/")
        print(f"jsonmodel={jsonmodel}")

        meshmodel = feelpp_config["cfpdes"]["mesh.filename"]
        meshmodel = meshmodel.replace(r"$cfgdir/", f"{basedir}/")
        print(f"meshmodel={meshmodel}")

    # Get Parameters from JSON model file
    parameters = {}
    with open(jsonmodel, "r") as jsonfile:
        dict_json = json.loads(jsonfile.read())
        parameters = dict_json["Parameters"]

    targets = {}
    postvalues = {}

    # args.mdata = currents:  {magnet.name: {'value': current.value, 'type': magnet.type, 'filter': '', 'flow_params': args.flow_params}}
    if args.mdata:
        for mname, values in args.mdata.items():
            print(f"mname={mname}, values={values}")
            filter = values["filter"]
            if values["type"] == "helix":
                # change rematch, params, control_params
                PowerM = {
                    "name": "PowerM",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_PowerM_{filter}\\w*integrate",
                    "post": {"type": "Statistics_PowerM", "math": "integrate"},
                    "unit": "W",
                }
                PowerH = {
                    "name": "PowerH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Power_{filter}H\\d+_integrate",
                    "post": {"type": "Statistics_Power", "math": "integrate"},
                    "unit": "W",
                }

                Flux = {
                    "name": "Flux",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Flux_{filter}Channel\\d+_integrate",
                    "post": {"type": "Statistics_Flux", "math": "integrate"},
                    "unit": "W",
                }

                HeatCoeff = {
                    "name": "HeatCoeff",
                    "params": [
                        ("Dh", f"{filter}Dh\\d+"),
                        ("Sh", f"{filter}Sh\\d+"),
                        ("hw", f"{filter}hw"),
                        ("h", f"{filter}hw\\d+"),
                    ],
                    "value": (getHeatCoeff),
                    "unit": "W/m2/K",
                }

                DT = {
                    "name": "DT",
                    "params": [
                        ("Tw", f"{filter}Tw"),
                        ("dTw", f"{filter}dTw"),
                        ("TwH", f"{filter}Tw\\d+"),
                        ("dTwH", f"{filter}dTw\\d+"),
                    ],
                    "value": (getDT),
                    "unit": "K",
                }
                MinT = {
                    "name": "MinT",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Stat_T_{filter}\\w*min",
                    "params": [],
                    "control_params": [],
                    "value": (getMinT),
                    "unit": "K",
                    "post": {"type": "Statistics_Stat_T", "math": "min"},
                }
                MeanT = {
                    "name": "MeanT",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Stat_T_{filter}\\w*mean",
                    "params": [],
                    "control_params": [],
                    "value": (getMeanT),
                    "unit": "K",
                    "post": {"type": "Statistics_Stat_T", "math": "mean"},
                }
                MaxT = {
                    "name": "MaxT",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Stat_T_{filter}\\w*max",
                    "params": [],
                    "control_params": [],
                    "value": (getMaxT),
                    "unit": "K",
                    "post": {"type": "Statistics_Stat_T", "math": "max"},
                }
                MinTH = {
                    "name": "MinTH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_T_{filter}\\w+\\d+_min",
                    "params": [],
                    "control_params": [],
                    "value": (getMinT),
                    "unit": "K",
                    "post": {"type": "Statistics_T", "math": "min"},
                }
                MeanTH = {
                    "name": "MeanTH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_T_{filter}\\w+\\d+_mean",
                    "params": [],
                    "control_params": [],
                    "value": (getMeanT),
                    "unit": "K",
                    "post": {"type": "Statistics_T", "math": "mean"},
                }
                MaxTH = {
                    "name": "MaxTH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_T_{filter}\\w+\\d+_max",
                    "params": [],
                    "control_params": [],
                    "value": (getMaxT),
                    "unit": "K",
                    "post": {"type": "Statistics_T", "math": "max"},
                }

                targets[f"{filter}I"] = {
                    "objectif": values["value"],
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Intensity_{filter}H\\w+_integrate",
                    "params": [("N", f"N_{filter}\\w+")],
                    "control_params": [(f"{filter}U", f"U_{filter}\\w+")],
                    "computed_params": [HeatCoeff, DT, Flux, PowerM, PowerH],
                    "unit": "A",
                    "name": f"Intensity_{filter}",
                    "post": {"type": "Statistics_Intensity", "math": "integrate"},
                    "waterflow": waterflow.flow_params(values["flow"]),
                }

            if values["type"] == "bitter":
                # change rematch, params, control_params
                PowerM = {
                    "name": "PowerM",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_PowerM_{filter}\\w*integrate",
                    "post": {"type": "Statistics_PowerM", "math": "integrate"},
                    "unit": "W",
                }
                PowerH = {
                    "name": "PowerH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Power_{filter}\\w+_B\\D_integrate",
                    "post": {"type": "Statistics_Power", "math": "integrate"},
                    "unit": "W",
                }

                Flux = {
                    "name": "Flux",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Flux_{filter}\\w+_Slit\\d+_integrate",
                    "post": {"type": "Statistics_Flux", "math": "integrate"},
                    "unit": "W",
                }

                HeatCoeff = {
                    "name": "HeatCoeff",
                    "params": [
                        ("Dh", f"{filter}\\w+Dh"),
                        ("Sh", f"{filter}\\w+Sh"),
                        ("hw", f"{filter}\\w+hw", "\\w+_Slit\\w+", False),
                        ("h", f"{filter}\\w+hw", "\\w+_Slit\\w+", True),
                    ],
                    "value": (getHeatCoeff),
                    "unit": "W/m2/K",
                }

                DT = {
                    "name": "DT",
                    "params": [
                        ("Tw", f"{filter}\\w+_Tw", "\\w+_Slit\\w+", False),
                        ("dTw", f"{filter}\\w+_dTw", "\\w+_Slit\\w+", False),
                        ("TwH", f"{filter}\\w+_Tw", "\\w+_Slit\\w+", True),
                        ("dTwH", f"{filter}\\w+_dTw", "\\w+_Slit\\w+", True),
                    ],
                    "value": (getDT),
                    "unit": "K",
                }
                MinT = {
                    "name": "MinT",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Stat_T_{filter}\\w*min",
                    "params": [],
                    "control_params": [],
                    "value": (getMinT),
                    "unit": "K",
                    "post": {"type": "Statistics_Stat_T", "math": "min"},
                }
                MeanT = {
                    "name": "MeanT",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Stat_T_{filter}\\w*mean",
                    "params": [],
                    "control_params": [],
                    "value": (getMeanT),
                    "unit": "K",
                    "post": {"type": "Statistics_Stat_T", "math": "mean"},
                }
                MaxT = {
                    "name": "MaxT",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Stat_T_{filter}\\w*max",
                    "params": [],
                    "control_params": [],
                    "value": (getMaxT),
                    "unit": "K",
                    "post": {"type": "Statistics_Stat_T", "math": "max"},
                }
                MinTH = {
                    "name": "MinTH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_T_{filter}\\w+_B\\D_min",
                    "params": [],
                    "control_params": [],
                    "value": (getMinT),
                    "unit": "K",
                    "post": {"type": "Statistics_T", "math": "min"},
                }
                MeanTH = {
                    "name": "MeanTH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_T_{filter}\\w+_B\\D_mean",
                    "params": [],
                    "control_params": [],
                    "value": (getMeanT),
                    "unit": "K",
                    "post": {"type": "Statistics_T", "math": "mean"},
                }
                MaxTH = {
                    "name": "MaxTH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_T_{filter}\\w+_B\\D_max",
                    "params": [],
                    "control_params": [],
                    "value": (getMaxT),
                    "unit": "K",
                    "post": {"type": "Statistics_T", "math": "max"},
                }

                targets[f"{filter}I"] = {
                    "objectif": values["value"],
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Intensity_{filter}\\w+_integrate",
                    "params": [("N", f"N_{filter}\\w+")],
                    "control_params": [(f"{filter}U", f"U_{filter}\\w+")],
                    "computed_params": [HeatCoeff, DT, Flux, PowerM, PowerH],
                    "unit": "A",
                    "name": f"Intensity{filter}",
                    "post": {"type": "Statistics_Intensity", "math": "integrate"},
                    "waterflow": waterflow.flow_params(values["flow"]),
                }

            postvalues[f"{filter}I"] = {
                "statsT": [MinT, MeanT, MaxT],
                "statsTH": [MinTH, MeanTH, MaxTH],
            }

    print(f"targets: {targets.keys()}")

    """
    # pvalues: dict for params key->actual target in targets
    pvalues = {
        "Flux" : "Flux",
        "PowerM" : "PowerM",
        "PowerH" : "PowerH",
        "HeatCoeff" : "HeatCoeff",
        "DT" : "DT"
    }
        
    postvalues = {
        "T": ['MinTH', 'MeanTH', 'MaxTH', 'PowerH'],
        "Stress": ['MinVonMises', 'MinHoop', 'MeanVonMises', 'MeanHoop', 'MaxVonMises', 'MaxHoop']
    }
    # pfields depends from method
    pfields = ['I', 'PowerH', 'MinTH', 'MeanTH', 'MaxTH', 'Flux',
                'MinVonMises', 'MeanVonMises', 'MaxVonMises',
                'MinHoop', 'MeanHoop', 'MaxHoop']
    
    """

    (table, dict_df) = oneconfig(
        feelpp_directory, jsonmodel, meshmodel, args, targets, postvalues, parameters
    )

    for target, values in dict_df.items():
        for key, df in values.items():
            if isinstance(df, pd.DataFrame):
                df = df.T
                outdir = f"{target[:-2]}_{key}.measures"
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
                df.to_csv(f"{outdir}/values.csv", index=True)
            if key in ["statsT", "statsTH"]:
                list_dfT = [dfT for keyT, dfT in df.items()]
                dfT = pd.concat(list_dfT).T
                outdir = f"{target[:-2]}_{key}.measures"
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
                dfT.to_csv(f"{outdir}/values.csv", index=True)

    outdir = f"U.measures"
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    table.to_csv(f"{outdir}/values.csv", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
