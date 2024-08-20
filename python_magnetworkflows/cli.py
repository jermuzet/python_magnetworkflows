"""
Run feelpp model
"""

import sys
import os
import argparse
import configparser

import pandas as pd
import json
import re
from warnings import simplefilter

simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

from natsort import natsorted
from tabulate import tabulate

from .waterflow import waterflow
from .cooling import getDT, getHeatCoeff

from .oneconfig import oneconfig
from .solver import init


def options(description: str, epilog: str):
    """
    define options
    """

    command_line = None
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("cfgfile", help="input cfg file (ex. HL-31.cfg)")

    # TODO: make a group: oneconfig
    # TODO make current a dict: magnet_name: value, targetkey: value
    parser.add_argument("--mdata", help="specify current data", type=json.loads)

    parser.add_argument(
        "--cooling",
        help="choose cooling type",
        type=str,
        choices=["mean", "grad", "meanH", "gradH", "gradHZ"],
        default="mean",
    )
    parser.add_argument(
        "--heatcorrelation",
        help="choose cooling model",
        type=str,
        choices=["Montgomery", "Dittus", "Colburn", "Silverberg"],
        default="Montgomery",
    )
    parser.add_argument(
        "--friction",
        help="choose friction method",
        type=str,
        choices=["Constant", "Blasius", "Filonenko", "Colebrook", "Swanee"],
        default="Constant",
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
    parser.add_argument("--reloadcfg", help="get feelpp config", action="store_true")
    parser.add_argument("--debug", help="activate debug", action="store_true")
    parser.add_argument("--verbose", help="activate verbose", action="store_true")

    return parser


def loadMdata(e, pwd: str, args, targets: dict, postvalues: dict):
    for mname, values in args.mdata.items():
        if e.isMasterRank():
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
                    ("Dh", f"Dh_{filter}\\w+"),
                    ("Sh", f"Sh_{filter}\\w+"),
                    ("hw", f"hw_{filter}Channel"),
                    ("hwH", f"hw_{filter}Channel\\d+"),
                    ("Zmax", f"Zmax_{filter}Channel"),
                    ("ZmaxH", f"Zmax_{filter}Channel\\d+"),
                ],
                "value": (getHeatCoeff),
                "unit": "W/m2/K",
            }

            DT = {
                "name": "DT",
                "params": [
                    ("Tw", f"Tw_{filter}Channel"),
                    ("dTw", f"dTw_{filter}Channel"),
                    ("TwH", f"Tw_{filter}Channel\\d+"),
                    ("dTwH", f"dTw_{filter}Channel\\d+"),
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
                "unit": "K",
                "post": {"type": "Statistics_Stat_T", "math": "min"},
            }
            MeanT = {
                "name": "MeanT",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_Stat_T_{filter}\\w*mean",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_Stat_T", "math": "mean"},
            }
            MaxT = {
                "name": "MaxT",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_Stat_T_{filter}\\w*max",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_Stat_T", "math": "max"},
            }
            MinTH = {
                "name": "MinTH",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_T_{filter}\\w+\\d+_min",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_T", "math": "min"},
            }
            MeanTH = {
                "name": "MeanTH",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_T_{filter}\\w+\\d+_mean",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_T", "math": "mean"},
            }
            MaxTH = {
                "name": "MaxTH",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_T_{filter}\\w+\\d+_max",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_T", "math": "max"},
            }

            targets[f"{filter}I"] = {
                "objectif": values["value"],
                "type": "helix",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_Intensity_{filter}H\\w+_integrate",
                "params": [("N", f"N_{filter}\\w+")],
                "control_params": [(f"{filter}U", f"U_{filter}\\w+")],
                "computed_params": [HeatCoeff, DT, Flux, PowerM, PowerH],
                "unit": "A",
                "name": f"Intensity_{filter}",
                "post": {"type": "Statistics_Intensity", "math": "integrate"},
                "waterflow": waterflow.flow_params(f'{pwd}/{values["flow"]}'),
            }
            if "Z" in args.cooling:
                if e.isMasterRank():
                    print("add FluxZ for Insert")
                FluxZ = {
                    "name": "FluxZ",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_FluxZ\\d+_{filter}Channel\\d+_integrate",
                    "post": {"type": "Statistics", "math": "integrate"},
                    "unit": "W",
                }
                targets[f"{filter}I"]["computed_params"].append(FluxZ)

            if "thmagel" in args.cfgfile:
                MinDispl = {
                    "name": "MinDispl",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Displ_{filter}\\w*min",
                    "params": [],
                    "control_params": [],
                    "unit": "m",
                    "post": {"type": "Statistics_Stat_Displ", "math": "min"},
                }
                # MeanDispl = {
                #     "name": "MeanDispl",
                #     "csv": "elastic.measures/values.csv",
                #     "rematch": f"Statistics_Stat_Displ_{filter}\\w*mean",
                #     "params": [],
                #     "control_params": [],
                #     "unit": "m",
                #     "post": {"type": "Statistics_Stat_Displ", "math": "mean"},
                # }
                MaxDispl = {
                    "name": "MaxDispl",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Displ_{filter}\\w*max",
                    "params": [],
                    "control_params": [],
                    "unit": "m",
                    "post": {"type": "Statistics_Stat_Displ", "math": "max"},
                }
                MinStress = {
                    "name": "MinStress",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Stress_{filter}\\w*min",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_Stress", "math": "min"},
                }
                MeanStress = {
                    "name": "MeanStress",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Stress_{filter}\\w*mean",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_Stress", "math": "mean"},
                }
                MaxStress = {
                    "name": "MaxStress",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Stress_{filter}\\w*max",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_Stress", "math": "max"},
                }
                MinVonMises = {
                    "name": "MinVonMises",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_VonMises_{filter}\\w*min",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_VonMises", "math": "min"},
                }
                MeanVonMises = {
                    "name": "MeanVonMises",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_VonMises_{filter}\\w*mean",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_VonMises", "math": "mean"},
                }
                MaxVonMises = {
                    "name": "MaxVonMises",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_VonMises_{filter}\\w*max",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_VonMises", "math": "max"},
                }
                MinDisplH = {
                    "name": "MinDisplH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Displ_{filter}\\w+\\d+_min",
                    "params": [],
                    "control_params": [],
                    "unit": "m",
                    "post": {"type": "Statistics_Displ", "math": "min"},
                }
                # MeanDisplH = {
                #     "name": "MeanDisplH",
                #     "csv": "elastic.measures/values.csv",
                #     "rematch": f"Statistics_Displ_{filter}\\w+\\d+_mean",
                #     "params": [],
                #     "control_params": [],
                #     "unit": "m",
                #     "post": {"type": "Statistics_Displ", "math": "mean"},
                # }
                MaxDisplH = {
                    "name": "MaxDisplH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Displ_{filter}\\w+\\d+_max",
                    "params": [],
                    "control_params": [],
                    "unit": "m",
                    "post": {"type": "Statistics_Displ", "math": "max"},
                }
                MinStressH = {
                    "name": "MinStressH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stress_{filter}\\w+\\d+_min",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stress", "math": "min"},
                }
                MeanStressH = {
                    "name": "MeanStressH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stress_{filter}\\w+\\d+_mean",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stress", "math": "mean"},
                }
                MaxStressH = {
                    "name": "MaxStressH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stress_{filter}\\w+\\d+_max",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stress", "math": "max"},
                }
                MinVonMisesH = {
                    "name": "MinVonMisesH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_VonMises_{filter}\\w+\\d+_min",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_VonMises", "math": "min"},
                }
                MeanVonMisesH = {
                    "name": "MeanVonMisesH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_VonMises_{filter}\\w+\\d+_mean",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_VonMises", "math": "mean"},
                }
                MaxVonMisesH = {
                    "name": "MaxVonMisesH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_VonMises_{filter}\\w+\\d+_max",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_VonMises", "math": "max"},
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
                    ("Dh", f"Dh_{filter}\\w+"),
                    ("Sh", f"Sh_{filter}\\w+"),
                    ("hw", f"hw_{filter}\\w+", "\\w+_Slit\\w+", False),
                    ("hwH", f"hw_{filter}\\w+", "\\w+_Slit\\w+", True),
                    ("Zmax", f"Zmax_{filter}\\w+", "\\w+_Slit\\w+", False),
                    ("ZmaxH", f"Zmax_{filter}\\w+", "\\w+_Slit\\w+", True),
                ],
                "value": (getHeatCoeff),
                "unit": "W/m2/K",
            }

            DT = {
                "name": "DT",
                "params": [
                    ("Tw", f"Tw_{filter}\\w+", "\\w+_Slit\\w+", False),
                    ("dTw", f"dTw_{filter}\\w+", "\\w+_Slit\\w+", False),
                    ("TwH", f"Tw_{filter}\\w+", "\\w+_Slit\\w+", True),
                    ("dTwH", f"dTw_{filter}\\w+", "\\w+_Slit\\w+", True),
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
                "unit": "K",
                "post": {"type": "Statistics_Stat_T", "math": "min"},
            }
            MeanT = {
                "name": "MeanT",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_Stat_T_{filter}\\w*mean",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_Stat_T", "math": "mean"},
            }
            MaxT = {
                "name": "MaxT",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_Stat_T_{filter}\\w*max",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_Stat_T", "math": "max"},
            }
            MinTH = {
                "name": "MinTH",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_T_{filter}\\w+_B\\D_min",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_T", "math": "min"},
            }
            MeanTH = {
                "name": "MeanTH",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_T_{filter}\\w+_B\\D_mean",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_T", "math": "mean"},
            }
            MaxTH = {
                "name": "MaxTH",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_T_{filter}\\w+_B\\D_max",
                "params": [],
                "control_params": [],
                "unit": "K",
                "post": {"type": "Statistics_T", "math": "max"},
            }

            targets[f"{filter}I"] = {
                "objectif": values["value"],
                "type": "bitter",
                "csv": "heat.measures/values.csv",
                "rematch": f"Statistics_Intensity_{filter}\\w+_integrate",
                "params": [("N", f"N_{filter}\\w+")],
                "control_params": [(f"{filter}U", f"U_{filter}\\w+")],
                "computed_params": [HeatCoeff, DT, Flux, PowerM, PowerH],
                "unit": "A",
                "name": f"Intensity{filter}",
                "post": {"type": "Statistics_Intensity", "math": "integrate"},
                "waterflow": waterflow.flow_params(f'{pwd}/{values["flow"]}'),
            }
            if "Z" in args.cooling:
                if e.isMasterRank():
                    print("add FluxZ for Bitter")
                FluxZ = {
                    "name": "FluxZ",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_FluxZ\\d+_{filter}\\w+_Slit\\d+_integrate",
                    "post": {"type": "Statistics", "math": "integrate"},
                    "unit": "W",
                }
                targets[f"{filter}I"]["computed_params"].append(FluxZ)

            if "thmagel" in args.cfgfile:
                MinDispl = {
                    "name": "MinDispl",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Displ_{filter}\\w*min",
                    "params": [],
                    "control_params": [],
                    "unit": "m",
                    "post": {"type": "Statistics_Stat_Displ", "math": "min"},
                }
                # MeanDispl = {
                #     "name": "MeanDispl",
                #     "csv": "elastic.measures/values.csv",
                #     "rematch": f"Statistics_Stat_Displ_{filter}\\w*mean",
                #     "params": [],
                #     "control_params": [],
                #     "unit": "m",
                #     "post": {"type": "Statistics_Stat_Displ", "math": "mean"},
                # }
                MaxDispl = {
                    "name": "MaxDispl",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Displ_{filter}\\w*max",
                    "params": [],
                    "control_params": [],
                    "unit": "m",
                    "post": {"type": "Statistics_Stat_Displ", "math": "max"},
                }
                MinStress = {
                    "name": "MinStress",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Stress_{filter}\\w*min",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_Stress", "math": "min"},
                }
                MeanStress = {
                    "name": "MeanStress",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Stress_{filter}\\w*mean",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_Stress", "math": "mean"},
                }
                MaxStress = {
                    "name": "MaxStress",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_Stress_{filter}\\w*max",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_Stress", "math": "max"},
                }
                MinVonMises = {
                    "name": "MinVonMises",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_VonMises_{filter}\\w*min",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_VonMises", "math": "min"},
                }
                MeanVonMises = {
                    "name": "MeanVonMises",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_VonMises_{filter}\\w*mean",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_VonMises", "math": "mean"},
                }
                MaxVonMises = {
                    "name": "MaxVonMises",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stat_VonMises_{filter}\\w*max",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stat_VonMises", "math": "max"},
                }
                MinDisplH = {
                    "name": "MinDisplH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Displ_{filter}\\w+_B\\D_min",
                    "params": [],
                    "control_params": [],
                    "unit": "m",
                    "post": {"type": "Statistics_Displ", "math": "min"},
                }
                # MeanDisplH = {
                #     "name": "MeanDisplH",
                #     "csv": "elastic.measures/values.csv",
                #     "rematch": f"Statistics_Displ_{filter}\\w+_B\\D_mean",
                #     "params": [],
                #     "control_params": [],
                #     "unit": "m",
                #     "post": {"type": "Statistics_Displ", "math": "mean"},
                # }
                MaxDisplH = {
                    "name": "MaxDisplH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Displ_{filter}\\w+_B\\D_max",
                    "params": [],
                    "control_params": [],
                    "unit": "m",
                    "post": {"type": "Statistics_Displ", "math": "max"},
                }
                MinStressH = {
                    "name": "MinStressH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stress_{filter}\\w+_B\\D_min",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stress", "math": "min"},
                }
                MeanStressH = {
                    "name": "MeanStressH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stress_{filter}\\w+_B\\D_mean",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stress", "math": "mean"},
                }
                MaxStressH = {
                    "name": "MaxStressH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_Stress_{filter}\\w+_B\\D_max",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_Stress", "math": "max"},
                }
                MinVonMisesH = {
                    "name": "MinVonMisesH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_VonMises_{filter}\\w+_B\\D_min",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_VonMises", "math": "min"},
                }
                MeanVonMisesH = {
                    "name": "MeanVonMisesH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_VonMises_{filter}\\w+_B\\D_mean",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_VonMises", "math": "mean"},
                }
                MaxVonMisesH = {
                    "name": "MaxVonMisesH",
                    "csv": "elastic.measures/values.csv",
                    "rematch": f"Statistics_VonMises_{filter}\\w+_B\\D_max",
                    "params": [],
                    "control_params": [],
                    "unit": "Pa",
                    "post": {"type": "Statistics_VonMises", "math": "max"},
                }

        targets[f"{filter}I"]["relax"] = 0
        if "relax" in values:
            targets[f"{filter}I"]["relax"] = values["relax"]

        targets[f"{filter}I"]["inductance"] = 0
        if "inductance" in values:
            targets[f"{filter}I"]["inductance"] = values["inductance"]

        targets[f"{filter}I"]["fuzzy"] = 1.0
        if "heatCorrelationFuzzyFactor" in values:
            targets[f"{filter}I"]["fuzzy"] = values["heatCorrelationFuzzyFactor"]
        elif values["type"] == "bitter":
            targets[f"{filter}I"]["fuzzy"] = 1.7

        postvalues[f"{filter}I"] = {
            "statsT": [MinT, MeanT, MaxT],
            "statsTH": [MinTH, MeanTH, MaxTH],
        }
        if "thmagel" in args.cfgfile:
            postvalues[f"{filter}I"]["statsDispl"] = [
                MinDispl,
                # MeanDispl,
                MaxDispl,
            ]
            postvalues[f"{filter}I"]["statsStress"] = [MinStress, MeanStress, MaxStress]
            postvalues[f"{filter}I"]["statsVonMises"] = [
                MinVonMises,
                MeanVonMises,
                MaxVonMises,
            ]
            postvalues[f"{filter}I"]["statsDisplH"] = [
                MinDisplH,
                # MeanDisplH,
                MaxDisplH,
            ]
            postvalues[f"{filter}I"]["statsStressH"] = [
                MinStressH,
                MeanStressH,
                MaxStressH,
            ]
            postvalues[f"{filter}I"]["statsVonMisesH"] = [
                MinVonMisesH,
                MeanVonMisesH,
                MaxVonMisesH,
            ]

    return (targets, postvalues)


def exportResults(
    args,
    parameters: dict,
    table,
    table_final,
    dict_df: dict,
    global_df=None,
    suffix: str = None,
):
    # Sum L*I², product I for Mutual
    sumLI2 = 0
    productI = 1
    for target, values in dict_df.items():
        mname = target[:-2]
        prefix = ""
        if mname:
            prefix = f"{mname}_"
        table_final[f"{prefix}I[A]"] = dict_df[target]["target"]
        table_final[f"{prefix}flow[l/s]"] = dict_df[target]["flow"] * 1e3
        table_final[f"{prefix}Tout[K]"] = dict_df[target]["Tout"]
        Tout_site = dict_df[target].get("MSite_Tout")

        sumLI2 += dict_df[target]["L"] * dict_df[target]["target"] ** 2
        productI *= dict_df[target]["target"]

        for key, df in values.items():
            if args.debug:
                print(f"dict key={key}", flush=True)
            if isinstance(df, pd.DataFrame):
                # df is dataframe -> set new index
                df["I"] = f'I={dict_df[target]["target"]}A'
                df.set_index("I", inplace=True)
                df_T = df.T
                if args.debug:
                    print(tabulate(df_T, headers="keys"), flush=True)
                df_T = df_T.reindex(index=natsorted(df_T.index))
                if global_df:  # if commissioning, store in global_df
                    global_df[mname][key] = pd.concat([global_df[mname][key], df])
                else:  # if cli, df to csv
                    outdir = f"{prefix}{key}.measures"
                    os.makedirs(outdir, exist_ok=True)
                    df_T.to_csv(f"{outdir}/values.csv", index=True)

                if key == "PowerM":
                    table_final[f"{prefix}PowerM[MW]"] = df_T.iloc[0, 0] * 1e-6
                elif key == "PowerH":
                    dfUcoil = df / dict_df[target]["target"]
                    for columnName, columnData in dfUcoil.items():
                        if "H" in columnName:  
                            # if helix, calculate Ucoil 2 helices by 2 
                            nH = int(columnName.split("H", 1)[1])

                            Uname = f"{prefix}Ucoil_H{nH-1}H{nH}[V]"
                            if nH % 2:
                                Uname = f"{prefix}Ucoil_H{nH}H{nH+1}[V]"

                            if Uname in table_final.columns:
                                table_final[Uname] += columnData.iloc[-1]
                            else:
                                table_final[Uname] = columnData.iloc[-1]

                        else:
                            table_final[f"{prefix}Ucoil_{columnName}[V]"] = (
                                columnData.iloc[-1]
                            )

            elif isinstance(df, dict):
                # df is a dict of df, change index and concat them in a single dataframe
                list_dfT = []
                for _key, _df in df.items():
                    if isinstance(_df, pd.DataFrame):
                        _df["I"] = f'{_key}_I={dict_df[target]["target"]}A'
                        _df.set_index("I", inplace=True)
                        list_dfT.append(_df)
                # list_dfT = [_df for keyT, _df in df.items() if isinstance(dfT, pd.DataFrame)]
                if args.debug:
                    print(f"keys={[keyT for keyT, dfT in df.items()]}", flush=True)
                    for _dft in list_dfT:
                        print(tabulate(_dft, headers="keys"), flush=True)
                dfT = pd.concat(list_dfT, sort=True)
                if args.debug:
                    print(f"dfT.keys={dfT.keys()}", flush=True)
                    print(f"dfT={dfT}", flush=True)
                dfT_T = dfT.T
                dfT_T = dfT_T.reindex(index=natsorted(dfT_T.index))
                if args.debug:
                    print(f"dfT_T.keys={dfT.keys()}", flush=True)
                    print(f"dfT_T={dfT_T}", flush=True)
                if global_df: # if commissioning, store in global_df
                    global_df[mname][key] = pd.concat([global_df[mname][key], dfT])
                else: # if cli, df to csv
                    outdir = f"{prefix}{key}.measures"
                    os.makedirs(outdir, exist_ok=True)
                    dfT_T.to_csv(f"{outdir}/values.csv", index=True)

                if key.endswith("H"):
                    # if TH, DisplH, StressH, VonMisesH: go also in table_final
                    symbol = key.replace("stats", "")
                    T_method = {
                        "Min": min,
                        "Max": max,
                    }
                    T_unit = {
                        "statsTH": "K",
                        "statsDisplH": "m",
                        "statsStressH": "Pa",
                        "statsVonMisesH": "Pa",
                    }
                    for columnName, columnData in dfT.items():
                        if args.debug:
                            print(f"columnName={columnName}", flush=True)
                            print(f"columnData={columnData.keys()}", flush=True)
                        for T in ["Min", "Max"]:
                            if "H" in columnName:
                                # if helices, calculate measures 2 helices by 2 
                                nH = int(columnName.split("H", 1)[1])

                                Tname = (
                                    f"{prefix}{T}{symbol}_H{nH-1}H{nH}[{T_unit[key]}]"
                                )
                                if nH % 2:
                                    Tname = f"{prefix}{T}{symbol}_H{nH}H{nH+1}[{T_unit[key]}]"

                                if Tname in table_final.columns:
                                    table_final[Tname] = T_method[T](
                                        table_final[Tname].iloc[-1],
                                        dfT.loc[
                                            f'{T}{symbol}_I={dict_df[target]["target"]}A'
                                        ][columnName],
                                    )
                                else:
                                    table_final[Tname] = dfT.loc[
                                        f'{T}{symbol}_I={dict_df[target]["target"]}A'
                                    ][columnName]

                            elif not re.search(r"_?R\d+", columnName):
                                table_final[
                                    f"{prefix}{T}{symbol}_{columnName}[{T_unit[key]}]"
                                ] = dfT.loc[
                                    f'{T}{symbol}_I={dict_df[target]["target"]}A'
                                ][
                                    columnName
                                ]
                        if symbol != "DisplH": #DisplH doesn't have mean <- fix !
                            if "H" in columnName:
                                # if helices, calculate mean 2 helices by 2 with area values 
                                nH = int(columnName.split("H", 1)[1])

                                Tname = (
                                    f"{prefix}Mean{symbol}_H{nH-1}H{nH}[{T_unit[key]}]"
                                )
                                if nH % 2:
                                    Tname = f"{prefix}Mean{symbol}_H{nH}H{nH+1}[{T_unit[key]}]"
                                    Area = (
                                        parameters[f"Area_{prefix}H{nH}"]
                                        + parameters[f"Area_{prefix}H{nH+1}"]
                                    )
                                else:
                                    Area = (
                                        parameters[f"Area_{prefix}H{nH-1}"]
                                        + parameters[f"Area_{prefix}H{nH}"]
                                    )

                                if Tname in table_final.columns:
                                    table_final[Tname] = (
                                        table_final[Tname].iloc[-1]
                                        + dfT.loc[
                                            f'Mean{symbol}_I={dict_df[target]["target"]}A'
                                        ][columnName]
                                        * parameters[f"Area_{prefix}H{nH}"]
                                    ) / Area
                                else:
                                    table_final[Tname] = (
                                        dfT.loc[
                                            f'Mean{symbol}_I={dict_df[target]["target"]}A'
                                        ][columnName]
                                        * parameters[f"Area_{prefix}H{nH}"]
                                    )

                            elif not re.search(r"_?R\d+", columnName):
                                table_final[
                                    f"{prefix}Mean{symbol}_{columnName}[{T_unit[key]}]"
                                ] = dfT.loc[
                                    f'Mean{symbol}_I={dict_df[target]["target"]}A'
                                ][
                                    columnName
                                ]

        for columnName, columnData in table_final.items():
            # Add R=Ucoil/I to table_final 
            if columnName.startswith(f"{prefix}Ucoil"):
                table_final[
                    columnName.replace("Ucoil", "R").replace("[V]", "[ohm]")
                ] = (columnData / dict_df[target]["target"])

        if dict_df[target]["L"] != 0:
            table_final[f"{prefix}L[H]"] = dict_df[target]["L"]

    if not global_df:
        outdir = f"U.measures"
        os.makedirs(outdir, exist_ok=True)
        table.to_csv(f"{outdir}/values.csv", index=False)

    if Tout_site:
        table_final["MSite_Tout[K]"] = Tout_site
    if "mag" in args.cfgfile:
        df = pd.read_csv("magnetic.measures/values.csv")
        table_final["B0[T]"] = df["Points_B0_expr_Bz"].iloc[-1]
        if len(dict_df) == 1:
            # if only one magnet calculate inductance
            table_final["L[H]"] = (
                2
                * df["Statistics_MagneticEnergy_integrate"].iloc[-1]
                / (dict_df[target]["target"] ** 2)
            )
        elif sumLI2 != 0:
            # if several magnets, calculate mutual
            table_final["M[H]"] = (
                df["Statistics_MagneticEnergy_integrate"].iloc[-1] - sumLI2 / 2
            ) / productI

    table_final.set_index("measures", inplace=True)
    if suffix is None:
        table_final.T.to_csv(f"measures.csv", index=True)
    else:
        table_final.T.to_csv(f"measures{suffix}.csv", index=True)

    print(table_final.T)
    return table_final, global_df


def main():
    fname = "cli"
    description = "Cfpdes model"
    epilog = (
        "Setup for Magnet or Site simulation\n"
        "Workflow: actually fix current and compute cooling BCs using selected heatcorrelation\n"
        "\n"
        "Before running you need a flow_params for each magnet\n"
    )
    # print(f"epilog: {epilog} (type={type(epilog)})", flush=True)

    parser = options(description, epilog)
    args = parser.parse_args()
    pwd = os.getcwd()

    # Load units:
    # TODO force millimeter when args.method == "HDG"
    # units = load_units('meter')

    # Load cfg as config
    dim = 0
    jsonmodel = ""
    meshmodel = ""
    feelpp_config = configparser.ConfigParser()
    basedir = None
    with open(args.cfgfile, "r") as inputcfg:
        feelpp_config.read_string("[DEFAULT]\n[main]\n" + inputcfg.read())
        if "case" in feelpp_config:
            dim = int(feelpp_config["case"]["dimension"])
        else:
            dim = int(feelpp_config["main"]["case.dimension"])
        feelpp_directory = feelpp_config["main"]["directory"]

        basedir = os.path.dirname(args.cfgfile)
        if not basedir:
            basedir = "."

        jsonmodel = feelpp_config["cfpdes"]["filename"]
        jsonmodel = jsonmodel.replace(r"$cfgdir/", f"{basedir}/")

        meshmodel = feelpp_config["cfpdes"]["mesh.filename"]
        meshmodel = meshmodel.replace(r"$cfgdir/", f"{basedir}/")

    # Get Parameters from JSON model file
    parameters = {}
    with open(jsonmodel, "r") as jsonfile:
        dict_json = json.loads(jsonfile.read())
        parameters = dict_json["Parameters"]

    e = None
    (e, f, fields) = init(
        fname,
        e,
        args,
        pwd,
        jsonmodel,
        meshmodel,
        directory=feelpp_directory,
        dimension=dim,
    )

    if e.isMasterRank():
        print(args)
        print(f"cwd: {pwd}")
        print(f"feelpp_directory={feelpp_directory}")
        print(f"dim={dim}")
        print(f"basedir={basedir}")
        print(f"jsonmodel={jsonmodel}")
        print(f"meshmodel={meshmodel}")

    targets = {}
    postvalues = {}

    # args.mdata =
    #  currents:  {magnet.name: {'value': current.value, 'type': magnet.type, 'filter': '', 'flow_params': args.flow_params}}
    if args.mdata:
        targets, postvalues = loadMdata(e, pwd, args, targets, postvalues)

    """
    postvalues = {
        "T": ['MinTH', 'MeanTH', 'MaxTH', 'PowerH'],
        "Stress": ['MinVonMises', 'MinHoop', 'MeanVonMises', 'MeanHoop', 'MaxVonMises', 'MaxHoop']
    }
    # pfields depends from method
    pfields = ['I', 'PowerH', 'MinTH', 'MeanTH', 'MaxTH', 'Flux',
                'MinVonMises', 'MeanVonMises', 'MaxVonMises',
                'MinHoop', 'MeanHoop', 'MaxHoop']
    
    """

    (table, dict_df, e) = oneconfig(
        fname,
        e,
        f,
        fields,
        feelpp_directory,
        f"{pwd}/{jsonmodel}",
        f"{pwd}/{meshmodel}",
        args,
        targets,
        postvalues,
        parameters,
    )
    if args.debug:
        print(f"oneconfig done, rank={e.worldCommPtr().localRank()}", flush=True)

    if e.isMasterRank():
        print("Export result", flush=True)
        table_final = pd.DataFrame(["values"], columns=["measures"])
        table_final, global_df = exportResults(
            args, parameters, table, table_final, dict_df
        )

    if args.debug:
        print(f"end of cli, rank={e.worldCommPtr().localRank()}", flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
