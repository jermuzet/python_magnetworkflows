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

from .oneconfig import oneconfig

from mpi4py import MPI


def main():
    # print(f'sys.argv({type(sys.argv)})={sys.argv}')

    epilog = (
        "Setup for cfpdes static formulations\n"
        "Support: only magnet of insert type\n"
        "Workflow: actually fix current and compute cooling BCs using Montgomery correlation with an average water velocity\n"
        "\n"
        "Before running adapt flow_params to your magnet setup \n"
    )

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

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
        "--heatcorrelation",
        help="choose cooling model",
        type=str,
        choices=["Montgomery", "Dittus", "Colburn", "Silverberg"],
        default="Montgomery",
    )
    # parser.add_argument(
    #     "--friction",
    #     help="choose friction method",
    #     type=str,
    #     choices=["Blasius", "Filonenko", "Colebrook", "Swanee"],
    #     default="Colebrook",
    # )
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
    if args.debug and rank == 0:
        print(args)

    pwd = os.getcwd()
    if rank == 0:
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
        if rank == 0:
            print(f"feelpp_directory={feelpp_directory}")

        basedir = os.path.dirname(args.cfgfile)
        if rank == 0:
            print(f"basedir={basedir}")
        if not basedir:
            basedir = "."

        jsonmodel = feelpp_config["cfpdes"]["filename"]
        jsonmodel = jsonmodel.replace(r"$cfgdir/", f"{basedir}/")
        if rank == 0:
            print(f"jsonmodel={jsonmodel}")

        meshmodel = feelpp_config["cfpdes"]["mesh.filename"]
        meshmodel = meshmodel.replace(r"$cfgdir/", f"{basedir}/")
        if rank == 0:
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
            if rank == 0:
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
                        ("Dh", f"{filter}Channel\\d+_Dh"),
                        ("Sh", f"{filter}Channel\\d+_Sh"),
                        ("hw", f"{filter}Channel_hw"),
                        ("hwH", f"{filter}Channel\\d+_hw"),
                        ("Zmax", f"{filter}Channel_Zmax"),
                        ("ZmaxH", f"{filter}Channel\\d+_Zmax"),
                    ],
                    "value": (getHeatCoeff),
                    "unit": "W/m2/K",
                }

                DT = {
                    "name": "DT",
                    "params": [
                        ("Tw", f"{filter}Channel_Tw"),
                        ("dTw", f"{filter}Channel_dTw"),
                        ("TwH", f"{filter}Channel\\d+_Tw"),
                        ("dTwH", f"{filter}Channel\\d+_dTw"),
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
                        ("hwH", f"{filter}\\w+hw", "\\w+_Slit\\w+", True),
                        ("Zmax", f"{filter}\\w+Zmax", "\\w+_Slit\\w+", False),
                        ("ZmaxH", f"{filter}\\w+Zmax", "\\w+_Slit\\w+", True),
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
                    "waterflow": waterflow.flow_params(values["flow"]),
                }

            postvalues[f"{filter}I"] = {
                "statsT": [MinT, MeanT, MaxT],
                "statsTH": [MinTH, MeanTH, MaxTH],
            }

    if rank == 0:
        print(f"targets: {targets.keys()}")

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

    (table, dict_df) = oneconfig(
        comm,
        feelpp_directory,
        jsonmodel,
        meshmodel,
        args,
        targets,
        postvalues,
        parameters,
    )

    if rank == 0:
        table_final = pd.DataFrame(["values"], columns=["measures"])

        for target, values in dict_df.items():
            mname = target[:-2]
            table_final[f"{mname}_flow[l/s]"] = dict_df[target]["flow"] * 1e3
            table_final[f"{mname}_Tout[K]"] = dict_df[target]["Tout"]
            print("\n")
            for key, df in values.items():
                if key in ["DT", "HeatCoeff"]:
                    outdir = f"{mname}_{key}.measures"
                    if not os.path.exists(outdir):
                        os.mkdir(outdir)
                    df.to_csv(f"{outdir}/values_noT.csv", index=True)
                if isinstance(df, pd.DataFrame):
                    df_T = df.T
                    outdir = f"{mname}_{key}.measures"
                    if not os.path.exists(outdir):
                        os.mkdir(outdir)
                    df_T.to_csv(f"{outdir}/values.csv", index=True)

                    if key == "PowerM":
                        table_final[f"{mname}_PowerM[MW]"] = df_T.iloc[0, 0] * 1e-6
                    elif key == "PowerH":
                        dfUcoil = df / dict_df[target]["target"]
                        for (columnName, columnData) in dfUcoil.iteritems():
                            if "H" in columnName:
                                nH = int(columnName.split("H", 1)[1])

                                Uname = f"{mname}_Ucoil_H{nH-1}H{nH}[V]"
                                if nH % 2:
                                    Uname = f"{mname}_Ucoil_H{nH}H{nH+1}[V]"

                                if Uname in table_final.columns:
                                    table_final[Uname] += columnData.iloc[-1]
                                else:
                                    table_final[Uname] = columnData.iloc[-1]

                            else:
                                table_final[
                                    f"{mname}_Ucoil_{columnName}[V]"
                                ] = columnData.iloc[-1]

                if key in ["statsT", "statsTH"]:
                    list_dfT = [dfT for keyT, dfT in df.items()]
                    dfT = pd.concat(list_dfT, sort=True).T
                    outdir = f"{mname}_{key}.measures"
                    if not os.path.exists(outdir):
                        os.mkdir(outdir)
                    dfT.to_csv(f"{outdir}/values.csv", index=True)

        outdir = f"U.measures"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        table.to_csv(f"{outdir}/values.csv", index=False)

        if "mag" in args.cfgfile:
            df = pd.read_csv("magnetic.measures/values.csv")
            table_final["B0[T]"] = df["Points_B0_expr_Bz"].iloc[-1]

        table_final.set_index("measures", inplace=True)
        table_final.T.to_csv(f"measures.csv", index=True)

        print(table_final.T)

    return 0


if __name__ == "__main__":
    sys.exit(main())
