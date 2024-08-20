"""
Run feelpp model
"""

import sys
import os
import configparser
import gc

import pandas as pd
import json

from .oneconfig import oneconfig
from .solver import init
from .cli import options, loadMdata, exportResults


def main():
    fname = "commissioning"
    description = "Cfpdes model for Commisionning"
    epilog = (
        "Setup for Magnet or Site Commissioning\n"
        "Workflow: actually fix current and compute cooling BCs using selected heatcorrelation\n"
        "\n"
        "Before running you need to have a flow_params for each magnet\n"
    )

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
        print("commissionning: load cfg", flush=True)

    if args.debug and e.isMasterRank():
        print(args)
        print(f"cwd: {pwd}", flush=True)
        print(f"feelpp_directory={feelpp_directory}", flush=True)
        print(f"dim={dim}", flush=True)
        print(f"basedir={basedir}", flush=True)
        print(f"jsonmodel={jsonmodel}", flush=True)
        print(f"meshmodel={meshmodel}", flush=True)

    Commissioning = True
    commissioning_df = pd.DataFrame()

    global_df = {}
    I = []
    step_i = {}
    for mname, values in args.mdata.items():
        filter = values["filter"]
        global_df[filter[:-1]] = {
            "PowerM": pd.DataFrame(),
            "PowerH": pd.DataFrame(),
            "Flux": pd.DataFrame(),
            "HeatCoeff": pd.DataFrame(),
            "DT": pd.DataFrame(),
            "statsT": pd.DataFrame(),
            "statsTH": pd.DataFrame(),
            "Uw": pd.DataFrame(),
        }
        if "Z" in args.cooling:
            global_df[filter[:-1]]["FluxZ"] = pd.DataFrame()
        if "H" in args.cooling:
            global_df[filter[:-1]]["cf"] = pd.DataFrame()
        if "thmagel" in args.cfgfile:
            global_df[filter[:-1]]["statsDispl"] = pd.DataFrame()
            global_df[filter[:-1]]["statsStress"] = pd.DataFrame()
            global_df[filter[:-1]]["statsVonMises"] = pd.DataFrame()
            global_df[filter[:-1]]["statsDisplH"] = pd.DataFrame()
            global_df[filter[:-1]]["statsStressH"] = pd.DataFrame()
            global_df[filter[:-1]]["statsVonMisesH"] = pd.DataFrame()

        if "steplist" in values:
            step_i[mname] = 0
            args.mdata[mname]["value"] = values["steplist"][step_i[mname]]
            I.append(values["steplist"][step_i[mname]])
        else:
            I.append(values["value"])

    global_df["MSite"] = {"U": pd.DataFrame()}

    targets = {}
    postvalues = {}
    # args.mdata = currents:  {magnet.name: {'value': current.value, 'type': magnet.type, 'filter': '', 'flow_params': args.flow_params}}
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

    nstep = 0
    while Commissioning:
        e.worldComm().barrier()
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

        post = ""
        for value in I:
            post += f"I={value}A-"

        if e.isMasterRank():
            print("Export results", flush=True)
            outdir = f"U_{post[:-1]}.measures"
            os.makedirs(outdir, exist_ok=True)
            table.to_csv(f"{outdir}/values.csv", index=False)

            table = table.iloc[-1:].drop(["it"], axis=1)
            table["I"] = post[:-1]
            table.set_index("I", inplace=True)
            global_df["MSite"]["U"] = pd.concat(
                [global_df["MSite"]["U"], table.iloc[-1:]]
            )
            table_final = pd.DataFrame([f"{I}"], columns=["measures"])
            table_final, global_df = exportResults(
                args,
                parameters,
                table,
                table_final,
                dict_df,
                global_df,
                f"_{post[:-1]}",
            )

            if commissioning_df.empty:
                commissioning_df = table_final.copy()
            else:
                commissioning_df = pd.concat([commissioning_df, table_final])

            del dict_df
            del table
            del table_final

        nstep += 1
        for i, (mname, values) in enumerate(args.mdata.items()):
            filter = values["filter"]
            if "steplist" in values:
                step_i[mname] += 1
                if len(values["steplist"]) == step_i[mname]:
                    Commissioning = False
                else:
                    targets[f"{filter}I"]["objectif"] = values["steplist"][
                        step_i[mname]
                    ]
                I[i] = targets[f"{filter}I"]["objectif"]
            else:
                targets[f"{filter}I"]["objectif"] -= values["step"]
                I[i] = targets[f"{filter}I"]["objectif"]
                if I[i] <= 0 or nstep >= values["stepmax"]:
                    Commissioning = False

        if Commissioning:
            del f
            gc.collect()
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

        e.worldComm().barrier()

    if e.isMasterRank():
        for target, values in global_df.items():
            mname = target
            if mname:
                mname = f"{mname}_"
            for key, df in values.items():
                if isinstance(df, pd.DataFrame):
                    df_T = df.T
                    outdir = f"{mname}{key}.measures"
                    os.makedirs(outdir, exist_ok=True)
                    df_T.to_csv(f"{outdir}/values.csv", index=True)

        commissioning_df.to_csv(f"measures.csv", index=False)

    if args.debug:
        print(f"end of commissioning, rank={e.worldCommPtr().localRank()}", flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
