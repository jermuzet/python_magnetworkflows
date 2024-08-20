from typing import Optional

import re
import os
import shutil
import numpy as np

import json
from tabulate import tabulate

import feelpp.core as fppc
import feelpp.toolboxes as fppt
import feelpp.toolboxes.core as fpptb
import feelpp.toolboxes.cfpdes as cfpdes

import pandas as pd
import gc

from .params import getTarget, getparam
from .waterflow import waterflow as w
from .cooling import rho, Cp, Uw, getDT, getHeatCoeff, getTout

from .error import compute_error


def create_field(
    feelpp_env,
    feel_pb,
    fields: dict,
    targets: dict,
    parameters: dict,
    wd: str,
    debug: bool = False,
):
    """
    create and save a field in h5 from parameters data
    """

    if feelpp_env.isMasterRank():
        print(f"create_field: fields={fields}", flush=True)

    params = []
    for target, values in targets.items():
        for p in values["control_params"]:
            tmp = getparam(p[0], parameters, p[1], debug)
            params += tmp

    # TODO: pass space and order in method params
    for field, data in fields.items():
        if field != "sigma" and field != "alpha":
            if feelpp_env.isMasterRank():
                print(f"fields[{field}]={data}", flush=True)
            Xh = fppc.functionSpace(
                space=data["space"], mesh=feel_pb.mesh(), order=data["order"]
            )
            usave = Xh.element()

            for param in params:
                if param.startswith(f"{field}_"):
                    marker = param.replace(f"{field}_", "")
                    # check if marker exists in mesh, if not skip the next lines
                    value = parameters[param]
                    # print(f"{param}={value}")
                    usave.on(
                        range=fppc.markedelements(feel_pb.mesh(), marker),
                        expr=fppc.expr(str(value)),
                    )

            usave.save(wd, name=field)
            if feelpp_env.isMasterRank():
                print(f"create_field: {field} done (wd={wd})", flush=True)


def init_field(e, jsonmodel: str, meshmodel: str, dimension: int):
    """
    create and save a field in h5"""
    if e.isMasterRank():
        print(f"init_field: jsonmodel={jsonmodel}, meshmodel={meshmodel}", flush=True)

    basedir = os.path.dirname(jsonmodel)
    with open(jsonmodel, "r") as file:
        data = json.load(file)

    if dimension == 3:
        m2d = fppc.load(
            fppc.mesh(dim=dimension, realdim=dimension), name=meshmodel, verbose=1
        )
    else:
        m2d = fppc.load(fppc.mesh(dim=dimension), name=meshmodel, verbose=1)

    fnames = {}
    # cfpdes: from data["Models"]
    if "Models" in data:
        for model in data["Models"]:
            if "Meshes" in data and model in data["Meshes"]:
                for field, value in data["Meshes"]["cfpdes"]["Fields"].items():
                    basis = value["basis"]
                    # split basis into Space + order:
                    split = re.split(r"(\d+)", basis)
                    Xh = fppc.functionSpace(
                        space=split[0], mesh=m2d, order=int(split[1])
                    )
                    usave = Xh.element()

                    for param in data["Parameters"]:
                        if param.startswith(f"{field}_"):
                            marker = param.replace(f"{field}_", "")
                            value = data["Parameters"][param]
                            # print(f"{param}={value}")
                            usave.on(
                                range=fppc.markedelements(m2d, marker),
                                expr=fppc.expr(str(value)),
                            )

                    usave.save(basedir, name=field)
                    fnames[field] = {
                        "model": model,
                        "space": split[0],
                        "order": int(split[1]),
                    }

    if e.isMasterRank():
        print(f"init_field: {fnames} done (basedir={basedir})", flush=True)

    return fnames


def update(e, jsonmodel: str, parameters: dict):
    """
    Update jsonmodel with parameters
    """

    dict_json = {}
    with open(jsonmodel, "r") as jsonfile:
        dict_json = json.loads(jsonfile.read())
        dict_json["Parameters"] = parameters

    if e.isMasterRank():
        with open(jsonmodel, "w+") as jsonfile:
            jsonfile.write(json.dumps(dict_json, indent=4))
            print(f"update {jsonmodel}", flush=True)
    e.worldComm().barrier()

    return 0


# TODO create toolboxes_options on the fly
def init(
    fname,
    e,
    args,
    pwd: str,
    jsonmodel: str,
    meshmodel: str,
    directory: Optional[str] = None,
    dimension: Optional[int] = 2,
):
    """
    init feelp env and feelpp problem
    """

    # pwd = os.getcwd()
    if not e:
        e = fppc.Environment(
            [f"{fname}.py"],
            opts=fpptb.toolboxes_options("coefficient-form-pdes", "cfpdes"),
            config=fppc.localRepository(directory),
        )
        if e.isMasterRank():
            print(
                f"Init: feelpp env created (pwd={pwd}, cwd={os.getcwd()})", flush=True
            )

    fields = init_field(e, f"{pwd}/{jsonmodel}", f"{pwd}/{meshmodel}", dimension)
    if e.isMasterRank():
        print(f"Init: fields={fields}", flush=True)

    e.setConfigFile(args.cfgfile)
    if e.isMasterRank():
        print(f"Init: setconfigfile", flush=True)

    f = cfpdes.cfpdes(dim=dimension)
    if e.isMasterRank():
        print("Create cfpdes done", flush=True)

    f.init()
    if e.isMasterRank():
        print("Init problem done", flush=True)
    # f.printAndSaveInfo()
    return (e, f, fields)


def solve(
    fname,
    e,
    f,
    fields,
    feelpp_directory: str,
    jsonmodel: str,
    meshmodel: str,
    args,
    targets: dict,
    postvalues: dict,
    params: dict,
    parameters: dict,
    dict_df: dict,
):
    """
    targets: dict of target
    params: dict(target, params:list of parameters name)
    """

    # suffix of tmp files
    post = ""
    for target, values in targets.items():
        post += f'{target}={values["objectif"]}{values["unit"]}-'

    basedir = os.path.dirname(jsonmodel)  # get absolute path instead??

    if e.isMasterRank():
        print(f"solve: jsonmodel={jsonmodel}, basedir={basedir}", flush=True)

    # save original jsonmodel
    save_json = f"{jsonmodel}.init"
    isFile = os.path.isfile(save_json)
    e.worldComm().barrier()
    if isFile:
        raise RuntimeError(
            f"solve: backup {jsonmodel} to {save_json} - fails since file already exists"
        )
    else:
        # Rename the file
        if e.isMasterRank():
            shutil.copy2(jsonmodel, save_json)
        e.worldComm().barrier()

    # save original U.h5 - fields
    for field in fields:
        save_h5 = f"{basedir}/{field}.h5.init"
        isFile = os.path.isfile(save_h5)
        e.worldComm().barrier()
        if isFile:
            raise RuntimeError(
                f"solve: backup {field}.h5 to {save_h5} - fails since file already exists"
            )
        else:
            # Rename the file
            if e.isMasterRank():
                shutil.copy2(f"{basedir}/{field}.h5", save_h5)
    e.worldComm().barrier()

    # how to do it properly in mpi ??
    # check on each proc, return an error??

    it = 0
    err_max = 10 * args.eps

    table = []
    headers = ["it"]
    for target in targets:
        for param in params[target]:
            headers.append(f"{param}")  # how to add unit??
        headers.append(f"err_max[{target}]")
    headers.append(f"err_max")

    # Xh = fppc.functionSpace(space="Pch", mesh=f.mesh(), order=1)
    # usave = Xh.element()

    while it < args.itermax:
        new_json = jsonmodel.replace(".json", f"-it{it}-{post[:-1]}.json")
        if e.isMasterRank():
            print(f"make a copy of files for it={it}, new_json={new_json}", flush=True)

        dict_json = {}
        with open(jsonmodel, "r") as jsonfile:
            dict_json = json.loads(jsonfile.read())
            dict_json["Meshes"]["cfpdes"]["Fields"]["U"][
                "filename"
            ] = f"$cfgdir/U-it{it}-{post[:-1]}.h5"
            # print(f'dict_json={dict_json}')
            csvfiles = [
                value["filename"]
                for p, value in dict_json["Parameters"].items()
                if isinstance(value, dict)
            ]
            # print(f'cvsfiles: {csvfiles}')
            for file in csvfiles:
                _file = file.replace("$cfgdir", basedir)
                if e.isMasterRank():
                    shutil.copy2(f"{_file}", f"{_file}-it{it}-{post[:-1]}.csv")
                e.worldComm().barrier()

        if e.isMasterRank():
            with open(new_json, "w+") as jsonfile:
                jsonfile.write(json.dumps(dict_json, indent=4))

            shutil.copy2(f"{basedir}/U.h5", f"{basedir}/U-it{it}-{post[:-1]}.h5")
            print(f"start calc for it={it}", flush=True)
            if "Z" not in args.cooling and args.debug and e.isMasterRank():
                print("Parameters:", f.modelProperties().parameters())
        e.worldComm().barrier()

        # Solve and export the simulation
        try:
            f.solve()
            f.exportResults()
        except:
            raise RuntimeError(
                "cfpdes solver or exportResults fails - check feelpp logs for more info"
            )

        # TODO: get csv to look for depends on cfpdes model used
        from .error import compute_error

        p_params = {}
        res = ()
        e.worldComm().barrier()
        if e.isMasterRank():
            res = compute_error(
                e,
                f,
                basedir,
                it,
                args,
                targets,
                postvalues,
                params,
                parameters,
                dict_df,
            )
        else:
            dict_df = None
            parameters = None
            p_params = None
            res = None

        # now broadcast: (err_max,..,table_), dict_df, parameters, p_params
        comm = e.worldCommPtr()
        dict_df = comm.localComm().bcast(dict_df, root=0)
        parameters = comm.localComm().bcast(parameters, root=0)
        p_params = comm.localComm().bcast(p_params, root=0)
        res = comm.localComm().bcast(res, root=0)

        # extract res
        assert isinstance(res, tuple), f"{res}"

        (err_max, err_max_dT, err_max_h, table_, p_params) = res

        table_.append(err_max)
        table.append(table_)
        if (
            err_max <= args.eps
            and err_max_dT <= max(args.eps, 1e-2)
            and err_max_h <= max(args.eps, 1e-2)
        ):
            break

        # make f.addParameterInModelProperties()
        if e.isMasterRank():
            print("update params", flush=True)

        for target in targets:
            for param in params[target]:
                if args.debug and e.isMasterRank():
                    print(
                        f"f.addParameterInModelProperties({param}, {parameters[param]}), rank={comm.localRank()}"
                    )
                f.addParameterInModelProperties(param, parameters[param])

        selected_params = ["hw", "dTw"]
        if "H" in args.cooling:
            selected_params = ["hwH", "dTwH"]
            if "Z" in args.cooling:
                selected_params = ["hwH"]

        if e.isMasterRank():
            print("update p_params", flush=True)
        for param, values in p_params.items():
            # print(
            #    f"param={param}, values={values} (type={type(values)}), rank={comm.localRank()}"
            # )
            if param in selected_params:
                # print(
                #    f"param={param}, values={values} (type={type(values)}), rank={comm.localRank()}, selected"
                # )
                if isinstance(values, list):
                    for val in values:
                        if args.debug and e.isMasterRank():
                            print(
                                f"f.addParameterInModelProperties({val}, {parameters[val]}, rank={comm.localRank()}, list"
                            )
                        f.addParameterInModelProperties(val, parameters[val])
                else:
                    if args.debug and e.isMasterRank():
                        print(
                            f"f.addParameterInModelProperties({values}, {parameters[values]}, rank={comm.localRank()}"
                        )
                    f.addParameterInModelProperties(values, parameters[values])

        # update Parameters
        f.updateParameterValues()

        if e.isMasterRank():
            print(
                f"it={it}, err_max={err_max:.3e}, err_max_dT={err_max_dT:.3e}, err_max_h={err_max_h:.3e}, eps={args.eps}, itmax={args.itermax}",
                flush=True,
            )

        create_field(e, f, fields, targets, parameters, basedir, args.debug)
        update(e, jsonmodel, parameters)

        # reload feelpp
        if args.reloadcfg:
            e.worldComm().barrier()
            if e.isMasterRank():
                print("reload cfg", flush=True)
            e.setConfigFile(args.cfgfile)
            dimension = f.mesh().dimension()
            f = cfpdes.cfpdes(dim=dimension)
            f.init()

        it += 1

    # Save table (need headers)
    if e.isMasterRank():
        print(tabulate(table, headers, tablefmt="simple"), flush=True)

    table_df = pd.DataFrame(table, columns=headers)

    if err_max > args.eps or it >= args.itermax:
        raise RuntimeError(f"Fail to solve {jsonmodel}: err_max={err_max}, it={it}")

    e.worldComm().barrier()
    if e.isMasterRank():
        print("remove temporary files", flush=True)
        for field in fields:
            save_h5 = f"{basedir}/{field}.h5.init"
            os.remove(save_h5)
        os.remove(save_json)
        for i in range(it):
            print(f"remove {basedir}/U-it{i}-{post[:-1]}.h5", flush=True)
            os.remove(f"{basedir}/U-it{i}-{post[:-1]}.h5")
            tmp_jsonmodel = jsonmodel.replace(".json", f"-it{i}-{post[:-1]}.json")
            print(f"remove {tmp_jsonmodel}", flush=True)
            os.remove(tmp_jsonmodel)
        if "Z" in args.cooling:
            for i in range(it):
                for file in csvfiles:
                    _file = file.replace("$cfgdir", basedir)
                    print(f"remove {_file}-it{i}-{post[:-1]}.csv", flush=True)
                    os.remove(f"{_file}-it{i}-{post[:-1]}.csv")
        print("end of solve")

    return (table_df, dict_df, e)
