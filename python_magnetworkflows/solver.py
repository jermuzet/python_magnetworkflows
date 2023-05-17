from typing import List

import re
import os
import shutil

import json
from tabulate import tabulate

import feelpp as fpp
import feelpp.toolboxes.core as tb
import feelpp.toolboxes.cfpdes as cfpdes

import pandas as pd

from .params import getTarget, getparam
from .waterflow import waterflow as w
from .real_methods import getDT, getHeatCoeff

def create_field(feel_pb, targets: dict, parameters: dict, wd: str, debug: bool = False):
    """
    create and save a field in h5"""
    print('create_field')

    # TODO: pass space and order in method params
    Xh = fpp.functionSpace(space="Pch", mesh=feel_pb.mesh(), order=1)
    usave = Xh.element()

    for target, values in targets.items():
        params = []
        for p in values['control_params']:
            tmp = getparam(p[0], parameters, p[1], debug)
            params += tmp

        for param in params:
            marker = param.replace('U_','') # 'U_': (values['control_params'][0])[0]
            # print(f'marker: {marker}, goal={values["goals"][marker].iloc[-1]} xx')
            value = float(parameters[param])
            usave.on(range=fpp.markedelements(feel_pb.mesh(), marker), expr=fpp.expr(str(value)))

    usave.save(wd, name='U')
    print('create_field: done')

def update(
    jsonmodel: str,
    parameters: dict,
):
    """
    Update jsonmodel with parameters
    """
    print(f'update {jsonmodel}')

    dict_json = {}
    with open(jsonmodel, "r") as jsonfile:
        dict_json = json.loads(jsonfile.read())
        dict_json["Parameters"] = parameters

    with open(jsonmodel, "w+") as jsonfile:
        jsonfile.write(json.dumps(dict_json, indent=4))

    return 0


# TODO create toolboxes_options on the fly
def init(args, directory: str = ""):
    print(f"init: pwd={os.getcwd()}")
    e = fpp.Environment(
        ["cli.py"],
        opts=tb.toolboxes_options("coefficient-form-pdes", "cfpdes"),
        config=fpp.localRepository(directory),
    )
    e.setConfigFile(args.cfgfile)

    if e.isMasterRank():
        print("Create cfpdes")
    f = cfpdes.cfpdes(dim=2)

    if e.isMasterRank():
        print("Init problem")
    f.init()
    # f.printAndSaveInfo()
    
    print(f"init: done")
    return (e, f)


def solve(feelpp_directory, jsonmodel, args, targets: dict, params: dict, parameters: dict):
    """

    targets: dict of target
    params: dict(target, params:list of parameters name)
    """

    # if args.debug:
    print(f"solve: jsonmodel={jsonmodel}, args={args}")
    
    # suffix of tmp files
    post = ""
    for target, values in targets.items():
        post += f'{target}={values["objectif"]}{values["unit"]}-'

    basedir = os.path.dirname(jsonmodel) # get absolute path instead??
    print(f"solve: jsonmodel={jsonmodel}, basedir={basedir}")
    # save original U.h5
    save_h5 = f'{basedir}/U.h5.init'
    if os.path.isfile(save_h5):
        RuntimeError(f'solve: backup U.h5 to {save_h5} - fails since file already exists')
    else:
        # Rename the file
        shutil.copy2(f'{basedir}/U.h5', save_h5)

    # save original jsonmodel
    save_json = f'{jsonmodel}.init'
    if os.path.isfile(save_json):
        RuntimeError(f'solve: backup jsonmodel to {save_json} - fails since file already exists')
    else:
        # Rename the file
        shutil.copy2(jsonmodel, save_json)

    it = 0
    err_max = 10 * args.eps

    table = []
    headers = ["it"]
    for target in targets:
        for param in params[target]:
            headers.append(f"{param}")  # how to add unit??
        headers.append(f"err_max[{target}]")
    headers.append(f"err_max")

    # Xh = fpp.functionSpace(space="Pch", mesh=f.mesh(), order=1)
    # usave = Xh.element()
    output_df = {}

    (e, f) = init(args, directory=feelpp_directory)
    while it < args.itermax:

        print(f'make a copy of files for it={it}')
        new_json = jsonmodel.replace('.json',f'-it{it}-{post[:-1]}.json')
        print(f'jsonmodel={jsonmodel}, new_json={new_json}')
        # shutil.copy2(jsonmodel, new_json)
        dict_json = {}
        with open(jsonmodel, "r") as jsonfile:
            dict_json = json.loads(jsonfile.read())
            dict_json["Meshes"]["cfpdes"]["Fields"]["U"]["filename"] = f'$cfgdir/U-it{it}.h5'
            # print(f'dict_json={dict_json}')
        with open(new_json, "w+") as jsonfile:
            jsonfile.write(json.dumps(dict_json, indent=4))

        shutil.copy2(f'{basedir}/U.h5', f'{basedir}/U-it{it}.h5')

        print(f'start calc for it={it}')
        if args.debug and e.isMasterRank():
            print("Parameters:", f.modelProperties().parameters())

        # Solve and export the simulation
        try:
            f.solve()
            f.exportResults()
        except:
            raise RuntimeError(
                "cfpdes solver or exportResults fails - check feelpp logs for more info"
            )

        # TODO: get csv to look for depends on cfpdes model used
        table_ = [it]
        err_max = 0
        for target, values in targets.items():
            objectif = -values["objectif"]  
            # multiply by -1 because of orientation of pseudo Axi domain Oy == -U_theta
            filtered_df = getTarget(targets, target, e, args.debug)

            # TODO: add stats for filtered_df to table_: mean, ecart type, min/max??

            error = filtered_df.div(-objectif).add(1)
            err_max_target = max(error.abs().max(axis=1))
            table_.append(err_max_target)
            err_max = max(err_max_target, err_max)

            print(f'{target}: objectif={objectif}')
            print(f'{target}: filtered_df={filtered_df}')
            print(f'{target}: error={error}')
            print(f'{target}: it={it}, err_max={err_max_target}, eps={args.eps}, itmax={args.itermax}')
            
            for param in params[target]:
                marker = param.replace(
                    "U_", ""
                )  # get name from values['control_params'] / change control_params to a list of dict?
                val = filtered_df[marker].iloc[-1]
                ovalue = float(parameters[param])
                table_.append(ovalue)
                nvalue = ovalue * objectif / val
                print(f'{it}: {marker}, goal={objectif}, val={val}, err={error[marker].iloc[-1]}, ovalue={ovalue}, nvalue={nvalue}')
                f.addParameterInModelProperties(param, nvalue)
                parameters[param] = nvalue

                # usave.on(
                #     range=fpp.markedelements(f.mesh(), marker),
                #     expr=fpp.expr(str(nvalue)),
                # )

            # update bcs
            print(f"{target}: it={it}, update BCs")
            p_params = {}
            # TODO upload p_df into a dict like {name: p_df} with name = target (aka key of targets)
            # this way we can get p_df per target as an output for solve
            # NB: pd_df[pname] is a pandas dataframe (pname is parameter name, eg Flux)
            p_df = {}

            for param in values["computed_params"]:
                name = param["name"]
                if e.isMasterRank():
                    print(f"{target}: computed_params {name}")

                if "csv" in param:
                    p_df[name] = getTarget({f"{name}": param}, name, e, args.debug)
                    if args.debug and e.isMasterRank():
                        print(f"{target}: {name}={p_df[name]}")
                else:
                    if args.debug and e.isMasterRank():
                        print(f"{target}: {name} computed")
                    for p in param["params"]:
                        pname = p[0]
                        # if args.debug and e.isMasterRank():
                        print(f"{name}: extract params for {p[0]} (len(p)={len(p)})")
                        tmp = getparam(p[0], parameters, p[1], args.debug)
                        if len(p) == 4:
                            regex_match = re.compile(p[2])
                            if p[3]:
                                tmp = [t for t in tmp if regex_match.fullmatch(t)]
                            else:
                                tmp = [t for t in tmp if not regex_match.fullmatch(t)]
                        tmp.sort()
                        # print(f'{name}: sorted tmp={tmp}')
                        if pname in p_params:
                            p_params[pname] += tmp
                        else:
                            p_params[pname] = tmp

            output_df[target] = p_df
            if args.debug and e.isMasterRank():
                print(f"p_df: {p_df.keys()}")
                print(f"p_params: {p_params.keys()}")
                print(f'p_params[Tw]={p_params["Tw"]}')

            Power = p_df["Power"].iloc[-1, 0]
            SPower_H = p_df["PowerH"].iloc[-1].sum()
            SFlux_H = p_df["Flux"].iloc[-1].sum()
            if e.isMasterRank():
                print(f'{target}: PowerH {p_df["PowerH"]}')
                print(
                    f'{target}: it={it} Power={Power} SPower_H={SPower_H} SFlux_H={SFlux_H} PowerH={p_df["PowerH"].iloc[-1]}'
                )

                if args.debug:
                    for key in p_params:
                        print(f"{target}: {key}={p_params[key]}")

            flow = values["waterflow"]
            Pressure = flow.pressure(abs(objectif))

            Dh = [float(parameters[p]) for p in p_params["Dh"]]
            Sh = [float(parameters[p]) for p in p_params["Sh"]]
            if args.debug and e.isMasterRank():
                print(f'Dh: {p_params["Dh"]}')
                print(f'hwH: {p_params["h"]}')

            Tw = [float(parameters[p]) for p in p_params["Tw"]]
            dTw = [float(parameters[p]) for p in p_params["dTw"]]
            hw = [float(parameters[p]) for p in p_params["hw"]]

            TwH = [float(parameters[p]) for p in p_params["TwH"]]
            dTwH = [float(parameters[p]) for p in p_params["dTwH"]]
            hwH = [float(parameters[p]) for p in p_params["h"]]

            # TODO verify if data are consistant??
            if e.isMasterRank():
                print(f'{target}: len(Dh)={len(Dh)}, len(TwH)={len(TwH)}, len(dTwH)={len(dTwH)}, len(Channels/Slits)={len(hwH)}')
            Umean = flow.umean(abs(objectif), sum(Sh))  # math.fsum(Sh)
            if e.isMasterRank():
                print(f"{target}: it={it}, objectif={abs(objectif)}, Umean={Umean}, Flow={flow}")

            # global:  what to do when len(Tw) != 1
            for i, T in enumerate(Tw):
                dTg = getDT(abs(objectif), flow, Power, Tw[i], Pressure)
                hg = getHeatCoeff(flow, sum(Dh) / float(len(Dh)), Umean, Tw[i])
                f.addParameterInModelProperties(p_params["dTw"][i], dTg)
                f.addParameterInModelProperties(p_params["hw"][i], hg)
                parameters[p_params["hw"][i]] = hg
                parameters[p_params["dTw"][i]] = dTg

            if e.isMasterRank():
                print(f'{target}: Tw={Tw[0]}, param={p_params["dTw"][0]}, umean={Umean}, Power={Power}, dTg={dTg}')

            # per Channel/Slit
            if args.debug and e.isMasterRank():
                print(f'{target} Flux: {p_df["Flux"]}')
            for i, (d, s) in enumerate(zip(Dh, Sh)):
                cname = p_params["Dh"][i].replace("_Dh", "")
                PowerCh = p_df["Flux"].iloc[-1, i]
                dTwi = getDT(abs(objectif), flow, PowerCh, TwH[i], Pressure)
                hi = getHeatCoeff(flow, d, Umean, TwH[i])
                f.addParameterInModelProperties(p_params["dTwH"][i], dTwi)
                f.addParameterInModelProperties(p_params["h"][i], hi)
                parameters[p_params["h"][i]] = hi
                parameters[p_params["dTwH"][i]] = dTwi
                if e.isMasterRank():
                    print(
                        f'{target} Channel{i}: cname={cname}, umean={Umean}, Dh={d}, Sh={s}, Power={PowerCh}, Tw={TwH[i]}, param={p_params["dTwH"][i]}, dTwi={dTwi}'
                    )

            # TODO: how to transform dTg, hg et DTwi, hi en dataframe??

        # update Parameters
        f.updateParameterValues()

        if e.isMasterRank():
            print(f"it={it}, err_max={err_max}, eps={args.eps}, itmax={args.itermax}")

        table_.append(err_max)
        table.append(table_)

        create_field(f, targets, parameters, basedir, args.debug)
        update(jsonmodel, parameters)

 
        if err_max <= args.eps:
            break

        # reload feelpp
        e.setConfigFile(args.cfgfile)
        f = cfpdes.cfpdes(dim=2)
        f.init()

        it += 1

    # Save table (need headers)
    if e.isMasterRank():
        print(tabulate(table, headers, tablefmt="simple"))

    if err_max > args.eps or it >= args.itermax:
        raise RuntimeError(f'Fail to solve {jsonmodel}: err_max={err_max}, it={it}')
    """
    results = {}
    legend = ""
    for target in targets:
        (name, val, paramsdict, params, bcs_params, targets, pfields, postvalues, pvalues) = objectif
        results[name] = bcparams
        legend = legend + f'{name}{val}A' 

    resfile = args.cfgfile.replace('.cfg', f'-{legend}.csv')
    if e.isMasterRank(): 
        print(f"Export result to csv: {os.getcwd()}/{resfile}")
    with open(resfile,"w+") as file:
        df = pd.DataFrame(table, columns = headers)
        df.to_csv(resfile, encoding='utf-8')
    """

    return output_df
