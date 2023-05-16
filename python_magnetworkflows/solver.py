from typing import List

import re
import os
import json
from tabulate import tabulate

import feelpp as fpp
import feelpp.toolboxes.core as tb
import feelpp.toolboxes.cfpdes as cfpdes

import pandas as pd

from .params import getTarget, getparam
from .waterflow import waterflow as w
from .real_methods import getDT, getHeatCoeff

def update(
    jsonmodel: str,
    parameters: dict,
    post: str,
    debug: bool = False,
):
    """
    Save a copy of jsonmodel with update parameters
    """

    with open(jsonmodel, "r") as jsonfile:
        dict_json = json.loads(jsonfile.read())
        dict_json["Parameters"] = parameters

    new_name_json = jsonmodel.replace(".json", f"-{post}.json")

    with open(new_name_json, "w+") as jsonfile:
        jsonfile.write(json.dumps(dict_json, indent=4))

    return 0

# TODO create toolboxes_options on the fly
def init(args, directory: str=""):
    print(f'init: pwd={os.getcwd()}')
    e = fpp.Environment(['cli.py'], opts=tb.toolboxes_options("coefficient-form-pdes","cfpdes"), config=fpp.localRepository(directory))
    e.setConfigFile(args.cfgfile)

    if e.isMasterRank(): print("Create cfpdes")
    f = cfpdes.cfpdes(dim=2)

    if e.isMasterRank(): print("Init problem")
    f.init()
    # f.printAndSaveInfo()

    return (e, f)

def solve(e, f, jsonmodel, args, targets: dict, params: dict, parameters: dict):
    """
    
    targets: dict of target
    params: dict(target, params:list of parameters name)
    """

    basedir = os.path.dirname(jsonmodel)

    if e.isMasterRank(): 
        print(f"solve: workingdir={ os.getcwd() }")
        print(f"args: {args}")
    
    it = 0
    err_max = 10 * args.eps

    table = []
    headers = ['it']
    for target in targets:
        for param in params[target]:
            headers.append(f'{param}') # how to add unit??
        headers.append(f'err_max[{target}]')
    headers.append(f'err_max')

    Xh = fpp.functionSpace(space="Pch", mesh=f.mesh(), order=1)
    usave = Xh.element()

    while err_max > args.eps and it < args.itermax :
        
        if args.debug and e.isMasterRank():
            print("Parameters:", f.modelProperties().parameters())

        # Solve and export the simulation
        try:
            f.solve()
            f.exportResults()
        except:
            raise RuntimeError("cfpdes solver or exportResults fails - check feelpp logs for more info")

        # TODO: get csv to look for depends on cfpdes model used
        table_ = [it]
        err_max = 0
        for target, values in targets.items():
            objectif = -values['objectif'] # multiply by -1 because of orientation of pseudo Axi domain Oy == -U_theta
            filtered_df = getTarget(targets, target, e, args.debug)

            # print(f'objectif={objectif}')
            # print(f'filtered_df={filtered_df}')
            # print(f'goals={values["goals"]}')
            for param in params[target]:
                marker = param.replace('U_','') # get name from values['control_params']
                val = filtered_df[marker].iloc[-1]
                # print(f'marker: {marker}, goal={objectif}, val={val}, err={error[marker].iloc[-1]}')
                ovalue = float(parameters[param])
                table_.append(ovalue)
                nvalue = ovalue * val / objectif 
                # print(f'marker: {marker}, ovalue={ovalue}, nvalue={nvalue}')
                f.addParameterInModelProperties(param, nvalue)
                parameters[param] = nvalue

                usave.on(range=fpp.markedelements(f.mesh(), marker), expr=fpp.expr(str(nvalue)))
                
            
            error = filtered_df.div(-objectif).add(1)
            # print(f'error: {error}')
            err_max_target = max(error.abs().max(axis=1))
            table_.append(err_max_target)
            # print(f'err_max[{target}]: {err_max_target}')
            err_max = max(err_max_target, err_max)
        
            # update bcs 
            print(f'{target}: it={it}, update BCs')
            p_params = {}
            p_df = {}
                    
            for param in values['computed_params']:
                name = param['name']
                if e.isMasterRank():
                    print(f'{target}: computed_params {name}')
                
                if "csv" in param:
                    p_df[name] = getTarget({f'{name}': param}, name, e, args.debug)
                    if args.debug and e.isMasterRank():
                        print(f'{target}: {name}={p_df[name]}')
                else:
                    if args.debug and e.isMasterRank():
                        print(f'{target}: {name} computed')
                    for p in param['params']:
                        pname = p[0]
                        # if args.debug and e.isMasterRank():
                        print(f"{name}: extract params for {p[0]}")
                        tmp = getparam(p[0], parameters, p[1], args.debug)
                        if len(p) == 3:
                            regex_match = re.compile(p[2])
                            if p[3]:
                                tmp = [t for t in tmp if regex_match.fullmatch(t)]
                            else:
                                tmp = [t for t in tmp if not regex_match.fullmatch(t)]
                        print(f'{name}: tmp={tmp}')
                        tmp.sort()
                        print(f'{name}: sorted tmp={tmp}')
                        if pname in p_params:
                            p_params[pname] += tmp
                        else:
                            p_params[pname] = tmp
                            
            print(f'p_df: {p_df.keys()}')
            print(f'p_params: {p_params.keys()}')
            print(f'p_params[Tw]={p_params["Tw"]}')

            Power = p_df['Power'].iloc[-1,0]
            SPower_H = p_df['PowerH'].iloc[-1].sum()
            if e.isMasterRank():
                print(f'{target}: it={it} Power={Power} SPower_H={SPower_H} PowerH={p_df["PowerH"].iloc[-1]}')
            
            for key in p_params:
                print(f'{target}: {key}={p_params[key]}') 
                
            flow = values['waterflow']
            Pressure = flow.pressure(abs(objectif))

            Dh = [float(parameters[p]) for p in p_params["Dh"]]
            Sh = [float(parameters[p]) for p in p_params["Sh"]]
            
            Tw = [float(parameters[p]) for p in p_params["Tw"]]
            dTw = [float(parameters[p]) for p in p_params["dTw"]]
            hw = [float(parameters[p]) for p in p_params["hw"]]

            TwH = [float(parameters[p]) for p in p_params["TwH"]]
            dTwH = [float(parameters[p]) for p in p_params["dTwH"]]
            hwH = [float(parameters[p]) for p in p_params["h"]]

            print(f'{target}: len(Dh)={len(Dh)}, len(Channels/Slits)={len(hwH)}')
            Umean = flow.umean(abs(objectif), sum(Sh)) #math.fsum(Sh)
            if args.debug and e.isMasterRank():
                print(f'it={it}, val={val}, Umean={Umean}, Flow={flow(val)}')

            # global
            dTg = getDT(abs(objectif), flow, Power, Tw[0], Pressure)
            hg = getHeatCoeff(flow, sum(Dh) / float(len(Dh)), Umean, Tw[0])
            f.addParameterInModelProperties(p_params["dTw"][0], dTg)
            f.addParameterInModelProperties(p_params["hw"][0], hg )
            print(f'{target}: Tw={Tw[0]}, param={p_params["dTw"][0]}, dTg={dTg}')

            # per Channel/Slit
            print(f'{target} Flux: {p_df["Flux"]}')
            for i,(d, s) in enumerate(zip(Dh, Sh)):
                cname =  p_params["Dh"][i].replace('_Dh','')
                print(f'cname={cname}')
                PowerCh = p_df['Flux'].iloc[-1, i]
                if e.isMasterRank():
                    print(f"Channel{i}: umean={Umean}, Dh={d}, Sh={s}, Power={PowerCh}")
                dTwi = getDT(abs(objectif), flow, PowerCh, TwH[i], Pressure)
                hi = getHeatCoeff(flow, d, Umean, TwH[i])
                f.addParameterInModelProperties(p_params["dTwH"][i], dTwi)
                f.addParameterInModelProperties(p_params["h"][i], hi)
                if e.isMasterRank():
                    print(f'{target} Channel{i}: Tw={TwH[i]}, param={p_params["dTwH"][i]}, dTwi={dTwi}')

        # update Parameters
        f.updateParameterValues()

        if e.isMasterRank():
            print(f'it={it}, err_max={err_max}')

        print(f'err_max: {err_max}')
        table_.append(err_max)
        table.append(table_)
        
        usave.save(basedir, name='U_new')

        # how to save jsonmodel to a tmp file
        post = ""
        for target, values in targets.items():
            post += f'{target}={values["objectif"]}{values["unit"]}-'
        
        update(jsonmodel, parameters, post[:-1], args.debug)

        it += 1
        break

    # Save table (need headers)
    print(tabulate(table, headers, tablefmt="simple"))
    
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
    
    results = {}
    df = pd.DataFrame()
    return (results, df)

