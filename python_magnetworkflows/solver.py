from typing import List

import sys
import os

import feelpp as fpp
import feelpp.toolboxes.core as tb
import feelpp.toolboxes.cfpdes as cfpdes

import pandas as pd

from .params import targetdefs, getTarget, getTargetUnit
from .real_methods import pressure, umean, flow, montgomery

# TODO create toolboxes_options on the fly
def init(args, directory: str=""):
    e = fpp.Environment(sys.argv, opts=tb.toolboxes_options("coefficient-form-pdes","cfpdes"), config=fpp.localRepository(directory))
    e.setConfigFile(args.cfgfile)

    if e.isMasterRank(): print("Create cfpdes")
    f = cfpdes.cfpdes(dim=2)

    if e.isMasterRank(): print("Init problem")
    f.init()
    # f.printAndSaveInfo()

    return (e, f)

def solve(e, f, args, objectifs: list):
    if e.isMasterRank(): 
        print(f"solve: workingdir={ os.getcwd() }")
    
    it = 0
    err_max = 10 * args.eps

    table = []
    headers = ['it']
    for objectif in objectifs:
        (name, val, paramsdict, params, bc_params, targets, pfields, postvalues, pvalues) = objectif
        for key in paramsdict:
            for p in params:
                headers.append(f'{p}_{key}')
        headers.append(f'err_max[{name}]')

    bcparams = {}
    while err_max > args.eps and it < args.itermax :
        
        # Update new value of U_Hi_Cuj on feelpp's senvironment
        for objectif in objectifs:
            (name, val, paramsdict, params, bc_params, targets, pfields, postvalues, pvalues) = objectif
            for key in paramsdict:
                for p in params:
                    entry = f'{p}_{key}'
                    val = float(paramsdict[key][p])
                    # print(f"{entry}: {val}")
                    f.addParameterInModelProperties(entry, val)
        f.updateParameterValues()

        if args.debug and e.isMasterRank():
            print("Parameters after change :", f.modelProperties().parameters())

        # Solve and export the simulation
        # try:
        f.solve()
        f.exportResults()
        # except:
        #    raise RuntimeError("cfpdes solver or exportResults fails - check feelpp logs for more info")

        # TODO: get csv to look for depends on cfpdes model used
        table_ = [it]
        for objectif in objectifs:
            (name, ref, paramsdict, params, bc_params, targets, pfields, postvalues, pvalues) = objectif
            filtered_df = getTarget(name, e, args.debug)

            err_max = 0

            # logging in table
            for key in targets:
                val = targetdefs[name]['value'][0](filtered_df, key)
                target = targets[key]
                """
                if e.isMasterRank():
                    print(f'key={key}, val={val}, target={target},  err={abs(1 - val/target)}')
                """
                err_max = max(abs(1 - val/target), err_max)

                # update val
                for p in targetdefs[name]['control_params']:
                    # print(f'p={p}')
                    table_.append(float(paramsdict[key][p[0]]))
                    paramsdict[key][p[0]] = p[2](paramsdict, key, target, val)
            if e.isMasterRank():
                print(f'{name}: it={it}, current={ref}, err_max={err_max}')
            table_.append(err_max)

        table.append(table_)

                
        # update bcs 
        for objectif in objectifs:
            (name, val, paramsdict, params, bcs_params, targets, pfields, postvalues, pvalues) = objectif
            if e.isMasterRank():
                print(f'{name}: it={it}, current={val}, update BCs')
            
            flux_df = getTarget(pvalues['Flux'], e, args.debug)
            power_df = getTarget(pvalues['PowerH'], e, args.debug)
            SPower_H = power_df.iloc[-1].sum()
            Power = getTarget("Power", e, args.debug)
            if args.debug and e.isMasterRank():
                print(f'it={it} Power={Power.iloc[-1,0]} SPower_H={SPower_H} PowerH={power_df.iloc[-1]}')
            Pressure = pressure(val) # dummy call to get a Pressure value

            Dh = []
            Sh = []
            for p in bcs_params:
                if "Dh" in p: Dh.append(float(bcs_params[p]['Dh']))
                if "Sh" in p: Sh.append(float(bcs_params[p]['Sh']))

            Umean = umean(val, sum(Sh))
            if args.debug and e.isMasterRank():
                print(f'it={it}, val={val}, Umean={Umean}, Flow={flow(val)}')
            for i,(d, s) in enumerate(zip(Dh, Sh)):
                # print(f"new={float(flux_df.iloc[-1, i])} orig={float(flux_df[f'Statistics_Flux_Channel{i}_integrate'].iloc[-1])}")
                PowerCh = flux_df.iloc[-1, i]
                if args.debug and e.isMasterRank():
                    print(f"Channel{i}: umean={Umean}, Dh={d}, Sh={s}, Power={PowerCh}")
                Tw = float(bcs_params[f'Tw{i}']['TwH'])
                dTwi = targetdefs[pvalues['DT']]['value'][0](val, PowerCh, Tw, Pressure)
                hi = targetdefs[pvalues['HeatCoeff']]['value'][0](d, Umean, Tw)
                f.addParameterInModelProperties(f'dTw{i}', dTwi)
                f.addParameterInModelProperties(f'h{i}', hi)        
                if args.debug and e.isMasterRank():
                    print(f'it={it} dTw{i}: {dTwi} hw{i}: {hi}')
                bcparams[f'dTw{i}'] = {'value': dTwi, 'unit': getTargetUnit('DT')}
                bcparams[f'h{i}'] = {'value': hi, 'unit': getTargetUnit('HeatCoeff')}

            Tw = float(bcs_params['Tw']['Tw'])
            dTw = targetdefs[pvalues['DT']]['value'][0](val, SPower_H, Tw, Pressure)
            hw = montgomery(Tw, Umean, sum(Dh)/len(Dh))
            f.addParameterInModelProperties("dTw", dTw)
            f.addParameterInModelProperties("hw", hw )        
            if args.debug and e.isMasterRank():
                print(f'it={it}: dTw={dTw} hw={hw}')
            bcparams['dTw'] = {'value': dTw, 'unit': getTargetUnit('DT')}
            bcparams['hw'] = {'value': hw, 'unit': getTargetUnit('HeatCoeff')}

            f.updateParameterValues()


        it += 1

    # Save table (need headers)
    # print(tabulate(table, headers, tablefmt="simple"))
    results = {}
    legend = ""
    for objectif in objectifs:
        (name, val, paramsdict, params, bcs_params, targets, pfields, postvalues, pvalues) = objectif
        results[name] = bcparams
        legend = legend + f'{name}{val}A' 

    resfile = args.cfgfile.replace('.cfg', f'-{legend}.csv')
    if e.isMasterRank(): 
        print(f"Export result to csv: {os.getcwd()}/{resfile}")
    with open(resfile,"w+") as file:
        df = pd.DataFrame(table, columns = headers)
        df.to_csv(resfile, encoding='utf-8')

    return (results, df)

