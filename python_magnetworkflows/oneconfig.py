"""
Run feelpp model for one config
"""
from typing import List

import os
from math import fabs 

import re
import pandas as pd
from tabulate import tabulate

from .params import targetdefs, setTarget, getTarget, getTargetUnit, update
from .solver import solve
# from ..units import load_units


def oneconfig(cwd, feelpp_env, feel_pb, jsonmodel, args, objectifs: list):
    """
    Run a simulation until currents are reached
    
    Right now only use vcurrents[0]
    
    # insert: IH, params N_H\* control_params U_H\* 'Statistics_Intensity_H\w+_integrate'
    # bitter: IB, other U_\*, extract name from U_* to get N_*
    # supra: IS params from I_\*, ............. I_\* to get N_*

    objectif: name, value, params, control_params, bc_params
    """

    table_values = []
    table_headers = []

    if feelpp_env.isMasterRank():
        print(f"oneconfig: workingdir={ os.getcwd() }")

    for objectif in objectifs:
        (name, value, params, control_params, bc_params, targets, pfields, postvalues, pvalues) = objectif
        if feelpp_env.isMasterRank():
            print(f"{name}: {value}")
        table_values.append(float(value))
        table_headers.append(f'{name}[{getTargetUnit(name)}]')    # print("targets:", targets)

    (bcparams, _Current_df) = solve(feelpp_env, feel_pb, args, objectifs)

    # update
    results = {}
    for objectif in objectifs:
        dict_df = {}
        
        (name, value, params, control_params, bc_params, targets, pfields, postvalues, pvalues) = objectif
        update(cwd, jsonmodel, params, control_params, bc_params, value, args.debug)
        if args.debug and feelpp_env.isMasterRank():
            print(f'bcparams[{name}]: {bcparams[name]}')
        # print(f'params: {params}')

        df = getTarget(pvalues['Power'], feelpp_env, args.debug)
        table_headers.append('R[Ohm]')
        table_values.append(float(df.iloc[-1][0]/(value*value)))
        if not df.empty and feelpp_env.isMasterRank():
            print(f"I: {value} [A]\tPower: {df.iloc[-1][0]} [W]")
   
        filename = os.path.join(os.getcwd(), 'magnetic.measures/values.csv')
        if os.path.isfile(filename):
            for field in ['Inductance', 'B0']:
                df = getTarget(field, feelpp_env, args.debug)
                table_headers.append(f'{field}[getTargetUnit{field}]')
                table_values.append(float(df.iloc[-1][0]))
                if feelpp_env.isMasterRank():
                    print(f'magnetic: {field}={df.iloc[-1][0]}')
        else:
            if feelpp_env.isMasterRank():
                print(f'{filename}: no such file')

        if feelpp_env.isMasterRank():
            print(f'Save tables to csv files: workingdir={os.getcwd()}')

        # Stats - to be displayed in results
        for p in pfields:
            df = getTarget(p, feelpp_env, args.debug)
            if not df.empty:
                # create new key dict
                keys = df.columns.values.tolist()
                nkeys = {}
                for item in keys:
                    nitem = item
                    if re.match(targetdefs[p]['rematch'], item):
                        regexp = re.split('_', targetdefs[p]['rematch'])
                        nitem = item.replace(regexp[0], '')
                        nitem = nitem.replace(regexp[1], '')
                        nitem = nitem.replace(regexp[3], '')
                        nitem = nitem.replace('__', '')
                        # remove _ if nitem end
                        if nitem.endswith('_'):
                            nitem = nitem[:-1]
                        # print(f"{p}: {nitem}")
                    nkeys[item] = nitem

                # remove all rows except latest row
                df.rename(columns=nkeys, inplace=True)
                if feelpp_env.isMasterRank():
                    print(f"{p} [{targetdefs[p]['unit']}]:\n{tabulate(df, headers='keys', tablefmt='simple')}\n")

                _df = df.drop(range(len(df.index)-1))
                # print(f'_df={_df.transpose()}')
                dict_df[p] = _df

        if args.debug and feelpp_env.isMasterRank():
            print(f'dict_df: {dict_df.keys()}')
            print("params:", params)
            print("bc_params:", bc_params)

        # Compute R per part
        _df = dict_df[pvalues['PowerH']]
        (nrows, ncolumns) = _df.shape
        for i in range(ncolumns):
            _value = float(_df.iloc[-1,i])/(value*value)
            table_headers.append(f"R{i}[Ohm]")
            table_values.append(_value)

        # save cooling params to csv
        for p in pfields:
            if p in dict_df:
                _df = dict_df[p]
                unit = f'{p}[{getTargetUnit(p)}]'
                _df.rename(index={1: unit}, inplace=True)

        _values = {}

        for param in bcparams[name]:
            table_headers.append(f"{param}[{bcparams[name][param]['unit']}]")
            _value = bcparams[name][param]['value']
            table_values.append(_value)
            if re.match('^h\d+', param):
                _name = f"Channel{param.replace('h','')}"
            if re.match('^dTw\d+', param):
                _name = f"Channel{param.replace('dTw','')}"
            _values[_name] = [_value]
            
        _h_df = pd.DataFrame(_values, index=['h[W/m2/K]','Tw[K]','dTw[K]'])
        _Cooling_df = pd.concat([_h_df, dict_df[pvalues['Flux']]])
        if args.debug and feelpp_env.isMasterRank():
            print(f'_Cooling_df: {_Cooling_df}')
        os.makedirs('cooling.measures', exist_ok=True)
        _Cooling_df.to_csv('cooling.measures/postvalues.csv', index=True)

        # Temp stats
        """
        /workspaces/python_magnetworkflows/python_magnetworkflows/cli.py:233: FutureWarning: Sorting because non-concatenation axis is not aligned.
        A future version of pandas will change to not sort by default.

        To accept the future behavior, pass 'sort=False'.

        To retain the current behavior and silence the warning, pass 'sort=True'.

        _T_df = pd.concat([dict_df['MinTH'],dict_df['MeanTH'], dict_df['MaxTH'], dict_df['PowerH']])
        /workspaces/python_magnetworkflows/python_magnetworkflows/cli.py:240: FutureWarning: Sorting because non-concatenation axis is not aligned.
        A future version of pandas will change to not sort by default.

        To accept the future behavior, pass 'sort=False'.

        To retain the current behavior and silence the warning, pass 'sort=True'.
        """

                    
        _T_df = pd.concat([dict_df[p] for p in postvalues['T']], sort=True)
        for p in postvalues['T']:
            _T_df.rename(index={f"{p}[{getTargetUnit(p)}]": f"{p.replace('H','')}[{getTargetUnit(p)}]"}, inplace=True)
        if args.debug and feelpp_env.isMasterRank():
            print(f'_T_df = {_T_df}')
        _T_df.to_csv('heat.measures/postvalues.csv', index=True)

        # Stress, VonMises stats
        if postvalues['Stress'][0] in dict_df:
            _Stress_df = pd.concat([dict_df[p] for p in postValues['Stress']], sort=True)
            for p in postvalues["Stress"]:
                _Stress_df.rename(index={f"{p}[{unit}]": f"{p}[{unit}]"}, inplace=True)
            if feelpp_env.isMasterRank():
                print(f'_Stress_df = {_Stress_df}')
            _Stress_df.to_csv('elastic.measures/postvalues.csv', index=True)

        if feelpp_env.isMasterRank():
            print(f"Commissionning:\n{tabulate([table_values], headers=table_headers, tablefmt='simple')}\n")
    
        results[name] = (table_headers, table_values)

    return results
