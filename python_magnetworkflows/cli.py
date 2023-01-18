"""
Run feelpp model 
"""

import sys
import os
import argparse
import configparser

import re
import json
import pandas as pd
from tabulate import tabulate

from .params import targetdefs, setTarget, getTarget, getTargetUnit, getparam, update, Merge
from .solver import init, solve
from .real_methods import flow_params
# from ..units import load_units

def main():

    epilog = "Setup for cfpdes static formulations\n" \
             "Support: only magnet of insert type\n" \
             "Workflow: actually fix current and compute cooling BCs using Montgomery correlation with an average water velocity\n" \
             "\n" \
             "Before running adapt flow_params to your magnet setup \n"

    command_line = None
    parser = argparse.ArgumentParser(description="Cfpdes HiFiMagnet Fully Coupled model")
    parser.add_argument("cfgfile", help="input cfg file (ex. HL-31.cfg)")
    parser.add_argument("--wd", help="set a working directory", type=str, default="")
    parser.add_argument("--current", help="specify requested current (default: 31kA)", nargs='+', metavar='Current', type=float, default=[31.e+3])
    # TODO: from_current, to_current, step_current to mimic commissionning
    # parser.add_argument("--current", help="specify requested current (default: 31kA)", nargs='+', metavar='Current', type=float, default=[31.e+3])
    parser.add_argument("--cooling", help="choose cooling type", type=str,
                    choices=['mean', 'grad', 'meanH', 'gradH'], default='mean')
    parser.add_argument("--eps", help="specify requested tolerance (default: 1.e-3)", type=float, default=1.e-3)
    parser.add_argument("--itermax", help="specify maximum iteration (default: 10)", type=int, default=10)
    parser.add_argument("--debug", help="activate debug", action='store_true')
    parser.add_argument("--verbose", help="activate verbose", action='store_true')
    parser.add_argument("--flow_params", help="select flow param json file", type=str, default="flow_params.json")
    # TODO add flow max, Vpmax, Imax, Vpump correlation

    args = parser.parse_args()
    if args.debug:
        print(args)

    cwd = os.getcwd()
    if args.wd:
        os.chdir(args.wd)

    # Load units:
    # TODO force millimeter when args.method == "HDG"
    # units = load_units('meter')

    # load flow params
    flow_params(args.flow_params)

    print(f"Current: {args.current}")
    currents = []
    if isinstance(args.current, list):
        currents = args.current
    else:
        currents.append(args.current)

    # Load cfg as config
    feelpp_config = configparser.ConfigParser()
    with open(args.cfgfile, 'r') as inputcfg:
        feelpp_config.read_string('[DEFAULT]\n[main]\n' + inputcfg.read())
        if args.debug:
            print(f"feelpp_cfg: {feelpp_config['main']}")

            for section in feelpp_config.sections():
                print(f"section: {section}")

        jsonmodel = feelpp_config['cfpdes']['filename']
        jsonmodel = jsonmodel.replace(r"$cfgdir/",'')
        if args.debug:
            print(f"jsonmodel={jsonmodel}")

    # Get Parameters from JSON model file
    parameters = {}
    with open(jsonmodel, 'r') as jsonfile:
        dict_json = json.loads(jsonfile.read())
        parameters = dict_json['Parameters']

    params = {}
    bc_params = {}
    control_params = []
    for p in targetdefs['I']['control_params']:
        if args.debug:
            print(f"extract control params for {p[0]}")
        control_params.append(p[0])
        tmp = getparam(p[0], parameters, p[1], args.debug)
        params = Merge(tmp, params, args.debug)

    for p in targetdefs['I']['params']:
        if args.debug: print(f"extract compute params for {p[0]}")
        tmp = getparam(p[0], parameters, p[1], args.debug)
        params = Merge(tmp, params, args.debug)

    # Get Bc params
    for p in [ 'HeatCoeff', 'DT' ]:
        if args.debug:
            print(f"extract bc params for {p}")
        for bc_p in targetdefs[p]['params']:
            if args.debug:
                print(f"{bc_p[0]}")
            tmp = getparam(bc_p[0], parameters, bc_p[1], args.debug)
            bc_params = Merge(tmp, bc_params, args.debug)
        if args.debug:
            print(f"extract bc control_params for {p}")
        for bc_p in targetdefs[p]['control_params']:
            if args.debug:
                print(f"{bc_p[0]}")
            tmp = getparam(bc_p[0], parameters, bc_p[1], args.debug)
            bc_params = Merge(tmp, bc_params, args.debug)

    # print("params:", params)
    if args.debug:
        print("bc_params:", bc_params)

    pfields = ['I', 'PowerH', 'MinTH', 'MeanTH', 'MaxTH', 'Flux']
    # TODO add MeanVonMisesH, MeanHoopH
    pfields.append('MinVonMisesH')
    pfields.append('MeanVonMisesH')
    pfields.append('MaxVonMisesH')
    pfields.append('MinHoopH')
    pfields.append('MeanHoopH')
    pfields.append('MaxHoopH')

    # define targets
    # TODO for insert + bitter [ + supra] ???
    # shall build targetdefs on the fly
    # insert: IH, params N_H\* control_params U_H\* 'Statistics_Intensity_H\w+_integrate'
    # bitter: IB, other U_\*, extract name from U_* to get N_*
    # supra: IS params from I_\*, ............. I_\* to get N_*
    targets = setTarget('I', params, currents[0], args.debug)
    # print("targets:", targets)

    # init feelpp env
    (feelpp_env, feel_pb) = init(args)

    # solve (output params contains both control_params and bc_params values )
    # TODO loop from_current to to_current by step current
    (params, bcparams) = solve(feelpp_env, feel_pb, args, 'I', params, control_params, bc_params, targets)

    # update
    update(cwd, jsonmodel, params, control_params, bcparams, currents[0], args.debug)

    # display csv results
    # TODO use units
    # get Power for Insert, Bitter
    if feelpp_env.isMasterRank():
        print(f"update: workingdir={ os.getcwd() }")
    df = getTarget('Power', feelpp_env, args.debug)
    if feelpp_env.isMasterRank():
        print(f"I: {currents[0]} [A]\tPower: {df.iloc[-1][0]} [W]")
        # TODO add bitter current and power if any Bitters

    # TODO when from_current:
    # create same cvs as below but tables are by current
    # mean that the following line are in the loop
    # and that the save is done afterwards

    if feelpp_env.isMasterRank():
        print(f'Save tables to csv files: workingdir={os.getcwd()}')
    
    # Stats - to be displayed in results
    # TODO: list of fields depends on method used
    dict_df = {}
    for p in pfields:
        df = getTarget(p, feelpp_env, args.debug)
        if not df is None:
            # TODO change keys (symbols+units)
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

            # save only last line
            # TODO remove all rows except latest row
            df.rename(columns=nkeys, inplace=True)
            if feelpp_env.isMasterRank():
                print(f"{p} [{targetdefs[p]['unit']}]:\n{tabulate(df, headers='keys', tablefmt='simple')}\n")
            
            if p in dict_df:
                if feelpp_env.isMasterRank():
                    print('copy last row from df to _df')
            else:
                _df = df.drop(range(len(df.index)-1))
                # print(f'_df={_df.transpose()}')
                dict_df[p] = _df
    
    # TODO
    # save cooling params to csv
    # save transposed version of the table
            
    for p in pfields:
        if p in dict_df:
            _df = dict_df[p]
            unit = f'{p}[{getTargetUnit(p)}]'
            _df.rename(index={1: unit}, inplace=True)
            # _df = _df.transpose()

    if args.debug and feelpp_env.isMasterRank():
        print("params:", params)
    if feelpp_env.isMasterRank():
        print("bc_params:", bc_params)
    _values = {}
    for param in bc_params:
        # print(f'{param}')
        if re.match('^h\d+', param):
            # print(f"bc[{param}] = {bc_params[param]['h']}")
            _h =  bc_params[param]['h']
            _name = f"Channel{param.replace('h','')}"
            _values[_name] = [_h]

    for param in bc_params:
        if re.match('^Tw\d+', param):
            _name = f"Channel{param.replace('Tw','')}"
            _Tw = bc_params[param]['TwH']
            _values[_name].append(_Tw)

    for param in bc_params:
        if re.match('^dTw\d+', param):
            _name = f"Channel{param.replace('dTw','')}"
            _dTw = bc_params[param]['dTwH']
            _values[_name].append(_dTw)

    _h_df = pd.DataFrame(_values, index=['h','Tw','dTw'])
    _Cooling_df = pd.concat([_h_df, dict_df['Flux']])
    # print(f'_Cooling_df: {_Cooling_df}')
    os.makedirs(f'cooling.measures', exist_ok=True)
    _Cooling_df.to_csv(f'cooling.measures/postvalues.csv', index=True)

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

    _T_df = pd.concat([dict_df['MinTH'],dict_df['MeanTH'], dict_df['MaxTH'], dict_df['PowerH']], sort=True)
    for p in ['PowerH', 'MinTH', 'MeanTH', 'MaxTH']:
        _T_df.rename(index={f"{p}[{getTargetUnit(p)}]": f"{p.replace('H','')}"}, inplace=True)
    if args.debug and feelpp_env.isMasterRank():
        print(f'_T_df = {_T_df}')
    _T_df.to_csv(f'heat.measures/postvalues.csv', index=True)

    # Stress, VonMises stats
    if 'MeanVonMisesH' in dict_df:
        _Stress_df = pd.concat([dict_df['MeanVonMisesH'],dict_df['MeanHoopH']], sort=True)
        for p in ['MeanVonMisesH', 'MeanHoopH']:
            _Stress_df.rename(index={f"{p}[{getTargetUnit(p)}]": f"{p.replace('H','')}"}, inplace=True)
        if args.debug and feelpp_env.isMasterRank():
            print(f'_Stress_df = {_Stress_df}')
        _Stress_df.to_csv(f'elastic.measures/postvalues.csv', index=True)

if __name__ == "__main__":
    sys.exit(main())
