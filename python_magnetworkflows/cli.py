import sys
import os
import argparse
import configparser

import re
import json
from tabulate import tabulate

from .params import targetdefs, setTarget, getTarget, getparam, update, Merge
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
        
    print("Current:", args.current)
    currents = []
    if isinstance(args.current, list):
        currents = args.current
    else:
        currents.append(args.current)
        

    # Load cfg as config
    import configparser
    feelpp_config = configparser.ConfigParser()
    with open(args.cfgfile, 'r') as inputcfg:
        feelpp_config.read_string('[DEFAULT]\n[main]\n' + inputcfg.read())
        if args.debug:
            print("feelpp_cfg:", feelpp_config['main'])
    
            for section in feelpp_config.sections():
                print("section:", section)

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
        if args.debug: print(f"extract control params for {p[0]}")
        control_params.append(p[0])
        tmp = getparam(p[0], parameters, p[1], args.debug)
        params = Merge(tmp, params, args.debug)

    for p in targetdefs['I']['params']:
        if args.debug: print(f"extract compute params for {p[0]}")
        tmp = getparam(p[0], parameters, p[1], args.debug)
        params = Merge(tmp, params, args.debug)

    # Get Bc params
    for p in [ 'HeatCoeff', 'DT' ]:
        if args.debug: print(f"extract bc params for {p}")
        for bc_p in targetdefs[p]['params']:
            if args.debug: print(f"{bc_p[0]}")
            tmp = getparam(bc_p[0], parameters, bc_p[1], args.debug)
            bc_params = Merge(tmp, bc_params, args.debug)
        if args.debug: print(f"extract bc control_params for {p}")
        for bc_p in targetdefs[p]['control_params']:
            if args.debug: print(f"{bc_p[0]}")
            tmp = getparam(bc_p[0], parameters, bc_p[1], args.debug)
            bc_params = Merge(tmp, bc_params, args.debug)
    
    # print("params:", params)
    if args.debug:
        print("bc_params:", bc_params)

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
    (params, bcparams) = solve(feelpp_env, feel_pb, args, 'I', params, control_params, bc_params, targets)
    
    # update
    update(cwd, jsonmodel, params, control_params, bcparams, currents[0], args.debug)

    # display csv results
    # TODO use units 
    # get Power for Insert, Bitter
    if feelpp_env.isMasterRank():
        print(f"update: workingdir={ os.getcwd() }")
    df = getTarget('Power', feelpp_env, args.debug)
    # df = getTarget('MeanT', feelpp_env, args.debug)
    # df = getTarget('MaxT', feelpp_env, args.debug)
    if feelpp_env.isMasterRank():
        print(f"I: {currents[0]} [A]\tPower: {df.iloc[-1][0]} [W]")
        # TODO add bitter current and power if any Bitters

    # stats by Helices
    for p in ['PowerH', 'MeanTH', 'MaxTH', 'Flux']:
        df = getTarget(p, feelpp_env, args.debug)

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
            
        df.rename(columns=nkeys, inplace=True)
        if feelpp_env.isMasterRank():
            print(f"{p} [{targetdefs[p]['unit']}]:\n{tabulate(df, headers='keys', tablefmt='simple')}\n")
            
    # Same for Bitters
    # Same for Supras

if __name__ == "__main__":
    sys.exit(main())
