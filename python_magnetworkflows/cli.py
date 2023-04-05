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
import numpy as np
from tabulate import tabulate

from .params import targetdefs, getparam, Merge, setTarget, update_U
from .solver import init
from .real_methods import flow_params
# from ..units import load_units

from .oneconfig import oneconfig

def main():
    # print(f'sys.argv({type(sys.argv)})={sys.argv}')

    epilog = "Setup for cfpdes static formulations\n" \
             "Support: only magnet of insert type\n" \
             "Workflow: actually fix current and compute cooling BCs using Montgomery correlation with an average water velocity\n" \
             "\n" \
             "Before running adapt flow_params to your magnet setup \n"

    command_line = None
    parser = argparse.ArgumentParser(description="Cfpdes model", epilog=epilog)
    parser.add_argument("cfgfile", help="input cfg file (ex. HL-31.cfg)")
    parser.add_argument("--wd", help="set a working directory", type=str, default="")
    
    # TODO: make a group: oneconfig
    # TODO make current a dict: magnet_name: value, targetkey: value -> '{ "magnet_name": {"value": current value, "type": "helix" or "bitter"} }'
    # ex dict :'{"Bitter_M10":{"value":11930,"type":"bitter"},"H6":{"value":11939,"type":"helix"}}'
    parser.add_argument("--mdata", help="specify requested current as dict ", type=json.loads)
    parser.add_argument("--current", help="specify requested current (default: 31kA)", nargs='+', metavar='Current', type=float, default=[31.e+3])

    # TODO: make a group: commissionning
    parser.add_argument("--current_from", help="specify starting currents", nargs='+', metavar='Current', type=float)
    parser.add_argument("--current_to", help="specify stopping current", nargs='+', metavar='Current', type=float)
    parser.add_argument("--current_step", help="specify requested number current steps", nargs='+', metavar='Current', type=int)

    parser.add_argument("--cooling", help="choose cooling type", type=str,
                    choices=['mean', 'grad', 'meanH', 'gradH'], default='mean')
    parser.add_argument("--eps", help="specify requested tolerance (default: 1.e-3)", type=float, default=1.e-3)
    parser.add_argument("--itermax", help="specify maximum iteration (default: 10)", type=int, default=10)
    parser.add_argument("--debug", help="activate debug", action='store_true')
    parser.add_argument("--verbose", help="activate verbose", action='store_true')

    # TODO make flow a dict: magnet_name: flow_params.json
    parser.add_argument("--flow_params", help="select flow param json file", type=str, default="flow_params.json")

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

    # Load cfg as config
    feelpp_config = configparser.ConfigParser()
    with open(args.cfgfile, 'r') as inputcfg:
        feelpp_config.read_string('[DEFAULT]\n[main]\n' + inputcfg.read())
        feelpp_directory = feelpp_config['main']['directory']
        """
        for section in feelpp_config['main']:
            print(f"feelpp_cfg/main/{section}: {feelpp_config['main'][section]}")
        if args.debug:
            print(f"feelpp_cfg: {feelpp_config['main']}")

            for section in feelpp_config.sections():
                print(f"section: {section}")
        """

        jsonmodel = feelpp_config['cfpdes']['filename']
        jsonmodel = jsonmodel.replace(r"$cfgdir/",'')
        meshfile = feelpp_config['cfpdes']['mesh.filename']
        meshfile = meshfile.replace(r"$cfgdir/",'')

    # args.mdata = currents:  {magnet.name: {'value': current.value, 'type': magnet.type}} ex:'{"Bitter_M10":{"value":11930,"type":"bitter"},"H6":{"value":11939,"type":"helix"}}'
    if args.mdata:
        import copy
        
        print(f'mdata: {args.mdata}')
        for key in args.mdata:
            data = args.mdata[key]
            key_type=data['type']
            if key_type == 'helix':
                # change rematch, params, control_params
                targetdefs[f'I_{key_type}'] = copy.deepcopy(targetdefs['I'])
            if key_type == 'bitter':
                # change rematch, params, control_params
                targetdefs[f'I_{key_type}'] = copy.deepcopy(targetdefs['I'])
                targetdefs[f'I_{key_type}']['rematch'] = f'Statistics_Intensity_B\w+_integrate'
                targetdefs[f'I_{key_type}']['params'] =  [('N',f'N_B\w+')]
                targetdefs[f'I_{key_type}']['control_params'] = [('U', f'U_B\w+', update_U)]

    # Get Parameters from JSON model file
    parameters = {}
    with open(jsonmodel, 'r') as jsonfile:
        dict_json = json.loads(jsonfile.read())
        parameters = dict_json['Parameters']

    if args.mdata:
        params_mdata=[]
        control_params_mdata=[]
        for key in args.mdata:
            data = args.mdata[key]
            key_type=data['type']
            params = {}
            control_params = []
            for p in targetdefs[f'I_{key_type}']['control_params']:
                if args.debug:
                    print(f"extract control params for {p[0]}")
                control_params.append(p[0])
                tmp = getparam(p[0], parameters, p[1], args.debug)
                params = Merge(tmp, params, args.debug)

            for p in targetdefs[f'I_{key_type}']['params']:
                if args.debug:
                    print(f"extract compute params for {p[0]}")
                tmp = getparam(p[0], parameters, p[1], args.debug)
                params = Merge(tmp, params, args.debug)

            params_mdata.append(params)
            control_params_mdata.append(control_params)
            
    else :
        params = {}
        control_params = []
        for p in targetdefs['I']['control_params']:
            if args.debug:
                print(f"extract control params for {p[0]}")
            control_params.append(p[0])
            tmp = getparam(p[0], parameters, p[1], args.debug)
            params = Merge(tmp, params, args.debug)

        for p in targetdefs['I']['params']:
            if args.debug:
                print(f"extract compute params for {p[0]}")
            tmp = getparam(p[0], parameters, p[1], args.debug)
            params = Merge(tmp, params, args.debug)
    
    # Get Bc params
    bc_params = {}
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

    # init feelpp env
    feelpp_env = init(args, feelpp_directory)
    if feelpp_env.isMasterRank():
        print(f'jsonmodel={jsonmodel}')

    # print("params:", params)
    if args.debug:
        print("bc_params:", bc_params)

    if args.current_from:
        from_currents = args.current_from
        if not args.current_to:
            if feelpp_env.isMasterRank():
                print('must have --to activated')
            return 1
        if len(args.current_to) != len(args.current_from):
            if feelpp_env.isMasterRank():
                print('--to must have the same number of args as --from')
            return 1

        to_currents = args.current_to
        if not args.current_step:
            if feelpp_env.isMasterRank():
                print('must have --step activated')
            return 1

        if len(args.current_step) != len(args.current_from):
            if feelpp_env.isMasterRank():
                print('--step must have the same number of args as --from')
            return 1
        steps_currents = args.current_step

        if len(from_currents) == 1:
            table_headers = None
            table_values = []

            current0 = from_currents[0]
            dcurrent0 = (to_currents[0]-from_currents[0])/float(steps_currents[0]-1)
            for i in range(steps_currents[0]):
                _values = []
                if feelpp_env.isMasterRank():
                    print(f'current[{i}]: {current0}')
                currents = [current0]
                
                # pvalues: dict for params key->actual target in targetdefs
                pvalues = {
                    "Flux" : "Flux",
                    "Power" : "Power",
                    "PowerH" : "PowerH",
                    "HeatCoeff" : "HeatCoeff",
                    "DT" : "DT"
                }
            
                postvalues = {
                    "T": ['MinTH', 'MeanTH', 'MaxTH', 'PowerH'],
                    "Stress": ['MinVonMises', 'MinHoop', 'MeanVonMises', 'MeanHoop', 'MaxVonMises', 'MaxHoop']
                }
                # pfields depends from method
                pfields = ['I', 'PowerH', 'MinTH', 'MeanTH', 'MaxTH', 'Flux',
                        'MinVonMises', 'MeanVonMises', 'MaxVonMises',
                        'MinHoop', 'MeanHoop', 'MaxHoop']
                targets = setTarget('I', params, current0, args.debug)
                objectifs = [('I', currents[0], params, control_params, bc_params, targets, pfields, postvalues, pvalues)]
                results = oneconfig(cwd, feelpp_env, jsonmodel, meshfile, args, objectifs)
                (table_headers, _values) = results['I']
                table_values.append(_values)
                current0 += dcurrent0

            if feelpp_env.isMasterRank():
                print(f"Commissionning:\n{tabulate(table_values, headers=table_headers, tablefmt='simple')}\n")
                with open('commissionning.csv', 'w+') as f:
                    f.write(tabulate(table_values, headers=table_headers, tablefmt='tsv'))
            
        elif len(from_currents) > 2:
            if feelpp_env.isMasterRank():
                print(f'magnetworkflow/cli: {len(from_currents)} magnets detected - not implement right now')
    else:
        if feelpp_env.isMasterRank():
            print(f"Current: {args.current}")
        currents = []
        if isinstance(args.current, list):
            currents = args.current
        else:
            currents.append(currents)

        # pvalues: dict for params key->actual target in targetdefs
        pvalues = {
            "Flux" : "Flux",
            "Power" : "Power",
            "PowerH" : "PowerH",
            "HeatCoeff" : "HeatCoeff",
            "DT" : "DT"
        }
            
        postvalues = {
            "T": ['MinTH', 'MeanTH', 'MaxTH', 'PowerH'],
            "Stress": ['MinVonMises', 'MinHoop', 'MeanVonMises', 'MeanHoop', 'MaxVonMises', 'MaxHoop']
        }

        if args.mdata:
            objectifs=[]
            mdata_i=0
            for key in args.mdata :
                data = args.mdata[key]
                key_type=data['type']
                # pfields depends from method
                pfields = [f'I_{key_type}', 'PowerH', 'MinTH', 'MeanTH', 'MaxTH', 'Flux',
                            'MinVonMises', 'MeanVonMises', 'MaxVonMises',
                            'MinHoop', 'MeanHoop', 'MaxHoop']
                targets = setTarget(f'I_{key_type}', params_mdata[mdata_i], data['value'], args.debug)
                objectif = (f'I_{key_type}', data['value'], params_mdata[mdata_i], control_params_mdata[mdata_i], bc_params, targets, pfields, postvalues, pvalues)
                
                objectifs.append(objectif)
                mdata_i+=1
        else:
            # pfields depends from method
            pfields = ['I', 'PowerH', 'MinTH', 'MeanTH', 'MaxTH', 'Flux',
                        'MinVonMises', 'MeanVonMises', 'MaxVonMises',
                        'MinHoop', 'MeanHoop', 'MaxHoop']
            targets = setTarget('I', params, currents[0], args.debug)
            objectifs = [('I', currents[0], params, control_params, bc_params, targets, pfields, postvalues, pvalues)]

        oneconfig(cwd, feelpp_env, jsonmodel, meshfile, args, objectifs)

    return 0

if __name__ == "__main__":
    sys.exit(main())
