"""
Run feelpp model
"""

import sys
import os
import argparse
import configparser

import json

from .waterflow import waterflow
from .real_methods import getDT, getHeatCoeff

from .oneconfig import oneconfig

# import feelpp as fpp
# from .params import getparam

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
    
    # TODO: make a group: oneconfig
    # TODO make current a dict: magnet_name: value, targetkey: value
    parser.add_argument("--mdata", help="specify current data", type=json.loads)

    parser.add_argument("--cooling", help="choose cooling type", type=str,
                    choices=['mean', 'grad', 'meanH', 'gradH'], default='mean')
    parser.add_argument("--eps", help="specify requested tolerance (default: 1.e-3)", type=float, default=1.e-3)
    parser.add_argument("--itermax", help="specify maximum iteration (default: 10)", type=int, default=10)
    parser.add_argument("--debug", help="activate debug", action='store_true')
    parser.add_argument("--verbose", help="activate verbose", action='store_true')

    args = parser.parse_args()
    if args.debug:
        print(args)

    pwd = os.getcwd()
    print(f'cwd: {pwd}')

    # Load units:
    # TODO force millimeter when args.method == "HDG"
    # units = load_units('meter')

    # Load cfg as config
    jsonmodel = ""
    feelpp_config = configparser.ConfigParser()
    with open(args.cfgfile, 'r') as inputcfg:
        feelpp_config.read_string('[DEFAULT]\n[main]\n' + inputcfg.read())
        feelpp_directory = feelpp_config['main']['directory']
        print(f'feelpp_directory={feelpp_directory}')

        jsonmodel = feelpp_config['cfpdes']['filename']
        basedir = os.path.dirname(args.cfgfile)
        print(f'basedir={basedir}')
        jsonmodel = jsonmodel.replace(r"$cfgdir/",f"{basedir}/")
        print(f'jsonmodel={jsonmodel}')

    # Get Parameters from JSON model file
    parameters = {}
    with open(jsonmodel, 'r') as jsonfile:
        dict_json = json.loads(jsonfile.read())
        parameters = dict_json['Parameters']

    targets = {}

    # args.mdata = currents:  {magnet.name: {'value': current.value, 'type': magnet.type, 'filter': '', 'flow_params': args.flow_params}}
    if args.mdata:
        for mname, values in args.mdata.items():
            print(f'mname={mname}, values={values}')
            filter = values['filter']
            if values['type'] == 'helix':
                # change rematch, params, control_params
                Power = {
                    "name": "Power",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Power_{filter}integrate",
                    "post": { "type": "Statistics_Power", "math": "integrate"},
                    "unit": "W",
                }
                PowerH = {
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Power_{filter}\\w+_integrate",
                    "name": "PowerH",
                    "post": { "type": "Statistics_Power", "math": "integrate"},
                    "unit": "W",
                }

                Flux = {
                    "name": "Flux",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Flux_{filter}Channel\\d+_integrate",
                    "post": { "type": "Statistics_Flux", "math": "integrate"},
                    "unit": "W",
                }

                HeatCoeff = {
                    "name": "HeatCoeff",
                    "params": [
                        ("Dh", f"{filter}Dh\\d+"), 
                        ("Sh", f"{filter}Sh\\d+"), 
                        ("hw", f"{filter}hw"), 
                        ("h", f"{filter}h\\d+")
                        ],
                    "value": (getHeatCoeff),
                    "unit": "W/m2/K",
                }
                
                DT = {
                    "name": "DT",
                    "params": [
                        ("Tw", f"{filter}Tw"), 
                        ("dTw", f"{filter}dTw"), 
                        ("TwH", f"{filter}Tw\\d+"), 
                        ("dTwH", f"{filter}dTw\\d+")
                        ],
                    "value": (getDT),
                    "unit": "K",
                }

                targets[f'{filter}I'] = {
                    "objectif": values['value'],
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Intensity_{filter}H\\w+_integrate",
                    "params": [("N", f"N_{filter}\\w+")],
                    "control_params": [(f'{filter}U', f"U_{filter}\\w+")],
                    "computed_params": [HeatCoeff, DT, Flux, Power, PowerH],
                    "unit": "A",
                    "name": f"Intensity_{filter}",
                    "post": { "type": "Statistics_Intensity", "math": "integrate"},
                    "waterflow": waterflow.flow_params(values['flow'])
                }

            if values['type'] == 'bitter':
                # change rematch, params, control_params
                Power = {
                    "name": "Power",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Power_{filter}integrate",
                    "post": { "type": "Statistics_Power", "math": "integrate"},
                    "unit": "W",
                }
                PowerH = {
                    "name": "PowerH",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Power_{filter}\\w+_integrate",
                    "post": { "type": "Statistics_Power", "math": "integrate"},
                    "unit": "W",
                }

                Flux = {
                    "name": "Flux",
                    "csv": "heat.measures/values.csv",
                    "rematch": f"Statistics_Flux_{filter}\\w+_Slit\\d+_integrate",
                    "post": { "type": "Statistics_Flux", "math": "integrate"},
                    "unit": "W",
                }

                HeatCoeff = {
                    "name": "HeatCoeff",
                    "params": [
                        ("Dh", f"{filter}\\w+Dh"), 
                        ("Sh", f"{filter}\\w+Sh"), 
                        ("hw", f"{filter}\\w+hw", "\\w+_Slit\\d+", False),
                        ("h", f"{filter}\\w+hw", "\\w+_Slit\\d+", True)
                        ],
                    "value": (getHeatCoeff),
                    "unit": "W/m2/K",
                }
                
                DT = {
                    "name": "DT",
                    "params": [
                        ("Tw", f"{filter}\\w+_Tw", "\\w+_Slit\\d+", False), 
                        ("dTw", f"{filter}\\w+_dTw", "\\w+_Slit\\d+", False),
                        ("TwH", f"{filter}\\w+_Tw", "\\w+_Slit\\d+", True), 
                        ("dTwH", f"{filter}\\w+_dTw", "\\w+_Slit\\d+", True)
                        ],
                    "value": (getDT),
                    "unit": "K",
                }

                targets[f'{filter}I'] = {
                    "objectif": values['value'],
                    "csv": "heat.measures/values.csv",
                    'rematch': f'Statistics_Intensity_{filter}\\w+_integrate',
                    'params': [('N',f'N_{filter}\\w+')],
                    'control_params': [(f'{filter}U', f'U_{filter}\\w+')],
                    "computed_params": [HeatCoeff, DT, Flux, Power, PowerH],
                    "unit": "A",
                    "name": f"Intensity{filter}",
                    "post": { "type": "Statistics_Intensity", "math": "integrate"},
                    "waterflow": waterflow.flow_params(values['flow'])
                }

    print(f'targets: {targets.keys()}')

    
    
    """
    # pvalues: dict for params key->actual target in targets
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
    
    """

    oneconfig(feelpp_directory, jsonmodel, args, targets, parameters)

    return 0

if __name__ == "__main__":
    sys.exit(main())
