"""
Run feelpp model for one config
"""

import os

from .params import getparam
from .solver import solve

# from ..units import load_units


def oneconfig(
    fname: str,
    e,
    f,
    fields,
    feelpp_directory: str,
    jsonmodel: str,
    meshmodel: str,
    args,
    targets: dict,
    postvalues: dict,
    parameters: dict,
):
    """
     Run a simulation until currents are reached

     e: feelpp environment
     jsonmodel:
     meshmodel:
     targets: dict of target (name: objectif , csv, ...)
     postvalues: dict for postprocessing quantities
     parameters: all jsonmodel parameters

    returns a tuple:
    table:
    dict_df:
    e:
    """
    comm = e.worldCommPtr()

    if e.isMasterRank():
        print("\n\nONECONFIG START", flush=True)

    # pwd = os.getcwd()
    basedir = os.path.dirname(args.cfgfile)
    if e.isMasterRank():
        print(f"oneconfig: jsonmodel={jsonmodel}, basedir={basedir}", flush=True)

    # capture actual params per target:
    params = {}
    for key, values in targets.items():
        params[key] = []
        for p in values["control_params"]:
            if args.debug and e.isMasterRank():
                print(f"extract control params for {p[0]}")
            tmp = getparam(p[0], parameters, p[1], args.debug)
            params[key] += tmp
            # print(f"params[{key}]={params[key]}, rank={comm.localRank()}", flush=True)

    # e.worldComm().barrier()
    (table, dict_df, e) = solve(
        fname,
        e,
        f,
        fields,
        feelpp_directory,
        jsonmodel,
        meshmodel,
        args,
        targets,
        postvalues,
        params,
        parameters,
    )

    if args.debug:
        print(f"end of oneconfig, rank={comm.localRank()}")

    e.worldComm().barrier()
    return (table, dict_df, e)
