from typing import List

import re
import os
import shutil
import copy

import json
from tabulate import tabulate

import feelpp as fpp
import feelpp.toolboxes.core as tb
import feelpp.toolboxes.cfpdes as cfpdes

import pandas as pd

from .params import getTarget, getparam
from .waterflow import waterflow as w
from .cooling import rho, Cp, Uw, getDT, getHeatCoeff, getTout


def create_field(
    feel_pb, targets: dict, parameters: dict, wd: str, debug: bool = False
):
    """
    create and save a field in h5"""
    print("create_field", flush=True)

    # TODO: pass space and order in method params
    Xh = fpp.functionSpace(space="Pdh", mesh=feel_pb.mesh(), order=0)
    usave = Xh.element()

    for target, values in targets.items():
        params = []
        for p in values["control_params"]:
            tmp = getparam(p[0], parameters, p[1], debug)
            params += tmp

        for param in params:
            marker = param.replace("U_", "")  # 'U_': (values['control_params'][0])[0]
            # print(f'marker: {marker}, goal={values["goals"][marker].iloc[-1]} xx')
            value = parameters[param]
            usave.on(
                range=fpp.markedelements(feel_pb.mesh(), marker),
                expr=fpp.expr(str(value)),
            )

    usave.save(wd, name="U")
    print("create_field: done")


def init_field(e, jsonmodel: str, meshmodel: str):
    """
    create and save a field in h5"""
    # print("init_field", flush=True)
    if e.isMasterRank():
        print(f"init_field: jsonmodel={jsonmodel}, meshmodel={meshmodel}")

    m2d = fpp.load(fpp.mesh(dim=2), name=meshmodel, verbose=1)
    Xh = fpp.functionSpace(space="Pdh", mesh=m2d, order=0)
    usave = Xh.element()

    basedir = os.path.dirname(jsonmodel)
    with open(jsonmodel, "r") as file:
        data = json.load(file)
        for param in data["Parameters"]:
            if param.startswith("U_"):
                value = data["Parameters"][param]
                usave.on(
                    range=fpp.markedelements(m2d, param[2:]), expr=fpp.expr(str(value))
                )

    usave.save(basedir, name="U")
    # print("init_field: done", flush=True)


def update(e, jsonmodel: str, parameters: dict):
    """
    Update jsonmodel with parameters
    """
    if e.isMasterRank():
        print(f"update {jsonmodel}", flush=True)

    dict_json = {}
    with open(jsonmodel, "r") as jsonfile:
        dict_json = json.loads(jsonfile.read())
        dict_json["Parameters"] = parameters

    if e.isMasterRank():
        with open(jsonmodel, "w+") as jsonfile:
            jsonfile.write(json.dumps(dict_json, indent=4))
    e.worldComm().barrier()

    return 0


# TODO create toolboxes_options on the fly
def init(fname, e, args, jsonmodel: str, meshmodel: str, directory: str = ""):
    pwd = os.getcwd()
    if not e:
        e = fpp.Environment(
            [f"{fname}.py"],
            opts=tb.toolboxes_options("coefficient-form-pdes", "cfpdes"),
            config=fpp.localRepository(directory),
        )
        if e.isMasterRank():
            print(f"init: feelpp env created (pwd={pwd}, cwd={os.getcwd()})")

    init_field(e, f"{pwd}/{jsonmodel}", f"{pwd}/{meshmodel}")

    e.setConfigFile(args.cfgfile)
    if e.isMasterRank():
        print(f"init: setconfigfile", flush=True)

    # if not comm is None:
    #    f = cfpdes.cfpdes(dim=2, worldComm=comm)
    # else:
    f = cfpdes.cfpdes(dim=2)
    if e.isMasterRank():
        # print(f"init: pwd={pwd}")
        print("Create cfpdes done", flush=True)

    f.init()
    if e.isMasterRank():
        print("Init problem done", flush=True)
    # f.printAndSaveInfo()

    return (e, f)


def solve(
    fname,
    e,
    f,
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

    if e is None:
        (e, f) = init(fname, e, args, jsonmodel, meshmodel, directory=feelpp_directory)
        if e.isMasterRank():
            print("solve: load cfg", flush=True)

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

    # save original U.h5
    save_h5 = f"{basedir}/U.h5.init"
    isFile = os.path.isfile(save_h5)
    e.worldComm().barrier()
    if isFile:
        raise RuntimeError(
            f"solve: backup U.h5 to {save_h5} - fails since file already exists"
        )
    else:
        # Rename the file
        if e.isMasterRank():
            shutil.copy2(f"{basedir}/U.h5", save_h5)
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

    # Xh = fpp.functionSpace(space="Pch", mesh=f.mesh(), order=1)
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
                    shutil.copy2(f"{_file}", f"{_file}-it{it}.csv")
                e.worldComm().barrier()

        if e.isMasterRank():
            with open(new_json, "w+") as jsonfile:
                jsonfile.write(json.dumps(dict_json, indent=4))

            shutil.copy2(f"{basedir}/U.h5", f"{basedir}/U-it{it}-{post[:-1]}.h5")
            # print(f"start calc for it={it}", flush=True)
            if args.debug:
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
        table_ = [it]
        err_max = 0
        err_max_dT = 0
        err_max_h = 0
        List_Tout = []
        List_VolMassout = []
        List_SpecHeatout = []
        List_Qout = []
        Tw0 = None
        for target, values in targets.items():
            objectif = -values["objectif"]
            # multiply by -1 because of orientation of pseudo Axi domain Oy == -U_theta
            filtered_df = getTarget(targets, target, e, args.debug)

            relax = values["relax"]

            # TODO: add stats for filtered_df to table_: mean, ecart type, min/max??

            error = filtered_df.div(-objectif).add(1)
            err_max_target = max(error.abs().max(axis=1))
            err_max = max(err_max_target, err_max)
            if e.isMasterRank():
                if args.debug:
                    print(
                        f"filtered_df: {filtered_df.columns.values.tolist()}",
                        flush=True,
                    )
                    print(f"{target}: objectif={objectif}")
                    print(f"{target}: relax={relax}")
                    print(f"{target}: filtered_df={filtered_df}")
                    print(f"{target}: error={error}")
                print(
                    f"{target}: it={it}, err_max={err_max_target:.3e}, eps={args.eps:.3e}, itmax={args.itermax}",
                    flush=True,
                )

            for param in params[target]:
                marker = param.replace(
                    "U_", ""
                )  # get name from values['control_params'] / change control_params to a list of dict?
                val = filtered_df[marker].iloc[-1]
                ovalue = parameters[param]
                table_.append(ovalue)
                nvalue = ovalue * objectif / val
                if e.isMasterRank():
                    if args.debug:
                        print(f"param={param}, marker={marker}")
                        print(
                            f"{it}: {marker}, goal={objectif:.3f}, val={val:.3f}, err={error[marker].iloc[-1]:.3e}, ovalue={ovalue:.3f}, nvalue={nvalue:.3f}",
                            flush=True,
                        )
                f.addParameterInModelProperties(param, nvalue)
                parameters[param] = nvalue

            table_.append(err_max_target)

            # update bcs
            p_params = {}
            # TODO upload p_df into a dict like {name: p_df} with name = target (aka key of targets)
            # this way we can get p_df per target as an output for solve
            # NB: pd_df[pname] is a pandas dataframe (pname is parameter name, eg Flux)

            for param in values["computed_params"]:
                name = param["name"]
                if e.isMasterRank():
                    print(f"{target}: computed_params {name}", flush=True)

                if "csv" in param:
                    dict_df[target][name] = getTarget(
                        {f"{name}": param}, name, e, args.debug
                    )
                    if args.debug and e.isMasterRank():
                        print(f"{target}: {name}={dict_df[target][name]}", flush=True)
                else:
                    if args.debug and e.isMasterRank():
                        print(f"{target}: {name}", flush=True)
                    for p in param["params"]:
                        pname = p[0]
                        if args.debug and e.isMasterRank():
                            print(
                                f"{name}: extract params for {p[0]} (len(p)={len(p)})",
                                flush=True,
                            )
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

            if args.debug and e.isMasterRank():
                print(f"p_df: {dict_df[target].keys()}")
                print(f"p_params: {p_params.keys()}")
                print(f'p_params[Tw]={p_params["Tw"]}', flush=True)

            for key in ["statsT", "statsTH"]:
                for param in postvalues[target][key]:
                    name = param["name"]
                    if args.debug and e.isMasterRank():
                        print(f"{target}: postvalues_params {name}", flush=True)

                    if "csv" in param:
                        dict_df[target][key][name] = getTarget(
                            {f"{name}": param}, name, e, args.debug
                        )

            if args.debug and e.isMasterRank():
                print(f"PowerM: {dict_df[target]['PowerM']}", flush=True)
                print(f"PowerH: {dict_df[target]['PowerH']}", flush=True)
                print(f"Flux: {dict_df[target]['Flux']}", flush=True)

            PowerM = dict_df[target]["PowerM"].iloc[-1, 0]
            SPower_H = dict_df[target]["PowerH"].iloc[-1].sum()
            SFlux_H = dict_df[target]["Flux"].iloc[-1].sum()

            Powers_Diff = abs(PowerM - SPower_H)
            PowerFlux_Diff = abs(PowerM - SFlux_H)
            if e.isMasterRank():
                print(
                    f"{target}: it={it} Power={PowerM:.3f} SPower_H={SPower_H:.3f} SFlux_H={SFlux_H:.3f}",
                    flush=True,
                )
                # print a table with key, power, ucoil
                t_headers = ["Part", "Power[W]", "U[V]"]
                t_parts = dict_df[target]["PowerH"].columns.values.tolist()
                t_power = dict_df[target]["PowerH"].iloc[-1]
                t_U = dict_df[target]["PowerH"].iloc[-1] / dict_df[target]["target"]
                print(
                    tabulate(list(zip(t_parts, t_power, t_U)), headers=t_headers),
                    flush=True,
                )
                # print(
                #    f'PowerH=\n{dict_df[target]["PowerH"].iloc[-1]} \nUcoil=\n{dict_df[target]["PowerH"].iloc[-1]/dict_df[target]["target"]}',
                #    flush=True,
                # )

                if args.debug:
                    print("PowerM : ", dict_df[target]["PowerM"], flush=True)
                    print(f'{target}: PowerH {dict_df[target]["PowerH"]}', flush=True)
                    for key in p_params:
                        print(f"{target}: {key}={p_params[key]}", flush=True)

            assert (
                Powers_Diff / PowerM < 1e-3
            ), f"Power!=SPower_H:{Powers_Diff/PowerM}    Power={PowerM:.3f} SPower_H={SPower_H:.3f}"
            assert (
                PowerFlux_Diff / PowerM < 1e-1
            ), f"Power!=SFlux_H:{PowerFlux_Diff/PowerM}    Power={PowerM:.3f} SFlux_H={SFlux_H:.3f}"

            # get dict_df[target]["Flux"] column names
            if e.isMasterRank():
                for i, cname in enumerate(
                    dict_df[target]["Flux"].columns.values.tolist()
                ):
                    flux = dict_df[target]["Flux"][cname].iloc[-1]
                    print(f"{target} Channel{i} Flux[cname={cname}]: {flux:.3f}")

            flow = values["waterflow"]
            Pressure = flow.pressure(abs(objectif))
            dPressure = flow.dpressure(abs(objectif))

            Dh = [parameters[p] for p in p_params["Dh"]]
            Sh = [parameters[p] for p in p_params["Sh"]]
            if args.debug and e.isMasterRank():
                i = 0
                for p in p_params["Dh"]:
                    print(f"Dh[{i}]: key={p}, value={parameters[p]}", flush=True)
                    i += 1
                print(f'Dh: {p_params["Dh"]}')

            Umean = flow.umean(abs(objectif), sum(Sh))  # math.fsum(Sh)
            dT_global = getDT(
                Umean * sum(Sh),
                PowerM,
                290.75,
                0,
                Pressure,
            )
            if e.isMasterRank():
                print(
                    f"{target}: it={it}, objectif={abs(objectif):.3f}, Umean={Umean:.3f}, Flow={flow.flow(abs(objectif)):.3f}, S={sum(Sh):.3f}, dT_global={dT_global:.3f}, Power={PowerM:.3f}",
                    flush=True,
                )
            dict_df[target]["flow"] = flow.flow(abs(objectif))

            error_dT = []
            error_h = []
            # per Channel/Slit
            if "H" in args.cooling:
                TwH = [parameters[p] for p in p_params["TwH"]]
                if e.isMasterRank():
                    for p in p_params["TwH"]:
                        print(f"TwH: parameters[{p}]={parameters[p]}")
                # print(f'TwH: {TwH} ({len(TwH)}/{len(Dh)})')
                dTwH = [parameters[p] for p in p_params["dTwH"]]
                # print(f"dTwH: {dTwH} ({len(dTwH)}/{len(Dh)})")
                hwH = [parameters[p] for p in p_params["hwH"]]
                # print(f'hwH: {hwH} ({len(hwH)}/{len(Dh)})')
                Lh = [
                    abs(parameters[p] - parameters[p.replace("max", "min")])
                    for p in p_params["ZmaxH"]
                ]
                if args.debug and e.isMasterRank():
                    for i, p in enumerate(p_params["hwH"]):
                        print(f"hwH[{i}]: key={p}, value={parameters[p]}", flush=True)

                    # TODO verify if data are consistant??
                    # assert len(Dh) == len(TwH) == len(dTwH) == len(hwH)
                    print(f'{target} Flux: {dict_df[target]["Flux"]}')

                if not dTwH:
                    dTwH = [0] * len(Dh)

                dTwi = [0] * len(Dh)
                Ti = [0] * len(Dh)
                hi = [0] * len(Dh)
                VolMass = [0] * len(Dh)
                SpecHeat = [0] * len(Dh)
                Q = [0] * len(Dh)

                if "Z" in args.cooling:
                    NCoolingCh = len(Dh)
                    FluxZ = dict_df[target]["FluxZ"]

                    for i, (d, s) in enumerate(zip(Dh, Sh)):
                        cname = p_params["Dh"][i].replace("Dh_", "")

                        U = Umean

                        # load T_z_old from csv
                        # print(f"TwH[{i}]: {TwH[i]}")
                        csvfile = TwH[i]["filename"].replace("$cfgdir", basedir)
                        # replace $cfgdir by actual value from feelpp environment e
                        # print(f"cwd={os.getcwd()}, csvfile={csvfile}")
                        Tw_data = pd.read_csv(csvfile, sep=",", engine="python")
                        nsections = len(Tw_data)
                        _old = Tw_data["Tw"].to_list()
                        dTwH[i] = _old[-1] - _old[0]
                        section = len(_old)

                        # get PowerCh_Z
                        key_dz = [
                            fkey
                            for fkey in FluxZ.columns.values.tolist()
                            if fkey.endswith(cname)
                        ]
                        assert (
                            nsections == len(key_dz) + 1
                        ), f"inconsistant data for Tw and FluxZ for {cname}"
                        if args.debug and e.isMasterRank():
                            print(f"key_dz={key_dz}")
                        FluxCh_dz = [
                            FluxZ.at[FluxZ.index[-1], f"FluxZ{i}_{cname}"]
                            for i in range(len(key_dz))
                        ]
                        PowerCh = sum(FluxCh_dz)
                        if e.isMasterRank():
                            print(
                                f'\nFlux/Sum(FluxZ)[{cname}]: {PowerCh:.3f} = {dict_df[target]["Flux"][cname].iloc[-1]:.3f}',
                                flush=True,
                            )
                            if args.debug:
                                for k, flux in enumerate(FluxCh_dz):
                                    print(f"{key_dz[k]} FluxZ[{k}]={flux:.3f}")
                        assert (
                            abs(1 - PowerCh / dict_df[target]["Flux"][cname].iloc[-1])
                            < 1e-1
                        ), f"Sum(FluxZ)!=Flux[{cname}]: PowerCh={PowerCh} Flux_H={dict_df[target]['Flux'][cname].iloc[-1]} (res={PowerCh / dict_df[target]['Flux'][cname].iloc[-1]})"

                        while True:
                            tmp_flow = U * s

                            dT_Ch = getDT(
                                tmp_flow,
                                PowerCh,
                                (_old[0] + _old[-1]) / 2.0,
                                0,
                                Pressure,
                            )
                            _new = copy.deepcopy(_old)
                            if args.debug and e.isMasterRank():
                                print(f"gradHZ: {cname} _old: {_old}")
                            for k, flux in enumerate(FluxCh_dz):
                                dT_old = _old[k + 1] - _old[k]
                                dT_new = getDT(
                                    tmp_flow,
                                    flux,
                                    (_old[k] + _old[k + 1]) / 2.0,
                                    0,
                                    Pressure,
                                )
                                if e.isMasterRank():
                                    print(
                                        f"gradHZ: {cname} {key_dz[k]} - flow={tmp_flow:.3f}, power={flux:.3f}, Tw={_old[k]:.3f}, dTw={dT_old:.3f}, P={Pressure:.3f}, dT={dT_new:.3f}, dTch={dT_Ch:.3f}"
                                    )
                                _new[k + 1] = _new[k] + dT_new

                            # print(f"_new: {_new}, relax={relax}")
                            for k in range(section):
                                # print(f's={s}, {1-relax} * {_new[s]} + {relax} * {_old[s]}')
                                _new[k] = (1 - relax) * _new[k] + relax * _old[k]
                            if e.isMasterRank():
                                print(
                                    f"gradHZ: {cname} _new (after relax={relax}): {_new}"
                                )

                            dTwi[i] = _new[-1] - _new[0]
                            tmp_hi = getHeatCoeff(
                                d,
                                Lh[i],
                                U,
                                _new[0] + dTwi[-1] / 2.0,
                                hwH[i],
                                Pressure,
                                dPressure,
                                model=args.heatcorrelation,
                                friction=args.friction,
                            )

                            tmp_U = Uw(
                                _new[0] + dTwi[-1] / 2.0,
                                Pressure,
                                dPressure,
                                d,
                                Lh[i],
                                friction=args.friction,
                                uguess=U,
                            )
                            n_tmp_flow = tmp_U * s
                            if args.debug and e.isMasterRank():
                                print(
                                    f"tmp_flow={tmp_flow:.3f}, n_tmp_flow={n_tmp_flow:.3f}, U={U:.3f}, tmp_U={tmp_U:.3f}, dTwH[{i}]={dTwH[i]:.3f}, tmp_dTwi={tmp_dTwi:.3f}"
                                )
                            U = tmp_U

                            if abs(1 - n_tmp_flow / tmp_flow) <= 1.0e-3:
                                break

                            _old = copy.deepcopy(_new)

                        Q[i] = tmp_flow
                        hi[i] = (1.0 - relax) * tmp_hi + relax * hwH[i]

                        # save back to csv: T_z.to_csv(f"Tw_{cname}.csv", index=False)
                        Tw_data["Tw"] = _new
                        # print(f'save _new={_new} to {csvfile}')
                        if e.isMasterRank:
                            Tw_data.to_csv(f"{csvfile}", index=False)
                        e.worldComm().barrier()

                        # f.addParameterInModelProperties(p_params["dTwH"][i], dTwi[-1])
                        f.addParameterInModelProperties(p_params["hwH"][i], hi[i])
                        parameters[p_params["hwH"][i]] = hi[i]
                        # parameters[p_params["dTwH"][i]] = dTwi[-1]
                        dict_df[target]["HeatCoeff"][p_params["hwH"][i]] = [hi[i]]
                        dict_df[target]["DT"][cname] = [dTwi[i]]
                        # !! export FluxZ = dict_df[target]["Flux"] with sections regrouped !!
                        # dict_df[target]["Flux"][cname] = PowerCh

                        error_dT.append(abs(1 - (dTwH[i] / dTwi[i])))
                        error_h.append(abs(1 - (hwH[i] / hi[i])))

                        Tw0 = Tw_data["Tw"].iloc[0]
                        if e.isMasterRank():
                            print(
                                f"{target} Cooling[{i}]: cname={cname}, u={U:.3f}, Dh={d:3e}, Sh={s:.3e}, Power={PowerCh:.3f}, TwH={Tw0:.3f}, dTwH={dTwH[i]:.3f}, hwH={hwH[i]:.3f}, dTwi={dTwi[i]:.3f}, hi={hi[i]:.3f}",
                                flush=True,
                            )

                        Ti[i] = Tw0 + dTwi[i]
                        VolMass[i] = rho(Tw0 + dTwi[i] / 2.0, Pressure)
                        SpecHeat[i] = Cp(Tw0 + dTwi[i] / 2.0, Pressure)

                    if e.isMasterRank():
                        print(
                            f"Flow={flow.flow(abs(objectif)):.3f}, Sum(Qh)={sum(Q):.3f}, Sum(Sh)={sum(Sh):.3e}",
                            flush=True,
                        )

                    # compute an estimate of dTg
                    Tout = getTout(Ti, VolMass, SpecHeat, Q)
                    VolMassout = rho(Tout, Pressure)
                    SpecHeatout = Cp(Tout, Pressure)
                    Qout = Umean * sum(Sh)

                    List_Tout.append(Tout)
                    List_VolMassout.append(VolMassout)
                    List_SpecHeatout.append(SpecHeatout)
                    List_Qout.append(Qout)

                    dTg = Tout - Tw0
                    if e.isMasterRank():
                        print(
                            f"{target}: Tout={Tout:.3f}  Tw={Tw0:.3f}, U={Umean:.3f}, Power={PowerM:.3f}, dTg={dTg:.3f} ({PowerM/(VolMass[0]*SpecHeat[0]*flow.flow(abs(objectif))):.3f})",
                            flush=True,
                        )
                    dict_df[target]["Tout"] = Tout
                    # exit(1)

                else:
                    for i, (d, s) in enumerate(zip(Dh, Sh)):
                        cname = p_params["Dh"][i].replace("Dh_", "")
                        PowerCh = dict_df[target]["Flux"][cname].iloc[-1]

                        U = Umean
                        tmp_dTwi = dTwH[i]
                        tmp_hi = hwH[i]
                        if e.isMasterRank():
                            print(
                                f"\ncname={cname}, i={i}, d={d:.5f}, s={s:.5f}, L={Lh[i]:.3f}, U={U:.3f}, PowerCh={PowerCh:.3f}, tmp_dTwi={tmp_dTwi:.3f}, tmp_hi={tmp_hi:.3f}",
                                flush=True,
                            )

                        while True:
                            tmp_flow = U * s
                            tmp_dTwi = getDT(
                                tmp_flow, PowerCh, TwH[i], dTwH[i], Pressure
                            )

                            tmp_hi = getHeatCoeff(
                                d,
                                Lh[i],
                                U,
                                TwH[i] + tmp_dTwi / 2.0,
                                hwH[i],
                                Pressure,
                                dPressure,
                                model=args.heatcorrelation,
                                friction=args.friction,
                            )

                            tmp_U = Uw(
                                TwH[i] + tmp_dTwi / 2.0,
                                Pressure,
                                dPressure,
                                d,
                                Lh[i],
                                friction=args.friction,
                                uguess=U,
                            )
                            n_tmp_flow = tmp_U * s
                            if args.debug and e.isMasterRank():
                                print(
                                    f"tmp_flow={tmp_flow:.3f}, n_tmp_flow={n_tmp_flow:.3f}, U={U:.3f}, tmp_U={tmp_U:.3f}, dTwH[{i}]={dTwH[i]:.3f}, tmp_dTwi={tmp_dTwi:.3f}"
                                )
                            U = tmp_U

                            if abs(1 - n_tmp_flow / tmp_flow) <= 1.0e-3:
                                break

                        Q[i] = tmp_flow
                        dTwi[i] = (1.0 - relax) * tmp_dTwi + relax * dTwH[i]
                        Ti[i] = TwH[i] + dTwi[i]
                        hi[i] = (1.0 - relax) * tmp_hi + relax * hwH[i]

                        f.addParameterInModelProperties(p_params["dTwH"][i], dTwi[i])
                        f.addParameterInModelProperties(p_params["hwH"][i], hi[i])
                        parameters[p_params["hwH"][i]] = hi[i]
                        parameters[p_params["dTwH"][i]] = dTwi[i]
                        dict_df[target]["HeatCoeff"][p_params["hwH"][i]] = [hi[i]]
                        dict_df[target]["DT"][p_params["dTwH"][i]] = [dTwi[i]]

                        error_dT.append(abs(1 - (dTwH[i] / dTwi[i])))
                        error_h.append(abs(1 - (hwH[i] / hi[i])))

                        if e.isMasterRank():
                            print(
                                f"{target} Cooling[{i}]: cname={cname}, u={U:.3f}, Dh={d:.3e}, Sh={s:.3e}, Power={PowerCh:.3f}, TwH={TwH[i]:.3f}, dTwH={dTwH[i]:.3f}, hwH={hwH[i]:.3f}, dTwi={dTwi[i]:.3f}, hi={hi[i]:.3f}",
                                flush=True,
                            )

                        VolMass[i] = rho(TwH[i] + dTwi[i] / 2.0, Pressure)
                        SpecHeat[i] = Cp(TwH[i] + dTwi[i] / 2.0, Pressure)

                    if e.isMasterRank():
                        print(
                            f"Flow={flow.flow(abs(objectif)):.3f}, Sum(Qh)={sum(Q):.3f}, Sum(Sh)={sum(Sh):.3e}",
                            flush=True,
                        )

                    # TODO compute an estimate of dTg
                    # Tout /= VolMass * SpecHeat * (Umean * sum(Sh))
                    Tout = getTout(Ti, VolMass, SpecHeat, Q)
                    VolMassout = rho(Tout, Pressure)
                    SpecHeatout = Cp(Tout, Pressure)
                    Qout = Umean * sum(Sh)

                    List_Tout.append(Tout)
                    List_VolMassout.append(VolMassout)
                    List_SpecHeatout.append(SpecHeatout)
                    List_Qout.append(Qout)

                    dTg = Tout - TwH[0]
                    if e.isMasterRank():
                        print(
                            f"{target}: Tout={Tout:.3f}  Tw={TwH[0]:.3f}, U={Umean:.3f}, Power={PowerM:.3f}, dTg={dTg:.3f} ({PowerM/(VolMass[0]*SpecHeat[0]*flow.flow(abs(objectif))):.3f})",
                            flush=True,
                        )
                    dict_df[target]["Tout"] = Tout

            # global:  what to do when len(Tw) != 1
            else:
                Tw = [parameters[p] for p in p_params["Tw"]]
                dTw = [parameters[p] for p in p_params["dTw"]]
                hw = [parameters[p] for p in p_params["hw"]]
                L = [
                    abs(parameters[p] - parameters[p.replace("max", "min")])
                    for p in p_params["Zmax"]
                ]

                for i, T in enumerate(Tw):
                    if args.debug and e.isMasterRank():
                        print(
                            f"T:{T} Tw:{Tw}  i:{i}  dTw:{dTw[i]}  hw:{hw[i]}",
                            flush=True,
                        )
                    dTg = getDT(
                        flow.flow(abs(objectif)), PowerM, Tw[i], dTw[i], Pressure
                    )
                    hg = getHeatCoeff(
                        Dh[i],
                        L[i],
                        Umean,
                        Tw[i] + dTg / 2.0,
                        hw[i],
                        Pressure,
                        dPressure,
                        model=args.heatcorrelation,
                        friction=args.friction,
                    )
                    f.addParameterInModelProperties(p_params["dTw"][i], dTg)
                    f.addParameterInModelProperties(p_params["hw"][i], hg)
                    parameters[p_params["hw"][i]] = hg
                    parameters[p_params["dTw"][i]] = dTg

                    error_dT.append(abs(1 - (dTw[i] / dTg)))
                    error_h.append(abs(1 - (hw[i] / hg)))

                    dict_df[target]["HeatCoeff"][p_params["hw"][i]] = [hg]
                    dict_df[target]["DT"][p_params["dTw"][i]] = [dTg]

                if args.debug and e.isMasterRank():
                    print(
                        f'{target}: Tw={Tw[0]}, param={p_params["dTw"][0]}, umean={Umean}, Power={PowerM}, dTg={dTg}, hg={hg}',
                        flush=True,
                    )
                dict_df[target]["Tout"] = Tw[0] + dTg

                VolMassout = rho(Tw[0] + dTg, Pressure)
                SpecHeatout = Cp(Tw[0] + dTg, Pressure)
                Qout = flow.flow(abs(objectif))

                List_Tout.append(Tw[0] + dTg)
                List_VolMassout.append(VolMassout)
                List_SpecHeatout.append(SpecHeatout)
                List_Qout.append(Qout)

            # TODO: how to transform dTg, hg et DTwi, hi en dataframe??

            err_max_dT = max(err_max_dT, max(error_dT))
            err_max_h = max(err_max_h, max(error_h))

        if len(List_Tout) > 1:
            Tout_site = getTout(List_Tout, List_VolMassout, List_SpecHeatout, List_Qout)

            if e.isMasterRank():
                print(f"MSITE Tout={Tout_site}", flush=True)

        # update Parameters
        f.updateParameterValues()

        if e.isMasterRank():
            print(
                f"it={it}, err_max={err_max:.3e}, err_max_dT={err_max_dT:.3e}, err_max_h={err_max_h:.3e}, eps={args.eps}, itmax={args.itermax}",
                flush=True,
            )

        table_.append(err_max)
        table.append(table_)
        create_field(f, targets, parameters, basedir, args.debug)
        if e.isMasterRank():
            # print("create_field : done", flush=True)
            update(e, jsonmodel, parameters)

        if (
            err_max <= args.eps
            and err_max_dT <= max(args.eps, 1e-2)
            and err_max_h <= max(args.eps, 1e-2)
        ):
            break

        # reload feelpp
        if args.reloadcfg:
            print("reload cfg")
            e.worldComm().barrier()
            e.setConfigFile(args.cfgfile)
            f = cfpdes.cfpdes(dim=2)
            f.init()

        it += 1

    # Save table (need headers)
    if e.isMasterRank():
        print(tabulate(table, headers, tablefmt="simple"), flush=True)

    table_df = pd.DataFrame(table, columns=headers)

    if err_max > args.eps or it >= args.itermax:
        raise RuntimeError(f"Fail to solve {jsonmodel}: err_max={err_max}, it={it}")

    if e.isMasterRank():
        os.remove(save_h5)
        os.remove(save_json)
        for i in range(it):
            print(f"remove {basedir}/U-it{i}-{post[:-1]}.h5")
            os.remove(f"{basedir}/U-it{i}-{post[:-1]}.h5")
            tmp_jsonmodel = jsonmodel.replace(".json", f"-it{i}-{post[:-1]}.json")
            print(f"remove {tmp_jsonmodel}")
            os.remove(tmp_jsonmodel)
        if "Z" in args.cooling:
            for i in range(1, it):
                for file in csvfiles:
                    _file = file.replace("$cfgdir", basedir)
                    print(f"remove {_file}-it{it}.csv")
                    os.remove(f"remove {_file}-it{it}.csv")

    return (table_df, dict_df, e)
