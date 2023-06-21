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
from .cooling import rho, Cp
from .real_methods import getDT, getHeatCoeff, getTout


def create_field(
    feel_pb, targets: dict, parameters: dict, wd: str, debug: bool = False
):
    """
    create and save a field in h5"""
    # print("create_field")

    # TODO: pass space and order in method params
    Xh = fpp.functionSpace(space="Pch", mesh=feel_pb.mesh(), order=1)
    usave = Xh.element()

    for target, values in targets.items():
        params = []
        for p in values["control_params"]:
            tmp = getparam(p[0], parameters, p[1], debug)
            params += tmp

        for param in params:
            marker = param.replace("U_", "")  # 'U_': (values['control_params'][0])[0]
            # print(f'marker: {marker}, goal={values["goals"][marker].iloc[-1]} xx')
            value = float(parameters[param])
            usave.on(
                range=fpp.markedelements(feel_pb.mesh(), marker),
                expr=fpp.expr(str(value)),
            )

    usave.save(wd, name="U")
    # print("create_field: done")


def create_field_init(jsonmodel: str, meshmodel: str):
    """
    create and save a field in h5"""
    # print("create_field")

    m2d = fpp.load(fpp.mesh(dim=2), name=meshmodel, verbose=1)
    Xh = fpp.functionSpace(space="Pch", mesh=m2d, order=1)
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
    # print("create_field: done")


def update(jsonmodel: str, parameters: dict):
    """
    Update jsonmodel with parameters
    """
    print(f"update {jsonmodel}")

    dict_json = {}
    with open(jsonmodel, "r") as jsonfile:
        dict_json = json.loads(jsonfile.read())
        dict_json["Parameters"] = parameters

    with open(jsonmodel, "w+") as jsonfile:
        jsonfile.write(json.dumps(dict_json, indent=4))

    return 0


# TODO create toolboxes_options on the fly
def init(e, args, jsonmodel: str, meshmodel: str, directory: str = ""):
    # pwd = os.getcwd()
    if not e:
        e = fpp.Environment(
            ["cli.py"],
            opts=tb.toolboxes_options("coefficient-form-pdes", "cfpdes"),
            config=fpp.localRepository(directory),
        )

    create_field_init(jsonmodel, meshmodel)

    e.setConfigFile(args.cfgfile)

    if e.isMasterRank():
        # print(f"init: pwd={pwd}")
        print("Create cfpdes")
    f = cfpdes.cfpdes(dim=2)

    if e.isMasterRank():
        print("Init problem")
    f.init()
    # f.printAndSaveInfo()

    if e.isMasterRank():
        print(f"init: done")
    return (e, f)


def solve(
    e,
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

    if args.debug:
        print(f"solve: jsonmodel={jsonmodel}, args={args}")

    # suffix of tmp files
    post = ""
    for target, values in targets.items():
        post += f'{target}={values["objectif"]}{values["unit"]}-'

    basedir = os.path.dirname(jsonmodel)  # get absolute path instead??

    (e, f) = init(e, args, jsonmodel, meshmodel, directory=feelpp_directory)

    if e.isMasterRank():
        print(f"solve: jsonmodel={jsonmodel}, basedir={basedir}")

        # save original jsonmodel
        save_json = f"{jsonmodel}.init"
        if os.path.isfile(save_json):
            raise RuntimeError(
                f"solve: backup jsonmodel to {save_json} - fails since file already exists"
            )
        else:
            # Rename the file
            shutil.copy2(jsonmodel, save_json)

    # save original U.h5
    if e.isMasterRank():
        save_h5 = f"{basedir}/U.h5.init"
        if os.path.isfile(save_h5):
            raise RuntimeError(
                f"solve: backup U.h5 to {save_h5} - fails since file already exists"
            )
        else:
            # Rename the file
            shutil.copy2(f"{basedir}/U.h5", save_h5)

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
            print(f"make a copy of files for it={it}")
            print(f"jsonmodel={jsonmodel}, new_json={new_json}")
        # shutil.copy2(jsonmodel, new_json)
        dict_json = {}
        with open(jsonmodel, "r") as jsonfile:
            dict_json = json.loads(jsonfile.read())
            dict_json["Meshes"]["cfpdes"]["Fields"]["U"][
                "filename"
            ] = f"$cfgdir/U-it{it}.h5"
            # print(f'dict_json={dict_json}')
        with open(new_json, "w+") as jsonfile:
            jsonfile.write(json.dumps(dict_json, indent=4))

        shutil.copy2(f"{basedir}/U.h5", f"{basedir}/U-it{it}.h5")

        if e.isMasterRank():
            print(f"start calc for it={it}")
            if args.debug:
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
        err_max_dT = 0
        err_max_h = 0
        List_Tout = []
        List_VolMassout = []
        List_SpecHeatout = []
        List_Qout = []
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
                    print(f"{target}: objectif={objectif}")
                    print(f"{target}: relax={relax}")
                    print(f"{target}: filtered_df={filtered_df}")
                    print(f"{target}: error={error}")
                print(
                    f"{target}: it={it}, err_max={err_max_target}, eps={args.eps}, itmax={args.itermax}"
                )

            for param in params[target]:
                marker = param.replace(
                    "U_", ""
                )  # get name from values['control_params'] / change control_params to a list of dict?
                val = filtered_df[marker].iloc[-1]
                ovalue = float(parameters[param])
                table_.append(ovalue)
                nvalue = ovalue * objectif / val
                if e.isMasterRank():
                    print(
                        f"{it}: {marker}, goal={objectif}, val={val}, err={error[marker].iloc[-1]}, ovalue={ovalue}, nvalue={nvalue}"
                    )
                f.addParameterInModelProperties(param, nvalue)
                parameters[param] = nvalue

                # usave.on(
                #     range=fpp.markedelements(f.mesh(), marker),
                #     expr=fpp.expr(str(nvalue)),
                # )

            table_.append(err_max_target)

            # update bcs
            if e.isMasterRank():
                print(f"{target}: it={it}, update BCs")
            p_params = {}
            # TODO upload p_df into a dict like {name: p_df} with name = target (aka key of targets)
            # this way we can get p_df per target as an output for solve
            # NB: pd_df[pname] is a pandas dataframe (pname is parameter name, eg Flux)

            for param in values["computed_params"]:
                name = param["name"]
                if args.debug and e.isMasterRank():
                    print(f"{target}: computed_params {name}")

                if "csv" in param:
                    dict_df[target][name] = getTarget(
                        {f"{name}": param}, name, e, args.debug
                    )
                    if args.debug and e.isMasterRank():
                        print(f"{target}: {name}={dict_df[target][name]}")
                else:
                    if args.debug and e.isMasterRank():
                        print(f"{target}: {name} computed")
                    for p in param["params"]:
                        pname = p[0]
                        if args.debug and e.isMasterRank():
                            print(
                                f"{name}: extract params for {p[0]} (len(p)={len(p)})"
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
                print(f'p_params[Tw]={p_params["Tw"]}')

            for key in ["statsT", "statsTH"]:
                for param in postvalues[target][key]:
                    name = param["name"]
                    if args.debug and e.isMasterRank():
                        print(f"{target}: postvalues_params {name}")

                    if "csv" in param:
                        dict_df[target][key][name] = getTarget(
                            {f"{name}": param}, name, e, args.debug
                        )

            PowerM = dict_df[target]["PowerM"].iloc[-1, 0]
            SPower_H = dict_df[target]["PowerH"].iloc[-1].sum()
            SFlux_H = dict_df[target]["Flux"].iloc[-1].sum()

            Powers_Diff = abs(PowerM - SPower_H)
            PowerFlux_Diff = abs(PowerM - SFlux_H)
            if e.isMasterRank():
                # Powers_Diff = abs(PowerM - SPower_H)
                # PowerFlux_Diff = abs(PowerM - SFlux_H)

                print(
                    f'{target}: it={it} Power={PowerM} SPower_H={SPower_H} SFlux_H={SFlux_H} \nPowerH=\n{dict_df[target]["PowerH"].iloc[-1]} \nUcoil=\n{dict_df[target]["PowerH"].iloc[-1]/dict_df[target]["target"]}'
                )

                print(
                    f"{target}: it={it} Power-SPower_H={Powers_Diff}     Power-SFlux_H={PowerFlux_Diff}"
                )

                if args.debug:
                    print("PowerM : ", dict_df[target]["PowerM"])
                    print(f'{target}: PowerH {dict_df[target]["PowerH"]}')
                    for key in p_params:
                        print(f"{target}: {key}={p_params[key]}")

            assert (
                Powers_Diff / PowerM < 1e-3
            ), f"Power!=SPower_H:{Powers_Diff/PowerM}    Power={PowerM} SPower_H={SPower_H}"
            assert (
                PowerFlux_Diff / PowerM < 1e-1
            ), f"Power!=SFlux_H:{PowerFlux_Diff/PowerM}    Power={PowerM} SFlux_H={SFlux_H}"

            flow = values["waterflow"]
            Pressure = flow.pressure(abs(objectif))
            dPressure = flow.dpressure(abs(objectif))

            Dh = [float(parameters[p]) for p in p_params["Dh"]]
            Sh = [float(parameters[p]) for p in p_params["Sh"]]
            if args.debug and e.isMasterRank():
                i = 0
                for p in p_params["Dh"]:
                    print(f"Dh[{i}]: key={p}, value={float(parameters[p])}")
                    i += 1
                print(f'Dh: {p_params["Dh"]}')

            Umean = flow.umean(abs(objectif), sum(Sh))  # math.fsum(Sh)
            if e.isMasterRank():
                print(
                    f"{target}: it={it}, objectif={abs(objectif)}, Umean={Umean}, Flow={flow.flow(abs(objectif))}"
                )
            dict_df[target]["flow"] = flow.flow(abs(objectif))

            error_dT = []
            error_h = []
            # per Channel/Slit
            if "H" in args.cooling:
                TwH = [float(parameters[p]) for p in p_params["TwH"]]
                dTwH = [float(parameters[p]) for p in p_params["dTwH"]]
                hwH = [float(parameters[p]) for p in p_params["hwH"]]
                Lh = [
                    abs(
                        float(parameters[p])
                        - float(parameters[p.replace("max", "min")])
                    )
                    for p in p_params["ZmaxH"]
                ]
                if args.debug and e.isMasterRank():
                    i = 0
                    for p in p_params["hwH"]:
                        print(f"hwH[{i}]: key={p}, value={float(parameters[p])}")
                        i += 1

                # TODO verify if data are consistant??
                if e.isMasterRank():
                    print(
                        f"{target}: len(Dh)={len(Dh)}, len(TwH)={len(TwH)}, len(dTwH)={len(dTwH)}, len(Channels/Slits)={len(hwH)}"
                    )
                    if args.debug:
                        print(f'{target} Flux: {dict_df[target]["Flux"]}')

                dTwi = []
                Ti = []
                hi = []
                VolMass = []
                SpecHeat = []
                Q = []
                for i, (d, s) in enumerate(zip(Dh, Sh)):
                    cname = p_params["Dh"][i].replace("Dh_", "")
                    PowerCh = dict_df[target]["Flux"][cname].iloc[-1]
                    Q.append(Umean * s)
                    dTwi.append(
                        getDT(Q[-1], PowerCh, TwH[i], dTwH[i], Pressure, relax=relax)
                    )
                    Ti.append(TwH[i] + dTwi[-1])
                    hi.append(
                        getHeatCoeff(
                            d,
                            Lh[i],
                            Umean,
                            TwH[i] + dTwi[-1] / 2.0,
                            hwH[i],
                            Pressure,
                            dPressure,
                            model=args.heatcorrelation,
                            relax=relax,
                        )
                    )
                    f.addParameterInModelProperties(p_params["dTwH"][i], dTwi[-1])
                    f.addParameterInModelProperties(p_params["hwH"][i], hi[-1])
                    parameters[p_params["hwH"][i]] = hi[-1]
                    parameters[p_params["dTwH"][i]] = dTwi[-1]
                    dict_df[target]["HeatCoeff"][p_params["hwH"][i]] = [hi[-1]]
                    dict_df[target]["DT"][p_params["dTwH"][i]] = [dTwi[-1]]

                    error_dT.append(abs(1 - (dTwH[i] / dTwi[-1])))
                    error_h.append(abs(1 - (hwH[i] / hi[-1])))

                    if e.isMasterRank():
                        print(
                            f"{target} {i}: cname={cname}, umean={Umean}, Dh={d}, Sh={s}, Power={PowerCh}, TwH={TwH[i]}, dTwH={dTwH[i]}, hwH={hwH[i]}, dTwi={dTwi[i]}, hi={hi[i]}"
                        )

                    VolMass.append(rho(TwH[i] + dTwi[-1] / 2.0, Pressure))
                    SpecHeat.append(Cp(TwH[i] + dTwi[-1] / 2.0, Pressure))

                if args.debug and e.isMasterRank():
                    print(f"Q-flow={sum(Q)-flow.flow(abs(objectif))}")
                    print(
                        f"len(Ti)={len(Ti)}  len(VolMass)={len(VolMass)} len(SpecHeat)={len(SpecHeat)} len(Dh)={len(Dh)}  len(Sh)={len(Sh)}"
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
                        f"{target}: Tout={Tout}  Tw={TwH[0]}, umean={Umean}, Power={PowerM}, dTg={dTg} ({PowerM/(VolMass[0]*SpecHeat[0]*flow.flow(abs(objectif)))})"
                    )
                dict_df[target]["Tout"] = Tout

            # global:  what to do when len(Tw) != 1
            else:
                Tw = [float(parameters[p]) for p in p_params["Tw"]]
                dTw = [float(parameters[p]) for p in p_params["dTw"]]
                hw = [float(parameters[p]) for p in p_params["hw"]]
                L = [
                    abs(
                        float(parameters[p])
                        - float(parameters[p.replace("max", "min")])
                    )
                    for p in p_params["Zmax"]
                ]

                for i, T in enumerate(Tw):
                    if e.isMasterRank():
                        print(f"T:{T} Tw:{Tw}  i:{i}  dTw:{dTw[i]}  hw:{hw[i]}")
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
                        f'{target}: Tw={Tw[0]}, param={p_params["dTw"][0]}, umean={Umean}, Power={PowerM}, dTg={dTg}, hg={hg}'
                    )
                dict_df[target]["Tout"] = Tw[0] + dTg

            # TODO: how to transform dTg, hg et DTwi, hi en dataframe??

            err_max_dT = max(err_max_dT, max(error_dT))
            err_max_h = max(err_max_h, max(error_h))

        if "H" in args.cooling and len(List_Tout) > 1:
            Tout_site = getTout(List_Tout, List_VolMassout, List_SpecHeatout, List_Qout)

            dTg = Tout_site - TwH[0]
            if e.isMasterRank():
                print(f"MSITE Tout={Tout_site}, Tw={TwH[0]}, dTg={dTg}")

        # update Parameters
        f.updateParameterValues()

        if e.isMasterRank():
            print(
                f"it={it}, err_max={err_max}, err_max_dT={err_max_dT}, err_max_h={err_max_h}, eps={args.eps}, itmax={args.itermax}"
            )

        table_.append(err_max)
        table.append(table_)
        if e.isMasterRank():
            print("create_field")
        create_field(f, targets, parameters, basedir, args.debug)
        if e.isMasterRank():
            print("create_field : done")
            update(jsonmodel, parameters)

        if (
            err_max <= args.eps
            and err_max_dT <= max(args.eps, 1e-2)
            and err_max_h <= max(args.eps, 1e-2)
        ):
            break

        # reload feelpp
        e.setConfigFile(args.cfgfile)
        f = cfpdes.cfpdes(dim=2)
        f.init()

        it += 1

    # Save table (need headers)
    if e.isMasterRank():
        print(tabulate(table, headers, tablefmt="simple"))

    table_df = pd.DataFrame(table, columns=headers)

    if err_max > args.eps or it >= args.itermax:
        raise RuntimeError(f"Fail to solve {jsonmodel}: err_max={err_max}, it={it}")
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
    if e.isMasterRank():
        os.remove(save_h5)
        os.remove(save_json)

    return (table_df, dict_df, e)

