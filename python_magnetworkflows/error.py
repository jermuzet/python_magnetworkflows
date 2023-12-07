import pandas as pd

from .params import getTarget, getparam
from .waterflow import waterflow as w
from .cooling import rho, Cp, Uw, getDT, getHeatCoeff, getTout

import re
from tabulate import tabulate

import copy


def compute_error(
    e,
    f,
    basedir: str,
    it: int,
    args,
    targets: dict,
    postvalues: dict,
    params: dict,
    parameters: dict,
    dict_df: dict,
):
    """
    it: actual iteration number
    args:
    e: feelpp env
    f: feelpp problem
    targets: dict of target
    params: dict(target, params:list of parameters name)
    parameters: all jsonmodel parameters
    dict_df:
    """
    print(f"compute_error: it={it}")

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
            if args.debug:
                print(f"param={param}, marker={marker}")
                print(
                    f"{it}: {marker}, goal={objectif:.3f}, val={val:.3f}, err={error[marker].iloc[-1]:.3e}, ovalue={ovalue:.3f}, nvalue={nvalue:.3f}",
                    flush=True,
                )
            # f.addParameterInModelProperties(param, nvalue)
            parameters[param] = nvalue

        table_.append(err_max_target)

        # update bcs
        p_params = {}
        # TODO upload p_df into a dict like {name: p_df} with name = target (aka key of targets)
        # this way we can get p_df per target as an output for solve
        # NB: pd_df[pname] is a pandas dataframe (pname is parameter name, eg Flux)

        for param in values["computed_params"]:
            name = param["name"]
            print(f"{target}: computed_params {name}", flush=True)

            if "csv" in param:
                dict_df[target][name] = getTarget(
                    {f"{name}": param}, name, e, args.debug
                )
                if args.debug:
                    print(f"{target}: {name}={dict_df[target][name]}", flush=True)
            else:
                if args.debug:
                    print(f"{target}: {name}", flush=True)
                for p in param["params"]:
                    pname = p[0]
                    if args.debug:
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

        if args.debug:
            print(f"p_df: {dict_df[target].keys()}")
            print(f"p_params: {p_params.keys()}")
            print(f'p_params[Tw]={p_params["Tw"]}', flush=True)

        for key in ["statsT", "statsTH"]:
            for param in postvalues[target][key]:
                name = param["name"]
                if args.debug:
                    print(f"{target}: postvalues_params {name}", flush=True)

                if "csv" in param:
                    dict_df[target][key][name] = getTarget(
                        {f"{name}": param}, name, e, args.debug
                    )

        if args.debug:
            print(f"PowerM: {dict_df[target]['PowerM']}", flush=True)
            print(f"PowerH: {dict_df[target]['PowerH']}", flush=True)
            print(f"Flux: {dict_df[target]['Flux']}", flush=True)

        PowerM = dict_df[target]["PowerM"].iloc[-1, 0]
        SPower_H = dict_df[target]["PowerH"].iloc[-1].sum()
        SFlux_H = dict_df[target]["Flux"].iloc[-1].sum()
        t_headers = ["Part", "Flux[W]"]
        t_parts = dict_df[target]["Flux"].columns.values.tolist()
        t_power = dict_df[target]["Flux"].iloc[-1]
        print(
            tabulate(list(zip(t_parts, t_power)), headers=t_headers),
            flush=True,
        )

        Powers_Diff = abs(PowerM - SPower_H)
        PowerFlux_Diff = abs(PowerM - SFlux_H)

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

        if Powers_Diff / PowerM > 1e-3:
            return f"Power!=SPower_H:{Powers_Diff/PowerM}    Power={PowerM:.3f} SPower_H={SPower_H:.3f}"
        if Powers_Diff / PowerM > 1e-3:
            return f"Power!=SPower_H:{Powers_Diff/PowerM}    Power={PowerM:.3f} SPower_H={SPower_H:.3f}"
        if PowerFlux_Diff / PowerM > 1e-1:
            return f"Power!=SFlux_H:{PowerFlux_Diff/PowerM}    Power={PowerM:.3f} SFlux_H={SFlux_H:.3f}"

        # assert (
        #    PowerFlux_Diff / PowerM < 1e-1
        # ), f"Power!=SFlux_H:{PowerFlux_Diff/PowerM}    Power={PowerM:.3f} SFlux_H={SFlux_H:.3f}"
        # assert (
        #    Powers_Diff / PowerM < 1e-3
        # ), f"Power!=SPower_H:{Powers_Diff/PowerM}    Power={PowerM:.3f} SPower_H={SPower_H:.3f}"
        # assert (
        #    PowerFlux_Diff / PowerM < 1e-1
        # ), f"Power!=SFlux_H:{PowerFlux_Diff/PowerM}    Power={PowerM:.3f} SFlux_H={SFlux_H:.3f}"

        # get dict_df[target]["Flux"] column names
        for i, cname in enumerate(dict_df[target]["Flux"].columns.values.tolist()):
            flux = dict_df[target]["Flux"][cname].iloc[-1]
            print(f"{target} Channel{i} Flux[cname={cname}]: {flux:.3f}")

        flow = values["waterflow"]
        Pressure = flow.pressure(abs(objectif))
        dPressure = flow.dpressure(abs(objectif))

        Dh = [parameters[p] for p in p_params["Dh"]]
        Sh = [parameters[p] for p in p_params["Sh"]]
        if args.debug:
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
            if args.debug:
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
                    if nsections != len(key_dz) + 1:
                        return f"inconsistant data for Tw and FluxZ for {cname}"
                    # assert (
                    #     nsections == len(key_dz) + 1
                    # ), f"inconsistant data for Tw and FluxZ for {cname}"
                    if args.debug and e.isMasterRank():
                        print(f"key_dz={key_dz}")
                    FluxCh_dz = [
                        FluxZ.at[FluxZ.index[-1], f"FluxZ{i}_{cname}"]
                        for i in range(len(key_dz))
                    ]
                    PowerCh = sum(FluxCh_dz)

                    print(
                        f'\nFlux/Sum(FluxZ)[{cname}]: {PowerCh:.3f} = {dict_df[target]["Flux"][cname].iloc[-1]:.3f}',
                        flush=True,
                    )
                    if args.debug:
                        for k, flux in enumerate(FluxCh_dz):
                            print(f"{key_dz[k]} FluxZ[{k}]={flux:.3f}")
                    if (
                        abs(1 - PowerCh / dict_df[target]["Flux"][cname].iloc[-1])
                        > 1e-1
                    ):
                        return f"Sum(FluxZ)!=Flux[{cname}]: PowerCh={PowerCh} Flux_H={dict_df[target]['Flux'][cname].iloc[-1]} (res={PowerCh / dict_df[target]['Flux'][cname].iloc[-1]})"
                    # assert (
                    #     abs(1 - PowerCh / dict_df[target]["Flux"][cname].iloc[-1])
                    #     < 1e-1
                    # ), f"Sum(FluxZ)!=Flux[{cname}]: PowerCh={PowerCh} Flux_H={dict_df[target]['Flux'][cname].iloc[-1]} (res={PowerCh / dict_df[target]['Flux'][cname].iloc[-1]})"

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
                        if args.debug:
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

                            print(
                                f"gradHZ: {cname} {key_dz[k]} - flow={tmp_flow:.3f}, power={flux:.3f}, Tw={_old[k]:.3f}, dTw={dT_old:.3f}, P={Pressure:.3f}, dT={dT_new:.3f}, dTch={dT_Ch:.3f}"
                            )
                            _new[k + 1] = _new[k] + dT_new

                        # print(f"_new: {_new}, relax={relax}")
                        for k in range(section):
                            # print(f's={s}, {1-relax} * {_new[s]} + {relax} * {_old[s]}')
                            _new[k] = (1 - relax) * _new[k] + relax * _old[k]
                        print(f"gradHZ: {cname} _new (after relax={relax}): {_new}")

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
                        if args.debug:
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
                    Tw_data.to_csv(f"{csvfile}", index=False)
                    e.worldComm().barrier()

                    # f.addParameterInModelProperties(p_params["dTwH"][i], dTwi[-1])
                    # f.addParameterInModelProperties(p_params["hwH"][i], hi[i])
                    parameters[p_params["hwH"][i]] = hi[i]
                    # parameters[p_params["dTwH"][i]] = dTwi[-1]
                    dict_df[target]["HeatCoeff"][p_params["hwH"][i]] = [round(hi[i], 3)]
                    dict_df[target]["DT"][cname] = [round(dTwi[i], 3)]
                    dict_df[target]["Uw"]["Uw_" + cname] = [round(U, 3)]
                    # !! export FluxZ = dict_df[target]["Flux"] with sections regrouped !!
                    # dict_df[target]["Flux"][cname] = PowerCh

                    error_dT.append(abs(1 - (dTwH[i] / dTwi[i])))
                    error_h.append(abs(1 - (hwH[i] / hi[i])))

                    Tw0 = Tw_data["Tw"].iloc[0]

                    print(
                        f"{target} Cooling[{i}]: cname={cname}, u={U:.3f}, Dh={d:3e}, Sh={s:.3e}, Power={PowerCh:.3f}, TwH={Tw0:.3f}, dTwH={dTwH[i]:.3f}, hwH={hwH[i]:.3f}, dTwi={dTwi[i]:.3f}, hi={hi[i]:.3f}",
                        flush=True,
                    )

                    Ti[i] = Tw0 + dTwi[i]
                    VolMass[i] = rho(Tw0 + dTwi[i] / 2.0, Pressure)
                    SpecHeat[i] = Cp(Tw0 + dTwi[i] / 2.0, Pressure)

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
                    print(
                        f"\ncname={cname}, i={i}, d={d:.5f}, s={s:.5f}, L={Lh[i]:.3f}, U={U:.3f}, PowerCh={PowerCh:.3f}, tmp_dTwi={tmp_dTwi:.3f}, tmp_hi={tmp_hi:.3f}",
                        flush=True,
                    )

                    while True:
                        tmp_flow = U * s
                        tmp_dTwi = getDT(tmp_flow, PowerCh, TwH[i], dTwH[i], Pressure)

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
                        if args.debug:
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

                    # f.addParameterInModelProperties(p_params["dTwH"][i], dTwi[i])
                    # f.addParameterInModelProperties(p_params["hwH"][i], hi[i])
                    parameters[p_params["hwH"][i]] = hi[i]
                    parameters[p_params["dTwH"][i]] = dTwi[i]
                    dict_df[target]["HeatCoeff"][p_params["hwH"][i]] = [round(hi[i], 3)]
                    dict_df[target]["DT"][p_params["dTwH"][i]] = [round(dTwi[i], 3)]
                    dict_df[target]["Uw"][p_params["dTwH"][i].replace("dTw", "Uw")] = [round(U, 3)]

                    error_dT.append(abs(1 - (dTwH[i] / dTwi[i])))
                    error_h.append(abs(1 - (hwH[i] / hi[i])))

                    print(
                        f"{target} Cooling[{i}]: cname={cname}, u={U:.3f}, Dh={d:.3e}, Sh={s:.3e}, Power={PowerCh:.3f}, TwH={TwH[i]:.3f}, dTwH={dTwH[i]:.3f}, hwH={hwH[i]:.3f}, dTwi={dTwi[i]:.3f}, hi={hi[i]:.3f}",
                        flush=True,
                    )

                    VolMass[i] = rho(TwH[i] + dTwi[i] / 2.0, Pressure)
                    SpecHeat[i] = Cp(TwH[i] + dTwi[i] / 2.0, Pressure)

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
                if args.debug:
                    print(
                        f"T:{T} Tw:{Tw}  i:{i}  dTw:{dTw[i]}  hw:{hw[i]}",
                        flush=True,
                    )
                dTg = getDT(flow.flow(abs(objectif)), PowerM, Tw[i], dTw[i], Pressure)
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
                # f.addParameterInModelProperties(p_params["dTw"][i], dTg)
                # f.addParameterInModelProperties(p_params["hw"][i], hg)
                parameters[p_params["hw"][i]] = hg
                parameters[p_params["dTw"][i]] = dTg

                error_dT.append(abs(1 - (dTw[i] / dTg)))
                error_h.append(abs(1 - (hw[i] / hg)))

                dict_df[target]["HeatCoeff"][p_params["hw"][i]] = [round(hg, 3)]
                dict_df[target]["DT"][p_params["dTw"][i]] = [round(dTg, 3)]
                dict_df[target]["Uw"][p_params["dTw"][i].replace("dTw", "Uw")] = [round(Umean, 3)]

            if args.debug:
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
        print(f"MSITE Tout={Tout_site}", flush=True)

    return (err_max, err_max_dT, err_max_h, table_, p_params)
