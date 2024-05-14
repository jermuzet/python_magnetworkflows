import sys
import os
import argparse
import configparser


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "type",
        help="choose cli or commissioning",
        type=str,
        choices=["cli", "commissioning"],
    )
    parser.add_argument(
        "--cfgfile",
        help="input cfg file",
        type=str,
    )
    parser.add_argument(
        "--mdata",
        help="mdata for cli",
        type=str,
    )
    parser.add_argument("--np", help="number of cores", type=int, default=1)
    parser.add_argument(
        "--coolings",
        help="what coolings do you want",
        nargs="+",
        metavar="coolings",
        type=str,
        choices=["mean", "meanH", "grad", "gradH", "gradHZ"],
        default=["mean", "meanH", "grad", "gradH", "gradHZ"],
    )
    parser.add_argument(
        "--hcorrelations",
        help="what heat correlations do you want",
        nargs="+",
        metavar="hcorrelations",
        type=str,
        choices=["Montgomery", "Dittus", "Colburn", "Silverberg"],
        default=["Montgomery"],
    )
    parser.add_argument(
        "--frictions",
        help="what frictions do you want",
        nargs="+",
        metavar="frictions",
        type=str,
        choices=["Constant", "Blasius", "Filonenko", "Colebrook", "Swanee"],
        default=["Constant"],
    )
    parser.add_argument(
        "--itermax",
        help="specify maximum iteration (default: 10)",
        type=int,
        default=100,
    )
    parser.add_argument("--debug", help="activate debug", action="store_true")
    args = parser.parse_args()

    feelpp_config = configparser.ConfigParser()
    basedir = None
    with open(args.cfgfile, "r") as inputcfg:
        feelpp_config.read_string("[DEFAULT]\n[main]\n" + inputcfg.read())
        feelpp_directory = feelpp_config["main"]["directory"]
        feelpp_directory = feelpp_directory.split("/")
        basedir = os.path.dirname(args.cfgfile)
        cfgfile = args.cfgfile.replace(basedir + "/", "")
        basedir = basedir.replace("/" + basedir.split("/")[-1], "")
        meshmodel = feelpp_config["cfpdes"]["mesh.filename"]
        if meshmodel.endswith(".msh"):
            meshmodel = meshmodel.replace(".msh", "")
        elif meshmodel.endswith(".json"):
            meshmodel = meshmodel.split("_p")[0]
        if args.debug:
            print("feelpp_directory=", feelpp_directory)
            print("basedir=", basedir)
            print("cfgfile=", cfgfile)
            print("meshmodel=", meshmodel)

    with open(f"{basedir}/run_{args.type}_matrix.sh", "w") as f:
        for cooling in args.coolings:
            f.write(f"# Cooling={cooling}\n")
            f.write(
                f"perl -pi -e 's|mesh.filename=.*|mesh.filename=\{meshmodel}_p{args.np}.json|' {basedir}/{cooling}/{cfgfile};\n\n"
            )

            heatcorrelations = ["Montgomery"]
            frictions = ["Constant"]
            if "H" in cooling:
                heatcorrelations = args.hcorrelations
                frictions = args.frictions

            for heatcorrelation in heatcorrelations:
                f.write(f"## Heat Correlation={heatcorrelation}\n")

                for friction in frictions:
                    f.write(f"### Friction={friction}\n")
                    f.write(
                        f"perl -pi -e 's|directory=.*|directory={feelpp_directory[0]}/{feelpp_directory[1]}/{cooling}/{heatcorrelation}/{friction}|' {basedir}/{cooling}/{cfgfile};\n"
                    )
                    commandline = f"mpirun -np {args.np} -bind-to core python -m python_magnetworkflows.{args.type} --mdata '{args.mdata}' {basedir}/{cooling}/{cfgfile} --reloadcfg --cooling {cooling} --heatcorrelation {heatcorrelation} --friction {friction} --itermax {args.itermax}"
                    if args.debug:
                        commandline = commandline + " --debug "
                    f.write(
                        f"{commandline} > {basedir}/{cooling}/{cooling}_{heatcorrelation}_{friction}.log 2>&1;\n\n"
                    )

    os.system(f"sh {basedir}/run_{args.type}_matrix.sh")
    if not args.debug:
        os.remove(f"{basedir}/run_{args.type}_matrix.sh")

    return 0


if __name__ == "__main__":
    sys.exit(main())
