import json
import os

import feelpp
from feelpp.toolboxes.core import *
from feelpp.toolboxes.cfpdes import *

import argparse
import configparser

parser = argparse.ArgumentParser()
parser.add_argument("cfgfile", help="input cfg file (ex. HL-31.cfg)")
args = parser.parse_args()

pwd = os.getcwd()
print(f"cwd: {pwd}")

jsonfile = ""
meshfile = ""
feelpp_config = configparser.ConfigParser()
with open(args.cfgfile, "r") as inputcfg:
    feelpp_config.read_string("[DEFAULT]\n[main]\n" + inputcfg.read())

    basedir = os.path.dirname(args.cfgfile)
    print(f"basedir={basedir}")
    if not basedir:
        basedir = "."

    jsonfile = feelpp_config["cfpdes"]["filename"]
    jsonfile = jsonfile.replace(r"$cfgdir/", f"{basedir}/")
    print(f"jsonfile={jsonfile}")

    meshfile = feelpp_config["cfpdes"]["mesh.filename"]
    meshfile = meshfile.replace(r"$cfgdir/", f"{basedir}/")
    print(f"meshfile={meshfile}")

app = feelpp.Environment(["myapp"], config=feelpp.localRepository(""))
comm = feelpp.Environment.worldCommPtr()

m2d = feelpp.load(feelpp.mesh(dim=2), name=meshfile, verbose=1)
Xh = feelpp.functionSpace(space="Pch", mesh=m2d, order=1)
usave = Xh.element()


with open(f"{pwd}/{jsonfile}", "r") as file:
    data = json.load(file)
    for param in data["Parameters"]:
        if param.startswith("U_"):
            expr = data["Parameters"][param]
            usave.on(
                range=feelpp.markedelements(m2d, param[2:]), expr=feelpp.expr(f"{expr}")
            )

usave.save(f"{pwd}/{basedir}", name="U")
print(f"odir={pwd}/{basedir}")
