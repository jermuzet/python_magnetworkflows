import json
import os

import feelpp
from feelpp.toolboxes.core import *
from feelpp.toolboxes.cfpdes import *

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--cfgfile", help="give cfg name", type=str)
parser.add_argument("--odir", help="give output directory", type=str)
args = parser.parse_args()

PWD = os.getcwd()
app = feelpp.Environment(["myapp"], config=feelpp.localRepository(""))
comm = feelpp.Environment.worldCommPtr()

with open(args.cfgfile, "r") as fp:
    for line in fp:
        if line.startswith("filename="):
            jsonfile = line.replace("filename=$cfgdir", PWD).rstrip()

        if line.startswith("mesh.filename="):
            meshfile = line.replace("mesh.filename=$cfgdir", PWD).rstrip()


m2d = feelpp.load(feelpp.mesh(dim=2), name=meshfile, verbose=1)
Xh = feelpp.functionSpace(space="Pch", mesh=m2d, order=1)
usave = Xh.element()


with open(jsonfile, "r") as file:
    data = json.load(file)
    for param in data["Parameters"]:
        if param.startswith("U_"):
            expr = data["Parameters"][param]
            usave.on(
                range=feelpp.markedelements(m2d, param[2:]), expr=feelpp.expr(f"{expr}")
            )

usave.save(
    args.odir,
    name="U",
)
