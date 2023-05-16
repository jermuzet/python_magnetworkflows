# Worflows

## Running in container

If python_magnetworkflows is installed in the container
just type:

```
singularity exec /home/singularity/feelpp-v0.110.2.sif \
    mpirun -np 2 python -m python_magnetworkflows.workflows.cli HL-test-cfpdes-thelec-Axi-sim.cfg --eps 1.e-5
```

[NOTE]
====
To check whether or not python_magnetworkflows is installed in the container
Run:

```
singularity exec /home/singularity/feelpp-v0.110.2.sif \
    dpkg -L python3-magnetworkflows
```

In devcontainer, you shall mount your data directory to /home/vscode to run singularity container,eg:

```
singularity exec -B /data/M9_M19020601-cfpdes-thmagel_hcurl-nonlinear-Axi-sim/data/geometries:/home/vscode \
    /home/singularity/hifimagnet-salome-9.8.4.sif \
    salome -w1 -t /opt/SALOME-9.8.0-UB20.04/INSTALL/HIFIMAGNET/bin/salome/HIFIMAGNET_Cmd.py args:M9_M19020601.yaml,--axi,--air,4,6
```

```
singularity exec -B /data/M9_M19020601-cfpdes-thmagel_hcurl-nonlinear-Axi-sim/data/geometries:/home/vscode \
    /home/singularity/hifimagnet-salome-9.8.4.sif \
    python3 -m python_magnetgeo.xao M9_M19020601-Axi_withAir.xao mesh  --group CoolingChannels --geo M9_M19020601.yaml []
```

```
mpirun -np 2 python -m python_magnetworkflows.workflows.cli HL-test-cfpdes-thelec-Axi-sim.cfg --eps 1.e-5
```

====

```
python -m python_magnetworkflows.cli \
  --mdata '{"HLtest":{"value":12000,"type":"helix","filter":"HLtest_","flow":"tmp/M9Bitters_HLtest/HLtest-flow_params.json"}, \
            "M9Bitters":{"value":31000,"type":"bitter","filter":"M9Bitters_","flow":"tmp/M9Bitters_HLtest/M9Bitters-flow_params.json"}}' \
  tmp/M9Bitters_HLtest/M9Bitters_HLtest-cfpdes-thelec-Axi-sim.cfg
```