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
    /home/singularity/hifimagnet-salome-9.8.3.sif \
    salome -w1 -t /opt/SALOME-9.8.0-UB20.04/INSTALL/HIFIMAGNET/bin/salome/HIFIMAGNET_Cmd.py args:M9_M19020601.yaml,--axi,--air,4,6
```

====
