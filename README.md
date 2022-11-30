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

====
