# Worflows

## Running in container

If python_magnetsetup is installed in the container
just type:

```
singularity exec /home/singularity/feelpp-toolboxes-v0.110.0-alpha.3.sif \
    mpirun -np 2 python -m python_magnetsetup.workflows.cli HL-test-cfpdes-thelec-Axi-sim.cfg --eps 1.e-5
```

Otherwise just copy workflows on machine in the simulation directory

```
scp -r python_magnetsetup/workflows machine
```

Then on machine:

```
singularity exec /home/singularity/feelpp-toolboxes-v0.110.0-alpha.3.sif \
    mpirun -np 2 python -m workflows.cli HL-test-cfpdes-thelec-Axi-sim.cfg --eps 1.e-5
```

[NOTE]
====
To check whether or not python_magnetsetup is installed in the container
Run:

```
singularity exec /home/singularity/feelpp-toolboxes-v0.110.0-alpha.3.sif \
    dpkg -L python3-magnetsetup
```

====
