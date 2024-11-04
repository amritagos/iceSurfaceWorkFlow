# Workflow for Ice Surface Calculations

## Notes for dropH2O

Now you need to create the input dictionary for the workflow in `2-dropH2O` before running the workflow.
Copy this into the `config` directory and then run the workflow as you would normally.

## Running the Snakemake with micromamba

In the top-level directory, first activate the micromamba environment
```bash
micromamba activate iceenv
```

Note: On the cluster it is a good idea to run Snakemake in a `tmux` shell (because if it is killed then the jobs die). One might do the following: 

```bash
tmux new -s my_session
micromamba activate iceenv
```
And then run Snakemake commands. To detach use <Ctrl-B then press D> and to re-attach

``` bash
tmux attach -t my_session
```

Regarding Snakemake: always do a dry run first, and use the `config.yml` inside `elja_profile`:
Go into the correct directory to run, for instance 1-dropH2O

```bash
snakemake -n --profile ../elja_profile/
```

**Caution:** Never run Snakemake on Elja without using the `config.yml` or else you will submit jobs to the login node.

To actually run the simulations, submit it like so: 

```bash
snakemake --profile ../elja_profile/
```