# Cas HIV Hits

This repository contains the code for the Nominate Diversify Narrow Filter (NDNF) pipeline for picking HIV gRNAs.
This pipeline is designed to pick gRNAs for Cas molecules based on the PAM and protospacer length.

## Organization

The code in the repository is organized as follows:

 - `code/` & `tests/` - Python library code and associated tests.
 - `data/` - Directory of data files for pipeline runs.
 - `notebooks/` - Example notebooks describing the NDNF pipeline and mutational data.
 - `results/` - Experimental results associated with the analysis performed for the Fronteers 2023 figures.
 - `scripts/` - Utility cli scripts that perform the nominate and filter steps.

Code requirements can be installed from the `requirements.conda` file provided.

Once DVC has been installed.
`make prepare` will download the remaining genomes.
`make test` will run all unit tests to ensure your installation is correct.


## Frontiers Manuscript

### Figure 1: Cas9 editor influences the number of safe, broad, and effective (SBE) targets. 
Generated using the notebooks in `results/pipeline_runs/`. These can also be used as templates for one's own exploration.

### Figure 2: Higher promiscuity leads to more broad targets but fewer safe ones.
Generated using the notebooks in `results/mismatch_effect/`, with `mismatch_maker.ipynb` responsible for generating the simulation data and `mismatch_figure.ipynb` for generating the visualization.

### Figure 3: The pattern of targetable regions is altered by PAM mutations.
Generated using the notebooks in `results/mutation_exploration/`, with `variant_maker.ipynb` responsible for generating the simulation data and `variant_figure.ipynb` for generating the visualization.

### Figure 4: Increase PAM specificity leads to a decrease in SBE gRNAs and loss of targetable sites.
Generated using the notebooks in `results/broadVsafe/`, with `additional_cas.ipynb` responsible for generating the simulation data and `broad_vs_safe.ipynb` for generating the visualization.

### Figure 5: Nomination and validation stages are unimpacted by random samplings.
Generated using the notebooks in `results/stability/`, with `pipeline_stability_exp.ipynb` responsible for generating the simulation data and `pipeline_stability_fig.ipynb` for generating the visualization.

### Figure 6: HIV-1 subtype distribution of full genomic sequences in the LANL dataset.
Generated using the notebook `notebooks/train_test_split.ipynb`.

### Figure 7: The effect of mutations on HIV-1 replication capacity.
Generated using the notebook `notebooks/mutation_scoring.ipynb`.

In order to relicate the Frontiers figures, the notebooks should be run in the following order to ensure intermediate files are properly created.

* `notebooks/train_test_split.ipynb` - Which will generate the random samplings for the other notebooks.
* `notebooks/mutation_scoring.ipynb` - Which will process the RC index file for other notebooks.
* Any notebooks in the `results/` folder will work.