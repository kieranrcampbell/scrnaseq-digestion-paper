import pandas as pd
import os
import json
import itertools

configfile: "private_config.yaml"

# Construct data frame of samples
df_inventory = pd.read_csv(config['sample_inventory_url']).dropna()

assert df_inventory.shape[0] == 45

df_inventory.index = df_inventory['id']
inventory_dict = df_inventory.to_dict('index')
ids = list(inventory_dict.keys())

cellranger_versions = config['cellranger_versions']

# all_figs and deliverables are added to by subsequent workflows
# Dict of images to populate for final output:
all_figs = {}

# Dictionary of deliverables (e.g. differential expression results, genesets, etc)
deliverables = {'murine-contamination': [config['murine_contamination_csv']]}

# Statistics
statistics = {s: 'data/statistics/{}.csv'.format(s) for s in config['statistics']}



include: 'pipeline/data-preparation/data-preparation.smk'

include: 'pipeline/murine-contamination/murine-contamination.smk'

include: 'pipeline/all-sample-overview/all-sample-overview.smk'

include: 'pipeline/live-dead-dying/live-dead-dying.smk'

include: 'pipeline/comparison-existing-10X/comparison-existing-10X.smk'

include: 'pipeline/differential-expression/differential-expression.smk'


rule all:
    input:
        sces_qc, # QC'd SingleCellExperiments
        list(itertools.chain(*all_figs.values())),
        list(itertools.chain(*deliverables.values())),
	list(statistics.values()),
	config['statfile']


rule collate_stats:
    input:
        list(statistics.values())
    output:
        config['statfile']
    shell:
        "Rscript scripts/create-latex-stats.R \
        --input_dir data/statistics \
        --output_latex {output}"