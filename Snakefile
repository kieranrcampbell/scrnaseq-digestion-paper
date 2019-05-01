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

# Dict of images to populate for final output:
all_figs = {}

# Dictionary of deliverables (e.g. differential expression results, genesets, etc)
deliverables = {}

include: 'pipeline/data-preparation/data-preparation.smk'

include: 'pipeline/murine-contamination/murine-contamination.smk'

include: 'pipeline/all-sample-overview/all-sample-overview.smk'

include: 'pipeline/live-dead-dying/live-dead-dying.smk'

include: 'pipeline/comparison-existing-10X/comparison-existing-10X.smk'

# TODO: includes for
# 1. All differential expression

print(all_figs)

rule all:
    input:
        sces_qc, # QC'd SingleCellExperiments
        config['murine_contamination_csv'],
        list(itertools.chain(*all_figs.values())),
        list(itertools.chain(*deliverables.values()))
