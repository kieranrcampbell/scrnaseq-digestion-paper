
configfile: "config.yaml"

files10X = ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']

import pandas as pd

df_inventory = pd.read_csv(config['sample_inventory_url']).dropna()
df_inventory.index = df_inventory['id']
inventory_dict = df_inventory.to_dict('index')

ids = list(inventory_dict.keys())

raw_files = expand("data/raw_10X/{id}/{file}", id=ids, file=files10X)
scesets_raw = expand("data/scesets/{id}_sceset_raw.rds", id=ids)

rule all:
    input:
        raw_files,
        scesets_raw

rule get_from_shahlab:
    params:
        shahlab_path=lambda wildcards: inventory_dict[wildcards.id]['shahlab_path'],
        genome=lambda wildcards: inventory_dict[wildcards.id]['genome']
    output:
        expand("data/raw_10X/{{id}}/{file}", file=files10X)
    shell:
        "scp thost:{params.shahlab_path}/outs/filtered_gene_bc_matrices/{params.genome}/* data/raw_10X/{wildcards.id}/"

rule convert_to_sce:
    params:

    input:
        expand("data/raw_10X/{{id}}/{file}", file=files10X)
    output:
        "data/scesets/{id}_sceset_raw.rds"
    shell:
        "Rscript pipeline/conversion_to_sceset/convert_to_sceset.R \
        --input_data_path data/raw10X/{wildcards.id} \
        --output_scepath {output} \
        --sample_id"


        
    





