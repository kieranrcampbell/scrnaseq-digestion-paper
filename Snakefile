import pandas as pd
import os
import json


configfile: "config.yaml"
files10X = ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']


df_inventory = pd.read_csv(config['sample_inventory_url']).dropna()

assert df_inventory.shape[0] == 42

df_inventory.index = df_inventory['id']
inventory_dict = df_inventory.to_dict('index')


ids = list(inventory_dict.keys())

raw_files = expand("data/raw_10X/{id}/{file}", id=ids, file=files10X)
scesets_raw = expand("data/scesets/{id}_sceset_raw.rds", id=ids)
scesets_qc = expand("data/scesets/{id}_sceset_qc.rds", id=ids)

mito_des = expand("data/mito_differential_expression/{id}_mito_de.csv",
id = ids)


rule all:
    input:
        raw_files,
        scesets_raw,
        scesets_qc

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
        metadata_json=lambda wildcards: json.dumps(inventory_dict[wildcards.id])
    input:
        expand("data/raw_10X/{{id}}/{file}", file=files10X)
    output:
        "data/scesets/{id}_sceset_raw.rds"
    shell:
        "Rscript pipeline/conversion_to_sceset/convert_to_sceset.R \
        --input_data_path data/raw_10X/{wildcards.id}/ \
        --output_scepath {output} \
        --metadata_json '{params.metadata_json}'"

rule qc_scesets:
    params:
        curr_dir = os.getcwd()
    input:
        "data/scesets/{id}_sceset_raw.rds"
    output:
        sce="data/scesets/{id}_sceset_qc.rds",
        report="reports/qc/qc_report_{id}.html"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/qc/sce_qc.Rmd', \
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}', \
        params=list(input_sce_path='{input}', output_sce_path='{output.sce}'))\" "

rule mito_de:
    input:
        "data/scesets/{id}_sceset_raw.rds"
    output:
        "data/mito_differential_expression/{id}_mito_de.csv"
    shell:
        "Rscript pipeline/all_sample_mito/all_sample_mito_de.R --sce_input_file {input} --output_csv {output}"
    


    





