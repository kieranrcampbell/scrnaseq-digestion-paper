import pandas as pd
import os
import json

# Config options
configfile: "private_config.yaml"



files10X = ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']


df_inventory = pd.read_csv(config['sample_inventory_url']).dropna()

assert df_inventory.shape[0] == 42

df_inventory.index = df_inventory['id']
inventory_dict = df_inventory.to_dict('index')

ids = list(inventory_dict.keys())

raw_10X_files = expand(os.path.join(config['raw_path_local'], "{id}/outs/filtered_feature_bc_matrix/{f}"), id=ids, f=files10X)

scesets_raw = expand("data/scesets/{id}_sceset_raw.rds", id=ids)
scesets_qc = expand("data/scesets/{id}_sceset_qc.rds", id=ids)

mito_des = expand("data/mito_differential_expression/{id}_mito_de.csv",
id = ids)


rule all:
    input:
        #raw_10X_files
        scesets_raw
        #scesets_qc

rule get_from_blob:
    params:
        blob_url = config['blob_url'],
        intermediate_dir = lambda wildcards: os.path.join(config['raw_path_local'], wildcards.id)
    output:
        expand("{c}/{{id}}/outs/filtered_feature_bc_matrix/{f}", f=files10X, c = config['raw_path_local'])
    shell:
        "wget -P {params.intermediate_dir} {params.blob_url}/{wildcards.id}.tar.gz && \
        tar -C {params.intermediate_dir} -xzf {params.intermediate_dir}/{wildcards.id}.tar.gz && \
        gunzip {params.intermediate_dir}/outs/filtered_feature_bc_matrix/* && \
        mv {params.intermediate_dir}/outs/filtered_feature_bc_matrix/features.tsv {params.intermediate_dir}/outs/filtered_feature_bc_matrix/genes.tsv"
        
rule convert_to_sce:
    params:
        metadata_json=lambda wildcards: json.dumps(inventory_dict[wildcards.id]),
        rp = config['raw_path_local']
    input:
        expand("{c}/{{id}}/outs/filtered_feature_bc_matrix/{f}", f=files10X, c = config['raw_path_local'])
    output:
        "data/scesets/{id}_sceset_raw.rds"
    shell:
        "Rscript pipeline/conversion_to_sceset/convert_to_sceset.R \
        --input_data_path {params.rp}/{wildcards.id}/outs/filtered_feature_bc_matrix \
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
    


    





