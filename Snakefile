import pandas as pd
import os
import json


configfile: "private_config.yaml"

files10X_v2 = ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']
files10X_v3 = ['barcodes.tsv', 'features.tsv', 'matrix.mtx']

df_inventory = pd.read_csv(config['sample_inventory_url']).dropna()

assert df_inventory.shape[0] == 42

df_inventory.index = df_inventory['id']
inventory_dict = df_inventory.to_dict('index')

ids = list(inventory_dict.keys())

cellranger_versions = ['v2','v3']

cellranger_outputs_v2 = expand("data/cellranger-v2-outputs/{id}/outs/filtered_gene_bc_matrices/{g}/{f}",
                               g=config['genome'], f=files10X_v2, id=ids)

cellranger_outputs_v3 = expand("data/cellranger-v3-outputs/{id}/outs/filtered_feature_bc_matrix/{f}",
                               f=files10X_v3, id=ids)

cellranger_outputs = cellranger_outputs_v2 + cellranger_outputs_v3

# SingleCellExperiments
sces_raw = expand("data/scesets/{cv}/{id}_sceset_{cv}_raw.rds",
                  id=ids, cv = cellranger_versions)



# raw_files_blob = expand("data/raw_scrna_blob/{id}/{file}", id=ids, file=files10X)
# scesets_raw = expand("data/scesets/{id}_sceset_raw.rds", id=ids)
# scesets_qc = expand("data/scesets/{id}_sceset_qc.rds", id=ids)

mito_des = expand("data/mito_differential_expression/{id}_mito_de.csv",
id = ids)


rule all:
    input:
        cellranger_outputs, sces_raw

rule download_from_blob:
    params:
        base_url = config['blob_base_url']
    output:
        "data/cellranger-{cv}-outputs/{id}.tar.gz"
    shell:
        "wget -P data/cellranger-{wildcards.cv}-outputs/ {params.base_url}{wildcards.cv}/{wildcards.id}.tar.gz"

rule extract_cellranger_v2:
    input:
        "data/cellranger-v2-outputs/{id}.tar.gz"
    output:
        expand("data/cellranger-v2-outputs/{{id}}/outs/filtered_gene_bc_matrices/{g}/{f}",
               g = config['genome'], f = files10X_v2)
    shell:
        # Snakemake has already created the directory!
        "tar -C data/cellranger-v2-outputs/{wildcards.id} -xzf {input}"

rule extract_cellranger_v3:
    input:
        "data/cellranger-v3-outputs/{id}.tar.gz"
    output:
        expand("data/cellranger-v3-outputs/{{id}}/outs/filtered_feature_bc_matrix/{f}",
               f = files10X_v3)
    shell:
        "tar -C data/cellranger-v3-outputs/{wildcards.id} -xzf {input} && \
        gunzip data/cellranger-v3-outputs/{wildcards.id}/outs/filtered_feature_bc_matrix/*"



rule convert_to_sce_v2:
    params:
        metadata_json=lambda wildcards: json.dumps(inventory_dict[wildcards.id]),
        genome=config['genome']
    input:
        expand("data/cellranger-v2-outputs/{{id}}/outs/filtered_gene_bc_matrices/{g}/{f}",
               g = config['genome'], f = files10X_v2)
    output:
        "data/scesets/v2/{id}_sceset_v2_raw.rds"
    shell:
        "Rscript pipeline/conversion_to_sceset/convert_to_sceset.R \
        --input_data_path data/cellranger-v2-outputs/{wildcards.id}/outs/filtered_gene_bc_matrices/{params.genome}/ \
        --output_scepath {output} \
        --metadata_json '{params.metadata_json}'"

rule convert_to_sce_v3:
    params:
        metadata_json=lambda wildcards: json.dumps(inventory_dict[wildcards.id])
    input:
        expand("data/cellranger-v2-outputs/{{id}}/outs/filtered_gene_bc_matrices/{g}/{f}",
               g = config['genome'], f = files10X_v2)
    output:
        "data/scesets/v3/{id}_sceset_v3_raw.rds"
    shell:
        "Rscript pipeline/conversion_to_sceset/convert_to_sceset.R \
        --input_data_path data/cellranger-v3-outputs/{wildcards.id}/outs/filtered_feature_bc_matrix/ \
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
    


    





