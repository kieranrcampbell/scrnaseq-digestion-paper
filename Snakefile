import pandas as pd
import os


configfile: "config.yaml"
files10X = ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']


df_inventory = pd.read_csv(config['sample_inventory_url']).dropna()
df_inventory.index = df_inventory['id']
inventory_dict = df_inventory.to_dict('index')

ids = list(inventory_dict.keys())

raw_files = expand("data/raw_10X/{id}/{file}", id=ids, file=files10X)
scesets_raw = expand("data/scesets/{id}_sceset_raw.rds", id=ids)
scesets_qc = expand("data/scesets/{id}_sceset_qc.rds", id=ids)

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
        sample_id=lambda wildcards: inventory_dict[wildcards.id]['sample_id'],
        batch_id=lambda wildcards: inventory_dict[wildcards.id]['batch_id'],
        sample_type=lambda wildcards: inventory_dict[wildcards.id]['sample_type'],
        cancer_type=lambda wildcards: inventory_dict[wildcards.id]['cancer_type'],
        digestion_temperature=lambda wildcards: inventory_dict[wildcards.id]['digestion_temperature'],
        tissue_state=lambda wildcards: inventory_dict[wildcards.id]['tissue_state'],
        enzyme_mix=lambda wildcards: inventory_dict[wildcards.id]['enzyme_mix'],
        jira_ticket=lambda wildcards: inventory_dict[wildcards.id]['jira_ticket'],
        cell_status=lambda wildcards: inventory_dict[wildcards.id]['cell_status'],
        genome=lambda wildcards: inventory_dict[wildcards.id]['genome'],
	filter_total_features=lambda wildcards: inventory_dict[wildcards.id]['filter_total_features'],
	filter_pct_counts_mito=lambda wildcards: inventory_dict[wildcards.id]['filter_pct_counts_mito']
    input:
        expand("data/raw_10X/{{id}}/{file}", file=files10X)
    output:
        "data/scesets/{id}_sceset_raw.rds"
    shell:
        "Rscript pipeline/conversion_to_sceset/convert_to_sceset.R \
        --input_data_path data/raw_10X/{wildcards.id}/ \
        --output_scepath {output} \
        --id {wildcards.id} \
        --sample_id {params.sample_id} \
        --batch_id {params.batch_id} \
        --sample_type {params.sample_type} \
        --cancer_type {params.cancer_type} \
        --digestion_temperature {params.digestion_temperature} \
        --tissue_state {params.tissue_state} \
        --enzyme_mix {params.enzyme_mix} \
        --jira_ticket '{params.jira_ticket}' \
        --cell_status {params.cell_status} \
        --genome {params.genome} \
	--filter_total_features {params.filter_total_features} \
	--filter_pct_counts_mito {params.filter_pct_counts_mito}"

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
        
        
    





