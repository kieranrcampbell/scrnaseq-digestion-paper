
# SingleCellExperiments
sces_raw = expand("data/scesets/{cv}/{id}_sceset_{cv}_raw.rds",
                  id=ids, cv = cellranger_versions)
sces_qc = expand("data/scesets/{cv}/{id}_sceset_{cv}_qc.rds",
                  id=ids, cv = cellranger_versions)

# TODO: move this to separate file
mito_des = expand("data/mito_differential_expression/{id}_mito_de.csv",
                  id = ids)



rule parse_nick_output:
    params:
        metadata_json=lambda wildcards: json.dumps(inventory_dict[wildcards.id])
    input:
        "data/temperature_rdata_v3/{id}.rdata" # <- file from nick's pipeline
    output:
        "data/scesets/v3/{id}_sceset_v3_raw.rds"
    shell:
        "Rscript pipeline/data-preparation/convert_to_sceset.R \
        --input_sce {input} \
        --output_sce {output} \
        --metadata_json '{params.metadata_json}'"
        

rule qc_scesets:
    params:
        curr_dir = os.getcwd()
    input:
        "data/scesets/{cv}/{id}_sceset_{cv}_raw.rds"
    output:
        sce="data/scesets/{cv}/{id}_sceset_{cv}_qc.rds",
        report="reports/qc/{cv}/qc_report_{id}_{cv}.html"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/data-preparation/sce_qc.Rmd', \
        output_file='{params.curr_dir}/{output.report}',  \
        knit_root_dir='{params.curr_dir}', \
        params=list(input_sce_path='{input}', output_sce_path='{output.sce}'))\" "

rule mito_de:
    input:
        "data/scesets/{id}_sceset_raw.rds"
    output:
        "data/mito_differential_expression/{id}_mito_de.csv"
    shell:
        "Rscript pipeline/all_sample_mito/all_sample_mito_de.R --sce_input_file {input} --output_csv {output}"
    


    





