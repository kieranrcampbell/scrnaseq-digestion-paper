
# Primary tumour fig
pt_fig = "figs/all-sample-overview/primary-tumour-cellprops.png"

# Final fig to be passed off
all_figs['all-sample-overview'] = expand("figs/all-sample-overview/all_sample_overview-{cv}.png",
                       cv=cellranger_versions) + \
                       [pt_fig]


umap_fig_rds = expand("figs/all-sample-overview/umap_all_{cv}.rds",
                      cv=cellranger_versions)
umap_fig_png = expand("figs/all-sample-overview/umap_all_{cv}.png",
                      cv=cellranger_versions)

pct_mito_fig = expand("figs/all-sample-overview/pct_counts_mito_all_samples-{cv}.png",
                      cv=cellranger_versions)


rule umap:
    input:
        sces_qc,
        config['murine_contamination_csv']
    output:
        rds="figs/all-sample-overview/umap_all_{cv}.rds",
        png="figs/all-sample-overview/umap_all_{cv}.png"
    shell:
        "Rscript pipeline/all-sample-overview/umap-all-samples.R \
        --cellranger_version {wildcards.cv} \
        --output_png {output.png} \
        --output_rds {output.rds} "

rule overview:
    params:
        curr_dir = os.getcwd()
    input:
        sces_qc,
        umap_fig="figs/all-sample-overview/umap_all_{cv}.rds"
    output:
        figure="figs/all-sample-overview/all_sample_overview-{cv}.png",
        pct_mito_fig="figs/all-sample-overview/pct_counts_mito_all_samples-{cv}.png",
        report="reports/all-sample-overview/all_sample_overview-{cv}.html",
        stats="data/statistics/all_cells_{cv}.csv"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/all-sample-overview/metric_overview.Rmd',\
        output_file='{params.curr_dir}/{output.report}', \
        knit_root_dir='{params.curr_dir}',\
        params=list(cellranger_version='{wildcards.cv}', \
        input_umap_rds='{input.umap_fig}', \
        output_figure='{output.figure}', \
        stats_file='{output.stats}', \
        pct_mito_fig='{output.pct_mito_fig}'))\" "


rule primary_tumour_fig:
    input:
        "data/primary_tumour_analysis/v5/sce_final_annotated/v3.rds"
    output:
        pt_fig
    shell:
        "Rscript pipeline/all-sample-overview/primary-tumour-plot.R \
        --input_sce {input} \
        --output_fig {output}"
