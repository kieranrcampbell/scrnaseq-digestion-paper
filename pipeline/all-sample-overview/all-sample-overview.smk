
# Final fig to be passed off
all_figs['all-sample-overview'] = expand("figs/all-sample-overview/all_sample_overview-{cv}.png",
                       cv=cellranger_versions)


umap_fig_rds = expand("figs/all-sample-overview/umap_all_{cv}.rds",
                      cv=cellranger_versions)
umap_fig_png = expand("figs/all-sample-overview/umap_all_{cv}.png",
                      cv=cellranger_versions)

pct_mito_fig = expand("figs/all-sample-overview/pct_counts_mito_all_samples-{cv}.png",
                      cv=cellranger_versions)


rule umap:
    input:
        sces_qc
    output:
        rds="figs/all-sample-overview/umap_all_{cv}.rds",
        png="figs/all-sample-overview/umap_all_{cv}.png"
    shell:
        "Rscript pipeline/all-sample-overview/umap-all-samples.R \
        --cellranger_version {wildcards.cv} \
        --output_png {output.png} \
        --output_rds {output.rds} "

rule overview:
    input:
        sces_qc,
        umap_fig="figs/all-sample-overview/umap_all_{cv}.rds"
    output:
        figure="figs/all-sample-overview/all_sample_overview-{cv}.png",
        pct_mito_fig="figs/all-sample-overview/pct_counts_mito_all_samples-{cv}.png",
        report="reports/all-sample-overview/all_sample_overview-{cv}.html"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/all-sample-overview/metric_overview.Rmd',\
        output_file='{output.report}', \
        params=list(cellranger_version='{wildcards.cv}', \
        input_umap_rds='{input.umap_fig}', \
        output_figure='{output.figure}', \
        pct_mito_fig='{output.pct_mito_fig}'))\" "
