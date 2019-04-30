import os


cellranger_versions = ['v2','v3']

#base_output_path = "../../data/primary_tumour_temp_de/{cv}/{cell_type}_"
base_figure_path = "../../figs/live_dead_dying/{cv}/"

figures = ['pct_mito', 'pathway_enrichment']
extensions = ['png', 'rds']

figs = expand(base_figure_path + "{f}.{e}",
              cv=cellranger_versions,
              f=figures,
              e=extensions)

figs = expand("../../figs/live_dead_dying/ldd_{cv}.rds",
              cv=cellranger_versions)


rule all:
    input:
        figs



rule make_figs:
    params:
        curr_dir = os.getcwd()
    output:
        rds="../../figs/live_dead_dying/ldd_{cv}.rds",
        report="../../reports/live_dying_dead/{cv}_live_dying_dead.html"
    shell:
        "Rscript -e \"rmarkdown::render('{params.curr_dir}/live_dead_dying.Rmd', \
        output_file='{params.curr_dir}/{output.report}', \
        knit_root_dir='{params.curr_dir}', \
        params=list(cellranger_version='{wildcards.cv}', \
        output_rds='{output.rds}'))\" "

