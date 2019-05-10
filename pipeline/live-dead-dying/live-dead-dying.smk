
all_figs['live-dead-dying'] = expand("figs/live-dead-dying/ldd_{cv}.png",
                                     cv=cellranger_versions)

deliverables['live-dead-dying-de-results'] = expand("data/deliverables/live_dead_dying_de-{cv}.csv",
                                                    cv=cellranger_versions)


rule live_dead_dying_analysis:
    params:
        curr_dir = os.getcwd()
    output:
        rds="figs/live-dead-dying/ldd_{cv}.rds",
        report="reports/live_dying_dead/{cv}_live_dying_dead.html"
    shell:
        "Rscript -e \"rmarkdown::render('{params.curr_dir}/pipeline/live-dead-dying/live-dead-dying.Rmd', \
        output_file='{params.curr_dir}/{output.report}', \
        knit_root_dir='{params.curr_dir}', \
        params=list(cellranger_version='{wildcards.cv}', \
        output_rds='{output.rds}'))\" "

rule collate_analysis:
    input:
        "figs/live-dead-dying/ldd_{cv}.rds",
    output:
        png="figs/live-dead-dying/ldd_{cv}.png",
        csv="data/deliverables/live-dead-dying_de-{cv}.csv",
	stats="data/statistics/live_dead_dying_{cv}.csv"
    shell:
        "Rscript pipeline/live-dead-dying/live-dead-dying-collate-for-paper.R \
        --results {input} \
	--output_stats {output.stats} \
        --output_png {output.png} \
        --output_csv {output.csv}"
        
        
    
