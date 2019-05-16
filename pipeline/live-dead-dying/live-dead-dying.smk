



de_results = {
    'cluster': expand('data/live-dead-dying/ldd-cluster-de-{cv}.rds', cv=cellranger_versions),
    'low-mito': expand('data/live-dead-dying/ldd-lowmito-de-{cv}.rds', cv=cellranger_versions)
    }



deliverables['tmp-ldd'] = expand('data/live-dead-dying/ldd-intermediate-results-{cv}.rds',
                                 cv=cellranger_versions) + \
                                 de_results['cluster'] + de_results['low-mito']

rule live_dead_dying_analysis:
    params:
        curr_dir = os.getcwd()
    output:
        rds="data/live-dead-dying/ldd-intermediate-results-{cv}.rds",
        oreport="reports/live_dying_dead/{cv}_live_dying_dead-analysis.html"
    shell:
        "Rscript -e \"rmarkdown::render('{params.curr_dir}/pipeline/live-dead-dying/live-dead-dying2.Rmd', \
        output_file='{params.curr_dir}/{output.oreport}', \
        knit_root_dir='{params.curr_dir}', \
        params=list(cellranger_version='{wildcards.cv}', \
        output_rds='{output.rds}'))\" "

rule cluster_de:
    input:
        "data/live-dead-dying/ldd-intermediate-results-{cv}.rds"
    output:
        'data/live-dead-dying/ldd-cluster-de-{cv}.rds'
    shell:
        "Rscript pipeline/live-dead-dying/cluster-de.R \
        --input_rds {input} --output_rds {output}"

rule lowmito_de:
    input:
        "data/live-dead-dying/ldd-intermediate-results-{cv}.rds"
    output:
        'data/live-dead-dying/ldd-lowmito-de-{cv}.rds'
    shell:
        "Rscript pipeline/live-dead-dying/low-mito-de.R \
        --input_rds {input} --output_rds {output}"

all_figs['live-dead-dying'] = expand("figs/live-dead-dying/ldd_{cv}.png",
                                     cv=cellranger_versions)

deliverables['live-dead-dying-de-results'] = ['data/deliverables/ldd-cluster-de-v3.csv','data/deliverables/ldd-lowmito-de-v3.csv']

rule collate_analysis:
    params:
        curr_dir=os.getcwd()
    input:
        cluster_de='data/live-dead-dying/ldd-cluster-de-{cv}.rds',
        lowmito_de='data/live-dead-dying/ldd-lowmito-de-{cv}.rds',
        intermediate_rds='data/live-dead-dying/ldd-intermediate-results-{cv}.rds'
    output:
        png='figs/live-dead-dying/ldd_{cv}.png',
        stats='data/statistics/live_dead_dying_{cv}.csv',
        cluster_de_results='data/deliverables/ldd-cluster-de-{cv}.csv',
        lowmito_de_results='data/deliverables/ldd-lowmito-de-{cv}.csv'
    shell:
        "Rscript -e \"rmarkdown::render('{params.curr_dir}/pipeline/live-dead-dying/live-dead-dying-collate-for-paper2.Rmd', \
        knit_root_dir='{params.curr_dir}', \
        params=list(cellranger_version='{wildcards.cv}', \
        output_png='{output.png}',\
        intermediate_rds='{input.intermediate_rds}',\
        lowmito_rds='{input.lowmito_de}',\
        cluster_rds='{input.cluster_de}',\
        output_stats='{output.stats}',\
        cluster_de_results='{output.cluster_de_results}',\
        lowmito_de_results='{output.lowmito_de_results}'))\" "


        
# rule collate_analysis:
#     input:
#         "figs/live-dead-dying/ldd_{cv}.rds",
#     output:
#         png="figs/live-dead-dying/ldd_{cv}.png",
#         csv="data/deliverables/live-dead-dying_de-{cv}.csv",
# 	stats="data/statistics/live_dead_dying_{cv}.csv"
#     shell:
#         "Rscript pipeline/live-dead-dying/live-dead-dying-collate-for-paper.R \
#         --results {input} \
# 	--output_stats {output.stats} \
#         --output_png {output.png} \
#         --output_csv {output.csv}"
        
        
    
