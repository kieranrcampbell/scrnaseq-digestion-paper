

comparisons = ['collagenase_2hvs30m', '2hr', '30min', 'coldprotease_2hvs30m']
time_des = expand('data/time/time-de-{comparison}.rds',
                              comparison=comparisons)

all_figs['time'] = ['figs/final/time.png']


rule parse_time_sce:
    input:
        sces_qc,
    output:
        'data/time/time-sce.rds'
    shell:
        'Rscript pipeline/time/parse-time-sce.R '
        '--output_sce {output}'

rule time_de:
    input:
        'data/time/time-sce.rds',
    output:
        'data/time/time-de-{comparison}.rds',
    shell:
        'Rscript pipeline/time/time-de.R '
        '--input_sce {input} '
        '--comparison {wildcards.comparison} '
        '--output_results {output} '

rule time_figs_stats:
    params:
        curr_dir=os.getcwd(),
    input:
        time_des,
        sce='data/time/time-sce.rds'
    output:
        png='figs/final/time.png',
        stats='data/statistics/time.csv',
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/time/explore-time.Rmd', \
        knit_root_dir='{params.curr_dir}',\
        params=list(output_fig='{output.png}',\
        output_stats='{output.stats}'))\" "
        
