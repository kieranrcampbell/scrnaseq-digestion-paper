

# From Nick's azure output
human_mouse_input_sces = expand(config['human_mouse_dir'] + "/{hm}_{id}.rdata",
                                hm = ['human','mouse'],
                                id = ids)




all_figs['murine-contamination'] = \
    expand('figs/murine-contamination/{p}.png',
           p=['pct_cells_mouse']) + \
	   ['figs/murine-contamination/murine-qc.png']



rule identify_human_mouse:
    params:
        input_dir=config['human_mouse_dir'],
        curr_dir = os.getcwd()
    input:
        human_mouse_input_sces
    output:
        all_figs['murine-contamination'][:-1],
        csv=config['murine_contamination_csv'],
        stats="data/statistics/murine_fp.csv"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/murine-contamination/identify-mouse-cells.Rmd', \
        knit_root_dir='{params.curr_dir}', \
        params=list(human_mouse_dir='{params.input_dir}',\
        stats='{output.stats}',\
        human_mouse_prop_file='{output.csv}'))\" "

rule collate_murine_stats:
    params:
        curr_dir = os.getcwd()
    input:
        sces_raw,
        config['murine_contamination_csv']
    output:
        fig='figs/murine-contamination/murine-qc.png',
	stat=statistics['murine_cell_count']
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/murine-contamination/mouse-qc-metrics.Rmd', \
        knit_root_dir='{params.curr_dir}', \
        params=list(mouse_qc_fig='{output.fig}'))\" "


