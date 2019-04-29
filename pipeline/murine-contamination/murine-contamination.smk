
# From Nick's azure output
human_mouse_input_sces = expand(config['human_mouse_dir'] + "/{hm}_{id}.rdata",
                                hm = ['human','mouse'],
                                id = ids)

all_figs['murine-contamination'] = \
    expand('figs/murine-contamination/{p}.png',
           p=['human_mouse_boxplot','human_mouse_histogram_pdxonly','pct_cells_mouse'])

rule identify_human_mouse:
    params:
        input_dir=config['human_mouse_dir'],
        curr_dir = os.getcwd()
    input:
        human_mouse_input_sces
    output:
        all_figs['murine-contamination'],
        csv=config['murine_contamination_csv']
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/murine-contamination/identify-mouse-cells.Rmd', \
        knit_root_dir='{params.curr_dir}', \
        params=list(human_mouse_dir='{params.input_dir}',\
        human_mouse_prop_file='{output.csv}'))\" "
