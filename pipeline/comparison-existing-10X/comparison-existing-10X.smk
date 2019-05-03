

e_ids = config['existing-10X'] # ID-genome dictionary for the existing data

files10X = ['barcodes.tsv',
           'genes.tsv',
           'matrix.mtx']

intermediate_csvs = expand('data/comparison-existing-10X/{id}.csv',id=e_ids.keys())

all_figs['comparison-existing-10X'] = ['figs/comparison-existing-10X/comparison-existing-10X.png']

rule get_metrics:
    params:
        genome = lambda wildcards: e_ids[wildcards.id]
    input:
        lambda wildcards: expand("data/external/existing-10X/{{id}}/filtered_gene_bc_matrices/{genome}/{f}",
                                 genome=e_ids[wildcards.id],
                                 f=files10X)
    output:
        "data/comparison-existing-10X/{id}.csv"
    shell:
        "Rscript pipeline/comparison-existing-10X/comparison-existing-10X.R \
        --id {wildcards.id} \
        --genome {params.genome} \
        --output_csv {output}"

rule make_plot:
    input:
        intermediate_csvs
    output:
        all_figs['comparison-existing-10X']
    shell:
        "Rscript pipeline/comparison-existing-10X/plot-existing.R \
        --output_png {output}"
