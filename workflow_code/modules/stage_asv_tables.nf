include { shellQuote; cleanedPrefix; assaySuffix } from './utils.nf'

process STAGE_ASV_TABLES {
    input:
    path(counts, stageAs: 'source/counts.tsv')
    path(taxonomy, stageAs: 'source/taxonomy.tsv')
    path(runsheet)
    output:
    path({ "${cleanedPrefix()}counts${assaySuffix()}.tsv" }), emit: counts
    path({ "${cleanedPrefix()}taxonomy${assaySuffix()}.tsv" }), emit: taxonomy
    script:
    """
    validate_asv_tables.py --counts ${shellQuote(counts)} --taxonomy ${shellQuote(taxonomy)} \\
        --runsheet ${shellQuote(runsheet)} --out-counts ${shellQuote(cleanedPrefix() + 'counts' + assaySuffix() + '.tsv')} \\
        --out-taxonomy ${shellQuote(cleanedPrefix() + 'taxonomy' + assaySuffix() + '.tsv')}
    """
}
