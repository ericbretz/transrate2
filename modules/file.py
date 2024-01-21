import pysam


def file(bamfile, bamdct, threads):
    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        references = bam.references

        for ref in references:
            bamdct[ref] = {
                'stats': {
                    'name'          : ref,
                    'p_seqtrue'     : 0.0,
                    'bridges'       : 0,
                    'length'        : bam.get_reference_length(ref),
                    'fragments'     : 0,
                    'both_mapped'   : 0,
                    'properpair'    : 0,
                    'good'          : 0,
                    'basesuncovered': 0,
                    'p_notsegmented': 0.0
                },
                'seq': {
                    'true'          : 0,
                    'count'         : 0
                },
                'good': {
                    'fragments'     : [],
                    'mean'          : 0.0,
                    'sd'            : 0.0,
                    'realistic'     : 0.0,
                },
                'coverage': {
                    'ranges_init'   : [],
                    'ranges_fnl'    : [],
                    'bases_count'   : 0,
                    'bases'         : [],
                    'length'        : 0,
                    'covered'       : [0] * bam.get_reference_length(ref)
                }
            }
    return bamdct    