from pysam import FastaFile

class IterContig:
    def __init__(self):
        return
    
    def mainRun(self, args):
        main    = args[0]
        i       = args[1]
        fasta   = FastaFile(main['assembly'])
        refs    = fasta.references[i::main['threads']]
        iterDct = {'assembly': {
                            'fragmentsMapped'        : 0,
                            'basesUncovered'         : 0,
                            'contigsUncovBase'       : 0,
        }}
        if main['mode'] == 2:
            iterDct['assembly'].update({
                            'goodMappings'           : 0,
                            'badMappings'            : 0,
                            'potentialBridges'       : 0,
                            'contigsSegmented'       : 0,
            })
        if main['mode'] >= 1:
            iterDct['assembly'].update({
                            'softclipped'            : 0,
                            'softclippedLength'      : 0,
            })

        filtered_df = main['contigDF'][main['contigDF']['name'].isin(refs)]

        for ref in refs:
            if ref == 'assembly':
                continue
            ref_data = filtered_df[filtered_df['name'] == ref]
            if not ref_data.empty:
                iterDct['assembly']['fragmentsMapped']  += ref_data['fragments'].iloc[0]
                iterDct['assembly']['basesUncovered']   += ref_data['basesUncovered'].iloc[0]
                iterDct['assembly']['contigsUncovBase'] += 1 if ref_data['basesUncovered'].iloc[0] > 0 else 0
                iterDct['assembly']['softclipped']      += 1 if ref_data['softclipped'].iloc[0] else 0

                if main['mode'] == 2:
                    iterDct['assembly']['goodMappings']     += ref_data['good'].iloc[0]
                    iterDct['assembly']['badMappings']      += ref_data['fragments'].iloc[0] - ref_data['good'].iloc[0]
                    iterDct['assembly']['potentialBridges'] += ref_data['bridges'].iloc[0]
                    iterDct['assembly']['contigsSegmented'] += 1 if ref_data['pNotSegmented'].iloc[0] < 0.5 else 0

        return iterDct
        