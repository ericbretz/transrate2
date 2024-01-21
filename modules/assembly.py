import pandas as pd
import os
import math
import numpy as np

def assembly(r_csv, fa_file, snapf, outdir, star):
    column_headers = ['assembly', 'n_seqs', 'smallest', 'largest', 'n_bases', 'mean_len',
        'median_len', 'std_len', 'n_under_200', 'n_over_1k', 'n_over_10k', 'n_with_orf',
        'mean_orf_percent', 'n90', 'n70', 'n50', 'n30', 'n10', 'gc', 'bases_n',
        'proportion_n', 'fragments', 'fragments_mapped', 'p_fragments_mapped',
        'good_mappings', 'p_good_mapping', 'bad_mappings', 'potential_bridges',
        'bases_uncovered', 'p_bases_uncovered', 'contigs_uncovbase',
        'p_contigs_uncovbase', 'contigs_uncovered', 'p_contigs_uncovered',
        'contigs_lowcovered', 'p_contigs_lowcovered', 'contigs_segmented',
        'p_contigs_segmented', 'score', 'optimal_score', 'cutoff', 'weighted']
    
    output_df  = pd.DataFrame(columns=column_headers)
    file_path  = fa_file
    file_name  = os.path.basename(file_path)

    tc = pd.read_csv(r_csv)

    contig_df               = tc.copy()
    contig_df['eff_length'] = 0
    contig_df['eff_count']  = 0
    contig_df['tpm']        = 0
    contig_df['coverage']   = 0
    contig_df['sCnuc']      = 0
    contig_df['sCcov']      = 0
    contig_df['sCord']      = 0
    contig_df['sCseg']      = 0

    def sc_metrics():
        for index, row in contig_df.iterrows():
            contig_df.loc[index, 'sCnuc'] = row['p_seqtrue']
            contig_df.loc[index, 'sCcov'] = round((row['length'] - row['basesuncovered']) / row['length'], 4)
            contig_df.loc[index, 'sCord'] = row['good']
            contig_df.loc[index, 'sCseg'] = row['basesuncovered']
            contig_df.loc[index, 'coverage'] = round(row['eff_count'] * row['length'] / row['eff_length'], 4)

    if star:
        for line in open(snapf, 'r'):
            if 'Number of input reads' in line:
                fragments = int(line.split('	')[1])
    else:
        with open(snapf, 'r') as file:
            lines = file.readlines()

        fragments = int(int(lines[3].strip().split(' ')[0].replace(',', '')) / 2)

    def longest_orf(seq):
        longest  = 0
        len_list = [0, 0, 0]
        sl       = len(seq)
        str      = seq
        for i in range(sl - 2):
            if str[i:i+3] == 'atg':
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 1
            elif str[i:i+3] in ['tag', 'taa', 'tga']:
                longest = max(longest, len_list[i % 3])
                len_list[i % 3] = -1
            else:
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 0
        longest  = max(longest, max(len_list))
        len_list = [0, 0, 0]
        for i in range(sl - 1, 1, -1):
            if str[i-2:i+1] == 'cat':
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 1
            elif str[i-2:i+1] in ['cta', 'tta', 'tct']:
                longest = max(longest, len_list[i % 3])
                len_list[i % 3] = -1
            else:
                len_list[i % 3] = len_list[i % 3] + 1 if len_list[i % 3] >= 0 else 0
        return max(longest, max(len_list))

    frame       = ''
    bases       = 0
    gc          = 0
    read_length = 0
    lengths     = []
    first       = True
    seqs        = 0
    n_count     = 0
    n_prop      = 0
    n_under_200 = 0
    n_over_1k   = 0
    n_over_10k  = 0
    n_with_orf  = 0
    orf_len_sum = 0

    with open(file_path, 'r') as file:
        for line in file:
            line = line.rstrip('\n').lower()
            if line.startswith('>'):
                seqs += 1
                if not first:
                    lengths.append(read_length)
                    n_under_200 += read_length < 200
                    n_over_1k += read_length > 1000
                    n_over_10k += read_length > 10000
                    orf = longest_orf(frame)
                    n_with_orf += orf > 149
                    orf_len_sum += orf
                    read_length = 0
                    frame = ''
                first = False
                continue
            read_length += len(line)
            bases += len(line)
            gc += line.count('g') + line.count('c')
            n_count += line.count('n')
            frame += line
        lengths.append(read_length)
        n_under_200 += read_length < 200
        n_over_1k += read_length > 1000
        n_over_10k += read_length > 10000
        orf = longest_orf(frame)
        n_with_orf += orf > 149
        orf_len_sum += orf
        read_length = 0
        frame = ''

    n_prop         = round(n_count/bases, 2)
    gc             = round(gc/bases, 2)
    smallest       = min(lengths)
    largest        = max(lengths)
    mean           = round(sum(lengths)/len(lengths), 2)
    median         = lengths[int(len(lengths)/2)]
    std            = round((sum([(x - mean)**2 for x in lengths]) / len(lengths))**0.5, 2)
    lengths_sorted = sorted(lengths)
    base_90        = int(bases * 0.9)
    base_70        = int(bases * 0.7)
    base_50        = int(bases * 0.5)
    base_30        = int(bases * 0.3)
    base_10        = int(bases * 0.1)
    n90            = 0
    n70            = 0
    n50            = 0
    n30            = 0
    n10            = 0
    b90            = False
    b70            = False
    b50            = False
    b30            = False
    b10            = False
    base_total     = 0

    mean_orf_percent = round((300 * orf_len_sum) / (seqs * mean), 2)

    for length in lengths_sorted[::-1]:
        base_total += length
        if base_total > base_90 and b90 == False:
            n90 = length
            b90 = True
        if base_total > base_70 and b70 == False:
            n70 = length
            b70 = True
        if base_total > base_50 and b50 == False:
            n50 = length
            b50 = True
        if base_total > base_30 and b30 == False:
            n30 = length
            b30 = True
        if base_total > base_10 and b10 == False:
            n10 = length
            b10 = True

    salmonf = pd.read_csv(f'{outdir}/salmon/quant.sf', sep='\t')
    contigs_uncovered  = 0
    contigs_lowcovered = 0
    tc['tpm'] = 0

    for index, row in salmonf.iterrows():
        if pd.notnull(row['EffectiveLength']):
            salmon_efflen = row['EffectiveLength']
        else:
            salmon_efflen = 0

        if pd.notnull(row['NumReads']):
            salmon_count  = row['NumReads']
        else:
            salmon_count  = 0

        if pd.notnull(row['TPM']):
            salmon_tpm    = row['TPM']
        else:
            salmon_tpm    = 0

        if pd.notnull(row['Length']):
            salmon_len    = row['Length']
        else:
            salmon_len    = 0

        if salmon_len    == 0:
            coverage      = 0
        else:
            coverage      = salmon_count * 150 / salmon_len

        contig_df.loc[contig_df['name'] == row['Name'], 'eff_length'] = salmon_efflen
        contig_df.loc[contig_df['name'] == row['Name'], 'eff_count']  = salmon_count
        contig_df.loc[contig_df['name'] == row['Name'], 'tpm']        = float(salmon_tpm)
        contigs_uncovered  += 1 if coverage < 1 else 0
        contigs_lowcovered += 1 if coverage < 10 else 0
    
    sc_metrics()

    assembly             = file_name
    fragments_mapped     = tc['fragments'].sum()
    p_fragments_mapped   = round(fragments_mapped / fragments, 2)
    good_mappings        = tc['good'].sum()
    p_good_mappings      = round(good_mappings / fragments, 2)
    bad_mappings         = round(fragments_mapped - good_mappings, 2)
    potential_bridges    = (tc['bridges'] > 0).sum()
    bases_uncovered      = tc['basesuncovered'].sum()
    p_bases_uncovered    = round(bases_uncovered / bases, 2)
    p_contigs_uncovered  = round(contigs_uncovered / tc.shape[0], 2)
    p_contigs_lowcovered = round(contigs_lowcovered / tc.shape[0], 2)
    contigs_segmented    = round((tc['p_notsegmented'] < 0.5).sum(), 2)
    p_contigs_segmented  = round(contigs_segmented / tc.shape[0], 2)
    contig_uncovbase     = round((tc['basesuncovered'] > 0).sum(), 2)
    p_contig_uncovbase   = round(contig_uncovbase / tc.shape[0], 2)

    def score():
        scores = []
        for index, row in tc.iterrows():
            prod = np.max([(row['length'] - row['basesuncovered'])/row['length'], 0.01]) * \
                    np.max([row['p_notsegmented'], 0.01]) * \
                    np.max([row['good'] / row['fragments'] if row['fragments'] > 1 else 0.01, 0.01]) * \
                    np.max([row['p_seqtrue'], 0.01])
            s = np.max([prod, 0.01])
            scores.append(s)
        return scores

    def raw_score(scores, good, total):
        contig_score = geomean(scores)
        return contig_score * (good / total)
    
    def geomean(x):
        sum = 0.0
        for v in x:
            sum += math.log(v)
        sum /= len(x)
        return math.exp(sum)

    def optimal_score():
        opt_product   = 0
        opt_good      = 0
        opt_product   = contig_df['scores'].apply(math.log).sum()
        opt_good      = tc['good'].sum()
        opt_count     = len(tc)
        cutoffscores  = {}
        scores_sorted = contig_df[['scores', 'good']].sort_values(by='scores')
        for index, row in scores_sorted.iterrows():
            opt_product -= math.log(row['scores'])
            opt_good    -= row['good']
            opt_count   -= 1 if opt_count > 1 else 0
            new_score    = math.exp(opt_product / opt_count) * (opt_good / fragments)
            cutoffscores[row['scores']] = new_score
        optimal = 0
        cutoff  = 0
        for c, score in cutoffscores.items():
            if score > optimal:
                optimal = score
                cutoff = c
        cutoffscores_df = pd.DataFrame(list(cutoffscores.items()), columns=['cutoff', 'score'])
        cutoffscores_df = cutoffscores_df.round(5)
        cuttofffile     = os.path.join(outdir, 'transrate', 'assembly_score_optimisation.csv')
        cutoffscores_df.to_csv(cuttofffile, index=False)
        return [round(optimal,6), round(cutoff,6)]

    def weighted_score():
        wscores = []
        for index, row in contig_df.iterrows():
            wscores.append(row['scores'] * row['tpm'])
        contig_weighted = sum(wscores) / len(wscores)
        return round(contig_weighted * (good_mappings / fragments), 6)

    scores = score()
    contig_df['scores'] = scores
    raw_score_val = round(raw_score(scores, good_mappings, fragments),6)
    optimal_score_val = optimal_score()

    #### CSV OUTPUTS ####
    output_df.loc[0, 'assembly']             = assembly
    output_df.loc[0, 'n_seqs']               = seqs
    output_df.loc[0, 'smallest']             = smallest
    output_df.loc[0, 'largest']              = largest
    output_df.loc[0, 'n_bases']              = bases
    output_df.loc[0, 'mean_len']             = mean
    output_df.loc[0, 'median_len']           = median
    output_df.loc[0, 'std_len']              = std
    output_df.loc[0, 'n_under_200']          = n_under_200
    output_df.loc[0, 'n_over_1k']            = n_over_1k
    output_df.loc[0, 'n_over_10k']           = n_over_10k
    output_df.loc[0, 'n_with_orf']           = n_with_orf
    output_df.loc[0, 'mean_orf_percent']     = mean_orf_percent
    output_df.loc[0, 'n90']                  = n90
    output_df.loc[0, 'n70']                  = n70
    output_df.loc[0, 'n50']                  = n50
    output_df.loc[0, 'n30']                  = n30
    output_df.loc[0, 'n10']                  = n10
    output_df.loc[0, 'gc']                   = gc
    output_df.loc[0, 'bases_n']              = n_count
    output_df.loc[0, 'proportion_n']         = n_prop
    output_df.loc[0, 'fragments']            = fragments
    output_df.loc[0, 'fragments_mapped']     = fragments_mapped
    output_df.loc[0, 'p_fragments_mapped']   = p_fragments_mapped
    output_df.loc[0, 'good_mappings']        = good_mappings
    output_df.loc[0, 'p_good_mapping']       = p_good_mappings
    output_df.loc[0, 'bad_mappings']         = bad_mappings
    output_df.loc[0, 'potential_bridges']    = potential_bridges
    output_df.loc[0, 'bases_uncovered']      = bases_uncovered
    output_df.loc[0, 'p_bases_uncovered']    = p_bases_uncovered
    output_df.loc[0, 'contigs_uncovbase']    = contig_uncovbase
    output_df.loc[0, 'p_contigs_uncovbase']  = p_contig_uncovbase
    output_df.loc[0, 'contigs_uncovered']    = contigs_uncovered
    output_df.loc[0, 'p_contigs_uncovered']  = p_contigs_uncovered
    output_df.loc[0, 'contigs_lowcovered']   = contigs_lowcovered
    output_df.loc[0, 'p_contigs_lowcovered'] = p_contigs_lowcovered
    output_df.loc[0, 'contigs_segmented']    = contigs_segmented
    output_df.loc[0, 'p_contigs_segmented']  = p_contigs_segmented
    output_df.loc[0, 'score']                = raw_score_val
    output_df.loc[0, 'optimal_score']        = optimal_score_val[0]
    output_df.loc[0, 'cutoff']               = optimal_score_val[1]
    output_df.loc[0, 'weighted']             = weighted_score()

    outfile    = os.path.join(outdir, 'transrate', 'assembly.csv')
    output_df.to_csv(outfile, index=False)

    contigfile = os.path.join(outdir, 'transrate', 'contigs.csv')
    contig_df.to_csv(contigfile, index=False)
