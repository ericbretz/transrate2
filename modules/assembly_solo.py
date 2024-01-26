import pandas as pd
import os

def assembly_solo(fa_file, outdir, multi=False):
    column_headers = ['assembly', 'n_seqs', 'smallest', 'largest', 'n_bases', 'mean_len',
                    'median_len', 'std_len', 'n_under_200', 'n_over_1k', 'n_over_10k', 
                    'n_with_orf', 'mean_orf_percent', 'n90', 'n70', 'n50', 'n30', 'n10',
                    'gc', 'bases_n', 'proportion_n']
    
    output_df  = pd.DataFrame(columns=column_headers)
    file_path  = fa_file
    file_name  = os.path.basename(file_path)

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

    assembly             = file_name

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
   
    outfile    = os.path.join(outdir, 'assembly.csv')
    if os.path.exists(outfile) and multi:
        output_df.to_csv(outfile, index=False, mode='a', header=False)
    else:
        output_df.to_csv(outfile, index=False)
