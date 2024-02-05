import pandas as pd
import os
import math
import numpy as np
import sys
import time

class Assembly:
    def __init__(self):
        #### Setup ####
        self.CSV      = ''
        self.FA       = ''
        self.SNAP     = ''
        self.OUTDIR   = ''
        self.BASE     = ''
        self.HEADERS  = []
        self.CHEADERS = []
        self.MULTI    = False
        self.SINGLE   = False
        self.SOLO     = False
        self.STAR     = False

        #### Dataframes ####
        self.OUTPUT_DF = None
        self.CONTIG_DF = None

        #### Metrics ####
        self.fragments   = 0

        #### Iter / ORF ####
        self.frame       = ''
        self.bases       = 0
        self.gc          = 0
        self.read_length = 0
        self.lengths     = []
        self.first       = True
        self.seqs        = 0
        self.n_count     = 0
        self.n_prop      = 0
        self.n_under_200 = 0
        self.n_over_1k   = 0
        self.n_over_10k  = 0
        self.n_with_orf  = 0
        self.orf_len_sum = 0
        self.lengths     = []

        #### Assembly Stats ####
        self.n_prop         = 0
        self.gc             = 0
        self.smallest       = 0
        self.largest        = 0
        self.mean           = 0
        self.median         = 0
        self.std            = 0
        self.lengths_sorted = []
        self.base_90        = 0
        self.base_70        = 0
        self.base_50        = 0
        self.base_30        = 0
        self.base_10        = 0
        self.n90            = 0
        self.n70            = 0
        self.n50            = 0
        self.n30            = 0
        self.n10            = 0
        self.b90            = False
        self.b70            = False
        self.b50            = False
        self.b30            = False
        self.b10            = False
        self.base_total     = 0

        #### Contig Stats ####
        self.prop_gc_contig      = 0
        self.orf_length_contig   = 0
        self.p_good_contig       = 0
        self.p_bases_covered_contig = 0
        self.good_dct_contigs = {}

        #### Salmon Stats ####
        self.SALMONF            = ''
        self.contigs_uncovered  = 0
        self.contigs_lowcovered = 0

        #### Assembly Stats ####
        self.assembly             = ''
        self.fragments_mapped     = 0
        self.p_fragments_mapped   = 0
        self.good_mappings        = 0
        self.p_good_mappings      = 0
        self.bad_mappings         = 0
        self.potential_bridges    = 0
        self.bases_uncovered      = 0
        self.p_bases_uncovered    = 0
        self.p_contigs_uncovered  = 0
        self.p_contigs_lowcovered = 0
        self.contigs_segmented    = 0
        self.p_contigs_segmented  = 0
        self.contig_uncovbase     = 0
        self.p_contig_uncovbase   = 0
        self.score               = []
        self.goodcontig           = 0
        self.goodcontiglist       = []



    def run(self):
        self.BASE = os.path.basename(self.FA)
        self.headers()
        self.starsnap() if not self.SOLO else None
        self.file_iter()
        self.basestats()
        self.nstats()
        self.salmonstats() if not self.SOLO else None
        self.assemblystats() if not self.SOLO else None
        self.score_run() if not self.SOLO and not self.SINGLE else None
        self.sc_metrics() if not self.SOLO else None
        self.assembly_output()
        self.contig_output() if not self.SOLO else None
        # t = time.perf_counter()
        # self.headers()
        # print('Headers:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.starsnap() if not self.SOLO else None
        # print('Snap:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.file_iter()
        # print('File Iter:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.basestats()
        # print('Base Stats:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.nstats()
        # print('N Stats:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.salmonstats() if not self.SOLO else None
        # print('Salmon Stats:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.assemblystats() if not self.SOLO else None
        # print('Assembly Stats:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.score_run() if not self.SOLO and not self.SINGLE else None
        # print('Score Run:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.sc_metrics() if not self.SOLO else None
        # print('SC Metrics:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.assembly_output()
        # print('Assembly Output:', time.perf_counter() - t)
        # t = time.perf_counter()
        # self.contig_output() if not self.SOLO else None
        # print('Contig Output:', time.perf_counter() - t)
        # return

    def headers(self):
        if self.SINGLE:
            self.HEADERS = ['assembly', 'n_seqs', 'smallest', 'largest', 'n_bases', 'mean_len',
                            'median_len', 'std_len', 'n_under_200', 'n_over_1k', 'n_over_10k', 
                            'n_with_orf', 'mean_orf_percent', 'n90', 'n70', 'n50', 'n30', 'n10', 
                            'gc', 'bases_n', 'proportion_n', 'fragments', 'fragments_mapped', 
                            'p_fragments_mapped', 'bases_uncovered', 'p_bases_uncovered', 
                            'contigs_uncovbase', 'p_contigs_uncovbase', 'contigs_uncovered', 
                            'p_contigs_uncovered', 'contigs_lowcovered', 'p_contigs_lowcovered']
        elif self.SOLO:
            self.HEADERS = ['assembly', 'n_seqs', 'smallest', 'largest', 'n_bases', 'mean_len',
                            'median_len', 'std_len', 'n_under_200', 'n_over_1k', 'n_over_10k', 
                            'n_with_orf', 'mean_orf_percent', 'n90', 'n70', 'n50', 'n30', 'n10',
                            'gc', 'bases_n', 'proportion_n']
        else:
            self.HEADERS = ['assembly', 'n_seqs', 'smallest', 'largest', 'n_bases', 'mean_len',
                            'median_len', 'std_len', 'n_under_200', 'n_over_1k', 'n_over_10k', 
                            'n_with_orf', 'mean_orf_percent', 'n90', 'n70', 'n50', 'n30', 'n10', 
                            'gc', 'bases_n', 'proportion_n', 'fragments', 'fragments_mapped', 
                            'p_fragments_mapped', 'good_mappings', 'p_good_mapping', 'bad_mappings', 
                            'potential_bridges', 'bases_uncovered', 'p_bases_uncovered', 
                            'contigs_uncovbase', 'p_contigs_uncovbase', 'contigs_uncovered', 
                            'p_contigs_uncovered', 'contigs_lowcovered', 'p_contigs_lowcovered', 
                            'contigs_segmented', 'p_contigs_segmented', 'score', 'optimal_score', 
                            'cutoff', 'weighted', 'goodcontig', 'sCnuc', 'sCcov', 'sCord', 'sCseg']
            
        if not self.SOLO:
            self.CHEADERS = [
            'name', 'length', 'prop_gc', 'orf_length', 'bridges', 'properpair',
            'fragments', 'both_mapped', 'good', 'p_good', 'basesuncovered', 
            'p_bases_covered', 'p_seqtrue', 'score', 'p_notsegmented', 'eff_length', 
            'eff_count', 'tpm', 'coverage', 'sCnuc', 'sCcov', 'sCord', 'sCseg'
            ]
            self.CONTIG_DF = pd.read_csv(self.CSV, sep=',')
            missing_columns = list(set(self.CHEADERS) - set(self.CONTIG_DF.columns))
            for column in missing_columns:
                self.CONTIG_DF[column] = 0

        self.OUTPUT_DF = pd.DataFrame(columns=self.HEADERS)
            

    def starsnap(self):
        if self.STAR:
            for line in open(self.SNAP, 'r'):
                if 'Number of input reads' in line:
                    self.fragments = int(line.split('	')[1])
        else:
            with open(self.SNAP, 'r') as file:
                lines = file.readlines()
                self.fragments = int(int(lines[3].strip().split(' ')[0].replace(',', '')) / 2)
   
    def longest_orf(self):
        longest  = 0
        len_list = [0, 0, 0]
        sl       = len(self.frame)
        str      = self.frame
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
    
    def orf(self):
        self.lengths.append(self.read_length)
        self.n_under_200 += self.read_length < 200
        self.n_over_1k   += self.read_length > 1000
        self.n_over_10k  += self.read_length > 10000
        orf               = self.longest_orf()
        self.orf_length_contig = orf
        self.n_with_orf  += orf > 149
        self.orf_len_sum += orf
        self.read_length  = 0

    def file_iter(self):
        with open(self.FA, 'r') as file:
            name = ''
            for line in file:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    self.seqs += 1
                    if not self.first:
                        self.orf()
                        self.good_dct_contigs[name] = {'orf_len': 0, 'prop_gc': 0}
                        self.prop_gc_contig = round((self.frame.lower().count('g') + self.frame.lower().count('c')) / len(self.frame), 2)
                        self.good_dct_contigs[name]['prop_gc'] = self.prop_gc_contig
                        self.good_dct_contigs[name]['orf_len'] = self.orf_length_contig
                        self.frame        = ''
                    name = line.split(' ')[0].strip('>')
                    self.first = False
                    continue
                self.read_length += len(line)
                self.bases       += len(line)
                self.gc          += line.lower().count('g') + line.lower().count('c')
                self.n_count     += line.lower().count('n')
                self.frame       += line

            self.orf()
            self.good_dct_contigs[name] = {'orf_len': 0, 'prop_gc': 0}
            self.good_dct_contigs[name]['orf_len'] = self.orf_length_contig
            self.prop_gc_contig = round((self.frame.lower().count('g') + self.frame.lower().count('c')) / len(self.frame), 2)
            self.good_dct_contigs[name]['prop_gc'] = self.prop_gc_contig

    def basestats(self):
        self.n_prop           = round(self.n_count/self.bases, 2)
        self.gc               = round(self.gc/self.bases, 2)
        self.smallest         = min(self.lengths)
        self.largest          = max(self.lengths)
        self.mean             = round(sum(self.lengths)/len(self.lengths), 2)
        self.median           = self.lengths[int(len(self.lengths)/2)]
        self.std              = round((sum([(x - self.mean)**2 for x in self.lengths]) / len(self.lengths))**0.5, 2)
        self.lengths_sorted   = sorted(self.lengths)
        self.base_90          = int(self.bases * 0.9)
        self.base_70          = int(self.bases * 0.7)
        self.base_50          = int(self.bases * 0.5)
        self.base_30          = int(self.bases * 0.3)
        self.base_10          = int(self.bases * 0.1)
        self.mean_orf_percent = round((300 * self.orf_len_sum) / (self.seqs * self.mean), 2)

    def nstats(self):
        for length in self.lengths_sorted[::-1]:
            self.base_total += length
            if self.base_total > self.base_90 and self.b90 == False:
                self.n90 = length
                self.b90 = True
            if self.base_total > self.base_70 and self.b70 == False:
                self.n70 = length
                self.b70 = True
            if self.base_total > self.base_50 and self.b50 == False:
                self.n50 = length
                self.b50 = True
            if self.base_total > self.base_30 and self.b30 == False:
                self.n30 = length
                self.b30 = True
            if self.base_total > self.base_10 and self.b10 == False:
                self.n10 = length
                self.b10 = True

    def salmonstats(self):
        self.SALMONF                          = pd.read_csv(os.path.join(self.OUTDIR,'salmon','quant.sf'), sep='\t')
        self.CONTIG_DF['eff_length']          = self.SALMONF['EffectiveLength'].fillna(0)
        self.CONTIG_DF['eff_count']           = self.SALMONF['NumReads'].fillna(0)
        self.CONTIG_DF['tpm']                 = self.SALMONF['TPM'].fillna(0).astype(float)
        self.CONTIG_DF['eff_length']          = self.CONTIG_DF['eff_length'].where(self.CONTIG_DF['eff_length'] != 0, 0)
        self.CONTIG_DF['eff_count']           = self.CONTIG_DF['eff_count'].where(self.CONTIG_DF['eff_count'] != 0, 0)
        self.CONTIG_DF['tpm']                 = self.CONTIG_DF['tpm'].where(self.CONTIG_DF['tpm'] != 0, 0)
        self.CONTIG_DF['eff_length']          = self.CONTIG_DF['eff_length'].where(self.SALMONF['Length'] != 0, 0)
        coverage                              = self.SALMONF['NumReads'] * 150 / self.SALMONF['Length']
        self.contigs_uncovered               += (coverage < 1).sum()
        self.contigs_lowcovered              += (coverage < 10).sum()

    def assemblystats(self):
        self.assembly             = self.BASE
        self.fragments_mapped     = self.CONTIG_DF['fragments'].sum()
        self.p_fragments_mapped   = round(self.fragments_mapped / self.fragments, 2)
        self.bases_uncovered      = self.CONTIG_DF['basesuncovered'].sum()
        self.p_bases_uncovered    = round(self.bases_uncovered / self.bases, 2)
        self.p_contigs_uncovered  = round(self.contigs_uncovered / self.CONTIG_DF.shape[0], 2)
        self.p_contigs_lowcovered = round(self.contigs_lowcovered / self.CONTIG_DF.shape[0], 2)
        self.contig_uncovbase     = round((self.CONTIG_DF['basesuncovered'] > 0).sum(), 2)
        self.p_contig_uncovbase   = round(self.contig_uncovbase / self.CONTIG_DF.shape[0], 2)
        if not self.SINGLE:
            self.good_mappings        = self.CONTIG_DF['good'].sum()
            self.p_good_mappings      = round(self.good_mappings / self.fragments, 2)
            self.bad_mappings         = round(self.fragments_mapped - self.good_mappings, 2)
            self.potential_bridges    = (self.CONTIG_DF['bridges'] > 0).sum()
            self.contigs_segmented    = round((self.CONTIG_DF['p_notsegmented'] < 0.5).sum(), 2)
            self.p_contigs_segmented  = round(self.contigs_segmented / self.CONTIG_DF.shape[0], 2)
    
    def scores(self):
        score = []
        for index, row in self.CONTIG_DF.iterrows():
            prod = np.max([(row['length'] - row['basesuncovered'])/row['length'], 0.01]) * \
                    np.max([row['p_notsegmented'], 0.01]) * \
                    np.max([row['good'] / row['fragments'] if row['fragments'] > 1 else 0.01, 0.01]) * \
                    np.max([row['p_seqtrue'], 0.01])
            s = np.max([prod, 0.01])
            score.append(s)
        return score
    
    def raw_score(self):
        contig_score = self.geomean(self.score)
        return contig_score * (self.good_mappings / self.fragments)
    
    def geomean(self, x):
        sum = 0.0
        for v in x:
            sum += math.log(v)
        sum /= len(x)
        return math.exp(sum)
    
    def optimal_score(self):
        opt_product   = 0
        opt_good      = 0
        opt_product   = self.CONTIG_DF['score'].apply(math.log).sum()
        opt_good      = self.CONTIG_DF['good'].sum()
        opt_count     = len(self.CONTIG_DF)
        cutoffscore  = {}
        score_sorted = self.CONTIG_DF[['score', 'good']].sort_values(by='score')
        for index, row in score_sorted.iterrows():
            opt_product -= math.log(row['score'])
            opt_good    -= row['good']
            opt_count   -= 1 if opt_count > 1 else 0
            new_score    = math.exp(opt_product / opt_count) * (opt_good / self.fragments)
            cutoffscore[row['score']] = new_score
        optimal = 0
        cutoff  = 0
        for c, score in cutoffscore.items():
            if score > optimal:
                optimal = score
                cutoff = c
        cutoffscore_df = pd.DataFrame(list(cutoffscore.items()), columns=['cutoff', 'score'])
        cutoffscore_df = cutoffscore_df.round(5)
        scorefile      = 'assembly_score_optimisation.csv' if not self.MULTI else f'{self.BASE}_assembly_score_optimisation.csv'
        cuttofffile     = os.path.join(self.OUTDIR, 'transrate', scorefile)
        cutoffscore_df.to_csv(cuttofffile, index=False)
        return [round(optimal,6), round(cutoff,6)]
    
    def weighted_score(self):
        wscore = []
        for index, row in self.CONTIG_DF.iterrows():
            wscore.append(row['score'] * row['tpm'])
        contig_weighted = sum(wscore) / len(wscore)
        return round(contig_weighted * (self.good_mappings / self.fragments), 6)
    
    def score_run(self):
        self.score = self.scores()
        self.CONTIG_DF['score'] = self.score
        self.raw_score_val = round(self.raw_score(),6)
        self.optimal_score_val = self.optimal_score()

    def sc_metrics(self):            
        for index, row in self.CONTIG_DF.iterrows():
            try:    
                self.CONTIG_DF.loc[index, 'p_good']          = round(row['good'] / row['fragments'], 4)
            except:
                self.CONTIG_DF.loc[index, 'p_good']          = 0.0
            try:
                self.CONTIG_DF.loc[index, 'p_bases_covered'] = round((row['length'] - row['basesuncovered']) / row['length'], 4)
            except:
                self.CONTIG_DF.loc[index, 'p_bases_covered'] = 0.0
            try:
                self.CONTIG_DF.loc[index, 'sCnuc']           = row['p_seqtrue']
            except:
                self.CONTIG_DF.loc[index, 'sCnuc']           = 0.0
            try:
                self.CONTIG_DF.loc[index, 'sCcov']           = round((row['length'] - row['basesuncovered']) / row['length'], 4)
            except:
                self.CONTIG_DF.loc[index, 'sCcov']           = 0.0
            try:
                self.CONTIG_DF.loc[index, 'sCseg']           = row['p_notsegmented']
            except:
                self.CONTIG_DF.loc[index, 'sCseg']           = 0.0
            try:
                self.CONTIG_DF.loc[index, 'coverage']        = round(row['eff_count'] * row['length'] / row['eff_length'], 4)
            except:
                self.CONTIG_DF.loc[index, 'coverage']        = 0.0
            try:
                self.CONTIG_DF.loc[index, 'orf_length']      = self.good_dct_contigs[row['name']]['orf_len']
            except:
                self.CONTIG_DF.loc[index, 'orf_length']      = 0
            try:
                self.CONTIG_DF.loc[index, 'prop_gc']         = self.good_dct_contigs[row['name']]['prop_gc']
            except:
                self.CONTIG_DF.loc[index, 'prop_gc']         = 0
            if not self.SINGLE:
                try:
                    self.CONTIG_DF.loc[index, 'sCord'] = round(row['good'] / row['fragments'], 4)
                except:
                    self.CONTIG_DF.loc[index, 'sCord'] = 0.0

    def good_contigs(self):
        for index, row in self.CONTIG_DF.iterrows():
            if self.CONTIG_DF.loc[index, 'score'] > self.optimal_score_val[1]:
                self.goodcontig += 1
                self.goodcontiglist.append(self.CONTIG_DF.loc[index, 'name'])
        try:
            with open(os.path.join(self.OUTDIR, 'transrate', f'good.{self.assembly}'), 'w') as good:
                good.write('')
            with open(os.path.join(self.OUTDIR, 'transrate', f'bad.{self.assembly}'), 'w') as bad:
                bad.write('')
        except:
            pass
        with open(self.FA, 'r') as file:
            fa = file.readlines()
            lineindex = 0
            for i, line in enumerate(fa):
                frame = []
                if line.startswith('>'):
                    if line.strip('\n').strip('>') in self.goodcontiglist:
                        frame.append(line.strip('\n'))
                        while True:
                            lineindex += 1
                            if lineindex == len(fa) or fa[lineindex].startswith('>'):
                                break
                            frame.append(fa[lineindex].upper().strip('\n'))
                        with open(os.path.join(self.OUTDIR, 'transrate', f'good.{self.assembly}'), 'a') as good:
                            good.write('\n'.join(frame))
                            good.write('\n')
                    else:
                        frame.append(line.strip('\n'))
                        while True:
                            lineindex += 1
                            if lineindex == len(fa) or fa[lineindex].startswith('>'):
                                break
                            frame.append(fa[lineindex].upper().strip('\n'))
                        with open(os.path.join(self.OUTDIR, 'transrate', f'bad.{self.assembly}'), 'a') as bad:
                            bad.write('\n'.join(frame))
                            bad.write('\n')


        if 'goodcontig' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'goodcontig'] = self.goodcontig

    def assembly_output(self):
        if 'assembly' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'assembly']             = self.BASE
        if 'n_seqs' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_seqs']               = self.seqs
        if 'smallest' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'smallest']             = self.smallest
        if 'largest' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'largest']              = self.largest
        if 'n_bases' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_bases']              = self.bases
        if 'mean_len' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'mean_len']             = self.mean
        if 'median_len' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'median_len']           = self.median
        if 'std_len' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'std_len']              = self.std
        if 'n_under_200' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_under_200']          = self.n_under_200
        if 'n_over_1k' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_over_1k']            = self.n_over_1k
        if 'n_over_10k' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_over_10k']           = self.n_over_10k
        if 'n_with_orf' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n_with_orf']           = self.n_with_orf
        if 'mean_orf_percent' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'mean_orf_percent']     = self.mean_orf_percent
        if 'n90' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n90']                  = self.n90
        if 'n70' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n70']                  = self.n70
        if 'n50' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n50']                  = self.n50
        if 'n30' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n30']                  = self.n30
        if 'n10' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'n10']                  = self.n10
        if 'gc' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'gc']                   = self.gc
        if 'bases_n' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'bases_n']              = self.n_count
        if 'proportion_n' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'proportion_n']         = self.n_prop
        if 'fragments' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'fragments']            = self.fragments
        if 'fragments_mapped' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'fragments_mapped']     = self.fragments_mapped
        if 'p_fragments_mapped' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_fragments_mapped']   = self.p_fragments_mapped
        if 'good_mappings' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'good_mappings']        = self.good_mappings
        if 'p_good_mapping' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_good_mapping']       = self.p_good_mappings
        if 'bad_mappings' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'bad_mappings']         = self.bad_mappings
        if 'potential_bridges' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'potential_bridges']    = self.potential_bridges
        if 'bases_uncovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'bases_uncovered']      = self.bases_uncovered
        if 'p_bases_uncovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_bases_uncovered']    = self.p_bases_uncovered
        if 'contigs_uncovbase' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'contigs_uncovbase']    = self.contig_uncovbase
        if 'p_contigs_uncovbase' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_contigs_uncovbase']  = self.p_contig_uncovbase
        if 'contigs_uncovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'contigs_uncovered']    = self.contigs_uncovered
        if 'p_contigs_uncovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_contigs_uncovered']  = self.p_contigs_uncovered
        if 'contigs_lowcovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'contigs_lowcovered']   = self.contigs_lowcovered
        if 'p_contigs_lowcovered' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_contigs_lowcovered'] = self.p_contigs_lowcovered
        if 'contigs_segmented' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'contigs_segmented']    = self.contigs_segmented
        if 'p_contigs_segmented' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'p_contigs_segmented']  = self.p_contigs_segmented
        if 'score' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'score']                = self.raw_score_val
        if 'optimal_score' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'optimal_score']        = self.optimal_score_val[0]
        if 'cutoff' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'cutoff']               = self.optimal_score_val[1]
        if 'weighted' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'weighted']             = self.weighted_score()
        if 'sCnuc' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCnuc']                = self.CONTIG_DF['sCnuc'].mean()
        if 'sCcov' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCcov']                = self.CONTIG_DF['sCcov'].mean()
        if 'sCord' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCord']                = self.CONTIG_DF['sCord'].mean()
        if 'sCseg' in self.OUTPUT_DF.columns:
            self.OUTPUT_DF.loc[0, 'sCseg']                = self.CONTIG_DF['sCseg'].mean()

        self.good_contigs()

        outfile    = os.path.join(self.OUTDIR, 'transrate', 'assembly.csv')
        if os.path.exists(outfile) and self.MULTI:
            self.OUTPUT_DF.to_csv(outfile, index=False, mode='a', header=False)
        else:
            self.OUTPUT_DF.to_csv(outfile, index=False)

    def contig_output(self):
        cfile = 'contigs.csv' if not self.MULTI else f'{self.assembly}.contigs.csv'
        contigfile = os.path.join(self.OUTDIR, 'transrate', f'{cfile}')
        self.CONTIG_DF = self.CONTIG_DF[self.CHEADERS]
        self.CONTIG_DF.to_csv(contigfile, index=False)
