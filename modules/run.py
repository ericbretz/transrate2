import time
import os
import subprocess
import sys
import threading
import pandas as pd
from modules import file,frag,seq,base,sgmt,good,assembly,assembly_solo,reference
import warnings
import shutil
warnings.filterwarnings('ignore')

class MAIN:
    def __init__(self):
        #### Parameters ####
        self.ASSEMBLYLIST = []
        self.ASSEMBLY   = ''
        self.LEFT       = ''
        self.RIGHT      = ''
        self.REFERENCE  = ''
        self.THREADS    = 1
        self.BASE       = ''
        self.TERM       = False
        self.STAR       = False

        #### Paths ####
        self.OUTDIR     = ''
        self.STARDIR    = ''
        self.STARINDEX  = ''
        self.SNAPDIR    = ''
        self.SNAPINDEX  = ''
        self.SALMONDIR  = ''
        self.RDIR       = ''
        self.LOGDIR     = ''

        #### Files ####
        self.BAM        = ''
        self.SORTEDBAM  = ''
        self.RDCT       = {}
        self.SNAPCOUNT  = ''
        self.CSVOUT     = ''
        self.ASSEMBLIES = ''
        self.GOODFA     = ''

        #### Logging ####
        self.RUNTIME    = time.strftime('%Y%m%d%H%M%S')
        self.LOGDCT     = {}
        self.LOGDIR     = ''
        self.LOGFILES   = {
            'snap_paired'   : [f'{self.RUNTIME}_snap_stdout.log',           f'{self.RUNTIME}_snap_stderr.log'          ],
            'snap_index'    : [f'{self.RUNTIME}_snap_index_stdout.log',     f'{self.RUNTIME}_snap_index_stderr.log'    ],
            'star_index'    : [f'{self.RUNTIME}_star_index_stdout.log',     f'{self.RUNTIME}_star_index_stderr.log'    ],
            'star'          : [f'{self.RUNTIME}_star_stdout.log',           f'{self.RUNTIME}_star_stderr.log'          ],
            'salmon_quant'  : [f'{self.RUNTIME}_salmon_stdout.log',         f'{self.RUNTIME}_salmon_stderr.log'        ],
            'samtools_sort' : [f'{self.RUNTIME}_samtools_stdout.log',       f'{self.RUNTIME}_samtools_stderr.log'      ],
            'samtools_index': [f'{self.RUNTIME}_samtools_index_stdout.log', f'{self.RUNTIME}_samtools_index_stderr.log'],
            'transrate'       : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'assembly'      : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'assembly_solo' : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'reference'     : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'file'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'frag'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'seq'           : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'good'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'base'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
            'sgmt'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'       ],
        }

        self.STDOUT     = ''
        self.STDERR     = ''
        self.PGID       = 0

        #### Output ####
        self.STAGE     = ''
        self.RSTAGE    = ''
        self.STAGEDONE = False
        self.RSTAGEDONE= False
        self.SPINCOUNT = 0
        self.STARTED   = False
        self.LOGOCOLOR = ''
        self.FINISHED  = False
        self.STATS     = False
        self.SOLOSTATS = False
        self.CLUTTER   = False
        self.REFFINISHED = False
        self.COLOR  = '\033[0;31m'
        self.COLORS = ['\033[0;31m', '\033[0;32m', '\033[0;33m']
        self.COLORCOUNT = 0

        self.DESC = {
                'file': ['• Preparing Database'], 
                'frag': ['• Scanning for fragmentation'], 
                'seq' : ['• Assessing sequence quality'], 
                'good': ['• Determining good/bad reads'], 
                'base': ['• Locating uncovered bases'], 
                'sgmt': ['• Calculating segmentation'],
                }
        
    def good_fa(self):
        goodlst = []
        goodseq = []

        with open(self.GOODFA, 'w') as good_f:
            good_f.write('')

        for k,v in self.RDCT.items():
            try:
                for g in v['good']['goodlst']:
                    if g not in goodlst:
                        goodlst.append(g)
            except:
                pass
        
        with open(self.ASSEMBLY, 'r') as assembly_f:
            for line in assembly_f:
                if line.startswith('>'):
                    if line.strip('>').strip('\n') in goodlst:
                        goodseq = [line]
                    else:
                        with open(self.GOODFA, 'a') as good_f:
                            for g in goodseq:
                                good_f.write(g)
                        goodseq = []
                else:
                    goodseq.append(line)
                

    def path_cleanup(self, path):
        clean_path = os.path.basename(path)
        return clean_path
    
    def path_check(self):
        if self.ASSEMBLY:
            if not os.path.exists(self.ASSEMBLY):
                print(f'Assembly file {self.ASSEMBLY} does not exist')
                sys.exit()
        if self.LEFT:
            if not os.path.exists(self.LEFT):
                print(f'Left reads file {self.LEFT} does not exist')
                sys.exit()
        if self.RIGHT:
            if not os.path.exists(self.RIGHT):
                print(f'Right reads file {self.RIGHT} does not exist')
                sys.exit()
        if self.REFERENCE:
            if not os.path.exists(self.REFERENCE):
                print(f'Reference file {self.REFERENCE} does not exist')
                sys.exit()

    def dir_set(self):
        #### Directories ####
        self.SNAPDIR    = os.path.join(self.OUTDIR, 'snap')
        self.SNAPINDEX  = os.path.join(self.SNAPDIR, 'snap_index')
        self.STARDIR    = os.path.join(self.OUTDIR, 'star')
        self.STARINDEX  = os.path.join(self.STARDIR, 'star_index')
        self.SALMONDIR  = os.path.join(self.OUTDIR, 'salmon')
        self.RDIR       = os.path.join(self.OUTDIR, 'transrate')
        self.LOGDIR     = os.path.join(self.OUTDIR, 'logs')

        #### Files ####
        if self.STAR:
            self.BAM    = os.path.join(self.STARDIR, f'{self.BASE}_Aligned.out.bam')
        else:
            self.BAM    = os.path.join(self.SNAPDIR, f'{self.BASE}.bam')
        self.SORTEDBAM  = os.path.join(self.SALMONDIR, 'postSample.sorted.bam')
        if self.STAR:
            self.SNAPCOUNT = os.path.join(self.STARDIR, self.BASE + '_Log.final.out')
        else:
            self.SNAPCOUNT  = os.path.join(self.RDIR, 'snapcount.txt')
        self.CSVOUT     = os.path.join(self.RDIR, 'transrate.csv')
        self.ASSEMBLIES = os.path.join(self.OUTDIR, 'assemblies.csv')
        self.GOODFA     = os.path.join(self.RDIR, f'good.{self.BASE}.fa')
    
    def log_time(self, process, type):
        if type == 'start':
            self.LOGDCT[process] = {}
            self.LOGDCT[process]['start'] = time.perf_counter()
        elif type == 'end':
            self.LOGDCT[process]['end'] = time.perf_counter()
        else:
            pass

    def log_set(self, process):
        self.STDOUT = os.path.join(self.LOGDIR, self.LOGFILES[process][0])
        self.STDERR = os.path.join(self.LOGDIR, self.LOGFILES[process][1])

    def log_write(self, stdout, stderr):
        with open(self.STDOUT, 'a') as stdout_f:
            stdout_f.write(stdout.decode('utf-8'))
        with open(self.STDERR, 'a') as stderr_f:
            stderr_f.write(stderr.decode('utf-8'))

    def star_index(self):
        self.log_time('star_index', 'start')
        self.log_set('star_index')
        self.STAGE = 'STAR Index'
        self.STARTED = True

        star_index_cmd  = ['STAR', '--runThreadN', f'{self.THREADS}', '--runMode', 'genomeGenerate', '--genomeDir', self.STARINDEX, '--genomeFastaFiles', self.ASSEMBLY, '--genomeSAindexNbases', '11', '--outFileNamePrefix', os.path.join(self.STARINDEX, 'star_')]
        star_index_run  = subprocess.Popen(star_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr  = star_index_run.communicate()
        returncode      = star_index_run.returncode
        self.PGID       = star_index_run.pid

        self.log_write(stdout, stderr)
        self.log_time('star_index', 'end')
        self.STAGEDONE = True

    def star(self):
        self.log_time('star', 'start')
        self.log_set('star')
        self.STAGE = 'STAR'
        self.SNAPCOUNT = os.path.join(self.STARDIR, self.BASE + '_Log.final.out')

        if self.LEFT.endswith('.gz') or self.RIGHT.endswith('.gz'):
            star_cmd = ['STAR', '--runThreadN', f'{self.THREADS}', '--genomeDir', self.STARINDEX, '--readFilesIn', self.LEFT, self.RIGHT, '--outFileNamePrefix', os.path.join(self.STARDIR, self.BASE + '_'), '--readFilesCommand', 'gunzip', '-c', '--outSAMtype BAM Unsorted']
        else:
            star_cmd = ['STAR', '--runThreadN', f'{self.THREADS}', '--genomeDir', self.STARINDEX, '--readFilesIn', self.LEFT, self.RIGHT, '--outFileNamePrefix', os.path.join(self.STARDIR, self.BASE + '_'), '--outSAMtype BAM Unsorted']
        star_run = subprocess.Popen(star_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = star_run.communicate()
        returncode = star_run.returncode

        self.log_write(stdout, stderr)
        self.log_time('star', 'end')
        self.STAGEDONE = True


    def snap_index(self):
        self.log_time('snap_index', 'start')
        self.log_set('snap_index')
        self.STAGE = 'Snap Index'
        self.STARTED = True

        snap_index_cmd  = ['snap-aligner', 'index', self.ASSEMBLY, self.SNAPINDEX, f'-t{self.THREADS}', '-bSpace', '-locationSize', '4']
        snap_index_run  = subprocess.Popen(snap_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr  = snap_index_run.communicate()
        returncode      = snap_index_run.returncode
        self.PGID       = snap_index_run.pid

        self.log_write(stdout, stderr)
        self.log_time('snap_index', 'end')
        self.STAGEDONE = True

    def snap(self):
        self.log_time('snap_paired', 'start')
        self.log_set('snap_paired')
        self.STAGE = 'Snap Paired'

        snap_cmd        = ['snap-aligner', 'paired', self.SNAPINDEX, self.LEFT, self.RIGHT, '-o', self.BAM, '-s', '0', '1000', '-H', '300000', '-h', '2000', '-d', '30', '-t', f'{self.THREADS}', '-b', '-M', '-D', '5', '-om', '5', '-omax', '10', '-mcp', '10000000']
        snap_run        = subprocess.Popen(snap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr  = snap_run.communicate()
        returncode      = snap_run.returncode
        self.PGID       = snap_run.pid
        self.SNAPCOUNT  = os.path.join(self.RDIR, 'snapcount.txt')
        with open(self.SNAPCOUNT, 'w') as snapout_f:
            snapout_f.write(stdout.decode('utf-8'))
        self.log_write(stdout, stderr)
        self.log_time('snap_paired', 'end')
        self.STAGEDONE = True

    def salmon(self):
        self.log_time('salmon_quant', 'start')
        self.log_set('salmon_quant')
        self.STAGE = 'Salmon Quant'

        salmon_cmd      = ['salmon', 'quant', '--libType', 'IU', '--alignments', self.BAM, '--targets', self.ASSEMBLY, f'--threads={self.THREADS}', '--sampleOut', '--sampleUnaligned', '--output', self.SALMONDIR, '--noEffectiveLengthCorrection']
        salmon_run      = subprocess.Popen(salmon_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr  = salmon_run.communicate()
        returncode      = salmon_run.returncode
        self.PGID       = salmon_run.pid

        self.log_write(stdout, stderr)
        self.log_time('salmon_quant', 'end')
        self.STAGEDONE = True

    def samtools_sort(self):
        self.log_time('samtools_sort', 'start')
        self.log_set('samtools_sort')
        self.STAGE = 'Samtools Sort'

        samtools_sort_cmd = ['samtools', 'sort', f'-@{self.THREADS}', self.BAM, '-o', self.SORTEDBAM]
        samtools_sort_run = subprocess.Popen(samtools_sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr    = samtools_sort_run.communicate()
        returncode        = samtools_sort_run.returncode
        self.PGID         = samtools_sort_run.pid

        self.log_write(stdout, stderr)
        self.log_time('samtools_sort', 'end')
        self.STAGEDONE = True

    def samtools_index(self):
        self.log_time('samtools_index', 'start')
        self.log_set('samtools_index')
        self.STAGE = 'Samtools Index'

        samtools_index_cmd = ['samtools', 'index', '-b', f'-@{self.THREADS}', self.SORTEDBAM]
        samtools_index_run = subprocess.Popen(samtools_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr     = samtools_index_run.communicate()

        self.log_write(stdout, stderr)
        self.log_time('samtools_index', 'end')
        self.STAGEDONE = True

    def transrate(self):
        self.log_time('transrate', 'start')
        self.log_set('transrate')
        self.STAGE = 'transrate'
        self.STARTED = True

        functions = [file.file, frag.frag, seq.seq, good.good, base.base, sgmt.sgmt]
        
        sys.stderr = open(os.path.join(self.LOGDIR, 'transrate.log'), 'w')

        for func in functions:
            self.RSTAGE = func.__name__.lower()
            self.DESC[func.__name__.lower()].append(time.perf_counter()) 
            self.log_time(func.__name__, 'start')
            self.log_set(func.__name__)
            
            self.RDCT = func(self.SORTEDBAM, self.RDCT, self.THREADS)
            self.log_time(func.__name__, 'end')
            self.RSTAGEDONE = True
 
        sys.stderr = sys.__stderr__

        self.log_time('transrate', 'end')
        self.STAGEDONE = True

    def csv(self):        
        self.STAGE = 'Assembly'
        csvout = pd.DataFrame(columns=['name', 'p_seqtrue', 'bridges', 'length', 'fragments',
                                           'both_mapped', 'properpair', 'good', 'basesuncovered',
                                           'p_notsegmented'])

        for key, value in self.RDCT.items():
            csvout.loc[len(csvout)] = value['stats']

        csvout.to_csv(self.CSVOUT, index=False)
        self.STATS = True

    def csv_solo(self):
        self.STAGE = 'Assembly'
        self.SOLOSTATS = True


    def assembly(self):
        self.log_time('assembly', 'start')
        self.log_set('assembly')
        self.STAGE = 'Assembly'

        assembly.assembly(self.CSVOUT, self.ASSEMBLY, self.SNAPCOUNT, self.OUTDIR, self.STAR)

        self.log_time('assembly', 'end')
        self.STAGEDONE = True

    def assembly_solo(self):
        self.log_time('assembly', 'start')
        self.log_set('assembly')
        self.STAGE = 'Assembly'
        self.STARTED = True

        assembly_solo.assembly_solo(self.ASSEMBLY, self.OUTDIR)
        
        self.log_time('assembly', 'end')
        self.STAGEDONE = True

    def reference(self):
        self.log_time('reference', 'start')
        self.log_set('reference')
        self.STAGE = 'Reference'
        if self.LEFT and self.RIGHT:
            assemblies = f'{self.RDIR}/assembly.csv'
        else:
            assemblies = f'{self.OUTDIR}/assemblies.csv'

        reference.reference(self.ASSEMBLY, self.REFERENCE, assemblies, self.RDIR, self.THREADS)

        self.log_time('reference', 'end')
        self.REFFINISHED = True
        self.STAGEDONE = True


    def run(self):
        # self.ASSEMBLY = assembly
        self.dir_set()
        self.BASE = os.path.basename(self.ASSEMBLY).split('.')[0]
        outputthread = threading.Thread(target=self.output)
        self.output_make()
        if self.LEFT and self.RIGHT:
            self.path_check()
            outputthread.start()
            if self.STAR:
                self.star_index()
                self.star()
            else:
                self.snap_index()
                self.snap()
            self.salmon()
            self.samtools_sort()
            self.samtools_index()
            self.transrate()
            self.csv()
            self.assembly()
            self.good_fa()
        else:
            self.path_check()
            outputthread.start()
            self.assembly_solo()
            self.csv_solo()
            
        if self.REFERENCE:
            self.reference()
        
        if self.CLUTTER:
            self.clutter()

        self.FINISHED = True

    def clutter(self):
        savefiles = ['transrate.csv', 'assemblies.csv', 'assembly.csv', 'contigs.csv']
        delfolders = ['snap', 'salmon']
        for root, dirs, files in os.walk(self.OUTDIR):
            for file in files:
                if file not in savefiles and 'logs' not in root:
                    os.remove(os.path.join(root, file))
            for dir in dirs:
                if dir in delfolders or os.path.getsize(os.path.join(root, dir)) == 64:
                    shutil.rmtree(os.path.join(root, dir))

                
    def output_make(self):
        if not os.path.exists(self.OUTDIR) and self.OUTDIR != '':
            os.makedirs(self.OUTDIR)
        elif self.OUTDIR == '':
            self.OUTDIR = os.path.join(os.getcwd())
        if not os.path.exists(self.LOGDIR):
            os.makedirs(self.LOGDIR)
        if not os.path.exists(self.SNAPDIR) and not self.STAR:
            os.makedirs(self.SNAPDIR)
            if not os.path.exists(self.SNAPINDEX):
                os.makedirs(self.SNAPINDEX)
        if not os.path.exists(self.STARDIR) and self.STAR:
            os.makedirs(self.STARDIR)
            if not os.path.exists(self.STARINDEX):
                os.makedirs(self.STARINDEX)
        if not os.path.exists(self.SALMONDIR):
            os.makedirs(self.SALMONDIR)
        if not os.path.exists(self.RDIR):
            os.makedirs(self.RDIR)


    def output(self):
        #[X] Current stage "running" and "finished"
        #[X] Assembly statistics
        while not self.STARTED:
            pass
        stage = self.STAGE

        def color_change():
            if self.STAGEDONE:
                return
            if self.COLORCOUNT == 3:
                self.COLORCOUNT = 0
            self.COLOR = self.COLORS[self.COLORCOUNT]
            self.COLORCOUNT += 1

        def spinner(stage):
            if self.STAGEDONE:
                return
            if self.SPINCOUNT == 0:
                color_change()

            if self.TERM:
                # spinner = ['◜ ', ' ◝', ' ◞', '◟ ']
                spinner = ['', '▏', '▎', '▍', '▌', '▋', '▊', '▉', '▊', '▋', '▌', '▍', '▎', '▏',]
            else:
                spinner = ['🮠', '🮧', '🮡', '🮥', '🮣', '🮦', '🮢', '🮤', '🮭', '🮮', '🮫', '🮤', '🮠', '🮡', '🮣', '🮢',]
            

            self.SPINCOUNT = self.SPINCOUNT + 1 if self.SPINCOUNT < (len(spinner) - 1) else 0
            try:
                description = self.DESC[stage.lower().replace(' ', '_')][0]
                out = f' {self.COLOR}█\033[m {description} {self.COLOR}  {spinner[self.SPINCOUNT]}\033[m'
            except Exception as e:
                description = stage + ' Running'
                out = f' {self.COLOR}█\033[m {description} {self.COLOR}  {spinner[self.SPINCOUNT]}\033[m'
            return out
        
        def running(stage):
            key = stage.lower().replace(' ', '_')
            self.LOGDCT[key] = {}
            self.LOGDCT[key]['start'] = time.perf_counter()
            while not self.STAGEDONE:
                print(spinner(stage), end='\r')
                if key == stage.lower().replace(' ', '_'):
                    time.sleep(0.2)
            runtime = time.perf_counter() - self.LOGDCT[key]['start']
            x = f' {self.LOGOCOLOR}░\033[m {key.capitalize().replace("_", " ")} Finished'
            xlen = 44 - len((str(x)))
            out = f'{x}{" " * xlen}({round(runtime, 2)}s)'
            clear = ' ' * 66
            print(clear, end='\r')
            print(out)


            self.STAGEDONE = False

        def rrunning():
            prevrstage = ''
            printed = False
            prevrdesc = ''
            self.LOGDCT['transrate'] = {}
            self.LOGDCT['transrate']['start'] = time.perf_counter()
            print(f' {self.LOGOCOLOR}░\033[m Transrate Running')
            while self.STAGE == 'transrate' or not printed:
                if self.RSTAGE != prevrstage:
                    out = f' {self.LOGOCOLOR}░\033[m {prevrdesc}'
                    xlen = 60 - len((str(out)))
                    if prevrdesc != '':
                        print(out + (' ' * xlen))
                    prevrstage = self.RSTAGE
                elif self.STAGE != 'transrate':
                    out = f' {self.LOGOCOLOR}░\033[m {prevrdesc}'
                    xlen = 60 - len((str(out)))
                    if prevrdesc != '':
                        print(out + (' ' * xlen))
                    printed = True
                else:
                    prevrdesc = self.DESC[self.RSTAGE][0]
                    print(spinner(self.RSTAGE), end='\r')
                    if prevrstage == self.RSTAGE:
                        time.sleep(0.2)

            runtime = time.perf_counter() - self.LOGDCT['transrate']['start']
            x = f' {self.LOGOCOLOR}░\033[m Transrate Finished'
            xlen = 44 - len((str(x)))
            out = f'{x}{" " * xlen}({round(runtime, 2)}s)'
            print(out)

            self.STAGEDONE = False


        def stats():
            clear = ' ' * 66
            print(clear)
            contiglabel = f'{self.LOGOCOLOR}  ┌{"─" * 28}\033[m  Contig Metrics  {self.LOGOCOLOR}{"─" * 28}┐\033[m'
            readmaplabel = f'{self.LOGOCOLOR}  ┌{"─" * 28}\033[m   Read Mapping   {self.LOGOCOLOR}{"─" * 28}┐\033[m'
            scorelabel = f'{self.LOGOCOLOR}  ┌{"─" * 28}\033[m  Quality Scores  {self.LOGOCOLOR}{"─" * 28}┐\033[m'
            
            contig = {
                        'n_seqs'              : '# Seqs',
                        'smallest'            : 'Smallest',
                        'largest'             : 'Largest',
                        'n_bases'             : '# Bases',
                        'mean_len'            : 'Mean Len',
                        'median_len'          : 'Median Len',
                        'std_len'             : 'Std Len',
                        'n_under_200'         : '# < 200',
                        'n_over_1k'           : '# > 1k',
                        'n_over_10k'          : '# > 10k',
                        'n_with_orf'          : '# ORF',
                        'mean_orf_percent'    : 'Mean ORF %',
                        'n90'                 : 'N90',
                        'n70'                 : 'N70',
                        'n50'                 : 'N50',
                        'n30'                 : 'N30',
                        'n10'                 : 'N10',
                        'gc'                  : 'GC',
                        'bases_n'             : '# N',
                        'proportion_n'        : 'p̂ N',
                        }
            readmap = {
                        'fragments'           : '# Fragments',
                        'fragments_mapped'    : '# Fragments Mapped',
                        'p_fragments_mapped'  : 'p̂ Fragments Mapped',
                        'good_mappings'       : '# Good Mappings',
                        'p_good_mapping'      : 'p̂ Good Mappings',
                        'bad_mappings'        : '# Bad Mappings',
                        'potential_bridges'   : '# Potential Bridges',
                        'bases_uncovered'     : '# Bases Uncovered',
                        'p_bases_uncovered'   : 'p̂ Bases Uncovered',
                        'contigs_uncovbase'   : '# Contigs Uncovbase',
                        'p_contigs_uncovbase' : 'p̂ Contigs Uncovbase',
                        'contigs_uncovered'   : '# Contigs Uncovered',
                        'p_contigs_uncovered' : 'p̂ Contigs Uncovered',
                        'contigs_lowcovered'  : '# Contigs Lowcovered',
                        'p_contigs_lowcovered': 'p̂ Contigs Lowcovered',
                        'contigs_segmented'   : '# Contigs Segmented',
                        'p_contigs_segmented' : 'p̂ Contigs Segmented',
                        }
            score = {
                        'score'               : 'Score',
                        'optimal_score'       : 'Optimal Score',
                        'cutoff'              : 'Cutoff',
                        'weighted'            : 'Weighted'
            }

            print(contiglabel)
            assemblies = pd.read_csv(self.RDIR + '/assembly.csv').to_dict()

            for k,v in assemblies.items():
                try:
                    if 'p̂' in k or 'p̂' in contig[k]:
                        print(f' {self.LOGOCOLOR}  •\033[m {contig[k]:<23}: {v[0]}')
                    else:
                        print(f' {self.LOGOCOLOR}  •\033[m {contig[k]:<22}: {v[0]}')
                except:
                    pass
            print(f'{self.LOGOCOLOR}  └{"─" * 74}┘\033[m')

            print(readmaplabel)
            for k,v in assemblies.items():
                try:
                    if 'p̂' in k or 'p̂' in readmap[k]:
                        print(f' {self.LOGOCOLOR}  •\033[m {readmap[k]:<23}: {v[0]}')
                    else:
                        print(f' {self.LOGOCOLOR}  •\033[m {readmap[k]:<22}: {v[0]}')
                except:
                    pass
            print(f'{self.LOGOCOLOR}  └{"─" * 74}┘\033[m')
            print(scorelabel)
            for k,v in assemblies.items():
                try:
                    if 'p̂' in k or 'p̂' in score[k]:
                        print(f' {self.LOGOCOLOR}  •\033[m {score[k]:<23}: {v[0]}')
                    else:
                        print(f' {self.LOGOCOLOR}  •\033[m {score[k]:<22}: {v[0]}')
                except:
                    pass

            print(f'{self.LOGOCOLOR}  └{"─" * 74}┘\033[m')

            try:
                contig = pd.read_csv(self.RDIR + '/contigs.csv')
                scdict = {
                    'scnuc_avg' : ['Avg. sCnuc', contig['sCnuc'].mean()],
                    'sccov_avg' : ['Avg. sCcov', contig['sCcov'].mean()],
                    'scord_avg' : ['Avg. sCord', contig['sCord'].mean()],
                    'scseg_avg' : ['Avg. sCseg', contig['sCseg'].mean()]
                }
                sclabel = f'{self.LOGOCOLOR}  ┌{"─" * 28}\033[m      Scores      {self.LOGOCOLOR}{"─" * 28}┐\033[m'
                print(sclabel)
                for k,v in scdict.items():
                    print(f' {self.LOGOCOLOR}  •\033[m {v[0]:<22}: {round(v[1],3)}')
                print(f'{self.LOGOCOLOR}  └{"─" * 74}┘\033[m')
            except Exception as e:
                print(e)

            print(clear)
            self.STATS = False

        def stats_solo():            
            clear = ' ' * 80
            print(clear)
            contiglabel = f'{self.LOGOCOLOR}  ┌{"─" * 28}\033[m  Contig Metrics  {self.LOGOCOLOR}{"─" * 28}┐\033[m'

            assemblies = pd.read_csv(self.OUTDIR + '/assemblies.csv').to_dict()
            contig = {
                        'n_seqs'              : '# Seqs',
                        'smallest'            : 'Smallest',
                        'largest'             : 'Largest',
                        'n_bases'             : '# Bases',
                        'mean_len'            : 'Mean Len',
                        'median_len'          : 'Median Len',
                        'std_len'             : 'Std Len',
                        'n_under_200'         : '# < 200',
                        'n_over_1k'           : '# > 1k',
                        'n_over_10k'          : '# > 10k',
                        'n_with_orf'          : '# ORF',
                        'mean_orf_percent'    : 'Mean ORF %',
                        'n90'                 : 'N90',
                        'n70'                 : 'N70',
                        'n50'                 : 'N50',
                        'n30'                 : 'N30',
                        'n10'                 : 'N10',
                        'gc'                  : 'GC',
                        'bases_n'             : '# N',
                        'proportion_n'        : 'p̂ N',
                        }
            
            print(contiglabel)

            for k,v in assemblies.items():
                try:
                    if 'p̂' in k or 'p̂' in contig[k]:
                        print(f' {self.LOGOCOLOR}  •\033[m {contig[k]:<23}: {v[0]}')
                    else:
                        print(f' {self.LOGOCOLOR}  •\033[m {contig[k]:<22}: {v[0]}')
                except:
                    pass
            print(f'{self.LOGOCOLOR}  └{"─" * 74}┘\033[m')
            print(clear)
            self.SOLOSTATS = False

        def ref():
            reflabel = f'{self.LOGOCOLOR}  ┌{"─" * 20}\033[m  Reference Statistics  {self.LOGOCOLOR}{"─" * 20}┐\033[m'
            clear = ' ' * 66
            print(clear)
            print(reflabel)
            if self.LEFT and self.RIGHT:
                assemblies = pd.read_csv(self.RDIR + '/assembly.csv').to_dict()
            else:
                assemblies = pd.read_csv(self.OUTDIR + '/assemblies.csv').to_dict()

            refs = {
                'CRBB_hits' : 'CRBB Hits',
                'n_contigs_with_CRBB' : '# Contigs with CRBB',
                'p_contigs_with_CRBB' : 'p̂ Contigs with CRBB',
                'rbh_per_reference' : 'RBH per Reference',
                'n_refs_with_CRBB' : '# Refs with CRBB',
                'p_refs_with_CRBB' : 'p̂ Refs with CRBB',
                'cov25' : 'Coverage 25',
                'p_cov25' : 'p̂ Coverage 25',
                'cov50' : 'Coverage 50',
                'p_cov50' : 'p̂ Coverage 50',
                'cov75' : 'Coverage 75',
                'p_cov75' : 'p̂ Coverage 75',
                'cov85' : 'Coverage 85',
                'p_cov85' : 'p̂ Coverage 85',
                'cov95' : 'Coverage 95',
                'p_cov95' : 'p̂ Coverage 95',
                'reference_coverage' : 'Reference Coverage',
            }

            for k,v in assemblies.items():
                try:
                    if 'p̂' in k or 'p̂' in refs[k]:
                        print(f' {self.LOGOCOLOR}  •\033[m {refs[k]:<23}: {v[0]}')
                    else:
                        print(f' {self.LOGOCOLOR}  •\033[m {refs[k]:<22}: {v[0]}')
                except:
                    pass
            print(f'{self.LOGOCOLOR}  └{"─" * 74}┘\033[m')

        while not self.FINISHED:
            if self.STAGE == 'transrate':
                rrunning()
            else:
                running(stage)
            stage = self.STAGE
            if self.STATS:
                stats()
                if not self.REFERENCE:
                    return
                else:
                    self.STATS = False
            elif self.SOLOSTATS:
                stats_solo()
                self.SOLOSTATS = False
                if not self.REFERENCE:
                    return
                else:
                    self.SOLOSTATS = False
            elif self.REFFINISHED:
                ref()
                return
            
