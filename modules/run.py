import time
import os
import subprocess
import sys
import threading
import pandas as pd
from modules import file,frag,seq,base,sgmt,good,assembly,reference
import warnings
import shutil
warnings.filterwarnings('ignore')

class MAIN:
    def __init__(self):
        #### Parameters ####
        self.ASSEMBLYLIST = []
        self.ASSEMBLY     = ''
        self.LEFT         = ''
        self.RIGHT        = ''
        self.REFERENCE    = ''
        self.THREADS      = 1
        self.BASE         = ''
        self.TERM         = False
        self.STAR         = True
        self.BT2          = False
        self.SNAP         = False
        self.MULTASSEMBLY = False
        self.QUIET        = False
        self.SINGLE       = ''
        self.SKIP         = False
        self.outputthread = None
        self.ERROR        = False
        self.STOPPED      = False

        #### Paths ####
        self.OUTDIR     = ''
        self.STARDIR    = ''
        self.STARINDEX  = ''
        self.BT2DIR     = ''
        self.BT2INDEX   = ''
        self.SNAPDIR    = ''
        self.SNAPINDEX  = ''
        self.SALMONDIR  = ''
        self.RDIR       = ''
        self.LOGDIR     = ''

        #### Files ####
        self.BAM        = ''
        self.SORTEDBAM  = ''
        self.SALMONBAM  = ''
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
            'snap_paired'   : [f'{self.RUNTIME}_snap_stdout.log',             f'{self.RUNTIME}_snap_stderr.log'],
            'snap_index'    : [f'{self.RUNTIME}_snap_index_stdout.log',       f'{self.RUNTIME}_snap_index_stderr.log'],
            'bowtie2_index' : [f'{self.RUNTIME}_bowtie2_index_stdout.log',    f'{self.RUNTIME}_bowtie2_index_stderr.log'],
            'bowtie2'       : [f'{self.RUNTIME}_bowtie2_stdout.log',          f'{self.RUNTIME}_bowtie2_stderr.log'],
            'star_index'    : [f'{self.RUNTIME}_star_index_stdout.log',       f'{self.RUNTIME}_star_index_stderr.log'],
            'star'          : [f'{self.RUNTIME}_star_stdout.log',             f'{self.RUNTIME}_star_stderr.log'],
            'salmon_quant'  : [f'{self.RUNTIME}_salmon_stdout.log',           f'{self.RUNTIME}_salmon_stderr.log'],
            'samtools_sort' : [f'{self.RUNTIME}_samtools_stdout.log',         f'{self.RUNTIME}_samtools_stderr.log'],
            'samtools_index': [f'{self.RUNTIME}_samtools_index_stdout.log',   f'{self.RUNTIME}_samtools_index_stderr.log'],
            'transrate'     : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'assembly'      : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'assembly_solo' : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'reference'     : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'file'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'frag'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'seq'           : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'good'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'base'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
            'sgmt'          : [f'{self.RUNTIME}_transrate_stdout.log',        f'{self.RUNTIME}_transrate_stderr.log'],
        }

        self.STDOUT     = ''
        self.STDERR     = ''
        self.PGID       = 0
        self.GOODCOUNT  = 0

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
        self.BT2DIR     = os.path.join(self.OUTDIR, 'bowtie2')
        self.BT2INDEX   = os.path.join(self.BT2DIR, 'bowtie2_index')

        #### Files ####
        if self.STAR:
            self.BAM    = os.path.join(self.STARDIR, f'{self.BASE}_Aligned.out.bam')
        elif self.BT2:
            self.BAM    = os.path.join(self.BT2DIR, f'{self.BASE}.bam')
        else:
            self.BAM    = os.path.join(self.SNAPDIR, f'{self.BASE}.bam')
        self.SORTEDBAM  = os.path.join(self.SALMONDIR, 'postSample.sorted.bam')
        # self.SORTEDBAM  = os.path.join(self.SALMONDIR, 'postSample.bam')

        if self.STAR:
            self.SNAPCOUNT = os.path.join(self.STARDIR, self.BASE + '_Log.final.out')
        else:
            self.SNAPCOUNT  = os.path.join(self.RDIR, 'snapcount.txt')
        csvout = os.path.join(self.RDIR, f'{self.BASE}.transrate.csv')
        self.CSVOUT     = csvout
        self.ASSEMBLIES = os.path.join(self.OUTDIR, 'assembly.csv')
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
        if self.STAGE == 'Bowtie2':
            self.SNAPCOUNT = self.STDERR

    def bowtie2_index(self):
        self.log_time('bowtie2_index', 'start')
        self.log_set('bowtie2_index')
        self.STAGE = 'Bowtie2 Index'
        self.STARTED = True

        bt2_index_cmd  = ['bowtie2-build', '--threads', str(self.THREADS), self.ASSEMBLY, self.BT2INDEX]
        bt2_index_run  = subprocess.Popen(bt2_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr  = bt2_index_run.communicate()
        returncode      = bt2_index_run.returncode
        self.PGID       = bt2_index_run.pid

        self.log_write(stdout, stderr)
        self.log_time('bowtie2_index', 'end')
        if bt2_index_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('Bowtie2 Index Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
        self.STAGEDONE = True        

    def bowtie2(self):
        self.log_time('bowtie2', 'start')
        self.log_set('bowtie2')
        self.STAGE = 'Bowtie2'
        self.STARTED = True

        # if self.LEFT.endswith('.gz') or self.RIGHT.endswith('.gz'):
        #     if self.SINGLE:
        #         al = '--al-gz'
        #         un = '--un-gz'
        #     else:
        #         al = '--al-conc-gz'
        #         un = '--un-conc-gz'
        # else:
        #         al = '--al'
        #         un = '--un'
        #         alconc = '--al-conc'
        #         unconc = '--un-conc'

        if not self.SINGLE:
            bt2_cmd = ['bowtie2', '--threads', f'{self.THREADS}', '--very-sensitive', '--phred33', '--no-unal', '--no-mixed', '--no-discordant', '--rdg', '1000,1000', '--rfg', '1000,1000', '-x', self.BT2INDEX, '-1', self.LEFT, '-2', self.RIGHT, '-S', self.BAM]
        else:
            bt2_cmd = ['bowtie2', '--threads', f'{self.THREADS}', '--very-sensitive', '--phred33', '--no-unal', '--no-mixed', '--no-discordant', '--rdg', '1000,1000', '--rfg', '1000,1000', '-x', self.BT2INDEX, '-U', self.SINGLE, '-S', self.BAM]
        
        bt2_run = subprocess.Popen(bt2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = bt2_run.communicate()
        returncode = bt2_run.returncode
        self.PGID = bt2_run.pid

        self.log_write(stdout, stderr)
        self.log_time('bowtie2', 'end')
        if bt2_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('Bowtie2 Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
        self.STAGEDONE = True

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
        if star_index_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('STAR Index Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
        self.STAGEDONE = True

    def star(self):
        self.log_time('star', 'start')
        self.log_set('star')
        self.STAGE = 'STAR'
        self.SNAPCOUNT = os.path.join(self.STARDIR, self.BASE + '_Log.final.out')
        reads = [self.SINGLE] if self.SINGLE else [self.LEFT, self.RIGHT]
        if self.LEFT.endswith('.gz') or self.RIGHT.endswith('.gz'):
            star_cmd = ['STAR', '--runThreadN', f'{self.THREADS}', '--genomeDir', self.STARINDEX, '--readFilesIn'] 
            star_cmd.extend(reads)
            star_cmd.extend([ '--outFileNamePrefix', os.path.join(self.STARDIR, self.BASE + '_'), '--readFilesCommand', 'gunzip', '-c', '--outSAMtype BAM Unsorted'])
        else:
            star_cmd = ['STAR', '--runThreadN', f'{self.THREADS}', '--genomeDir', self.STARINDEX, '--readFilesIn']
            star_cmd.extend(reads)
            star_cmd.extend(['--outFileNamePrefix', os.path.join(self.STARDIR, self.BASE + '_'), '--outSAMtype BAM Unsorted'])
        star_run = subprocess.Popen(star_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr = star_run.communicate()
        returncode = star_run.returncode

        self.log_write(stdout, stderr)
        self.log_time('star', 'end')
        if star_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('STAR Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
        self.STAGEDONE = True


    def snap_index(self):
        self.log_time('snap_index', 'start')
        self.log_set('snap_index')
        self.STAGE = 'Snap Index'
        self.STARTED = True

        snap_index_cmd  = ['snap-aligner', 'index', self.ASSEMBLY, self.SNAPINDEX, '-s', '23', f'-t{self.THREADS}', '-bSpace', '-locationSize', '4']
        snap_index_run  = subprocess.Popen(snap_index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr  = snap_index_run.communicate()
        returncode      = snap_index_run.returncode
        self.PGID       = snap_index_run.pid

        self.log_write(stdout, stderr)
        self.log_time('snap_index', 'end')
        if snap_index_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('Snap Index Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
        self.STAGEDONE = True

    def snap(self):
        self.log_time('snap_paired', 'start')
        self.log_set('snap_paired')
        self.STAGE = 'Snap Paired'
        if not self.SINGLE:
            snap_cmd        = ['snap-aligner', 'paired', self.SNAPINDEX, self.LEFT, self.RIGHT, '-o', self.BAM, '-s', '0', '1000', '-H', '300000',  '-h', '2000', '-d', '30', '-t', f'{self.THREADS}', '-b', '-M', '-D', '5', '-om', '5', '-omax', '10']
        else:
            snap_cmd        = ['snap-aligner', 'single', self.SNAPINDEX, self.SINGLE, '-o', self.BAM, '-h', '2000', '-d', '30', '-t', f'{self.THREADS}', '-b', '-D', '5', '-om', '5', '-omax', '10']
        snap_run        = subprocess.Popen(snap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)

        stdout, stderr  = snap_run.communicate()
        returncode      = snap_run.returncode
        self.PGID       = snap_run.pid
        self.SNAPCOUNT  = os.path.join(self.RDIR, 'snapcount.txt')
        with open(self.SNAPCOUNT, 'w') as snapout_f:
            snapout_f.write(stdout.decode('utf-8'))
        self.log_write(stdout, stderr)
        self.log_time('snap_paired', 'end')
        if snap_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('Snap Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
        self.STAGEDONE = True

    def salmon(self):
        self.log_time('salmon_quant', 'start')
        self.log_set('salmon_quant')
        self.STAGE = 'Salmon Quant'
        self.SALMONBAM = os.path.join(self.SALMONDIR, 'postSample.bam')
        salmon_cmd      = ['salmon', 'quant', '--libType', 'IU', '--alignments', self.BAM, '--targets', self.ASSEMBLY, f'--threads={self.THREADS}', '--sampleOut', '--sampleUnaligned', '--output', self.SALMONDIR]
        salmon_run      = subprocess.Popen(salmon_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr  = salmon_run.communicate()
        returncode      = salmon_run.returncode
        self.PGID       = salmon_run.pid

        self.log_write(stdout, stderr)
        self.log_time('salmon_quant', 'end')
        if salmon_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('Salmon Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
        self.STAGEDONE = True

    def samtools_sort(self):
        self.log_time('samtools_sort', 'start')
        self.log_set('samtools_sort')
        self.STAGE = 'Samtools Sort'

        samtools_sort_cmd = ['samtools', 'sort', f'-@{self.THREADS}', self.SALMONBAM, '-o', self.SORTEDBAM]
        samtools_sort_run = subprocess.Popen(samtools_sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid, shell=False)
        stdout, stderr    = samtools_sort_run.communicate()
        returncode        = samtools_sort_run.returncode
        self.PGID         = samtools_sort_run.pid

        self.log_write(stdout, stderr)
        self.log_time('samtools_sort', 'end')
        if samtools_sort_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('Samtools Sort Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
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
        if samtools_index_run.returncode != 0:
            self.ERROR = True
            while not self.STOPPED:
                pass
            print('\n')
            print('Samtools Index Failed')
            print(f'Check Logs at {self.LOGDIR}')
            sys.exit()
        self.STAGEDONE = True

    def transrate(self):
        self.log_time('transrate', 'start')
        self.log_set('transrate')
        self.STAGE = 'transrate'
        self.STARTED = True

        functions = [file.file, frag.frag, seq.seq, good.good, base.base, sgmt.sgmt]
        
        sys.stderr = open(os.path.join(self.LOGDIR, 'transrate.log'), 'w')

        for func in functions:
            if self.SINGLE:
                if func.__name__.lower() in ['good', 'sgmt']:
                    break
            self.RSTAGE = func.__name__.lower()
            self.DESC[func.__name__.lower()].append(time.perf_counter()) 
            self.log_time(func.__name__, 'start')
            self.log_set(func.__name__)
            self.RDCT = func(self.SORTEDBAM, self.RDCT, self.THREADS, bool(self.SINGLE))
            self.log_time(func.__name__, 'end')
            self.RSTAGEDONE = True
 
        sys.stderr = sys.__stderr__

        self.log_time('transrate', 'end')
        self.STAGEDONE = True

    def csv(self):        
        self.STAGE = 'Assembly'
        columns = ['name', 'p_seqtrue', 'bridges', 'length', 'fragments',
                    'both_mapped', 'properpair', 'good', 'basesuncovered', 'p_notsegmented'] if not self.SINGLE else \
                    ['name', 'p_seqtrue', 'length', 'fragments', 'basesuncovered']
        csvout = pd.DataFrame(columns=columns)

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
        self.STARTED = True
        self.STAGE = 'Assembly'

        assembly_main        = assembly.Assembly()
        assembly_main.CSV    = self.CSVOUT
        assembly_main.FA     = self.ASSEMBLY
        assembly_main.SNAP   = self.SNAPCOUNT
        assembly_main.OUTDIR = self.OUTDIR
        assembly_main.STAR   = self.STAR
        assembly_main.BT2    = self.BT2
        assembly_main.MULTI  = self.MULTASSEMBLY
        assembly_main.SINGLE = self.SINGLE
        assembly_main.SOLO   = True if not self.LEFT and not self.RIGHT else False
        assembly_main.run()

        self.log_time('assembly', 'end')
        self.STAGEDONE = True

    def reference(self):
        self.log_time('reference', 'start')
        self.log_set('reference')
        self.STAGE = 'Reference'
        assemblies = f'{self.RDIR}/assembly.csv'

        reference.reference(self.ASSEMBLY, self.REFERENCE, assemblies, self.RDIR, self.THREADS)

        self.log_time('reference', 'end')
        self.REFFINISHED = True
        self.STAGEDONE = True

    def pair_check(self):
        if bool(self.LEFT) != bool(self.RIGHT):
            self.SINGLE = self.LEFT if self.LEFT else self.RIGHT

    def run(self):
        self.dir_set()
        self.BASE = os.path.basename(self.ASSEMBLY).split('.')[0]
        if not self.QUIET:
            self.outputthread = threading.Thread(target=self.output)
        self.output_make()
        
        if self.LEFT or self.RIGHT:
            self.path_check()
            self.pair_check()
            self.outputthread.start() if not self.QUIET else None
            if not self.SKIP:
                if self.STAR:
                    self.star_index()
                    self.star()
                elif self.BT2:
                    self.bowtie2_index()
                    self.bowtie2()
                else:
                    self.snap_index()
                    self.snap()
                self.salmon()
                self.samtools_sort()
                self.samtools_index()
            else:
                self.RSTAGE = 'file'
            self.transrate()
            self.csv()
            # outputthread.join()
            self.assembly()
        else:
            self.path_check()
            self.outputthread.start() if not self.QUIET else None
            self.assembly()
            self.csv_solo()
            
        if self.REFERENCE:
            self.reference()
        
        if self.CLUTTER:
            self.clutter()

        self.FINISHED = True

    def clutter(self):
        delfolders = ['snap', 'salmon'] # ! Add STAR and BOWTIE2 directories
        for root, dirs, files in os.walk(self.OUTDIR):
            for file in files:
                if not file.endswith('.csv') and 'logs' not in root:
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
        if not os.path.exists(self.SNAPDIR) and self.SNAP:
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
        if not os.path.exists(self.BT2DIR) and self.BT2:
            os.makedirs(self.BT2DIR)
            if not os.path.exists(self.BT2INDEX):
                os.makedirs(self.BT2INDEX)

    def output(self):
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
                if self.ERROR:
                    self.STOPPED = True
                    sys.exit()
                    
                print(spinner(stage), end='\r')
                if key == stage.lower().replace(' ', '_'):
                    time.sleep(0.2)
            runtime = time.perf_counter() - self.LOGDCT[key]['start']
            if 'star' in key:
                x = f' {self.LOGOCOLOR}░\033[m {key.upper().replace("_", " ")} Finished'
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
                if self.ERROR:
                    self.STOPPED = True
                    sys.exit()
                    
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
                        'weighted'            : 'Weighted',
                        'goodcontig'          : 'Good Contig',
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
            if not self.SINGLE:
                print(scorelabel)
                for k,v in assemblies.items():
                    try:
                        if 'p̂' in k or 'p̂' in score[k]:
                            print(f' {self.LOGOCOLOR}  •\033[m {score[k]:<23}: {v[0]}')
                        else:
                            print(f' {self.LOGOCOLOR}  •\033[m {score[k]:<22}: {v[0]}')
                    except:
                        pass
                try:
                    if self.GOODCOUNT:
                        print(f' {self.LOGOCOLOR}  •\033[m {"Good Contigs":<22}: {self.GOODCOUNT}')
                except:
                    pass
                    

                print(f'{self.LOGOCOLOR}  └{"─" * 74}┘\033[m')

            try:
                contigfile = os.path.join(self.RDIR,'contigs.csv') if not self.MULTASSEMBLY else os.path.join(self.RDIR, os.path.basename(self.ASSEMBLY) + '.contigs.csv')
                contig = pd.read_csv(contigfile)
                scdict = {
                    'scnuc_avg' : ['Avg. sCnuc', contig['sCnuc'].mean()],
                    'sccov_avg' : ['Avg. sCcov', contig['sCcov'].mean()],
                    'scord_avg' : ['Avg. sCord', contig['sCord'].mean()],
                    'scseg_avg' : ['Avg. sCseg', contig['sCseg'].mean()]
                } if not self.SINGLE else {
                    'scnuc_avg' : ['Avg. sCnuc', contig['sCnuc'].mean()],
                    'sccov_avg' : ['Avg. sCcov', contig['sCcov'].mean()],
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

            assemblies = pd.read_csv(self.OUTDIR + '/transrate/assembly.csv').to_dict()
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
            reflabel = f'{self.LOGOCOLOR}  ┌{"─" * 25}\033[m  Reference Statistics  {self.LOGOCOLOR}{"─" * 25}┐\033[m'
            clear = ' ' * 66
            print(clear)
            print(reflabel)
            assemblies = pd.read_csv(self.RDIR + '/assembly.csv').to_dict()

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
            
