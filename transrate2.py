import __init__
import sys
import argparse
import time
import psutil
import requests
import os
from modules.run import MAIN
from modules.reqs import req_check

if __name__ == '__main__':
    class PARSER:
        def __init__(self):
            self.color         = ''
            self.assemblycount = 0
            self.clean         = {}
            self.colors        = ['\033[0;33m', '\033[0;31m', '\033[0;32m']
            self.colorcount    = 0
            self.quiet         = False
            self.reqs          = False
            self.options       = {}
            
            try:
                url = 'https://api.github.com/repos/ericbretz/transrate2/releases'
                header = {'Accept': 'application/vnd.github+json'}
                response = requests.get(url, headers=header)
                self.latest = response.json()[0]['tag_name']
            except:
                self.latest = ''

        def logoprint(self, star = False, bowtie2 = False):
            if self.quiet:
                return
            
            C = self.colors[0]
            H = '\033[m'
            self.color = C
            W = '\033[37m'
            
            rversion = f'v{__init__.__version__}'
            
            if self.latest:
                if self.latest != __init__.__version__:
                    rversion = f'{rversion} (Update Available: {self.latest})'
            else:
                pass

            B = '\033[0;32m'
            A = '\033[0;33m'
            C = '\033[0;31m'
            W = '\033[m'
            transrate_c = f'''
 {B}██████{W}┐{B}██████{W}┐  {B}█████{W}┐ {B}███{W}┐  {B}██{W}┐ {B}██████{W}┐{A}██████{W}┐  {A}█████{W}┐ {A}██████{W}┐{A}███████{W}┐{C}██████{W}┐ 
 └─{B}██{W}┌─┘{B}██{W}┌──{B}██{W}┐{B}██{W}┌──{B}██{W}┐{B}████{W}┐ {B}██{W}│{B}██{W}┌────┘{A}██{W}┌──{A}██{W}┐{A}██{W}┌──{A}██{W}┐└─{A}██{W}┌─┘{A}██{W}┌────┘└────{C}██{W}┐
   {B}▓▓{W}│  {B}▓▓▓▓▓▓{W}┌┘{B}▓▓▓▓▓▓▓{W}│{B}▓▓{W}┌{B}▓▓{W}┐{B}▓▓{W}│└{B}▓▓▓▓▓{W}┐ {A}▓▓▓▓▓▓{W}┌┘{A}▓▓▓▓▓▓▓{W}│  {A}▓▓{W}│  {A}▓▓▓▓▓{W}┐    {C}▓▓▓{W}┌─┘
   {B}▒▒{W}│  {B}▒▒{W}┌──{B}▒▒{W}┐{B}▒▒{W}┌──{B}▒▒{W}│{B}▒▒{W}│└{B}▒▒▒▒{W}│ └───{B}▒▒{W}┐{A}▒▒{W}┌──{A}▒▒{W}┐{A}▒▒{W}┌──{A}▒▒{W}│  {A}▒▒{W}│  {A}▒▒{W}┌──┘  {C}▒▒{W}┌──┘  
   {B}░░{W}│  {B}░░{W}│  {B}░░{W}│{B}░░{W}│  {B}░░{W}│{B}░░{W}│ └{B}░░░{W}│{B}░░░░░░{W}┌┘{A}░░{W}│  {A}░░{W}│{A}░░{W}│  {A}░░{W}│  {A}░░{W}│  {A}░░░░░░░{W}┐{C}░░░░░░░{W}┐
 {W}  └─┘  └─┘  └─┘└─┘  └─┘└─┘  └──┘└─────┘ └─┘  └─┘└─┘  └─┘  └─┘  └──────┘└──────┘
{W}{"Quality analysis for de-novo transcriptome assemblies":^80}
{W}             {B}░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ {A}░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ {C}░▓▓▓▓▓▓▓◠▓▓▓▓▓▓▓░ {W}   EC Bretz
{W}{rversion:>78}\033[0m
'''

            print(transrate_c)
            self.reqs = req_check(star, bowtie2)

        def aligners(self):
            print(f'{self.color}  ┌{"─" * 27}\033[m    Aligner Error   {self.color}{"─" * 27}┐\033[m')
            print(f'{self.color}  │ \033[m Please choose only one aligner to use.{self.color}{" " * 34}│\033[m')
            print(f'{self.color}  │  •\033[m STAR    (-s) is the default aligner.{self.color}{" " * 34}│\033[m')
            print(f'{self.color}  │  •\033[m Snap    (-p){" " * 58}{self.color}│\033[m')
            print(f'{self.color}  │  •\033[m Bowtie2 (-b){" " * 58}{self.color}│\033[m')
            print(f'{self.color}  │ \033[m Or use (--help) for more options.{self.color}{" " * 39}│\033[m')
            print(f'{self.color}  └{"─" * 74}┘\033[m')

        def helpoptions(self):

            topbar = f'{self.color}  ┌{"─" * 28}\033[m   Help Options   {self.color}{"─" * 28}┐\033[m'
            bottombar = f'{self.color}  └{"─" * 74}┘\033[m'

            modebar = f'{self.color}  ┌{"─" * 28}\033[m    Mode Types    {self.color}{"─" * 28}┐\033[m'
            extrabar = f'{self.color}  ┌{"─" * 28}\033[m     #Threads     {self.color}{"─" * 28}┐\033[m'
            self.options = {
                'Assembly'        : ['--assembly', '-a', 'Path to assembly file (FASTA)'],
                'Left Reads'      : ['--left', '-l', 'Path to left reads file (FASTQ)'],
                'Right Reads'     : ['--right', '-r', 'Path to right reads file (FASTQ)'],
                'Reference'       : ['--reference', '-f', 'Path to reference file (FASTA)'],
                'Output Directory': ['--outdir', '-o', 'Path to output directory'],
                'Threads'         : ['--threads', '-t', 'Number of threads to use'],
                'STAR'            : ['--star', '-s', 'Use STAR aligner (default)'],
                'Snap'            : ['--snap', '-p', 'Use Snap aligner'],
                'Bowtie2'         : ['--bowtie2', '-b', 'Use Bowtie2 aligner'],
                'Clutter'         : ['--clutter', '-c', 'Remove intermediate files'],
                'Quiet'           : ['--quiet', '-q', 'Supress terminal output'],
                'Help'            : ['--help', '-h', 'Display this help message']
            }
            modes = {
                'Assembly' : ['-a', 'Run assembly analysis only.'],
                'Reads'    : ['-a -l -r', 'Run assembly with paired-end reads analysis'],
                'All'      : ['-a -l -r -f', 'Run assembly with paired-end reads and reference'],
                'Single'   : ['-a -l', 'Run assembly with single-end reads analysis'],
                'Reference': ['-a -f', 'Run assembly with reference analysis'],
            }


            print(topbar)
            for k,v in self.options.items():
                print(f'{self.color}  │\033[m {v[0]:<20}{v[1]:<15}{v[2]:<38}{self.color}│\033[m')
            print(bottombar)
            print(modebar)
            for k,v in modes.items():
                body = f'{self.color}  │\033[m {v[0]:<20}{v[1]:<5}{self.color}'
                length = 94 - len(str(body))
                print(f'{body}{" " * (length)}│\033[m')
            print(bottombar)

            print(extrabar)
            logical = psutil.cpu_count(logical=True)
            physical = psutil.cpu_count(logical=False)
            threadslabel = f'{self.color}  │\033[m Threads: {self.color}'
            logbody = f'{self.color}  │\033[m     Logical:  '
            physbody = f'{self.color}  │\033[m     Physical: '
            threalen = 94 - len(str(threadslabel))
            loglen = 87 - len(str(logbody)) - len(str(logical))
            physlen = 87 - len(str(physbody)) - len(str(physical))
            print(f'{threadslabel}{self.color}{" " * threalen}│\033[m')
            print(f'{logbody}{logical}{" " * loglen}{self.color}│\033[m')
            print(f'{physbody}{physical}{" " * physlen}{self.color}│\033[m')
            print(bottombar)

        def parser(self):
            parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.RawTextHelpFormatter, description='Transrate2')
            parser.add_argument('--assembly', '-a', type=str , help='Assembly file')
            parser.add_argument('--left', '-l', type=str , help='Left reads file')
            parser.add_argument('--right', '-r', type=str , help='Right reads file')
            parser.add_argument('--reference', '-f', type=str , help='Reference file')
            parser.add_argument('--output', '-o', type=str , help='Output directory')
            parser.add_argument('--threads', '-t', type=int , help='Number of threads')
            parser.add_argument('--star', '-s', action='store_true', help='Use STAR aligner (default)')
            parser.add_argument('--snap', '-p', action='store_true', help='Use Snap aligner')
            parser.add_argument('--bowtie2', '-b', action='store_true', help='Use Bowtie2 aligner')
            parser.add_argument('--clutter', '-c', action='store_true', help='Remove clutter from output directory')
            parser.add_argument('--help', '-h', action='store_true', help='Display this help message')
            parser.add_argument('--quiet', '-q', action='store_true', help='Supress terminal output')
            parser.add_argument('--skip', '-k', action='store_true', help='Skip to transrate')

            args = parser.parse_args()
            

            def get_term(self):
                term = os.environ.get('TERM')
                if 'xterm' in term:
                    return True
                else:
                    return False
                
            def parameters():
                if self.quiet:
                    return
                topbar = f'{self.color}  ┌{"─" * 74}┐\033[m'
                bottombar = f'{self.color}  └{"─" * 74}┘\033[m'
                print(topbar)
                if self.assemblytotal > 1:
                    y = 38 - len(str(self.assemblycount + 1)) - len(str(self.assemblytotal))
                    print(f'{self.color}  │\033[m Assembly #{" " * 24}{self.assemblycount + 1}/{self.assemblytotal}{self.color}{" " * y}│\033[m')
                for k,v in args.__dict__.items():
                    if v:
                        if k == 'assembly':
                            v = v.strip(' ').split(',')
                            v = v[self.assemblycount]
                            self.assemblycount += 1
                        elif k == 'left' and not args.right:
                            k = 'reads'
                        elif k == 'right' and not args.left:
                            k = 'reads'
                        xlen = 74 - len(str(k)) - 40
                        x = ' ' * xlen
                        y = 64 - len(str(v)[-30:]) - 25
                        if len(str(v)) >= 30:
                            v = f'...{str(v)[-27:]}'
                        line =  f'{self.color}  │\033[m {k.capitalize()}{x}{str(v)[-30:]}{self.color}{" "*y}│\033[m'
                        print(line)
                        
                print(bottombar)
                print('')

            if args.help:
                self.logoprint(args.star, args.bowtie2)
                self.helpoptions()
                sys.exit()

            elif len(sys.argv) > 2:

                self.quiet = args.quiet if args.quiet else False

                transrate_start              = MAIN()
                transrate_start.TERM         = get_term(self)
                transrate_start.ASSEMBLYLIST = str(args.assembly).strip(' ').split(',') if args.assembly else ''
                transrate_start.ASSEMBLYLIST = [x for x in transrate_start.ASSEMBLYLIST if x.strip()]
                transrate_start.SKIP         = args.skip if args.skip else False
                self.assemblytotal           = len(transrate_start.ASSEMBLYLIST)

                if len(transrate_start.ASSEMBLYLIST) > 1:
                    transrate_start.MULTASSEMBLY = True
                self.clean                   = transrate_start.__dict__

                if args.star + args.bowtie2 + args.snap > 1:
                    self.logoprint(args.star, args.bowtie2)
                    self.aligners()
                    sys.exit()

                for x in transrate_start.ASSEMBLYLIST:
                    self.color = self.colors[self.colorcount]
                    transrate_start.__dict__  = self.clean
                    transrate_start.FINISHED  = False
                    transrate_start.ASSEMBLY  = x
                    transrate_start.BASE      = os.path.basename(transrate_start.ASSEMBLY).split('.')[0]
                    transrate_start.LEFT      = args.left if args.left else ''
                    transrate_start.RIGHT     = args.right if args.right else ''
                    transrate_start.REFERENCE = args.reference if args.reference else ''
                    transrate_start.THREADS   = args.threads if args.threads else 1
                    transrate_start.OUTDIR    = args.output if args.output else ''
                    transrate_start.CLUTTER   = args.clutter if args.clutter else False
                    transrate_start.STAR      = False if any([args.snap, args.bowtie2]) else True
                    transrate_start.BT2       = args.bowtie2 if args.bowtie2 else False
                    transrate_start.SNAP      = args.snap if args.snap else False
                    transrate_start.QUIET     = args.quiet if args.quiet else False
                    transrate_start.LOGOCOLOR = self.color

                    self.colorcount += 1
                    if self.colorcount == 3:
                        self.colorcount = 0
                    self.logoprint(args.star, args.bowtie2)
                    if self.reqs:
                        sys.exit()
                    parameters()
                    transrate_start.run()
                    while not transrate_start.FINISHED:
                        pass
                    time.sleep(1)
            else:
                self.logoprint(args.star, args.bowtie2)
                self.helpoptions()
                sys.exit()

    
    transrateparser = PARSER()
    transrateparser.parser()

