import __init__
import sys
import argparse
import threading
import time
import random
import psutil
import requests
import json
import os
from modules.run import MAIN

if __name__ == '__main__':
    class PARSER:
        def __init__(self):
            self.color = ''
            try:
                url = 'https://api.github.com/repos/ericbretz/transrate/releases'
                header = {'Accept': 'application/vnd.github+json'}
                response = requests.get(url, headers=header)
                self.latest = response.json()[0]['tag_name']
            except:
                self.latest = ''

        def logoprint(self):
            colors = {
                'red': '\033[0;31m',
                'green': '\033[0;32m',
                'yellow': '\033[0;33m',
            }
            C = colors['yellow']
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

        def helpoptions(self):

            topbar = f'{self.color}  ┌{"─" * 28}\033[m   Help Options   {self.color}{"─" * 28}┐\033[m'
            bottombar = f'{self.color}  └{"─" * 74}┘\033[m'

            modebar = f'{self.color}  ┌{"─" * 28}\033[m    Mode Types    {self.color}{"─" * 28}┐\033[m'
            extrabar = f'{self.color}  ┌{"─" * 28}\033[m     #Threads     {self.color}{"─" * 28}┐\033[m'
            options = {
                'Assembly': ['--assembly', '-a', 'Path to assembly file (FASTA)'],
                'Left Reads': ['--left', '-l', 'Path to left reads file (FASTQ)'],
                'Right Reads': ['--right', '-r', 'Path to right reads file (FASTQ)'],
                'Reference': ['--reference', '-f', 'Path to reference file (FASTA)'],
                'Output Directory': ['--outdir', '-o', 'Path to output directory'],
                'Threads': ['--threads', '-t', 'Number of threads to use'],
                'STAR': ['--STAR', '-s', 'Use STAR aligner'],
                'SNAP': ['--SNAP', '-p', 'Use SNAP aligner (default)'],
                'Clutter': ['--clutter', '-c', 'Remove intermediate files'],
                'Help': ['--help', '-h', 'Display this help message']
            }
            modes = {
                'Assembly': ['-a', 'Only run assembly analysis'],
                'Reads': ['-a -l -r', 'Run assembly with reads analysis'],
                'All': ['-a -l -r -f', 'Run assembly with reads and reference'],
                'Reference': ['-a -f', 'Run assembly with reference analysis'],
            }

            print(topbar)
            for k,v in options.items():
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
            parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.RawTextHelpFormatter)
            parser.add_argument('--assembly', '-a', type=str , help='Assembly file')
            parser.add_argument('--left', '-l', type=str , help='Left reads file')
            parser.add_argument('--right', '-r', type=str , help='Right reads file')
            parser.add_argument('--reference', '-f', type=str , help='Reference file')
            parser.add_argument('--output', '-o', type=str , help='Output directory')
            parser.add_argument('--threads', '-t', type=int , help='Number of threads')
            parser.add_argument('--STAR', '-s', action='store_true', help='Use STAR aligner')
            parser.add_argument('--SNAP', '-p', action='store_true', help='Use SNAP aligner (default)')
            parser.add_argument('--clutter', '-c', action='store_true', help='Remove clutter from output directory')
            parser.add_argument('--help', '-h', action='store_true', help='Display this help message')

            args = parser.parse_args()

            def get_term(self):
                term = os.environ.get('TERM')
                if 'xterm' in term:
                    return True
                else:
                    return False
                
            def parameters():
                topbar = f'{self.color}  ┌{"─" * 74}┐\033[m'
                bottombar = f'{self.color}  └{"─" * 74}┘\033[m'
                print(topbar)
                for k,v in args.__dict__.items():
                    if v:
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
                self.helpoptions()
                sys.exit()

            elif len(sys.argv) != 1:

                transrate_start           = MAIN()
                transrate_start.TERM      = get_term(self)
                transrate_start.BASE      = os.path.basename(args.assembly).split('.')[0]
                transrate_start.ASSEMBLY  = args.assembly if args.assembly else ''
                transrate_start.LEFT      = args.left if args.left else ''
                transrate_start.RIGHT     = args.right if args.right else ''
                transrate_start.REFERENCE = args.reference if args.reference else ''
                transrate_start.THREADS   = args.threads if args.threads else 1
                transrate_start.OUTDIR    = args.output if args.output else ''
                transrate_start.CLUTTER   = args.clutter if args.clutter else False
                transrate_start.STAR      = args.STAR if args.STAR else False
                parameters()
                transrate_start.LOGOCOLOR = self.color
                transrate_start.run()
            else:
                self.helpoptions()
                sys.exit()

    
    transrateparser = PARSER()
    transrateparser.logoprint()
    transrateparser.parser()

