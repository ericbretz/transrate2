<p align="center">
<picture><img src="https://i.imgur.com/ksFvFqp.png"
     alt="TransRate 2.0.5"/><br></picture>
Quality analysis for de-novo transcriptome assemblies</p>
<p align="center"><a href="#"><img alt="GitHub Tag" src="https://img.shields.io/github/v/tag/ericbretz/transrate2"></a> <a href="#"><img alt='GitHub Clones' src='https://img.shields.io/badge/dynamic/json?color=success&label=Clone&query=count&url=https://gist.githubusercontent.com/ericbretz/b15efc70dafd2fa92db05edc588c1fae/raw/clone.json&logo=github'></a> <a href="#"><img alt="GitHub Downloads (all assets, all releases)" src="https://img.shields.io/github/downloads/ericbretz/transrate2/total"></a></p>

TransRate2 is a reimplementation of the original bioinformatics tool, TransRate (Smith-Unna et al. 2015).<br>
https://github.com/blahah/transrate

## Usage/Examples

<code>transrate2</code> <code>-a</code> ASSEMBLY.fa <code>-l</code> LEFT.fq <code>-r</code> RIGHT.fq <code>-f</code> REFERENCE.fa <code>-o</code> OUTDIR <code>-t</code> THREADS <br>
<code>transrate2</code> <code>--help</code>

$$\Huge\textcolor{yellow}{\textsf{Under Development}}$$

## Requirements
<b>Salmon</b>&emsp;&emsp;&emsp;&emsp;https://github.com/COMBINE-lab/salmon<br>

<b>Snap-aligner</b>&emsp;&nbsp;https://github.com/amplab/snap<br>

<b>STAR</b> (optional)&nbsp;https://github.com/alexdobin/STAR<br>

<b>Samtools</b>&emsp;&emsp;&emsp;https://github.com/samtools/samtools<br>

<b>Diamond</b>&emsp;&emsp;&emsp;https://github.com/bbuchfink/diamond

## Bug Reporting
Please report bugs in the <a href="https://github.com/ericbretz/transrate2/issues">Issue Tracker</a>
