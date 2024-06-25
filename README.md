<p align="center">
<picture><img src="https://i.imgur.com/1tJvioD.png"
     alt="TransRate 2.1.3"/><br></picture>
Quality analysis for de-novo transcriptome assemblies</p>
<p align="center"><a href="#"><img alt="GitHub Tag" src="https://img.shields.io/github/v/tag/ericbretz/transrate2?style=flat-square"></a></p>

TransRate2 is a reimplementation of the original bioinformatics tool, TransRate (Smith-Unna et al. 2015).<br>
https://github.com/blahah/transrate

> [!CAUTION]
> $\Huge\textcolor[RGB]{248, 82, 73}{\textsf{TransRate2 is currently under development}}$

## Usage/Examples

<code>transrate2</code> <code>-a</code> ASSEMBLY.fa <code>-l</code> LEFT.fq <code>-r</code> RIGHT.fq <code>-f</code> REFERENCE.fa <code>-o</code> OUTDIR <code>-t</code> THREADS <br>
<code>transrate2</code> <code>--help</code>

## Requirements
| Package | Repo <img width=400px></img>|
| :--- | :--- |
| **Salmon** | https://github.com/COMBINE-lab/salmon |
| **STAR** | https://github.com/alexdobin/STAR |
| **Bowtie2** _<sub>(optional)</sub>_ | https://github.com/BenLangmead/bowtie2 |
| **Snap** _<sub>(optional)</sub>_ | https://github.com/amplab/snap |
| **Samtools** | https://github.com/samtools/samtools |
| **Diamond** | https://github.com/bbuchfink/diamond |

## Bug Reporting
Please report bugs in the <a href="https://github.com/ericbretz/transrate2/issues">Issue Tracker</a>
