<p align="center">
<picture><img src="https://i.imgur.com/ksFvFqp.png"
     alt="TransRate 2.0.4"/><br></picture>
<b>Version 2.0.4</b><br>
Quality analysis for de-novo transcriptome assemblies</p>

Transrate2 is a reimplementation of the original bioinformatics tool, TransRate (Smith-Unna et al. 2015).<br>
https://github.com/blahah/transrate

## Usage/Examples

<code>transrate2</code> <code>-a</code> ASSEMBLY.fa <code>-l</code> LEFT.fq <code>-r</code> RIGHT.fq <code>-f</code> REFERENCE.fa <code>-o</code> OUTDIR <code>-t</code> THREADS <br>
<code>transrate2</code> <code>--help</code>

$$\Huge\textcolor{yellow}{\textsf{Under Development}}$$

## Requirements
<b>Salmon</b>&emsp;&emsp;&emsp;&emsp;https://github.com/COMBINE-lab/salmon<br>
```
apt install salmon
```
<b>Snap-aligner</b>&emsp;&nbsp;https://github.com/amplab/snap<br>
```
apt install snap-aligner
```
<b>STAR</b> (optional)&nbsp;https://github.com/alexdobin/STAR<br>
```
apt install rna-star
```
<b>Samtools</b>&emsp;&emsp;&emsp;https://github.com/samtools/samtools<br>
```
apt install samtools
```
<b>Diamond</b>&emsp;&emsp;&emsp;https://github.com/bbuchfink/diamond
```
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
```
