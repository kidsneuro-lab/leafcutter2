# LeafCutter2

## Introduction

This script takes junction files (e.g., produced by regtools) as input and constructs intron clusters from them. Then, it processes the intron clusters to identify rarely spliced introns based on certain filtering cut-offs. The output is a text file that follows the same format as the standard leafcutter tool. The first column of each row shows the genome coordinates of introns and labels them as **UP** (unproductive), **PR** (productive/protein-coding), **NE** (ambiguous in their functional effect) or **IN** (intergenic).


Main input files:
    - junction files, e.g., processed using `regtools extract junctions`. 


Main output files:

- `{out_prefix}_perind.counts.noise.gz`: output functional introns (intact), and  noisy introns. Note the start and end coordinates of noisy introns are merged to the min(starts) and max(ends) of all functional introns within cluster.

- `{out_prefix}_perind_numers.counts.noise.gz`: same as above, except write numerators.

- `{out_prefix}_perind.counts.noise_by_intron.gz`: same as the first output, except here noisy introns' coordinates are kept as their original coordinates. This is useful for diagnosis purposes. 


## Clustering, classifying and quantifying splice junctions

A basic LeafCutter2 run work works as follows:

```
python scripts/leafcutter2.py \
    -j junction_files.txt \
    -r output_dir \
    -o leafcutter2 \
    -A gtf_file.gtf \
    -G genome.fa
```

This mode first generate intron clusters based on the junction files. Then it counts junction reads towards each classified
introns.
-    `junction_files.txt` should be a text file listing path to each junction file, one path per line.
-    `-r` specify the directory of output. 
-    `-o` specify the prefix of output file names (not including directory name).
-    `-A` a GTF annotation with genes, start codons and stop codons.
-    `-G` a genome assembly FASTA file. It must correspond to the same assembly as the GTF file.


**Note:** LeafCutter2 will still run if you skip the `-A` and `-G` parameters, but it will not classify the splice junctions.


## Pre-clustering splice junctions

Generating intron clusters first can save time form multiple subsequent runs. The script `scripts/leafcutter_make_clusters.py`
makes intron clusters separately that can be later used as input for `scripts/leafcutter2.py`. You can generate clusters by running: 

```
python scripts/leafcutter_make_clusters.py \
    -j junction_files.txt \
    -r output_dir \
    -o leafcutter2 
```
This will generate a file named `output_dir/leafcutter2_refined_noisy` that can be later used as an input for LeafCutter2. This will skip the clustering step:

```
python scripts/leafcutter2.py \
    -j junction_files.txt \
    -r output_dir \
    -o leafcutter2 \
    -A gtf_file.gtf \
    -G genome.fa \
    -c output_dir/leafcutter2_refined_noisy
```


### Parameters

```
python scripts/leafcutter2.py -h

usage: leafcutter2.py [-h] -j JUNCFILES 
                              [-o OUTPREFIX] [-q] [-r RUNDIR]
                              [-l MAXINTRONLEN] [-m MINCLUREADS]
                              [-M MINREADS] [-p MINCLURATIO]
                              [-c CLUSTER] [-k] [-C]
                              [-N NOISECLASS] [-f OFFSET] [-T]

required arguments:
  -j JUNCFILES, --juncfiles JUNCFILES
                        text file with all junction files to be processed
                        
optional arguments:
  -h, --help            show this help message and exit
  -r RUNDIR, --rundir RUNDIR
                        write to directory (default ./)
  -o OUTPREFIX, --outprefix OUTPREFIX
                        output prefix (default leafcutter)
  -q, --quiet           don't print status messages to stdout
  -A ANNOTATION, --annot ANNOTATION
                        GTF annotation file
  -G GENOME, --genome GENOME
                        Genome fasta file
  -N MAXJUNCS, --max_juncs MAXJUNCS
                        skip solveNMD function if gene contains more than N juncs. 
                        Juncs in skipped genes are assigned Coding=False. Default 10000
  -l MAXINTRONLEN, --maxintronlen MAXINTRONLEN
                        maximum intron length in bp (default 100,000bp)
  -m MINCLUREADS, --minclureads MINCLUREADS
                        minimum reads in a cluster (default 30 reads)
  -M MINREADS, --minreads MINREADS
                        minimum reads for a junction to be considered for
                        clustering(default 5 reads)
  -p MINCLURATIO, --mincluratio MINCLURATIO
                        minimum fraction of reads in a cluster that support a
                        junction (default 0.001)
  -D MINREADSTD, --minreadstd MINREADSTD
                        minimum standard deviation of reads across samples for a \
                        junction to be included in output (default 0.5)
  -p MINCLURATIO, --mincluratio MINCLURATIO
                        minimum fraction of reads in a cluster that support a junction \
                        (default 0.001)
  -c CLUSTER, --cluster CLUSTER
                        refined cluster file when clusters are already made
  -k, --checkchrom      check that the chromosomes are well formated e.g. chr1,
                        chr2, ..., or 1, 2, ...
  -C, --includeconst    also include constitutive introns
  -f OFFSET, --offset OFFSET
                        Offset sometimes useful for off by 1 annotations. (default 0)
  -T, --keeptemp        keep temporary files. (default false)
```


## Prerequisites

- Minimum python version - `python v3.6`
- It is recommended to use `regtools` to junction files from BAM files. E.g.: `regtools junctions extract -a 8 -i 50 -I 500000 bamfile.bam -o outfile.junc` . See detailed regtools documentations [here](https://regtools.readthedocs.io/en/latest/commands/junctions-extract/).

## NOTE

- Make sure intron junction annotation files are BED formatted (0 based left close right open). This is different from the standard leafcutter.
- leafcutter2 outputs BED formatted coordinates.
- key differences from leafcutter:
    - leafcutter2 use the same set of filters to construct intron clusters.
    - However, when counting junction reads towards predefined or on-demand-run intron clusters, no read filter is applied, essentially all junction reads are counted towards introns.
    - The filtering options (--MINCLUREADS, --MINREADS, --MINCLURATIO) are only used if you are not providing a pre-defined set of intron clusters. They are not used against counting junction reads towards introns.

