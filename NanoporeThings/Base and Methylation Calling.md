# Processing Raw Nanopore Sequence Data
Here is my code for going from raw sequence outputs from Nanopore to basecalled and methylation sequeences

## Basecalling
`module load guppy/6.5.7-cuda11.4
guppy_basecaller -x "auto" -i ./pod5_al/barcode01 -s ./basecall.out/barcode01 -c dna_r10.4.1_e8.2_400bps_hac_prom.cfg --min_qscore 7 --num_callers 1`

## Methylation Calling
### Pre-processing
Some files needed to be quaility checked and/or converted to the correct file format before methylation calling

#### POD5 to blow5
POD5 output from Promethion need to be converted to slow5/blow5 files. I used [blue-crab](https://github.com/Psy-Fer/blue-crab)
Used -o flag for outputing all in single file. 

`source /scratch/arutle14/blue-crab-venv/bin/activate
blue-crab p2s /scratch/arutle14/West_v_East/RawData/WestRaw/pod5_al/barcode02 -o ../methylation/barcode02.meth/barcode02.all.blow5`

#### Fast5 to blow5
Fast5 output needed to be converted to slow5/blow5 files and then merged. Used [slow5tools](https://github.com/hasindu2008/slow5tools)
`./slow5tools-v1.2.0/slow5tools f2s ../../RawData/EastRaw/fast5_all/barcode01/ -d ./barcode01.blow5`
`./slow5tools-v1.2.0/slow5tools merge ./barcode01.blow5 -o bc01.all.blow5`

#### Quaility Check
`slow5tools/slow5tools quickcheck XXXX.blow5`

### Call Methylation
I used [f5c](https://github.com/hasindu2008/f5c) for calling methylation. First,  fastq reads need to be indexed... 
`module load f5c
f5c index --slow5 ../slow5/bc01.all.blow5 ./barcode01.meth/NH_1_bc01.fastq`

...and then aligned to the genome. 
`module load minimap2
module load samtools
minimap2 -a -x ava-ont uk_jaNemVect1.1_genomic.fna barcode01.meth/CoosBay_A_bc01.fastq | samtools sort -T tmp -o merge.test1.ukbc01.sorted.bam
samtools index merge.test1.ukbc01.sorted.bam`

Then can call methylation and get methylation frequency
`f5c call-methylation --slow5 ../../slow5/bc02.all.blow5 -b merge.NH2.ukbc02.sorted.bam -g ../../../West/methylation/uk_jaNemVect1.1_genomic.fna -r NH_2_bc02.fastq > ./NH2_bc02_methcall.tsv
f5c meth-freq -i barcode02.meth/NH2_bc02_methcall.tsv > barcode02.meth/NH2_bc02_methfreq.tsv`
