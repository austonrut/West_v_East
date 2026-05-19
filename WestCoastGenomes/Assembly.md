# Genome Assembly for West Coast Populations
I assembled scaffold level draft genomes for each West Coast population.
I used the fastq files from each indvidual for a given population and used [RedBean](https://github.com/ruanjue/wtdbg2) to assemble

`module load wtdbg`

`wtdbg2 -x ont -g 290m -t 16 -i ../../../fastqs/All_CoosBay.fastq -fo CoosAll.red`

`wtpoa-cns -t16 -i CoosAll.red.ctg.lay.gz -fo CoosAll.red.ctg.fa`
