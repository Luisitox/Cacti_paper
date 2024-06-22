source ~/anaconda3/etc/profile.d/conda.sh;

### genome assembly
conda activate hifiasm-env
hifiasm -o peyote_z20.asm -t 32 -z20 m64175e_211210_162607.hifi_reads.fq.gz 64175e_211215_124958.hifi_reads.fq.gz

conda deactivate;

### pre-processing of illumina RNAseq fastq files

# Rcorrector
module load perl

RCOR_SCRIPT="~/bin/rcorrector/run_rcorrector.pl"
FASTQ_FOLDER="~/transcriptome_assembly_May21/raw_data"

for file in $FASTQ_FOLDER/*R1.fastq.gz
do
        bn=$(basename $file | cut -d_ -f1,2)
        #echo $bn
        cmd="perl $RCOR_SCRIPT -maxcorK 3 -t 6 -1 $FASTQ_FOLDER/${bn}_R1.fastq.gz -2 $FASTQ_FOLDER/${bn}_R2.fastq.gz"
        cmd_2="bsub -q "new-short" -n 6 -R "rusage[mem=5000]" -J rcorr_${bn} -o rcorr_${bn}_%J.out -e rcorr_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        # eval $cmd_eval
        #sleep 1
done

# TrimGalore
conda activate trimgalore-env

FASTQ_FOLDER="~/transcriptome_assembly_May21/raw_data"

for file in $FASTQ_FOLDER/LW*R1.fastq.gz
do
        bn=$(basename $file | cut -d_ -f1,2)
        #echo $bn
	cmd="trim_galore --paired --retain_unpaired --phred33 --output_dir . --length 36 -q 5 --stringency 1 -e 0.1 -j 6 $FASTQ_FOLDER/${bn}_R1.fastq.gz $FASTQ_FOLDER/${bn}_R2.fastq.gz"
        cmd_2="bsub -q "new-short" -n 6 -R "rusage[mem=3000]" -J trimgalore_${bn} -o trimgalore_${bn}_%J.out -e trimgalore_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        eval $cmd_eval
        #sleep 1
done

conda deactivate;

### mapping of RNAseq fastq files to the genome
conda activate star-env

FASTA_FILE="~/genomics_jan_2022/peyote/assembly/peyote_z20.asm.bp.p_ctg.gfa.fa"
GENOME_DIR="~/stardb/peyote_jan2022"
FASTQ_FOLDER="~/cacti/trueseq_reads/LW"
BAM_FOLDER="~/genomics_jan_2022/peyote/bam_files"

cmd_bsub="bsub -q 'new-short' -n 20 -R 'span[hosts=1] rusage[mem=3000]' -J genome_peyote -o genome_peyote_%J.out -e genome_peyote_%J.err"
cmd_genome="$cmd_bsub STAR --runThreadN 20 --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $FASTA_FILE"

echo $cmd_genome
eval $cmd_genome

for file in $FASTQ_FOLDER/*1.fq.gz
do
        bn=$(basename $file | cut -d_ -f1,2)
        cmd="STAR --genomeDir $GENOME_DIR --readFilesIn $FASTQ_FOLDER/${bn}_R1.cor_val_1.fq.gz $FASTQ_FOLDER/${bn}_R2.cor_val_2.fq.gz --readFilesCommand zcat --runThreadN 20 --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $BAM_FOLDER/${bn}_STAR"
        cmd_2="bsub  -w 'done(genome_peyote)' -q 'new-short' -n 22 -R 'span[hosts=1] rusage[mem=3000]' -J STAR_${bn} -o STAR_${bn}_%J.out -e STAR_${bn}_%J.err"
        cmd_eval="${cmd_2} ${cmd}"
        echo $cmd_eval
        eval $cmd_eval
done

conda deactivate;

### Genome-Independent Trinity-assembly
conda activate trinity-env;

for d in */
do
		bn=$(echo $d | cut -d"/" -f1)
		# echo $bn
		cmd_bsub="bsub -m 'public_2019_hosts public_2017_hosts' -q new-medium -n 24 -R 'span[hosts=1] rusage[mem=10000]' -J trinity_${bn} -o trinity_${bn}.out -e trinity_${bn}.err"
		cmd="Trinity --max_memory 100G --CPU 24 --seqType fq --left ${bn}/*R1*.fq.gz --right ${bn}/*R2*.fq.gz --output ${bn}/trinity_${bn}"
		cmd_eval="${cmd_bsub} ${cmd}"
		echo $cmd_eval
		eval $cmd_eval
done

conda deactivate;

### Trinity-assembly mapping to genome
conda activate minimap2-env

FASTA_FILE="~/genomics_jan_2022/peyote/assembly/peyote_z20.asm.bp.p_ctg.gfa.fa"
TRINITY_FASTQ="~/genomics_jan_2022/peyote/trinity_assembly/LW_trinity.fasta"
BAM_FOLDER="~/genomics_jan_2022/peyote/bam_files"

cmd_bsub="bsub -q 'new-short' -n 30 -R 'span[hosts=1] rusage[mem=2000]' -J minimap_trinity_peyote -o $BAM_FOLDER/minimap_trinity_peyote_%J.out -e $BAM_FOLDER/minimap_trinity_peyote_%J.err"

cmd="$cmd_bsub 'minimap2 -t 30 -ax splice:hq -uf --secondary=no -C5 -O6,24 -B4 $FASTA_FILE $TRINITY_FASTQ > $BAM_FOLDER/peyote_trinity.sam 2> $BAM_FOLDER/peyote_trinity.sam.log'"

echo $cmd
eval $cmd

### Isoseq pre-processing
conda activate isoseq-env

BAM_FILE="~/genomics_jan_2022/isoseq_peyote_salvia_wildpotatoes/cacti-salvia/19133_Aharoni_Asaph/raw_data/20220208/Peyote_RNA/m64175e_220119_170937.hifi_reads.bam"
PRIMERS_FILE="~/genomics_jan_2022/salvia/isoseq_processing/primers.fasta"

bn="Peyote_isoseq"

cmd_bsub="bsub -q 'new-short' -n 30 -R 'span[hosts=1] rusage[mem=2000]' -J lima_${bn} -o lima_${bn}_%J.out -e lima_${bn}_%J.err"
cmd="$cmd_bsub lima --isoseq --dump-clips --peek-guess -j 30 $BAM_FILE $PRIMERS_FILE ${bn}_lima.bam"

echo $cmd
eval $cmd

cmd_bsub="bsub -w 'done(lima_${bn})' -q 'new-short' -n 30 -R 'span[hosts=1] rusage[mem=2000]' -J  refine_${bn} -o refine_${bn}_%J.out -e refine_${bn}_%J.err"
cmd="$cmd_bsub isoseq3 refine -j 30 --require-polya ${bn}_lima.NEB_5p--NEB_Clontech_3p.bam $PRIMERS_FILE ${bn}_refine.bam"

echo $cmd
eval $cmd

cmd_bsub="bsub -w 'done(refine_${bn})' -q 'new-short' -n 30 -R 'span[hosts=1] rusage[mem=2000]' -J  cluster_${bn} -o cluster_${bn}_%J.out -e cluster_${bn}_%J.err"
cmd="$cmd_bsub isoseq3 cluster -j 30 ${bn}_refine.bam ${bn}_polished.bam --verbose --use-qvs"

echo $cmd
eval $cmd

conda deactivate;

### Isoseq mapping to assembly
conda activate isoseq-env

FASTA_FILE="~/genomics_jan_2022/peyote/assembly/peyote_z20.asm.bp.p_ctg.gfa.fa"
MMI_FILE="~/genomics_jan_2022/peyote/assembly/peyote_z20.asm.bp.p_ctg.gfa.mmi"
ISOSEQ_FASTQ="~/genomics_jan_2022/peyote/isoseq_processing/Peyote_isoseq_polished.hq.bam"
BAM_FOLDER="~/genomics_jan_2022/peyote/bam_files"

cmd_bsub="bsub -q 'new-short' -n 10 -R 'span[hosts=1] rusage[mem=2000]' -J pbmm2_index -o pbmm2_index_%J.out -e pbmm2_index_%J.err"

cmd="$cmd_bsub pbmm2 index $FASTA_FILE $MMI_FILE"

echo $cmd
eval $cmd

cmd="$cmd_bsub pbmm2 align $MMI_FILE $ISOSEQ_FASTQ $BAM_FOLDER/peyote_isoseq_polished_pbmm2.bam --preset ISOSEQ --sort -j 30"

cmd_bsub="bsub -w 'done(pbmm2_index)' -q 'new-short' -n 30 -R 'span[hosts=1] rusage[mem=2000]' -J pbmm_peyote_polished -o $BAM_FOLDER/pbmm_peyote_polished_%J.out -e $BAM_FOLDER/pbmm_peyote_polished_%J.err"

cmd="$cmd_bsub pbmm2 align $MMI_FILE $ISOSEQ_FASTQ $BAM_FOLDER/peyote_isoseq_polished_pbmm2.bam --preset ISOSEQ --sort -j 30"

echo $cmd
eval $cmd

### collapse mapped transcripts to get gene annotations

cmd_bsub="bsub -w 'done(pbmm_peyote_polished)' -q 'new-short' -n 30 -R 'span[hosts=1] rusage[mem=2000]' -J peyote_isoseq_polished_pbmm2 -o $BAM_FOLDER/peyote_isoseq_polished_pbmm_collapse_%J.out -e $BAM_FOLDER/peyote_isoseq_polished_pbmm_collapse_%J.err"

cmd="$cmd_bsub isoseq3 collapse $BAM_FOLDER/peyote_isoseq_polished_pbmm2.bam $BAM_FOLDER/peyote_isoseq_polished_pbmm2.bam.gff"

echo $cmd
eval $cmd
