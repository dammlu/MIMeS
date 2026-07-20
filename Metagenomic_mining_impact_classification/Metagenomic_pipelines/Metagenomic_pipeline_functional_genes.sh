Final workflow functional genes from BGR samples SO295


# Workflow DNA MIMes
================================================================================================================================================================================= 

# concatenate single sequencing file to one file using cat file
# make sample.txt

for file in *_R[12].fastq.gz; do
  sample=$(echo "$file" | sed -E 's/_(R[12])?\.fastq\.gz//')
  echo "$sample"
done | sort -u > samples.txt


=================================================================================================================================================================================

# QC using FastQC/MultiQC
# select files like so: fastqc /path/to/sample/*DNA_[1]{8,9}_* -o /path/to/output/

module load bioconda/latest
sbatch -N 1 --account=RPTU-MIMeS --tasks=24 --mem=100000 --time=168:00:00 --mail-type=END --wrap="fastqc /work/RPTU-MIMeS/concat_rawfiles_SO295/*fastq.gz \
-o /work/RPTU-MIMeS/DNA_processing/quality_control/fastqc/"

# multiqc in folder where the .zip files from fastqc are outputted
module load bioconda/latest
multiqc .

=================================================================================================================================================================================

# bbmap quality trimming using bbduk
# here: we do not have any trimming adapters which should be trimmed.therefore: only quality timming

for i in `cat /work/RPTU-MIMeS/concat_rawfiles_SO295/samples.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i})
  mkdir ${SAMPLE}
  module load bioconda/latest
  sbatch -N 1 --account=RPTU-MIMeS --tasks=12 --mem=10000 --time=03:00:00 --wrap="bbduk.sh in1=/work/RPTU-MIMeS/concat_rawfiles_SO295/${SAMPLE}_R1.fastq.gz \
  in2=/work/RPTU-MIMeS/concat_rawfiles_SO295/${SAMPLE}_R2.fastq.gz out1=${SAMPLE}/${SAMPLE}_cleaned_R1.fastq.gz out2=${SAMPLE}/${SAMPLE}_cleaned_R2.fastq.gz \
  threads=24 qtrim=rl trimq=30 minlength=50"
done


=================================================================================================================================================================================

# removal of optical duplicates: clumpify

for i in `cat /work/RPTU-MIMeS/concat_rawfiles_SO295/samples.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i})
  mkdir ${SAMPLE}
  module load bioconda/latest
  sbatch -N 1 --partition=bio-eco --tasks=24 --mem=90GB --time=01:30:00 --wrap="clumpify.sh in=/work/RPTU-MIMeS/DNA_processing/bbduk_q30_len50/${SAMPLE}/${SAMPLE}_cleaned_R1.fastq.gz \
  in2=/work/RPTU-MIMeS/DNA_processing/bbduk_q30_len50/${SAMPLE}/${SAMPLE}_cleaned_R2.fastq.gz out=${SAMPLE}/${SAMPLE}_dedup_R1.fastq.gz out2=${SAMPLE}/${SAMPLE}_dedup_R2.fastq.gz \
  dedupe optical dist=40"
done


=================================================================================================================================================================================

#MEGAHIT assembly 

# paired end mode

for i in `cat /work/RPTU-MIMeS/DNA_processing/Megahit_q30_len50/samples.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i})
  module load bioconda/latest
  sbatch -N 1 --account=RPTU-MIMeS --constraint=SWAP:600 --partition=bio-eco --tasks=12 --mem=90000 --time=15:00:00 --wrap="megahit -1 /work/RPTU-MIMeS/DNA_processing/clumpify_q30_50/${SAMPLE}/${SAMPLE}_dedup_R1.fastq.gz -2 /work/RPTU-MIMeS/DNA_processing/clumpify_q30_50/${SAMPLE}/${SAMPLE}_dedup_R2.fastq.gz -t 32 -o /work/RPTU-MIMeS/DNA_processing/Megahit_q30_len50/${SAMPLE}"
done


=================================================================================================================================================================================

## quast for assembly statistics
# start jobs from QUAST folder!
for i in `cat /work/RPTU-MIMeS/DNA_processing/Megahit_q30_len50/quast_Meg_q30_50/samples1.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i}) 
  echo ${SAMPLE}
  sbatch -N 1 --account=RPTU-MIMeS --partition=bio-eco --constraint=SWAP:600 --tasks=6 --mem=10GB --time=08:00:00 --wrap="/work/RPTU-MIMeS/tools/quast_files/metaquast.py /work/RPTU-MIMeS/DNA_processing/Megahit_q30_len50/${SAMPLE}/final.contigs.fa --max-ref-number 0 -o ${SAMPLE}"
done



##extract quast statitsics using .sh script

#!/bin/bash

# Base directory
BASE_DIR="/work/RPTU-MIMeS/DNA_processing/Megahit/quast"

# Output CSV file
OUTPUT="${BASE_DIR}/contig_summary.csv"

# Write header
echo "Sample,# contigs,# contigs (>= 1000 bp),# contigs (>= 10000 bp),N50" > "$OUTPUT"

# Loop through sample directories
for SAMPLE_DIR in "$BASE_DIR"/SO295_DNA_*; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    REPORT="$SAMPLE_DIR/report.tsv"

    if [[ -f "$REPORT" ]]; then
        CONTIGS=$(awk -F'\t' '$1=="# contigs" {print $2}' "$REPORT")
        CONTIGS_1K=$(awk -F'\t' '$1=="# contigs (>= 1000 bp)" {print $2}' "$REPORT")
        CONTIGS_10K=$(awk -F'\t' '$1=="# contigs (>= 10000 bp)" {print $2}' "$REPORT")
        N50=$(awk -F'\t' '$1=="N50" {print $2}' "$REPORT")

        echo "$SAMPLE,$CONTIGS,$CONTIGS_1K,$CONTIGS_10K,$N50" >> "$OUTPUT"
    else
        echo "Warning: No report.tsv in $SAMPLE_DIR" >&2
    fi
done

=================================================================================================================================================================================

##whokaryote on megahit output with loop for all samples new megahit q30l50

module load bioconda/latest
for i in `cat /work/RPTU-MIMeS/functional_genes_pipeline/whokaryote_output/megahit_new_q30_l50/sampleAO.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i}) 
  mkdir ${SAMPLE}
  cd ${SAMPLE}
rz-singularity --nonv whokaryote_1.1.2--pyhdfd78af_0.simg "whokaryote.py --contigs /work/RPTU-MIMeS/DNA_processing/Megahit_q30_len50/${SAMPLE}/final.contigs.fa --outdir /work/RPTU-MIMeS/functional_genes_pipeline/whokaryote_output/megahit_new_q30_l50/${SAMPLE} --f"
  cd /work/RPTU-MIMeS/functional_genes_pipeline/whokaryote_output/megahit_new_q30_l50
done


=================================================================================================================================================================================

##prodigal on whokaryote (megahit) q30l50

module load bioconda/latest
for i in `cat /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50/samples1.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i}) 
  mkdir ${SAMPLE}
  cd ${SAMPLE}
  sbatch -N 1 --account=RPTU-MIMes --partition=bio-eco --constraint=SWAP:600 --tasks=24 --mem=40000 --time=16:00:00 --wrap="prodigal -i /work/RPTU-MIMeS/functional_genes_pipeline/whokaryote_output/megahit_new_q30_l50/${SAMPLE}/prokaryotes.fasta -a /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50/${SAMPLE}/orfs_${SAMPLE}.faa -o /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50/${SAMPLE}/orfs_${SAMPLE}.gff -d /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50/${SAMPLE}/orfs_${SAMPLE}.fna -p meta"
  cd /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50
done


#linearize orfs

for i in `cat /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50/samples.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i}) 
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < //work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50/${SAMPLE}/orfs_${SAMPLE}.fna | sed 's/\t/\n/g' > /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50/${SAMPLE}/orfs_${SAMPLE}_lin.fna
done

## count ORFs
grep -c "^>" **/*lin.faa > count_orfs.txt

## count complete ORFs
grep -c "^>" **/*orfs_complete_***.faa > count_orfs_complete.txt

##combine.faa files for later annotations
cd /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50
find . -type f -name "*.faa" -exec cat {} + > all_BGR_samples_sortme_notaligned_aminoacid.faa


=================================================================================================================================================================================

#filter out rRNA using SortmeRNA  .sh script
#!/bin/bash

# Load modules once
module load bioconda/latest

# Sample file and base directory
SAMPLES_FILE="/work/RPTU-MIMeS/functional_genes_pipeline/SortmeRNA_q30/samples.txt"
BASE_DIR="/work/RPTU-MIMeS/functional_genes_pipeline/SortmeRNA_q30"

cd "$BASE_DIR"

for SAMPLE in $(cat "$SAMPLES_FILE"); do
  echo "Submitting SortMeRNA job for $SAMPLE"
  mkdir -p "${BASE_DIR}/${SAMPLE}"

  sbatch -N 1 --partition=bio-eco \
    --ntasks=16 --mem=300G --time=96:00:00 \
    --output=${BASE_DIR}/${SAMPLE}/sortme_%j.out \
    --error=${BASE_DIR}/${SAMPLE}/sortme_%j.err \
    --wrap="sortmerna \
    --ref /software/bioconda/sortmerna-data-4.3.6/rRNA_databases/rfam-5.8s-database-id98.fasta \
    --ref /software/bioconda/sortmerna-data-4.3.6/rRNA_databases/rfam-5s-database-id98.fasta --ref /software/bioconda/sortmerna-data-4.3.6/rRNA_databases/silva-arc-16s-id95.fasta \
    --ref /software/bioconda/sortmerna-data-4.3.6/rRNA_databases/silva-arc-23s-id98.fasta --ref /software/bioconda/sortmerna-data-4.3.6/rRNA_databases/silva-bac-16s-id90.fasta \
    --ref /software/bioconda/sortmerna-data-4.3.6/rRNA_databases/silva-bac-23s-id98.fasta --ref /software/bioconda/sortmerna-data-4.3.6/rRNA_databases/silva-euk-18s-id95.fasta \
    --ref /software/bioconda/sortmerna-data-4.3.6/rRNA_databases/silva-euk-28s-id98.fasta \
      --reads /work/RPTU-MIMeS/functional_genes_pipeline/prodigal/megahit_q30_l50/${SAMPLE}/orfs_${SAMPLE}_lin.fna \
      --fastx \
      --other ${BASE_DIR}/${SAMPLE}/sortme_notaligned \
      --aligned ${BASE_DIR}/${SAMPLE}/sortme_aligned \
      -a 16 \
      --workdir ${BASE_DIR}/${SAMPLE}"
done

=================================================================================================================================================================================

##combine sortmeRNA output from all samples and use mmseqs2 to cluster, but keep BGR and GSR samples separated


#BGR
cd /work/RPTU-MIMeS/functional_genes_pipeline/SortmeRNA_q30/comb_count_table_BGR
find . -type f -name "sortme_notaligned.fa" -exec cat {} + > combined_func_catalog/all_BGR_samples_sortme_notaligned.fa



#BGR
module load bioconda/latest
mmseqs easy-cluster /work/RPTU-MIMeS/functional_genes_pipeline/SortmeRNA_q30/comb_count_table_BGR/combined_func_catalog/all_BGR_samples_sortme_notaligned.fa non_redundant_catalog tmp --min-seq-id 0.99 --cov-mode 1 -c 0.9

output: non_redundant_catalog_rep_seq.fasta

#BGR 95%sim
module load bioconda/latest
mmseqs easy-cluster /work/RPTU-MIMeS/functional_genes_pipeline/SortmeRNA_q30/comb_count_table_BGR/combined_func_catalog/cluster_95/all_BGR_samples_sortme_notaligned.fa non_redundant_catalog_95 tmp --min-seq-id 0.95 --cov-mode 1 -c 0.9

output: non_redundant_catalog_95_rep_seq.fasta


=================================================================================================================================================================================

#use kallisto to get read counts, TPM etc 

#kallisto index

#BGR  99% sim
kallisto index --index=STRING /work/RPTU-MIMeS/functional_genes_pipeline/kallisto_new/BGR_func_genes/non_redundant_catalog_rep_seq.fasta --make-unique

#BGR 95% sim
kallisto index --index=STRING /work/RPTU-MIMeS/functional_genes_pipeline/kallisto_new/BGR_func_genes/Kallisto_func_95/non_redundant_catalog_95_rep_seq.fasta --make-unique


## run kallisto using .sh script

#!/bin/bash

#SBATCH -N 1
#SBATCH --job-name=kallisto
#SBATCH --partition=bio-eco    
#SBATCH --ntasks=24
#SBATCH --mem=450G
#SBATCH --time=120:00:00 
#SBATCH --output=kallisto%j.out           
#SBATCH --error=kallisto%j.err                
#SBATCH --mail-type=END

# Define input and output directories
INPUT_DIR="/work/RPTU-MIMeS/DNA_processing/clumpify"
OUTPUT_DIR="/work/RPTU-MIMeS/functional_genes_pipeline/kallisto_new/BGR_func_genes/kallisto_results"
KALLISTO_INDEX="/work/RPTU-MIMeS/functional_genes_pipeline/kallisto_new/BGR_func_genes/kallisto_results/STRING"
SAMPLE_FILE="/work/RPTU-MIMeS/functional_genes_pipeline/kallisto_new/BGR_func_genes/samples_BGR.txt"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Read sample names from sample_file.txt into an array
mapfile -t SAMPLES < "$SAMPLE_FILE"

# Loop through each sample listed in the file
for SAMPLE_NAME in "${SAMPLES[@]}"; do
    # Define input FASTQ files
    SAMPLE_DIR="${INPUT_DIR}/${SAMPLE_NAME}"
    READ1="${SAMPLE_DIR}/${SAMPLE_NAME}_dedup_R1.fastq.gz"
    READ2="${SAMPLE_DIR}/${SAMPLE_NAME}_dedup_R2.fastq.gz"

    # Check if both FASTQ files exist
    if [[ -f "$READ1" && -f "$READ2" ]]; then
        echo "Processing sample: $SAMPLE_NAME"

        # Define sample output directory
        SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE_NAME}"
        mkdir -p "$SAMPLE_OUTPUT"

        # Run Kallisto quant
        kallisto quant --index="$KALLISTO_INDEX" -o "$SAMPLE_OUTPUT" -t 4 "$READ1" "$READ2"

        echo "Finished processing: $SAMPLE_NAME"
    else
        echo "FASTQ files missing for $SAMPLE_NAME, skipping..."
    fi
done

echo "All specified samples processed!"



=================================================================================================================================================================================



#extract TPM values per sample from kallisto result and create sample-to-feature matrix (TPM) values 

#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name=combBGR_fast
#SBATCH --ntasks=8
#SBATCH --partition=bio-eco
#SBATCH --account=RPTU-MIMeS
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --output=combBGR_fast%j.out
#SBATCH --error=combBGR_fast%j.err
#SBATCH --mail-type=END

KALLISTO_DIR="/work/RPTU-MIMeS/functional_genes_pipeline/kallisto_new/BGR_func_genes/kallisto_results"
OUTPUT_FILE="merged_tpm_matrix_BGR.tsv"

# get sample list
samples=($(ls "$KALLISTO_DIR"))
n_samples=${#samples[@]}

echo "Found ${n_samples} samples"

# step 1 — extract target_id + tpm from each sample
# add sample name as column header
mkdir -p tmp_merge

for sample in "${samples[@]}"; do
    abundance_file="${KALLISTO_DIR}/${sample}/abundance.tsv"
    if [[ -f "$abundance_file" ]]; then
        # extract id and tpm, add header
        awk -v s="$sample" \
            'NR==1{print "target_id\t"s} NR>1{print $1"\t"$5}' \
            "$abundance_file" \
            > "tmp_merge/${sample}.tsv"
        echo "Extracted: $sample"
    else
        echo "WARNING: missing $abundance_file"
    fi
done

echo "All samples extracted — now merging with join..."

# step 2 — sort all files by target_id (required for join)
for f in tmp_merge/*.tsv; do
    # sort keeping header on top
    header=$(head -1 "$f")
    tail -n +2 "$f" | sort -k1,1 > "${f}.sorted"
    echo "$header" | cat - "${f}.sorted" > "$f"
    rm "${f}.sorted"
done

echo "All files sorted"

# step 3 — iterative join across all samples
# start with first sample
first_sample="${samples[0]}"
cp "tmp_merge/${first_sample}.tsv" tmp_merge/merged_tmp.tsv

for sample in "${samples[@]:1}"; do
    join \
        -t $'\t' \
        -a 1 -a 2 \
        -e "0" \
        -o auto \
        tmp_merge/merged_tmp.tsv \
        "tmp_merge/${sample}.tsv" \
        > tmp_merge/merged_new.tsv
    mv tmp_merge/merged_new.tsv tmp_merge/merged_tmp.tsv
    echo "Joined: $sample"
done

mv tmp_merge/merged_tmp.tsv "$OUTPUT_FILE"

# cleanup
rm -rf tmp_merge

echo "Done — matrix saved to: $OUTPUT_FILE"
wc -l "$OUTPUT_FILE"

