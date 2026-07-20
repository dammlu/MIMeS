Final workflow taxonomic marker genes (16S) from BGR samples SO295



 
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

##use bbmap to separate 16S and 18S genes (map against SILVA v138.1 database)

#final 16S
for i in `cat /work/RPTU-MIMeS/taxonomic_genes_pipeline/bbmap_16S/samples.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i})
  mkdir ${SAMPLE}
  cd ${SAMPLE}
module load bioconda/latest
sbatch -N 1 --tasks=16 --mem=100000 --time=24:00:00 --wrap="bbmap.sh ref=/software/bioconda/sortmerna-data-4.3.6/rRNA_databases/silva-bac-16s-id90.fasta \
         in=/work/RPTU-MIMeS/DNA_processing/clumpify_q30_50/${SAMPLE}/${SAMPLE}_dedup_R1.fastq.gz in2=/work/RPTU-MIMeS/DNA_processing/clumpify_q30_50/${SAMPLE}/${SAMPLE}_dedup_R2.fastq.gz \
         outm=aligned_reads.fastq.gz \
         outu=nonaligned_reads.fastq \
         minid=0.85 threads=16"
cd /work/RPTU-MIMeS/taxonomic_genes_pipeline/bbmap_16S
done



=================================================================================================================================================================================

##pass /bbmap/final16S/${SAMPLE}/aligned_reads.fastq.gz  to the phyloflash pipeline

###Phyloflash
#it finally worked! M.Hillenbrand changed storage location of phyloflash and its dependencies to bioconda!


#test phyloflash loop --> worked! (but not as slurm job) 
for i in `cat /work/RPTU-MIMeS/taxonomic_genes_pipeline/phyloFlash/samples.txt`
do
  echo ${i}
  SAMPLE=$(echo ${i})
  mkdir ${SAMPLE}
  cd ${SAMPLE}
module load usearch/latest
module load bioconda/latest
source /software/bioconda/latest/etc/profile.d/conda.sh
conda activate phyloflash
phyloFlash.pl -lib ${SAMPLE}_16S_emirge_fastq -emirge -dbhome /groups/RPTU-MIMeS/phyloFlash/138.1 -read1 /work/RPTU-MIMeS/taxonomic_genes_pipeline/bbmap_16S/${SAMPLE}/aligned_reads.fastq.gz -interleaved
conda deactivate
cd /work/RPTU-MIMeS/taxonomic_genes_pipeline/phyloFlash
done


#get number of assembled genes
grep -c ">" filename.fasta

#get mean length of assembled genes
awk '/^>/ {if (seq) {len+=length(seq); count++; seq=""}} !/^>/ {seq=seq$0} END {len+=length(seq); count++; print len/count}' genes.fasta



=================================================================================================================================================================================

##extract complete 16S tax marker genes info using .sh script

#!/bin/bash
#SBATCH --job-name=extract16Sinfo
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --account=RPTU-MIMeS
#SBATCH --mail-type=END


# Output CSV file
output_csv="summary_16S_marker_genes.csv"

# Write the header for the CSV file
echo "Sample,Num_Marker_Genes,Mean_Length" > $output_csv

# Loop through all directories matching the pattern SO295_DNA_*
for dir in SO295_DNA_*; do
    # Check if the directory exists and contains the target file
    fasta_file="$dir/${dir}_16S_emirge_fastq.SSU.collection.fasta"
    if [[ -f $fasta_file ]]; then
        echo "Processing $fasta_file"

        # Count the number of sequences (marker genes) in the FASTA file
        num_marker_genes=$(grep -c '^>' "$fasta_file")

        # Calculate the mean sequence length
        mean_length=$(awk '/^>/ { if (seq) { print length(seq); seq="" } next } { seq=seq $0 } END { if (seq) print length(seq) }' "$fasta_file" | \
                      awk '{ sum += $1; count++ } END { if (count > 0) print sum / count; else print 0 }')

        # Append the results to the CSV file
        echo "$dir,$num_marker_genes,$mean_length" >> $output_csv
    else
        echo "Warning: $fasta_file not found!"
    fi
done

echo "Summary saved to $output_csv"


=================================================================================================================================================================================

#concatenate all 16S marker genes (.fasta) in one file -> create gene catalogue

cat /work/RPTU-MIMeS/taxonomic_genes_pipeline/phyloFlash_BGR/*/*_16S_emirge_fastq.SSU.collection.fasta > new_combined_16S_marker_genes_BGR.fasta
cat /work/RPTU-MIMeS/taxonomic_genes_pipeline/phyloFlash_GSR/*/*_16S_emirge_fastq.SSU.collection.fasta > new_combined_16S_marker_genes_GSR.fasta


#clustering at 97% sim
#BGR
module load bioconda/latest
mmseqs easy-cluster new_combined_16S_marker_genes_BGR.fasta combined_16S_cluster_BGR_97 tmp --min-seq-id 0.97 --cov-mode 1 -c 0.9


#GSR
module load bioconda/latest
mmseqs easy-cluster new_combined_16S_marker_genes_GSR.fasta combined_16S_cluster_GSR_97 tmp --min-seq-id 0.97 --cov-mode 1 -c 0.9


output: /work/RPTU-MIMeS/taxonomic_genes_pipeline/mmseqs/combined_16S_cluster_BGR_97_rep_seq.fasta


=================================================================================================================================================================================

#use kallisto to get read counts, TPM etc 

#kallisto index

kallisto index --index=STRING /work/RPTU-MIMeS/taxonomic_genes_pipeline/Kallisto/taxo_BGR/taxo_BGR_97/combined_16S_cluster_BGR_97_rep_seq.fasta --make-unique
kallisto index --index=STRING /work/RPTU-MIMeS/taxonomic_genes_pipeline/Kallisto/taxo_GSR/taxo_GSR_97/combined_16S_cluster_GSR_97_rep_seq.fasta --make-unique



## run kallisto using .sh script

#!/bin/bash

#SBATCH -N 1
#SBATCH --job-name=kallisto
#SBATCH --partition=bio-eco    
#SBATCH --ntasks=24
#SBATCH --mem=200G
#SBATCH --time=96:00:00 
#SBATCH --output=kallisto%j.out           
#SBATCH --error=kallisto%j.err                
#SBATCH --mail-type=END

# Define input and output directories
INPUT_DIR="/work/RPTU-MIMeS/DNA_processing/clumpify"
OUTPUT_DIR="/work/RPTU-MIMeS/First_Paper_SML/Kallisto/kallisto_tax/kallisto_results"
KALLISTO_INDEX="/work/RPTU-MIMeS/First_Paper_SML/Kallisto/kallisto_tax/STRING"
SAMPLE_FILE="/work/RPTU-MIMeS/First_Paper_SML/Kallisto/kallisto_tax/sample_file.txt"

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

extract TPM values per sample from kallisto result and create sample-to-feature matrix (TPM) values 

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

KALLISTO_DIR="/work/RPTU-MIMeS/First_Paper_SML/Kallisto/kallisto_tax/kallisto_results"
OUTPUT_FILE="merged_tpm_matrix_BGR_tax.tsv"

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
