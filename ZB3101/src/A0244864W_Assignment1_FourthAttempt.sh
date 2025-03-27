#!/bin/bash

# NGS Analysis for ZB3101 Assignment Saccharomyces
# cerevisiae detection of insertion of transposon (Ty5-6p)
# Name: A0244864W E0878509 Choi Jae Bin
# March 12, 2025 Y2024/2025 S2

# Reference/Files used for Data Analysis
Seq_Dir="/Users/saranganddasom/Documents/ZB3101/data"
INPUT_FASTQ_1="$Seq_Dir/A0244864W_1.fq"        # Update with your FASTQ file names
INPUT_FASTQ_2="$Seq_Dir/A0244864W_2.fq"
REF_GENOME="$Seq_Dir/sacCer3.fa"       # Yeast reference genome
TRANSPOSON_REF="$Seq_Dir/ty5_6p.fa"         # Transposon reference
THREADS=4                                         # Number of CPU threads
OUTPUT_DIR="/Users/saranganddasom/Documents/ZB3101/analysis_results"      # Output directory

#Create File Output Directories and Make it organized
#makes subdirectores fastqc, trimmed, alignment, variants, bed_files
mkdir -p $OUTPUT_DIR/{fastqc,trimmed,alignment,bed_files}

# 1. Processing the Data
# 1. A. Quality Control via fastqc
echo -e "Running FastQC on Genome"
fastqc $INPUT_FASTQ_1 $INPUT_FASTQ_2 -o $OUTPUT_DIR/fastqc

# 1. B. Trimming the Sequences Via Trimmomatic
# Threads = 4, Use Fastq1 and Fastq2, output into trimmed data subdirectory
echo -e "Running Trimmomatic to Trim Low Quality Bases"
trimmomatic PE -threads $THREADS \
  $INPUT_FASTQ_1 $INPUT_FASTQ_2 \
  $OUTPUT_DIR/trimmed/trimmed_1.fq.gz $OUTPUT_DIR/trimmed/unpaired_1.fq.gz \
  $OUTPUT_DIR/trimmed/trimmed_2.fq.gz $OUTPUT_DIR/trimmed/unpaired_2.fq.gz \
  ILLUMINACLIP:/opt/homebrew/opt/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# 2. GENOME ALIGNMENT 1 (BOWTIE2)
echo -e "Creating Index of Insert for BowTie"
# 2. A. Produce index for yeast alignment
bowtie2-build $REF_GENOME $OUTPUT_DIR/alignment/yeast_index

# 2. B. Aligning the trimmed data to the reference yeast_index
echo -e "Aligning Fastqc data to Ref index"
bowtie2 -x $OUTPUT_DIR/alignment/yeast_index \
  -1 $OUTPUT_DIR/trimmed/trimmed_1.fq.gz \
  -2 $OUTPUT_DIR/trimmed/trimmed_2.fq.gz \
  -S $OUTPUT_DIR/alignment/aligned.sam \
  --very-sensitive-local\
  --threads $THREADS 

# 2. SAM to BAM
echo -e "Converting SAM TO BAM and sort"
samtools view -bS $OUTPUT_DIR/alignment/aligned.sam > $OUTPUT_DIR/alignment/aligned.bam

# 2. C. Extracting soft clipping
echo -e "Extracting Soft Clippings"
samtools view -h $OUTPUT_DIR/alignment/aligned.bam | \
awk '($6 ~ /S/) {print}' > $OUTPUT_DIR/alignment/softclipped.sam

#add the header back to file then change to BAM
echo -e "Add the heading back intpo the soft clipping file"
samtools view -H $OUTPUT_DIR/alignment/aligned.sam > $OUTPUT_DIR/alignment/header.sam
cat $OUTPUT_DIR/alignment/header.sam $OUTPUT_DIR/alignment/softclipped.sam > $OUTPUT_DIR/alignment/softclipped_withheader.sam
echo -e "Change the softclip file to BAM"
samtools view -bS $OUTPUT_DIR/alignment/softclipped_withheader.sam > $OUTPUT_DIR/alignment/softclipped_withheader.bam
echo -e "Change bam file into fq for next alignment"
samtools fastq $OUTPUT_DIR/alignment/softclipped_withheader.bam > $OUTPUT_DIR/alignment/softclipped.fq

# 3. Align soft clips to transposon index via Bowtie
echo -e "Aligning Soft clipping data to Transposon index via bowtie"

# 3. A. Prepare the transposon database
echo -e "Create index of transposon"
bowtie2-build $TRANSPOSON_REF $OUTPUT_DIR/alignment/transposon_index

# 3. B. Align soft clip to transposon via bowtie
echo -e "Aligning soft clip to transposon"
bowtie2 -x $OUTPUT_DIR/alignment/transposon_index \
        -U $OUTPUT_DIR/alignment/softclipped.fq \
        -S $OUTPUT_DIR/alignment/transposon_alignment.sam \
        --very-sensitive-local \
        --threads $THREADS

# 3. C. Changing Sam to Bam then sorting
echo -e "Changing Sam to Bam"
samtools view -bS $OUTPUT_DIR/alignment/transposon_alignment.sam > $OUTPUT_DIR/alignment/transposon_alignment.bam
echo -e "Sort the confirmed transposon sequences"
samtools sort $OUTPUT_DIR/alignment/transposon_alignment.bam -o $OUTPUT_DIR/alignment/transposon_alignment_sorted.bam
samtools index $OUTPUT_DIR/alignment/transposon_alignment_sorted.bam

# 3. D. Extract reads that match the transposon
echo -e "Extract reads that match transposon and is paired"
samtools view -F 4 $OUTPUT_DIR/alignment/transposon_alignment_sorted.bam | \
awk '{print $1}' | sort | uniq > $OUTPUT_DIR/alignment/transposon_matching_reads.txt

# 3. E. Extract the matching transposon sequences located in the initial trim
echo -e "Aligning transposon matching reads back to the genome"
seqtk subseq $OUTPUT_DIR/trimmed/trimmed_1.fq.gz $OUTPUT_DIR/alignment/transposon_matching_reads.txt > $OUTPUT_DIR/alignment/matching_1.fq
seqtk subseq $OUTPUT_DIR/trimmed/trimmed_2.fq.gz $OUTPUT_DIR/alignment/transposon_matching_reads.txt > $OUTPUT_DIR/alignment/matching_2.fq

# 3. F. Ensure that the reads have pairs to each other, are matching
# Why dont you workkkkkkkkk.......
# This ensures that we can performed pair-end reads
#Read number of reads for each
reads_1=$(tr -d '\r' < $OUTPUT_DIR/alignment/matching_1.fq | grep -c "^@")
reads_2=$(tr -d '\r' < $OUTPUT_DIR/alignment/matching_2.fq | grep -c "^@")
# Print the number of reads in each file to check
echo "Number of reads in matching_1.fq: $reads_1"
echo "Number of reads in matching_2.fq: $reads_2"

# 3. F. Step 1: Extract read IDs from both FASTQ files and normalize them (remove trailing /1 or /2)
grep "^@" $OUTPUT_DIR/alignment/matching_1.fq | sed 's/\/1//; s/ .*//' | sort > $OUTPUT_DIR/alignment/r1_ids.txt
grep "^@" $OUTPUT_DIR/alignment/matching_2.fq | sed 's/\/2//; s/ .*//' | sort > $OUTPUT_DIR/alignment/r2_ids.txt

# 3. F. Step 2: Find the common read IDs between the two files
comm -12 $OUTPUT_DIR/alignment/r1_ids.txt $OUTPUT_DIR/alignment/r2_ids.txt > $OUTPUT_DIR/alignment/common_ids.txt

# 3. F. Step 3: Add the /1 /2 to compensate for them existing in matchings (forward and reverse)
# Append /1 to common IDs for R1
sed 's/$/\/1/' $OUTPUT_DIR/alignment/common_ids.txt > $OUTPUT_DIR/alignment/common_ids_r1.txt
# Append /2 to common IDs for R2
sed 's/$/\/2/' $OUTPUT_DIR/alignment/common_ids.txt > $OUTPUT_DIR/alignment/common_ids_r2.txt

# 3. G. Extract the matching end_pairs and create a new file
# For matching_1_common_fixed.fq
grep -A 3 -Ff $OUTPUT_DIR/alignment/common_ids_r1.txt $OUTPUT_DIR/alignment/matching_1.fq | grep -v '^--$' > $OUTPUT_DIR/alignment/matching_1_common_fixed.fq
# For matching_2_common_fixed.fq
grep -A 3 -Ff $OUTPUT_DIR/alignment/common_ids_r2.txt $OUTPUT_DIR/alignment/matching_2.fq | grep -v '^--$' > $OUTPUT_DIR/alignment/matching_2_common_fixed.fq
# ============== Now we finally have 2 pair end matchings so file has same number of sequences==========

# 4. Perform alignment of mapped insert sequences against yeast_genome
echo -e "Final alignment of transposon to yeast_genome"
bowtie2 -x $OUTPUT_DIR/alignment/yeast_index \
        -1 $OUTPUT_DIR/alignment/matching_1_common_fixed.fq \
        -2 $OUTPUT_DIR/alignment/matching_2_common_fixed.fq \
        -S $OUTPUT_DIR/alignment/insert_genome_final_alignment.sam \
        --very-sensitive-local \
        --threads $THREADS

# 4. A. Convert to BAM
echo -e "Converting SAM to BAM"
samtools view -bS $OUTPUT_DIR/alignment/insert_genome_final_alignment.sam > $OUTPUT_DIR/alignment/insert_genome_final_alignment.bam
samtools sort $OUTPUT_DIR/alignment/insert_genome_final_alignment.bam -o $OUTPUT_DIR/alignment/transposon_alignment_sorted.bam
samtools rmdup $OUTPUT_DIR/alignment/transposon_alignment_sorted.bam $OUTPUT_DIR/alignment/transposon_alignment_dedup.bam
samtools view -b -f 2 "$OUTPUT_DIR/alignment/transposon_alignment_dedup.bam" > "$OUTPUT_DIR/alignment/properly_mapped.bam"
samtools index $OUTPUT_DIR/alignment/properly_mapped.bam


#5 Getting the answers

#5 A. # Extract chromosome (RNAME), position (POS), and orientation (FLAG) from the aligned BAM file
echo -e "Creating BED file with insert locations and orientations"

samtools view "$OUTPUT_DIR/alignment/properly_mapped.bam" | while read -r line; do
    chrom=$(echo "$line" | cut -f3)
    pos=$(echo "$line" | cut -f4)
    flag=$(echo "$line" | cut -f2)

    # Skip unmapped reads
    if [[ "$chrom" == "*" || "$pos" -le 0 ]]; then
        continue
    fi
    # Correct strand detection using bitwise AND
    if (( $flag & 16 )); then
        strand="-"
    else
        strand="+"
    fi
    # Convert 1-based POS to 0-based start for BED format
    start=$((pos - 1))
    end=$((start + 1))

    echo -e "$chrom\t$start\t$end\tinsertion_$((++count))\t$flag\t$strand"
done > "$OUTPUT_DIR/bed_files/transposon_insertions.bed"

# Step 2: Merge overlapping or nearby insertions
bedtools merge -i $OUTPUT_DIR/bed_files/transposon_insertions.bed \
               -c 4,5,6 \
               -o distinct,distinct,distinct \
               > $OUTPUT_DIR/bed_files/transposon_insertions_merged.bed

# Step 3: Count the number of unique insertion sites
num_unique_inserts=$(wc -l < $OUTPUT_DIR/bed_files/transposon_insertions_merged.bed)
echo -e "Total number of unique insertion sites: $num_unique_inserts"