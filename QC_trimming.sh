# Current directory is /home/user/Documents/wd

# Make sure you have installed fastqc and multiqc using:
# sudo apt install fastqc
# sudo apt install multiqc

##### Merge fastq files from different lanes by running './fq.files/fastq_lane_merging.sh' script
# After merging, move the resulting files to a new folder for better organization
mkdir data
mv ./fq.files/*_ME_* ./data

##### Perform quality control (QC) before trimming
# Replace '/home/user/Documents/wd' with your actual working directory in the code below
mkdir /home/user/Documents/wd/QC/fastqc_results

# Run fastqc in parallel on all *.fastq.gz files, outputting results to QC/fastqc_results
find *.fastq.gz | parallel -j 24 "fastqc {} --outdir /home/user/Documents/wd/QC/fastqc_results"
cd /home/shego/Documents/wd/QC/fastqc_results

# Generate a summary report of the fastqc results using multiqc
multiqc .

##### Trim both adapters (since we are uncertain of the adapters used)
# Replace '/home/user/Documents/wd' with your actual working directory in the code below
# But since Nextera adapters are suspected, run trimming for those instead
mkdir /home/user/Documents/wd/trimmed_data_illu

# Nextera adapter sequences (CTGTCTCTTATACACATCT, AGATGTGTATAAGAGACAG) will be trimmed in the next section
mkdir /home/user/Documents/wd/trimmed_data

# Trim Nextera adapters from the raw fastq files
for i in `ls -1 /home/user/Documents/wd/data/*_ME_L001_R1_001.fastq.gz`; \
    do dname=$(dirname ${i}); name=$(basename ${i} _ME_L001_R1_001.fastq.gz); \
    bbduk.sh -Xmx100g in1=${dname}/${name}_ME_L001_R1_001.fastq.gz in2=${dname}/${name}_ME_L001_R2_001.fastq.gz \
    out1=/home/user/Documents/wd/trimmed_data/${name}_trimmed_1.fq.gz out2=/home/user/Documents/wd/trimmed_data/${name}_trimmed_2.fq.gz \
    qin=33 literal=CTGTCTCTTATACACATCT,AGATGTGTATAAGAGACAG restrictright=20 rcomp=f k=19 mink=8 ktrim=r hdist=1 hdist2=0 qtrim=rl trimq=20 maq=20 minlen=25 >> QC/trim_log 2>&1 ; 
done

##### Perform quality control (QC) on the trimmed data
# Replace '/home/user/Documents/wd' with your actual working directory in the code below
mkdir /home/user/Documents/wd/QC/trimmed_fastqc_results

# Run fastqc on the trimmed fastq files and output results to QC/trimmed_fastqc_results
find *.fq.gz | parallel -j 24 "fastqc {} --outdir /home/user/Documents/wd/QC/trimmed_fastqc_results"
cd /home/user/Documents/wd/QC/trimmed_fastqc_results

# Generate a summary report of the trimmed fastqc results using multiqc
multiqc .

