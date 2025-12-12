# Yleaf v3.2 (Enhanced Batch Processing & T2T Support)

**This is a modified version of Yleaf featuring robust batch processing, T2T reference support, and crash protection.**

*   **Original Authors:** Arwin Ralf, Diego Montiel Gonzalez, Kaiyin Zhong, Manfred Kayser (Erasmus MC)
*   **Enhanced By:** [Your Name/Organization]
*   **Key Features:** 
    *   **Batch Processing:** Process hundreds of BAM/CRAM files automatically (`batch_process.py`).
    *   **Fast Fail:** Automatically detects and skips corrupted files or broken indices without hanging.
    *   **Auto-Healing:** Attempts to fix missing or outdated indices (`samtools index`) on the fly.
    *   **T2T Support:** Full support for the T2T-CHM13 reference genome.

---

# Original Yleaf Description

#### Arwin Ralf, Diego Montiel Gonzalez, Kaiyin Zhong and Manfred Kayser

### Department of Genetic Identification 
#### Erasmus MC University Medical Centre Rotterdam, The Netherlands

## Requirements

    Operating system: Linux only. 
    Internet connection: when running for the first time for downloading the reference genome. Alternatively you 
                         can configure your own references.
    Data storage: For installation we recommend a storage capacity of > 8 GB. 

## Installation

The easiest way to get Yleaf up and running is by using a conda environment. 

```bash
# first clone this repository to get the environment_yleaf.yaml
git clone https://github.com/genid/Yleaf.git
cd Yleaf
# create the conda environment from the .yaml the environment will be called yleaf
conda env create --file environment_yleaf.yaml
# activate the environment
conda activate yleaf
# pip install the cloned yleaf into your environment. Using the -e flag allows you to modify the config file in your cloned folder
pip install -e .

# verify that Yleaf is installed correctly. You can call this command from any directory on your system
Yleaf -h 
```      
or manually install everything
```bash
# install python and libraries
apt-get install python3.6
pip3 install pandas
pip3 install numpy
# install Burrows-Wheeler Aligner for FASTQ files
sudo apt-get install minimap2 
# install SAMtools
wget https://github.com/samtools/samtools/releases/download/1.4.1/
samtools-1.4.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2
cd samtools-1.4.1/
./configure make
make install
# clone the yleaf repository
git clone https://github.com/genid/Yleaf.git
# pip install the yleaf repository
cd Yleaf
pip install -e .

# verify that Yleaf is installed correctly. You can call this command from any directory on your system
Yleaf -h 
```
After installation you can navigate to the yleaf/config.txt folder and add custom paths for the files listed there. This will make sure that Yleaf does not download the files on the first go or downloads the files in the provided location. This allows you to use a custom reference if you want. Please keep in mind that custom reference files might cause other issues or give problems in combination with already existing data files. Positions are based on either hg38 or hg19.

## Usage and examples
Here follow some minimal working examples of how to use Yleaf with different input files. There are additional options
that can be used to tune how strict Yleaf is as well as options to get private mutations as well as a graph showing 
the positioning of predicted haplogroups of all your samples in the Haplogroup tree.

_Note: In version 3.0 we switched to using YFull (v10.01) for the underlying tree structure of the haplogroups.
 This also means that predictions are a bit different compared to earlier versions._
### Yleaf: FASTQ (raw reads)

    Yleaf -fastq raw_reads.fastq -o fastq_output --reference_genome hg38
        
### Yleaf: BAM or CRAM format
    Yleaf -bam file.bam -o bam_output --reference_genome hg19 
    Yleaf -cram file.bam -o cram_output --reference_genome hg38 

### With drawing predicted haplogroups in a tree and showing all private mutations

    Yleaf -bam file.bam -o bam_output --reference_genome hg19 -dh -p

### Batch Processing (New in v3.2)

Yleaf now includes a robust batch processing script `batch_process.py` located in the root directory. This tool is designed to handle large collections of BAM/CRAM files efficiently.

**Features:**
*   **Directory Scanning:** Recursively finds all `.bam` and `.cram` files in a folder.
*   **Multi-threading:** Process multiple files in parallel using the `-t` argument.
*   **Fast Fail & Auto-Healing:** Automatically checks for valid indices (`.bai`/`.crai`).
    *   If an index is missing or invalid, it attempts to generate one using `samtools index` (with a timeout).
    *   If indexing fails or the file is corrupted, it skips the file immediately to prevent hanging.
    *   Supports auto-detection of CRAM references to prevent hangs on older samtools versions.
*   **Improved Reporting:** Detects samples with very low Y-DNA coverage/quality (e.g., female samples) and flags them as `Low_Y_Signal` instead of generic failure.
*   **Summary Reporting:** Generates a `summary_table.csv` with results for all samples (Haplogroup, QC Score, Status).
*   **Resume Capability:** Logs are saved per sample; you can re-run the script to retry failed samples or process new ones.

**Usage:**

```bash
# Process all BAM/CRAM files in a directory (recursively) with 16 threads
python3 batch_process.py /path/to/your/bam_folder -d -o output_directory -t 16

# Process a list of files from a text file
python3 batch_process.py file_list.txt -o output_directory
```

**Output:**
*   `summary_table.csv`: A CSV file containing the haplogroup prediction, quality scores, and links to YFull for each processed sample.
*   `output_directory/`: Contains individual result folders for each sample.

### NAS / Docker Indexing Helper (`smart_index.sh`)

If your data is stored on a NAS (e.g., TrueNAS Scale) or remote server, indexing files over the network (SMB/NFS) can be slow. Use the `smart_index.sh` script to fix indices directly on the storage server.

1.  Copy `smart_index.sh` to your data folder on the NAS.
2.  Edit the reference genome paths in the script (`REF_T2T`, etc.) to match the NAS file structure.
3.  Run it inside the NAS shell/container:
    ```bash
    chmod +x smart_index.sh
    ./smart_index.sh .
    ```
4.  Once finished, run `batch_process.py` on your PC. It will detect the valid indices and skip re-indexing.

### Optimization Recommendation (WGSExtract / Samtools)

Processing full Whole Genome Sequencing (WGS) BAM files (often 30GB - 100GB+) can be slow. Yleaf only needs reads mapped to the Y chromosome.

**Recommendation:** Extract the Y chromosome reads before running Yleaf. This reduces file size to ~50-200MB and reduces processing time from hours to seconds per sample.

**Linux (samtools):**
```bash
samtools view -b input.bam chrY > input_chrY.bam
samtools index input_chrY.bam
```
*Note: Check if your reference uses 'chrY' or 'Y'.*

**Windows (WGSExtract):**
If you are preparing files on Windows, use **WGSExtract** (Beta) to extract the Y chromosome into a smaller BAM file, then transfer it to your Linux environment for Yleaf.

### Known Issues & Modifications (v3.2 Enhanced)

1.  **Marker G-PF3346 Removal:**
    *   The marker `G-PF3346` (pos 23141720 in T2T) has been observed to show ancestral state (False Negative) in some high-quality WGS samples, blocking detection of deeper clades like `G-L1264`.
    *   **Action:** This marker was removed from `yleaf/data/t2t/new_positions.txt`.
    *   **Note:** If updating the database, ensure this marker remains excluded or verified.

2.  **Base Majority Threshold:**
    *   The default `base_majority` threshold has been lowered from **90% to 75%**.
    *   **Reason:** WGS data often has moderate coverage (e.g. 8x). A single sequencing error (1 vs 7 reads = 87.5%) would previously cause valid markers to be rejected. 75% allows for 1 error in 4 reads.

## Additional information

For a more comprehensive manual please have a look at the [yleaf_manual](yleaf_manual.pdf).

If you have a bug to report or a question about installation consider sending an email to 
 a.ralf at erasmusmc.nl or create an issue on GitHub.


### References and Supporting Information
A. Ralf, et al., Yleaf: software for human Y-chromosomal haplogroup inference from next generation sequencing data (2018).

https://academic.oup.com/mbe/article/35/5/1291/4922696

