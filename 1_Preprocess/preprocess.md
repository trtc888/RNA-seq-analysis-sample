# Data preprocessing of breast cancer RNA req raw results
RNA seq results are usually saved as .fastq formats. They will need to be cleaned up and aligned before we can analyze them. An example of a fastq block is shown below.
```
@A01619:77:HCLWCDSX3:3:1101:10158:1235 1:N:0:TCGGATTC+ACAAGCTC
CTCGTCTGTTGGGGTGGTGCTCAGGGTAAAGGGCTATGGGCAACAGGGGACCAGACCAGGGATGAGTGGGGAGGGCACAAGGACCATTTGCCA
+
!''*((((***+))%%%++)(%%%%).1***-+*''5CCF>>>>>>CCCCCCC65:FFFFFF:FFFFFFFFFF:FFFFFFFFFFFFFFFFFFF
```
In the block above:
- The first line following @ is the identifier and the optional description.
- the second line is the sequence.
- The thrid line is the "+" followed by any optional identifier and description.
- the fourth line is the quality values of the matched sequence code. The quality value is usually the ASCII code of teh quality value plus 33 or 64 and displayed as a ASCII symbol, namely phred33 or phred64, which will be used later on for trimming the sequence. 

In this document, raw files of RNA seq from the paper below will be used.

> Publication Pommerenke C, Nagel S, Haake J, Koelz AL, Christgen M, Steenpass L, Eberth S. (2024). Molecular Characterization and Subtyping of Breast Cancer Cell Lines Provide Novel Insights into Cancer Relevant Genes. Volume 13. Issue 4.
pmid 38391914 doi 10.3390/cells13040301

## Step 1: download the files.
1. Before downloading, create a working folder and subfolders. The subfolders I created are called "00ref", "01raw", "02clean","03bam","04counts".
2. The raw data are stored at <https://www.ebi.ac.uk/biostudies/studies/S-BSST1200?query=S-BSST1200>. Select the interesting cell types and follow the guide to download. Each cell type contains two sequencing documents. Both need to be dowloaded. In the current setting, the downloaded files are kept under "./01raw/" folder. In current case, files for MDA-MB-231, MDA-MB-453, and MDA-MB-468 cell lines are downloaded.
3. Download the reference files. For simplicity, directly download the pre-compiled index file for later Hisat2 use. The file is kept down "00ref" foler.
```sh
$ wget ftp://anonymous@ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
```
4. Download file*Homo_sapiens.GRCh38.113.gtf.gz* from <https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/> to "00ref", and unzip it.
5. Download the method package of the paper from <https://zenodo.org/records/6401600> and store the adapter file *illu_ad_sel.fa* under "00ref" folder.

## Step 2: Quality check using FastQC MultiQC
Once all raw files are kept in the "01raw" folder, run FastQC and MultiQC sequentially.
```
$ fastqc *
```
then
```
$ multiqc
```
FastQC will check the quality of the fastq files and generate the report. MultiQC will aggregte the reports for easier reading.

## Step 3: Trim the sequence using Trimmomatic.
This step is to remove adapter sequence before further processing using Trimmomatic. Use the command below under the "01raw" folder.
```sh
$ java -jar /Path/To/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 NG-29643_MDA_MB_231_lib581314_8005_3_1.fastq.gz NG-29643_MDA_MB_231_lib581314_8005_3_2.fastq.gz /Path/To/breast_cancer/02clean/231_1_paired.fq.gz /Path/To/breast_cancer/02clean/231_1_unpaired.fq.gz /Path/To/02clean/231_2_paired.fq.gz /Path/To/breast_cancer/02clean/231_2_unpaired.fq.gz -summary /Path/To/breast_cancer/02clean/231_summary.txt ILLUMINACLIP:/Path/To/breast_cancer/05paper_method/illu_ad_sel.fa:2:30:10:1:True LEADING:3 TRAILING:3 MINLEN:36 -phred33
```

A brief explain of this long command:
- java -jar /Path/To/Trimmomatic-0.39/trimmomatic-0.39.jar: call the trimmomatic program downloaded. Make sure to change "/Path/To" to the real directory.
- PE: let the program know that this is paired mode. Another option is SE which will not be used here.
- -threads 4: tell the program to use 4 cores. Change the number based on the computer capacity.
- NG-29643_MDA_MB_231_lib581314_8005_3_1.fastq.gz NG-29643_MDA_MB_231_lib581314_8005_3_2.fastq.gz: Two input files. If the command is not run under "01raw" folder, add path before the file names.
- Path/To/breast_cancer/02clean/231_1_paired.fq.gz /Path/To/breast_cancer/02clean/231_1_unpaired.fq.gz /Path/To/02clean/231_2_paired.fq.gz /Path/To/breast_cancer/02clean/231_2_unpaired.fq.gz: four output files with path. Make sure the order is correct.
- -summary /Path/To/breast_cancer/02clean/231_summary.txt: this saves the summary after the trimming is down.
- ILLUMINACLIP:/Path/To/breast_cancer/00ref/illu_ad_sel.fa:2:30:10:1:True LEADING:3 TRAILING:3 MINLEN:36 -phred33: just keep them. Don't have to change the numbers at beginning. 

Once the trimming is done, there should be four output files in the "02clean" folder.

Trim all raw files before next step.

## Step 4: generate alignment file using Hisat2
Because we are using pre-compiled index file, we can directly align the results. Otherwise we should download the reference sequence and create the index file by ourselves.
1. Unzip the *grch38.tar.gz* file under "00ref" folder. The command below can be used.
```
$ tar -xvzf grch38.tar.gz
```
2. Use Hisat2 to aligned the file one by one. Below is an example.
```
$ hisat2 --new-summary --thread 4 -x 00ref_hisat/grch38/genome -1 02clean/231_1_paired.fq -2 02clean/231_2_paired.fq 2>231.ht2.log | samtools sort > 03aligned/231.bam
```
A brief explantion:
- --new-summary: this tells the program what it is going to do.
- --thread 4: use four cores.
- -x 00ref_hisat/grch38/genome: location of the reference file. There will be multiple files starting with "genome".
- -1 02clean/231_1_paired.fq -2 02clean/231_2_paired.fq 2>231.ht2.log: two input files and generate a log file.
- | samtools sort > 03bam/231.bam: pass the results to samtools to make bam file.

After this command, the bam file should be generated.

## Step 5: Count the readings using FeatureCounts
Use the code below to generate the count file. Note that gene_id will be used for each count.
```
$ featureCounts -a 00ref/Homo_sapiens.GRCh38.113.gtf -F GTF -g gene_id -p -T 4 -o 04counts/expression.txt 03bam/*.bam
```
After that, a count file *expression.txt" is generated. It can be analyzed using DESeq2 package with R, which will be discussed next.