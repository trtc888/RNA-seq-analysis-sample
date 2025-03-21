# Computer configurration

In this part, I am going to briefly describe how to configure the working environment of the computer for RNA seq analysis.

There are many ways to configure the computer, a workibng laptop, an old computer, a server, a node of HPC, or a cloud server, such as an Amazon EC2. The configurations on different operation systems also varies. Users may have their preferences in which program to use, which will also impact what hardware and software to choose. For example, if you use STAR for sample preprocess, you may need 32 GB RAM with numerous cores for human transcriptome. For bench scientists who are not performing RNA seq analysis everyday, a more convenient environment may be favored.

Here, I will describe two ways of setting up a working environment, local and remote computer.

If you happen to have an old computer, either a desktop or a laptop, with multiple cores on your CPU, and about 8 GB or more RAM, which is pretty common with 10-year-old computers, you may use it specifically to practice the RNA-seq analysis. If you prefer to use your current working computer, you may use it too.

## Step 1. System preparation (optional)
### 1) Install the Ubuntu system locally
- For a dedicated computer, download the Ubuntu (better LTS) iamge from <https://ubuntu.com/download/desktop#how-to-install-NobleNumbat>, and burn it into a 8 GB (or larger) USB drive following the offical instructions. This USB drive will be your installation media
- Boot your computer from the USB drive and install a clean system following its instructions.

### 2) Use a remote system
- If you need to set up a cloud computer, such as an Amazon EC2 instance, choose Ubuntu as the operation system.

### 3) Use a docker on your computer
If you want to use your own computer but prefer a clean environment for bioinformatics, you may use Docker to mount a Ubuntu on it. Briefly, there are three steps.
1. Install docker. Get the installation file from <https://docs.docker.com/get-started/get-docker/> and install.
2. In the terminal, run
```
$ docker pull ubuntu
```
This will pull the latest ubuntu image.

3. Mount the image
```
$ docker run -it ubuntu /bin/bash
```
This will start the image.

There are many online articles talking about how to use docker, suc as <https://www.makeuseof.com/run-ubuntu-as-docker-container/>

### 4) Use the current system as is.
I check all programs I used for this guide, and none of them really require a dedicated operation system. Although I did not test, I think it way work to use whatever major system your are currently working with.

## Step 2. Install anaconda/miniconda
Miniconda may be a better option for a cean environment, although I myself am using anaconda.

- Pick up either one to install after registration (or skip) at <https://www.anaconda.com/download/>.
- By default, a few tools, such as python, have been installed alone with conda.
- Once the conda is installed, there should be a (base) in the terminal.
- Create a conda environment
```
(base) .. $ conda create -n myenv
```
Here, myenv can be replaced with a more preferred name.
After that, activate myenv.
```
(base) .. $ conda activate myenv
(myenv) .. $
```

## Step 3. Install java runtime environment (JRE)
On ubuntu, use the command below to instlall JRE if not already installed.
```sh
$ sudo apt install openjdk-21-jre
```
On other systems, an appropriate JDK may be downloaded from <https://www.oracle.com/java/technologies/downloads/> and installed.

## Step 4. Install needed programs
1. Fastqc and multiqc for result quality check.
```
$ sudo apt install fastqc multiqc
```
2. Trimmomatic for sequence trimming.
```
$ sudo apt install trimmomatic
```
3. Set up bioconda following <https://bioconda.github.io/>. Below are the commands I borrowed from their website.
```
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
$ conda config --set channel_priority strict
```
This will add bioconda repositories as the resources to install other programs.
4. Install Hisat2 following <https://anaconda.org/bioconda/hisat2>. Aso,here is the command:
```
$ conda install bioconda::hisat2
```
5. Install R and R Studio.
- Install R by using the command below:
```
$ sudo apt install r-base
```
- Download R Studio from <https://posit.co/download/rstudio-desktop/> and install. Ussally, a double click on the package should automatically initiate the installation.

6. Install bioconductor in R Studio following <https://www.bioconductor.org/install/>
Currently, the code is below:
```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
```
7. Install needed R or bioconductor packages
```r
install.packages("dplyr")
install.packages("stringr")
install.packages("jsonlite")
install.packages("pheatmap")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("enrichplot")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("ReactomePA")
```
If the dependent package $fgsea$ cannot be properly called likely due to an outdated version installed, an error will show up. In this case, install the right version of this package mammually first using the codes below.
```r
install.packages("remotes") #Dependent pacakge "fgsea" cannot be installed by other methods. Will have to install from Github.
remotes::install_github("ctlab/fgsea")
```
8. Install ssh service for remote access (optional)
HPC or cloud computers may have already provided this service. The ssh service will allow the user to access the computer remotely. For example, a Windows user can use Putty to access the Linux computer through ssh and upload download files through sftp.

```
$ sudo apt install openssh-server
```
After rebooting, the ssh service should be automatically running. the GUI can be turned off and on through the following two commands. A dedicated command line terminal environment may be more favorable for computation purposes.
```
$ sudo systemctl isolate multi-user.target
$ sudo systemctl isolate graphical.target
```
To set the command line as the default interface, use the command below:
```
$ sudo systemctl set-default multi-user.target
```


After these steps, a working environment should be pretty much all set.