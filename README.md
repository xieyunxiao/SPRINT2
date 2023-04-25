
# SPRINT2

SPRINT2 is an enhanced version for SPRINT to detect RNA editing sites from RNA-seq data.



## 1. Old version.

**Binary sources and details at** https://github.com/jumphone/SPRINT/

**Publication** https://doi.org/10.1093/bioinformatics/btx473


## 2. Installation.

**You can install the package via pip in Linux:**

```
pip install sprint2
```
**or install via Git clone:** 
```
git clone https://github.com/xieyunxiao/SPRINT2/
```

## 3. Download annotation files.
**Repeat annotation files can download at:** https://github.com/xieyunxiao/SPRINT2/repeat_data/

**Double-stranded RNA annotation files can download at:** https://github.com/xieyunxiao/SPRINT2/dsRNA_data/


## 4. Requirements.
```
Unix & Python == 3.8
```
**If you had installed Anaconda, please run:**

``
bash require.bash
``

**If you didn't install Anaconda, please install these tools:**
```
SAMTOOLS == 1.2 (using htslib == 1.2.1) 
BEDTOOLS == v2.30.0
BWA ==  0.7.12 
BLAT == v.36x7
```

**Please install some python packages via pip in Linux:** 

``
pip install -r requirements
``



## 5. Usage.

**Step1.Mask A with G on reference FASTA file.**
```
usage: MaskAGref [-h] [-r REFERENCE] [-o OUTPUT] [-t GTF] [-b BWA] [-p CPU]

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        path to reference FASTA file (default:~/reference.fa)
  -o OUTPUT, --output OUTPUT
                        path to output directory (default:~/)
  -t GTF, --gtf GTF     path to transcript annotation GTF file #Optional
  -b BWA, --bwa BWA     path to BWA (default:~/bwa)
  -p CPU, --cpu CPU     CPU number (default:CPU=1)
```
**Step2.Get candidate double-stranded RNAs.**
```
usage: getDsRNA [-h] [-o OUT_DIR] [-r REFERENCE_FILE] [-t TRANSCRIPT_FILE] [-b BLAT_FILE] [-B BEDTOOLS_PATH] [-lc LC_PATH] [-p CPU]

optional arguments:
  -h, --help            show this help message and exit
  -o OUT_DIR, --OUT_DIR OUT_DIR
                        path to output directory (default:~/)
  -r REFERENCE_FILE, --Reference_file REFERENCE_FILE
                        path to reference FASTA file (default:~/reference.fa)
  -t TRANSCRIPT_FILE, --Transcript_file TRANSCRIPT_FILE
                        path to transcript annotation bed file (default:~/transcript.bed)
  -b BLAT_FILE, --Blat_file BLAT_FILE
                        path to BLAT (default:~/blat)
  -B BEDTOOLS_PATH, --Bedtools_path BEDTOOLS_PATH
                        path to BEDTOOLS (default:~/bedtools)
  -lc LC_PATH, --lc_path LC_PATH
                        path to Low_complexity.txt (default:~/Low_complexity.txt)
  -p CPU, --cpu CPU     CPU number (default：CPU=1)
```
**Step3.Get RNA editing sites.**
*(Note: Do not use .GZ files. Please "gunzip" the ".fastq.gz" files.)*
```
usage: getRES [-h] [-s SNV_PATH] [-o RES_PATH] [-r1 READ1] [-r2 READ2] [-r REFERENCE_FILE] [-rp REPEAT_FILE] [-b BWA] [-B BEDTOOLS] [-S SAMTOOLS] [-ds CANDIDATE_DSRNA] [-p CPU]

Get candidate double-stranded RNA

optional arguments:
  -h, --help            show this help message and exit
  -s SNV_PATH, --SNV_PATH SNV_PATH
                        path to SNV directory (default:~/1_SNV_calling/)
  -o RES_PATH, --RES_PATH RES_PATH
                        path to RES directory (default:~/2_RES_calling/)
  -r1 READ1, --read1 READ1
                        path to FASTQ_1 file (default:~/test_1.fastq.fa)
  -r2 READ2, --read2 READ2
                        path to FASTQ_2 file #Optional
  -r REFERENCE_FILE, --Reference_file REFERENCE_FILE
                        path to reference FASTA file (default:~/reference.fa)
  -rp REPEAT_FILE, --Repeat_file REPEAT_FILE
                        path to repeat annotation BED file #Optional
  -b BWA, --bwa BWA     path to BWA (default:~/bwa)
  -B BEDTOOLS, --bedtools BEDTOOLS
                        path to BEDTOOLS (default:~/bedtools)
  -S SAMTOOLS, --samtools SAMTOOLS
                        path to SAMTOOLS (default:~/samtools)
  -ds CANDIDATE_DSRNA, --Candidate_dsRNA CANDIDATE_DSRNA
                        path to candidate dsRNA directory #Optional
  -p CPU, --cpu CPU     CPU number (default：CPU=1)
```

**Here is an example of how to use sprint2:**
```
cd test/
MaskAGref -o ./ -r reference.fa -t annotation.gtf -b ./bwa -p 10
getDsRNA -o ./ -r reference.fa -t transcript.bed -b ./blat -B /local/bin/bedtools -lc Low_complexity.txt -p 10
getRES -s ./1_SNV_calling/ -o ./2_RES_calling/ -r1 test_1.fastq -r2 test_2.fastq -r reference.fa -rp repeat.bed -b ./bwa -B /local/bin/bedtools -S ./samtools -ds ./dsRNA_file/ -p 10
```

## 6. Get repeat annotation file.
**Step1.Install RepeatMasker:**
Download the latest version of RepeatMasker installation package from the RepeatMasker website (http://www.repeatmasker.org/RMDownload.html) and install it on a Linux system using the following commands:
```
wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.5.tar.gz
tar zxvf RepeatMasker-4.1.5.tar.gz
cp RepeatMasker-4.1.5.tar.gz /usr/local/
cd /usr/local/
tar zxvf RepeatMasker-4.1.5.tar.gz
cd RepeatMasker
perl ./configure
```

During the installation, other software and tools required by RepeatMasker need to be downloaded and installed.

**Step2.Download genome sequence:**
Download the fasta format file of species genome sequence from UCSC Genome Browser, for example:
```
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
```
**Step3.Download the sequence databases required by RepeatMasker:**
Please download the required sequence databases from the RepeatMasker website, for example:
```
wget https://www.dfam.org/releases/Dfam_3.7/families/Dfam.h5.gz
gunzip Dfam.h5.gz
mv Dfam.h5 /usr/local/RepeatMasker/Libraries
tar zxvf RepeatMaskerGenomeAnnotations_LATEST.tar.gz
```
**Step4.Please use the following command to run RepeatMasker:**
```
RepeatMasker -species human -dir output_dir hg19.fa
```
The -species parameter specifies the sequence database used by RepeatMasker, the -dir parameter specifies the output directory, and hg19.fa is the input sequence file. After the analysis is finished, RepeatMasker will generate multiple result files in the output directory, including the masked sequence file and the statistics file.

**Step5.Get repeat annotation file and low complex annotation file**
```
mv hg19.out repeat.bed
getLowComplexity.py repeat.bed Low_complexity.txt
```
Note that the sequence databases downloaded in the above steps are protected by copyright and require a license to use. Also, the above code is for reference only and the actual operation steps may vary depending on the situation.

## 7.Input data files.
 - reference.fa: Reference genome FASTA file.
 - annotation.gtf: Gene annotation GTF file.
 - transcript.bed: Transcripts BED file (Optional file).
 
```
______________________________________________________________________________________________
| Chrom | Start | End | transcript | Strand | Gene Symbol | Transcript ID |  Transcript type |
______________________________________________________________________________________________
```
 - read1.fastq: Required!!!
 - read2.fastq: If single-read sequencing, ignore this parameter.
 - repeat.bed: Repeat annotation file. If this parameter is missing, without repeat-based regular RES.
 - Low_complexity.txt: Low complexity region file.
 - dsRNA_file/: Candidate double-stranded RNA. If this parameter is missing, without dsRNA-based regular RES.
 
```
____________________________________________________
| Chrom1 | Start1 | End1 |  Chrom2 | Start2 | End2 | 
____________________________________________________
```
 - bwa: Path to BWA.
 - blat: Path to BLAT.
 - bedtools: Path to BEDTOOLS.
 - samtools: Path to SAMTOOLS.

## 8.Output files.
 - reference_mskAG.fa: Mask A with G on reference FASTA file.
 - dsRNA_file/: Candidate double-stranded RNA.
 - 1_SNV_calling/: SNV files.
  regular.snv  | hyper_mskAG.snv  | hyper_mskTC.snv 
```
__________________________________________________________________
| Chrom | Start(0-base) | End(1-base) |  Type  | Supporting_reads | 
__________________________________________________________________
```
 - 2_RES_calling/: RES files.
```
__________________________________________________________________________
| Chrom | Start(0-base) | End(1-base) |  Type | Supporting_reads | Depth |
__________________________________________________________________________
```
## 9.License.

This project is licensed under the MIT License - see the LICENSE file for details.








