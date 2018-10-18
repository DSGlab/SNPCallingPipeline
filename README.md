# SNPCallingPipeline
This is a SNP calling tool, which implements read mapping from three different algorithms ([BWA](http://bio-bwa.sourceforge.net), [LAST](http://last.cbrc.jp/doc/last.html) and [Novoalign](http://www.novocraft.com/products/novoalign/)) and evaluates their output. This pipeline currently is only suitable for bacterial genomes.

SNP filtering pipeline by julio.diaz@mail.utoronto.ca @ [Guttman lab](https://guttman.csb.utoronto.ca).<br>


# Install
* Download the [SNPCallingPipeline](https://github.com/DSGlab/SNPCallingPipeline/raw/master/SNPCallingPipeline.jar)<br>
`$ wget https://github.com/DSGlab/SNPCallingPipeline/raw/master/SNPCallingPipeline.jar`
* Download the sample [configuration file](https://github.com/DSGlab/SNPCallingPipeline/raw/master/conf.txt)<br>
`$ wget https://github.com/DSGlab/SNPCallingPipeline/raw/master/conf.txt`
* Test
```sh
$ java -jar SNPCalling Pipeline.jar
SNP Calling Pipeline v. 1.12
Questions? julio.diaz@mail.utoronto.ca

Usage:	java -jar SNPCallingPipeline ANALYSIS_NAME [config file]

Step 0:		scinetjobcreator	Create jobs to be used in scinet's niagara
Step 1:		gethqsnps		Get list of High quuality (confidence) SNPs
Optional step:	getintraclonalsnps	Remove SNPs that are fixed between the isolates and the reference
Step 2:		snpchecker		Gets raw calls at the SNP positions in all the isolates
Step 3:		snpfilter		Filters raw SNP calls
Step 4:		createalignment		Creates alignment based on filtered SNP calls
```

# Setup:
### Required Files
* Paired-End sequencing reads in the format: `\<id>_1.fq \<id>_2.fq` (id is identical as described in the id list).<br>
* Reference headers in the format: `\<refName>\_\<repliconType>\_\<repliconNum>` (e.g. pao1_chromosome_1).
* File including list of strain names (One strain name per line).

### Index Reference files
##### Index BWA
Example: `bwa index reference.fa`
##### Index LAST
Example: `lastdb reference_LAST reference.fa`
##### Index NOVOALIGN
Example: `novoindex reference_NOVOALIGN reference.fa`

### Strain List File
One strain per line.
```
strainA
strainB
strainC
strainD
```

### Configuration File:
##### SETTTINGS FOR scinetojobcreator STEP
* `REF_FILE` - The reference file in fasta format. Must include complete path to file.
* `ID_LIST_FILE` - The file including the IDs of all the strains to be taken into acccount (one per line). Must include complete path to file.
* `ALIGNER` - The name of the aligner to be taken into account for the second step (options: "BWA","LAST","NOVOALIGN").
* `INPUT_DIR` - The directory which includes the sequencing reads.
* `ALIGNMENT_DIR` - The directory where the alignment files will be saved.
* `JOBS_DIR` - The directory where the bash jobs will be saved.
* `INCLUDE_QUAKE` - TRUE if pipeline is set to include non-quality controlled reads (FALSE otherwise).
* `INCLUDE_BWA` - TRUE if 'BWA' alignment is required (FALSE not supported at the moment).
* `INCLUDE_LAST` - TRUE if 'LAST' alignment is required (FALSE not supported at the moment).
* `INCLUDE_NOVOALIGN` - TRUE if 'NOVOALIGN' alignment is required (FALSE not supported at the moment).
* `INCLUDE_BOWTIE2` - TRUE if 'BWA' alignment is required (FALSE is the only option supported at the moment).
* `BWA_REF` - The 
* `NOVOALIGN_REF` - 
* `LAST_REF` - 
* `BOWTIE_REF` - 
* `BWA_PATH` - Path to the 'BWA' executable.
* `LAST_PATH` - Path to the 'LAST' executable.
* `NOVOALIGN_PATH` - Path to the 'NOVOALIGN' executable.
* `SAMTOOLS_PATH` - Path to the 'SAMTOOLS' executable.
* `BCFTOOLS_PATH` - Path to the 'BCFTOOLS' executable.
##### SETTINGS FOR gethqsnps STEP
* `HQ_SNP_LIST_OUTPUT` - The output file of this step. Must include complete path to file.
* `HQ_SNP_WORKING_DIR` - The directory where the alignment results are stored.
* `HQ_SNP_QUALITY` - The minimum quality score required to be considered a HQ SNP.
* `HQ_SNP_DEPTH` - The minimum depth required to be considered a HQ SNP.
* `HQ_SNP_DIST_TO_CONTIG_END` - The minimum distance to the end of the replicon required to be considered a HQ SNP.
* `HQ_SNP_READ_BALANCE` - The minimum number of reads required from both reverse and forward directions to be considered a HQ SNP.
* `HQ_SNP_REF_TO_ALT_BALANCE` - The minimum reference to alternative call ratio required to be considered a HQ SNP.
* `HQ_SNP_CLUSTER_SIZE` - The required distance between other SNPs to be considered a HQ SNP.
* `HQ_SNP_REQUIRED_NUM` - The minimum number of aligners that call the SNP to be considered a HQ SNP (Currently only '3' is supported).
* `INCLUDE_PRE_QC` - TRUE if pipeline is set to include non-quality controlled reads (FALSE otherwise).
##### SETTINGS FOR getintraclonalsnps STEP
* `INTRA_SNP_INPUT` - The output file from the gethqsnps STEP. Must include complete path to file.
* `INTRA_SNP_WORKING_DIR` - The directory where the alignment results are stored
* `INTRA_SNP_OUTPUT` - The output file of this step. 
##### SETTINGS FOR snpchecker STEP
* `SNP_CHECKER_INPUT` - The output file from the gethqsnps (or getintraclonalsnps) STEP. Must include complete path to file.
* `SNP_CHECKER_WORKING_DIR` - The directory where the alignment results are stored
* `SNP_CHECKER_OUTPUT` - The output file of this step. 
##### SETTINGS FOR snpfilter STEP
* `SNP_FILTER_INPUT` - The output file from the snpchecker STEP. Must include complete path to file.
* `SNP_FILTER_OUTPUT` - The output file of this step. 
* `SNP_FILTER_MIN_GOOD_BALANCE` - The minimum number of reads required from both reverse and forward directions to be considered a SNP.
* `SNP_FILTER_MIN_GOOD_CALL_RATIO` - The minimum ratio required between a reference and alternative call to be considered a SNP.
* `SNP_FILTER_MIN_GOOD_QUALITY` - The minimum quality score required to be considered a SNP.
##### SETTINGS FOR createalignment STEP
* `CREATE_ALIGN_INPUT` - The output file from the snpfilter STEP. Must include complete path to file.
* `CREATE_ALIGN_OUTPUT` - The output alignment file in fasta.
* `CREATE_ALIGN_LIST_OUTPUT` - The output file including the positions where SNPs were detected.

# Easy Guide
### Setup
1. Create [strain list file](###Strain List File).
2. Have a reference sequences and paired-end sequencing reads in the [required format](###Required Files).
3. [Index](###Index Reference files) the reference.
4. Modify [configuration file](###Configuration File) to match your criteria.
5. Check that setup was completed
```sh
$ java -jar SNPCallingPipeline.jar datachecker conf.txt

SNP Calling Pipeline v. 1.12
Questions? julio.diaz@mail.utoronto.ca

Usage:	java -jar SNPCallingPipeline ANALYSIS_NAME [config file]

Checking input files.

0 error(s) found in the reference or seq read files
0 error(s) found in the path to the software
```
6. Create jobs
```sh
$ java -jar SNPCallingPipeline.jar scinetjobcreator conf.txt
SNP Calling Pipeline v. 1.12
Questions? julio.diaz@mail.utoronto.ca

Usage:	java -jar SNPCallingPipeline ANALYSIS_NAME [config file]

Starting Create_Alignment analysis.
	Printing task for isolate:	strainA
	Printing task for isolate:	strainB
	Printing task for isolate:	strainC
	Printing task for isolate:	strainD
```
This created a number of bash scripts in the jobs directory
```sh
$ ls jobsDirectory
aligner0.sh     bwa-strainA.sh       bwa-strainB.sh       bwa-strainC.sh       bwa-strainD.sh
last-strainA.sh      last-strainB.sh      last-strainC.sh       last-strainD.sh
novoalign-strainA.sh novoalign-strainB.sh novoalign-strainC.sh novoalign-strainD.sh

```
The `bwa*.sh`, `novoalign*.sh`, and `last*.sh` files can be executed in the terminal as usual
```sh
$ bash bwa-strainA.sh
$ bash bwa-strainB.sh
$ bash bwa-strainC.sh
$ bash bwa-strainD.sh
$ bash last-strainA.sh
$ bash last-strainB.sh
$ bash last-strainC.sh
$ bash last-strainD.sh
$ bash novoalign-strainA.sh
$ bash novoalign-strainB.sh
$ bash novoalign-strainC.sh
$ bash novoalign-strainD.sh
``` 
OR the `aligner?.sh` bash scripts can be submitted to niagara in scinet.
```
$ sbatch aligner0.sh
```
Wait for the scripts to finsh.<br>
6. Run `gethqns` step
```sh
$ java -jar SNPCallingPipeline.jar gethqsnps conf.txt
```
7. Run `getintraclonalsnps` step
```sh
$ java -jar SNPCallingPipeline.jar getintraclonalsnps conf.txt
```
8. Run `snpchecker` step
```sh
$ java -jar SNPCallingPipeline.jar snpchecker conf.txt
```
9. Run `snpfilter` step
```sh
$ java -jar SNPCallingPipeline.jar snpfilter conf.txt
```
10. Run `createalignment` step
```sh
$ java -jar SNPCallingPipeline.jar createalignment conf.txt
```
11. Results
```sh
$ls outputDirectory
snp_alignment.fa
snp_list.fa
```
* `snp_alignment.fa` - all positions with HQ SNPs with the proper call
* `snp_list.fa` - info about every polymorphic position