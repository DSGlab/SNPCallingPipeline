# SNPCallingPipeline
SNP filtering pipeline by julio.diaz@mail.utoronto.ca @ Guttman lab.<br>
Download jar file: https://github.com/DSGlab/SNPCallingPipeline/raw/master/SNPCallingPipeline.jar<br>
Sample configuration file: https://github.com/DSGlab/SNPCallingPipeline/raw/master/conf.txt<br>

To get details on how this pipeline works, run:<br>
<code>java -jar SNPCallingPipeline.jar</code>

# Requirements:
- Paired-End sequencing reads in the format: \<id>_1.fq \<id>_2.fq (id is identical as described in the id list)<br>
- Reference headers in the format: \<refName>\_\<repliconType>\_\<repliconNum> (e.g. pao1_chromosome_1)
