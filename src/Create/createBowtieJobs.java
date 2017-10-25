package Create;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.File;

public class createBowtieJobs{
    public static void main(String [] args) throws FileNotFoundException{

		String INDIR = "/scratch/d/dguttman/jdiaz/unfilteredData/";

		String ALREF = "/scratch/d/dguttman/jdiaz/reference/C1Q_BOWTIE2";

		String REF = "/scratch/d/dguttman/jdiaz/reference/C1Q.fa";

		String OUTDIR = "/scratch/d/dguttman/jdiaz/CF67_C1Q/BOWTIE_QUAKE_C1Q/";
	
		String JOBS_DIR = "";

		String isoList[] = {"64B"};
	
		boolean quaked = true;

		for(String id : isoList){
	    	run(INDIR, ALREF, REF, OUTDIR, id, quaked, JOBS_DIR);
		}
    }

    public static void run(String INPUT_DIR, String ALIGNER_REF, String REFERENCE, String OUTPUT_DIR, String ISOLATE_ID, boolean isQuakeCorrected, String JOBS_DIR) throws FileNotFoundException{
	    
	/* Ensure directory String is appropriate */
	INPUT_DIR = checkDirString(INPUT_DIR);
	OUTPUT_DIR = checkDirString(OUTPUT_DIR);

	/* Set corrected or non corrected reads */
	String COR_EXTENSION = "";
	String COR_ID = "";
	if(isQuakeCorrected){
	    COR_EXTENSION = ".cor";
	    COR_ID = "COR_";
	}

        String ISOLATE_DIR = "BOWTIE2-"+COR_ID+ISOLATE_ID+"/";
        PrintWriter out = new PrintWriter(JOBS_DIR+"/"+"bowtie2-"+COR_ID+ISOLATE_ID+".sh");

        /*Create Folder where the data will be saved*/
        new File(OUTPUT_DIR+ISOLATE_DIR).mkdirs();
	
	/* Headers required by scinet to run the job */
	out.println("#!/bin/bash");
	out.println("#MOAB/Torque submission script for Multiple Serial Jobs");
	out.println("#PBS -l nodes=2:ppn=8,walltime=6:00:00");
	out.println("#PBS -N serialx8");
	out.println("module load gnu-parallel\n");
	out.println("module load bowtie2");
	out.println("module load samtools");
	out.println("module load bcftools");
	
	/* The actual alignment */                                                        
	out.println("bowtie2 --rg-id "+ISOLATE_ID+" --rg LB:library --rg PL:MiSeq --rg SM:"+ISOLATE_ID+" --rg CN:CAGEF -1 "+INPUT_DIR+ISOLATE_ID+"_1"+COR_EXTENSION+".fq -2 "+INPUT_DIR+ISOLATE_ID+"_2"+COR_EXTENSION+".fq -x "+ALIGNER_REF+" > "+OUTPUT_DIR+ISOLATE_DIR+"r.sam\n");
	
	out.println("## SAM FILE IS CONVERTED INTO BINARY FORM ##");                                                  
	out.println("samtools view -bS "+OUTPUT_DIR+ISOLATE_DIR+"r.sam > "+OUTPUT_DIR+ISOLATE_DIR+"r.bam");
	out.println("## SORT BAM FILE FOR ##");                                                       
	out.println("samtools sort "+OUTPUT_DIR+ISOLATE_DIR+"r.bam "+OUTPUT_DIR+ISOLATE_DIR+"r_sorted");
	out.println("## REMOVE r.sam FILE ##");
	//out.println("rm "+MAINDIR+ISOLATE_DIR+"r.sam");
	//out.println("rm "+MAINDIR+ISOLATE_DIR+"r.bam\n");
	
	out.println("## CREATE LIST OF POTENTIAL SNP OR INDEL ##");                                                      
	out.println("samtools mpileup -uf "+REFERENCE+" "+OUTPUT_DIR+ISOLATE_DIR+"r_sorted.bam > "+OUTPUT_DIR+ISOLATE_DIR+"r.bcf");
	out.println("## PARSE POTENTIAL SNP OR INDEL USING BAYESIAN INFERENCE ##");                         
	out.println("bcftools view -bvcg "+OUTPUT_DIR+ISOLATE_DIR+"r.bcf > "+OUTPUT_DIR+ISOLATE_DIR+"r2.bcf");
	out.println("## BCF FILE IS CONVERTED VIEWABLE FORM ##");                                          
	out.println("bcftools view "+OUTPUT_DIR+ISOLATE_DIR+"r2.bcf > "+OUTPUT_DIR+ISOLATE_DIR+"r.vcf");
	out.close();
    }

    private static String checkDirString(String inputDir){
	char lastChar = inputDir.charAt(inputDir.length()-1);
	if(lastChar=='/'){
	    return inputDir;
	}else{
	    return inputDir+"/";
	}
    }
}