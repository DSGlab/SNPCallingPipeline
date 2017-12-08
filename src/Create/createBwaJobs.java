package Create;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.File;

public class createBwaJobs{
    public static void main(String [] args) throws FileNotFoundException{

		String INDIR = "/scratch/d/dguttman/cizydor/fastq_files/";//need to make sure this has a '/' at the end

		String ALREF = "/home/d/dguttman/jdiaz/reference/LESB58.fa";

		String REF = "/home/d/dguttman/jdiaz/reference/LESB58.fa";

		String OUTDIR = "/scratch/d/dguttman/jdiaz/CONRAD_CF/";//need to make sure this has a '/' at the end

		String JOBS_DIR = "";

		String isoList[] = {"10","130","14","152","160","17","180","196","19","200","217","21","225","235","23","251","257","258","259","263","265","270","273","275","276","278","280","282","288","290","293","298C","303","306","30","310","315","317","318","323","325","326","327","330","331","332","335","336","337","340","341","342","358","359","360","366","367","369","36","375","379","37","380","385","388","390","392","393","395","400","401","402","404C","405C","406C","409","410C","411C","412","417C","418C","419C","420N2","427","428","429","430","431","432","433","446","448","449","450","453","455","457","459C","465","471","472","473","475","476","479","480","487","501","503","504","505","506","507","508","509","510","511","512","513","514","521","525","526","527","528","529","530","531","532","538","539","53","540","541-2","542","544","548","549","550","551","552","553C","557","559-2","562","563","566","567","568","569","570","571","572-1","573-2","573-3","575","577","578","580","581","583","6","7","8","97"};

		boolean quaked = false;

		for(String id:isoList){
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

        String ISOLATE_DIR = "BWA-"+COR_ID+ISOLATE_ID+"/";
		PrintWriter out = new PrintWriter(JOBS_DIR+"/"+"bwa-"+COR_ID+ISOLATE_ID+".sh");

		/*Create Folder where the data will be saved*/
		new File(OUTPUT_DIR+ISOLATE_DIR).mkdirs();

		/* Header required for scinet */
		out.println("#!/bin/bash");
		out.println("#MOAB/Torque submission script for Multiple Serial Jobs");
		out.println("#PBS -l nodes=2:ppn=8,walltime=18:00:00");
		out.println("#PBS -N serialx8");
		out.println("module load gnu-parallel\n");
		out.println("module load bwa");
		out.println("module load samtools");
		out.println("module load bcftools");

			/* The actual alignment */
		out.println("bwa aln "+ALIGNER_REF+" "+INPUT_DIR+ISOLATE_ID+"_1"+COR_EXTENSION+".fq > "+OUTPUT_DIR+ISOLATE_DIR+"rOne.sai");
		out.println("bwa aln "+ALIGNER_REF+" "+INPUT_DIR+ISOLATE_ID+"_2"+COR_EXTENSION+".fq > "+OUTPUT_DIR+ISOLATE_DIR+"rTwo.sai");
		out.println("bwa sampe "+ALIGNER_REF+" "+OUTPUT_DIR+ISOLATE_DIR+"rOne.sai "+OUTPUT_DIR+ISOLATE_DIR+"rTwo.sai "+INPUT_DIR+ISOLATE_ID+"_1"+COR_EXTENSION+".fq "+INPUT_DIR+ISOLATE_ID+"_2"+COR_EXTENSION+".fq > "+OUTPUT_DIR+ISOLATE_DIR+"r.sam" );
		out.println("rm "+OUTPUT_DIR+ISOLATE_DIR+"rOne.sai");
		out.println("rm "+OUTPUT_DIR+ISOLATE_DIR+"rTwo.sai\n");

		out.println("## SAM FILE IS CONVERTED INTO BINARY FORM ##");
		out.println("samtools view -bS "+OUTPUT_DIR+ISOLATE_DIR+"r.sam > "+OUTPUT_DIR+ISOLATE_DIR+"r.bam");
		out.println("## SORT BAM FILE FOR ##");
		out.println("samtools sort "+OUTPUT_DIR+ISOLATE_DIR+"r.bam "+OUTPUT_DIR+ISOLATE_DIR+"r_sorted");
		out.println("## REMOVE r.sam FILE ##");
		//out.println("rm "+MAINDIR+ISOLATE+"/r.sam");
		//out.println("rm "+MAINDIR+ISOLATE+"/r.bam\n");

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