package Create;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.File;

public class createNovoalignJobs{
    public static void main(String [] args) throws FileNotFoundException{

		String INDIR = "/scratch/d/dguttman/cizydor/fastq_files/";//need to make sure this has a '/' at the end

		String ALREF = "/home/d/dguttman/jdiaz/reference/LESB58_NOVOALIGN";

		String REF = "/home/d/dguttman/jdiaz/reference/LESB58.fa";

        String OUTDIR = "/scratch/d/dguttman/jdiaz/CONRAD_LAST/";//need to make sure this has a '/' at the end

		String JOBS_DIR = "";

		String isoList[] = {"10","130","14","152","160","17","180","196","19","200","217","21","225","235","23","251","257","258","259","263","265","270","273","275","276","278","280","282","288","290","293","298C","303","306","30","310","315","317","318","323","325","326","327","330","331","332","335","336","337","340","341","342","358","359","360","366","367","369","36","375","379","37","380","385","388","390","392","393","395","400","401","402","404C","405C","406C","409","410C","411C","412","417C","418C","419C","420N2","427","428","429","430","431","432","433","446","448","449","450","453","455","457","459C","465","471","472","473","475","476","479","480","487","501","503","504","505","506","507","508","509","510","511","512","513","514","521","525","526","527","528","529","530","531","532","538","539","53","540","541-2","542","544","548","549","550","551","552","553C","557","559-2","562","563","566","567","568","569","570","571","572-1","573-2","573-3","575","577","578","580","581","583","6","7","8","97"};

		boolean quaked = false;

		for(String id:isoList){
		    run(INDIR, ALREF, REF, OUTDIR, id, quaked, JOBS_DIR,"","","");
		}
    }

    public static void run(String INPUT_DIR, String ALIGNER_REF, String REFERENCE, String OUTPUT_DIR, String ISOLATE_ID, boolean isQuakeCorrected, String JOBS_DIR, String NOVOALIGN_PATH, String SAMTOOLS_PATH, String BCFTOOLS_PATH) throws FileNotFoundException{
  
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

        String ISOLATE_DIR = "NOVOALIGN-"+COR_ID+ISOLATE_ID+"/";
		PrintWriter out = new PrintWriter(JOBS_DIR+"/"+"novoalign-"+COR_ID+ISOLATE_ID+".sh");

		/*Create Folder where the data will be saved*/
        new File(OUTPUT_DIR+ISOLATE_DIR).mkdirs();
	
		/* Header required by scinet */
		out.println("#!/bin/bash");

		/* The actual alignment */
		out.println(NOVOALIGN_PATH+"/novoalign -o SAM -f "+INPUT_DIR+ISOLATE_ID+"_1"+COR_EXTENSION+".fq "+INPUT_DIR+ISOLATE_ID+"_2"+COR_EXTENSION+".fq -d "+ALIGNER_REF+" > "+OUTPUT_DIR+ISOLATE_DIR+"r.sam\n");

		out.println(SAMTOOLS_PATH+"/samtools fixmate -O BAM "+OUTPUT_DIR+ISOLATE_DIR+"r.sam "+OUTPUT_DIR+ISOLATE_DIR+"r.bam");
		out.println(SAMTOOLS_PATH+"/samtools sort -O BAM -o "+OUTPUT_DIR+ISOLATE_DIR+"r.sorted.bam "+OUTPUT_DIR+ISOLATE_DIR+"r.bam\n");

		out.println(BCFTOOLS_PATH+"/bcftools mpileup -Ob -o "+OUTPUT_DIR+ISOLATE_DIR+"r.bcf -f "+REFERENCE+" "+OUTPUT_DIR+ISOLATE_DIR+"r.sorted.bam");
		out.println(BCFTOOLS_PATH+"/bcftools call --ploidy 1 -vmO v -o "+OUTPUT_DIR+ISOLATE_DIR+"r.vcf "+OUTPUT_DIR+ISOLATE_DIR+"r.bcf\n");

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

	private static String getParentDirectory(String fullPath){
		return fullPath.split("/")[fullPath.split("/").length-1];
	}
}
