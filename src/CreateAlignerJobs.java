import java.io.FileNotFoundException;
import java.io.File;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.Hashtable;
import java.util.ArrayList;

public class CreateAlignerJobs {
    public static void main(String [] args) throws FileNotFoundException, InterruptedException {

	    //String CONF_FILE = "/Users/juliofdiaz/Dropbox/CF/snp_calling/TEST/alignconf.txt";
	    //Hashtable<String,String> conf = getConfiguration(CONF_FILE);
	    Hashtable<String,String> conf = getConfiguration( args[0] );

	    String [] ISOLATE_LIST = getIsolates(conf);;

	    String INPUT_DIR = getInputDir(conf);
	    isStringDir(INPUT_DIR);
	    String REFERENCE = getReference(conf);
	    isStringFile(REFERENCE);
	    String OUTPUT_DIR = getOutputDir(conf);
	    isStringDir(OUTPUT_DIR);
        String JOBS_DIR = getJobsDir(conf);
        isStringDir(JOBS_DIR);

        String BWA_PATH = getBwaPath(conf);
        isStringDir(BWA_PATH);
        String NOVOALIGN_PATH = getNovoAlignPath(conf);
		isStringDir(NOVOALIGN_PATH);
        String LAST_PATH = getLastPath(conf);
		isStringDir(LAST_PATH);
        String SAMTOOLS_PATH = getSamtoolsPath(conf);
		isStringDir(SAMTOOLS_PATH);
        String BCFTOOLS_PATH = getBcftoolsPath(conf);
		isStringDir(BCFTOOLS_PATH);

	    boolean includeQuake = isIncludeQuake(conf);
    	boolean includeBowtie = isIncludeBowtie(conf);
    	boolean includeBwa = isIncludeBwa(conf);
    	boolean includeLast = isIncludeLast(conf);
    	boolean includeNovoalign = isIncludeNovoalign(conf);

    	String BOWTIE_REF = getBowtieRef(conf);
    	if(includeBowtie){
    	    isStringFile(BOWTIE_REF);
    	}
    	String BWA_REF = REFERENCE;
    	if(includeBwa){
    	    isStringFile(BWA_REF);
    	}
    	String LAST_REF = getLastRef(conf);
	    //if(includeLast){
    	//    isStringFile(LAST_REF);
    	//}
    	String NOVO_REF = getNovoalignRef(conf);
	    if(includeNovoalign){
    	    isStringFile(NOVO_REF);
    	}

    	ArrayList<String> tasks = new ArrayList<String>();
    	for(String ISOLATE : ISOLATE_LIST){
	        System.out.println("\tPrinting task for isolate:\t"+ISOLATE);
            if(includeBowtie){
	    		Create.createBowtieJobs.run(INPUT_DIR, BOWTIE_REF, REFERENCE, OUTPUT_DIR, ISOLATE, false, JOBS_DIR);
	    		tasks.add("bash bowtie2-"+ISOLATE+".sh 2>bowtie2-"+ISOLATE+".sh.err >bowtie2-"+ISOLATE+".sh.out &");
	        }
	        if(includeBwa){
                Create.createBwaJobs.run(INPUT_DIR, BWA_REF, REFERENCE, OUTPUT_DIR, ISOLATE, false, JOBS_DIR, BWA_PATH, SAMTOOLS_PATH, BCFTOOLS_PATH);
				tasks.add("bash bwa-"+ISOLATE+".sh 2>bwa-"+ISOLATE+".sh.err >bwa-"+ISOLATE+".sh.out &");
	        }
            if(includeLast){
                Create.createLastJobs.run(INPUT_DIR, LAST_REF, REFERENCE, OUTPUT_DIR, ISOLATE, false, JOBS_DIR, LAST_PATH, SAMTOOLS_PATH, BCFTOOLS_PATH);
				tasks.add("bash last-"+ISOLATE+".sh 2>last-"+ISOLATE+".sh.err >last-"+ISOLATE+".sh.out &");
			}
	        if(includeNovoalign){
                Create.createNovoalignJobs.run(INPUT_DIR, NOVO_REF, REFERENCE, OUTPUT_DIR, ISOLATE, false, JOBS_DIR, NOVOALIGN_PATH, SAMTOOLS_PATH, BCFTOOLS_PATH);
				tasks.add("bash novoalign-"+ISOLATE+".sh 2>novoalign-"+ISOLATE+".sh.err >novoalign-"+ISOLATE+".sh.out &");
			}

	        if(includeQuake){
		        if(includeBowtie){
                    Create.createBowtieJobs.run(INPUT_DIR, BOWTIE_REF, REFERENCE, OUTPUT_DIR, ISOLATE, true, JOBS_DIR);
					tasks.add("bash bowtie2-COR_"+ISOLATE+".sh 2>bowtie2-COR_"+ISOLATE+".sh.err >bowtie2-COR_"+ISOLATE+".sh.out &");
				}
		        if(includeBwa){
                    Create.createBwaJobs.run(INPUT_DIR, BWA_REF, REFERENCE, OUTPUT_DIR, ISOLATE, true, JOBS_DIR, BWA_PATH, SAMTOOLS_PATH, BCFTOOLS_PATH);
					tasks.add("bash bwa-COR_"+ISOLATE+".sh 2>bwa-COR_"+ISOLATE+".sh.err >bwa-COR_"+ISOLATE+".sh.out &");
				}
		        if(includeLast){
                    Create.createLastJobs.run(INPUT_DIR, LAST_REF, REFERENCE, OUTPUT_DIR, ISOLATE, true, JOBS_DIR, LAST_PATH, SAMTOOLS_PATH, BCFTOOLS_PATH);
					tasks.add("bash last-COR_"+ISOLATE+".sh 2>last-COR_"+ISOLATE+".sh.err >last-COR_"+ISOLATE+".sh.out &");
		        }
		        if(includeNovoalign){
                    Create.createNovoalignJobs.run(INPUT_DIR, NOVO_REF, REFERENCE, OUTPUT_DIR, ISOLATE, true, JOBS_DIR, NOVOALIGN_PATH, SAMTOOLS_PATH, BCFTOOLS_PATH);
					tasks.add("bash novoalign-COR_"+ISOLATE+".sh 2>novoalign-COR_"+ISOLATE+".sh.err >novoalign-COR_"+ISOLATE+".sh.out &");
		        }
	        }
	    }

		createSubmitters(tasks, JOBS_DIR);



	    System.out.println();
    }



    private static void createSubmitters(ArrayList<String> tasks, String dir) throws FileNotFoundException {
    	int count = 0;
		PrintWriter out = null;// = new PrintWriter(dir+"/submitter"+1+".sh");
		//System.out.println(dir+"/submitter"+1+".sh");
		for(String task:tasks){
			count++;
			if(count%40==1){
				out = new PrintWriter(dir+"/aligner"+count/40+".sh");

				out.println("#!/bin/bash");
				out.println("#SBATCH --nodes=1");
				out.println("#SBATCH --cpus-per-task=40");
				out.println("#SBATCH --time=24:00:00");
				out.println("#SBATCH --job-name=aligner_jobs-"+count/40);
				out.println("#SBATCH --output=aligner_jobs-"+count/40+".txt\n");

				out.println(task);
			}else if(count==tasks.size()){
				out.println(task);
				out.println("wait");
				out.close();
				//System.out.println("close printwriter");
			}else if(count%40==0){
				out.println(task);
				out.println("wait");
				out.close();
				//System.out.println("close printwriter");
			}else{
				out.println(task);
			}
		}
	}

    /*
     * WILL NOT DETECT IT IF NAME IS TOO SHORT 
     */
    private static void isStringFile(String inString){
	    boolean isFileFound = false;
        //System.out.println(inString);
	    File testFile = new File(inString);
	    String testFileName = testFile.getName();
	    File parentFile = testFile.getParentFile();
	    File [] listOfFiles = parentFile.listFiles();

    	for(File curFile:listOfFiles){
	        String curFileName = curFile.getName();
	        if( curFileName.length() >= testFileName.length() ){
		        if(testFileName.equals(curFileName.substring(0,testFileName.length()))){
		            //System.out.print("same\t");
		            isFileFound = true;
		        }
	        }
	    }

	    if(!isFileFound){
            System.out.println("ERROR: "+inString+" IS NOT A FILE");
            System.exit(1);
	    }
    }

    private static void isStringDir(String inString){
	    File testFile = new File(inString);
	    if(!testFile.isDirectory()){
	        System.out.println("ERROR: "+inString+" IS NOT A DIRECTORY");
	        System.exit(1);
	    }
    }

    private static Hashtable<String,String> getConfiguration(String confFile)throws FileNotFoundException{
    
	    Scanner input = new Scanner(new File(confFile));
	
	    Hashtable<String, String> confTable = new Hashtable<String, String>();

	    while(input.hasNextLine()){
	        String curLine = input.nextLine();
	        if(!curLine.equals("")){
		        char firstChar =  curLine.charAt(0);
		        if(firstChar != '#'){
		            String[] curItem = curLine.split("=");
		            if(curItem.length == 2){
			            confTable.put(curItem[0].trim(), curItem[1].trim());
		            }else{
			            System.out.println("CONFIGURATION FILE ERROR");
			            System.exit(1);
		            }
		        }
	        }
	    }

	    return confTable;
    }

    private static String getBwaPath(Hashtable<String, String> confTable){
		return getParameter(confTable, "BWA_PATH");
	}

	private static String getLastPath(Hashtable<String, String> confTable){
		return getParameter(confTable,"LAST_PATH");
	}

	private static String getNovoAlignPath(Hashtable<String, String> confTable){
		return getParameter(confTable,"NOVOALIGN_PATH");
	}

	private static String getSamtoolsPath(Hashtable<String, String> confTable){
		return getParameter(confTable,"SAMTOOLS_PATH");
	}

	private static String getBcftoolsPath(Hashtable<String, String> confTable){
		return getParameter(confTable,"BCFTOOLS_PATH");
	}

    private static String getInputDir(Hashtable<String, String> confTable){
	    return getParameter(confTable,"INPUT_DIR");
    }

    private static String getOutputDir(Hashtable<String, String> confTable){
	    return getParameter(confTable,"ALIGNMENT_DIR");
    }

    private static String getJobsDir(Hashtable<String, String> confTable){
        return getParameter(confTable,"JOBS_DIR");
    }

    private static String getReference(Hashtable<String, String> confTable){
	    return getParameter(confTable,"REF_FILE");
    }

    private static boolean isIncludeQuake(Hashtable<String, String> confTable){
	    String result= getParameter(confTable,"INCLUDE_QUAKE").toLowerCase();

        return result.equals("true");
    }

    private static boolean isIncludeBwa(Hashtable<String, String> confTable){
	
        String result= getParameter(confTable,"INCLUDE_BWA").toLowerCase();

        return result.equals("true");
    }
    
    private static boolean isIncludeLast(Hashtable<String, String> confTable){

        String result= getParameter(confTable,"INCLUDE_LAST").toLowerCase();

        return result.equals("true");
    }

    private static boolean isIncludeNovoalign(Hashtable<String, String> confTable){
        String result= getParameter(confTable,"INCLUDE_NOVOALIGN").toLowerCase();
        //return Boolean.valueOf(result);
        return result.equals("true");
    }
    
    private static boolean isIncludeBowtie(Hashtable<String, String> confTable){

        String result= getParameter(confTable,"INCLUDE_BOWTIE2").toLowerCase();

        return result.equals("true");
    }

    private static String getBowtieRef(Hashtable<String, String> confTable){
	    return getParameter(confTable,"BOWTIE_REF");
    }

    private static String getBwaRef(Hashtable<String, String> confTable){
	    return getParameter(confTable,"BWA_REF");
    }

    private static String getLastRef(Hashtable<String, String> confTable){
	    return getParameter(confTable,"LAST_REF");
    }

    private static String getNovoalignRef(Hashtable<String, String> confTable){
	    return confTable.get("NOVOALIGN_REF");
    }

    private static String[] getIsolates(Hashtable<String, String> confTable)throws FileNotFoundException{
	    String isolatesFileName = getParameter(confTable,"ID_LIST_FILE");

	    Scanner isFile = new Scanner(new File(isolatesFileName));

        ArrayList<String> isolateList = new ArrayList<String>();

        while(isFile.hasNextLine()){
            String curLine = isFile.nextLine();
            String [] curLineArray = curLine.split("\\s+");
            for(String curItem : curLineArray){
                   isolateList.add(curItem);
            }
        }

        String [] result = new String[isolateList.size()];
        for(int i=0; i<isolateList.size(); i++){
            result[i] = isolateList.get(i);
        }

        return result;
    }

	private static String getParameter(Hashtable<String, String> confTable, String parameterKey){
		if(confTable.keySet().contains(parameterKey)){
			return confTable.get(parameterKey);
		}else{
			System.out.println("ERROR: "+parameterKey+" IS NOT SET");
			System.exit(1);
		}
		return "ERROR!!!!";
	}

}