import java.io.FileNotFoundException;
import java.io.File;
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

    	for(String ISOLATE : ISOLATE_LIST){
	        System.out.println("\tPrinting jobs for isolate:\t"+ISOLATE);
            if(includeBowtie){
	    		Create.createBowtieJobs.run(INPUT_DIR, BOWTIE_REF, REFERENCE, OUTPUT_DIR, ISOLATE, false, JOBS_DIR);
	        }
	        if(includeBwa){
                Create.createBwaJobs.run(INPUT_DIR, BWA_REF, REFERENCE, OUTPUT_DIR, ISOLATE, false, JOBS_DIR);
	        }
            if(includeLast){
                Create.createLastJobs.run(INPUT_DIR, LAST_REF, REFERENCE, OUTPUT_DIR, ISOLATE, false, JOBS_DIR);
	        }
	        if(includeNovoalign){
                Create.createNovoalignJobs.run(INPUT_DIR, NOVO_REF, REFERENCE, OUTPUT_DIR, ISOLATE, false, JOBS_DIR);
	        }

	        if(includeQuake){
		        if(includeBowtie){
                    Create.createBowtieJobs.run(INPUT_DIR, BOWTIE_REF, REFERENCE, OUTPUT_DIR, ISOLATE, true, JOBS_DIR);
		        }
		        if(includeBwa){
                    Create.createBwaJobs.run(INPUT_DIR, BWA_REF, REFERENCE, OUTPUT_DIR, ISOLATE, true, JOBS_DIR);
		        }
		        if(includeLast){
                    Create.createLastJobs.run(INPUT_DIR, LAST_REF, REFERENCE, OUTPUT_DIR, ISOLATE, true, JOBS_DIR);
		        }
		        if(includeNovoalign){
                    Create.createNovoalignJobs.run(INPUT_DIR, NOVO_REF, REFERENCE, OUTPUT_DIR, ISOLATE, true, JOBS_DIR);
		        }
	        }
	    }
	    System.out.println();
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
		        //System.out.println(testFileName+"\t"+testFileName.length()+"\t"+curFile.getName().substring(0,testFileName.length()));
		        //System.out.println( curFile.getName()/*.substring(0,testFileName.length()-1)*/ );
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

    private static String getInputDir(Hashtable<String, String> confTable){
	    return confTable.get("INPUT_DIR");
    }

    private static String getOutputDir(Hashtable<String, String> confTable){
	    return confTable.get("ALIGNMENT_DIR");
    }

    private static String getJobsDir(Hashtable<String, String> confTable){
        return confTable.get("JOBS_DIR");
    }

    private static String getReference(Hashtable<String, String> confTable){
	    return confTable.get("REF_FILE");
    }

    private static boolean isIncludeQuake(Hashtable<String, String> confTable){

	    String result= confTable.get("INCLUDE_QUAKE").toLowerCase();

        return result.equals("true");
    }

    private static boolean isIncludeBwa(Hashtable<String, String> confTable){
	
        String result= confTable.get("INCLUDE_BWA").toLowerCase();

        return result.equals("true");
    }
    
    private static boolean isIncludeLast(Hashtable<String, String> confTable){

        String result= confTable.get("INCLUDE_LAST").toLowerCase();

        return result.equals("true");
    }

    private static boolean isIncludeNovoalign(Hashtable<String, String> confTable){
        String result= confTable.get("INCLUDE_NOVOALIGN").toLowerCase();
        return Boolean.valueOf(result);
        //return result.equals("true");
    }
    
    private static boolean isIncludeBowtie(Hashtable<String, String> confTable){

        String result= confTable.get("INCLUDE_BOWTIE2").toLowerCase();

        return result.equals("true");
    }

    private static String getBowtieRef(Hashtable<String, String> confTable){
	    return confTable.get("BOWTIE_REF");
    }

    private static String getBwaRef(Hashtable<String, String> confTable){
	    return confTable.get("BWA_REF");
    }

    private static String getLastRef(Hashtable<String, String> confTable){
	    return confTable.get("LAST_REF");
    }

    private static String getNovoalignRef(Hashtable<String, String> confTable){
	    return confTable.get("NOVOALIGN_REF");
    }

    private static String[] getIsolates(Hashtable<String, String> confTable)throws FileNotFoundException{
	    String isolatesFileName = confTable.get("ID_LIST_FILE");

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

}