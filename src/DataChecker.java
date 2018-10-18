import JavaBrew.FastaPrinter;
import JavaBrew.Utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.LinkedHashMap;

public class DataChecker {

    private static final String REF_FILE_STRING = "REF_FILE";
    private static final String ID_LIST_FILE_STRING = "ID_LIST_FILE";
    private static final String INPUT_DIR_STRING = "INPUT_DIR";

    private static final String BWA_PATH_STRING = "BWA_PATH";
    private static final String LAST_PATH_STRING = "LAST_PATH";
    private static final String NOVOALIGN_PATH_STRING = "NOVOALIGN_PATH";
    private static final String SAMTOOLS_PATH_STRING = "SAMTOOLS_PATH";
    private static final String BCFTOOLS_PATH_STRING = "BCFTOOLS_PATH";

    public static void main(String[] args) throws FileNotFoundException {
        Hashtable<String,String> options = Utilities.getConfigurationTable(args[0]);
        /* The reference genome in fasta file */
        String REF_FILE = options.get( REF_FILE_STRING );
        String INPUT_DIR = options.get( INPUT_DIR_STRING );

        /* Software */
        String BWA_PATH =  options.get( BWA_PATH_STRING );
        String LAST_PATH =  options.get( LAST_PATH_STRING );
        String NOVOALIGN_PATH =  options.get( NOVOALIGN_PATH_STRING );
        String SAMTOOLS_PATH =  options.get( SAMTOOLS_PATH_STRING );
        String BCFTOOLS_PATH =  options.get( BCFTOOLS_PATH_STRING );

        ArrayList<String> result = new ArrayList<String>();
        File file = new File( REF_FILE );
        FastaPrinter ref = new FastaPrinter( file );

        /* Check reference format */
        Integer refErrorCount = checkReference(ref);

        /* Check input directory */
        Integer inputDirErrorCount = isDirReal(INPUT_DIR);

        /* Check Sequencing Data is labeled in the correct format */
        Integer readErrorCount = 0;
        if(inputDirErrorCount==0) {
            String[] IDS = Utilities.getIds(options.get(ID_LIST_FILE_STRING));
            readErrorCount = checkReads(INPUT_DIR,IDS);
        }

        /* Checks software */
        Integer bwaError = checkSoftware(BWA_PATH,"bwa");
        Integer lastError = checkSoftware(LAST_PATH,"src/lastal");
        Integer lastPPError = checkSoftware(LAST_PATH,"src/last-pair-probs");
        Integer lastMCError = checkSoftware(LAST_PATH,"scripts/maf-convert");
        Integer novoalignError = checkSoftware(NOVOALIGN_PATH,"novoalign");
        Integer bcftoolsError = checkSoftware(BCFTOOLS_PATH,"bcftools");
        Integer samtoolsError = checkSoftware(SAMTOOLS_PATH,"samtools");

        /* NEEDS TO CHECK REFERENCE */


        System.out.println("\n"+ (readErrorCount + refErrorCount+inputDirErrorCount)+" error(s) found in the reference or seq read files");
        System.out.println( lastMCError+lastPPError+bwaError + lastError+novoalignError+samtoolsError+bcftoolsError+" error(s) found in the path to the software");
        System.out.println();

    }

    private static Integer checkSoftware(String bin, String label){
        File file = new File(bin+"/"+label);
        if(!file.exists()){
            System.err.println("ERROR\tThe binaries for "+label+" cannot be found.");
            return 1;
        }
        return 0;
    }

    /**
     *
     * @param dir
     * @param ids
     * @return
     */
    private static Integer checkReads(String dir, String[] ids){
        Integer errorCount = 0;
        for(String id: ids) {
            File rOne = new File(dir + "/" + id + "_1.fq");
            File rTwo = new File(dir + "/" + id + "_2.fq");

            if (!rOne.exists()) {
                System.err.println("ERROR\tRead file '"+id+"_1.fq' does not exist");
                errorCount++;
            }
            if (!rTwo.exists()) {
                System.err.println("ERROR\tRead file '"+id+"_2.fq' does not exist");
                errorCount++;
            }

        }
        return errorCount;
    }

    /**
     *
     *
     * @param dirName
     * @return
     */
    private static Integer isDirReal(String dirName){
       File dir = new File(dirName);
       if(!dir.isDirectory()){
           System.err.println("ERROR\tSequencing reads directory does not exist");
           return 1;
       }
       return 0;
    }

    /**
     * This method identifies and reports problems with the reference
     * file
     *
     * @param ref reference file
     * @return the number of errors encountered
     */
    private static Integer checkReference(FastaPrinter ref){
        Integer errorCount = 0;

        LinkedHashMap<String,String> seqs = ref.getSequences();
        ArrayList<String> headerList = new ArrayList<String>();
        for(String header:seqs.keySet()){
            String[] parts = header.split("_");
            headerList.add(parts[0]);

            /* Check that the third element of header is a number */
            /* Check proper format of header of reference fasta */
            try{
                Integer.parseInt(parts[2]);
            }catch(NumberFormatException e){
                System.err.println("ERROR\tSequence '"+header+"' must include the replicon number");
                errorCount++;
            }catch(ArrayIndexOutOfBoundsException e2){
                System.err.println("ERROR\tSequence '"+header+"' is not in the format <REFERENCE-NAME>_<REPLICON-TYPE>_<REPLICON-NUMBER>");
                errorCount++;
            }
        }
        if( !isListAllSame(headerList) ){
            System.err.println("ERROR\tNot all the headers have the same replicon name");
            errorCount++;
        }
        return errorCount;
    }

    /**
     * This method tests whether all the elements in a list are the
     * same
     *
     * @param list a list to be tested
     * @return true if all the elements are identical
     */
    private static Boolean isListAllSame(ArrayList<String> list){
        String first = list.get(0);
        int sameCount = 0;
        for(String item:list){
            if(item.equals(first)){
                sameCount = sameCount+1;
            }
        }
        return sameCount == list.size();
    }
}
