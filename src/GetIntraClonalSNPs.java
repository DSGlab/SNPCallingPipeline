import JavaBrew.Utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Scanner;

/**
 * This method complements GetHQSNPs when the reference employed is not from an
 * organism from the same population. Tnis method takes al the SNPs and finds
 * which ones do not differ in at least one isolate and returns a list of only
 * the SNPs that segregate in the population.
 *
 * This is STEP ONE.ONE from the SNP calling pipeline. It takes the output from
 * GetHQSNPs and it returns only intraclonal SNPs. The output should be reviewed
 * in STEP TWO (SNPChecker).
 *
 * @author juliofdiazc
 */
public class GetIntraClonalSNPs {
    //private static final String CONF_FILE = Utilities.CONF_FILE;
    private static final String DIR_SEPARATOR = Utilities.DIR_SEPARATOR;
    private static final String VCF_NAME = Utilities.VCF_NAME;

    private static final String ALIGNER_STRING = "ALIGNER";

    private static final String ID_LABEL_SEPARATOR = Utilities.ID_LABEL_SEPARATOR;
    private static final String REPLICON_POS_SEPARATOR = Utilities.REPLICON_POS_SEPARATOR;

    private static final String ID_LIST_FILE_STRING = "ID_LIST_FILE";
    private static final String INTRA_SNP_WORKING_DIR_STRING = "INTRA_SNP_WORKING_DIR";
    private static final String INTRA_SNP_INPUT_STRING = "INTRA_SNP_INPUT";
    private static final String INTRA_SNP_OUTPUT_STRING = "INTRA_SNP_OUTPUT";

    public static void main ( String[] args ) throws Exception {
        //Hashtable<String,String> options = Utilities.getConfigurationTable(CONF_FILE);
        Hashtable<String,String> options = Utilities.getConfigurationTable(args[0]);

        /* We look for SNPs only among the isolates with the following ids */
        String[] IDS = Utilities.getIds(options.get( ID_LIST_FILE_STRING ));

        /* The directory where the folders containing the vcf files are found */
        String WORKING_DIR = options.get( INTRA_SNP_WORKING_DIR_STRING );

        /* The complete list of HQ SNPs */
        String INPUT = options.get( INTRA_SNP_INPUT_STRING );

        /* The new list of intra samples HQ SNPs */
        String OUTPUT = options.get( INTRA_SNP_OUTPUT_STRING );

        /* The name of the aligner to use */
        String ALIGNER = options.get( ALIGNER_STRING );

        System.out.println( "1. Importing list of HQ SNPs." );
        ArrayList<String> HQsnpList = getSNPs( INPUT );
        System.out.println( "2. Importing raw SNP calls." );
        LinkedHashMap<String, VCFHolder> snpCalls = getAllSNPCalls( WORKING_DIR, VCF_NAME, IDS, ALIGNER );

        System.out.println("3. Listing only intraclonal HQ SNPs.");

        PrintWriter out = new PrintWriter( OUTPUT );
        for(String snp: HQsnpList){
            String[] temp = snp.split( REPLICON_POS_SEPARATOR );
            String tempContig = temp[0];
            int tempPos = Integer.parseInt( temp[1] );

            int count = 0;
            for(String id:snpCalls.keySet()){

                if( snpCalls.get(id).isVariant( tempContig, tempPos ) ){
                    count++;
                }
            }
            if ( count!= IDS.length && count!=0 ) {
                out.println( snp );
            }
        }
        out.close();

        System.out.println( "4. List of intraclonal HQ SNPs was written in: " + OUTPUT );
        System.out.println();
    }

    /**
     * This method reads the vcf file of each isolate in the list of ids. Each vcf file
     * needs to be in the following format:
     * {working directory}{id}_{aligner}/{name for vcf files}
     *
     * @param prefix beginning of file
     * @param suffix end of file
     * @param ids list of isolate ids
     * @return list of all the SNP calls
     * @throws FileNotFoundException if file not found
     */
    private static LinkedHashMap<String, VCFHolder> getAllSNPCalls(String prefix,
                                                                  String suffix, String[] ids, String aligner)
            throws FileNotFoundException {

        LinkedHashMap<String, VCFHolder> snpCalls = new LinkedHashMap<String, VCFHolder>();

        for ( String id: ids ) {
            VCFHolder vcfh = new VCFHolder(prefix + id + ID_LABEL_SEPARATOR + aligner + DIR_SEPARATOR + suffix);
            snpCalls.put(id, vcfh);
        }
        return snpCalls;
    }

    /**
     * This method reads a file containing positions that will be reviewed and it returns
     * an ArrayList with them. The file needs to be a single text file where each line is
     * a variant position in the following format:
     * {contig/chromosome}-{position in contig/chromosome}
     *
     * The contig/chromosome identifier must be identical to the CHROM column from the vcf
     * format and the position identifier must be identical to the POS column from the vcf
     * format.
     *
     * @param fileName the name of the file containing the SNP positions that need to
     *                 be reviewed.
     * @return a list of all positions that need to be reviewed
     * @throws FileNotFoundException if file is not found
     */
    private static ArrayList<String> getSNPs ( String fileName )
            throws FileNotFoundException {
        Scanner in = new Scanner( new File( fileName ));

        ArrayList<String> result = new ArrayList<String>();
        while ( in.hasNextLine() ) {
            result.add(in.nextLine());
        }

        return result;
    }
}
