import JavaBrew.Utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 2/19/15.
 *
 * This class inspects every loci provided by a SNP list and returns all the read
 * information for each isolate from the cohort at those positions.
 *
 * This is STEP TWO of the SNP calling Pipeline. This program takes the HQ SNP positions
 * identified in STEP ONE (GetHQSNPs) [If the reference is not a bacteria within the same
 * population then get intra conal SNPs at GetIntraClonalSNPs] and reviews them in all
 * the isolates. The resulting reviewed SNP positions should be assessed to find out if
 * they may be ambiguous calls. Then, The reviewed SNP positions, can then be aligned
 * using CreateAlignmentFromSNPChecked.
 *
 * @author julio.diaz@mail.utoronto.ca
 *
 */
public class SNPChecker {
    private static final String CONF_FILE = Utilities.CONF_FILE;
    private static final String DIR_SEPARATOR = Utilities.DIR_SEPARATOR;
    private static final String VCF_NAME = Utilities.VCF_NAME;

    private static final String ALIGNER_STRING = "ALIGNER";

    private static final String ID_LABEL_SEPARATOR = Utilities.ID_LABEL_SEPARATOR;
    private static final String SEPARATOR = Utilities.SEPARATOR;

    private static final String ID_LIST_FILE_STRING = "ID_LIST_FILE";
    private static final String SNP_CHECKER_WORKING_DIR_STRING = "SNP_CHECKER_WORKING_DIR";
    private static final String SNP_CHECKER_INPUT_STRING = "SNP_CHECKER_INPUT";
    private static final String SNP_CHECKER_OUTPUT_STRING = "SNP_CHECKER_OUTPUT";

    public static void main ( String[] args) throws FileNotFoundException {
        Hashtable<String,String> options = Utilities.getConfigurationTable(CONF_FILE);

        //String[] isolates = {"1","2","3","4","5","6C","7","8","9","10","11","12C","13","14","15","16","17","18","19","20"};
        //String[] isolates = {"20"};
        //String[] isolates = {"1a","1b","1c","1d","1e","1f","1g","1h","1i","1j","2a","2b","2c","2d","2e","2f","2g","2h","2i","2j","3a","3b","3c","3d","3e","3f","3g","3h","3i","3j","4a","4b","4c","4d","4e","4f","4g","4h","4i","4j","5a","5b","5c","5d","5e","5f","5g","5h","5i","5j","6a","6b","6c","6d","6e","6f","6g","6h","6i","6j","7a","7b","7c","7d","7e","7f","7g","7h","7i","7j","8a","8b","8c","8d","8e","8f","8g","8h","8i","8j"};
        //String[] isolates = {"0a","1a","1b","1c","1d","1e","1f","1g","1h","1i","1j","2a","2b","2c","2d","2e","2f","2g","2h","2i","2j"};//,"3a","3b","3c","3d","3e","3f","3g","3h","3i","3j","4a","4b","4c","4d","4e","4f","4g","4h","4i","4j","5a","5b","5c","5d","5e","5f","5g","5h","5i","5j","6a","6b","6c","6d","6e","6f","6g","6h","6i","6j","7a","7b","7c","7d","7e","7f","7g","7h","7i","7j","8a","8b","8c","8d","8e","8f","8g","8h","8i","8j","9a","9b","9c","9d","9e","9f","9g","9h","9i","9j","10a","10b","10c","10d","10e","10f","10g","10h","10i","10j","11a","11b","11c","11d","11e","11f","11g","11h","11i","11j"};
        //String[] isolates = {"10","130","14","152","160","17","180","196","19","200","217","21","225","235","23","251","257","258","259","263","265","270","273","275","276","278","280","282","288","290","293","298C","303","306","30","310","315","317","318","323","325","326","327","330","331","332","335","336","337","340","341","342","358","359","360","366","367","369","36","375","379","37","380","385","388","390","392","393","395","400","401","402","404C","405C","406C","409","410C","411C","412","417C","418C","419C","420N2","427","428","429","430","431","432","433","446","448","449","450","453","455","457","459C","465","471","472","473","475","476","479","480","487","501","503","504","505","506","507","508","509","510","511","512","513","514","521","525","526","527","528","529","530","531","532","538","539","53","540","541-2","542","544","548","549","550","551","552","553C","557","559-2","562","563","566","567","568","569","570","571","572-1","573-2","573-3","575","577","578","580","581","583","6","7","8","97"};

        String[] IDS = Utilities.getIds(options.get( ID_LIST_FILE_STRING ));

        String WORKING_DIR = options.get( SNP_CHECKER_WORKING_DIR_STRING );

        String ALIGNER = options.get( ALIGNER_STRING );

        String INPUT = options.get( SNP_CHECKER_INPUT_STRING );

        String OUTPUT = options.get( SNP_CHECKER_OUTPUT_STRING );

        System.out.println("1. Retrieving raw SNP data.");
        PrintWriter out = new PrintWriter( OUTPUT );
        for( String i : IDS ) {
            System.out.println("\tReading:\t" + i);

            //String snpList = "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/snplist2.txt";
            //String vcf = "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/BWA-COR_" + i + "/r.vcf";

            //String snpList = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_B6CQ/snplist.txt";
            //String vcf = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_B6CQ/NOVOALIGN-COR_" + i + "/r.vcf";

            //String snpList = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_NEW/SNPList.txt";
            //String vcf = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_NEW/" + i + "_BWA/r.vcf";

            //String snpList ="/Users/juliofdiaz/Dropbox/CF/snp_calling/DOLOSA_V21-INTRA/SNPList.txt";
            //String vcf = "/Users/juliofdiaz/Dropbox/CF/snp_calling/DOLOSA_V21-INTRA/" + i + "_BWA/r.vcf";

            String vcf = WORKING_DIR + i + ID_LABEL_SEPARATOR + ALIGNER + DIR_SEPARATOR + VCF_NAME;

            ArrayList<String> checkedList = getVCFVariantSubset( INPUT, vcf );
            for(String item: checkedList){
                out.println(item);
            }
        }
        out.close();

        System.out.println( "2. Information of HQ SNPs was written in: " + OUTPUT );
        System.out.println();
    }

    /**
     * This method checks every position provided by the SNP list in the vcf file provided
     * and it returns information (if found) about that positions
     *
     * @param SNPList the path and name of the file containing the list of variants
     * @param vcfFile the path and name of the vcf file
     * @throws FileNotFoundException if vcf file or variant list is not found.
     */
    private static ArrayList<String> getVCFVariantSubset ( String SNPList, String vcfFile )
            throws FileNotFoundException {
        ArrayList<String> SNPPositions = getSNPs( SNPList );
        VCFHolder vcfh = new VCFHolder( vcfFile );

        ArrayList<String> checkedList = new ArrayList<String>();
        for( VCFVariant vcfv : vcfh.getVariants() ){
            if ( SNPPositions.contains( vcfv.getChromosome()+"-"+vcfv.getPosition() ) ){
                checkedList.add( vcfFile + SEPARATOR +
                        vcfv.getChromosome() + SEPARATOR + vcfv.getPosition() + SEPARATOR +
                        vcfv.getReference() + SEPARATOR + vcfv.getAlternative() + SEPARATOR +
                        vcfv.getQuality() + SEPARATOR + vcfv.getDepth() + SEPARATOR +
                        vcfv.getQualityDepth() + SEPARATOR + SEPARATOR +
                        vcfv.getDp4ReferenceForward() + SEPARATOR + vcfv.getDp4ReferenceReverse() +
                        SEPARATOR + SEPARATOR + vcfv.getDp4AlternativeForward() + SEPARATOR +
                        vcfv.getDp4AlternativeReverse() + SEPARATOR + vcfv.getReferenceToAlternativeRatio() );
            }
        }
        return checkedList;
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
