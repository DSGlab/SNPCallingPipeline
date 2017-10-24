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
    //private static final String CONF_FILE = Utilities.CONF_FILE;
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
        //Hashtable<String,String> options = Utilities.getConfigurationTable(CONF_FILE);
        Hashtable<String,String> options = Utilities.getConfigurationTable(args[0]);

        /* */
        String[] IDS = Utilities.getIds(options.get( ID_LIST_FILE_STRING ));

        /* */
        String WORKING_DIR = options.get( SNP_CHECKER_WORKING_DIR_STRING );

        /* */
        String ALIGNER = options.get( ALIGNER_STRING );

        /* */
        String INPUT = options.get( SNP_CHECKER_INPUT_STRING );

        /* */
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
