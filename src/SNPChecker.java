import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 2/19/15.
 *
 * This class inspects every loci provided by a SNP list and returns all the read information
 * for each isolate from the cohort at those positions.
 *
 * @author julio.diaz@mail.utoronto.ca
 *
 */
public class SNPChecker {


    public static void main ( String[] args) throws FileNotFoundException {

        //String[] isolates = {"1","2","3","4","5","6C","7","8","9","10","11","12C","13","14","15","16","17","18","19","20"};
        //String[] isolates = {"20"};
        String[] isolates = {"1a","1b","1c","1d","1e","1f","1g","1h","1i","1j","2a","2b","2c","2d","2e","2f","2g","2h","2i","2j","3a","3b","3c","3d","3e","3f","3g","3h","3i","3j","4a","4b","4c","4d","4e","4f","4g","4h","4i","4j","5a","5b","5c","5d","5e","5f","5g","5h","5i","5j","6a","6b","6c","6d","6e","6f","6g","6h","6i","6j","7a","7b","7c","7d","7e","7f","7g","7h","7i","7j","8a","8b","8c","8d","8e","8f","8g","8h","8i","8j"};

        for( String i : isolates ) {
            //System.out.println("ITEM " + i);

            //String snpList = "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/snplist2.txt";
            //String vcf = "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/BWA-COR_" + i + "/r.vcf";

            //String snpList = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_B6CQ/snplist.txt";
            //String vcf = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_B6CQ/NOVOALIGN-COR_" + i + "/r.vcf";

            String snpList ="/Users/juliofdiaz/Dropbox/CF/snp_calling/DOLOSA_AU0158/intraSNPList.txt";
            String vcf = "/Users/juliofdiaz/Dropbox/CF/snp_calling/DOLOSA_AU0158/" + i + "-COR_BWA/r.vcf";

            getVCFVariantSubset(snpList, vcf);
        }
    }

    /**
     * This method checks every position provided by the SNP list in the vcf file provided
     * and it returns information (if found) about that positions
     *
     * @param SNPList the path and name of the file containing the list of variants
     * @param vcfFile the path and name of the vcf file
     * @throws FileNotFoundException if vcf file or variant list is not found.
     */
    public static void getVCFVariantSubset ( String SNPList, String vcfFile )
            throws FileNotFoundException {
        ArrayList<String> SNPPositions = getSNPs( SNPList );
        VCFHolder vcfh = new VCFHolder( vcfFile );

        for( VCFVariant vcfv : vcfh.getVariants() ){
            if ( SNPPositions.contains( vcfv.getChromosome()+"-"+vcfv.getPosition() ) ){
                System.out.println( vcfv.getChromosome() + "\t" + vcfv.getPosition() + "\t" +
                        vcfv.getReference() + "\t" + vcfv.getAlternative() + "\t" + vcfv.getQuality() +
                        "\t" + vcfv.getDepth() + "\t" + vcfv.getQualityDepth() + "\t" +
                        "\t" + vcfv.getDp4ReferenceForward() + "\t" + vcfv.getDp4ReferenceReverse() + "\t" +
                        "\t" + vcfv.getDp4AlternativeForward() + "\t" + vcfv.getDp4AlternativeReverse() + "\t" +
                        vcfv.getReferenceToAlternativeRatio() );
            }
        }
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
    public static ArrayList<String> getSNPs ( String fileName )
            throws FileNotFoundException {
        Scanner in = new Scanner( new File( fileName ));

        ArrayList<String> result = new ArrayList<String>();
        while ( in.hasNextLine() ) {
            result.add(in.nextLine());
        }

        return result;
    }

}
