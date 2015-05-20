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

        String[] isolates = {"1","2","3","4","5","6C","7","8","9","10","11","12C","13","14","15","16","17","18","19","20"};
        //String[] isolates = {"20"};

        for( String i : isolates ) {
            System.out.println("ITEM " + i);

            //String snpList = "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/snplist2.txt";
            //String vcf = "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/BWA-COR_" + i + "/r.vcf";
            String snpList = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_B6CQ/snplist.txt";
            String vcf = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_B6CQ/NOVOALIGN-COR_" + i + "/r.vcf";

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

    private static ArrayList<String> getIds () {
        ArrayList<String> result = new ArrayList<String>();
        File f = new File("/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71");

        for ( File cur : f.listFiles() ) {
            if ( cur.isDirectory() ) {
                String tempo = cur.getName().split("-")[1].split("_")[0];
                if( !tempo.equals( "COR" ) ) {
                    if ( !result.contains(tempo) ) {
                        result.add( tempo );
                    }
                }
            }
        }
        //System.out.println(result.size());
        return result;
    }
}
