import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 2/19/15.
 *
 *
 */
public class SNPChecker {

    public static void main ( String[] args) throws FileNotFoundException {


        for( String i : getIds() ) {
            System.out.println("ITEM " + i);

            String snpList = "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/snplist2.txt";
            String vcf = "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/BWA-COR_" + i + "/r.vcf";

            getVCFVariantSubset(snpList, vcf);
        }
    }

    public static void getVCFVariantSubset ( String SNPList, String vcfFile )
            throws FileNotFoundException {
        ArrayList<String> SNPPositions = getSNPs( SNPList );
        VCFHolder vcfh = new VCFHolder( vcfFile );

        for( VCFVariant vcfv : vcfh.getVariants() ){

            if ( SNPPositions.contains( vcfv.getChromosome()+"-"+vcfv.getPosition() ) ){

                if ( vcfv.getQuality() < 25 || vcfv.getReferenceToAlternativeRatio() >0.2 ){//meetsLowerQualityThreshold( vcfv ) ) {
                    System.out.println( "WARNING "+"\t"+vcfv.getChromosome()+"\t"+vcfv.getPosition()+"\t"+vcfv.getReference()+"\t"+vcfv.getAlternative()+"\t"+vcfv.getQuality()+"\t"+vcfv.getDepth()+"\t"+vcfv.getQualityDepth()+"\t"+vcfv.getReferenceToAlternativeRatio() );
                }else{
                    System.out.println( vcfv.getChromosome()+"\t"+vcfv.getPosition()+"\t"+vcfv.getReference()+"\t"+vcfv.getAlternative()+"\t"+vcfv.getQuality()+"\t"+vcfv.getDepth()+"\t"+vcfv.getQualityDepth()+"\t"+vcfv.getReferenceToAlternativeRatio() );
                }
            }
        }
    }

    private static ArrayList<String> getSNPs ( String fileName )
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
                        //System.out.println(tempo);
                        result.add( tempo );
                    }
                }
            }
        }
        //System.out.println(result.size());
        return result;
    }
}
