import JavaBrew.FastaPrinter;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 * Modified by juliofdiaz on 2/26/15.
 *
 *
 * @author julio.diaz@mail.utoronto.ca
 *
 */
public class GetHQSNPs {
    public static void main ( String[] args ) throws Exception {


        String SUFFIX = "/r.vcf";
        String WORKING_DIR = checkDir("/Users/juliofdiaz/Dropbox/CF/snp_calling/CF67_C71/");
        String REF_FILE = "/Users/juliofdiaz/Dropbox/CF/references/C71.fa";
        String id = "1";

        System.out.println( "ITEM: " + id );
        run( REF_FILE, WORKING_DIR, id, SUFFIX );
    }

    /**
     * This code is used to extract positions where high quality (HQ) SNPs were found. It employs
     * variant calls from six different read mapping methods ( i.e. bwa, last and novoalign
     * using quality filtered reads and raw reads ). The loci where HQ SNPs were found are labeled
     * as {reference contig}-{position}
     *
     * @param refFile the path and name of the file containing the reference in fasta file
     * @param workngDir the working directory where the vcf file are located
     * @param id the identifier of the isolate
     * @param suffix the suffix of the file of each vcf file ( e.g. ".vcf" )
     * @throws Exception
     */
    private static void run ( String refFile, String workngDir, String id, String suffix )
            throws Exception {

        File file = new File( refFile );
        FastaPrinter ref = new FastaPrinter( file );

        String[] varFiles = { workngDir + "BWA-" + id + suffix,
                              workngDir + "BWA-COR_" + id + suffix,
                              workngDir + "LAST-" + id + suffix,
                              workngDir + "LAST-COR_" + id + suffix,
                              workngDir + "NOVOALIGN-" + id + suffix,
                              workngDir + "NOVOALIGN-COR_" + id + suffix };

        LinkedHashMap<String, ArrayList<VCFVariant>> main = VCFHolderMap.getVariantsMap(varFiles);

        /* Here we filter out the variants whose quality is lower than threshold */
        main = VCFHolderMap.qualityFilterVariantsMap(main, 30);

        /* Here we filter out the variants whose sequencing depth is lower than threshold */
        main = VCFHolderMap.qualityDepthFilterVariantsMap(main, 20);

        /* Here we filter out the variants that are closer to the contig ends than threshold */
        main = VCFHolderMap.contigEndFilterVariantsMap(main, ref, 250);

        /* Here we filter out the variants that are covered unevenly by the forward and reverse reads */
        main = VCFHolderMap.readBalanceFilterVariantsMap(main, 5);

        /* Here we filter out the variants where sufficient read support the reference */
        main = VCFHolderMap.refToAltRatioFilterVariantsMap(main, 0.1);

        /*  Here we filter out the variants that are clustered too close*/
        main = VCFHolderMap.clusterFilterVariantsMap(main, 15);


        /* Loops through each position that reports a variant */
        for ( String key : main.keySet()) {

            /* Ensures that position actually reports variant. IOW, variants were not pruned by filtering steps */
            if ( main.get(key).size() != 0 ) {

                /* Filters out positions reporting variants where the reference includes "N" */
                if ( !main.get(key).get(0).isReferenceN() ) {

                     /* Filters out indels */
                    if ( !main.get(key).get(0).isIndel() ) {

                        /* Filters out positins that are not reported by all the methods */
                        if ( main.get(key).size() == 6 ) {
                            System.out.print( key + "\t" );
                            System.out.println();
                        }
                    }
                }
            }
        }

    }

    /**
     * This method analyses the string entered as a directory and ensures its in fact a
     * directory. Otherwise it throws an error
     *
     * @param dirString the string referring to the directory
     * @return the String that is formatted as a directory
     */
    private static String checkDir( String dirString ) {
        File dir = new File ( dirString );
        if ( dir.isDirectory() ) {
            return dir.getPath()+"/";
        } else {
            return "";
        }
    }

}
