import JavaBrew.FastaPrinter;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 * Modified by juliofdiaz on 2/26/15.
 *
 * This class implements the VCFHolderMap class. It creates a VCFHolderMap employing vcf files
 * from six different read mapping methods ( i.e. bwa, last and novoalign using quality
 * filtered reads and raw reads ) and it filters the variant calls to meet HQ status. This is
 * performed for all isolates in our sampling cohort to identify all positions with signal for
 * HQ SNPs.
 *
 * @author julio.diaz@mail.utoronto.ca
 *
 */
public class GetHQSNPs {
    private static final String DIR_SEPARATOR = "/";

    public static void main ( String[] args ) throws Exception {
        ArrayList<String> snpPositions = new ArrayList<String>();

        //First Set of data from CF170
        //String[] ids = {"1","2","3","4","5","6C","7","8","9","10","11","12C","13","14","15","16","17","18","19","20"};
        String[] ids = {"1a","1b","1c","1d","1e","1f","1g","1h","1i","1j","2a","2b","2c","2d","2e","2f","2g","2h","2i","2j","3a","3b","3c","3d","3e","3f","3g","3h","3i","3j","4a","4b","4c","4d","4e","4f","4g","4h","4i","4j","5a","5b","5c","5d","5e","5f","5g","5h","5i","5j","6a","6b","6c","6d","6e","6f","6g","6h","6i","6j","7a","7b","7c","7d","7e","7f","7g","7h","7i","7j","8a","8b","8c","8d","8e","8f","8g","8h","8i","8j"};
        //String[] ids = {"1c","1d","1e","1f","1g","1h","1i","1j","2c","2d","2e","2f","2g","2h","2i","2j","3c","3d","3e","3f","3g","3h","3i","3j","4c","4d","4e","4f","4g","4h","4i","4j","5c","5d","5e","5f","5g","5h","5i","5j","6c","6d","6e","6f","6g","6h","6i","6j","7c","7d","7e","7f","7g","7h","7i","7j","8c","8d","8e","8f","8g","8h","8i","8j","9c","9d","9e","9f","9g","9h","9i","9j","10c","10d","10e","10f","10g","10h","10i","10j"};

        for( String id : ids ) {
            String SUFFIX = "/r.vcf";
            //System.out.println("\t"+id);
            //String WORKING_DIR = checkDir("/Users/juliofdiaz/Dropbox/CF/snp_calling/CF67_C71/");
            //String REF_FILE = "/Users/juliofdiaz/Dropbox/CF/references/C71.fa";

            //String WORKING_DIR = checkDir("/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_B6CQ/");
            //String REF_FILE = "/Users/juliofdiaz/Dropbox/CF/references/B6CQ.fa";

            String WORKING_DIR = checkDir("/Users/juliofdiaz/Dropbox/CF/snp_calling/DOLOSA_V21-INTRA/");
            String REF_FILE = "/Users/juliofdiaz/Dropbox/CF/references/8f.fa";

            //String WORKING_DIR = checkDir("/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_NEW/");
            //String REF_FILE = "/Users/juliofdiaz/Dropbox/CF/references/B6CQ.fa";


            ArrayList<String> temp = run(REF_FILE, WORKING_DIR, id, SUFFIX);
            //System.out.println( temp );

            for ( String t : temp ) {
                if ( !snpPositions.contains(t) ) {
                    snpPositions.add(t);
                    System.out.println( t );
                }
            }
        }

        System.out.println( snpPositions.size() );
    }

    /**
     *
     * @param refFile
     * @param workingDir
     * @param id
     * @param suffix
     * @throws Exception
     */
    public static void printOnScreen ( String refFile, String workingDir, String id, String suffix )
            throws Exception {
        ArrayList<String> snps = run(refFile, workingDir, id, suffix);
        for ( String s : snps ) {
            System.out.println(s);
        }
    }

    /**
     * This code is used to extract positions where high quality (HQ) SNPs were found. It employs
     * variant calls from six different read mapping methods ( i.e. bwa, last and novoalign
     * using quality filtered reads and raw reads ). The loci where HQ SNPs were found are labeled
     * as {reference contig}-{position}
     *
     * @param refFile the path and name of the file containing the reference in fasta file
     * @param workingDir the working directory where the vcf file are located
     * @param id the identifier of the isolate
     * @param suffix the suffix of the file of each vcf file ( e.g. ".vcf" )
     * @throws Exception
     * @note The ArrayList should eventually hold VCFVariants
     */
    private static ArrayList<String> run ( String refFile, String workingDir, String id, String suffix )
            throws Exception {
        ArrayList<String> result = new ArrayList<String>();
        File file = new File( refFile );
        FastaPrinter ref = new FastaPrinter( file );

        String[] varFiles = {
                              workingDir + id + "_BWA" + suffix,
                              //workingDir + id + "-COR_BWA" + suffix,
                              workingDir + id + "_LAST" + suffix,
                              //workingDir + id + "-COR_LAST" + suffix,
                              workingDir + id + "_NOVOALIGN" + suffix,
                              //workingDir + id + "-COR_NOVOALIGN" + suffix
                            };
        LinkedHashMap<String, ArrayList<VCFVariant>> main = VCFHolderMap.getVariantsMap(varFiles);

        /* Here we filter out the variants whose quality is lower than threshold */
        main = VCFHolderMap.qualityFilterVariantsMap(main, 30);


        for ( String key : main.keySet() ) {
            for(VCFVariant temp: main.get(key)){
                if( temp.getQuality() < 30 ) {
                    System.out.println(temp.getQuality());
                }
            }
        }


        /* Here we filter out the variants whose sequencing depth is lower than threshold */
        main = VCFHolderMap.qualityDepthFilterVariantsMap(main, 20);

        /* Here we filter out the variants that are closer to the contig ends than threshold */
        main = VCFHolderMap.contigEndFilterVariantsMap(main, ref, 250);

        /* Here we filter out the variants that are covered unevenly by the forward and reverse reads */
        main = VCFHolderMap.readBalanceFilterVariantsMap(main, 3);

        /* Here we filter out the variants where sufficient reads support the reference */
        main = VCFHolderMap.refToAltRatioFilterVariantsMap(main, 0.2);

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
                        if ( main.get(key).size() == 3 ) {
                            result.add(key);
                        }
                    }
                }
            }
        }
        return result;
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
            return dir.getPath() + DIR_SEPARATOR;
        } else {
            return "";
        }
    }

}
