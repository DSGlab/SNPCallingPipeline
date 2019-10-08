import JavaBrew.FastaPrinter;
import JavaBrew.Utilities;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
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
 * This is the STEP ONE of the SNP calling pipeline. If the reference is the genome of a bacteria
 * in the same population check the HQ SNP positions in all the isolates with SNPChecker; otherwise,
 * extract the intra clonal SNPs with GetIntraClonalSNPs.
 *
 * @author julio.diaz@mail.utoronto.ca
 *
 */
public class GetHQSNPs {
    //private static final String CONF_FILE = Utilities.CONF_FILE;
    private static final String DIR_SEPARATOR = Utilities.DIR_SEPARATOR;
    private static final String VCF_NAME = Utilities.VCF_NAME;

    private static final String BWA_STRING = "BWA";
    private static final String NOVOALIGN_STRING = "NOVOALIGN";
    private static final String LAST_STRING = "LAST";
    private static final String ID_COR_MODIFIER = "-COR";
    private static final String ID_LABEL_SEPARATOR = Utilities.ID_LABEL_SEPARATOR;

    private static final String ID_LIST_FILE_STRING = "ID_LIST_FILE";
    private static final String HQ_SNP_LIST_OUTPUT_STRING = "HQ_SNP_LIST_OUTPUT";
    private static final String HQ_SNP_WORKING_DIR_STRING = "HQ_SNP_WORKING_DIR";
    private static final String REF_FILE_STRING = "REF_FILE";
    private static final String INCLUDE_PRE_QC_STRING = "INCLUDE_PRE_QC";

    private static final String HQ_SNP_QUALITY_STRING = "HQ_SNP_QUALITY";
    private static final String HQ_SNP_DEPTH_STRING = "HQ_SNP_DEPTH";
    private static final String HQ_SNP_DIST_TO_CONTIG_END_STRING = "HQ_SNP_DIST_TO_CONTIG_END";
    private static final String HQ_SNP_READ_BALANCE_STRING = "HQ_SNP_READ_BALANCE";
    private static final String HQ_SNP_REF_TO_ALT_BALANCE_STRING = "HQ_SNP_REF_TO_ALT_BALANCE";
    private static final String HQ_SNP_CLUSTER_SIZE_STRING = "HQ_SNP_CLUSTER_SIZE";
    private static final String HQ_SNP_REQUIRED_NUM_STRING = "HQ_SNP_REQUIRED_NUM";

    public static void main ( String[] args ) throws Exception {
        //Hashtable<String,String> options = Utilities.getConfigurationTable(CONF_FILE);
        Hashtable<String,String> options = Utilities.getConfigurationTable(args[0]);
        ArrayList<String> snpPositions = new ArrayList<String>();

        /* We look for HQ SNPs in the isolates with the following ids */
        String[] IDS = Utilities.getIds(options.get( ID_LIST_FILE_STRING ));
        /* The list of HQ SNPs will be found in the following file */
        String OUTPUT = options.get( HQ_SNP_LIST_OUTPUT_STRING );
        /* The directory where the folders containing the vcf files are found */
        String WORKING_DIR = checkDir(options.get( HQ_SNP_WORKING_DIR_STRING ));
        /* The reference genome in fasta file */
        String REF_FILE = options.get( REF_FILE_STRING );
        /* True if results from non-qc reads are to be included  */
        Boolean INCLUDE_PRE_QC = Boolean.valueOf(options.get( INCLUDE_PRE_QC_STRING ));

        System.out.println("1. Acquiring list of HQ SNPs.");

        PrintWriter out = new PrintWriter( OUTPUT );
        for( String id : IDS ) {
            System.out.println("\tTesting:\t"+id);

             ArrayList<String> temp = run(REF_FILE, WORKING_DIR, id, VCF_NAME, options, INCLUDE_PRE_QC);

            for ( String t : temp ) {
                if ( !snpPositions.contains(t) ) {
                    snpPositions.add(t);
                    out.println( t );
                }
            }
        }
        out.close();

        System.out.println( "2. Num. of HQ SNP:\t"+snpPositions.size() );
        System.out.println( "3. List of HQ SNPs was written in: " + OUTPUT );
        System.out.println();
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
     * @throws Exception issues with file reading or writing
     * @note The ArrayList should eventually hold VCFVariants
     */
    private static ArrayList<String> run ( String refFile, String workingDir, String id, String suffix,
                                           Hashtable<String,String> options, Boolean includePreQc)
            throws Exception {
        /* Obtain settings for HQ SNP calling */
        Integer MIN_QUALITY = Integer.parseInt( options.get( HQ_SNP_QUALITY_STRING ) );
        Integer MIN_DEPTH = Integer.parseInt( options.get( HQ_SNP_DEPTH_STRING ) );
        Integer MIN_EDGE_DISTANCE = Integer.parseInt( options.get( HQ_SNP_DIST_TO_CONTIG_END_STRING ) );
        Integer MIN_READ_BALANCE = Integer.parseInt( options.get( HQ_SNP_READ_BALANCE_STRING ) );
        Double MIN_REF_TO_ALT_RATIO = Double.parseDouble( options.get( HQ_SNP_REF_TO_ALT_BALANCE_STRING ) );
        Integer CLUSTER_SIZE = Integer.parseInt( options.get( HQ_SNP_CLUSTER_SIZE_STRING ) );
        Integer ALIGNERS_NUM = Integer.parseInt( options.get( HQ_SNP_REQUIRED_NUM_STRING ) );

        ArrayList<String> result = new ArrayList<String>();
        File file = new File( refFile );
        FastaPrinter ref = new FastaPrinter( file );

        /* Setting list of vcf files to include in the analysis */
        String[] varFiles;
        if(!includePreQc) {
            varFiles = new String[]{
                    workingDir + BWA_STRING + "-" + id + DIR_SEPARATOR + suffix,
                    ////workingDir + LAST_STRING + "-" + id + DIR_SEPARATOR + suffix,
                    ////workingDir + NOVOALIGN_STRING + "-" + id + DIR_SEPARATOR + suffix
            };
        }else {
            varFiles = new String[]{
                    workingDir + BWA_STRING + "-" + id + DIR_SEPARATOR + suffix,
                    workingDir + BWA_STRING + "-" + "COR" + "_" + id + DIR_SEPARATOR + suffix,
                    workingDir + LAST_STRING + "-" + id + DIR_SEPARATOR + suffix,
                    workingDir + LAST_STRING + "=" + "COR" + "_" + id + DIR_SEPARATOR + suffix,
                    workingDir + NOVOALIGN_STRING + "-" + id + DIR_SEPARATOR + suffix,
                    workingDir + NOVOALIGN_STRING + "-" + "COR" + "_" + id + DIR_SEPARATOR + suffix
            };
        }

        LinkedHashMap<String, ArrayList<VCFVariant>> main = VCFHolderMap.getVariantsMap(varFiles);

        /* Here we filter out the variants whose quality is lower than threshold */
        main = VCFHolderMap.qualityFilterVariantsMap(main, MIN_QUALITY);

        /* Sanity check */
        for ( String key : main.keySet() ) {
            for(VCFVariant temp: main.get(key)){
                if( temp.getQuality() < MIN_QUALITY ) {
                    System.out.println(temp.getQuality()+"what is this?");
                }
            }
        }

        /* Here we filter out the variants whose sequencing depth is lower than threshold */
        main = VCFHolderMap.qualityDepthFilterVariantsMap(main, MIN_DEPTH);
        /* Here we filter out the variants that are closer to the contig ends than threshold */
        main = VCFHolderMap.contigEndFilterVariantsMap(main, ref, MIN_EDGE_DISTANCE);
        /* Here we filter out the variants that are covered unevenly by the forward and reverse reads */
        main = VCFHolderMap.readBalanceFilterVariantsMap(main, MIN_READ_BALANCE);
        /* Here we filter out the variants where sufficient reads support the reference */
        main = VCFHolderMap.refToAltRatioFilterVariantsMap(main, MIN_REF_TO_ALT_RATIO);
        /*  Here we filter out the variants that are clustered too close*/
        main = VCFHolderMap.clusterFilterVariantsMap(main, CLUSTER_SIZE);

        /* Loops through each position that reports a variant */
        for ( String key : main.keySet()) {
            /* Ensures that position actually reports variant. IOW, variants were not pruned by filtering steps */
            if ( main.get(key).size() != 0 ) {
                /* Filters out positions reporting variants where the reference includes "N" */
                if ( !main.get(key).get(0).isReferenceN() ) {
                     /* Filters out indels */
                    if ( !main.get(key).get(0).isIndel() ) {
                        /* Filters out positions that are not reported by all the methods */
                        ////if ( main.get(key).size() == ALIGNERS_NUM ) {
                            result.add(key);
                        ////}
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