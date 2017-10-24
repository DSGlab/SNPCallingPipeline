import JavaBrew.Utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by juliofdiaz on 5/20/15.
 *
 * This program takes a list of SNPs in the format specified in
 * doc/Sample-SNP_FILTERED_FILE.txt; then it returns a sorted list
 * and a fasta formatted alignment of the SNPs identified.
 *
 * This is STEP THREE of the SNP calling pipeline. It takes a reviewed version
 * of the output of STEP TWO (SNPChecker).
 *c
 *  @author juliofdiazc
 */
public class CreateAlignmentFromSNPChecked {
    //private static final String CONF_FILE = Utilities.CONF_FILE;

    private static final String FASTA_START = Utilities.FASTA_START;
    private static final String REF_HEADER = FASTA_START + Utilities.REFERENCE_DEFAULT_NAME;
    private static final String REPLICON_POS_SEPARATOR = Utilities.REPLICON_POS_SEPARATOR;
    private static final String ID_LABEL_SEPARATOR = Utilities.ID_LABEL_SEPARATOR;
    private static final String SEPARATOR = Utilities.SEPARATOR;

    private static final String CREATE_ALIGN_INPUT_STRING= "CREATE_ALIGN_INPUT";
    private static final String CREATE_ALIGN_OUTPUT_STRING = "CREATE_ALIGN_OUTPUT";
    private static final String CREATE_ALIGN_LIST_OUTPUT_STRING = "CREATE_ALIGN_LIST_OUTPUT";

    public static void main (String[] args) throws FileNotFoundException {
        //Hashtable<String,String> options = Utilities.getConfigurationTable(CONF_FILE);
        Hashtable<String,String> options = Utilities.getConfigurationTable(args[0]);

        //String SNP_CHECKED_FILE =  "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_NEW/SNPChecked_BWA.txt";
        //String SNP_CHECKED_FILE =  "/Users/juliofdiaz/Dropbox/CF/indel_calling/CF170_NEW/IndelChecked_BWA.txt";
        /* */
        String INPUT =  options.get( CREATE_ALIGN_INPUT_STRING );//"/Users/juliofdiaz/Dropbox/CF/snp_calling/TEST/SNPFiltered_BWA.txt";

        /* */
        String OUTPUT = options.get( CREATE_ALIGN_OUTPUT_STRING );//"/Users/juliofdiaz/Dropbox/CF/snp_calling/TEST/SNPAlignment.fa";

        /* */
        String OUTPUT_LIST = options.get( CREATE_ALIGN_LIST_OUTPUT_STRING );//"/Users/juliofdiaz/Dropbox/CF/snp_calling/TEST/FinalSNPList.txt";

        /* */
        System.out.println( "1. Retrieving SNP calls" );
        LinkedHashMap<String, ArrayList<VCFVariant>> variants = getVariantsByParentId(INPUT);


        /* */
        System.out.println( "2. Collecting reference sequence" );
        LinkedHashMap <String, String> ref = getReference( getAllVariants(INPUT), OUTPUT_LIST );

        /* */
        System.out.println( "3. Recording SNP alignment in: " + OUTPUT );
        PrintWriter out = new PrintWriter(OUTPUT);
        for( String isol : variants.keySet() ){
            System.out.println("\tWriting:\t"+isol);
            out.println( FASTA_START + isol );
            for( String r : ref.keySet() ){
                VCFVariant tempVar = contains(  variants.get( isol ), r );
                if( tempVar != null ) {
                    out.print( tempVar.getAlternative() );

                }else{
                    out.print( ref.get(r) );
                }
            }
            out.println();
        }
        out.println( REF_HEADER );
        for ( String s : ref.keySet() ) {
            out.print( ref.get( s ) );
        }
        out.close();

    }

    private static VCFVariant contains( ArrayList<VCFVariant> haystack, String needle ){
        for ( VCFVariant var : haystack ) {
            if ( needle.equals( var.getChromosome() + REPLICON_POS_SEPARATOR + var.getPosition() ) ){
                return var;
            }
        }
        return null;
    }

    /**
     *
     *
     * @return ArrayList
     * @throws FileNotFoundException if input file of getAllVariants is not found.
     */
    private static LinkedHashMap<String, ArrayList<VCFVariant>> getVariantsByParentId ( String SNPCheckedFile )
            throws FileNotFoundException {
        LinkedHashMap<String,ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        ArrayList<VCFVariant> allVariants = getAllVariants( SNPCheckedFile );

        for ( VCFVariant var : allVariants ){
            if ( result.keySet().contains( var.getParentId() ) ) {
                ArrayList<VCFVariant> temp = result.get( var.getParentId() );
                temp.add( var );
                result.put( var.getParentId(), temp );
            } else {
                ArrayList<VCFVariant> temp = new ArrayList<VCFVariant>();
                temp.add( var );
                result.put( var.getParentId(), temp);
            }
        }

        return result;
    }

    /**
     * The input needs to be specified
     * @return ArrayList
     *
     */
    private static ArrayList<VCFVariant> getAllVariants ( String SNPCheckedFile ) throws FileNotFoundException {
        ArrayList<VCFVariant> result = new ArrayList<VCFVariant>();
        //Scanner in = new Scanner( new File( "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_NEW/SNPChecked_BWA.txt" ) );
        Scanner in = new Scanner( new File( SNPCheckedFile ) );

        //If the inout file has a header, turn the next line on
        //in.nextLine();
        while ( in.hasNext() ) {
            VCFVariant tempVar = new VCFVariant();
            String[] temp = in.nextLine().split( SEPARATOR );
            tempVar.setParentId(temp[1]);
            tempVar.setChromosome(temp[2]);
            tempVar.setPosition(Integer.parseInt(temp[3]));
            tempVar.setReference(temp[4]);
            tempVar.setAlternative( temp[5] );
            result.add(tempVar);
        }
        return result;
    }

    /**
     * From a list of VCFVariants, this method retrieves all variant positions and their
     * respective reference base.The list of variants should be called from the same reference.
     *
     * @param variants list of variants
     * @return LinkedHashMap map of reference identity of each SNP
     */
    private static LinkedHashMap<String, String> getReference ( ArrayList<VCFVariant> variants, String listOutput ) throws FileNotFoundException {
        ArrayList<String> temp = new ArrayList<String>();
        LinkedHashMap<String, String> result = new LinkedHashMap<String, String>();

        for ( VCFVariant s : variants ) {
            if ( !temp.contains( s.getChromosome() + REPLICON_POS_SEPARATOR + s.getPosition() )  ) {
                temp.add( s.getChromosome() + REPLICON_POS_SEPARATOR + s.getPosition() );
                result.put( s.getChromosome() + REPLICON_POS_SEPARATOR + s.getPosition(), s.getReference() );
            }
        }

        if( !isSingleReplicon(result) ) {
            return sortReferenceItems( result, listOutput );
        } else {
            return sortCompleteReferenceItems(result, listOutput);
        }
    }

    /**
     * Analyze list of items to see how many replicons are part of the analysis.
     *
     * @param result list of SNPs
     * @return true if a single replicon has SNP
     */
    private static Boolean isSingleReplicon(LinkedHashMap<String, String> result){
        ArrayList<String> replicons = new ArrayList<String>();
        for(String key:result.keySet()){
            String repl = key.split("-")[0];
            if(!replicons.contains(repl)){
                replicons.add(repl);
            }
        }
        return replicons.size() == 1;
    }

    /**
     * This method sorts SNP items according to its position in the replicon.
     *
     * @param ref the reference
     * @param output the ordered list
     * @return the ordered list
     * @throws FileNotFoundException file cannot be written
     */
    private static LinkedHashMap<String,String> sortCompleteReferenceItems ( LinkedHashMap<String,String> ref, String output ) throws FileNotFoundException {
        LinkedHashMap<String,String> result = new LinkedHashMap<String, String>();

        Comparator<String> comparator = new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                String[] temp = o1.split( REPLICON_POS_SEPARATOR );
                Integer pos1 = Integer.parseInt( temp[1] );
                String[] temp2 = o2.split( REPLICON_POS_SEPARATOR );
                Integer pos2 = Integer.parseInt( temp2[1] );

                return comparatorHelper(pos1,pos2);
            }
        };
        SortedSet<String> keys = new TreeSet<String>(comparator);

        keys.addAll( ref.keySet() );

        PrintWriter out = new PrintWriter( output );
        for(String s: keys){
            result.put(s, ref.get(s));
            out.println(s);
        }
        out.close();
        return result;
    }

    /**
     * This method sorts items first by contig and then by position in contig. This
     * method needs to be better annotated.
     *
     * @param ref reference
     * @return LinkedHashMap the list of sorted reference items
     */
    private static LinkedHashMap<String,String> sortReferenceItems ( LinkedHashMap<String,String> ref, String output ) throws FileNotFoundException {
        LinkedHashMap<String,String> result = new LinkedHashMap<String, String>();

        Comparator<String> comparator = new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                String[] temp = o1.split( REPLICON_POS_SEPARATOR );
                String[] temp2 = temp[0].split( ID_LABEL_SEPARATOR );
                Integer chrom1 = Integer.parseInt(temp2[2]);
                Integer pos1 = Integer.parseInt( temp[1] );

                String[] temp22 = o2.split( REPLICON_POS_SEPARATOR );
                String[] temp222 = temp22[0].split( ID_LABEL_SEPARATOR );
                Integer chrom2 = Integer.parseInt(temp222[2]);
                Integer pos2 = Integer.parseInt( temp22[1] );

                if( chrom1 > chrom2 ) {
                    return 1;
                } else if (chrom1.equals(chrom2)) {
                    return comparatorHelper(pos1,pos2);
                } else {
                    return -1;
                }


            }
        };
        SortedSet<String> keys = new TreeSet<String>(comparator);

        keys.addAll( ref.keySet() );

        PrintWriter out = new PrintWriter( output );
        for(String s: keys){
            result.put(s, ref.get(s));
            out.println(s);
        }
        out.close();

        return result;
    }

    /**
     * This method just helps a instantiation of comparator by returning the appropriate
     * value when giving two numbers to compare.
     *
     * @param pos1 the first number to compare
     * @param pos2 the second number to compare
     * @return 1 if pos1 is greater than pos2, -1 if pos2 is greater than pos1 and 0 if equal
     */
    private static Integer comparatorHelper(Integer pos1, Integer pos2){
        if ( pos1 > pos2 ) {
            return 1;
        } else if (pos1 < pos2) {
            return -1;
        } else {
            return 0;
        }
    }

}
