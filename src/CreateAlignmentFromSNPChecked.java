import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by juliofdiaz on 5/20/15.
 *
 * This program takes a list of SNPs in the format specified in
 * doc/Sample-SNP_CHECKED_FILE.tx
 *
 *  @author juliofdiazc
 */
public class CreateAlignmentFromSNPChecked {
    private static final String FASTA_START = ">";
    private static final String REF_HEADER = ">REF";
    private static final String CONTIG_POS_SEPARATOR = "-";

    public static void main (String[] args) throws FileNotFoundException {
        String SNP_CHECKED_FILE =  "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_NEW/SNPChecked_BWA.txt";
        LinkedHashMap<String, ArrayList<VCFVariant>> variants = getVariantsByParentId(SNP_CHECKED_FILE);
        LinkedHashMap <String, String> ref = getReference( getAllVariants(SNP_CHECKED_FILE) );
        for( String isol : variants.keySet() ){
            System.out.println( FASTA_START + isol );
            for( String r : ref.keySet() ){
                VCFVariant tempVar = contains(  variants.get( isol ), r );
                if( tempVar != null ) {
                    System.out.print( tempVar.getAlternative() );

                }else{
                    System.out.print( ref.get(r) );
                }
            }
            System.out.println();
        }
        System.out.println( REF_HEADER );
        for ( String s : ref.keySet() ) {
            System.out.print( ref.get( s ) );
        }


    }

    private static VCFVariant contains( ArrayList<VCFVariant> haystack, String needle ){
        for ( VCFVariant var : haystack ) {
            if ( needle.equals( var.getChromosome() + CONTIG_POS_SEPARATOR + var.getPosition() ) ){
                return var;
            }
        }
        return null;
    }

    /**
     *
     * @return ArrayList
     * @throws FileNotFoundException
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

        //this may not be necessary. it may just be needed if file has header
        in.nextLine();
        while ( in.hasNext() ) {
            VCFVariant tempVar = new VCFVariant();
            String[] temp = in.nextLine().split( "\t" );
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
     * @param variants a l
     * @return LinkedHashMap
     */
    private static LinkedHashMap<String, String> getReference ( ArrayList<VCFVariant> variants ){
        ArrayList<String> temp = new ArrayList<String>();
        LinkedHashMap<String, String> result = new LinkedHashMap<String, String>();

        for ( VCFVariant s : variants ) {
            if ( !temp.contains( s.getChromosome() + CONTIG_POS_SEPARATOR + s.getPosition() )  ) {
                temp.add( s.getChromosome() + CONTIG_POS_SEPARATOR + s.getPosition() );
                result.put( s.getChromosome() + CONTIG_POS_SEPARATOR + s.getPosition(), s.getReference() );
            }
        }

        return sortReferenceItems( result );
    }

    /**
     * This method needs to be better annotated
     *
     * @param ref               reference
     * @return LinkedHashMap
     */
    private static LinkedHashMap<String,String> sortReferenceItems ( LinkedHashMap<String,String> ref ) {
        LinkedHashMap<String,String> result = new LinkedHashMap<String, String>();

        Comparator<String> comparator = new Comparator<String>() {
            @Override
            public int compare(String o1, String o2) {
                String[] temp = o1.split("-");
                String[] temp2 = temp[0].split("_");
                Integer chrom1 = Integer.parseInt(temp2[2]);
                Integer pos1 = Integer.parseInt( temp[1] );

                String[] temp22 = o2.split("-");
                String[] temp222 = temp22[0].split("_");
                Integer chrom2 = Integer.parseInt(temp222[2]);
                Integer pos2 = Integer.parseInt( temp22[1] );

                if( chrom1 > chrom2 ) {
                    return 1;
                } else if (chrom1.equals(chrom2)) {
                    if ( pos1 > pos2 ) {
                        return 1;
                    } else if (pos1 < pos2) {
                        return -1;
                    } else {
                        return 0;
                    }
                } else {
                    return -1;
                }


            }
        };
        SortedSet<String> keys = new TreeSet<String>(comparator);


        keys.addAll( ref.keySet() );

        for(String s: keys){
            result.put(s, ref.get(s));
            System.out.println(s);
        }

        return result;
    }



}
