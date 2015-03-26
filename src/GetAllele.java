import JavaBrew.FastaPrinter;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 * Created by juliofdiaz on 3/23/15.
 */
public class GetAllele {
    public static void main ( String[] args ) {

        String DIR = System.getProperty( "user.home" );
        FastaPrinter fp = new FastaPrinter( new File( DIR+"/Dropbox/CF67_snpAlignment_03102015-exclude(11r,6k).fa") );
        LinkedHashMap<String, String> seqs = fp.getSequences();

        for ( String s : seqs.keySet() ) {
            System.out.print( s + "\t" );
            String allele = ""+complement( run( seqs.get(s), 38 ) ) + complement( run( seqs.get(s), 37 ) ) + complement( run( seqs.get(s), 36 ) ) ;
            System.out.print(allele + "\t");

            if ( allele.equals( "GAC" ) ) {
                System.out.print( "AMaA\tancestral major allele" );
            } else if ( allele.equals( "GCC" ) ) {
                System.out.print( "AMiA\tancestral minor allele" );
            } else if ( allele.equals( "GAT" ) ) {
                System.out.print("SMiA\tsweep minor allele");
            } else if ( allele.equals( "AAC" ) ) {
                System.out.print( "SMiB\tsweep minor allele B" );
            } else {
                System.out.print( "NNnN\tambiguous" );
            }

            ArrayList<String> cladeA = getCladeA();
            if ( cladeA.contains( s ) ) {
                System.out.println( "\tcladeA" );
            } else {
                System.out.println( "\tcladeB" );
            }
        }

    }

    public static char run ( String seq, int pos ) {
        return seq.charAt( pos-1 );
    }

    private static char complement ( char c ) {
        if ( c == 'A' ) {
            return 'T';
        } else if ( c == 'T' ) {
            return 'A';
        } else if ( c == 'G' ) {
            return 'C';
        } else if ( c == 'C') {
            return 'G';
        } else {
            return 'N';
        }
    }

    private static ArrayList<String> getCladeA () {
        String[] members = new String[] {"1a","1b","1c","1d","1e","1f","1g","1h","1i","1j","1k","1l","1m","1n","1o","1p","1q","1r","1s","2b","2c","2d","2e","2r","3b","3c","3g","3h","3i","3j","3k","3l","3n","3q","3r","4a","4b","4c","4d","4e","4f","4g","4h","4i","4j","4k","4l","4p","4t","5n","6l","6m","7k","7l","7p","7r","7s","7t","8o","8p","9o","9q","9s","9t","10i","10j","10k","10l","11a","11f","11g","11h","11i","11j","11k","11l","11m","11n","11t","12c","12i"};
        ArrayList<String> result = new ArrayList<String>();
        for ( int i=0;i<members.length; i++ ) {
            result.add( members[i] );
        }
        return result;
    }
}
