package JavaBrew;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;

/**
 * This program takes a list of sequences in the fasta format
 * and prints each sequence in a separate fasta file
 *
 * @author juliofdiazc
 */
public class FastaSplitter {
    public static void main ( String [] args ) {
        String FASTA_FILE = "/Users/juliofdiaz/Dropbox/CF/snp_annotation/SNP_annotNucleotideSeq.fa";
        run( FASTA_FILE );
    }

    public static void run ( String fileName ) {
        File fastaFile = new File( fileName );
        LinkedHashMap<String, String> lhm = new FastaPrinter(fastaFile).getSequences();
        String path = getPath( fileName );

        for ( String key : lhm.keySet() ) {
            PrintWriter out = null;
            try {
                out = new PrintWriter( path + key + ".fa" );
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
            assert out != null;
            out.println( ">"+key );
            out.println( lhm.get( key ) );
            out.close();
        }
    }


    private static String getPath( String fileName ){
        File f = new File(fileName);
        return f.getParent()+"/";
    }
}
