package JavaBrew;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;

/**
 * Created by juliofdiaz on 30/09/14.
 */
public class FastaPrinter {
    private LinkedHashMap<String, String> sequences;
    private int seqLength;

    public FastaPrinter () {
        this.seqLength = 60;
        this.sequences = new LinkedHashMap<String, String>();
    }

    public FastaPrinter ( LinkedHashMap<String,String> in ) {
        this.sequences = in;
        this.seqLength = 60;
    }

    public FastaPrinter ( File file ) {
        this();

        LinkedHashMap<String, DNASequence> contigs = null;
        try {
            contigs = FastaReaderHelper.readFastaDNASequence( file );
        } catch (Exception e) {
            System.err.println( "ERROR\t" + file +  " IS NOT IN FASTA FORMAT OR IT DOES NOT EXIST." );
            //e.printStackTrace();
        }
        //THIS SHOULD GO AWAY BECAUSE VALUES SHOULD BE DNASequence
        if ( contigs != null ) {
            for ( String key : contigs.keySet() ) {
                this.sequences.put( key, contigs.get( key ).getSequenceAsString() );
            }
        }


    }


    public void setSequences ( LinkedHashMap<String,String> newSequences ) {
        this.sequences = newSequences;
    }

    public LinkedHashMap<String,String> getSequences () {
        return this.sequences;
    }

    public void setSeqLength ( int newSeqLength ) {
        this.seqLength = newSeqLength;
    }

    public int getSeqLength () {
        return this.seqLength;
    }

    /**
     * This method allows the sequences in the current Fasta to be printed in the standard Fasta format
     *
     * @param outFile       The name of the newly generated Fasta File
     */
    public void simplePrint ( String outFile ) {
        PrintWriter out = null;
        try {
            out = new PrintWriter( outFile );
        } catch ( FileNotFoundException e ) {
            System.err.println( outFile+" does not exist." );
            e.printStackTrace();
        }

        if (out != null) {

            for (String key : this.sequences.keySet()) {
                out.println(">" + key);

                String tempSequence = this.sequences.get(key);
                for (int i = 0; i < tempSequence.length(); i += this.seqLength) {
                    if (i + this.seqLength < tempSequence.length()) {
                        out.println(tempSequence.substring(i, (i + this.seqLength)));
                    } else {
                        out.println(tempSequence.substring(i, tempSequence.length()));
                    }
                }
            }
            out.close();
        }
    }

}