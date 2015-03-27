import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Modified by juliofdiaz on 06/11/14.
 *
 * This class reads a VCF file and records every variant line
 * as a VCFVariant object.
 *
 * @author julio.diaz@mail.utoronto.ca
 *
 */
public class VCFHolder {
    private ArrayList<VCFVariant> variants;

    public ArrayList<VCFVariant> getVariants () {
        return this.variants;
    }

    public VCFHolder ( String s ) throws FileNotFoundException{
        this( new File ( s ) );
    }

    public VCFHolder ( File f ) throws FileNotFoundException {
        this.variants = new ArrayList<VCFVariant>();
        Scanner in = new Scanner( f );

        while( in.hasNextLine() ){
            String tempLine = in.nextLine();
            if( tempLine.charAt(0) != '#' ){
                this.variants.add(new VCFVariant(tempLine));
            }
        }
    }

    /**
     * This method checks if a position in a contig (or chromosome) is a variant
     *
     * @param contig the contig or chromosome to check
     * @param pos the position in the contig (or chomosome) to check
     * @return true if position is declared a variant
     */
    public Boolean isVariant ( String contig, int pos ){
        for ( VCFVariant vv : this.variants ) {
            if ( vv.getChromosome().equals(contig) && vv.getPosition()==pos ) {
                return true;
            }
        }
        return false;
    }

    /**
     * This method goes through all the variants in the holder and only returns the
     * indels
     *
     * @return ArrayList of VCFVariants that are "INDELS"
     */
    public ArrayList<VCFVariant> getIndels () {
        ArrayList<VCFVariant> result = new ArrayList<VCFVariant>();
        for(VCFVariant var: this.getVariants()){
            if ( var.isIndel() ) {
                result.add( var );
            }

        }
        return result;
    }


}
