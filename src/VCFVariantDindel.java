import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

/**
 * Modified by juliofdiaz on 06/11/14.
 *
 * This class acquires the data from each variant in a VCF file. The vcf
 * file will be read by VCFHolder and each line become this object.
 *
 * @author julio.diaz@mail.utoronto.ca
 *
 */
public class VCFVariantDindel {
    private String chromosome;
    private Integer position;
    private String reference;
    private String alternative;
    private String filter;
    private Double quality;
    private String type;
    private Integer depth;
    private Integer nfs;
    private Integer nrs;
    private Integer nf;
    private Integer nr;
    private Double af;

    public VCFVariantDindel ( String VCFLine ) {
        String []data = VCFLine.split("\t");
        this.chromosome = data[0];
        this.position = Integer.parseInt( data[1] );
        this.reference = data[3];
        this.alternative = data[4];
        this.quality = Double.parseDouble( data[5] );
        this.filter = data[6];
        this.type = "INDEL";
        setInfoStep1(data[7]);
    }

    public String getChromosome () {
        return this.chromosome;
    }

    public Integer getPosition () {
        return this.position;
    }

    public String getReference () {
        return this.reference;
    }

    public String getAlternative () {
        return this.alternative;
    }

    public Double getQuality () {
        return this.quality;
    }

    public String getType () {
        return this.type;
    }

    public Integer getDepth () {
        return this.depth;
    }

    public String getFilter () {
        return this.filter;
    }

    public int getNfs () {
        return this.nfs;
    }

    public int getNrs () {
        return  this.nrs;
    }

    public int getNf () {
        return this.nf;
    }

    public int getNr () {
        return this.nr;
    }

    public Double getAf () {
        return this.af;
    }


    /**
     * This method is ONLY used by the class constructor and it sets the type of variant
     * and variant read information.
     *
     * @param infoLine the line of the variant in the vcf file
     */
    private void setInfoStep1( String infoLine ){
        ArrayList<String> info = new ArrayList<String>( Arrays.asList( infoLine.split(";") ) );

        Hashtable<String,String> ht = new Hashtable<String, String>();
        for ( String s : info ) {
            ht.put( s.split("=")[0], s.split("=")[1]);
        }

        this.depth = Integer.parseInt( ht.get( "DP" ) );
        this.nfs = Integer.parseInt( ht.get( "NFS" ) );
        this.nrs = Integer.parseInt( ht.get("NRS") );
        this.nf = Integer.parseInt( ht.get("NF") );
        this.nr = Integer.parseInt( ht.get("NR") );
        this.af = Double.parseDouble( ht.get("AF") );
    }

    /**
     * This method tests if the position of the variant is "N" in the reference
     *
     * @return true is reference is "N"
     */
    public Boolean isReferenceN(){

        ArrayList<Character> chars= new ArrayList<Character>();
        for ( int i=0; i<this.reference.length(); i++ ) {
            chars.add( this.reference.charAt( i ) );
        }

        return chars.contains('N');
    }

    /**
     * This method tests if variant is labeled as "INDEL"
     *
     * @return true if variant is an indel
     */
    public Boolean isIndel(){
        return this.type.equals( "INDEL" );
    }
}