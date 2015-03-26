import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;

/**
 * Modified by juliofdiaz on 06/11/14.
 *
 * This object acquires the data from each variant in a VCF file. The vcf
 * file will be read by VCFHolder and each line become this object.
 *
 */
public class VCFVariant {
    private String chromosome;
    private Integer position;
    private String reference;
    private String alternative;
    private Double quality;
    private String type;
    private Integer depth;
    private Integer qualityDepth;
    private Integer dp4ReferenceForward;
    private Integer dp4ReferenceReverse;
    private Integer dp4AlternativeForward;
    private Integer dp4AlternativeReverse;

    public VCFVariant ( String VCFLine ) {
        String []data = VCFLine.split("\t");
        this.chromosome = data[0];
        this.position = Integer.parseInt( data[1] );
        this.reference = data[3];
        this.alternative = data[4];
        this.quality = Double.parseDouble( data[5] );
        setInfoStep1(data[7]);
        this.qualityDepth = dp4AlternativeForward + dp4AlternativeReverse +
                dp4ReferenceForward + dp4ReferenceReverse;
    }

    public Integer getQualityDepth () {
        return this.qualityDepth;
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

    public Integer getDp4ReferenceForward () {
        return this.dp4ReferenceForward;
    }

    public Integer getDp4ReferenceReverse() {
        return this.dp4ReferenceReverse;
    }

    public Integer getDp4AlternativeForward () {
        return this.dp4AlternativeForward;
    }

    public Integer getDp4AlternativeReverse () {
        return this.dp4AlternativeReverse;
    }

    /**
     * This method is ONLY used by the class constructor and it sets the type of variant
     * and variant read information.
     *
     * @param infoLine the line of the variant in the vcf file
     */
    private void setInfoStep1( String infoLine ){
        ArrayList<String> info = new ArrayList<String>( Arrays.asList( infoLine.split(";") ) );

        if( info.get(0).equals("INDEL") ){
            this.type = "INDEL";
            info.remove(0);
            setInfoStep2( info );
        }else{
            this.type = "SNP";
            setInfoStep2( info );
        }
    }

    /**
     * This method is ONLY used by the setInfoStep1() in the class constructor and it
     * extracts variant read information.
     *
     * @param infoArray the portion of the variant line that contains the read information
     */
    private void setInfoStep2 ( List<String> infoArray ) {
        Hashtable<String,String> ht = new Hashtable<String, String>();
        for ( String s : infoArray ) {
            ht.put( s.split("=")[0], s.split("=")[1]);
        }

        this.depth = Integer.parseInt( ht.get( "DP" ) );
        this.dp4ReferenceForward = Integer.parseInt( ht.get( "DP4" ).split(",")[0] );
        this.dp4ReferenceReverse = Integer.parseInt(ht.get("DP4").split(",")[1]);
        this.dp4AlternativeForward = Integer.parseInt( ht.get( "DP4" ).split(",")[2] );
        this.dp4AlternativeReverse = Integer.parseInt(ht.get("DP4").split(",")[3]);
    }

    /**
     * This method calculates the ratio of the number of reads supporting the reference
     * to the number of reads supporting the variant.
     *
     * @return reference reads divided by variant reads
     */
    public Double getReferenceToAlternativeRatio () {
        Double ref = this.dp4ReferenceReverse + this.dp4ReferenceForward + 0.0;
        Double alt = this.dp4AlternativeForward + this.dp4AlternativeReverse + 0.0;
        return ref / alt;
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