import JavaBrew.FastaPrinter;
import JavaBrew.Utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by juliofdiaz on 3/28/16.
 *
 * @author juliofdiazc
 */
public class SNPIdentityFromBlast {

    public static void main(String[] args) throws FileNotFoundException {

        /*
         * The file containing the variant file.
         */
        String VARIANT_FILE = "/Users/juliofdiaz/Dropbox/CF/snp_annotation/variants.txt";

        /*
         * The prefix of the mview files. It may include the directory.
         */
        String MVIEW_FILE_PREFIX = "/Users/juliofdiaz/Dropbox/CF/snp_annotation/";

        /*
         * The suffix of the mview files. This could be changed but for
         * simplification purposes should remain the same.
         */
        String MVIEW_FILE_SUFFIX = "_mview.fa";

        /*
         * The list of ids to be assesed. If the ids are in a numerical
         * range, then the static getRange() method can be used
         */
        ArrayList<String> IDS_ARRAY = Utilities.getRange(1,1936);


        LinkedHashMap<String, Variant> variants = getVariants( VARIANT_FILE );
        for(String id : IDS_ARRAY){
            String curFile = MVIEW_FILE_PREFIX + id + MVIEW_FILE_SUFFIX;
            BlastMatch result = new BlastMatch();

            /* THE result OBJECT IS THE ONE HOLDING THE DESIRED INFORMATION*/
            result = work( variants.get( id ), getMViewObject( curFile ) );

            /* START OF PRINTING INFORMATION OF BLAST RESULT */
            System.out.print(result.getVariantId()+"\t");
            System.out.print(result.getVariantContig()+"\t");
            System.out.print(result.getVariantPos()+"\t");
            System.out.print(result.getVariantStrand()+"\t");
            System.out.print(result.getHitQueryStart()+"\t");
            System.out.print(result.getHitQueryStop()+"\t");
            System.out.print(result.getVariantPosInNucleotideSeq()+"\t");
            System.out.print(result.getVariantRef()+"\t");
            System.out.print(result.getHitContig()+"\t");


            //GOOD BUT 1 INDEXED. APPROXIMATE
            if(result.getHitStrand()=='+') {
                System.out.print( ( result.getHitStart() + result.getVariantHomologousPos() ) + "\t");
            }else {
                System.out.print( ( result.getHitStart() - result.getVariantHomologousPos() )+"\t");
            }

            System.out.print(result.getHitStrand()+"\t");
            System.out.print(result.getHitStart()+"\t");
            System.out.print(result.getHitStop()+"\t");
            System.out.print(result.getVariantHomologousPos()+"\t");
            System.out.print(result.getVariantHomologousBase()+"\t");
            System.out.println(result.getNote());

            /* END OF PRINTING INFORMATION */
        }


    }

    /**
     *
     *
     * @param var the variant
     * @param mViewObject the mview object
     * @return
     */
    private static BlastMatch work( Variant var, MViewObject mViewObject  ){
        BlastMatch main = new BlastMatch();
        main.setVariantId( var.getId() );
        main.setVariantContig( var.getContig() );
        main.setVariantPos( var.getPos() );
        main.setVariantStrand( var.getStrand() );
        main.setVariantRef( var.getRef() );
        main.setVariantAlt( var.getAlt() );
        main.setVariantBaseInNucleotideSeq( var.getRefInNucleotideSeq() );
        main.setVariantPosInNucleotideSeq( var.getPosInNucleotideSeq() );
        main.setVariantNucleotideSeq( var.getNucleotideSeq() );

        if(mViewObject!=null) {

            if(var.getPosInNucleotideSeq()<mViewObject.getRefStop() && var.getPosInNucleotideSeq()>=(mViewObject.getRefStart()-1)) {
                Integer newPos = (var.getPosInNucleotideSeq() - mViewObject.getRefStart() + 1) ;
                ArrayList<MViewHits> curHits = mViewObject.getHits();
                if( mViewObject.getHits().size()==1 ){
                    MViewHits firstHit = curHits.get(0);
                    String firstNucleotideSeq = firstHit.getNucleotide();
                    Character newBase = properStrand(firstNucleotideSeq.charAt(newPos), var.getStrand());

                    main.setHitContig( firstHit.getContig() );
                    main.setHitNucleotideSeq( firstHit.getNucleotide() );
                    main.setHitStart( firstHit.getStart() );
                    main.setHitStop( firstHit.getStop() );
                    main.setHitEValue( firstHit.geteValue() );
                    main.setHitBitscore( firstHit.getBitscore() );
                    main.setHitStrand( firstHit.getStrand() );
                    main.setHitQueryStart( firstHit.getRefStart() );
                    main.setHitQueryStop( firstHit.getRefStop() );

                    main.setVariantHomologousPos( newPos );
                    main.setVariantHomologousBase( newBase );

                }else{
                    // This portion needs improvement. This portion is set
                    // for when the query has multiple hits. Currently, it
                    // returns empty hit information and for the homologous
                    // base it returns '-' for gap. The actual information
                    // is in the notes.

                    MViewHits bestHit = getBestHit( curHits );

                    main.setHitContig( "-" );
                    main.setHitNucleotideSeq( "-" );
                    main.setHitStart( -1 );
                    main.setHitStop( -1 );
                    main.setHitEValue( -1.0 );
                    main.setHitBitscore( -1.0 );
                    main.setHitStrand( 'n' );
                    main.setHitQueryStart( curHits.get(0).getRefStart() );
                    main.setHitQueryStop( curHits.get(0).getRefStop() );

                    main.setVariantHomologousPos( -1 );
                    main.setVariantHomologousBase( '-' );

                    main.setNote( "MULTIPLE_HITS:"+ printAllHits(curHits, newPos) );
                }

            }else{
                // If the blast result does not include the region where the
                // the SNP is called, the return all empty hit information.
                // Homologous base returns '-' in order to represent a gap.
                main.setHitContig( "-" );
                main.setHitNucleotideSeq( "-" );
                main.setHitStart( -1 );
                main.setHitStop( -1 );
                main.setHitEValue( -1.0 );
                main.setHitBitscore( -1.0 );
                main.setHitStrand( 'n' );
                main.setHitQueryStart( -1 );
                main.setHitQueryStop( -1 );

                main.setVariantHomologousPos( -1 );
                main.setVariantHomologousBase( '-' );

                main.setNote(" ");
            }
        }else{
            // If there was no blast result, then return all empty hit
            // information. Also homologous base returns '-' to represent
            // gap.
            main.setHitContig( "-" );
            main.setHitNucleotideSeq( "-" );
            main.setHitStart( -1 );
            main.setHitStop( -1 );
            main.setHitEValue( -1.0 );
            main.setHitBitscore( -1.0 );
            main.setHitStrand( 'n' );
            main.setHitQueryStart( -1 );
            main.setHitQueryStop( -1 );

            main.setVariantHomologousPos( -1 );
            main.setVariantHomologousBase( '-' );

            main.setNote("NO_BLAST_RESULT");
        }
        return main;
    }


    private static char properStrand(char base, char strand){
        if(strand=='+' || strand=='n'){
            return base;
        }else{
            if(base=='A'){
                return 'T';
            }else if(base=='T'){
                return 'A';
            }else if (base=='G'){
                return 'C';
            }else if(base=='C'){
                return 'G';
            }else if(base=='-'){
                return '-';
            }else{
                return 'N';
            }
        }
    }

    private static int getPosInVariantContig(MViewHits hit, int pos){
        int start;
        if( hit.getStart()>hit.getStop() ){
            start = hit.getStop();
        }else{
            start = hit.getStart();
        }
        return start+pos;
    }

    private static String printAllHits(ArrayList<MViewHits>hits, int pos){
        String result = "";
        for(MViewHits tempHit:hits){
            result = result+tempHit.getContig()+","+(tempHit.getStart()+pos)+","+tempHit.getStart()+","+tempHit.getStop()+","+pos+","+tempHit.getNucleotide().charAt(pos)+";";
        }
        return result;
    }

    private static MViewHits getBestHit( ArrayList<MViewHits> hits ){
        Double highestBitScore = 0.0;
        MViewHits bitScoreHit = new MViewHits();
        Double lowestEValue = 1000.0;
        MViewHits eValueHit = new MViewHits();

        for(MViewHits curhit:hits){
            if(curhit.getBitscore()>highestBitScore){
                highestBitScore = curhit.getBitscore();
                bitScoreHit = curhit;
            }

            if(curhit.geteValue()<lowestEValue){
                lowestEValue = curhit.geteValue();
                eValueHit = curhit;
            }
        }

        if( !eValueHit.equals(bitScoreHit) ){
            System.out.println("fux");
        }

        return null;
    }

    private static class BlastMatch{
        private String variantId;
        private String variantContig;
        private Integer variantPos;
        private Character variantStrand;
        private Character variantRef;
        private Character variantAlt;
        private Character variantBaseInNucleotideSeq;
        private Integer variantPosInNucleotideSeq;
        private String variantNucleotideSeq;
        private String hitContig;
        private String hitNucleotideSeq;
        private Integer hitStart;
        private Integer hitStop;
        private Double hitEValue;
        private Double hitBitscore;
        private Character hitStrand;
        private Integer hitQueryStart;
        private Integer hitQueryStop;
        private Integer variantHomologousPos;
        private Character variantHomologousBase;
        private String note;

        public BlastMatch(){
            this.variantId = null;
            this.variantContig = null;
            this.variantPos = null;
            this.variantStrand = null;
            this.variantRef = null;
            this.variantAlt = null;
            this.variantBaseInNucleotideSeq = null;
            this.variantPosInNucleotideSeq = null;
            this.variantNucleotideSeq = null;
            this.hitContig = null;
            this.hitNucleotideSeq = null;
            this.hitStart = null;
            this.hitStop = null;
            this.hitEValue = null;
            this.hitBitscore = null;
            this.hitStrand = null;
            this.hitQueryStart = null;
            this.hitQueryStop = null;
            this.variantHomologousPos = null;
            this.variantHomologousBase = null;
            this.note = null;
        }

        public String getVariantContig() {
            return variantContig;
        }

        public void setVariantContig(String variantContig) {
            this.variantContig = variantContig;
        }

        public String getVariantId() {
            return variantId;
        }

        public void setVariantId(String variantId) {
            this.variantId = variantId;
        }

        public Integer getVariantPos() {
            return variantPos;
        }

        public void setVariantPos(Integer variantPos) {
            this.variantPos = variantPos;
        }

        public Character getVariantStrand() {
            return variantStrand;
        }

        public void setVariantStrand(Character variantStrand) {
            this.variantStrand = variantStrand;
        }

        public Character getVariantRef() {
            return variantRef;
        }

        public void setVariantRef(Character variantRef) {
            this.variantRef = variantRef;
        }

        public Character getVariantAlt() {
            return variantAlt;
        }

        public void setVariantAlt(Character variantAlt) {
            this.variantAlt = variantAlt;
        }

        public Character getVariantBaseInNucleotideSeq() {
            return variantBaseInNucleotideSeq;
        }

        public void setVariantBaseInNucleotideSeq(Character variantBaseInNucleotideSeq) {
            this.variantBaseInNucleotideSeq = variantBaseInNucleotideSeq;
        }

        public Integer getVariantPosInNucleotideSeq() {
            return variantPosInNucleotideSeq;
        }

        public void setVariantPosInNucleotideSeq(Integer variantPosInNucleotideSeq) {
            this.variantPosInNucleotideSeq = variantPosInNucleotideSeq;
        }

        public String getVariantNucleotideSeq() {
            return variantNucleotideSeq;
        }

        public void setVariantNucleotideSeq(String variantNucleotideSeq) {
            this.variantNucleotideSeq = variantNucleotideSeq;
        }

        public String getHitContig() {
            return hitContig;
        }

        public void setHitContig(String hitContig) {
            this.hitContig = hitContig;
        }

        public String getHitNucleotideSeq() {
            return hitNucleotideSeq;
        }

        public void setHitNucleotideSeq(String hitNucleotideSeq) {
            this.hitNucleotideSeq = hitNucleotideSeq;
        }

        public Integer getHitStart() {
            return hitStart;
        }

        public void setHitStart(Integer hitStart) {
            this.hitStart = hitStart;
        }

        public Integer getHitStop() {
            return hitStop;
        }

        public void setHitStop(Integer hitStop) {
            this.hitStop = hitStop;
        }

        public Double getHitEValue() {
            return hitEValue;
        }

        public void setHitEValue(Double hitEValue) {
            this.hitEValue = hitEValue;
        }

        public Double getHitBitscore() {
            return hitBitscore;
        }

        public void setHitBitscore(Double hitBitscore) {
            this.hitBitscore = hitBitscore;
        }

        public Character getHitStrand() {
            return hitStrand;
        }

        public void setHitStrand(Character hitStrand) {
            this.hitStrand = hitStrand;
        }

        public Integer getHitQueryStart() {
            return hitQueryStart;
        }

        public void setHitQueryStart(Integer hitQueryStart) {
            this.hitQueryStart = hitQueryStart;
        }

        public Integer getHitQueryStop() {
            return hitQueryStop;
        }

        public void setHitQueryStop(Integer hitQueryStop) {
            this.hitQueryStop = hitQueryStop;
        }

        public Integer getVariantHomologousPos() {
            return variantHomologousPos;
        }

        public void setVariantHomologousPos(Integer variantHomologousPos) {
            this.variantHomologousPos = variantHomologousPos;
        }

        public Character getVariantHomologousBase() {
            return variantHomologousBase;
        }

        public void setVariantHomologousBase(Character variantHomologousBase) {
            this.variantHomologousBase = variantHomologousBase;
        }

        public String getNote() {
            return note;
        }

        public void setNote(String note) {
            this.note = note;
        }
    }

    /**
     * This method retrieves the blast result formatted in fasta with MView.
     * See Notebook 201604.
     *
     * @param fileName  The name of the fasta formatted blast result. This file
     *                  name should be in the format $id"-mview.fa where $i is
     *                  the isolate id.
     * @return
     */
    private static MViewObject getMViewObject( String fileName ){
        MViewObject result = new MViewObject();
        ArrayList<MViewHits> hits = new ArrayList<MViewHits>();

        File file = new File(fileName);
        LinkedHashMap<String, String> fasta = new FastaPrinter(file, false).getSequences();
        if( fasta.size()==0 ){
            return null;
        }
        if( fasta.size()==0 ){
            System.out.println( "error beep boop" );
            return null;
        }

        //get first item of LinkedHashMap
        String key = fasta.keySet().iterator().next();
        String value = fasta.get(key);
        fasta.remove(key);

        String[] curRefPos = key.split("\\s+")[6].split(":");
        result.setRefStart( Integer.parseInt( curRefPos[0] ) );
        result.setRefStop( Integer.parseInt( curRefPos[1] ) );
        result.setRefNucleotide( value );

        for( String curKey : fasta.keySet() ){
            MViewHits curHit = new MViewHits();
            String[] curInfo = curKey.split("\\s+");
            curHit.setContig( curInfo[0] );
            curHit.setBitscore( Double.parseDouble( curInfo[1] ) );
            curHit.seteValue( Double.parseDouble( curInfo[2] ) );
            curHit.setRefStrand( curInfo[4].charAt(0) );
            curHit.setStrand( curInfo[5].charAt(0) );
            String[] curInternalRefPos = curInfo[6].split(":");
            curHit.setRefStart( Integer.parseInt(curInternalRefPos[0]) );
            curHit.setRefStop( Integer.parseInt(curInternalRefPos[1]) );
            String[] curPos = curInfo[7].split(":");
            curHit.setStart( Integer.parseInt( curPos[0] ) );
            curHit.setStop( Integer.parseInt( curPos[1] ) );
            curHit.setNucleotide( fasta.get(curKey) );
            /*System.out.println( curKey );
            System.out.println("contig: "+curHit.getContig());
            System.out.println("bitscore: "+curHit.getBitscore());
            System.out.println("evalue: "+curHit.geteValue());
            System.out.println("ref strand: "+curHit.getRefStrand());
            System.out.println("strand: "+curHit.getStrand());
            System.out.println("ref Start: "+curHit.getRefStart());
            System.out.println("ref Stop: "+curHit.getRefStop());
            System.out.println("start: "+curHit.getStart());
            System.out.println("stop: "+curHit.getStop());
            System.out.println("fasta array: "+fasta);
            System.out.println("nucleotide: "+fasta.get(curKey));
            System.out.println("nucleotide: "+curHit.getNucleotide());*/
            hits.add(curHit);
        }
        result.setHits(hits);

        return result;
    }

    /**
     * This class holds information from the blast result formatted in
     * fasta alignment format and created with MView.
     */
    private static class MViewObject{
        private Integer refStop;
        private Integer refStart;
        private String refNucleotide;
        private ArrayList<MViewHits> hits;

        public  MViewObject(){}

        public Integer getRefStart() {
            return refStart;
        }

        public void setRefStart(Integer refStart) {
            this.refStart = refStart;
        }

        public String getRefNucleotide() {
            return refNucleotide;
        }

        public void setRefNucleotide(String refNucleotide) {
            this.refNucleotide = refNucleotide;
        }

        public ArrayList<MViewHits> getHits() {
            return hits;
        }

        public void setHits(ArrayList<MViewHits> hits) {
            this.hits = hits;
        }

        public Integer getRefStop() {
            return refStop;
        }

        public void setRefStop(Integer refStop) {
            this.refStop = refStop;
        }
    }

    /**
     * This class holds information about the hits from the blast result
     * formatted in fasta alignment and created with MView.
     */
    private static class MViewHits{
        private String contig;
        private String nucleotide;
        private Integer start;
        private Integer stop;
        private Double eValue;
        private Double bitscore;
        private Character refStrand;
        private Character strand;
        private Integer refStart;
        private Integer refStop;

        public MViewHits(){}


        public String getContig() {
            return contig;
        }

        public void setContig(String contig) {
            this.contig = contig;
        }

        public Integer getStart() {
            return start;
        }

        public void setStart(Integer start) {
            this.start = start;
        }

        public Double geteValue() {
            return eValue;
        }

        public void seteValue(Double eValue) {
            this.eValue = eValue;
        }

        public Character getRefStrand() {
            return refStrand;
        }

        public void setRefStrand(Character refStrand) {
            this.refStrand = refStrand;
        }

        public String getNucleotide() {
            return nucleotide;
        }

        public void setNucleotide(String nucleotide) {
            this.nucleotide = nucleotide;
        }


        public Double getBitscore() {
            return bitscore;
        }

        public void setBitscore(Double bitscore) {
            this.bitscore = bitscore;
        }

        public Integer getStop() {
            return stop;
        }

        public void setStop(Integer stop) {
            this.stop = stop;
        }

        public Character getStrand() {
            return strand;
        }

        public void setStrand(Character strand) {
            this.strand = strand;
        }

        public Integer getRefStart() {
            return refStart;
        }

        public void setRefStart(Integer refStart) {
            this.refStart = refStart;
        }

        public Integer getRefStop() {
            return refStop;
        }

        public void setRefStop(Integer refStop) {
            this.refStop = refStop;
        }

        public boolean equals(MViewHits other){
            if(this.bitscore==other.getBitscore()){
                if(this.eValue==other.geteValue()){
                    if(this.contig==other.getContig()) {
                        if(this.strand==other.getStrand()) {
                            if(this.start==other.getStart()) {
                                if(this.stop==other.getStop()) {
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
            return false;
        }

        public String printHitInfo(){
            return this.contig+" "+this.bitscore+" "+this.eValue+" "+
                    "#"+" "+this.refStrand+" "+this.strand+" "+this.refStart+":"+
                    this.refStop+" "+this.start+":"+this.stop;
        }
    }

    /**
     * This method extracts the variants in a variant file created by the user.
     * See Notebook 201604.
     *
     * @param fileName  The name of the file containing the information about the
     *                  variants.
     * @return          A HashMap containing variants. The key is the id of the
     *                  variants and the target is the variant information.
     * @throws FileNotFoundException
     */
    private static LinkedHashMap<String, Variant> getVariants(String fileName)
            throws FileNotFoundException {
        LinkedHashMap<String, Variant> result = new LinkedHashMap<String, Variant>();
        Scanner in = new Scanner(new File(fileName));

        in.nextLine();
        while( in.hasNextLine() ){
            Variant curVar = new Variant();
            String line = in.nextLine();
            String[] curInfo = line.split("\\s+");
            curVar.setId( curInfo[0] );
            curVar.setContig( curInfo[1] );
            curVar.setPos( Integer.parseInt( curInfo[2] ) );
            curVar.setStrand( curInfo[3].charAt(0) );
            curVar.setRef( curInfo[4].charAt(0) );
            curVar.setAlt( curInfo[5].charAt(0) );
            curVar.setRefInNucleotideSeq( curInfo[6].charAt(0) );
            curVar.setPosInNucleotideSeq( Integer.parseInt( curInfo[7] ) ) ;
            curVar.setNucleotideSeq( curInfo[8] );

            result.put( curInfo[0], curVar);
        }
        return result;
    }

    /**
     * This class hold the information about variants.
     */
    private static class Variant{
        private String id;
        private String contig;
        private Integer pos;
        private Character strand;
        private Character ref;
        private Character alt;
        private Character refInNucleotideSeq;
        private Integer posInNucleotideSeq;
        private  String nucleotideSeq;

        public Variant(){}

        public String getContig() {
            return contig;
        }

        public void setContig(String contig) {
            this.contig = contig;
        }

        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
        }

        public Integer getPos() {
            return pos;
        }

        public void setPos(Integer pos) {
            this.pos = pos;
        }

        public Character getStrand() {
            return strand;
        }

        public void setStrand(Character strand) {
            this.strand = strand;
        }

        public Character getRef() {
            return ref;
        }

        public void setRef(Character ref) {
            this.ref = ref;
        }

        public Character getAlt() {
            return alt;
        }

        public void setAlt(Character alt) {
            this.alt = alt;
        }

        public Character getRefInNucleotideSeq() {
            return refInNucleotideSeq;
        }

        public void setRefInNucleotideSeq(Character refInNucleotideSeq) {
            this.refInNucleotideSeq = refInNucleotideSeq;
        }

        public Integer getPosInNucleotideSeq() {
            return posInNucleotideSeq;
        }

        public void setPosInNucleotideSeq(Integer posInNucleotideSeq) {
            this.posInNucleotideSeq = posInNucleotideSeq;
        }

        public String getNucleotideSeq() {
            return nucleotideSeq;
        }

        public void setNucleotideSeq(String nucleotideSeq) {
            this.nucleotideSeq = nucleotideSeq;
        }
    }

}
