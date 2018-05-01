import JavaBrew.FastaPrinter;
import JavaBrew.Utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
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
        String VARIANT_FILE = Utilities.HOME_PATH+"/Dropbox/CF/snp_annotation/DOLOSA_3j/variants.txt";

        /*
         * The prefix of the mview files. It may include the directory.
         */
        String MVIEW_FILE_PREFIX = Utilities.HOME_PATH+"/Dropbox/CF/snp_annotation/DOLOSA_3j/";

        /*
         * The suffix of the mview files. This could be changed but for
         * simplification purposes should remain the same.
         */
        String MVIEW_FILE_SUFFIX = "_mview.fa";

        /*
         * The list of ids to be assessed. If the ids are in a numerical
         * range, then the static getRange() method can be used
         */
        ArrayList<String> IDS_ARRAY = Utilities.getRange(1,517);
        //ArrayList<String> IDS_ARRAY = new ArrayList<String>(Arrays.asList( new String[]{"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283","284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299","300","301","302","303","304","305","306","307","308","309","310","311","312","313","314","315","316","317","318","319","320","321","322","323","324","325","326","327","328","329","330","331","332","333","334","335","336","337","338","339","340","341","342","343","344","345","346","347","348","349","350","351","352","353","355","356","357","358","359","360","361","362","363","364","365","366","367","368","369","370","371","372","373","374","375","376","377","378","379","380","381","382","383","384","386","387","388","389","390","391","392","393","394","395","396","397","398","399","400","401","402","403","404","405","406","407","408","409","410","411","412","413","414","415","416","417","418","419","420","422","423","424","425","426","427","428","429","430","431","432","433","434","435","436","437","438","439","440","441","442","443","444","445","446","447","448","449","450","451","453","454","455","456","457","458","459","460","461","462","463","464","465","466","467","468","469","470","471","472","473","474","475","476","477","478","479","480","481","482","483","484","485","486","487","488","489","490","491","492","493","494","495","496","497","498","499","501","503","504","505","506","507","508","509","510","515","516","517"} ));


        /*
         * The location and name of the file where the output will be printed.
         */
        String OUT_FILE = Utilities.HOME_PATH+"/Dropbox/CF/snp_annotation/DOLOSA_3j/tbd.txt";


        process(VARIANT_FILE,MVIEW_FILE_PREFIX,MVIEW_FILE_SUFFIX,IDS_ARRAY,OUT_FILE);

    }

    /**
     *
     *
     * @param varFile           The location of the file with the variant information.
     * @param mviewFilePrefix   The directory where the mview files is located.
     * @param mviewFileSuffix   The text after the isolate id of the mview file.
     * @param idsArray          The list of variant ids in this analysis.
     * @param outFile           The location and file name of the output file.
     * @throws FileNotFoundException
     */
    public static void process(String varFile, String mviewFilePrefix, String mviewFileSuffix,
                               ArrayList<String> idsArray, String outFile ) throws FileNotFoundException {

        PrintWriter out = new PrintWriter(outFile);

        out.println("MULTIPLE_HITS:hit contig, hit pos, hit start in contig, hit end in contig," +
                " homologous pos in hit, homologous base in hit, bitscore, evalue;NEXT HIT");
        out.println("ID\tCONTIG\tPOS\tVAR STRAND\tQUERY START\tQUERY END\tPOS (0 index)\t" +
                "REF BASE\tHIT CONTIG\tHIT POS (1 indexed, approx)\tHIT STRAND\tHIT START IN CONTIG\t" +
                "HIT END IN CONTIG\tHOMOLOGOUS POS IN HIT (0 index)\tHOMOLOGOUS BASE IN HIT\tNOTE");

        LinkedHashMap<String, Variant> variants = getVariants( varFile );
        for (String id : idsArray) {
            String curFile = mviewFilePrefix + id + mviewFileSuffix;
            BlastMatch result = new BlastMatch();
            System.out.println(id);
            /* THE result OBJECT IS THE ONE HOLDING THE DESIRED INFORMATION*/
            result = work(variants.get(id), getMViewObject(curFile));

            /* START OF PRINTING INFORMATION OF BLAST RESULT */
            out.print(result.getVariantId() + "\t");
            out.print(result.getVariantContig() + "\t");
            out.print(result.getVariantPos() + "\t");
            out.print(result.getVariantStrand() + "\t");
            out.print(result.getHitQueryStart() + "\t");
            out.print(result.getHitQueryStop() + "\t");
            out.print(result.getVariantPosInNucleotideSeq() + "\t");
            out.print(result.getVariantRef() + "\t");
            out.print(result.getHitContig() + "\t");


            //GOOD BUT 1 INDEXED. APPROXIMATE
            if (result.getHitStrand() == '+') {
                out.print((result.getHitStart() + result.getVariantHomologousPos()) + "\t");
            } else {
                out.print((result.getHitStart() - result.getVariantHomologousPos()) + "\t");
            }

            out.print(result.getHitStrand() + "\t");
            out.print(result.getHitStart() + "\t");
            out.print(result.getHitStop() + "\t");
            out.print(result.getVariantHomologousPos() + "\t");
            out.print(result.getVariantHomologousBase() + "\t");
            out.println(result.getNote());

            /* END OF PRINTING INFORMATION */
        }

        out.close();
    }

    /**
     * This program takes a variant and the respective blast hit of its annotation
     * and it returns a BlastMatch object with the proper information about the hit
     * and the respective information about the homologous base and position of the variant.
     *
     * @param var           The variant object.
     * @param mViewObject   The mview object.
     * @return              The resulting object encapsulating the information about
     *                      the blast hit an the respective homologous positions to the
     *                      variant.
     */
    private static BlastMatch work( Variant var, MViewObject mViewObject  ){
        BlastMatch main = new BlastMatch();
        //System.out.println(var.getId());
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

    /**
     * This method complements a base if the reference is in the '-' strand or
     * leaves it as is if the strand is '+' or if the location is intergenic ('n').
     *
     * @param base      The input base to be tested.
     * @param strand    The direction of the strand
     * @return          If the strand is '+' or 'n', the it returns the same
     *                  base. If the base is '-', then it returns the complement.
     */
    private static char properStrand(char base, char strand){
        if(strand=='+' || strand=='n'){
            return base;
        }else{
            return Utilities.getComplementaryBase(base);
        }
    }

    /**
     * This method is aimed at find the position of the homologous position
     * of the variant in the contig that holds the hit.
     *
     * @param hit   The hit where we find the homologous position from the
     *              variant.
     * @param pos   The position where we find the variant in the annotation.
     * @return      The position where we find the homologous position in
     *              the contig of the hit.
     */
    private static int getPosInVariantContig(MViewHits hit, int pos){
        int start;
        if( hit.getStart()>hit.getStop() ){
            start = hit.getStop();
        }else{
            start = hit.getStart();
        }
        return start+pos;
    }

    /**
     * This method turns the information of a list of Hits in a String.
     *
     * @param hits  The hits whose information is to be returned.
     * @param pos   The position of the variant in the hits.
     * @return      The information of al the input hits, including the base at
     *              the variant position in the hit.
     */
    private static String printAllHits(ArrayList<MViewHits>hits, int pos){
        String result = "";
        for(MViewHits tempHit:hits){
            result = result+tempHit.getContig()+",";

            if (tempHit.getStrand() == '+') {
                result = result+(tempHit.getStart() + pos);
            } else {
                result = result+(tempHit.getStart() - pos);
            }

            result = result+","+tempHit.getStart()+","+tempHit.getStop()+","+pos+","+
                     tempHit.getNucleotide().charAt(pos)+","+tempHit.getBitscore()+
                     ","+tempHit.geteValue()+";";
        }
        return result;
    }

    /**
     * Ideally this method should compare all the hits and return the MViewHits
     * that is the best hit for the query.
     *
     * @param hits  A list of all the hits to be compared.
     * @return      The best hit.
     */
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

    /**
     * The class holds all the output information.
     */
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
     * @return          The MView object retrieved from the mview fasta file.
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
