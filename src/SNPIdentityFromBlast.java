import JavaBrew.FastaPrinter;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.logging.FileHandler;

/**
 * Created by juliofdiaz on 3/28/16.
 *
 * @author juliofdiazc
 */
public class SNPIdentityFromBlast {

    public static void main(String[] args){
        for(int i=1; i<=1936; i++){
        //for( int i=1; i<=1; i++ ) {
                System.out.println(i);
                String curFile = "/Users/juliofdiaz/Dropbox/CF/snp_annotation/" + i + "_mview.fa";
                MViewObject curMviewObj = getMViewObject( curFile );

        }
    }


    private static MViewObject getMViewObject( String fileName ){
        MViewObject result = new MViewObject();
        ArrayList<MViewHits> hits = new ArrayList<MViewHits>();

        File file = new File(fileName);
        LinkedHashMap<String, String> fasta = new FastaPrinter(file).getSequences();
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
            curHit.setNucleotide( fasta.get(key) );
            System.out.println( curKey );
            System.out.println("contig: "+curHit.getContig());
            System.out.println("bitscore: "+curHit.getBitscore());
            System.out.println("evalue: "+curHit.geteValue());
            System.out.println("ref strand: "+curHit.getRefStrand());
            System.out.println("strand: "+curHit.getStrand());
            System.out.println("ref Start: "+curHit.getRefStart());
            System.out.println("ref Stop: "+curHit.getRefStop());
            System.out.println("start: "+curHit.getStart());
            System.out.println("stop: "+curHit.getStop());

        }

        return null;
    }

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
    }
}
