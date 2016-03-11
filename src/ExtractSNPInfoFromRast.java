import java.util.ArrayList;

/**
 * Created by juliofdiaz on 3/11/16.
 *
 * @author juliofdiazc
 */
public class ExtractSNPInfoFromRast {

    public static void main (String[] args){

    }

    public static ArrayList<Annotation> getAnnotations(){
        return null;
    }

    private class Annotation{
        private String contig;
        private String feature;
        private String type;
        private String location;
        private Integer start;
        private Integer stop;
        private Character strand;
        private String function;
        private String aliases;
        private String figfam;
        private String evidenceCodes;
        private String nucleotideSeq;
        private String aminoAcidSeq;

        public Annotation (){ }

        public String getContig() {
            return contig;
        }

        public void setContig(String contig) {
            this.contig = contig;
        }

        public String getFeature() {
            return feature;
        }

        public void setFeature(String feature) {
            this.feature = feature;
        }

        public String getFunction() {
            return function;
        }

        public void setFunction(String function) {
            this.function = function;
        }

        public Character getStrand() {
            return strand;
        }

        public void setStrand(Character strand) {
            this.strand = strand;
        }

        public Integer getStart() {
            return start;
        }

        public void setStart(Integer start) {
            this.start = start;
        }

        public String getType() {
            return type;
        }

        public void setType(String type) {
            this.type = type;
        }

        public Integer getStop() {
            return stop;
        }

        public void setStop(Integer stop) {
            this.stop = stop;
        }

        public String getLocation() {
            return location;
        }

        public void setLocation(String location) {
            this.location = location;
        }

        public String getAliases() {
            return aliases;
        }

        public void setAliases(String aliases) {
            this.aliases = aliases;
        }

        public String getEvidenceCodes() {
            return evidenceCodes;
        }

        public void setEvidenceCodes(String evidenceCodes) {
            this.evidenceCodes = evidenceCodes;
        }

        public String getAminoAcidSeq() {
            return aminoAcidSeq;
        }

        public void setAminoAcidSeq(String aminoAcidSeq) {
            this.aminoAcidSeq = aminoAcidSeq;
        }

        public String getNucleotideSeq() {
            return nucleotideSeq;
        }

        public void setNucleotideSeq(String nucleotideSeq) {
            this.nucleotideSeq = nucleotideSeq;
        }

        public String getFigfam() {
            return figfam;
        }

        public void setFigfam(String figfam) {
            this.figfam = figfam;
        }
    }
}
