import JavaBrew.FastaPrinter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 2017-01-24.
 *
 */
public class CountIntergenicReg {
    public static void main(String[] args) throws FileNotFoundException {

        //Integer REG_REGION = Integer.parseInt("143");
        //String ANNOTATION_FILE = "/Users/juliofdiaz//Dropbox/CF/references/B6CQ.txt";
        //String REF_CONTIGS = "/Users/juliofdiaz/Dropbox/CF/references/B6CQ.fa";

        Integer REG_REGION = Integer.parseInt(args[0]);
        String ANNOTATION_FILE = args[1];
        String REF_CONTIGS = args[2];
        String OUT_FILE = args[3];

        System.out.print("Reading annotation file...");
        ArrayList<Annotation> annotList = getAnnotations(ANNOTATION_FILE);
        System.out.println("\tDONE");

        System.out.print("Reading contigs file...");
        LinkedHashMap<String, String> reference = new FastaPrinter(new File(REF_CONTIGS)).getSequences();
        System.out.println("\tDONE");

        int intCount = 0;
        int regCount = 0;

        PrintWriter out = new PrintWriter(OUT_FILE);
        for (String contig : reference.keySet()) {
            System.out.print("Analyzing contig "+contig+"...");
            for(int i=1;i<=reference.get(contig).length();i++ ){

                Variant curVariant = new Variant();
                curVariant.setContig(contig);
                curVariant.setPosition(i);

                if (isIntergenic(curVariant,annotList)) {
                    //System.out.print(contig+"\t"+i);

                Annotation prevAnnot = new Annotation();
                prevAnnot.setContig("NA");
                prevAnnot.setStrand('+');

                for (Annotation curAnnot : annotList) {
                    if (curAnnot.getContig().equals(curVariant.getContig())) {
                        if (curAnnot.getStart() > curVariant.getPosition()) {

                            int startSeq;
                            if (!prevAnnot.getContig().equals(curAnnot.getContig())) {
                                startSeq = 0;
                            } else {
                                startSeq = getFwdDirectionStart(prevAnnot) - 1;
                            }
                            int endSeq = getFwdDirectionEnd(curAnnot);

                            int distanceToPrev;
                            if (!prevAnnot.getContig().equals(curAnnot.getContig())) {
                                distanceToPrev = curVariant.getPosition();
                            } else {
                                distanceToPrev = curVariant.getPosition() - getFwdDirectionEnd(prevAnnot);
                            }
                            int distanceToNext = getFwdDirectionStart(curAnnot) - curVariant.getPosition();

                            String type;
                            if ((distanceToPrev <= REG_REGION && prevAnnot.getStrand() == '-') && (distanceToNext <= REG_REGION && curAnnot.getStrand() == '+')) {
                                type = "REG";
                            } else if (distanceToPrev <= REG_REGION && prevAnnot.getStrand() == '-') {
                                type = "REG";
                            } else if (distanceToNext <= REG_REGION && curAnnot.getStrand() == '+') {
                                type = "REG";
                            } else {
                                type = "INT";
                            }

                            Integer pos = curVariant.getPosition() - startSeq - 1;

                            //System.out.print( "\t" + (startSeq + 1) + "\t" + endSeq);
                            //System.out.print( "\t" + pos );
                            //System.out.println( "\t" + type);

                            if(type.equals("INT")){
                                intCount++;
                            }else{
                                regCount++;
                            }


                            //System.out.print("\t" + distanceToPrev +"\t"+ prevAnnot.getStrand());
                            //System.out.print("\t" + distanceToNext +"\t"+ curAnnot.getStrand());
                            break;
                        }
                    }
                    prevAnnot = curAnnot;
                }
            }
        }
        System.out.println("\tDONE");
    }
    out.println("intergenic regulatory");
    out.println(intCount+" "+regCount);
    out.close();
    }

    /**
     *
     * @param curVar    the variant to be tested
     * @param annotList the list of annotations
     * @return          true is the variant is intergenic
     */
    private static Boolean isIntergenic(Variant curVar, ArrayList<Annotation> annotList){
        for (Annotation curAnnot : annotList) {
            if (isVariantInAnnotation(curVar, curAnnot)) {
                return false;
            }
        }
        return true;
    }

    /**
     *
     * @param var
     * @param annot
     * @return
     */
    private static Boolean isVariantInAnnotation(Variant var, Annotation annot) {
        if (annot.getStrand() == '+') {
            if (annot.getStart() <= var.getPosition() && annot.getStop() >= var.getPosition()) {
                if (annot.getContig().equals(var.getContig())) {
                    return true;
                }
            }
        } else {
            if (annot.getStart() >= var.getPosition() && annot.getStop() <= var.getPosition()) {
                if (annot.getContig().equals(var.getContig())) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     *
     * @param annotFile
     * @return
     * @throws FileNotFoundException
     */
    private static ArrayList<Annotation> getAnnotations(String annotFile)
            throws FileNotFoundException {
        ArrayList<Annotation> result = new ArrayList<Annotation>();
        Scanner in = new Scanner(new File(annotFile));
        in.nextLine();
        while (in.hasNextLine()) {
            Annotation newAnnot = new Annotation();
            String[] temp = in.nextLine().split("\t");
            newAnnot.setContig(temp[0]);
            newAnnot.setFeature(temp[1]);
            newAnnot.setType(temp[2]);
            newAnnot.setLocation(temp[3]);
            newAnnot.setStart(temp[4]);
            newAnnot.setStop(temp[5]);
            newAnnot.setStrand(temp[6]);
            newAnnot.setFunction(temp[7]);
            newAnnot.setAliases(temp[8]);
            newAnnot.setFigfam(temp[9]);
            newAnnot.setEvidenceCodes(temp[10]);
            newAnnot.setNucleotideSeq(temp[11]);
            if (temp.length == 13) {
                newAnnot.setAminoAcidSeq(temp[12]);
            } else {
                newAnnot.setAminoAcidSeq("-");
            }
            result.add(newAnnot);
        }
        return result;
    }

    private static int getFwdDirectionStart(Annotation previous) {
        if (previous.getStart() > previous.getStop()) {
            return previous.getStop();
        } else {
            return previous.getStart();
        }
    }


    private static int getFwdDirectionEnd(Annotation next) {
        if (next.getStart() > next.getStop()) {
            return next.getStart();
        } else {
            return next.getStop();
        }
    }

    /**
     * This class can hold information about annotations.
     */
    static class Annotation {
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

        public Annotation() {
        }

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

        public void setStrand(String strand) {
            this.strand = strand.charAt(0);
        }

        public Integer getStart() {
            return start;
        }

        public void setStart(Integer start) {
            this.start = start;
        }

        public void setStart(String start) {
            this.start = Integer.parseInt(start);
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

        public void setStop(String stop) {
            this.stop = Integer.parseInt(stop);
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

    /**
     * This class can hold information about variants.
     */
    static class Variant {
        private String id;
        private String contig;
        private Integer position;
        private Character reference;
        private Character alternative;

        public Variant() {
        }

        public String getContig() {
            return contig;
        }

        public void setContig(String contig) {
            this.contig = contig;
        }

        public Integer getPosition() {
            return position;
        }

        public void setPosition(Integer position) {
            this.position = position;
        }

        public void setPosition(String position) {
            this.position = Integer.parseInt(position);
        }

        public Character getAlternative() {
            return alternative;
        }

        public void setAlternative(Character alternative) {
            this.alternative = alternative;
        }

        public void setAlternative(String alternative) {
            this.alternative = alternative.charAt(0);
        }

        public Character getReference() {
            return reference;
        }

        public void setReference(Character reference) {
            this.reference = reference;
        }

        public void setReference(String reference) {
            this.reference = reference.charAt(0);
        }

        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
        }
    }

}
