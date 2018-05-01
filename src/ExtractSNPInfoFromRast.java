import JavaBrew.FastaPrinter;
import JavaBrew.Utilities;
import org.jetbrains.annotations.Contract;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Scanner;

/**
 * This program returns information about SNPs given an
 * annotated reference. The annotated reference should be in
 * "Spreadsheet (tab-separated text format)" style from RAST. The
 * SNPs should be in the format specified in doc/Sample-VAR_FILE.txt
 * In order to get information about the SNPs found in intragenic
 * regions, include the reference in fasta format.
 *
 * Type of SNPS:    SYN synonymous
 *                  NSY non-synonymous
 *                  INT intergenic
 *                  RGN Regulatory region of next gene
 *                  RGP Regulatory region of previous gene
 *                  RGB Regulatory region of both adjacent genes
 *
 * Regulatory regions are specified in REG_REGION.
 *
 * @author juliofdiazc
 */
public class ExtractSNPInfoFromRast {

    private static final Integer REG_REGION = 150;

    public static void main(String[] args) throws FileNotFoundException {
        /*
         *
         */
        String ANNOT_FILE = Utilities.HOME_PATH+"/Dropbox/CF/references/3j.txt";

        /*
         * The file containing information about the variants
         */
        String VAR_FILE = Utilities.HOME_PATH+"/Dropbox/CF/snp_calling/DOLOSA_3j/variants.txt";

        /*
         * The reference assembly contigs in a fasta format
         */
        String REF_CONTIGS = Utilities.HOME_PATH+"/Dropbox/CF/references/3j.fa";

        /*
         * The name of the output file containing the information of the
         * annotations with variants.
         */
        String OUT_FILE = Utilities.HOME_PATH+"/Dropbox/CF/snp_annotation/DOLOSA_3j/SNP_annot.txt";

        /*
         * Here, the files are processed and the output is created
         */
        process(ANNOT_FILE,VAR_FILE,REF_CONTIGS,OUT_FILE);
    }

    /**
     * This method takes the main variables of this program and runs the main
     * method work().
     *
     * @param annotFile     The location of the annotation file
     * @param varFile       The location of the variant file
     * @param refContigs    The location of the contings in a fasta format
     * @param outFile       The location where the output file will be saved.
     * @throws FileNotFoundException
     */
    public static void process(String annotFile, String varFile, String refContigs, String outFile)
            throws FileNotFoundException {
        ArrayList<Annotation> annotList = getAnnotations( annotFile );
        ArrayList<Variant> varList = getVariants( varFile );
        LinkedHashMap<String, String> reference = new FastaPrinter(new File( refContigs )).getSequences();

        PrintWriter out = new PrintWriter( outFile );
        out.println("CONTIG\tPOSITION\tSTRAND\tSTART ANNOT\tEND ANNOT\tREF\tALT\tREF (in annot)\tPOS IN CODON (0 index)\t" +
                "POS IN ANNOT (0 index)\tREF CODON\tALT CODON\tREF AA\tALT AA\tTYPE\tDNA");

        for (Variant curVar : varList) {
            work(curVar, annotList, reference, out);
        }

        out.close();
    }


    /**
     * This method does most of the work. Here we find the annotation related to each
     * variant. And it als retrieves the nucleotide sequences for each annotation. Indels
     * are also checked to see if they are found in a regulatory region. Based on the
     *
     * @param curVar    variant information
     * @param annotList annotation information
     * @param reference reference contigs in fasta file
     * @param out       out file
     * @return          The modified out file
     */
    private static PrintWriter work(Variant curVar, ArrayList<Annotation> annotList,
                                    LinkedHashMap<String, String> reference, PrintWriter out) {

        Boolean isIntragenic = false;
        out.print(curVar.getContig() + "\t" + curVar.getPosition());
        for (Annotation curAnnot : annotList) {
            if (isVariantInAnnotation(curVar, curAnnot)) {
                isIntragenic = true;
                if (curAnnot.getStrand() == '+') {
                    out.print("\t" + curAnnot.getStrand() + "\t" + curAnnot.getStart() + "\t" + curAnnot.getStop());
                    out.print("\t" + curVar.getReference() + "\t" + curVar.getAlternative() + "\t");
                    int posInGene = curVar.getPosition() - curAnnot.getStart();
                    int posInCodon = posInGene % 3;
                    int codonStart = posInGene - posInCodon;
                    String codon = curAnnot.getNucleotideSeq().substring(codonStart, codonStart + 3).toUpperCase();
                    String altCodon = getAlternativeCodon(codon, posInCodon, curVar.getAlternative(), curAnnot.getStrand()).toUpperCase();
                    char refAA = Utilities.getAA(codon.toUpperCase());
                    char altAA = Utilities.getAA(altCodon.toUpperCase());
                    String mutType = getMutationType(refAA, altAA);
                    out.print(curAnnot.getNucleotideSeq().charAt(posInGene) + "\t" + posInCodon + "\t" + posInGene + "\t" + codon + "\t" + altCodon + "\t" + refAA + "\t" + altAA + "\t" + mutType + "\t" + curAnnot.getFeature() + "\t" + curAnnot.getFunction() + "\t" + curAnnot.getNucleotideSeq());

                } else {
                    out.print("\t" + curAnnot.getStrand() + "\t" + curAnnot.getStart() + "\t" + curAnnot.getStop());
                    out.print("\t" + curVar.getReference() + "\t" + curVar.getAlternative() + "\t");
                    int posInGene = curAnnot.getStart() - curVar.getPosition();
                    int posInCodon = posInGene % 3;
                    int codonStart = posInGene - posInCodon;
                    String codon = curAnnot.getNucleotideSeq().substring(codonStart, codonStart + 3).toUpperCase();
                    String altCodon = getAlternativeCodon(codon, posInCodon, curVar.getAlternative(), curAnnot.getStrand()).toUpperCase();
                    char refAA = Utilities.getAA(codon);
                    char altAA = Utilities.getAA(altCodon);
                    String mutType = getMutationType(refAA, altAA);
                    out.print(curAnnot.getNucleotideSeq().charAt(posInGene) + "\t" + posInCodon + "\t" + posInGene + "\t" + codon + "\t" + altCodon + "\t" + refAA + "\t" + altAA + "\t" + mutType + "\t" + curAnnot.getFeature() + "\t" + curAnnot.getFunction() + "\t" + curAnnot.getNucleotideSeq());
                }
            }
        }

        if (!isIntragenic) {
            Annotation prevAnnot = null;
            for (Annotation curAnnot : annotList) {
                if (curAnnot.getContig().equals(curVar.getContig())) {
                    if (curAnnot.getStart() > curVar.getPosition()) {
                        out.print("\t" + "na");

                        int startSeq;
                        if (!prevAnnot.getContig().equals(curAnnot.getContig())) {
                            startSeq = 0;
                        } else {
                            startSeq = getFwdDirectionStart(prevAnnot) - 1;
                        }
                        int endSeq = getFwdDirectionEnd(curAnnot);

                        int distanceToPrev;
                        if (!prevAnnot.getContig().equals(curAnnot.getContig())) {
                            distanceToPrev = curVar.getPosition();
                        } else {
                            distanceToPrev = curVar.getPosition() - getFwdDirectionEnd(prevAnnot);
                        }
                        int distanceToNext = getFwdDirectionStart(curAnnot) - curVar.getPosition();

                        String type;
                        if ((distanceToPrev <= REG_REGION && prevAnnot.getStrand() == '-') && (distanceToNext <= REG_REGION && curAnnot.getStrand() == '+')) {
                            type = "RGB";
                        } else if (distanceToPrev <= REG_REGION && prevAnnot.getStrand() == '-') {
                            type = "RGP";
                        } else if (distanceToNext <= REG_REGION && curAnnot.getStrand() == '+') {
                            type = "RGN";
                        } else {
                            type = "INT";
                        }

                        String nucleotideSeq = reference.get(curVar.getContig()).substring(startSeq, endSeq);
                        Integer pos = curVar.getPosition() - startSeq - 1;

                        out.print("\t" + (startSeq + 1) + "\t" + endSeq);
                        out.print("\t" + curVar.getReference() + "\t" + curVar.getAlternative() + "\t" + "na" + "\t" + "na");
                        out.print("\t" + pos + "\t" + "na" + "\t" + "na" + "\t" + "na" + "\t" + "na" + "\t" + type);
                        out.print("\t" + nucleotideSeq);

                        if (type.equals("RGP") || type.equals("RGB")) {
                            out.print("\t" + prevAnnot.getFeature() + "\t" + prevAnnot.getFunction());
                        }
                        if (type.equals("RGN") || type.equals("RGB")) {
                            out.print("\t" + curAnnot.getFeature() + "\t" + curAnnot.getFunction());
                        }

                        //System.out.print("\t" + distanceToPrev +"\t"+ prevAnnot.getStrand());
                        //System.out.print("\t" + distanceToNext +"\t"+ curAnnot.getStrand());
                        break;
                    }
                }
                prevAnnot = curAnnot;
            }
        }
        out.println();


        return out;
    }

    /**
     * This method evaluates an annotation and it returns the start position
     * of the annotation as if the annotation was in the '+' strand.
     *
     * @param previous  The annotation to be assessed.
     * @return          The start position of the annotation in its contig.
     */
    private static int getFwdDirectionStart(Annotation previous) {
        if (previous.getStart() > previous.getStop()) {
            return previous.getStop();
        } else {
            return previous.getStart();
        }
    }

    /**
     *  This method evaluates an annotation and it returns the end position
     *  of the annotation as if the annotation was in the '+' strand.
     *
     * @param next  The annotation to be assessed.
     * @return      The end position of the annotation in its contig.
     */
    private static int getFwdDirectionEnd(Annotation next) {
        if (next.getStart() > next.getStop()) {
            return next.getStart();
        } else {
            return next.getStop();
        }
    }

    /**
     * This method evaluates mutations in coding sequences
     *
     * @param refAA One-letter representation of the reference amino acid.
     * @param altAA One-letter representation of the alternate amino acid.
     * @return      "SYN" for synonymous mutations and "NSY" for non synonymous
     *              mutations.
     */
    @NotNull
    @Contract(pure = true)
    private static String getMutationType(char refAA, char altAA) {
        if (refAA == altAA) {
            return "SYN";
        } else {
            return "NSY";
        }
    }

    /**
     * This method allows the alteration of a codon given a position that
     * needs to be modified and the base that needs to be modified for.
     *
     * @param codon         A String of three letters representing a codon.
     *                      The only bases allowed in the codons are 'A',
     *                      'T', 'G' and 'C'.
     * @param posInCodon    The position in the codon that will be altered
     * @param alternative   The base to be inserted instead of the one in the
     *                      posInCodon.
     * @param strand        The strand of the annotation.
     * @return              The altered codon with the new alternative base at
     *                      the posInCodon in the codon.
     */
    private static String getAlternativeCodon(String codon, int posInCodon,
                                              char alternative, char strand) {
        if (strand == '+') {
            return codon.substring(0, posInCodon) + "" +
                    alternative + "" +
                    codon.substring(posInCodon + 1, 3);

        } else {
            return codon.substring(0, posInCodon) + "" +
                    Utilities.getComplementaryBase(alternative) + "" +
                    codon.substring(posInCodon + 1, 3);
        }
    }

    /**
     * This method tests if a variant can be found in a given annotation.
     *
     * @param var   The variant to be tested.
     * @param annot The annotation where the variant would be searched for.
     * @return      True if the annotation encompasses the variant. False if
     *              the annotation does not emcompass the variant.
     */
    @NotNull
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
     * This method retrieves the variants in a variant file created by the user
     * (see notebook).
     *
     * @param varFile The location of the file containing the information of
     *                the variants to be analyzed.
     * @return        A list of the variants found in the varFile.
     * @throws FileNotFoundException
     */
    private static ArrayList<Variant> getVariants(String varFile)
            throws FileNotFoundException {
        ArrayList<Variant> result = new ArrayList<Variant>();
        Scanner in = new Scanner(new File(varFile));
        while (in.hasNextLine()) {
            Variant newVariant = new Variant();
            String[] temp = in.nextLine().split("\t");
            newVariant.setId(temp[0]);
            newVariant.setContig(temp[1]);
            newVariant.setPosition(temp[2]);
            newVariant.setReference(temp[3]);
            newVariant.setAlternative(temp[4]);
            result.add(newVariant);
        }
        return result;
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

    /**
     * This method evaluates the tab-separated values file containing the annotation
     * information obtained from RAST.
     *
     * @param annotFile The location of the tab-separated values file downloaded from
     *                  RAST.
     * @return          A list of the annotations found in the annotFile.
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
}