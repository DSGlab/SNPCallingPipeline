import JavaBrew.FastaPrinter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
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

    public static void main (String[] args) throws FileNotFoundException {
        String ANNOT_FILE = "/Users/juliofdiaz/Dropbox/CF/references/B6CQ.txt";
        String VAR_FILE = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_NEW/variants.txt";
        String REF_CONTIGS = "/Users/juliofdiaz/Dropbox/CF/references/B6CQ.fa";
        String OUT_FILE = "/Users/juliofdiaz/Dropbox/CF/snp_calling/CF170_NEW/SNP_annot.txt";

        Integer REG_REGION = 150;

        ArrayList<Annotation> annotList = getAnnotations(ANNOT_FILE);
        ArrayList<Variant> varList = getVariants(VAR_FILE);
        LinkedHashMap<String, String> reference = new FastaPrinter(new File(REF_CONTIGS)).getSequences();

        PrintWriter out = new PrintWriter( OUT_FILE );
        out.println("CONTIG\tPOSITION\tSTRAND\tSTART_ANNOT\tEND_ANNOT\t\tREF\tALT\tREF (in annot)\tCODON_POS (0 index)\t" +
                "POS_IN_GENE (o index)\tREF_CODON\tALT_CODON\tREF_AA\tALT_AA\tTYPE\tDNA");
        for( Variant curVar : varList ){
            Boolean isIntragenic = false;
            out.print( curVar.getContig() + "\t" + curVar.getPosition() );
            for( Annotation curAnnot : annotList ){
                if ( isVariantInAnnotation( curVar, curAnnot ) ) {
                    isIntragenic = true;
                    if (curAnnot.getStrand() == '+') {
                        out.print("\t"+curAnnot.getStrand() + "\t" + curAnnot.getStart() + "\t" + curAnnot.getStop() + "\t");
                        out.print( "\t" + curVar.getReference() + "\t" + curVar.getAlternative() + "\t" );
                        int posInGene = curVar.getPosition()-curAnnot.getStart();
                        int posInCodon = posInGene%3;
                        int codonStart = posInGene-posInCodon;
                        String codon = curAnnot.getNucleotideSeq().substring(codonStart,codonStart+3).toUpperCase();
                        String altCodon = getAlternativeCodon(codon, posInCodon, curVar.getAlternative(), curAnnot.getStrand()).toUpperCase();
                        char refAA = getAA(codon.toUpperCase());
                        char altAA = getAA(altCodon.toUpperCase());
                        String mutType = getMutationType(refAA,altAA);
                        out.print( curAnnot.getNucleotideSeq().charAt( posInGene ) + "\t"+posInCodon+"\t"+codonStart+"\t"+codon+"\t"+altCodon+"\t"+refAA+"\t"+altAA+"\t"+mutType+"\t"+curAnnot.getNucleotideSeq() );

                    } else {
                        out.print("\t"+curAnnot.getStrand() + "\t" + curAnnot.getStart() + "\t" + curAnnot.getStop() + "\t");
                        out.print( "\t" + curVar.getReference() + "\t" + curVar.getAlternative() + "\t");
                        int posInGene = curAnnot.getStart()-curVar.getPosition();
                        int posInCodon = posInGene%3;
                        int codonStart = posInGene-posInCodon;
                        String codon = curAnnot.getNucleotideSeq().substring(codonStart,codonStart+3).toUpperCase();
                        String altCodon = getAlternativeCodon(codon, posInCodon, curVar.getAlternative(), curAnnot.getStrand()).toUpperCase();
                        char refAA = getAA(codon);
                        char altAA = getAA(altCodon);
                        String mutType = getMutationType(refAA,altAA);
                        out.print( curAnnot.getNucleotideSeq().charAt( posInGene ) + "\t"+posInCodon+"\t"+codonStart+"\t"+codon+"\t"+altCodon+"\t"+refAA+"\t"+altAA+"\t"+mutType+"\t"+curAnnot.getNucleotideSeq() );
                    }
                }
            }

            if(!isIntragenic){
                Annotation prevAnnot = null;
                for(Annotation curAnnot: annotList) {
                    if ( curAnnot.getContig().equals( curVar.getContig() ) ) {
                        if (curAnnot.getStart() > curVar.getPosition()) {
                            out.print("\t"+"na");

                            int startSeq;
                            if( !prevAnnot.getContig().equals(curAnnot.getContig()) ){
                                startSeq = 0;
                            }else{
                                startSeq = getFwdDirectionStart(prevAnnot)-1;
                            }
                            int endSeq = getFwdDirectionEnd(curAnnot);

                            int distanceToPrev;
                            if(!prevAnnot.getContig().equals(curAnnot.getContig())) {
                                distanceToPrev = curVar.getPosition();
                            }else{
                                distanceToPrev = curVar.getPosition() - getFwdDirectionEnd(prevAnnot);
                            }
                            int distanceToNext = getFwdDirectionStart(curAnnot)-curVar.getPosition();

                            String type;
                            if(( distanceToPrev<=REG_REGION && prevAnnot.getStrand()=='-' )&&( distanceToNext<=REG_REGION && curAnnot.getStrand()=='+' )){
                                type = "RGB";
                            }else if( distanceToPrev<=REG_REGION && prevAnnot.getStrand()=='-' ){
                                type = "RGP";
                            }else if( distanceToNext<=REG_REGION && curAnnot.getStrand()=='+' ) {
                                type = "RGN";
                            }else{
                                type = "INT";
                            }

                            String nucleotideSeq = reference.get( curVar.getContig()).substring( startSeq, endSeq);
                            Integer pos = curVar.getPosition()-startSeq-1;

                            out.print("\t" + (startSeq+1) +"\t"+endSeq+"\t");
                            out.print("\t" + curVar.getReference()+"\t"+curVar.getAlternative()+"\t"+"na"+"\t"+"na");
                            out.print("\t" + pos +"\t"+"na"+"\t"+"na"+"\t"+"na"+"\t"+"na"+"\t"+type);
                            out.print("\t" + nucleotideSeq );
                            out.print("\t" + distanceToPrev +"\t"+ prevAnnot.getStrand());
                            out.print("\t" + distanceToNext +"\t"+ curAnnot.getStrand());
                            break;
                        }
                    }
                    prevAnnot = curAnnot;
                }
            }
            out.println();
        }

        out.close();


    }

    public static int getFwdDirectionStart(Annotation previous){
        if ( previous.getStart()>previous.getStop() ) {
            return previous.getStop();
        }else{
            return previous.getStart();
        }
    }

    public static int getFwdDirectionEnd(Annotation next){
        if ( next.getStart()>next.getStop() ) {
            return next.getStart();
        }else{
            return next.getStop();
        }
    }

    public static String getMutationType( char refAA, char altAA ){
        if( refAA == altAA ){
            return "SYN";
        }else{
            return "NSY";
        }
    }

    public static String getAlternativeCodon(String codon, int posInCodon, char alternative, char strand ){
        if(strand=='+'){
            return codon.substring(0,posInCodon)+""+alternative+""+codon.substring(posInCodon+1,3);

        }else{
            return codon.substring(0,posInCodon)+""+getComplementaryBase( alternative )+""+codon.substring(posInCodon+1,3);
        }
    }

    public static Boolean isVariantInAnnotation( Variant var, Annotation annot ){
        if( annot.getStrand() == '+' ){
            if ( annot.getStart()<= var.getPosition() && annot.getStop()>=var.getPosition() ) {
                if ( annot.getContig().equals( var.getContig() ) ) {
                    return true;
                }
            }
        }else{
            if ( annot.getStart()>= var.getPosition() && annot.getStop()<=var.getPosition() ) {
                if ( annot.getContig().equals( var.getContig() ) ) {
                    return true;
                }
            }
        }
        return false;
    }

    public static ArrayList<Variant> getVariants(String varFile) throws FileNotFoundException {
        ArrayList<Variant> result = new ArrayList<Variant>();
        Scanner in = new Scanner(new File(varFile));
        while(in.hasNextLine()){
            Variant newVariant = new Variant();
            String[] temp = in.nextLine().split("\t");
            newVariant.setId( temp[0] );
            newVariant.setContig( temp[1] );
            newVariant.setPosition( temp[2] );
            newVariant.setReference( temp[3] );
            newVariant.setAlternative( temp[4] );
            result.add(newVariant);
        }
        return result;
    }

    private static class Variant{
        private String id;
        private String contig;
        private Integer position;
        private Character reference;
        private Character alternative;
        public Variant (){ }

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

    public static ArrayList<Annotation> getAnnotations(String annotFile) throws FileNotFoundException {
        ArrayList<Annotation> result = new ArrayList<Annotation>();
        Scanner in = new Scanner(new File(annotFile));
        in.nextLine();
        while( in.hasNextLine() ){
            Annotation newAnnot = new Annotation();
            String [] temp = in.nextLine().split("\t");
            newAnnot.setContig( temp[0] );
            newAnnot.setFeature( temp[1] );
            newAnnot.setType( temp[2] );
            newAnnot.setLocation( temp[3] );
            newAnnot.setStart( temp[4] );
            newAnnot.setStop( temp[5] );
            newAnnot.setStrand( temp[6] );
            newAnnot.setFunction( temp[7] );
            newAnnot.setAliases( temp[8] );
            newAnnot.setFigfam( temp[9] );
            newAnnot.setEvidenceCodes( temp[10] );
            newAnnot.setNucleotideSeq( temp[11] );
            if(temp.length==13) {
                newAnnot.setAminoAcidSeq(temp[12]);
            }else{
                newAnnot.setAminoAcidSeq("-");
            }
            result.add(newAnnot);
        }
        return result;
    }

    private static class Annotation{
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

    private static char getAA ( String inputCodon ) {
        Hashtable<String, Character> table = new Hashtable<String, Character>();
        table.put ("TTT", 'F');
        table.put ("TTC", 'F');
        table.put ("TTA", 'L');
        table.put ("TTG", 'L');
        table.put ("TCT", 'S');
        table.put ("TCC", 'S');
        table.put ("TCA", 'S');
        table.put ("TCG", 'S');
        table.put ("TAT", 'Y');
        table.put ("TAC", 'Y');
        table.put ("TAA", '*');
        table.put ("TAG", '*');
        table.put ("TGT", 'C');
        table.put ("TGC", 'C');
        table.put ("TGA", '*');
        table.put ("TGG", 'W');
        table.put ("CTT", 'L');
        table.put ("CTC", 'L');
        table.put ("CTA", 'L');
        table.put ("CTG", 'L');
        table.put ("CCT", 'P');
        table.put ("CCC", 'P');
        table.put ("CCA", 'P');
        table.put ("CCG", 'P');
        table.put ("CAT", 'H');
        table.put ("CAC", 'H');
        table.put ("CAA", 'Q');
        table.put ("CAG", 'Q');
        table.put ("CGT", 'R');
        table.put ("CGC", 'R');
        table.put ("CGA", 'R');
        table.put ("CGG", 'R');
        table.put ("ATT", 'I');
        table.put ("ATC", 'I');
        table.put ("ATA", 'I');
        table.put ("ATG", 'M');
        table.put ("ACT", 'T');
        table.put ("ACC", 'T');
        table.put ("ACA", 'T');
        table.put ("ACG", 'T');
        table.put ("AAT", 'N');
        table.put ("AAC", 'N');
        table.put ("AAA", 'K');
        table.put ("AAG", 'K');
        table.put ("AGT", 'S');
        table.put ("AGC", 'S');
        table.put ("AGA", 'R');
        table.put ("AGG", 'R');
        table.put ("GTT", 'V');
        table.put ("GTC", 'V');
        table.put ("GTA", 'V');
        table.put ("GTG", 'V');
        table.put ("GCT", 'A');
        table.put ("GCC", 'A');
        table.put ("GCA", 'A');
        table.put ("GCG", 'A');
        table.put ("GAT", 'D');
        table.put ("GAC", 'D');
        table.put ("GAA", 'E');
        table.put ("GAG", 'E');
        table.put ("GGT", 'G');
        table.put ("GGC", 'G');
        table.put ("GGA", 'G');
        table.put ("GGG", 'G');
        return table.get(inputCodon);
    }

    private static char getComplementaryBase ( char base ) {
        if ( base=='A' || base=='a' ) {
            return 'T';
        } else if ( base=='T' || base=='t' ) {
            return 'A';
        } else if ( base=='C' || base=='c' ) {
            return 'G';
        } else if ( base=='G' || base=='g' ) {
            return 'C';
        } else {
            return 'N';
        }
    }
}
