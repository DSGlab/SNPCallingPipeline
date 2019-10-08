import org.biojava.bio.BioException;
import org.biojava.bio.seq.Feature;
import org.biojavax.Note;
import org.biojavax.RichAnnotation;
import org.biojavax.RichObjectFactory;
import org.biojavax.bio.seq.RichFeature;
import org.biojavax.bio.seq.RichLocation;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;

//import org.biojava.utils.bytecode.CodeException;
import org.biojavax.ontology.ComparableTerm;



import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;

public class GbkReader {

    public static void main(String[] args){
        String file = "/Users/juliofdiaz/Documents/results/prokka/p190523-dr-tm-sp-pi/EP8/EP8.gbk";
        ArrayList<ExtractSNPInfoFromRast.Annotation> annot = process(file);
        System.out.println(annot);
    }

    public static ArrayList<ExtractSNPInfoFromRast.Annotation> process(String annotFile) {
        String file = annotFile;

        ComparableTerm locusTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("locus_tag");
        ComparableTerm inferenceTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("inference");
        ComparableTerm productTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("product");
        ComparableTerm aaTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("translation");
        ComparableTerm geneTerm = RichObjectFactory.getDefaultOntology().getOrCreateTerm("gene");

        RichSequenceIterator richSeq;

        ArrayList<ExtractSNPInfoFromRast.Annotation> result = new ArrayList<ExtractSNPInfoFromRast.Annotation>();
        try {
            richSeq = RichSequence.IOTools.readGenbankDNA(new BufferedReader(new FileReader(file)), null);


            while (richSeq.hasNext()) {
                RichSequence temp_rs = richSeq.nextRichSequence();
                String contig = temp_rs.getName();
                Iterator<Feature> fh = temp_rs.features();
                while(fh.hasNext()){
                    RichFeature feat = (RichFeature) fh.next();
                    RichLocation loc = (RichLocation) feat.getLocation();

                    char strand = feat.getStrand().getToken();
                    int start = (strand=='+')?loc.getMin():loc.getMax();
                    int stop = (strand=='+')?loc.getMax():loc.getMin();

                    String seq = feat.getSequence().seqString().substring(loc.getMin()-1,loc.getMax());
                    String type = feat.getType();

                    RichAnnotation annot = (RichAnnotation) feat.getAnnotation();

                    String feature_id = "";
                    String location = contig + "_" + start + "_" + stop;
                    String inference = "";
                    String product = "";
                    String aa = "-";
                    String gene = "";
                    String figfam = "";

                    Set<Note> set= annot.getNoteSet();
                    for(Note n : set){

                        if(n.getTerm().equals(locusTerm)){ feature_id = n.getValue(); }
                        if(n.getTerm().equals(inferenceTerm)){ inference = n.getValue(); }
                        if(n.getTerm().equals(productTerm)){ product = n.getValue(); }
                        if(n.getTerm().equals(aaTerm)){ aa = n.getValue(); }
                        if(n.getTerm().equals(geneTerm)){ gene = n.getValue(); }
                    }
                    if(!type.equals("source")) {
                        ExtractSNPInfoFromRast.Annotation newAnnot = new ExtractSNPInfoFromRast.Annotation();

                        newAnnot.setContig( contig );
                        newAnnot.setFeature( feature_id );
                        newAnnot.setType( type );
                        newAnnot.setLocation( location );
                        newAnnot.setStart( start );
                        newAnnot.setStop( stop );
                        newAnnot.setStrand( strand );
                        newAnnot.setFunction( product );
                        newAnnot.setAliases( gene );
                        newAnnot.setFigfam( figfam );
                        newAnnot.setEvidenceCodes( inference );
                        newAnnot.setNucleotideSeq( seq );
                        newAnnot.setAminoAcidSeq( aa );

                        result.add(newAnnot);

                        //System.out.println(contig + "\t" + feature_id + "\t" + type + "\t" + location + "\t" + start + "\t" + stop + "\t" + strand + "\t" + product + "\t" + gene + "\t" + figfam + "\t" + inference + "\t" + seq + "\t" + aa);
                    }
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (BioException e) {
            e.printStackTrace();
        }
        return result;
    }
}

