package JavaBrew;

import org.jetbrains.annotations.Contract;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 4/28/16.
 *
 * Class only filled with useful methods and variables
 */
public class Utilities {
    public static final String CONF_FILE = "/Users/juliofdiaz/Dropbox/CF/snp_calling/TEST/conf.txt";
    public static final String REPLICON_POS_SEPARATOR = "-";
    public static final String ID_LABEL_SEPARATOR = "_";
    public static final String DIR_SEPARATOR = "/";
    public static final String VCF_NAME = "r.vcf";
    public static final String SEPARATOR = "\t";

    public static final String FASTA_START = ">";
    public static final String REFERENCE_DEFAULT_NAME = "reference";

    /**
     *
     */
    public static final String HOME_PATH = System.getProperty("user.home");

    /**
     * This method takes a base ('A', 'T', 'G' or 'C') and it returns the
     * complement base.
     *
     * @param base  The base to be complemented. Only valid for inputs 'A',
     *              'T', 'G', 'C' or 'N'
     * @return      The complemented base of the input base. Input 'A' would
     *              return 'T'. Input 'T' would return 'A'. Input 'G' would
     *              return 'C'. Input 'C' would return 'G'. Any other input
     *              would return 'N'
     */
    @Contract(pure = true)
    public static char getComplementaryBase (char base ) {
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

    /**
     * This method takes a codon and it returns its respective amino acid
     * one-letter code.
     *
     * @param inputCodon    A String of three letters representing a codon.
     *                      The only bases allowed in the codons are 'A',
     *                      'T', 'G' and 'C'.
     * @return              The one-letter representation of the amino acid
     *                      respective to the input codon. '*' represents
     *                      stop codons.
     */
    public static char getAA ( String inputCodon ) {
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

    /**
     * This method return a range of numbers as a list of Strings. This could be useful
     * when looping through isolates with a numerical ID.
     *
     * @param start The beginning of the range.
     * @param end   The end of the range.
     * @return      The range of numbers as a list of Strings.
     */
    public static ArrayList<String> getRange(int start, int end){
        ArrayList<String> result = new ArrayList<String>();

        for(int i=start; i<=end; i++){
            result.add(i+"");
        }
        return result;
    }

    /**
     *
     * @return
     * @throws FileNotFoundException
     */
    public static Hashtable<String, String> getConfigurationTable(String optionTable)
            throws FileNotFoundException {
        Scanner in = new Scanner(new File(optionTable));
        Hashtable<String,String> result = new Hashtable<String, String>();
        while(in.hasNextLine()){
            String line = in.nextLine();
            if(!line.isEmpty()) {
                if (line.charAt(0) != '#') {
                    String[] line2 = line.split("=");
                    result.put(line2[0].trim(), line2[1].trim());
                }
            }
        }
        return result;
    }

    /**
     *
     * @param input
     * @return
     */
    public static String[] getIds(String input) throws FileNotFoundException {
        Scanner in = new Scanner(new File(input));
        ArrayList<String> preList= new ArrayList<String>();
        while(in.hasNextLine()){
            preList.add(in.nextLine().trim());
        }

        String [] result = preList.toArray(new String[preList.size()]);
        return result;
    }

}
