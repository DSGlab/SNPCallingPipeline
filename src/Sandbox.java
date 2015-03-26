import JavaBrew.FastaPrinter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 2/26/15.
 *
 *
 */
public class Sandbox {
    public static void main ( String[] args ) throws FileNotFoundException {
        LinkedHashMap<String, String> seqs = getReads();

        Scanner in = new Scanner( new File( "/Users/juliofdiaz/Desktop/newlistsnp.txt" ) );
        while( in.hasNextLine() ){
            String[] line = in.nextLine().split( "\t" );
            String iso = line[0];
            int pos = Integer.parseInt( line[3] );
            String base = line[5];
            //System.out.println( iso+"\t"+pos+"\t"+base );
            String tempSeq = seqs.get(iso);
            //System.out.println( tempSeq );
            String newtemp = tempSeq.substring( 0, (pos-1) )+base+tempSeq.substring( pos );
            seqs.put(iso, newtemp);
            //System.out.println( seqs.get( iso ) );
        }

        FastaPrinter newOne = new FastaPrinter();
        newOne.setSequences( seqs );
        newOne.simplePrint( "/Users/juliofdiaz/Desktop/CF67-snp-02252015.fasta" );
    }

    private static LinkedHashMap<String, String> getReads() throws FileNotFoundException {
        File f = new File( "/Users/juliofdiaz/Dropbox/CF67_final.fa" );
        FastaPrinter pw = new FastaPrinter( f );
        return pw.getSequences();
    }


/*
//  This Cde allowed me to link snpnames to snp number and isolate id to seq number

    public static void main ( String[] args ) throws FileNotFoundException {
        Scanner in = new Scanner( new File( "/Users/juliofdiaz/Desktop/newlistsnp.txt" ) );

        while ( in.hasNextLine() ) {
            String[] line = in.nextLine().split("\t");
            System.out.println( getId( line[0] )+"\t"+line[0]+"\t"+line[1]+"\t"+getSNPNum( line[1] )+"\t"+line[2] );
        }
    }

    private static String getId ( String seqNum ) throws FileNotFoundException {
        LinkedHashMap<String,String> result = new LinkedHashMap<String, String>();

        String[] id = new String[]{
                "1a", "2a", "3a", "4a", "5a", "6a", "7a", "8a", "9a", "10a", "11a", "12a",
                "1b", "2b", "3b", "4b", "5b", "6b", "7b", "8b", "9b", "10b", "11b", "12b",
                "1c", "2c", "3c", "4c", "5c", "6c", "7c", "8c", "9c", "10c", "11c", "12c",
                "1d", "2d", "3d", "4d", "5d", "6d", "7d", "8d", "9d", "10d", "11d", "12d",
                "1e", "2e", "3e", "4e", "5e", "6e", "7e", "8e", "9e", "10e", "11e", "12e",
                "1f", "2f", "3f", "4f", "5f", "6f", "7f", "8f", "9f", "10f", "11f", "12f",
                "1g", "2g", "3g", "4g", "5g", "6g", "7g", "8g", "9g", "10g", "11g", "12g",
                "1h", "2h", "3h", "4h", "5h", "6h", "7h", "8h", "9h", "10h", "11h", "12h",
                "1i", "2i", "3i", "4i", "5i", "6i", "7i", "8i", "9i", "10i", "11i", "12i",
                "1j", "2j", "3j", "4j", "5j", "6j", "7j", "8j", "9j", "10j", "11j", "12j",
                "1k", "2k", "3k", "4k", "5k", "6k", "7k", "8k", "9k", "10k", "11k", "12k",
                "1l", "2l", "3l", "4l", "5l", "6l", "7l", "8l", "9l", "10l", "11l", "12l",
                "1m", "2m", "3m", "4m", "5m", "6m", "7m", "8m", "9m", "10m", "11m", "12m",
                "1n", "2n", "3n", "4n", "5n", "6n", "7n", "8n", "9n", "10n", "11n", "12n",
                "1o", "2o", "3o", "4o", "5o", "6o", "7o", "8o", "9o", "10o", "11o", "12o",
                "1p", "2p", "3p", "4p", "5p", "6p", "7p", "8p", "9p", "10p", "11p", "12p",
                "1q", "2q", "3q", "4q", "5q", "6q", "7q", "8q", "9q", "10q", "11q", "12q",
                "1r", "2r", "3r", "4r", "5r", "6r", "7r", "8r", "9r", "10r", "11r", "12r",
                "1s", "2s", "3s", "4s", "5s", "6s", "7s", "8s", "9s", "10s", "11s",
                "1t", "2t", "3t", "4t", "5t", "6t", "7t", "8t", "9t", "10t", "11t" };

        String[] num = new String[]{
                "51", "29", "89", "109AB", "129", "1", "149", "21", "169", "26B", "31", "33",
                "52", "30", "90AB", "110", "130ACD", "2", "150AB", "22", "170", "27B", "32", "34",
                "53", "71", "91", "111", "131", "3", "151", "23B", "171", "28B", "35", "223",
                "54", "72", "92AB", "112AB", "132", "4", "152", "24", "172", "189", "206", "224",
                "55", "73", "93", "113AB", "133", "5", "153", "25B", "173", "190", "207", "225AB",
                "56", "74", "94", "114", "134", "6", "154", "36", "174", "191", "208AB", "226",
                "57", "75", "95", "115", "135", "7", "155", "37", "175AB", "192AB", "209", "227",
                "58", "76", "96", "116", "136", "8", "156", "38", "176", "193", "210", "228AB",
                "59", "77", "97", "117", "137", "9", "157AB", "39", "177", "194", "211", "229",
                "60", "78", "98AB", "118", "138", "10", "158AB", "40", "178AB", "195AB", "212AB", "230AB",
                "61B", "79", "99", "119", "139", "11", "159", "41", "179", "196", "213", "231",
                "62B", "80", "100AB", "120AB", "140", "12", "160AB", "42", "180AB", "197", "214", "232",
                "63", "81", "101", "121", "141", "13", "161", "43", "181", "198ACD", "215AB", "233",
                "64AB", "82", "102AB", "122AB", "142AB", "14", "162AB", "44", "182AB", "199", "216", "234",
                "65AB", "83", "103", "123", "143", "15", "163", "45", "183", "200", "217", "235AB",
                "66", "84", "104", "124", "144AB", "16", "164", "46", "184", "201", "218AB", "236",
                "67", "85", "105", "125", "145AB", "17", "165AB", "47", "185AB", "202AB", "219", "237",
                "68", "86", "106", "126", "146", "18", "166", "48", "186", "203", "220", "238",
                "69", "87", "107", "127", "147AB", "19", "167", "49", "187", "204", "221",
                "70", "88", "108", "128", "148AB", "20", "168AB", "50", "188AB", "205", "222AB" };

        for ( int i=0; i<238; i++) {
            result.put(num[i], id[i]);
        }

        return result.get( seqNum );
    }

    private static String getSNPNum ( String snpName ) throws FileNotFoundException {
        Scanner in = new Scanner( new File( "/Users/juliofdiaz/Documents/CF67/snp_calling/CF67_C71/snplist2.txt" ) );
        Hashtable<String, String> table = new Hashtable<String, String>();

        int count=1;
        while ( in.hasNextLine() ) {
            //System.out.println( in.nextLine() );
            table.put(in.nextLine().trim(), count + "");

            count++;
        }
        //System.out.println( snpName );
        return table.get( snpName );
    }*/

}
