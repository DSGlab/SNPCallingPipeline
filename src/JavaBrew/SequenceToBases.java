package JavaBrew;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedHashMap;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 5/9/16.
 */
public class SequenceToBases {

    public static void main(String[]main) throws FileNotFoundException {
        String INPUT = "/Users/juliofdiaz/Desktop/tbd.txt";

        File inputFile = new File(INPUT);
        //LinkedHashMap<String,String> seqs = new FastaPrinter(inputFile).getSequences();
        LinkedHashMap<String,String> seqs = getSeqs(inputFile);

        for(String key: seqs.keySet()){
            System.out.print(key+"\t");
            for(int i=0; i<seqs.get(key).length(); i++){
                System.out.print(seqs.get(key).charAt(i)+"\t");
            }
            System.out.println();
        }
    }

    private static LinkedHashMap<String,String> getSeqs(File input) throws FileNotFoundException {
        LinkedHashMap<String,String> result = new LinkedHashMap<String, String>();

        Scanner in = new Scanner(input);
        while(in.hasNextLine()){

            result.put( in.nextLine(), in.nextLine());
        }
        return result;
    }
}
