package IndelAnalysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Hashtable;
import java.util.Scanner;

public class IndelInfoExtracter {
    public static void main(String[] args) throws FileNotFoundException {

        for(String i : new String[]{"1","2","3","4","5","6","7","8"}) {
            for (String j : new String[]{"a","b","c","d","e","f","g","h","i","j"}) {
                String inFile = "/Users/juliofdiaz/Dropbox/CF/indel_calling/DOLOSA_2e/"+i+j+".ind.vcf";
                run(inFile);
            }
        }
    }

    public static void run(String inFile) throws FileNotFoundException {
        Scanner in = new Scanner(new File(inFile));
        while(in.hasNextLine()){
            String line = in.nextLine();
            if(line.charAt(0)!='#'){
                String[] info = line.split("\t");
                System.out.print(inFile+"\t");
                System.out.print(info[0]+"\t");
                System.out.print(info[1]+"\t");
                //System.out.println(info[2]);
                System.out.print(info[3]+"\t");
                System.out.print(info[4]+"\t");
                System.out.print(info[5]+"\t");
                //System.out.println(info[6]);
                //System.out.println(info[7]);
                System.out.print(getDP4(info[7])+"\t");
                System.out.println(isSTR(info[7]));
                //System.out.println(info[9]);
            }
        }
    }

    private static String isSTR(String cigar){
        String[] info = cigar.split(";");
        for(String item: info){
            if(item.equals("STR")){
                return "STR";
            }
        }
        return "";
    }

    private static String getDP4(String cigar){
        Hashtable<String, String> table = new Hashtable<String, String>();
        String[] info = cigar.split(";");
        for(String item : info){
            String[] curItem = item.split("=");
            if(curItem.length>1) {
                table.put(curItem[0], curItem[1]);
            }
        }
        return table.get("DP");
    }
}
