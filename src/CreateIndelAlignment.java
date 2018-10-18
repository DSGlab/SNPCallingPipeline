import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Scanner;

public class CreateIndelAlignment {
    public static void main(String[]args) throws FileNotFoundException {
        String indelListFile = "/Users/juliofdiaz/Dropbox/CF/indel_calling/DOLOSA_2e/final_indel_list_noSTR.txt";
        String indelInfoFile = "/Users/juliofdiaz/Dropbox/CF/indel_calling/DOLOSA_2e/checked_indels_noSTR.txt";

        ArrayList<String> indelList = getIndelList(indelListFile);
        Hashtable<String, ArrayList<String>> indelCalls = getIndelCalls(indelInfoFile);
        System.out.println(indelList);
        System.out.println(indelCalls);
        System.out.println(indelCalls.get("1a").get(0));


        for(String i : new String[]{"1","2","3","4","5","6","7","8"}) {
            for (String j : new String[]{"a","b","c","d","e","f","g","h","i","j"}) {
                System.out.print(i+j+"\t");
                ArrayList<String> temp = indelCalls.get(i+j);
                for (String snp:indelList){
                    if(temp.contains(snp)){
                        System.out.print("1\t");
                    }else{
                        System.out.print("0\t");
                    }
                }
                System.out.println();
            }
        }
    }


    private static ArrayList<String> getIndelList(String indelListFile) throws FileNotFoundException {
        Scanner in = new Scanner(new File(indelListFile));

        ArrayList<String> indelList = new ArrayList<String>();
        while(in.hasNextLine()){
            String line = in.nextLine();
            indelList.add(line);
        }
        return indelList;
    }

    private static Hashtable<String, ArrayList<String>> getIndelCalls(String indelInfoFile) throws FileNotFoundException {
        Scanner in = new Scanner(new File(indelInfoFile));

        Hashtable<String, ArrayList<String>> result = new Hashtable<String, ArrayList<String>>();
        while(in.hasNextLine()){
            String[] line = in.nextLine().split("\t");
            String id = line[0].split("/")[7].split("\\.")[0];
            String curSNPId = line[1]+"-"+line[2];

            if(result.keySet().contains(id)){
                ArrayList<String> temp = result.get(id);
                temp.add(curSNPId);
                result.put(id,temp);
            }else{
                ArrayList<String> newList = new ArrayList<String>();
                newList.add(curSNPId);
                result.put(id,newList);
            }
        }
        return result;
    }
}
