import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Scanner;

public class Tester {
    public static void main(String[] args) throws FileNotFoundException {
        //Scanner in = new Scanner(new File("/Users/juliofdiaz/Downloads/Psy0524.vcf"));
        //Scanner in = new Scanner(new BufferedReader(new FileReader("/Users/juliofdiaz/Downloads/Psy0524.vcf")));
        //while(in.hasNextLine()){
        //    //String line = ;
        //    System.out.println(in.nextLine());

//        }


        System.out.println("1. Testing VCF File for Renanzinho.");
        VCFHolder test = new VCFHolder(new FileReader("/Users/juliofdiaz/Downloads/Psy0524.vcf"));
        System.out.println("2. Seeing inside VCF File");
        for(VCFVariant var : test.getVariants()){
            System.out.println(var.getChromosome()+"-"+var.getPosition()+" "+var.getReference()+">"+var.getAlternative());
        }
    }
}
