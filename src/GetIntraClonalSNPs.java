import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 1/11/16.
 *
 * @author juliofdiazc
 */
public class GetIntraClonalSNPs {
    public static void main ( String[] args ) throws Exception {
        ArrayList<String> HQsnpList = getSNPs("/Users/juliofdiaz/Dropbox/CF/snp_calling/DOLOSA_NEW/AllSNPList.txt");


        String[] ids = {"1a","1b","1c","1d","1e","1f","1g","1h","1i","1j","2a","2b","2c","2d","2e","2f","2g","2h","2i","2j","3a","3b","3c","3d","3e","3f","3g","3h","3i","3j","4a","4b","4c","4d","4e","4f","4g","4h","4i","4j","5a","5b","5c","5d","5e","5f","5g","5h","5i","5j","6a","6b","6c","6d","6e","6f","6g","6h","6i","6j","7a","7b","7c","7d","7e","7f","7g","7h","7i","7j","8a","8b","8c","8d","8e","8f","8g","8h","8i","8j"};
        String suffix = "/r.vcf";
        String workingDir = "/Users/juliofdiaz/Dropbox/CF/snp_calling/DOLOSA_V2/";

        LinkedHashMap<String, VCFHolder> snpCalls = getAllSNPCalls(workingDir, suffix, ids);

        //System.out.print(snpCalls.size());

        for(String snp: HQsnpList){
            String[] temp = snp.split("-");
            String tempContig = temp[0];
            int tempPos = Integer.parseInt( temp[1] );

            int count = 0;
            for(String id:snpCalls.keySet()){

                if( snpCalls.get(id).isVariant( tempContig, tempPos ) ){
                    count++;
                }
            }
            if ( count!= 80 && count!=0 ) {
                System.out.println( snp );
            }

        }


    }

    /**
     *
     * @param prefix
     * @param suffix
     * @param ids
     * @return
     * @throws FileNotFoundException
     */
    public static LinkedHashMap<String, VCFHolder> getAllSNPCalls(String prefix, String suffix, String[] ids)
            throws FileNotFoundException {

        LinkedHashMap<String, VCFHolder> snpCalls = new LinkedHashMap<String, VCFHolder>();

        for ( String id: ids ) {
            VCFHolder vcfh = new VCFHolder(prefix + id + "_NOVOALIGN" + suffix);
            snpCalls.put(id, vcfh);
        }
        return snpCalls;
    }

    /**
     * This method reads a file containing positions that will be reviewed and it returns
     * an ArrayList with them. The file needs to be a single text file where each line is
     * a variant position in the following format:
     * {contig/chromosome}-{position in contig/chromosome}
     *
     * The contig/chromosome identifier must be identical to the CHROM column from the vcf
     * format and the position identifier must be identical to the POS column from the vcf
     * format.
     *
     * @param fileName the name of the file containing the SNP positions that need to
     *                 be reviewed.
     * @return a list of all positions that need to be reviewed
     * @throws FileNotFoundException if file is not found
     */
    public static ArrayList<String> getSNPs ( String fileName )
            throws FileNotFoundException {
        Scanner in = new Scanner( new File( fileName ));

        ArrayList<String> result = new ArrayList<String>();
        while ( in.hasNextLine() ) {
            result.add(in.nextLine());
        }

        return result;
    }
}
