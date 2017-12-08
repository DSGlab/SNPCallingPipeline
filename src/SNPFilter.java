import JavaBrew.Utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.Hashtable;
import java.util.Scanner;

/**
 * Created by juliofdiaz on 2017-09-22.
 *
 */
public class SNPFilter {
    //private static final String CONF_FILE = Utilities.CONF_FILE;

    private static final String DIR_SEPARATOR = Utilities.DIR_SEPARATOR;
    private static final String ID_LABEL_SEPARATOR= Utilities.ID_LABEL_SEPARATOR;
    private static final String SEPARATOR = Utilities.SEPARATOR;

    private static final String SNP_FILTER_INPUT_STRING = "SNP_FILTER_INPUT";
    private static final String SNP_FILTER_OUTPUT_STRING = "SNP_FILTER_OUTPUT";

    private static final String SNP_FILTER_MIN_GOOD_BALANCE_STRING = "SNP_FILTER_MIN_GOOD_BALANCE";
    private static final String SNP_FILTER_MIN_GOOD_CALL_RATIO_STRING = "SNP_FILTER_MIN_GOOD_CALL_RATIO";
    private static final String SNP_FILTER_MIN_GOOD_QUALITY_STRING = "SNP_FILTER_MIN_GOOD_QUALITY";

    public static void main(String[] args) throws FileNotFoundException {
        //Hashtable<String,String> options = Utilities.getConfigurationTable(CONF_FILE);
        Hashtable<String,String> options = Utilities.getConfigurationTable(args[0]);

        /* The inout file for this step */
        String INPUT = options.get( SNP_FILTER_INPUT_STRING );

        /* The output file for this step */
        String OUTPUT = options.get( SNP_FILTER_OUTPUT_STRING );

        /* Minimum required balance between forward and reverse reads */
        Integer SNP_FILTER_MIN_GOOD_BALANCE = Integer.parseInt(options.get(SNP_FILTER_MIN_GOOD_BALANCE_STRING));

        /*  Minimum required ratio of reads supporting the reference of the alternative */
        BigDecimal SNP_FILTER_MIN_GOOD_CALL_RATIO = new BigDecimal(options.get(SNP_FILTER_MIN_GOOD_CALL_RATIO_STRING));

        /* Minimum required quality for the SNP call */
        Integer SNP_FILTER_MIN_GOOD_QUALITY = Integer.parseInt(options.get(SNP_FILTER_MIN_GOOD_QUALITY_STRING));

        System.out.println("1. Reviewing raw SNP calls");
        Scanner in = new Scanner(new File(INPUT));
        PrintWriter out = new PrintWriter(OUTPUT);
        while(in.hasNextLine()){
            String[] data = in.nextLine().split(SEPARATOR);
            String isolate = data[0].split(DIR_SEPARATOR)[data[0].split(DIR_SEPARATOR).length-2].split("-")[1];
            String replicon = data[1];
            String position = data[2];
            String reference = data[3];
            String alternative = data[4];
            Double quality = Double.parseDouble( data[5] );
            //Integer depth = Integer.parseInt( data[6] );
            //Integer qualityDepth = Integer.parseInt( data[7] );
            Integer refForward = Integer.parseInt( data[9] );
            Integer refReverse = Integer.parseInt( data[10] );
            Integer altForward = Integer.parseInt( data[12] );
            Integer altReverse = Integer.parseInt( data[13] );

            Integer refTotal = refForward+refReverse;
            Integer altTotal = altForward+altReverse;

            Double refRatio = (0.0+refTotal)/(altTotal+refTotal);

            Double refRequiredRatio = SNP_FILTER_MIN_GOOD_CALL_RATIO.doubleValue();
            Double altRequiredRatio = new BigDecimal(1).subtract(SNP_FILTER_MIN_GOOD_CALL_RATIO).doubleValue();

            String type = (refRatio>=( refRequiredRatio ))? "ref" : ((refRatio<=( altRequiredRatio ))? "alt":"amb");

            Boolean goodQuality = (quality >= SNP_FILTER_MIN_GOOD_QUALITY);
            Boolean goodBalance = isGoodBalance(type,refForward,refReverse,altForward,altReverse,SNP_FILTER_MIN_GOOD_BALANCE);
            Boolean goodBoth = (goodQuality && goodBalance);

            String call = getCall(type,reference,alternative,goodBoth);
            String callFinal = (call.length()==1)&&(reference.length()==1)?call:"N";
            reference = (reference.length()==1)?reference:reference.charAt(0)+"";

            out.print(isolate + SEPARATOR + isolate + SEPARATOR +
                    replicon + SEPARATOR + position + SEPARATOR +
                    reference + SEPARATOR + callFinal );
            out.println();
        }
        out.close();

        System.out.println("2. The reviewed SNP calls are written in: " + OUTPUT);
        System.out.println();
    }


    /**
     * This method evaluates the SNP call and it returns the appropriate base. If the SNP
     * call is ambiguous or if the quality of the call is not good enough, then an N is
     * returned. If the call is strong enough to support the reference or the alternative
     * call, the the appropriate base is returned.
     *
     * @param type the type of mutation
     * @param ref the reference base
     * @param alt the alt base
     * @param goodBooth does it pass the quality controls
     * @return the final call
     */
    private static String getCall(String type, String ref,
                                  String alt, Boolean goodBooth){
        if(type.equals("amb")){
            return "N";
        }else if(type.equals("ref")){
            return ref;
        }else if(type.equals("alt")){
            return (goodBooth)? alt:"N";
        }else{
            return "?";
        }
    }

    /**
     * This method evaluates the support for the balance between reverse and forward reads. It
     * returns true if there is a minimum number of reads on both the reverse and forward
     * reads. This is evaluated only for the call that is accepted by the filtering.
     *
     * @param type the type of mutation
     * @param refForward the number of reference bases in the forward direction
     * @param refReverse the number of reference bases in the reverse direction
     * @param altForward the number of alt bases in the forward direction
     * @param altReverse the number of alt bases in the reverse direction
     * @return true if well balance between alternative and reference reads
     */
    private static Boolean isGoodBalance(String type, Integer refForward,
                                         Integer refReverse, Integer altForward,
                                         Integer altReverse, Integer minGoodBalance){
        if (type.equals("alt")) {
            return altForward >= minGoodBalance && altReverse >= minGoodBalance;
        } else
            return type.equals("ref") && refForward >= minGoodBalance && refReverse >= minGoodBalance;
    }
}
