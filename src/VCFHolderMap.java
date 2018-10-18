import JavaBrew.FastaPrinter;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 * Modified by juliofdiaz on 11/11/14.
 *
 * This class is set to be implemented by another class. The user would instantiate a variant map
 * using the getVariantsMap() method. This method requires vcf files obtained from different
 * mapping methods using the same reads and same reference. The getVariantsMap() method returns a
 * map as a LinkedHashMap where the key is the variant identifier and the target is a list of the
 * variant calls from the different mapping methods.
 *
 * @author julio.diaz@mail.utoronto.ca
 *
 */
public class VCFHolderMap {
    public static void main(String[] args) throws Exception {


    }
    /**
     * This method surveys multiple vcf files from a single isolate and different mapping methods
     * (the reference needs to be the same). It creates a map where each variant holds the calls
     * from the different mapping methods.
     *
     * @param variantFiles list of the vcf files
     * @return a LinkedHashMap where the each key is any observed variant (Each variant is labeled as
     *         {contig/chromosome}-{position}). The target is a list of variants corresponding to said
     *         variant.
     * @throws Exception is thrown if any of the files is not found
     */
    public static LinkedHashMap<String,ArrayList<VCFVariant>> getVariantsMap ( String... variantFiles )
            throws Exception{
        LinkedHashMap<String, ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        ArrayList<VCFHolder> holders  = getVcfHolders( variantFiles );

        for ( VCFHolder v : holders ) {
            ArrayList<VCFVariant> vv = v.getVariants();
            for ( VCFVariant var : vv ) {
                if ( result.keySet().contains( var.getChromosome()+"-"+var.getPosition() ) ) {
                    result.get( var.getChromosome()+"-"+var.getPosition() ).add(var);
                }else{
                    ArrayList<VCFVariant> temp = new ArrayList<VCFVariant>();
                    temp.add( var );
                    result.put(var.getChromosome() + "-" + var.getPosition(), temp);
                }
            }
        }
        return result;
    }

    /**
     * This method removes variant calls that are clustered next to other variant calls
     * bellow a set minimum
     *
     * @param map maps of variants from all the sources
     * @param min minimum distance between variant calls
     * @return the variant map without variant calls bellow the threshold
     */
    public static LinkedHashMap<String, ArrayList<VCFVariant>> clusterFilterVariantsMap
            ( LinkedHashMap<String, ArrayList<VCFVariant>> map, Integer min ) {
        LinkedHashMap<String, ArrayList<Integer>> contigs = new LinkedHashMap<String, ArrayList<Integer>>();


        for ( String s : map.keySet() ) {
            String s2 = s.split("-")[0];
            if ( contigs.keySet().contains(s2) ) {
                //ADD IT TO ARRAYLIST OF THIS KEY
                ArrayList<Integer> tempIf = contigs.get( s2 );
                tempIf.add( Integer.parseInt( s.split( "-" )[1] ) );

                contigs.put( s2, tempIf );

            }else {
                //CREATE NEW KEY
                ArrayList <Integer> tempElse = new ArrayList<Integer>();
                tempElse.add( Integer.parseInt( s.split( "-" )[1] ) );

                contigs.put( s2,  tempElse);
            }
        }

        ArrayList<String> badList = new ArrayList<String>();
        for ( String c : contigs.keySet() ) {
            ArrayList<Integer> list = contigs.get( c );
            for(Integer i : list){
                for ( int k=i-min; k<=i+min ;k++ ) {
                    if ( list.contains(k) && k!=i) {
                        badList.add( c+"-"+k );
                    }
                }
            }
        }

        LinkedHashMap<String, ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        for ( String st : map.keySet() ) {
            if ( !badList.contains( st ) ) {
                result.put( st, map.get( st ) );
            }
        }

        return result;
    }

    /**
     * This method removes variant calls that are not supported by both forward and reverse
     * reads above the set minimum
     *
     * @param map maps of variants from all the sources
     * @param min minimum number of reads each direction
     * @return the variant map without variant calls bellow the threshold
     */
    public static LinkedHashMap<String, ArrayList<VCFVariant>> readBalanceFilterVariantsMap
            ( LinkedHashMap<String, ArrayList<VCFVariant>>  map, Integer min ) {
        LinkedHashMap<String,ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        for ( String key : map.keySet() ) {
            ArrayList<VCFVariant> pass = new ArrayList<VCFVariant>();

            for ( VCFVariant var : map.get( key ) ) {
                if ( var.getDp4AlternativeForward() > min && var.getDp4AlternativeReverse() > min ) {
                    pass.add( var );
                }
            }
            map.put(key,pass);

            if ( !pass.isEmpty() ) {
                result.put(key,pass);
            }
        }
        return result;
    }

    /**
     * This method removes variant calls that are too close to the edge of contigs/chromosomes.
     *
     * @param map the variant map
     * @param reference the reference used to map the reads and create the
     * @param min the minimum distance between the edge of a contig or chromosome
     * @return the variant map without variant calls below the threshold.
     */
    public static LinkedHashMap<String, ArrayList<VCFVariant>> contigEndFilterVariantsMap
            ( LinkedHashMap<String, ArrayList<VCFVariant>>  map, FastaPrinter reference, Integer min) {
        LinkedHashMap<String, Integer> ref = new LinkedHashMap<String, Integer>();
        for ( String k : reference.getSequences().keySet() ) {
            ref.put( k, reference.getSequences().get(k).length() );
        }

        LinkedHashMap<String,ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        for ( String key : map.keySet() ) {
            ArrayList<VCFVariant> pass = new ArrayList<VCFVariant>();

            for ( VCFVariant var : map.get( key ) ) {
                //System.out.println("1 "+var.getPosition());
                //System.out.println("2 "+var.getChromosome());
                //System.out.println("3 "+ref.get( var.getChromosome() ));
                //System.out.println("4 "+ref.keySet());

                if ( var.getPosition() > min && var.getPosition()<ref.get( var.getChromosome() )-min ) {
                    pass.add( var );
                }
            }
            map.put(key,pass);

            if ( !pass.isEmpty() ) {
                result.put(key,pass);
            }
        }
        return result;
    }

    /**
     * This method removes variant calls where the ratio of reads supporting the reference  to the
     * reads supporting the variant is under a given threshold.
     *
     * @param map the variant map
     * @param max the maximum ratio of reads supporting the reference to the reads supporting
     *            the alternative
     * @return the variant map without the variant calls bellow the threshold
     */
    public static LinkedHashMap<String, ArrayList<VCFVariant>> refToAltRatioFilterVariantsMap
            ( LinkedHashMap<String, ArrayList<VCFVariant>>  map, Double max) {
        LinkedHashMap<String,ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        for ( String key : map.keySet() ) {
            ArrayList<VCFVariant> pass = new ArrayList<VCFVariant>();

            for ( VCFVariant var : map.get( key ) ) {
                if ( var.getReferenceToAlternativeRatio() < max ) {
                    pass.add( var );
                }
            }

            if ( !pass.isEmpty() ) {
                result.put(key,pass);
            }
        }
        return result;
    }

    /**
     * This method removes variant calls that do not meet the reads depth threshold from
     * the variant map.
     *
     * @param map the variant map
     * @param min the minimum number of reads at the variant position
     * @return the variant map without the variant calls bellow the threshold
     */
    public static LinkedHashMap<String, ArrayList<VCFVariant>> qualityDepthFilterVariantsMap
            ( LinkedHashMap<String, ArrayList<VCFVariant>>  map, Integer min) {
        LinkedHashMap<String,ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        for ( String key : map.keySet() ) {
            ArrayList<VCFVariant> pass = new ArrayList<VCFVariant>();

            for ( VCFVariant var : map.get( key ) ) {
                if ( var.getQualityDepth() >= min ) {
                    pass.add( var );
                }
            }
            map.put(key,pass);

            if ( !pass.isEmpty() ) {
                result.put(key,pass);
            }
        }
        return result;
    }

    /**
     * This method removes variant calls that do not meet the raw reads depth threshold from
     * the variant map.
     *
     * @param map the variant map
     * @param min the minimum number of raw reads at the variant position
     * @return the variant map without the variant calls bellow the threshold
     */
    public static LinkedHashMap<String, ArrayList<VCFVariant>> depthFilterVariantsMap
            ( LinkedHashMap<String, ArrayList<VCFVariant>>  map, Integer min) {
        LinkedHashMap<String,ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        for ( String key : map.keySet() ) {
            ArrayList<VCFVariant> pass = new ArrayList<VCFVariant>();

            for ( VCFVariant var : map.get( key ) ) {
                if ( var.getDepth() >= min ) {
                    pass.add( var );
                }
            }
            map.put(key,pass);

            if ( !pass.isEmpty() ) {
                result.put(key,pass);
            }
        }
        return result;
    }

    /**
     * This method removes variant calls that do not meet the quality threshold from the variant
     * map.
     *
     * @param map the variant map
     * @param min the minimum variant phred score
     * @return the variant map without the variant calls bellow the threshold
     */
    public static LinkedHashMap<String, ArrayList<VCFVariant>> qualityFilterVariantsMap
            ( LinkedHashMap<String, ArrayList<VCFVariant>>  map, Integer min) {
        LinkedHashMap<String,ArrayList<VCFVariant>> result = new LinkedHashMap<String, ArrayList<VCFVariant>>();
        for ( String key : map.keySet() ) {
            ArrayList<VCFVariant> pass = new ArrayList<VCFVariant>();

            for ( VCFVariant var : map.get( key ) ) {
                if ( var.getQuality() >= min ) {
                    pass.add( var );
                }
            }
            map.put(key,pass);

            if ( !pass.isEmpty() ) {
                result.put(key,pass);
            }
        }
        return result;
    }

    /**
     * This method gathers the vcf files from the same isolate but from different
     * read mapping methods
     *
     * @param variantFiles list of the vcf files
     * @return a list of the vcf files as VCFHolders
     * @throws Exception is thrown if any of the files is not found
     */
    private static ArrayList<VCFHolder> getVcfHolders ( String... variantFiles )
            throws Exception {
        File file;
        VCFHolder vh;
        ArrayList<VCFHolder> vcfHolders = new ArrayList<VCFHolder>();

        for ( String temp : variantFiles ) {
            //file = new File( temp );
            vh = new VCFHolder(temp);
            vcfHolders.add(vh);
        }

        return vcfHolders;
    }

}