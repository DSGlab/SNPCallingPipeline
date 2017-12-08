import JavaBrew.Utilities;

public class main {
    public static void main(String[] args) throws Exception {
        System.out.println("SNP Calling Pipeline v. 1.07");
        System.out.println("Questions? julio.diaz@mail.utoronto.ca");
        System.out.println();
        System.out.println("Usage:\tjava -jar SNPCallingPipeline ANALYSIS_NAME [config file]");
        System.out.println();

        String confFile = System.getProperty("user.dir")+Utilities.DIR_SEPARATOR+Utilities.DEFAULT_CONF_NAME;

        switch (args.length){
            case 0:
                System.out.println("Step 1:\t\tgethqsnps\t\tGet list of High quuality (confidence) SNPs");
                System.out.println("Optional step:\tgetintraclonalsnps\tRemove SNPs that are fixed between the isolates and the reference");
                System.out.println("Step 2:\t\tsnpchecker\t\tGets raw calls at the SNP positions in all the isolates");
                System.out.println("Step 3:\t\tsnpfilter\t\tFilters raw SNP calls");
                System.out.println("Step 4:\t\tcreatealignment\t\tCreates alignment based on filtered SNP calls");
                break;
            case 1:
                runAnalysis(args[0], confFile);
                break;
            case 2:
                runAnalysis(args[0],args[1]);
                break;
            default:
                System.out.println("max two arguments needed");
                break;
        }

    }

    private static void runAnalysis(String analysis, String configFile) throws Exception {
        String[] input = {configFile};
        try {
            Analysis analysisName = Analysis.valueOf( analysis );

            switch (analysisName) {
                case gethqsnps:
                    System.out.println("Starting Get_HQ-SNPs analysis.");
                    GetHQSNPs.main(input);
                    break;
                case getintraclonalsnps:
                    System.out.println("Starting Get_IntraClonal-SNPs analysis.");
                    GetIntraClonalSNPs.main(input);
                    break;
                case snpchecker:
                    System.out.println("Starting SNP_Checker analysis.");
                    SNPChecker.main(input);
                    break;
                case snpfilter:
                    System.out.println("Starting SNP_Filter analysis.");
                    SNPFilter.main(input);
                    break;
                case createalignment:
                    System.out.println("Starting Create_Alignment analysis.");
                    CreateAlignmentFromSNPChecked.main(input);
                    break;
                case scinetjobcreator:
                    System.out.println("Starting Create_Alignment analysis.");
                    CreateAlignerJobs.main(input);
                    break;
                default:
                    System.out.println("\nwhat are you doing here?\n");
                    break;
            }
        }catch (IllegalArgumentException e){
            System.out.println( e.toString() );
            System.out.println("\"" + analysis + "\" is not an analysis.");
        }
    }

    public enum Analysis{
        gethqsnps, getintraclonalsnps, snpchecker, snpfilter, createalignment, scinetjobcreator
    }
}
