package antibody;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;

import basicUtils.Utils;
import errorUtils.ErrorThrower;

public class GenerateAllPRMs {

	public static String UsageInfo = "Usage: java -jar GenerateAllPRMs.jar version 2011.10.18\n"

			+ "[REQUIRED]:\n"
			+ "-r [Directory] Directory containing the spectrum files\n"
			+ "-w [Directory] Directory to write the prm files\n"
			+ "-x [FILE] Path to PepNovo executable\n"
			+ "-m [DIR] Directory containing PepNovo models\n"
			+ "[OPTIONAL]:\n"
			+ "-y [NUM] Specifies the parent mass tolerance (in Daltons).  Default is "
			+ AntibodyUtils.DefaultPMTolerance
			+ " Da.\n"
			+ "-z [NUM] Specifies the fragment mass tolerance (in Daltons).  Default is "
			+ AntibodyUtils.DefaultFragTolerance
			+ " Da.\n"
			+ "-f The instrument is high accuracy (e.g. FT)\n"
			+ "-t [FILE] File containing which spectrum files and spectra to use\n"
			+ "-c [NUM] Mass in Da of the cysteine protecting group (Default is 0)\n"
			+ "-p Tryptic digest was used on the samples.  Default is to guess from the file names.\n"
			+ "-a [NUM] Activation type.  0:CID (default), 1:ETD, 2:HCD, 3:Multi/Paired, 4:Unpaired\n"
			+ "-q [NUM] If merging pairs/triplets, number of consecutive spectra to merged";

	/*
	 * public static String PepNovoDir; public static String ModelDir; public
	 * static String InputDir; public static String OutputDir; public static
	 * String Digest; public static int CysteineMass; public static int
	 * Instrument = 1; // 0 = LTQ 1 = ORBI
	 * 
	 * public static String InputFile = null; public static ArrayList
	 * InputFiles;
	 */

	public static void RunPepNovoOnDir(String inputDir, String outputDir,
			String pepNovoExec, String modelDir, int cysteineMass,
			String inputFile, String digest, double pMTol, double fragTol,
			int actType, int numConsec, boolean isFT) {

		String[] specFiles = null;
		if (inputFile != null)
			specFiles = AntibodyUtils.LoadLinesFromFile(inputFile);
		else {
			specFiles = AntibodyUtils.ListSpectrumDir(inputDir);
		}
		if (specFiles.length == 0) {
			ErrorThrower.ThrowError(4,
					"Unable to load any files from spectrum directory "
							+ inputDir);
		}

		System.out.println("Total spectrum files: " + specFiles.length);

		basicUtils.Utils.activationType aType = basicUtils.Utils.activationType.CID;
		if (actType == 0)
			aType = basicUtils.Utils.activationType.CID;
		else if (actType == 1)
			aType = basicUtils.Utils.activationType.ETD;
		else if (actType == 2)
			aType = basicUtils.Utils.activationType.HCD;
		else if (actType == 3) {
			aType = basicUtils.Utils.activationType.MultiPaired;
			if (numConsec <= 0) {
				ErrorThrower
						.ThrowError(4,
								"Number of consecutive spectra (-q) to merge must be 1 or more");
			}
		} else if (actType == 4)
			aType = basicUtils.Utils.activationType.Unpaired;
		else {
			ErrorThrower
					.ThrowError(4,
							"Activation type (-a) must be 0:CID, 1:ETD, 2:HCD, 3:Multi/Paired, 4:Unpaired");
		}

		for (int i = 0; i < specFiles.length; ++i) {

			String OutputFile = outputDir + File.separator
					+ Utils.GetBaseNameNoExtension(specFiles[i]) + ".prm";
			String errFileName = outputDir + File.separator
					+ Utils.GetBaseNameNoExtension(specFiles[i]) + ".err";
			String cMass = null;
			if (cysteineMass > 0)
				cMass = cysteineMass + "";
			RunPepNovo(specFiles[i], OutputFile, errFileName, pepNovoExec,
					pMTol, fragTol, modelDir, cMass, digest, aType, isFT);
		}
	}

	public static void RunPepNovo(String SpecFileName, String OutputFile,
			String errFileName, String PepNovoExec, double PMTolerance,
			double FragTolerance, String PepNovoModelDir, String CysteineMass,
			String Digest, basicUtils.Utils.activationType actType, boolean isFT) {
		String FullInputFile = SpecFileName;
		String Command = "";

		if (actType == basicUtils.Utils.activationType.Unpaired) {
			/*
			 * RunPepNovoSeparate(SpecFileName, OutputFile, errFileName,
			 * PepNovoExec, PMTolerance, FragTolerance, PepNovoModelDir,
			 * CysteineMass, Digest);
			 */
			System.out.println("This parameters set is under construction");
			return;
		} else if (actType == basicUtils.Utils.activationType.MultiPaired) {
			/*
			 * RunPepNovoTogether(SpecFileName, OutputFile, errFileName,
			 * PepNovoExec, PMTolerance, FragTolerance, PepNovoModelDir,
			 * CysteineMass, Digest, numConsec);
			 */
			System.out.println("This parameters set is under construction");
			return;
		}
		if (Digest != null) {
			if (Digest.toLowerCase().compareTo("trypsin") == 0)
				Digest = "TRYPSIN";
			else
				Digest = "NON_SPECIFIC";
		}

		Command = PepNovoExec + " ";
		if (actType == basicUtils.Utils.activationType.CID && isFT) {
			Command += " -model CID_FT_Tryp_z";
		} else if (actType == basicUtils.Utils.activationType.CID) {
			Command += " -model CID_IT_TRYP";
		} else {
			System.out.println("This parameters set is under construction");
			return;
		}

		Command += " -prm_norm -pm_tolerance " + PMTolerance
				+ " -fragment_tolerance " + FragTolerance + " -model_dir "
				+ PepNovoModelDir;

		if (CysteineMass != null)
			Command += " -PTMs C+" + CysteineMass;
		if (Digest == null)
			Command += " -digest "
					+ AntibodyUtils.DetermineDigestTypePepNovo(FullInputFile)
					+ " -file " + FullInputFile;
		// + " > " + OutputFile;
		else
			Command += " -digest " + Digest + " -file " + FullInputFile;
		// + " > " + OutputFile;
		
		try {
			Utils.RunCommandWithStdOutToFile(Command, OutputFile,
					null, Utils.GetFilePath(PepNovoExec));
		} catch (Exception E) {
			E.printStackTrace();
		}
	}

	/*
	 * + "Usage: java -jar GenerateAllPRMs.jar\n" + "[REQUIRED]:" +
	 * "-r [Directory] Directory containing the spectrum files\n" +
	 * "-w [Directory] Directory to write the prm files\n" +
	 * "-x [FILE] Path to PepNovo executable\n" +
	 * "-m [DIR] Directory containing PepNovo models\n" + "[OPTIONAL]:\n" +
	 * "-y [NUM] Specifies the parent mass tolerance (in Daltons).  Default is 3.0 Da.\n"
	 * +
	 * "-z [NUM] Specifies the fragment mass tolerance (in Daltson).  Default is 0.5 Da.\n"
	 * + "-t [FILE] File containing which spectrum files and spectra to use\n" +
	 * "-c [NUM] Mass in Da of the cysteine protecting group.  Only integral masses are accepted.\n"
	 * +
	 * "-p Tryptic digest was used on the samples.  Default is to guess from the file names.\n"
	 * ;
	 */
	public static void main(String[] args) {
		String[] Commands = { "-r", "-w", "-x", "-m", "-y", "-z", "-t", "-c",
				"-p", "-a", "-q", "-f" };
		boolean[] Values = { true, true, true, true, true, true, true, true,
				false, true, true, false };
		Hashtable Options = Utils.ParseCommandLine(args, Commands, Values);

		if (!Options.containsKey("-r") || !Options.containsKey("-w")
				|| !Options.containsKey("-x") || !Options.containsKey("-m")) {
			System.err
					.println("ERROR: Must specify an input directory, output directory, and path to PepNovo executable and model dir");
			System.out.println(UsageInfo);
			return;
		}

		String InputDir = (String) (Options.get("-r"));
		String OutputDir = (String) (Options.get("-w"));
		String PepNovoDir = (String) (Options.get("-x"));
		String ModelDir = (String) (Options.get("-m"));
		int CysteineMass = 0;
		if (Options.containsKey("-c"))
			CysteineMass = Integer.parseInt((String) (Options.get("-c")));

		String InputFile = null;
		if (Options.containsKey("-t"))
			InputFile = (String) (Options.get("-t"));

		String Digest = null;
		if (Options.containsKey("-p"))
			Digest = ((String) (Options.get("-p"))).toUpperCase();

		double PMTol = AntibodyUtils.DefaultPMTolerance;
		double FragTol = AntibodyUtils.DefaultFragTolerance;
		if (Options.containsKey("-y"))
			PMTol = Double.parseDouble((String) (Options.get("-y")));
		if (Options.containsKey("-z"))
			FragTol = Double.parseDouble((String) (Options.get("-z")));

		int actType = 0;
		int numConsec = 1;
		if (Options.containsKey("-a"))
			actType = Integer.parseInt((String) (Options.get("-a")));
		if (Options.containsKey("-q"))
			numConsec = Integer.parseInt((String) (Options.get("-q")));

		boolean isFT = false;
		if (Options.containsKey("-f"))
			isFT = true;

		GenerateAllPRMs.RunPepNovoOnDir(InputDir, OutputDir, PepNovoDir,
				ModelDir, CysteineMass, InputFile, Digest, PMTol, FragTol,
				actType, numConsec, isFT);

	}

}
