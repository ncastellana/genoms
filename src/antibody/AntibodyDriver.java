package antibody;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;

import specLibUtils.SpectrumLibraryAnnotation;
import trieUtils.TrieDB;
import basicUtils.MSGFDBAnnotation;
import basicUtils.PeptideSpectrumMatch;
import basicUtils.Utils;
import edu.ucsd.msjava.ui.MSGFDB;
import errorUtils.ErrorThrower;
//import basicUtils.InspectAnnotation;

//import uk.ac.ebi.jmzidml.model.mzidml.AbstractParam;

/**
 * This is the main driver class for the antibody project. In reality it can be
 * used for any set of ordered sequences.
 * 
 * Requirements: Inspect v2009.01.22 or higher
 * 
 * Updates
 * 
 * 2012-05-15 Updated load PRMs to accept MGF files in addition to PepNovo
 * output files 2012-07-16 Make input and output file absolute before running.
 * 2012-07-25 Removed 'O' as a possible amino acid 2012-0825 Added direct call
 * to MSGFDB, instead of system call. 2012-09-07 Fixed file presence checks for
 * DB files in constructor. Also fixed dealing with empty PRM spectra.
 * 2012-10-02 Added cygwin compatibility 2013-01-07 Added activation type
 * capability for CID, HCD, ETD, and pairs/triplets 2013-05-09 Added ability to
 * interpret spectra when PrmClust is used 2013-12-02 Updated PRMClust
 * interpretation to be based on parent mass (not enforcing strict num
 * consecutive) 2013-12-16 Updated to parse new ClusterSet data for
 * PRMClustering 2013-12-19 Updated to include spectral library search results
 * 2015-08-25 Updated constraint file format to include the class type
 * 2015-10-10 Added PepMass to sps_seqs.mgf to fix offset in reports
 * 2015-11-05 Added C+58, C+105 to cysteine mods
 */
public class AntibodyDriver {

	public static boolean Debug = true;

	public static final String UsageInfo = "---- GenoMS Version 1.0 "
			+ AntibodyUtils.VersionInfo
			+ " ----\n"
			+ "Usage: java -jar GenoMS.jar\n"
			+ "[REQUIRED]:\n"
			+ "-i [FILE] The input configuration file\n"
			+ "-o [FILE] The output file\n"
			+ "[OPTIONAL]:\n"
			+ "-r [DIR] The resource directory to look for things like DBs (default .)\n"
			+ "-k [DIR] The project directory to write all output stuff to (default .)\n"
			+ "-e [DIR] The execution directory where all executables are found (default .)\n"
			+ "-x Force a rerun DB searches\n"
			+ "-p Use a missing peak penalty when scoring HMM Match states\n"
			+ "-s Require all2all similarity for Multiple Spectrum Alignment\n"
			+ "-a Continue to add eligible match states, even after all spectra are added to HMM\n"
			+ "-w [NUM] FDR cutoff for DB search results (Default is 0.01)\n"
			+ "-q Generate files for report generating (NOTE: Only use this for spec net integration.  \n"
			+ "   It makes the assumption that the input spectra are in a single file.\n"
			+ "-l [FILE] Write output to a log file (default is to stdout)\n";

	/*
	 * These parameters are specified on the command line
	 */
	// Config file, see documentation
	private String configFileName;
	// Output file
	private String outputFileName;
	private String detailsOutputFileName;

	private String reportFileDir = null;

	// private int MaxXLSWidth = 30;

	private boolean AddMatchStatesAtEndFlag = false;
	private double fdrCutoff = 0.01;

	private String PepNovoModelDir = null;
	private String PepNovoExec = null;

	// Data parsed from the config file
	private String[] SpectraFileNames; // The files/directories containing the
	// spectra for this experiment

	// to the spectra (omitting the DB)
	private String[] dbSearchOutputFiles;
	private String[] DatabaseRootName; // The root name(s) for the databases,
	// see documentation

	// database
	private String CombinedShuffledDBName = null; // The trie file of the
	// combined shuffled
	// database
	private boolean runDBSearchFlag = false;
	private String[] PRMFileNames; // The file/directory containing the
	// PRMspectra (produced from PepNovo) for
	// the spectra in this experiment.

	private boolean ClusterFlag = false; // Require Spectra to have similarity
	// to previously aligned spectra in
	// HMM

	private TrieDB trieDatabase = null;
	// private InspectAnnotation[] AllAnnotations;
	private int annCount = 0; // THe number of spectrum identified

	private SpectrumNode[] allPRMSpectra; // One SpectrumNode per PRM spectrum

	private Hashtable<String, SpectrumNode> spectrumNodeIndex; // SpectrumFile->ScanNumber->Index
																// into
	// AllSpectra
	private double[] FixedMods = null;

	// New stuf added to get rid of Template classes
	private String TemplateConstraintFile;

	// Each SelectedTemplates[i] ->Template->Template...
	private Template[] SelectedTemplates;
	private Template[] allTemplates;
	private int MinPeptidesPerPath = 2;

	private double PMTolerance = AntibodyUtils.DefaultPMTolerance;
	private double FragTolerance = AntibodyUtils.DefaultFragTolerance;
	private double PRMTolerance = AntibodyUtils.DefaultFragTolerance;

	private String Digest = "other";

	private String FASTADatabase = null;

	private String projectDir;

	private String exeDir;

	private String ResourceDir;

	// By default we treat each spectrum individually, as CID spectra
	private basicUtils.Utils.activationType actType = basicUtils.Utils.activationType.CID;;
	// private int numConsec = 1; //Number of consecutive spectra from the same
	// precursor

	private boolean isFT = false; // Is high accuracy instrument
	private int clusterMinSize = -1; // Min size clusters to allow
	private String clusterTool = null; // cluster tool

	private boolean unRavelPRMCluster = false; // Do we need to unravel the prms
												// from the same peptide
	private boolean mergeSamePrec = false;
	private Hashtable<Integer, int[]> PRM2SpectrumMap = null; // from PrmClust
	// private Hashtable<Integer, Integer> SpecFile2SpectrumMap = null;
	private String inst = "LTQ";

	private ArrayList<antibody.Peptide> dbPepObj;

	// Spectral library params
	private boolean runLibrarySearch; // Was a library search run prior to
										// running GenoMS?
	private String specLibSearchResultsFile; // Where were the library search
												// results written?

	
	private boolean applyMatchPenalty = false;
	private double matchPenalty;
	
	// private int[] SpecNode2DBSearchResults = null;
	// private ArrayList[] DBResult2SpecNodeKey = null;
	// private Hashtable SpecNodeKey2DBResult = null;

	// private int HistogramNum = 0;

	// private String InstrumentType;

	/**
	 * 
	 * @return Returns true if the config file contains the required parameters
	 *         (see documentation).
	 */
	private boolean ParseConfigFile() {
		ArrayList<String[]> Values;

		String[] Options = { "exe_dir", "spectra", "dbcombined", "dbrootname",
				"prms", "fixedmod", "templateconstraintfile", "pepnovoexec",
				"pepnovomodeldir", "tolerance_pm", "tolerance_peak", "digest",
				"activation", "intrument_type", "cluster_min_size",
				"cluster_tool", "merge_same_prec", "genoms_instrument",
				"results_dir", "run_library_search" };

		int[] NumValues = { 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1 };

		Hashtable<String, ArrayList<String[]>> Parameters = null;
		if (this.reportFileDir != null)
			Parameters = Utils.ParseCSPSInputFile(this.configFileName, Options,
					NumValues);
		else
			Parameters = Utils.ParseInputFile(this.configFileName, Options,
					NumValues);
		boolean runPRMGenerator = false;

		// ////////////////////////////////////////////////////////
		// / Load params related to spectral library search
		// ////////////////////////////////////////////////////////

		if (Parameters.containsKey("run_library_search")) {
			Values = (ArrayList<String[]>) Parameters.get("run_library_search");
			int runSearch = Integer.parseInt(((String[]) (Values.get(0)))[0]);

			if (runSearch == 1) {
				this.runLibrarySearch = true;
				if (Parameters.containsKey("results_dir")) {
					Values = (ArrayList<String[]>) Parameters
							.get("results_dir");
					this.specLibSearchResultsFile = ((String[]) (Values.get(0)))[0];
				} else {
					ErrorThrower
							.ThrowError(4,
									"Must specify 'results_dir' when running library search");
				}
			} else
				this.runLibrarySearch = false;

		}

		// ////////////////////////////////////////////////////////
		// / Load params related to clustering
		// ////////////////////////////////////////////////////////
		if (Parameters.containsKey("cluster_min_size")) {
			Values = (ArrayList<String[]>) Parameters.get("cluster_min_size");
			this.clusterMinSize = Integer
					.parseInt(((String[]) (Values.get(0)))[0]);
		}
		if (Parameters.containsKey("genoms_instrument")) {
			Values = (ArrayList<String[]>) Parameters.get("genoms_instrument");
			this.inst = ((String[]) (Values.get(0)))[0].toUpperCase();

			if (!this.inst.equals("LTQ") && !this.inst.equals("FT")
					&& !this.inst.equals("TOF")) {
				ErrorThrower
						.ThrowError(
								4,
								"'"
										+ this.inst
										+ "''.  GenoMS_Instrument must be 'LTQ', 'FT', or 'TOF'");
			}
		}

		if (Parameters.containsKey("cluster_tool")) {
			Values = (ArrayList<String[]>) Parameters.get("cluster_tool");
			this.clusterTool = ((String[]) (Values.get(0)))[0].toUpperCase();
			if (this.clusterTool.equals("PRMCLUST") && this.clusterMinSize > 0) {
				// ErrorThrower.ThrowError(22,
				// "CLUSTER_TOOL:PRMCLUST is not implemented in this version");
				this.unRavelPRMCluster = true;
			}
		}

		if (Parameters.containsKey("merge_same_prec")) {
			Values = (ArrayList<String[]>) Parameters.get("merge_same_prec");
			if (Integer.parseInt(((String[]) (Values.get(0)))[0]) == 1) {
				this.mergeSamePrec = true;
				this.unRavelPRMCluster = true;
			} else
				this.mergeSamePrec = false;

		}

		if (Parameters.containsKey("instrument_type")) {
			Values = (ArrayList<String[]>) Parameters.get("instrument_type");
			String type = ((String[]) (Values.get(0)))[0].toUpperCase();
			if (type.equals("FT"))
				this.isFT = true;
		}

		/*
		 * if(Parameters.containsKey("num_consecutive")) { Values =
		 * (ArrayList<String[]>) Parameters.get("num_consecutive");
		 * this.numConsec = Integer.parseInt(((String[])(Values.get(0)))[0]);
		 * if(this.numConsec < 1) ErrorThrower.ThrowError(4, "'" +
		 * this.numConsec +
		 * "'.  Number of consecutive spectra to merge must be greater than 0");
		 * }
		 */

		// If we are running with SPS, then activation will not be explicitly
		// set, so lets figure it out
		if (this.mergeSamePrec)
			this.actType = basicUtils.Utils.activationType.MultiPaired;
		else
			this.actType = basicUtils.Utils.activationType.Unpaired;

		if (Parameters.containsKey("activation")) {
			Values = (ArrayList<String[]>) Parameters.get("activation");
			String act = ((String[]) (Values.get(0)))[0].toUpperCase();
			if (act.equals("CID"))
				this.actType = basicUtils.Utils.activationType.CID;
			else if (act.equals("ETD"))
				this.actType = basicUtils.Utils.activationType.ETD;
			else if (act.equals("HCD"))
				this.actType = basicUtils.Utils.activationType.HCD;
			else if (act.equals("MULTIPAIRED"))
				this.actType = basicUtils.Utils.activationType.MultiPaired;
			else if (act.equals("UNPAIRED"))
				this.actType = basicUtils.Utils.activationType.Unpaired;
			else {
				ErrorThrower
						.ThrowError(
								4,
								"'"
										+ act
										+ "'.  Activation type must be 'CID','ETD','HCD','MultiPaired', or 'Unpaired'");
			}
		} else
			System.out.println("Using default activiation: " + this.actType);
		/**
		 * Specify a fixed modification on an amino acid
		 */
		if (Parameters.containsKey("fixedmod")) {

			Values = (ArrayList<String[]>) Parameters.get("fixedmod");
			char AA = '0';
			for (int i = 0; i < Values.size(); ++i) {
				try {
					AA = (((String[]) Values.get(i))[0]).toUpperCase()
							.charAt(0);
				} catch (Exception E) {
					ErrorThrower.ThrowErrorCustum(
							ErrorThrower.CUSTOM_ERROR,
							"FixedMod is ill formed: "
									+ ((String[]) Values.get(i))[0]);
				}
				double Mass = 0;
				try {
					Mass = Double.parseDouble(((String[]) Values.get(i))[1]);
				} catch (Exception E) {
					System.err.println("ERROR: FixedMod is ill formed: "
							+ ((String[]) Values.get(i))[1]);
				}

				if (AA < 65 || AA >= 91)
					ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
							"FixedMod on AA " + AA + " is invalid");
				if (this.FixedMods != null && this.FixedMods[AA - 65] != 0)
					continue;
				// System.out.println("FIXEDMOD: " + AA + " + " + Mass);
				AntibodyUtils.AAMasses[AA - 65] += Mass;
				AntibodyUtils.AAMassesScaled[AA - 65] += (int) (Mass * AntibodyUtils.MASS_SCALE);
				if (this.FixedMods == null)
					this.FixedMods = new double[AntibodyUtils.AAMasses.length];
				this.FixedMods[AA - 65] = Mass;
			}
		}

		/**
		 * Find a database in whatever form
		 */
		if (Parameters.containsKey("templateconstraintfile")) {
			Values = (ArrayList<String[]>) Parameters
					.get("templateconstraintfile");
			String FileName = ((String[]) Values.get(0))[0];
			this.TemplateConstraintFile = AntibodyUtils.makePathAbsolute(
					FileName, this.projectDir);
			File test = new File(this.TemplateConstraintFile);
			if (!test.exists()) {
				this.TemplateConstraintFile = AntibodyUtils.makePathAbsolute(
						FileName, this.ResourceDir);
				test = new File(this.TemplateConstraintFile);
				if (!test.exists())
					ErrorThrower.ThrowError(1, this.TemplateConstraintFile);
			}

		}
		// Look for a DB, must have specified either dbcombined or dbrootname
		if (Parameters.containsKey("dbcombined")) {
			if (!Parameters.containsKey("templateconstraintfile")) {
				ErrorThrower
						.ThrowError(3,
								"Cannot specify a combined db without a template contraint file!!");

			}

			Values = (ArrayList<String[]>) Parameters.get("dbcombined");
			String FileName = ((String[]) Values.get(0))[0];
			if (Utils.GetFileExtension(FileName).compareTo("trie") == 0) {
				// System.out.println("Great!DB seems to be formatted!");
				this.CombinedShuffledDBName = AntibodyUtils.makePathAbsolute(
						FileName, this.projectDir);
				File test = new File(this.CombinedShuffledDBName);
				if (!test.exists()) {
					this.CombinedShuffledDBName = AntibodyUtils
							.makePathAbsolute(FileName, this.ResourceDir);
					test = new File(this.FASTADatabase);
					if (!test.exists())
						ErrorThrower.ThrowError(1, this.CombinedShuffledDBName);
				}
				this.trieDatabase = new TrieDB(this.CombinedShuffledDBName);

			} else {
				// CombinedDBName = FileName;
				this.FASTADatabase = AntibodyUtils.makePathAbsolute(FileName,
						this.projectDir);
				File test = new File(this.FASTADatabase);
				if (!test.exists()) {
					this.FASTADatabase = AntibodyUtils.makePathAbsolute(
							FileName, this.ResourceDir);
					test = new File(this.FASTADatabase);
					if (!test.exists())
						ErrorThrower.ThrowError(1, this.FASTADatabase);
				}
			}
		}

		else if (Parameters.containsKey("dbrootname")) {
			Values = (ArrayList<String[]>) Parameters.get("dbrootname");
			ArrayList<String> allNames = new ArrayList<String>();
			for (int i = 0; i < Values.size(); ++i) {
				String fullVal = ((String[]) (Values.get(i)))[0];
				// System.out.println(fullVal);
				String[] vals = (fullVal).split(";");
				// System.out.println(vals.length);
				for (int j = 0; j < vals.length; ++j)
					allNames.add(vals[j]);
			}
			this.DatabaseRootName = new String[allNames.size()];
			for (int i = 0; i < this.DatabaseRootName.length; ++i) {
				String FileName = (String) (allNames.get(i));
				this.DatabaseRootName[i] = AntibodyUtils.makePathAbsolute(
						FileName, this.projectDir);

			}
		} else {
			ErrorThrower.ThrowError(3, "Must specify either DBCombined or "
					+ "DBRootName in Config file");
		}

		/**
		 * Determine if there are PRMS
		 */
		if (!Parameters.containsKey("prms")) {
			if (!Parameters.containsKey("pepnovoexec")
					|| !Parameters.containsKey("pepnovomodeldir")) {
				ErrorThrower.ThrowError(3,
						"Must specify either prms or path to PepNovo");
			}
			Values = (ArrayList) Parameters.get("pepnovoexec");
			String FileName = ((String[]) (Values.get(0)))[0];
			this.PepNovoExec = AntibodyUtils.makePathAbsolute(FileName,
					this.projectDir);
			File test = new File(this.PepNovoExec);
			if (!test.exists())
				this.PepNovoExec = AntibodyUtils.makePathAbsolute(FileName,
						this.exeDir);

			Values = (ArrayList) Parameters.get("pepnovomodeldir");
			FileName = ((String[]) (Values.get(0)))[0];
			this.PepNovoModelDir = AntibodyUtils.makePathAbsolute(FileName,
					this.projectDir);
			test = new File(this.PepNovoModelDir);
			if (!test.exists())
				this.PepNovoModelDir = AntibodyUtils.makePathAbsolute(FileName,
						this.exeDir);

			runPRMGenerator = true;
		} else {
			Values = (ArrayList<String[]>) Parameters.get("prms");

			this.PRMFileNames = new String[Values.size()];
			for (int i = 0; i < Values.size(); ++i) {
				this.PRMFileNames[i] = AntibodyUtils.makePathAbsolute(
						((String[]) (Values.get(i)))[0], this.projectDir);
			}
		}

		/*
		 * Look for spectra
		 */
		if (!Parameters.containsKey("spectra")) {
			ErrorThrower.ThrowError(3, "Must specify Spectra in Config file");
		}

		Values = (ArrayList<String[]>) Parameters.get("spectra");
		this.SpectraFileNames = new String[Values.size()];
		if (this.SpectraFileNames.length > 1) {
			ErrorThrower.ThrowError(4,
					"MSGFDB can only be run on a single spectrum file!");
		}
		for (int i = 0; i < Values.size(); ++i) {
			this.SpectraFileNames[i] = AntibodyUtils.makePathAbsolute(
					((String[]) (Values.get(i)))[0], this.projectDir);
		}

		if (Parameters.containsKey("tolerance_pm")) {
			Values = (ArrayList<String[]>) Parameters.get("tolerance_pm");
			this.PMTolerance = Double
					.parseDouble(((String[]) (Values.get(0)))[0]);
		}

		if (Parameters.containsKey("tolerance_peak")) {
			Values = (ArrayList<String[]>) Parameters.get("tolerance_peak");
			this.FragTolerance = Double
					.parseDouble(((String[]) (Values.get(0)))[0]);
		}

		if (Parameters.containsKey("digest")) {
			Values = (ArrayList<String[]>) Parameters.get("digest");
			this.Digest = ((String[]) (Values.get(0)))[0];
		}

		if (runPRMGenerator) {
			System.out
					.println("\n[GenoMS Preprocessing]: Running PepNovo+ to generate PRM Spectra");
			String CysteineMass = null;
			this.PRMFileNames = new String[this.SpectraFileNames.length];
			if (this.FixedMods != null && this.FixedMods['C' - 65] != 0) {

				int Mass = (int) (this.FixedMods['C' - 65]);
				CysteineMass = (Mass) + "";
			}
			for (int i = 0; i < this.SpectraFileNames.length; ++i) {
				String OutputFile = Utils
						.GetFileNameNoExtension(this.SpectraFileNames[i])
						+ ".prm";
				String errFile = Utils
						.GetFileNameNoExtension(this.SpectraFileNames[i])
						+ ".err";
				this.PRMFileNames[i] = OutputFile;
				GenerateAllPRMs.RunPepNovo(this.SpectraFileNames[i],
						OutputFile, errFile, this.PepNovoExec,
						this.PMTolerance, this.FragTolerance,
						this.PepNovoModelDir, CysteineMass, this.Digest,
						this.actType, this.isFT);
			}
		}

		if (Debug) {
			if (this.DatabaseRootName != null) {
				for (int k = 0; k < this.DatabaseRootName.length; ++k)
					System.out.println("DBRoot: " + this.DatabaseRootName[k]);
			}
			if (this.CombinedShuffledDBName != null)
				System.out
						.println("CombinedDB: " + this.CombinedShuffledDBName);

			System.out.println("RunDBSearchFlag: " + this.runDBSearchFlag);

			// for(int i = 0; i < this.SpectraFileNames.length; ++i)

			// System.out.println("Spectra: " + this.SpectraFileNames[i] +
			// " Config: " + this.InspectConfigFiles[i]);
			// for(int i = 0; i < this.PRMFileNames.length; ++i)
			// System.out.println("PRMs: " + this.PRMFileNames[i]);
			if (this.ResourceDir != null)
				System.out.println("ResourceDir: " + this.ResourceDir);
			if (this.projectDir != null)
				System.out.println("ProjectDir: " + this.projectDir);
			if (this.exeDir != null)
				System.out.println("Exe Dir: " + this.exeDir);

		}
		return true;
	}

	private PeptideSpectrumMatch[] ParseMSGFDBResults(String FileName) {
		PeptideSpectrumMatch[] AllAnnotations = MSGFDBAnnotation
				.LoadResultsFileLess(FileName, MSGFDBAnnotation.qValue, 0.01);

		// We have to update the spectrum indexes since MSGFDB uses 1-based
		// numbering for the index
		// and SPS uses 0-based numbering
		for (int i = 0; i < AllAnnotations.length; ++i) {
			int[] indexes = AllAnnotations[i].getSpectrumIndex();
			for (int j = 0; j < indexes.length; ++j)
				indexes[j] -= 1;
			AllAnnotations[i].setSpectrumIndex(indexes);
		}

		return AllAnnotations;
	}

	private boolean runMSGFDB(boolean runInspectFlag) {

		System.out
				.println("\n[GenoMS Part 2]: Run Spectrum Identification in normal mode");
		/**
		 * If SLGF was run before GenoMS, then we need to do the following 1)
		 * Load the SLGF results 2) Generate a set of unannotated spectra to run
		 * the database search
		 */
		String[] searchSpectraFileNames = new String[this.SpectraFileNames.length];
		String[] mapFileNames = new String[this.SpectraFileNames.length];
		PeptideSpectrumMatch[] allSpecLibAnnotations = null;
		if (this.runLibrarySearch) {
			allSpecLibAnnotations = SpectrumLibraryAnnotation
					.LoadResultsFileLess(this.specLibSearchResultsFile,
							SpectrumLibraryAnnotation.qValue, 0.01);

			for (int i = 0; i < allSpecLibAnnotations.length; ++i) {
				int[] sNums = allSpecLibAnnotations[i].getScanNumber();
				for (int j = 0; j < sNums.length; ++j)
					sNums[j] -= 1;
				allSpecLibAnnotations[i].setSpectrumIndex(sNums);
				// allSpecLibAnnotations[i].debugPrint();

			}
			for (int i = 0; i < this.SpectraFileNames.length; ++i) {
				searchSpectraFileNames[i] = Utils
						.GetFileNameNoExtension(this.SpectraFileNames[i])
						+ ".unannotated"
						+ Utils.GetFileExtension(this.SpectraFileNames[i]);
				String mFileName = SpectrumLibraryAnnotation
						.writeUnIdentifiedSpectra(allSpecLibAnnotations,
								this.SpectraFileNames[i],
								searchSpectraFileNames[i]);

				mapFileNames[i] = mFileName;
			}
			System.out.println("Loaded " + allSpecLibAnnotations.length
					+ " annotations from spectral library search");
		}

		else {
			// If we didn't do a library search, then all spectra are searched
			// by the DBsearch
			for (int i = 0; i < searchSpectraFileNames.length; ++i) {
				searchSpectraFileNames[i] = this.SpectraFileNames[i];
				mapFileNames[i] = null;
			}
		}

		if (this.SpectraFileNames.length > 1) {
			ErrorThrower
					.ThrowWarningCustom(ErrorThrower.CUSTOM_WARNING,
							"More than one spectrum file was specified! We are only expecting one!");
		}

		this.dbSearchOutputFiles = new String[this.SpectraFileNames.length];

		// for (int i = 0; i < this.SpectraFileNames.length; ++i) {
		for (int i = 0; i < searchSpectraFileNames.length; ++i) {
			this.dbSearchOutputFiles[i] = Utils
					.GetFileNameNoExtension(searchSpectraFileNames[i]) + ".tsv";
			File Test = new File(this.dbSearchOutputFiles[i]);
			if (runInspectFlag && Test.exists())
				Utils.DeleteFile(this.dbSearchOutputFiles[i]);
			if (runInspectFlag || !Test.exists()) {

				ArrayList<String> tempList = new ArrayList<String>();

				tempList.add("-s");
				tempList.add(searchSpectraFileNames[i]);

				tempList.add("-d");

				tempList.add(this.FASTADatabase);

				tempList.add("-t");
				tempList.add(this.PMTolerance + "Da");

				tempList.add("-thread");
				tempList.add("1");

				tempList.add("-o");
				tempList.add(dbSearchOutputFiles[i]);

				tempList.add("-inst");
				if (this.inst.equals("LTQ")) {
					tempList.add("0");
				} else if (this.inst.equals("FT")) {
					tempList.add("1");
				} else if (this.inst.equals("TOF")) {
					tempList.add("2");
				}

				tempList.add("-e");
				// Figure out the digest
				if (this.Digest.toLowerCase().compareTo("trypsin") == 0)
					tempList.add("1");
				else if (this.Digest.toLowerCase().compareTo("chymotrypsin") == 0)
					tempList.add("2");
				else
					tempList.add("0");

				// Figure out the fixed mods, yo
				if (this.FixedMods == null) {
					tempList.add("-mod");
					tempList.add(AntibodyUtils.makePathAbsolute(
							this.ResourceDir, "AASetStandard.txt"));
				} else {
					int CDelta = (int) (this.FixedMods['C' - 65]);
					if (CDelta == 99) {
						tempList.add("-mod");
						tempList.add(AntibodyUtils.makePathAbsolute(
								"AASetC99.txt", this.ResourceDir));
					} else if (CDelta == 58) {
						tempList.add("-mod");
						tempList.add(AntibodyUtils.makePathAbsolute(
								"AASetC58.txt", this.ResourceDir));
						
					} else if (CDelta == 105) {
						tempList.add("-mod");
						tempList.add(AntibodyUtils.makePathAbsolute(
								"AASetC105.txt", this.ResourceDir));
					
					}else if (CDelta != 57) {
						tempList.add("-mod");
						tempList.add(AntibodyUtils.makePathAbsolute(
								"AASetStandard.txt", this.ResourceDir));
					}
				}
				tempList.add("-tda");
				tempList.add("1");

				// Add the fragmentation type
				tempList.add("-m");
				if (this.actType == basicUtils.Utils.activationType.CID)
					tempList.add("1");
				else if (this.actType == basicUtils.Utils.activationType.ETD)
					tempList.add("2");
				else if (this.actType == basicUtils.Utils.activationType.HCD)
					tempList.add("3");
				else if (this.actType == basicUtils.Utils.activationType.MultiPaired)
					tempList.add("4");
				else if (this.actType == basicUtils.Utils.activationType.Unpaired)
					tempList.add("0");

				String[] argv = Utils.ConvertArraylistToStringArray(tempList);
				System.out.println("MSGFDB args: "
						+ Utils.JoinStringArray(argv, " "));
				try {
					MSGFDB.main(argv);
				} catch (Exception e) {
					e.printStackTrace();
					ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
							"Unable to run MSGFDB!");
				}

			}
		}
		// This works because we only allow MSGFDB to analyze one file
		PeptideSpectrumMatch[] AllAnnotations = this
				.ParseMSGFDBResults(this.dbSearchOutputFiles[0]);
		if (AllAnnotations == null) {
			ErrorThrower.ThrowError(10, this.dbSearchOutputFiles[0]);

		}
		System.out.println("Loaded " + AllAnnotations.length
				+ " MSGFDB annotations");
		this.annCount = AllAnnotations.length;

		// If we ran the Library search, then we need to update the
		// spectrumIndexes for the DBsearch results
		if (this.runLibrarySearch) {
			PeptideSpectrumMatch[] combinedAnns = new PeptideSpectrumMatch[AllAnnotations.length
					+ allSpecLibAnnotations.length];
			for (int i = 0; i < searchSpectraFileNames.length; ++i) {
				AllAnnotations = SpectrumLibraryAnnotation
						.updateAnnotationIndexes(mapFileNames[i],
								AllAnnotations);

			}
			for (int i = 0; i < AllAnnotations.length; ++i)
				combinedAnns[i] = AllAnnotations[i];
			// Augment Annotations with the speclib results
			for (int i = 0; i < allSpecLibAnnotations.length; ++i) {
				combinedAnns[i + AllAnnotations.length] = allSpecLibAnnotations[i];
			}
			AllAnnotations = combinedAnns;
		}

		// Create the map from MSGFDB results to PRMSpectra
		ArrayList<Integer> currScansIndexes = new ArrayList<Integer>();
		ArrayList<String> currPeps = new ArrayList<String>();

		for (int i = 0; i < this.allPRMSpectra.length; ++i) {
			int prmIndex = this.allPRMSpectra[i].GetSpecIndex()[0];
			boolean localDebug = false;

			// What raw spectral scans match this PRM scan?
			if (this.unRavelPRMCluster) {
				int[] specInfo = this.PRM2SpectrumMap
						.get(new Integer(prmIndex));
				// currScans = new int[specInfo.length];
				currScansIndexes.clear();
				currPeps.clear();
				// currPeps = new String[specInfo.length];
				for (int j = 0; j < specInfo.length; ++j) {
					currScansIndexes.add(new Integer(specInfo[j]));
					// if(specInfo[j] == 9721)
					// localDebug = true;
					currPeps.add(null);
				}
			}
			// If we are not unraveling which MS/MS indexes went with what PRM
			// spectrum, then its just one to one.
			else {
				// currScans = new int[this.numConsec];
				// currPeps = new String[this.numConsec];
				currScansIndexes.clear();
				currPeps.clear();
				currScansIndexes.add(new Integer(prmIndex));
				currPeps.add(null);
			}
			if (localDebug) {
				System.out.println("Found PRM spectrum " + i
						+ " contains our raw spectrum 9721");
				for (int j = 0; j < currScansIndexes.size(); ++j)
					System.out.println(" - " + currScansIndexes.get(j));
			}
			// Any ID's for these raw spectral scans?
			for (int j = 0; j < AllAnnotations.length; ++j) {
				int[] rawScanIndexes = AllAnnotations[j].getSpectrumIndex();
				for (int l = 0; l < currScansIndexes.size(); ++l) {
					for (int k = 0; k < rawScanIndexes.length; ++k) {
						if (currScansIndexes.get(l).intValue() == rawScanIndexes[k]) {
							currPeps.set(l,
									AllAnnotations[j].getModifiedPeptide());
						}
					}
				}
			}

			// Determine a consensus for the SpectrumNode
			String bestVote = AntibodyUtils.getMostPopular(currPeps);

			if (localDebug) {
				System.out.println("Peptide chosen: " + bestVote);
				for (int j = 0; j < currPeps.size(); ++j)
					System.out.println(" - " + currPeps.get(j));

			}
			// System.out.println(Utils.JoinStringArray(currPeps, ","));
			// System.out.println("Best: " + bestVote);
			if (bestVote != null)
				this.allPRMSpectra[i].setDBAnnotation(bestVote);
		}
		return true;
	}

	public void DeterminePValue(String RawOutputFile, String FinalOutputFile,
			double PValueCutoff, int ScoreMode) {
		FDRCalculator.performFDRCutoffWithDecoys(RawOutputFile, FinalOutputFile,
				PValueCutoff, ScoreMode);
	}

	/**
	 * We only create SpectrumNodes for spectra with PRMs precomputed using
	 * PepNovo. Creates the list of spectrumNodes and the indexes.
	 * 
	 * @return Returns true if load is successful
	 */
	private boolean LoadPRMs() {
		System.out.println("\n[GenoMS Part 1]: Load PRM Spectra");
		PRMSpectrum[] allPRMs = null;

		// An ArrayList in which to create a single list of PRMSpectrum objects
		// Simply concatenates the PRMSpectrum objects loaded from each file.
		ArrayList<PRMSpectrum> tempList = new ArrayList<PRMSpectrum>();

		int totalPRMs = 0; // The total number of PRM peaks in all PRM spectra
		double totalPRMScore = 0; // The total PRM peak scores in all PRM
									// spectra

		for (int i = 0; i < this.PRMFileNames.length; ++i) {
			System.out.println("Loading PRMs from : " + this.PRMFileNames[i]);
			try {
				allPRMs = PRMSpectrum.LoadAllPRMSpectrum(this.PRMFileNames[i]);
			} catch (Exception E) {
				ErrorThrower.ThrowError(10, this.PRMFileNames[i]);
			}

			// Add the new PRMSpecturm objects to our list
			for (int j = 0; j < allPRMs.length; ++j) {

				tempList.add(allPRMs[j]);
			}
		}

		// Convert the ArrayList to a PRMSpectrum array
		allPRMs = new PRMSpectrum[tempList.size()];
		for (int i = 0; i < tempList.size(); ++i) {

			allPRMs[i] = (PRMSpectrum) (tempList.get(i));
			totalPRMs += allPRMs[i].PRMs.length;
			totalPRMScore += allPRMs[i].PRMs.length * allPRMs[i].MeanScore;
		}

		this.spectrumNodeIndex = new Hashtable<String, SpectrumNode>();

		// Create the link between SpectrumNode and PRM Spectrum
		this.allPRMSpectra = new SpectrumNode[allPRMs.length];

		// If PrmCluster was used, then we need to create a mapping between the
		// PRM spectrum index and the
		// raw MS/MS spectrum index. Normally there would be 1 PRM spectrum per
		// MS/MS spectrum
		if (this.unRavelPRMCluster) {
			this.PRM2SpectrumMap = this.loadPRM2SpectrumMapping();
		}

		for (int i = 0; i < this.allPRMSpectra.length; ++i) {
			this.allPRMSpectra[i] = new SpectrumNode(allPRMs[i]);

			int[] specDetails = null;
			if (this.unRavelPRMCluster)
				specDetails = this.PRM2SpectrumMap.get(new Integer(i));

			// Set the number of MS/MS spectra that are merged into this single
			// PRM spectrum
			if (this.unRavelPRMCluster)
				this.allPRMSpectra[i].setNumSpectraMerged(specDetails.length);

			// ASSUMPTION: If we are not unraveling PRMs, then the number of
			// consecutive ids merged is 1
			else
				this.allPRMSpectra[i].setNumSpectraMerged(1);

			String FileName = Utils.GetBaseNameNoExtension(allPRMs[i].FileName)
					+ ".mgf";

			// Create the key to lookup the SpectrumNode in the
			// spectrumNodeIndex
			String hashKey = AntibodyUtils.CreateSpectrumKey(FileName,
					allPRMs[i].SpecIndex);

			// If we are writing to reports, then we have a single input file,
			// so we only need
			// the PRM spectrum index as a key.
			if (this.reportFileDir != null)
				hashKey = allPRMs[i].SpecIndex[0] + "";

			// Check for a collision in this key (it should be unique)
			if (this.spectrumNodeIndex.containsKey(hashKey)) {
				ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
						"LoadPRMs found two PRMs for " + hashKey);
			}

			// Add the key and SpectrumNode to our index
			this.spectrumNodeIndex.put(hashKey, this.allPRMSpectra[i]);
		}
		System.out.println("Loaded " + this.allPRMSpectra.length
				+ " PRMSpectra");

		if (Debug) {

			System.out.println("Mean PRMScore: " + (totalPRMScore / totalPRMs));
		}

		// In our alignment of PRM spectra to the HMM, we penalize for missing
		// peaks
		// by subtracting the mean PRM score.
		if (this.applyMatchPenalty)
			this.matchPenalty = (totalPRMScore / totalPRMs);

		return true;
	}

	/**
	 * When PrmClust is the clustering tool, we must unravel which raw input
	 * spectra are merged together. Looks in spectra/specs_scored.clust. Loads
	 * the spectrum indexes. These are 0-based!
	 * 
	 * @return
	 */
	private Hashtable<Integer, int[]> loadPRM2SpectrumMapping() {

		boolean localDebug = false;

		DataInputStream inFile = null;
		/*
		 * String clusterFileName = this.reportFileDir + File.separator +
		 * "spectra" + File.separator + "clusterData.bin"; String mappingFile =
		 * this.reportFileDir + File.separator + "spectra" + File.separator +
		 * "input_mapping.bin";
		 */
		String mapFileName = this.projectDir + File.separator + "spectra"
				+ File.separator + "specs_scored.clust";

		// Maps the PRM spectrum index to the original spectrum index
		Hashtable<Integer, int[]> ret = new Hashtable<Integer, int[]>();

		try {
			inFile = new DataInputStream(new FileInputStream(mapFileName));
		} catch (Exception E) {
			ErrorThrower.ThrowError(5, mapFileName);
		}

		// First read the version info
		try {

			// Read the strings
			int numStrings = Integer.reverseBytes(inFile.readInt());
			String[] keys = new String[numStrings];
			short[] values = new short[numStrings];

			// First read the string keys (len:int key:string)
			for (int i = 0; i < numStrings; ++i) {
				int strLen = Integer.reverseBytes(inFile.readInt());
				if (strLen > 0) {
					byte[] b = new byte[strLen];
					inFile.read(b);
					keys[i] = new String(b);
				} else
					keys[i] = "";
			}

			// Read the values (short)
			for (int i = 0; i < numStrings; ++i) {
				values[i] = Short.reverseBytes(inFile.readShort());
				if (localDebug)
					System.out.println(keys[i] + " = " + values[i]);
			}
		} catch (IOException E) {
			ErrorThrower.ThrowError(6, mapFileName);
		}

		// Now read the clusters
		try {
			int numClusters = Integer.reverseBytes(inFile.readInt());

			if (localDebug)
				System.out.println("LOCALDEBUG: Total clusters = "
						+ numClusters);
			// Now we read in each cluster
			for (int clusterIndex = 0; clusterIndex < numClusters; ++clusterIndex) {
				int cScan = Integer.reverseBytes(inFile.readInt());
				int cIndex = Integer.reverseBytes(inFile.readInt());
				int fileIndex = Integer.reverseBytes(inFile.readInt());

				if (localDebug)
					System.out.println("LOCALDEBUG: scan=" + cScan + ", index="
							+ cIndex + ", fileIndex=" + fileIndex);

				// Read in the fileNames
				int numStrings = Integer.reverseBytes(inFile.readInt());
				String[] clusterFileName = new String[numStrings];
				for (int i = 0; i < numStrings; ++i) {
					int strLen = Integer.reverseBytes(inFile.readInt());
					if (strLen > 0) {
						byte[] b = new byte[strLen];
						inFile.read(b);
						clusterFileName[i] = new String(b);
					} else
						clusterFileName[i] = "";
				}
				if (localDebug)
					System.out.println("LOCALDEBUG: file name="
							+ clusterFileName[0]);
				int numSpectra = Integer.reverseBytes(inFile.readInt());

				if (localDebug)
					System.out.println("LOCALDEBUG: spectra=" + numSpectra);
				int[] spectraIndexes = new int[numSpectra];
				// int[] fileIndexes = new int[numSpectra];

				// Read scanNumbers (for PRMClust these are set to the cluster
				// index)
				for (int i = 0; i < numSpectra; ++i)
					Integer.reverseBytes(inFile.readInt());

				// Read indexes
				for (int i = 0; i < numSpectra; ++i)
					spectraIndexes[i] = Integer.reverseBytes(inFile.readInt());

				// Read file indexes (since this is from PRMClust, there should
				// only be one input file)
				for (int i = 0; i < numSpectra; ++i)
					Integer.reverseBytes(inFile.readInt());

				// Read spectrum file names
				numStrings = Integer.reverseBytes(inFile.readInt());
				String[] specFileNames = new String[numStrings];
				for (int i = 0; i < numStrings; ++i) {
					int strLen = Integer.reverseBytes(inFile.readInt());
					if (strLen > 0) {
						byte[] b = new byte[strLen];
						inFile.read(b);
						specFileNames[i] = new String(b);
					} else
						specFileNames[i] = "";
				}

				// Put this guy into our map
				if (localDebug) {
					System.out.println("LOCALDEBUG: Cluster index " + cIndex
							+ " maps to ");
					for (int i = 0; i < spectraIndexes.length; ++i) {
						System.out.print(spectraIndexes[i] + " ");
					}
					System.out.print("\n");
				}
				ret.put(new Integer(cIndex), spectraIndexes);
			}
		} catch (IOException E) {
			ErrorThrower.ThrowError(6, mapFileName);
		}
		try {
			inFile.close();
		} catch (Exception E) {
			ErrorThrower.ThrowError(8, mapFileName);
		}
		return ret;

		// load the input_mapping (maps the global MS/MS spectrum index to the
		// file number and file spectrum index)
		/*
		 * try{ inFile = new DataInputStream(new FileInputStream(mappingFile));
		 * } catch(Exception E) { ErrorThrower.ThrowError(5, mappingFile); }
		 * 
		 * int numRows = -1; int numCols = -1; try { numRows =
		 * Integer.reverseBytes(inFile.readInt()); numCols =
		 * Integer.reverseBytes(inFile.readInt());
		 * 
		 * } catch(Exception E) { ErrorThrower.ThrowError(6, mappingFile); }
		 * 
		 * if(numCols != 2) {
		 * ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
		 * "It seems that input_mapping.bin is not in the correct 2-column format"
		 * ); } if(localDebug) System.out.println("Loading " + numRows +
		 * " total raw spectra!"); Hashtable<Integer,Integer> map = new
		 * Hashtable<Integer,Integer>();
		 * 
		 * for(int i = 0; i < numRows; ++i) { int fileNum = -1; int SpecNum =
		 * -1;
		 * 
		 * try { fileNum = Integer.reverseBytes(inFile.readInt()); SpecNum =
		 * Integer.reverseBytes(inFile.readInt());
		 * 
		 * } catch(Exception E) { ErrorThrower.ThrowError(6, mappingFile); }
		 * 
		 * 
		 * map.put(new Integer(fileNum*1000000+SpecNum), new Integer(i+1)); }
		 * 
		 * try{ inFile.close(); } catch(Exception E) {
		 * ErrorThrower.ThrowError(8, mappingFile); }
		 * 
		 * //Load the clusterData (maps the PRM spectra to the MS/MS file number
		 * and file indexes) try{ inFile = new DataInputStream(new
		 * FileInputStream(clusterFileName)); } catch(Exception E) {
		 * ErrorThrower.ThrowError(5, clusterFileName); }
		 * 
		 * int numPRMs = 0; try { numPRMs =
		 * Integer.reverseBytes(inFile.readInt());
		 * 
		 * } catch(Exception E) { ErrorThrower.ThrowError(6, clusterFileName); }
		 * if(numPRMs != this.allPRMSpectra.length) {
		 * ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,"Num " +
		 * "PRM spectra listed in spectra/clusterData.bin (" + numPRMs +
		 * ") does " + "not match the number loaded (" +
		 * this.allPRMSpectra.length + ")!"); }
		 * 
		 * Hashtable<Integer,int[]> ret = new Hashtable<Integer,int[]>();
		 * 
		 * for(int i = 0; i < numPRMs; ++i) { int prmNum = -1; int numSpecs =
		 * -1;
		 * 
		 * try { prmNum = Integer.reverseBytes(inFile.readInt()); numSpecs =
		 * Integer.reverseBytes(inFile.readInt());
		 * 
		 * } catch(Exception E) { ErrorThrower.ThrowError(6, clusterFileName); }
		 * if(localDebug) System.out.println(" " + i + " - " + prmNum + " has "
		 * + numSpecs + " spectra");
		 * 
		 * int[] specIdxs = new int[numSpecs]; int fileIdx = -1; int specIdx =
		 * -1;
		 * 
		 * for(int j = 0; j < numSpecs; ++j) { try { fileIdx =
		 * Integer.reverseBytes(inFile.readInt()); specIdx =
		 * Integer.reverseBytes(inFile.readInt());
		 * 
		 * } catch(Exception E) { ErrorThrower.ThrowError(6, clusterFileName); }
		 * specIdxs[j] = (map.get(new
		 * Integer(fileIdx*1000000+specIdx))).intValue(); if(localDebug)
		 * System.out.println("  " + j + " - file " + fileIdx + " spectrum " +
		 * specIdx + "=" + specIdxs[j]); }
		 * 
		 * ret.put(new Integer(prmNum), specIdxs); }
		 * 
		 * try{ inFile.close(); } catch(Exception E) {
		 * ErrorThrower.ThrowError(8, clusterFileName); } return ret;
		 */

	}

	/**
	 * Driver method for anchor identification and template selection. Performs
	 * the following: 1. Greedy assignment of peptides to templates 2. Scores
	 * templates 3. Greedily finds set of templates 4. Creates anchors on
	 * selected templates (populates this.SelectedTemplates)
	 * 
	 * @return
	 */
	private boolean createSequenceSeeds() {
		System.out.println("\n[GenoMS Part 4]:Create sequence seeds");

		boolean LocalDebug = false; // [ProteinID]

		// We use the TrieDB utilities to probe the DB, so convert the FASTA
		// Database to a trie database
		// if we don't already have a trie database
		if (this.trieDatabase == null) {
			String trieDB = Utils.GetFileNameNoExtension(this.FASTADatabase)
					+ ".trie";
			TrieDB.prepDB(this.FASTADatabase, trieDB);
			this.trieDatabase = new TrieDB(trieDB);

		}

		/*
		 * Load all DB search Annotations, populating ProteinIDs and Peptides
		 * and Peptide2ProteinIDs
		 */
		System.out.println("Loading db search results...");
		Object[] deets = AntibodyUtils.createSeedsFromDBSearchResults(
				this.dbPepObj, this.trieDatabase, false, LocalDebug);

		// TemplateID/ProteinID -> ArrayList of Peptide sequences
		Hashtable<Integer, ArrayList<String>> proteinIDs2PepSeq = (Hashtable<Integer, ArrayList<String>>) (deets[0]);

		// Peptide -> ArrayList of Peptide Object indexes
		Hashtable<String, ArrayList<Integer>> pepSeq2PepObj = (Hashtable) (deets[1]);

		// Peptide -> ArrayList of protein IDs
		Hashtable<String, ArrayList<Object[]>> pepSeq2ProteinIDs = (Hashtable) (deets[2]);

		/*
		 * Score each template according to coverage. Forbid templates that
		 * share AntibodyUtils.MAX_SHARED_PEPTIDES or more peptides
		 */
		if (this.allTemplates == null) {
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"No templates loaded when trying to create seeds");

		}

		System.out.println("Determining template scores...");

		// Maps the ProteinID to the number of peptides shared between current
		// template and the ProteinID template.
		Hashtable<Integer, Integer> RedundantPeptides = new Hashtable<Integer, Integer>();

		for (int i = 0; i < this.allTemplates.length; ++i) {

			// Skip thsis template if we have no peptides matching it
			if (!proteinIDs2PepSeq.containsKey(new Integer(
					this.allTemplates[i].ProteinID)))
				continue;

			System.out.println("Scoring " + i + "/" + this.allTemplates.length);
			if (LocalDebug)
				System.out.println("Considering template " + i + ":"
						+ this.allTemplates[i].TemplateName + ":"
						+ this.allTemplates[i].ProteinID);
			RedundantPeptides.clear();

			// ArrayList of Peptide objects supporting this template
			ArrayList<Peptide> peptideObjs = new ArrayList<Peptide>();

			// ArrayList of peptide sequences supporting this template
			ArrayList<String> peptideList = (ArrayList<String>) (proteinIDs2PepSeq
					.get(new Integer(this.allTemplates[i].ProteinID)));

			if (LocalDebug)
				System.out.println("  contains " + peptideList + " peptides");

			// Iterate over each peptide matching this template, and add the
			// Peptides to the list 'peptideObjs'
			for (int p = 0; p < peptideList.size(); ++p) {
				String peptideSeq = (String) (peptideList.get(p));

				// Add those Peptides Objects to the template
				ArrayList<Integer> pepObjIndexes = (ArrayList<Integer>) pepSeq2PepObj
						.get(peptideSeq);
				for (int j = 0; j < pepObjIndexes.size(); ++j) {
					int Index = ((Integer) (pepObjIndexes.get(j))).intValue();
					Peptide obj = (antibody.Peptide) this.dbPepObj.get(Index);
					peptideObjs.add(obj);
				}

				// Count the other Proteins that share this peptide
				ArrayList<Object[]> OtherProteins = (ArrayList<Object[]>) (pepSeq2ProteinIDs
						.get(peptideSeq));
				if (LocalDebug)
					System.out.println("  '" + peptideSeq + "' has "
							+ OtherProteins.size() + " proteins");

				for (int j = 0; j < OtherProteins.size(); ++j) {
					Object[] otherLocation = ((Object[]) (OtherProteins.get(j)));

					int otherProteinID = ((Integer) (otherLocation[TrieDB.TrieLocationColumns.ProteinID]))
							.intValue();
					if (otherProteinID == this.allTemplates[i].ProteinID)
						continue;

					if (!RedundantPeptides.containsKey(new Integer(
							otherProteinID))) {
						RedundantPeptides.put(new Integer(otherProteinID),
								new Integer(0));
					}
					int Count = ((Integer) (RedundantPeptides.get(new Integer(
							otherProteinID)))).intValue();
					RedundantPeptides.put(new Integer(otherProteinID),
							new Integer(Count + 1));
				}
			}

			// Create anchors on that template
			System.out.println("Creating anchors on template with "
					+ peptideObjs.size() + " peptides ...");
			Template.createAnchorsOnTemplate(this.allTemplates,
					this.allTemplates[i].ProteinID, peptideObjs, false,
					this.trieDatabase);

			Enumeration<Integer> CompetingTemplates = RedundantPeptides.keys();

			System.out.println("Determining competing templates...");
			while (CompetingTemplates.hasMoreElements()) {

				int offendingID = ((Integer) (CompetingTemplates.nextElement()))
						.intValue();
				int pepCount = ((Integer) (RedundantPeptides.get(new Integer(
						offendingID)))).intValue();
				if (pepCount > AntibodyUtils.MAX_SHARED_PEPTIDES
						|| ((double) pepCount) / (peptideList.size()) >= 0.5) {
					for (int k = 0; k < this.allTemplates.length; ++k) {
						if (this.allTemplates[k].ProteinID == offendingID) {
							if (this.allTemplates[k].TemplateName == null) {
								{
									this.allTemplates[k].TemplateName = this.trieDatabase
											.getProteinName(this.allTemplates[k].ProteinID);

								}
							}
							this.allTemplates[i]
									.AddForbiddenEdge(this.allTemplates[k]);
							break;
						}
					}
				}
			}

			// if (LocalDebug)
			// AntibodyUtils.//();
		}

		System.out.println("\n[GenoMS Part 5]: Select best templates");
		this.SelectedTemplates = Template.chooseBestTemplateSetMultiPath(
				this.allTemplates, this.MinPeptidesPerPath);
		// if (LocalDebug || this.Debug)
		// this.DebugPrintSelectedTemplates();
		// AntibodyUtils.WaitForEnter();
		if (this.SelectedTemplates != null)
			LearnQFromInspect();
		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			Template t = this.SelectedTemplates[i];
			System.out.print("[" + i + "]: ");
			while (t != null) {
				System.out.print(t.TemplateName + " -> ");
				t = t.NextInList;
			}
			System.out.println("");
		}
		return (this.SelectedTemplates != null);
	}

	public void DebugPrintSelectedTemplates() {
		System.out.println("Selected " + this.SelectedTemplates.length
				+ " paths from the graphs");
		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			System.out.println("Template list [" + i + "]");
			Template Ptr = this.SelectedTemplates[i];
			while (Ptr != null) {
				if (Ptr.Sequence != null)
					Ptr.DebugPrint();
				Ptr = Ptr.NextInList;
			}
		}
		// AntibodyUtils.WaitForEnter();

	}

	/**
	 * Once we have selected templates and created anchors, determine the
	 * average, min, and max alignment scores for some guys
	 */
	private void LearnQFromInspect() {
		boolean LocalDebug = false;
		if (this.SelectedTemplates == null)
			return;

		System.out
				.println("\n[GenoMS Part 6]: Learn spectrum alignment parameters");
		double TrueSpectraScores = 0.0;
		int TrueSpectraCount = 0;
		double FalseSpectraScores = 0.0;
		int FalseSpectraCount = 0;

		double MaxTrueScore = 0.0;
		double MinTrueScore = Double.MAX_VALUE;
		double MaxFalseScore = 0.0;
		double MinFalseScore = Double.MAX_VALUE;

		// double[] TrueCounts = new double[200];
		// double[] FalseCounts = new double[200];
		ArrayList<Double> TrueCounts = new ArrayList<Double>();
		ArrayList<Double> FalseCounts = new ArrayList<Double>();
		int[] FalseDistro = new int[200];
		int[] TrueDistro = new int[200];

		// Iterate over each template set (path)
		for (int i = 0; i < this.SelectedTemplates.length; ++i) {

			Template Ptr = this.SelectedTemplates[i];

			// Iterate over all templates in the set
			while (Ptr != null) {

				// Skip templates that don't have any anchors
				if (Ptr.anchors == null) {
					Ptr = Ptr.NextInList;
					continue;
				}

				Seed CurrAnchor = Ptr.anchors;
				while (CurrAnchor != null) {
					if (LocalDebug)
						System.out.println("Considering Anchor: "
								+ CurrAnchor.AnchoredSeedSequence);
					ArrayList<SpectrumNode> SpectraList = CurrAnchor.SupportingSpectra;
					int[][] Starts = new int[SpectraList.size()][2];
					SpectrumNode[] CurrSpectra = new SpectrumNode[Starts.length];

					for (int j = 0; j < SpectraList.size(); ++j) {
						SpectrumNode CurrAnn = (SpectrumNode) (SpectraList
								.get(j));
						String CurrUnModdedAnn = CurrAnn
								.getUnModdedDBAnnotation();
						// If this node is supporting an anchor, it should have
						// a seqence!
						if (CurrUnModdedAnn == null) {
							ErrorThrower
									.ThrowErrorCustum(
											ErrorThrower.CUSTOM_ERROR,
											"Spectrum ("
													+ AntibodyUtils.CreateSpectrumKey(
															CurrAnn.GetSourceFileName(),
															CurrAnn.GetSpecIndex())
													+ ") is attached to anchor '"
													+ CurrAnchor.AnchoredSeedSequence
													+ "' but has no seq!");
						}

						// Find a spectrum node for this annotation
						CurrSpectra[j] = CurrAnn;
						Starts[j][0] = CurrAnchor.AnchoredSeedSequence
								.indexOf(CurrUnModdedAnn);
						Starts[j][1] = Starts[j][0] + CurrUnModdedAnn.length();
					}

					for (int j = AntibodyUtils.MIN_OVERLAPPING_PEAKS; j < CurrAnchor.AnchoredSeedSequence
							.length(); ++j) {
						if (LocalDebug)
							System.out.println("Considering anchor ending at "
									+ j);

						// Get the spectra which overlap this anchor as the
						// 'true' set
						for (int k = 0; k < Starts.length; ++k) {
							if (Starts[k][0] < 0)
								continue;
							if (Starts[k][0] <= j
									- AntibodyUtils.MIN_OVERLAPPING_PEAKS + 1
									&& Starts[k][1] >= j) // This guy overlaps
							// and is a
							// candidate
							{
								if (LocalDebug) {
									System.out
											.println("Found a spectrum with overlaps this, starts at "
													+ Starts[k][0]
													+ " and ends at "
													+ Starts[k][1]);
									CurrSpectra[j].DebugPrint();
								}

								String TestSeq = CurrAnchor.AnchoredSeedSequence
										.substring(Math.max(0, j - 10), j);
								double CurrScore = AntibodyUtils
										.alignmentScoreSpecAlign(
												TestSeq,
												CurrSpectra[k].GetPRMSpectrum(),
												this.PRMTolerance);
								if (LocalDebug) {
									System.out.println("Score: " + CurrScore);
									Utils.WaitForEnter();
								}
								if (CurrScore < 0)
									continue;
								if ((int) (CurrScore) >= TrueDistro.length)
									CurrScore = TrueDistro.length - 1;
								TrueSpectraCount += 1;
								TrueCounts.add(new Double(CurrScore));
								TrueDistro[(int) (CurrScore)] += 1;

								if (CurrScore > MaxTrueScore)
									MaxTrueScore = CurrScore;
								if (CurrScore < MinTrueScore)
									MinTrueScore = CurrScore;
								TrueSpectraScores += CurrScore;
							} else if (j < Starts[k][0] || j > Starts[k][1] + 4) {

								String TestSeq = CurrAnchor.AnchoredSeedSequence
										.substring(Math.max(0, j - 10), j);
								double Score = AntibodyUtils
										.alignmentScoreSpecAlign(
												TestSeq,
												CurrSpectra[k].GetPRMSpectrum(),
												this.PRMTolerance);
								if (LocalDebug) {
									System.out
											.println("Found a spectrum with does NOT overlaps this, starts at "
													+ Starts[k][0]
													+ " and ends at "
													+ Starts[k][1]);
									CurrSpectra[k].DebugPrint();
									System.out.println("False Score: " + Score);
									Utils.WaitForEnter();
								}
								if (Score < 0)
									continue;
								FalseSpectraCount += 1;
								FalseCounts.add(new Double(Math.min(
										FalseDistro.length - 1, Score)));
								if (Score > MaxFalseScore)
									MaxFalseScore = Score;
								if (Score < MinFalseScore)
									MinFalseScore = Score;
								FalseSpectraScores += Score;
								FalseDistro[Math.min(FalseDistro.length - 1,
										(int) (Score))] += 1;
							}
						}

					}
					CurrAnchor = CurrAnchor.NextSeed;

				}
				Ptr = Ptr.NextInList;
			}

		}
		if (TrueCounts.size() == 0)
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Unable to load any true alignments!");
		if (FalseCounts.size() == 0)
			ErrorThrower.ThrowWarningCustom(ErrorThrower.CUSTOM_WARNING,
					"Unable to load any false alignments!");

		double MedianTrue = Utils.getMedian(Utils
				.ConvertArraylistToDoubleArray(TrueCounts));
		double MedianFalse = Double.NaN;
		if (FalseCounts != null && FalseCounts.size() > 0)
			Utils.getMedian(Utils.ConvertArraylistToDoubleArray(FalseCounts));

		double AverageTrue = TrueSpectraScores / TrueSpectraCount;
		double AverageFalse = 0;
		if (FalseSpectraCount != 0)
			AverageFalse = FalseSpectraScores / FalseSpectraCount;
		else {
			MinFalseScore = Double.NaN;
		}

		AntibodyUtils.MSALIGN_SCORE_CUTOFF = AverageTrue;

		System.out.println("TotalTrueSpectra: " + TrueSpectraCount);
		System.out.println("AverageTrueScore: " + AverageTrue);
		System.out.println("MedianTrueScore: " + MedianTrue);
		System.out.println("MinTrue: " + MinTrueScore);
		System.out.println("MaxTrue: " + MaxTrueScore);
		System.out.println("TotalFalseSpectra: " + FalseSpectraCount);
		System.out.println("AverageFalseScore: " + AverageFalse);
		System.out.println("MedianFalseScore: " + MedianFalse);
		System.out.println("MinFalse: " + MinFalseScore);
		System.out.println("MaxFalse: " + MaxFalseScore);

		double CurrCount = 0;
		double CurrFalse = 0;

		AntibodyUtils.MSALIGN_SCORE_CUTOFF = MinTrueScore;

		// Adjust the false so taht it is roughly equal to the number of true.
		if (FalseSpectraCount > 0) {
			double ValueOfFalse = TrueSpectraCount
					/ ((double) FalseSpectraCount);
			if (LocalDebug || !LocalDebug)
				System.out.println("ValueofFalse: " + ValueOfFalse);
			for (int i = FalseDistro.length - 1; i >= 0; --i) {
				CurrCount += FalseDistro[i] + TrueDistro[i];
				if (CurrCount == 0)
					continue;
				CurrFalse += FalseDistro[i] * ValueOfFalse;
				// TrueCount += TrueDistro[i];

				if (CurrFalse / CurrCount <= .05) {
					// System.out.println("Score " + i + " has FDR "
					// + (CurrFalse / CurrCount));
					AntibodyUtils.MSALIGN_SCORE_CUTOFF = i;
					// break;
				}
			}
		}

		System.out.println("Score cutoff: "
				+ AntibodyUtils.MSALIGN_SCORE_CUTOFF);
		// Utils.WaitForEnter();
		// AntibodyUtils.MSALIGN_SCORE_CUTOFF = 62;
		// AntibodyUtils.MSALIGN_SCORE_CUTOFF = 72;

		// if(LocalDebug || ! LocalDebug)
		// AntibodyUtils.WaitForEnter();

	}

	private MSAlignmentType[] GatherSpectraRight(Seed Sequence,
			SpectrumNode[] Spectra) {

		SpecAlign Gatherer = new SpecAlign(Sequence, Spectra, false,
				this.trieDatabase);
		MSAlignmentType[] Alignments = null;

		Alignments = Gatherer.RunAlignmentRight(this.SelectedTemplates,
				this.PRMTolerance);

		MSAlignmentType[] Ret = MSAlignmentType.GetBestScoreCutOff(
				AntibodyUtils.MSALIGN_SCORE_CUTOFF, Alignments,
				AntibodyUtils.MAX_ALIGNED_SPECTRA);
		// MSAlignmentType[] Ret = MSAlignmentType.GetBestNoCutOff(Alignments,
		// Utils.MAX_ALIGNED_SPECTRA);
		if (Debug) {
			System.out.println("Aligned " + Alignments.length + " spectra to "
					+ Sequence.GetSeedSequence());
			// for(int i = 0; i < Alignments.length; ++i)
			// System.out.println("[" + i + "]: " +
			// Alignments[i].toStringDebug(this.OracleHash));
			System.out.println("Reduced to best " + Ret.length);
			// for(int i = 0; i < Ret.length; ++i)
			// Ret[i].GetSpectrumNode().DebugPrint();
			// AntibodyUtils.WaitForEnter();

		}

		if (Ret.length == 0)
			return null;
		return Ret;
	}

	/**
	 * Does the same thing as GatherSpectr, but attempts to gather spectra that
	 * overlap/extend the Sequence to the left
	 * 
	 * @param Sequence
	 *            The Seed to extend
	 * @param Spectra
	 *            A list of spectra to align to the seed
	 * @return
	 */
	private MSAlignmentType[] GatherSpectraLeft(Seed Sequence,
			SpectrumNode[] Spectra) {

		SpecAlign Gatherer = new SpecAlign(Sequence, Spectra, false,
				this.trieDatabase);
		MSAlignmentType[] Alignments = Gatherer.RunAlignmentLeft(null,
				this.SelectedTemplates, this.PRMTolerance);
		// MSAlignmentType[] Ret =
		// MSAlignmentType.GetBestNonRedundantOdds(Utils.MAX_ALIGNED_SPECTRA,
		// Alignments);
		// MSAlignmentType[] Ret =
		// MSAlignmentType.GetBestOddsCutOff(Utils.MSALIGN_PVALUE_CUTOFF,
		// Alignments, Utils.MAX_ALIGNED_SPECTRA);
		MSAlignmentType[] Ret = MSAlignmentType.GetBestScoreCutOff(
				AntibodyUtils.MSALIGN_SCORE_CUTOFF, Alignments,
				AntibodyUtils.MAX_ALIGNED_SPECTRA);
		// MSAlignmentType[] Ret = MSAlignmentType.GetBestNoCutOff(Alignments,
		// Utils.MAX_ALIGNED_SPECTRA);

		boolean LocalDebug = false;
		if (Debug) {
			System.out.println("Aligned " + Alignments.length + " spectra to "
					+ Sequence.GetSeedSequence());
			System.out.println("Reduced to best " + Ret.length);
			// for(int i = 0; i < Ret.length; ++i)
			// Ret[i].GetSpectrumNode().DebugPrint();
			// AntibodyUtils.WaitForEnter();
			/*
			 * for(int i = 0; i < Ret.length; ++i) { //Ret[i].Print();
			 * Ret[i].PrintPeaks(); }
			 */
		}

		return Ret;
	}

	private ConsensusPRMSpectrum PerformAlignmentWithScoresRight(
			MSAlignmentType[] Spectra, Seed SeedSequence) {
		FileWriter Writer = null;

		GatherConsensus ModelRunner = new GatherConsensus(this.applyMatchPenalty,this.matchPenalty);
		ConsensusPRMSpectrum Consensus = ModelRunner
				.RunGatherAlignWithScoresRight(Spectra, "Temp.txt",
						SeedSequence, null, this.ClusterFlag,
						this.AddMatchStatesAtEndFlag, this.PRMTolerance);

		try {
			Writer = new FileWriter(this.detailsOutputFileName, true);
			Writer.write(ModelRunner.GetModelMatchString());
			Writer.close();
		} catch (Exception E) {
			ErrorThrower.ThrowError(7, this.detailsOutputFileName);
		}

		if (Consensus == null) {

			try {
				Writer = new FileWriter(this.detailsOutputFileName, true);
				Writer.write("ConsensusPRM:**No consensus Produced!!");
				Writer.close();
			} catch (Exception E) {
				ErrorThrower.ThrowError(7, this.detailsOutputFileName);
			}
			return Consensus;
		}
		String PRMs = "";
		String Scores = "";
		for (int i = 0; i < Consensus.ScaledPeaks.length; ++i) {
			PRMs += Consensus.ScaledPeaks[i] + " ";
			Scores += Consensus.PeakScores[i] + " ";
		}

		try {
			Writer = new FileWriter(this.detailsOutputFileName, true);
			Writer.write("ConsensusPRM: " + PRMs);
			Writer.close();
		} catch (Exception E) {
			ErrorThrower.ThrowError(7, this.detailsOutputFileName);
		}
		return Consensus;
	}

	private ConsensusPRMSpectrum PerformAlignmentWithScoresLeft(
			MSAlignmentType[] Spectra, Seed SeedSequence) {
		FileWriter Writer = null;

		GatherConsensus ModelRunner = new GatherConsensus(this.applyMatchPenalty, this.matchPenalty);
		ConsensusPRMSpectrum Consensus = ModelRunner
				.RunGatherAlignWithScoresLeft(Spectra, "Temp.txt",
						SeedSequence, null, this.ClusterFlag,
						this.AddMatchStatesAtEndFlag, this.PRMTolerance);

		try {
			Writer = new FileWriter(this.detailsOutputFileName, true);
			Writer.write(ModelRunner.GetModelMatchString());
			Writer.close();
		} catch (Exception E) {
			ErrorThrower.ThrowError(7, this.detailsOutputFileName);
		}

		if (Consensus == null) {

			try {
				Writer = new FileWriter(this.detailsOutputFileName, true);
				Writer.write("ConsensusPRM:**No consensus Produced!!");
				Writer.close();
			} catch (Exception E) {
				ErrorThrower.ThrowError(7, this.detailsOutputFileName);
			}
			return Consensus;
		}
		String PRMs = "";
		String Scores = "";
		for (int i = 0; i < Consensus.ScaledPeaks.length; ++i) {
			PRMs += Consensus.ScaledPeaks[i] + " ";
			Scores += Consensus.PeakScores[i] + " ";
		}

		try {
			Writer = new FileWriter(this.detailsOutputFileName, true);
			Writer.write("ConsensusPRM: " + PRMs);
			Writer.close();
		} catch (Exception E) {
			ErrorThrower.ThrowError(7, this.detailsOutputFileName);
		}
		return Consensus;
	}

	// Accessors and Setters
	public void SetConfigFileName(String FileName) {
		this.configFileName = FileName;
	}

	public String GetConfigFileName() {
		return this.configFileName;
	}

	public void SetOutputFileName(String FileName) {
		this.outputFileName = FileName;
	}

	public String GetOutputFileName() {
		return this.outputFileName;
	}

	public void SetResourceDir(String FileName) {
		this.ResourceDir = FileName;
	}

	public String GetResourceDir() {
		return this.ResourceDir;
	}

	private void EvaluateDBSearchAndPRMsNoOracle() {
		if (this.allPRMSpectra == null)
			return;

		int TotalSpectra = 0;
		double TotalScore = 0;
		int TotalLength = 0;
		double MaxScore = Double.MIN_VALUE;
		double MinScore = Double.MAX_VALUE;
		int MaxLength = 0;
		int MinLength = Integer.MAX_VALUE;
		int TotalPRMs = 0;
		for (int i = 0; i < this.allPRMSpectra.length; ++i) {
			if (this.allPRMSpectra[i].getModdedDBAnnotation() == null)
				continue;

			SpectrumNode Spectrum = this.allPRMSpectra[i];
			PRMSpectrum CurrSpectrum = Spectrum.GetPRMSpectrum();

			String Peptide = Spectrum.getUnModdedDBAnnotation();

			int[] PRMs = AntibodyUtils.GetSeqPRMScaled(Peptide);
			TotalPRMs += PRMs.length;
			double Score = 0.0;
			int Length = 0;
			for (int k = 0; k < PRMs.length; ++k) {
				for (int l = 0; l < CurrSpectrum.PRMs.length; ++l) {
					if (AntibodyUtils.EquivalentScaled(PRMs[k],
							CurrSpectrum.PRMs[l], this.PRMTolerance)) {
						Score += CurrSpectrum.PRMScores[l];
						Length += 1;
					}
				}
			}
			TotalSpectra += 1;
			TotalScore += Score;
			TotalLength += Length;
			if (Score > MaxScore)
				MaxScore = Score;
			if (Score > 0 && Score < MinScore)
				MinScore = Score;
			if (Length > MaxLength)
				MaxLength = Length;
			if (Length < MinLength)
				MinLength = Length;
		}

		System.out.println("Total PRM spectra with DB annotation: "
				+ TotalSpectra);
		System.out.println("Avg Score: " + TotalScore / TotalSpectra);
		System.out.println("Avg PRMScore: " + TotalScore / TotalLength);
		System.out.println("Max Score: " + MaxScore);
		System.out.println("Min Score: " + MinScore);

		// AntibodyUtils.MIN_PEAK_GROUP_SCORE =
		// AntibodyUtils.MIN_PEAK_GROUP_SCORE;
		AntibodyUtils.MIN_PEAK_GROUP_SCORE = Math.min(MinScore, TotalScore
				/ TotalLength);

		System.out.println("MIN PEAK GROUP SCORE: "
				+ AntibodyUtils.MIN_PEAK_GROUP_SCORE);
		// AntibodyUtils.WaitForEnter();

		System.out.println("Avg Length: " + TotalLength / TotalLength);
		System.out.println("Max Length: " + MaxLength);
		System.out.println("Min Length: " + MinLength);
		// AntibodyUtils.WaitForEnter();

	}

	/**
	 * Testing methods to evaluate extension length and accuracy
	 * 
	 * @throws IOException
	 */

	private int RunDriver() throws IOException {

		if (!this.ParseConfigFile())
			ErrorThrower.ThrowError(10, this.configFileName);
		this.validateParams();

		// Load the PRMs, populating SpectrumNodes
		if (!LoadPRMs())
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Unable to load any PRMs!");

		if (this.DatabaseRootName != null) {
			this.FASTADatabase = FragmentDB
					.createFASTADBandConstraintsFromRootName(this.projectDir,
							this.DatabaseRootName, this.Debug);
			if (this.trieDatabase == null) {
				String trieDB = Utils
						.GetFileNameNoExtension(this.FASTADatabase) + ".trie";
				TrieDB.prepDB(this.FASTADatabase, trieDB);
				this.trieDatabase = new TrieDB(trieDB);
			}
			// this.CombinedDBName = this.FASTADatabase;
			this.DatabaseRootName = null;
			//this.TemplateConstraintFile = Utils
				//	.GetFileNameNoExtension(this.FASTADatabase)
					//+ "_Constraints.bin";
			this.TemplateConstraintFile = FragmentDB.getConstraintFileName(this.trieDatabase.GetDBFileName());
			
			// System.out.println("FASTADATABASE: " + this.FASTADatabase);
			this.allTemplates = Template.CreateConstraintsFromDB(
					this.TemplateConstraintFile,
					this.trieDatabase.GetDBFileName());

			if (this.allTemplates == null) {
				ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
						"Unable to create constraint file from DBRoots!");
			}
			// Template.WriteConstraintsToFile(this.AllTemplates,
			// this.TemplateConstraintFile);
			FragmentDB.WriteConstraintsToFileBinary(this.allTemplates,
					this.TemplateConstraintFile);

		}

		if (!runDBSearchTool(this.runDBSearchFlag)) {
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Unable to run database search");
		}

		// if(this.OracleHash != null)
		// this.EvaluateInspectAndPRMs();
		// else
		this.EvaluateDBSearchAndPRMsNoOracle();
		/*
		 * Whata if we added in the PRMs that supported this peptide!!
		 */
		this.AugmentPRMSpectraForInspectResults();

		System.out.println("\n[GenoMS Part 3]: Load templates and constraints");
		if (this.allTemplates == null) {
			this.allTemplates = FragmentDB
					.LoadAllTemplateFromConstraints(this.TemplateConstraintFile);
			// this.AllTemplates = Template
			// .LoadAllTemplateConstraints(this.TemplateConstraintFile);
		}
		this.mergeSpectraFromSamePeptide();
		this.createSequenceSeeds();

		// for(int i = 0; i < this.SeedSequences.length; ++i)
		// this.PrintSeeds(i);

		// Utils.WaitForEnter();

		this.WriteOutputHeaderInfo(this.detailsOutputFileName);
		this.BuildAndExtendSeeds(this.detailsOutputFileName, 50);

		// System.out.println("Final seeds: ");
		/*
		 * for (int i = 0; i < this.SelectedTemplates.length; ++i) { Template
		 * Ptr = this.SelectedTemplates[i]; while (Ptr != null) {
		 * Ptr.DebugPrint(); Ptr = Ptr.NextInList; } }
		 */
		if (this.reportFileDir == null)
			this.WriteSeedInfoToFile(this.detailsOutputFileName);
		System.out
				.println("\n[GenoMS Part 7]: Reconstruct full protein sequences");

		FileWriter Writer = null;
		try {
			Writer = new FileWriter(this.outputFileName, false);
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, this.outputFileName);
		}
		System.out.println(this.SelectedTemplates.length
				+ " total proteins sequences sequenced");
		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			System.out.println("Details for sequence " + (i + 1) + "/"
					+ this.SelectedTemplates.length);
			String[] seqs = this.ReconstructPath(i);
			System.out.println("Final Sequence: " + seqs[0]);

			System.out.println("Support: " + seqs[1]);
			System.out.println("Overlap: " + seqs[2]);
			Writer.write("# ----FINAL SEQUENCE INFO----\n");
			Writer.write("# Final Seq: " + seqs[0] + "\n");
			Writer.write("# Sitewise Support: " + seqs[1] + "\n");
			Writer.write("# Sitewise Overlap: " + seqs[2] + "\n");
			String Temp = seqs[0].replace(" ... ", "");
			String[] AAs = AntibodyUtils.SplitStringtoAA(Temp);

			Writer.write("# FinalLen: " + AAs.length + "\n");
		}
		try {
			Writer.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(7, this.outputFileName);
		}
		// String XLSFileName =
		// Utils.GetFileNameNoExtension(this.OutputFileName) + ".csv";
		// this.WriteSpreadsheetFile(XLSFileName);
		// If Report Generating is turned on, then we also need to produce a
		// number of extra files!!
		if (this.reportFileDir != null) {
			this.writeCSPSReportFiles();
		}
		return 0;

	}

	/**
	 * Once we've parsed the config file and the command line make sure there
	 * isn't any funny business
	 */
	private void validateParams() {

		// Check 1: If we are doing Spectral Library search, then we can't be
		// merging spectra together into PRMs
		if (this.runLibrarySearch && this.unRavelPRMCluster) {
			ErrorThrower
					.ThrowError(4,
							"Cannot use library search and PRM clustering together with GenoMS!");
		}

	}

	/**
	 * Creates a Peptide object for each unique peptide sequence
	 */
	private void mergeSpectraFromSamePeptide() {

		boolean localDebug = true;
		this.dbPepObj = new ArrayList<Peptide>();
		HashSet<String> peps = new HashSet<String>();

		// Iterate over all PRM spectra
		for (int i = 0; i < this.allPRMSpectra.length; ++i) {

			// Check to see if it has a DB annotation
			String pep1 = this.allPRMSpectra[i].getModdedDBAnnotation();
			if (pep1 == null || pep1.length() == 0)
				continue;

			//pep1 = Utils.GetUnModded(pep1);
			if (peps.contains(pep1))
				continue;

			Peptide newPep = new Peptide(this.allPRMSpectra[i]);
			peps.add(pep1);

			// Find all other PRM Spectra that support this peptide
			for (int j = i + 1; j < this.allPRMSpectra.length; ++j) {
				String pep2 = this.allPRMSpectra[j].getModdedDBAnnotation();
				if (pep2 == null || pep2.length() == 0)
					continue;
				//pep2 = Utils.GetUnModded(pep2);
				if (pep1.equals(pep2))
					newPep.addNewSpectrum(this.allPRMSpectra[j]);
			}
			this.dbPepObj.add(newPep);
		}
		System.out.println("Created " + this.dbPepObj.size()
				+ " peptide objects");

	}

	private boolean runDBSearchTool(boolean runInspectFlag2) {
		return runMSGFDB(this.runDBSearchFlag);
	}

	/**
	 * This might be hacky, but it adds the PRMs for the identified peptide to
	 * the PRM spectrum
	 */
	private void AugmentPRMSpectraForInspectResults() {
		for (int i = 0; i < this.allPRMSpectra.length; ++i) {
			if (this.allPRMSpectra[i].getModdedDBAnnotation() == null)
				continue;

			PRMSpectrum CurrSpectrum = this.allPRMSpectra[i].GetPRMSpectrum();
			// InspectAnnotation Ann =
			// this.AllSpectra[i].GetInspectAnnotation();
			String Peptide = Utils.GetUnModded(Utils
					.ConvertMutations(this.allPRMSpectra[i].getModdedDBAnnotation()));

			int[] pepPRMs = AntibodyUtils.GetSeqPRMScaled(Peptide);
			// System.out.println("Augmenting: " +
			// this.AllSpectra[i].GetScanNumber() + " with peptide " +
			// Ann.Annotation + "(" + Peptide + ")");
			// System.out.println(" - pepPRMs: " + Utils.JoinIntArray(pepPRMs,
			// " "));
			// System.out.println(" - currPRMS: " +
			// Utils.JoinIntArray(CurrSpectrum.PRMs, " "));
			// CurrSpectrum.addPRMs(pepPRMs,AntibodyUtils.MIN_PEAK_GROUP_SCORE,
			// this.PRMTolerance);
			CurrSpectrum.setPRMs(pepPRMs, AntibodyUtils.MIN_PEAK_GROUP_SCORE,
					this.PRMTolerance);

			// CurrSpectrum.PrecursorCharge = Ann.Charge;
			CurrSpectrum.PrecursorMass = AntibodyUtils
					.GetSeqMassScaled(Peptide);
			// CurrSpectrum.PrecursorMass = (int)(Utils.GetMFromMZ(Ann.ParentMZ,
			// Ann.Charge)*AntibodyUtils.MASS_SCALE);
			// System.out.println(" - newPRMS: " +
			// Utils.JoinIntArray(CurrSpectrum.PRMs, " "));
			// Utils.WaitForEnter();
			this.allPRMSpectra[i].SetPRMSpectrum(CurrSpectrum);

		}

	}

	private boolean writeCSPSReportFiles() {
		boolean LocalDebug = true;

		System.out.println("\n[GenoMS PostProcessing]: Writing report files");
		// Make subdirectories called spectra, assembly, and homology
		Utils.MakeDir(this.reportFileDir + File.separator + "spectra");
		Utils.MakeDir(this.reportFileDir + File.separator + "assembly");
		Utils.MakeDir(this.reportFileDir + File.separator + "homology");
		if (LocalDebug)
			System.out.println("Created directories for output in "
					+ this.reportFileDir);

		int contigCount = this.countContigs();

		// Create file of star spectra
		// For now this is just the PRM spectra
		// (5)
		// String starSpectrumFile = this.reportFileDir + File.separator
		// + "spectra" + File.separator + "stars.pklbin";
		String starSpectrumFile = this.reportFileDir + File.separator
				+ "spectra" + File.separator + "stars.mgf";
		if (LocalDebug)
			System.out.println("Creating the star spectrum file '"
					+ starSpectrumFile + "'...(" + this.allPRMSpectra.length + ")");
		// this.writePRMSpectraAsPKLBIN(starSpectrumFile);
		this.writePRMSpectraAsMGF(starSpectrumFile);

		// Create a file of contig spectra
		// (6)
		// String contigSpectrumFile = this.reportFileDir + File.separator
		// + "assembly" + File.separator + "sps_seqs.pklbin";
		String contigSpectrumFile = this.reportFileDir + File.separator
				+ "assembly" + File.separator + "sps_seqs.mgf";
		if (LocalDebug)
			System.out.println("Creating the contig spectrum file '"
					+ contigSpectrumFile + "'...");
		// this.writeSeedSpectraAsPKLBIN(contigSpectrumFile);
		this.writeSeedSpectraAsMGF(contigSpectrumFile);

		// Create file of a bruijn graphs
		// (7)
		String ABFileName = this.reportFileDir + File.separator + "assembly"
				+ File.separator + "component_info.bin";
		if (LocalDebug)
			System.out.println("Creating the AB file '" + ABFileName + "'...");
		this.writeABFile(ABFileName);

		// Write reference protein to file
		// prot_id.fasta
		// (8)
		String proteinFile = this.reportFileDir + File.separator
				+ "protid.fasta";
		if (LocalDebug)
			System.out.println("Creating the protein file '" + proteinFile
					+ "'...");
		// correct contig sequences to be all letters, no numbers

		String[] finalSeqs = this.createProteinFromTemplates(proteinFile);

		// Create file with contig-protein mapping
		// homology/contigs_mp_all.bin
		/*
		 * Reuses the binary array format (*.bin), unsigned int (4-byte) values,
		 * where the i-th row indicates the match information for the i-th
		 * spectrum (3 columns): 1. Matched protein index (0-based) 2. Number of
		 * modifications in the best-scoring spectrum/protein match 3. Direct
		 * (0) or reversed (1) spectrum/protein match
		 */
		// (9)
		// String otherMapFileName = this.reportFileDir + File.separator
		// + "homology" + File.separator + "contigs_mp.bin";
		String mapFileName = this.reportFileDir + File.separator + "homology"
				+ File.separator + "contigs_mp.bin";
		if (LocalDebug)
			System.out.println("Creating the contig-protein map '"
					+ mapFileName + "'...");
		this.createContigProteinMap(mapFileName, contigCount);
		// Utils.copyFile(mapFileName, otherMapFileName);

		// Create file with contig-protein alignment
		// homology/contigs_midx_all.pklbin
		/*
		 * Reuses the spectrum file format (.pklbin): the i-th spectrum has the
		 * set of matched masses for the i-th spectrum against the matched
		 * protein the k-th mass/intensity pair (a,b) indicates the spectrum
		 * peak index a and the protein mass index b for the k-th matched
		 * spectrum peak (a protein of length N has N+1 mass indexes: 0 to N)
		 */
		// (10)
		// String otherAlnIndexFileName = this.reportFileDir + File.separator
		// + "homology" + File.separator + "contigs_midx.mgf";
		String alnIndexFileName = this.reportFileDir + File.separator
				+ "homology" + File.separator + "contigs_midx.mgf";
		// String alnMassFileName = this.reportFileDir + File.separator
		// + "homology" + File.separator + "contigs_matches.mgf";
		/*
		 * String otherAlnIndexFileName = this.reportFileDir + File.separator +
		 * "homology" + File.separator + "contigs_midx.pklbin"; String
		 * alnIndexFileName = this.reportFileDir + File.separator + "homology" +
		 * File.separator + "contigs_midx_all.pklbin"; String alnMassFileName =
		 * this.reportFileDir + File.separator + "homology" + File.separator +
		 * "contigs_matches_all.pklbin";
		 */

		if (LocalDebug)
			System.out.println("Creating the contig protein alignment file '"
					+ alnIndexFileName + "'...");
		// this.createContigProteinAlignment(alnIndexFileName, alnMassFileName,
		// finalSeqs);
		this.createContigProteinAlignmentMGF(alnIndexFileName, finalSeqs);
		// Utils.copyFile(alnIndexFileName, otherAlnIndexFileName);
		// Create a file containing all contigs that match a protein (same as
		// sps_seqs.pklbin)
		// (11)
		String contigFileName = this.reportFileDir + File.separator + "spectra"
				+ File.separator + "contigs.mgf";
		// String contigFileName = this.reportFileDir + File.separator +
		// "spectra"
		// + File.separator + "contigs.pklbin";
		if (LocalDebug)
			System.out.println("Creating the mapped contig file '"
					+ contigFileName + "'...");
		Utils.copyFile(contigSpectrumFile, contigFileName);

		// Create files for mapping indices of contigs that matched a protein.
		// Since contigs.pklbin is the same as sps_seqs.pklbin, this is trivial
		/*
		 * Two-dimensional arrays of values: Number of rows (unsigned int, 4
		 * bytes) Number of columns (unsigned int, 4 bytes) N_rows * N_cols
		 * values, one row at a time (ie., concatenated rows) Size of values
		 * written to / read from file depends on template type used for Load /
		 * Save methods
		 */
		// (12)
		String contigIndexMap = this.reportFileDir + File.separator + "spectra"
				+ File.separator + "contigs_indices.bin";
		if (LocalDebug)
			System.out.println("Creating the contig index map '"
					+ contigIndexMap + "'...");
		this.createContigMap(contigIndexMap, contigCount);

		// Create a file of contig names
		// (18)
		String contigNames = this.reportFileDir + File.separator + "homology"
				+ File.separator + "ref_sps_names.txt";
		if (LocalDebug)
			System.out.println("Creating the contig name file '" + contigNames
					+ "'...");
		this.createContigNameFile(contigNames);

		// Copy protien to fasta
		// (19)
		// String inputFastaName = this.ReportFileDir + File.separator +
		// "homology"
		String temp = this.reportFileDir + File.separator + "homology"
				+ File.separator + "homglue_ref_mp.bin";
		Utils.copyFile(mapFileName, temp);

		temp = this.reportFileDir + File.separator + "homology"
				+ File.separator + "homglue_ref_midx.mgf";
		// temp = this.reportFileDir + File.separator + "homology"
		// + File.separator + "homglue_ref_midx.pklbin";
		Utils.copyFile(alnIndexFileName, temp);

		temp = this.reportFileDir + File.separator + "homology"
				+ File.separator + "homglue_matches.mgf";
		// temp = this.reportFileDir + File.separator + "homology"
		// + File.separator + "homglue_matches.pklbin";
		// Utils.createEmptyFile(temp);
		// Utils.copyFile(alnMassFileName, temp);
		Utils.copyFile(contigFileName, temp);

		temp = this.reportFileDir + File.separator + "homology"
				+ File.separator + "homglue_matches_mp.bin";
		// Utils.createEmptyFile(temp);
		// Utils.copyFile(mapFileName, temp);
		this.createCSPSContigProteinMatch(temp, contigCount);

		temp = this.reportFileDir + File.separator + "homology"
				+ File.separator + "homglue_matches_midx.mgf";
		// temp = this.reportFileDir + File.separator + "homology"
		// + File.separator + "homglue_matches_midx.pklbin";
		// Utils.createEmptyFile(temp);
		// this.createCSPSContigMap(temp);
		Utils.copyFile(alnIndexFileName, temp);

		return true;
	}

	/*
	 * private boolean createCSPSContigMap(String contigSpectrumFile) { return
	 * Template.WriteAllSeedSpectraToBlankPKLBIN(contigSpectrumFile,
	 * this.SelectedTemplates);
	 * 
	 * }
	 */

	/**
	 * We don't want our GenoMS guys to be CSPS contigs, so we match them to a
	 * -1 protein
	 * 
	 * @param mapFileName
	 * @param contigCount
	 * @return
	 */
	private boolean createCSPSContigProteinMatch(String mapFileName,
			int contigCount) {
		DataOutputStream f = null;

		try {
			f = new DataOutputStream(new FileOutputStream(mapFileName));
			f.writeInt(Integer.reverseBytes(contigCount));
			f.writeInt(Integer.reverseBytes(3));

		} catch (IOException E) {
			ErrorThrower.ThrowError(5, mapFileName);
		}
		int proteinIndex = 0;
		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			if (this.SelectedTemplates[i] == null)
				continue;
			Template currTemplate = this.SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					try {
						f.writeInt(Integer.reverseBytes(proteinIndex)); // Protein
																		// is
																		// the
																		// first
																		// one
																		// in
																		// the
																		// file
						f.writeInt(Integer.reverseBytes(0)); // Modifications,
																// no
																// information
																// here
						f.writeInt(Integer.reverseBytes(0)); // All 0 since all
																// contigs are
																// direct
																// alignments,
																// not reversed
					} catch (IOException E) {
						ErrorThrower.ThrowError(7, mapFileName);
					}

					currSeed = currSeed.NextSeed;
				}
				currTemplate = currTemplate.NextInList;
			}
			proteinIndex += 1;
		}

		try {
			f.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(8, mapFileName);
		}

		return true;

	}

	/*
	 * private void correctContigSeqs() { for(int i = 0; i <
	 * this.SelectedTemplates.length; ++i) { if(this.SelectedTemplates[i] ==
	 * null) continue; Template currTemplate = this.SelectedTemplates[i];
	 * while(currTemplate != null) { Seed currSeed = currTemplate.Anchors;
	 * while(currSeed != null) {
	 * 
	 * String[] seq = AntibodyUtils.SplitStringtoAA(currSeed.SeedSequence);
	 * for(int j = 0; j < seq.length; ++j) { if(seq[j].charAt(0) == '[') //this
	 * is a number not a char!! { int mass = (int)
	 * (Double.parseDouble(seq[j].substring
	 * (1,seq[j].length()-1))*AntibodyUtils.MASS_SCALE); String s =
	 * AntibodyUtils.FindClosestDoubleAAScaledSeqs(mass, this.PRMTolerance);
	 * seq[j] = s; } } //System.out.println("Old Seed:" +
	 * currSeed.SeedSequence); currSeed.SeedSequence =
	 * Utils.JoinStringArray(seq, ""); //System.out.println("New Seed:" +
	 * currSeed.SeedSequence); //Utils.WaitForEnter(); currSeed =
	 * currSeed.NextSeed;
	 * 
	 * } currTemplate = currTemplate.NextInList; } }
	 * 
	 * }
	 */
	private int countContigs() {
		int contigCount = 0;

		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			if (this.SelectedTemplates[i] == null)
				continue;
			Template currTemplate = this.SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					currSeed = currSeed.NextSeed;
					contigCount += 1;
				}
				currTemplate = currTemplate.NextInList;
			}
		}
		return contigCount;
	}

	/**
	 * //Create files for mapping indices of contigs that matched a protein.
	 * //Since contigs.pklbin is the same as sps_seqs.pklbin, this is trivial /*
	 * Two-dimensional arrays of values: Number of rows (unsigned int, 4 bytes)
	 * Number of columns (unsigned int, 4 bytes) N_rows * N_cols values, one row
	 * at a time (ie., concatenated rows) Size of values written to / read from
	 * file depends on template type used for Load / Save methods
	 * 
	 * @param contigIndexMap
	 */

	private boolean createContigMap(String contigIndexMap, int contigCount) {
		DataOutputStream f = null;

		try {
			f = new DataOutputStream(new FileOutputStream(contigIndexMap));
			f.writeInt(Integer.reverseBytes(contigCount));
			f.writeInt(Integer.reverseBytes(1));
			for (int i = 0; i < contigCount; ++i) {
				f.writeInt(Integer.reverseBytes(i));
			}
			f.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, contigIndexMap);
		}
		return true;

	}

	/**
	 * Writes a file of contig names Row i is the name for contig i
	 * 
	 * @param contigNames
	 */
	private boolean createContigNameFile(String contigNames) {
		FileWriter f = null;

		try {
			f = new FileWriter(contigNames);
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, contigNames);
		}
		int contigCount = 1;
		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			if (this.SelectedTemplates[i] == null)
				continue;
			Template currTemplate = this.SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					try {
						f.write("GenoMS:" + contigCount + "\n");
					} catch (IOException E) {
						ErrorThrower.ThrowError(7, contigNames);
					}
					currSeed = currSeed.NextSeed;
					contigCount += 1;
				}
				currTemplate = currTemplate.NextInList;
			}
		}
		try {
			f.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(8, contigNames);
		}
		return true;

	}

	/**
	 * Creates a protein from the concatenated templates used by GenoMS. May
	 * include multiple sequences if running heavy and light chain together
	 * 
	 * @param proteinFile
	 *            The file to write to (FASTA format)
	 * @return Returns true if file is successfully created, false otherwise
	 */
	private String[] createProteinFromTemplates(String proteinFile) {

		boolean LocalDebug = false;

		String[] proteinSeqs = new String[this.SelectedTemplates.length];
		String[] proteinNames = new String[this.SelectedTemplates.length];

		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			if (this.SelectedTemplates[i] == null)
				continue;
			// We also construct the predicted protein sequences here.
			proteinSeqs[i] = "";
			proteinNames[i] = ">";

			// Iterate through each template
			Template currTemplate = this.SelectedTemplates[i];
			String updatedTemplateSeq = null;
			// What we should do here is check that the current template does
			// not align with the prevous sequence, but if it does, then we can
			// truncate or whatever.
			while (currTemplate != null) {
				if (currTemplate.Sequence == null) {
					if (LocalDebug) {
						System.out.println("Next template ("
								+ currTemplate.TemplateName
								+ " has no sequence (and no anchors)");
					}
					// proteinSeqs[i] += "X";
					currTemplate = currTemplate.NextInList;
					continue;
				}
				Seed prevSeed = null;
				Seed currSeed = currTemplate.anchors;

				if (LocalDebug)
					System.out.println("Considering template: "
							+ currTemplate.TemplateName);

				if (currSeed == null) {

					proteinSeqs[i] += Utils
							.CleanUpPeptideString(currTemplate.Sequence);
					proteinNames[i] += currTemplate.TemplateName.split("[|]")[0]
							+ "+";
					currTemplate = currTemplate.NextInList;
					if (LocalDebug) {
						System.out
								.println("This template has no anchors, so we just add the sequence as is");
						Utils.WaitForEnter();
					}
					continue;
				}

				// The updatedTemplateSeq is null if this is the first template,
				// or if the previous seed did not overlap with this template
				if (updatedTemplateSeq == null)
					updatedTemplateSeq = Utils
							.CleanUpPeptideString(currTemplate.Sequence);

				while (currSeed != null) {
					if (LocalDebug)
						System.out.println("Aligning " + currSeed.SeedSequence);
					String[] seq = AntibodyUtils
							.SplitStringtoAA(currSeed.SeedSequence);
					for (int j = 0; j < seq.length; ++j) {
						if (seq[j].charAt(0) == '[') // this is a number not a
														// char!!
						{
							int mass = (int) (Double.parseDouble(seq[j]
									.substring(1, seq[j].length() - 1)) * AntibodyUtils.MASS_SCALE);
							String s = AntibodyUtils
									.FindClosestDoubleAAScaledSeqs(mass,
											this.PRMTolerance);
							seq[j] = s;
						}
					}
					// System.out.println("Old Seed:" + currSeed.SeedSequence);
					String fixedSequence = Utils.JoinStringArray(seq, "");

					// FIrst check, if it is unchanged from template, then just
					// continue
					if (updatedTemplateSeq.indexOf(fixedSequence) >= 0) {
						if (LocalDebug)
							System.out
									.println("Seed is unchanged from template, so continue");
						currSeed = currSeed.NextSeed;
						continue;
					}
					Object[] seedTemplateAlignment = AntibodyUtils
							.GetLongestAlignmentDetailedSeqs(fixedSequence,
									updatedTemplateSeq, this.PRMTolerance, 2, 2);
					// Object[] test
					// AntibodyUtils.GetLongestAlignmentStrictPRMs(PRMsA, PRMsB,
					// Tolerance)
					int[][] alignmentScores = (int[][]) (seedTemplateAlignment[1]);
					int[][] prevCell = (int[][]) (seedTemplateAlignment[0]);

					if (LocalDebug)
						System.out.println("Finished alignment");
					int bestScore = 0;
					int[] endIndices = new int[2];
					int[] startIndices = new int[2];

					// Find the best alignment, and it's start and end
					for (int j = 0; j < alignmentScores.length; ++j) {
						for (int k = 0; k < alignmentScores[0].length; ++k) {
							if (alignmentScores[j][k] > bestScore) {
								bestScore = alignmentScores[j][k];
								endIndices[0] = j; // PRM Index of last seed prm
													// in alignment
								endIndices[1] = k; // PRM Index of last template
													// prm in alignment
							}
						}
					}

					if (LocalDebug) {
						System.out.println("Alignment Score: " + bestScore);
						System.out.println("prevCell array: " + prevCell.length
								+ "x" + prevCell[0].length);
					}
					int currRow = endIndices[0];
					int currCol = endIndices[1];

					while (prevCell[currRow][currCol] >= 0) {
						// if(LocalDebug)
						// {
						// System.out.println("[" + currRow + "," + currCol +
						// "]");
						// Utils.WaitForEnter();
						// }
						int currCell = prevCell[currRow][currCol];
						currRow = currCell / prevCell[0].length;
						currCol = currCell % prevCell[0].length;
						// if(LocalDebug)
						// {
						// // System.out.println("new: [" + currRow + "," +
						// currCol + "]");
						// Utils.WaitForEnter();
						// }

					}

					startIndices[0] = currRow; // PRM Index of first seed prm in
												// alignment
					startIndices[1] = currCol; // PRM Index of first template
												// prm in alignment

					// Sub out the sequence of the template for the seed in the
					// aligned regions
					int startOfExcision = Math.max(0, startIndices[1]
							- startIndices[0]);
					int endOfExcision = Math
							.min(updatedTemplateSeq.length(),
									endIndices[1]
											+ (alignmentScores.length - 1 - endIndices[0]));

					String newSeq = updatedTemplateSeq.substring(0,
							startOfExcision)
							+ fixedSequence
							+ updatedTemplateSeq.substring(endOfExcision);

					if (LocalDebug) {
						System.out.println("oldSeq: " + updatedTemplateSeq);
						System.out.println("newSeq: " + newSeq);
						Utils.WaitForEnter();
					}
					updatedTemplateSeq = newSeq;

					prevSeed = currSeed;
					currSeed = currSeed.NextSeed;
				}

				proteinSeqs[i] += updatedTemplateSeq;
				proteinNames[i] += currTemplate.TemplateName.split("[|]")[0]
						+ "+";
				currTemplate = currTemplate.NextInList;

				updatedTemplateSeq = null;

				/**
				 * CHECK TO SEE IF WE OVERLAP WITH NEXT TEMPLATE< WOOOOOO!!!
				 */
				if (prevSeed != null && prevSeed.hasMergedAcrossTemplates) {
					if (LocalDebug)
						System.out
								.println("This seed has overlap with next template!!");

					// Find the next template with sequence!
					while (currTemplate != null
							&& currTemplate.Sequence == null)
						currTemplate = currTemplate.NextInList;

					// If we found the next one
					if (currTemplate != null) {
						updatedTemplateSeq = Utils
								.CleanUpPeptideString(currTemplate.Sequence);
						if (LocalDebug) {
							System.out.println(prevSeed.SeedSequence);
							System.out.println(updatedTemplateSeq);
						}
						String[] seq = AntibodyUtils
								.SplitStringtoAA(prevSeed.SeedSequence);
						for (int j = 0; j < seq.length; ++j) {
							if (seq[j].charAt(0) == '[') // this is a number not
															// a char!!
							{
								int mass = (int) (Double.parseDouble(seq[j]
										.substring(1, seq[j].length() - 1)) * AntibodyUtils.MASS_SCALE);
								String s = AntibodyUtils
										.FindClosestDoubleAAScaledSeqs(mass,
												this.PRMTolerance);
								seq[j] = s;
							}
						}
						// System.out.println("Old Seed:" +
						// currSeed.SeedSequence);
						String fixedSequence = Utils.JoinStringArray(seq, "");

						Object[] seedTemplateAlignment = AntibodyUtils
								.GetLongestAlignmentDetailedSeqs(fixedSequence,
										updatedTemplateSeq, this.PRMTolerance,
										2, 2);
						int[][] alignmentScores = (int[][]) (seedTemplateAlignment[1]);
						int[][] prevCell = (int[][]) (seedTemplateAlignment[0]);

						if (LocalDebug)
							System.out.println("Finished alignment");
						int bestScore = 0;
						int[] endIndices = new int[2];
						int[] startIndices = new int[2];

						// Find the best alignment, and it's start and end
						for (int j = 0; j < alignmentScores.length; ++j) {
							for (int k = 0; k < alignmentScores[0].length; ++k) {
								if (alignmentScores[j][k] > bestScore) {
									bestScore = alignmentScores[j][k];
									endIndices[0] = j; // PRM Index of last seed
														// prm in alignment
									endIndices[1] = k; // PRM Index of last
														// template prm in
														// alignment
								}
							}
						}

						if (LocalDebug) {
							System.out.println("Alignment Score: " + bestScore);
							System.out.println("prevCell array: "
									+ prevCell.length + "x"
									+ prevCell[0].length);
						}

						// If the alignment is good, wOOOWOOOW
						if (bestScore >= 5) {

							int currRow = endIndices[0];
							int currCol = endIndices[1];

							while (prevCell[currRow][currCol] >= 0) {
								// if(LocalDebug)
								// {
								// System.out.println("[" + currRow + "," +
								// currCol + "]");
								// Utils.WaitForEnter();
								// }
								int currCell = prevCell[currRow][currCol];
								currRow = currCell / prevCell[0].length;
								currCol = currCell % prevCell[0].length;
								// if(LocalDebug)
								// {
								// System.out.println("new: [" + currRow + "," +
								// currCol + "]");
								// Utils.WaitForEnter();
								// }

							}

							startIndices[0] = currRow; // PRM Index of first
														// seed prm in alignment
							startIndices[1] = currCol; // PRM Index of first
														// template prm in
														// alignment

							// Sub out the sequence of the template for the seed
							// in the aligned regions
							int startOfExcision = Math.max(0, startIndices[1]
									- startIndices[0]);
							int endOfExcision = Math
									.min(updatedTemplateSeq.length(),
											endIndices[1]
													+ (alignmentScores.length - 1 - endIndices[0]));

							// We want it to start close to the beginning of the
							// template
							if (startOfExcision <= 3) {
								String newSeq = updatedTemplateSeq
										.substring(endOfExcision);

								if (LocalDebug) {
									System.out.println("oldSeq: "
											+ updatedTemplateSeq);
									System.out.println("newSeq: " + newSeq);
									Utils.WaitForEnter();
								}

								updatedTemplateSeq = newSeq;
							}
						}
					}
				}

			}
		}

		FileWriter f = null;

		try {
			f = new FileWriter(proteinFile);
			for (int i = 0; i < proteinSeqs.length; ++i) {
				if (proteinSeqs[i].length() == 0)
					continue;
				f.write(proteinNames[i] + " GenoMS Predicted Sequence\n");
				f.write(proteinSeqs[i] + "\n");

				if (LocalDebug) {
					System.out.println("[" + i + "]:" + proteinNames[i]);
					System.out.println(proteinSeqs[i]);
				}
			}
			f.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, proteinFile);
		}

		return proteinSeqs;

	}

	/**
	 * Reuses the spectrum file format (.pklbin): the i-th spectrum has the set
	 * of matched masses for the i-th spectrum against the matched protein the
	 * k-th mass/intensity pair (a,b) indicates the spectrum peak index a and
	 * the protein mass index b for the k-th matched spectrum peak (a protein of
	 * length N has N+1 mass indexes: 0 to N)
	 * 
	 * @param alnFileName
	 * @param alnMassFileName
	 * @param finalSeqs
	 * @return Returns true if file is created correctly, false otherwise
	 */
	private boolean createContigProteinAlignment(String alnFileName,
			String alnMassFileName, String[] finalSeqs) {

		boolean LocalDebug = false;
		DataOutputStream f = null;
		DataOutputStream f2 = null;

		// FileWriter f = Utils.openFileWriter(alnFileName);
		// FileWriter f2 = Utils.openFileWriter(alnMassFileName);

		ArrayList<ArrayList<int[]>> peakAlignments = new ArrayList<ArrayList<int[]>>();
		ArrayList<ArrayList<float[]>> massAlignments = new ArrayList<ArrayList<float[]>>();
		ArrayList<Double> parentMasses = new ArrayList<Double>();

		int finalSeqIndex = 0;
		for (int i = 0; i < SelectedTemplates.length; ++i) {
			if (SelectedTemplates[i] == null)
				continue;
			Template currTemplate = SelectedTemplates[i];
			int[] seqPeaks = AntibodyUtils
					.GetSeqPRMScaled(finalSeqs[finalSeqIndex]);
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {

					int[] contigPeaks = AntibodyUtils
							.GetSeqPRMScaled(currSeed.SeedSequence);

					Object[] seedTemplateAlignment = AntibodyUtils
							.GetLongestAlignmentDetailedPRMs(contigPeaks,
									seqPeaks, this.PRMTolerance, 2, 2);
					ArrayList<int[]> alnPeaks = new ArrayList<int[]>();
					ArrayList<float[]> alnMasses = new ArrayList<float[]>();
					int[][] alignmentScores = (int[][]) (seedTemplateAlignment[1]);
					int[][] prevCell = (int[][]) (seedTemplateAlignment[0]);

					int bestScore = 0;
					int[] endIndices = new int[2];
					int[] startIndices = new int[2];
					float[] endMasses = new float[2];
					// Find the best alignment, and it's start and end
					for (int j = 0; j < alignmentScores.length; ++j) {
						for (int k = 0; k < alignmentScores[0].length; ++k) {
							if (alignmentScores[j][k] > bestScore) {
								bestScore = alignmentScores[j][k];
								endIndices[0] = j; // PRM Index of last seed prm
													// in alignment
								endIndices[1] = k; // PRM Index of last template
													// prm in alignment
								endMasses[0] = (float) (((double) contigPeaks[j]) / (double) AntibodyUtils.MASS_SCALE);
								endMasses[1] = (float) (((double) seqPeaks[k]) / (double) AntibodyUtils.MASS_SCALE);
							}
						}
					}
					alnPeaks.add(endIndices);
					alnMasses.add(endMasses);

					int currRow = endIndices[0];
					int currCol = endIndices[1];

					while (prevCell[currRow][currCol] >= 0) {
						int currCell = prevCell[currRow][currCol];
						currRow = currCell / prevCell[0].length;
						currCol = currCell % prevCell[0].length;
						int[] newAln = new int[2];
						float[] newMasses = new float[2];
						newAln[0] = currRow;
						newAln[1] = currCol;
						newMasses[0] = (float) (((double) contigPeaks[currRow]) / (double) AntibodyUtils.MASS_SCALE);
						newMasses[1] = (float) (((double) seqPeaks[currCol]) / (double) AntibodyUtils.MASS_SCALE);

						alnPeaks.add(0, newAln);
						alnMasses.add(0, newMasses);
					}

					startIndices[0] = currRow; // PRM Index of first seed prm in
												// alignment
					startIndices[1] = currCol; // PRM Index of first template
												// prm in alignment
					// alnPeaks.add(0,startIndices);

					peakAlignments.add(alnPeaks);
					massAlignments.add(alnMasses);
					parentMasses.add(new Double(Utils
							.GetSeqMass(currSeed.SeedSequence)));
					currSeed = currSeed.NextSeed;
				}
				currTemplate = currTemplate.NextInList;

			}
			finalSeqIndex += 1;
		}

		/**
		 * Write that shit
		 */
		try {
			f = new DataOutputStream(new FileOutputStream(alnFileName));
			f.writeInt(Integer.reverseBytes(peakAlignments.size()));

			f2 = new DataOutputStream(new FileOutputStream(alnMassFileName));
			f2.writeInt(Integer.reverseBytes(peakAlignments.size()));

			if (LocalDebug)
				System.out.println("Total spectra: " + peakAlignments.size());

			// Next is an array of spectrum sizes (shorts)
			for (int i = 0; i < peakAlignments.size(); ++i) {
				ArrayList currSpec = (ArrayList) (peakAlignments.get(i));
				f.writeShort(Short.reverseBytes((short) (currSpec.size())));
				f2.writeShort(Short.reverseBytes((short) (currSpec.size())));
				if (LocalDebug)
					System.out.println("peakList[" + i + "]: "
							+ (short) (currSpec.size()));

			}
			// Finally, each spectrum has a float value for the precursor mass,
			// precursor charge, and then the peak list (mass, log Odds).
			for (int i = 0; i < peakAlignments.size(); ++i) {
				ArrayList currSpec = (ArrayList) (peakAlignments.get(i));
				ArrayList currMasses = (ArrayList) (massAlignments.get(i));

				double currMass = ((Double) (parentMasses.get(i)))
						.doubleValue();
				f.writeInt(Integer.reverseBytes(Float
						.floatToIntBits((float) currMass)));
				f.writeInt(Integer.reverseBytes(Float
						.floatToIntBits((float) 1.0)));

				f2.writeInt(Integer.reverseBytes(Float
						.floatToIntBits((float) currMass)));
				f2.writeInt(Integer.reverseBytes(Float
						.floatToIntBits((float) 1.0)));
				if (LocalDebug)
					System.out.println("S" + i + ":" + currMass);
				for (int j = 0; j < currSpec.size(); ++j) {
					int[] currPeaks = (int[]) (currSpec.get(j));
					float[] currMassAln = (float[]) (currMasses.get(j));
					if (LocalDebug)
						System.out.println((float) (currPeaks[0]) + " "
								+ currPeaks[1]);
					f.writeInt(Integer.reverseBytes(Float
							.floatToIntBits((float) (currPeaks[0]))));
					f.writeInt(Integer.reverseBytes(Float
							.floatToIntBits((float) (currPeaks[1]))));

					f2.writeInt(Integer.reverseBytes(Float
							.floatToIntBits(currMassAln[0])));
					f2.writeInt(Integer.reverseBytes(Float
							.floatToIntBits(currMassAln[1])));
					// f.writeFloat((float)(currPeaks[0]));
					// f.writeFloat((float)(currPeaks[1]));
				}
			}
			f.close();
			f2.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, alnFileName + " and " + alnMassFileName);
		}

		return true;

	}

	/**
	 * Reuses the spectrum file format (.mgf): the i-th spectrum has the set of
	 * matched masses for the i-th spectrum against the matched protein the k-th
	 * mass/intensity pair (a,b) indicates the spectrum peak index a and the
	 * protein mass index b for the k-th matched spectrum peak (a protein of
	 * length N has N+1 mass indexes: 0 to N)
	 * 
	 * @param alnFileName
	 * @param alnMassFileName
	 * @param finalSeqs
	 * @return Returns true if file is created correctly, false otherwise
	 */
	private boolean createContigProteinAlignmentMGF(String alnFileName,
			String[] finalSeqs) {

		FileWriter f = Utils.openFileWriter(alnFileName);
		// FileWriter f2 = Utils.openFileWriter(alnMassFileName);

		// ArrayList peakAlignments = new ArrayList();
		// ArrayList massAlignments = new ArrayList();
		// ArrayList parentMasses = new ArrayList();

		int finalSeqIndex = 0;
		for (int i = 0; i < SelectedTemplates.length; ++i) {
			if (SelectedTemplates[i] == null)
				continue;
			Template currTemplate = SelectedTemplates[i];
			int[] seqPeaks = AntibodyUtils
					.GetSeqPRMScaled(finalSeqs[finalSeqIndex]);
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {

					int[] contigPeaks = AntibodyUtils
							.GetSeqPRMScaled(currSeed.SeedSequence);

					Object[] seedTemplateAlignment = AntibodyUtils
							.GetLongestAlignmentDetailedPRMs(contigPeaks,
									seqPeaks, this.PRMTolerance, 2, 2);
					ArrayList<int[]> alnPeaks = new ArrayList<int[]>();
					// ArrayList<float[]> alnMasses = new ArrayList<float[]>();
					int[][] alignmentScores = (int[][]) (seedTemplateAlignment[1]);
					int[][] prevCell = (int[][]) (seedTemplateAlignment[0]);

					int bestScore = 0;
					int[] endIndices = new int[2];
					int[] startIndices = new int[2];
					// float[] endMasses = new float[2];
					// Find the best alignment, and it's start and end
					for (int j = 0; j < alignmentScores.length; ++j) {
						for (int k = 0; k < alignmentScores[0].length; ++k) {
							if (alignmentScores[j][k] > bestScore) {
								bestScore = alignmentScores[j][k];
								endIndices[0] = j; // PRM Index of last seed prm
													// in alignment
								endIndices[1] = k; // PRM Index of last template
													// prm in alignment
								// endMasses[0] = (float) (((double)
								// contigPeaks[j]) / (double)
								// AntibodyUtils.MASS_SCALE);
								// endMasses[1] = (float) (((double)
								// seqPeaks[k]) / (double)
								// AntibodyUtils.MASS_SCALE);
							}
						}
					}
					alnPeaks.add(endIndices);
					// alnMasses.add(endMasses);

					int currRow = endIndices[0];
					int currCol = endIndices[1];

					while (prevCell[currRow][currCol] >= 0) {
						int currCell = prevCell[currRow][currCol];
						currRow = currCell / prevCell[0].length;
						currCol = currCell % prevCell[0].length;
						int[] newAln = new int[2];
						// float[] newMasses = new float[2];
						newAln[0] = currRow;
						newAln[1] = currCol;
						// newMasses[0] = (float) (((double)
						// contigPeaks[currRow]) / (double)
						// AntibodyUtils.MASS_SCALE);
						// newMasses[1] = (float) (((double) seqPeaks[currCol])
						// / (double) AntibodyUtils.MASS_SCALE);

						alnPeaks.add(0, newAln);
						// alnMasses.add(0, newMasses);
					}

					startIndices[0] = currRow; // PRM Index of first seed prm in
												// alignment
					startIndices[1] = currCol; // PRM Index of first template
												// prm in alignment
					// alnPeaks.add(0,startIndices);

					// Write to MGF
					Utils.writeLine(f, alnFileName, "BEGIN IONS\n");
					Utils.writeLine(f, alnFileName, "PEPMASS=0.00000\n");
					Utils.writeLine(f, alnFileName, "CHARGE=0+\n");
					Utils.writeLine(f, alnFileName, "ACTIVATION=CID\n");
					Utils.writeLine(f, alnFileName, "TITLE=Scan Number: 0\n");
					Utils.writeLine(f, alnFileName, "SCANS=0\n");

					/*
					 * Utils.writeLine(f2, alnMassFileName, "BEGIN IONS");
					 * Utils.writeLine(f2, alnMassFileName,
					 * "PEPMASS=0.00000\n"); Utils.writeLine(f2,
					 * alnMassFileName, "CHARGE=0+\n"); Utils.writeLine(f2,
					 * alnMassFileName, "ACTIVATION=CID\n"); Utils.writeLine(f2,
					 * alnMassFileName, "TITLE=Scan Number: 0\n");
					 * Utils.writeLine(f2, alnMassFileName, "SCANS=0\n");
					 */
					for (int j = 0; j < alnPeaks.size(); ++j) {
						int[] currPeaks = alnPeaks.get(j);
						// float[] currMassAln = alnMasses.get(j);

						Utils.writeLine(f, alnFileName, currPeaks[0] + " "
								+ currPeaks[1] + "\n");
						// Utils.writeLine(f2, alnMassFileName, currMassAln[0] +
						// " " + currMassAln[1] + "\n");

					}

					Utils.writeLine(f, alnFileName, "END IONS\n");
					// Utils.writeLine(f2, alnMassFileName, "END IONS");
					// peakAlignments.add(alnPeaks);
					// massAlignments.add(alnMasses);
					// parentMasses.add(new Double(Utils
					// .GetSeqMass(currSeed.SeedSequence)));
					currSeed = currSeed.NextSeed;
				}
				currTemplate = currTemplate.NextInList;

			}
			finalSeqIndex += 1;
		}

		Utils.closeFileWriter(f, alnFileName);
		// Utils.closeFileWriter(f2, alnMassFileName);
		return true;

	}

	/**
	 * Reuses the binary array format (*.bin), unsigned int (4-byte) values,
	 * where the i-th row indicates the match information for the i-th spectrum
	 * (3 columns): 1. Matched protein index (0-based) 2. Number of
	 * modifications in the best-scoring spectrum/protein match (not filled in,
	 * but left blank) 3. Direct (0) or reversed (1) spectrum/protein match
	 * (always 0)
	 * 
	 * @param mapFileName
	 */
	private boolean createContigProteinMap(String mapFileName, int contigCount) {

		DataOutputStream f = null;

		try {
			f = new DataOutputStream(new FileOutputStream(mapFileName));
			f.writeInt(Integer.reverseBytes(contigCount));
			f.writeInt(Integer.reverseBytes(3));

		} catch (IOException E) {
			ErrorThrower.ThrowError(5, mapFileName);
		}
		int proteinIndex = 0;
		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			if (this.SelectedTemplates[i] == null)
				continue;
			Template currTemplate = this.SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					try {
						f.writeInt(Integer.reverseBytes(proteinIndex)); // Protein
																		// is
																		// the
																		// first
																		// one
																		// in
																		// the
																		// file
						f.writeInt(Integer.reverseBytes(0)); // Modifications,
																// no
																// information
																// here
						f.writeInt(Integer.reverseBytes(0)); // All 0 since all
																// contigs are
																// direct
																// alignments,
																// not reversed
					} catch (IOException E) {
						ErrorThrower.ThrowError(7, mapFileName);
					}

					currSeed = currSeed.NextSeed;
				}
				currTemplate = currTemplate.NextInList;
			}
			proteinIndex += 1;
		}

		try {
			f.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(8, mapFileName);
		}

		return true;

	}

	/**
	 * Create the A-Bruijn graph for each component/extended seed File Format
	 * nNumComponents [PROJ:int] per component: num_specs [PROJ:short] per spec:
	 * specIndex [PROJ:int, 0-based], specFlipped [PROJ:short 0/1]
	 * num_ABruijn_vertices [PROJ:short] per vertex: num_peaks_in_vertex
	 * [PROJ:short] per peak: specIndex [PROJ:int, 0-based], peakMass
	 * [PROJ:float]
	 * 
	 * THIS IS WHAT GETS WRITTEN (NOT THE ABOVE STUFF) numUsedSpectra [PROJ:int]
	 * (a) numUsedSpectra-by-1 unsigned int array of all used spectrum indices
	 * (b) numUsedSpectra-by-2 unsigned short array of [PROJ:component
	 * index,specFlipped] (c) numComponents [PROJ:int] (d) numComponents-by-1
	 * unsigned short array of number of ABruijn vertices per component (e)
	 * numABVertices-by-1 unsigned short array of number of spectrum peaks per
	 * ABruijn vertex (f) totNumSpecPeaks-by-1 unsigned int array of spectrum
	 * index per peak (g) totNumSpecPeaks-by-1 float array of peak masses per
	 * spectrum peak (h)
	 * 
	 * 
	 * This program makes several HUGE assumptions 1. The input file is MGF, so
	 * the spectrum scan number is the same as the order of the spectrum in the
	 * input file 2. All input spectra are in one file, so the index is
	 * sufficient information to distinguish them.
	 * 
	 * @param aBFileName
	 *            The file name to write the AB info
	 */
	private boolean writeABFile(String aBFileName) {

		boolean LocalDebug = false;
		if (LocalDebug)
			System.out.println("Creating AB -graph in " + aBFileName);

		DataOutputStream f = null;

		// The number of spectra used across all contigs. (a)
		int spectraCount = 0;

		// List of allSpectra used in the contigs (b)
		ArrayList<Integer> allSpectra = new ArrayList<Integer>();
		// List of spectra info, same indices as for allSpectra. spectraInfo[i]
		// = [componentIndex,isFlipped] (c)
		ArrayList<int[]> spectraInfo = new ArrayList<int[]>();

		// The number of contigs (d)
		int contigCount = 0;

		// peakCounts[i] = The number of vertices in contig i (e)
		ArrayList<Short> peakCounts = new ArrayList<Short>();

		// ABVertexCounts[i] = the number of spectra supporting vertex i (f)
		ArrayList<Short> ABVertexCounts = new ArrayList<Short>();

		// In order of contig/vertex, the info for each spectrum
		ArrayList<Integer> specIndex = new ArrayList<Integer>();
		ArrayList<Float> peakMasses = new ArrayList<Float>();

		for (int i = 0; i < SelectedTemplates.length; ++i) {
			if (SelectedTemplates[i] == null)
				continue;
			Template currTemplate = SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {

					// Add the number of vertices for this contig
					peakCounts
							.add(new Short((short) (currSeed.SeedPRMs.length)));

					if (LocalDebug) {
						System.out
								.println("New Seed: " + currSeed.SeedSequence);
						currSeed.DebugPrint();
					}
					// collect all of the supporting spectra
					for (int j = 0; j < currSeed.PeakSupporting.length; ++j) {
						if (LocalDebug)
							System.out.println("-Investigating vertex " + j);

						// The list of spectra supporting this vertex
						// (FileName_ScanNum)
						String[] supportingPeaks = (String[]) (currSeed.PeakSupporting[j]);
						ABVertexCounts.add(new Short(
								(short) (supportingPeaks.length)));

						if (LocalDebug)
							System.out.println(" - There are "
									+ supportingPeaks.length
									+ " spectra supporting");

						// Add the spectrum index to set of all Spectra, and its
						// info
						for (int k = 0; k < supportingPeaks.length; ++k) {

							if (supportingPeaks[k] == null) {
								System.out
										.println("WTF!!, supporting peaks[k] is null!");
								Utils.WaitForEnter();
								continue;
							}
							String[] keyInfo = supportingPeaks[k].split("_");

							Integer ScanNumber = (new Integer(
									keyInfo[keyInfo.length - 1]));

							if (LocalDebug) {
								System.out.println("  - Spectrum: "
										+ supportingPeaks[k]);

							}
							if (!allSpectra.contains(ScanNumber)) {
								allSpectra.add(ScanNumber);
								int[] specInfo = new int[2];
								specInfo[0] = contigCount;
								specInfo[1] = 0;
								spectraInfo.add(specInfo);
							}
							// Add this spectrum to the list of details for the
							// vertex
							specIndex.add(ScanNumber);

							String key = ScanNumber.intValue() + "";
							// Determine what the mass is for this spectrum and
							// whats up
							SpectrumNode specNode = (SpectrumNode) (this.spectrumNodeIndex
									.get(key));

							if (specNode == null) {
								System.out.println("WTF!! specNode for scan "
										+ ScanNumber + " is null!!");
							}

							if (LocalDebug) {
								if (specNode.getModdedDBAnnotation() != null)
									System.out.println("  -  Peptide: "
											+ specNode.getModdedDBAnnotation());
								if (specNode.HMMAnnotation != null) {
									System.out.println("ANCHOR SEQ: "
											+ specNode.GetMSAnnotation()
													.GetAnchorSequence());
									System.out.println("  -  HMMAlingment: "
											+ specNode.HMMAnnotation
													.GetMassShift());
									System.out.println("  -  HMMAlignment: "
											+ Utils.JoinIntArray(
													specNode.HMMAnnotation
															.GetShiftedPRMs(),
													" "));
									System.out.println("  -  HMMAlignment: "
											+ Utils.JoinIntArray(
													specNode.HMMAnnotation
															.getAlignedPeaks(),
													" "));
									System.out.println("  -  PathString: "
											+ specNode.HMMAnnotation
													.GetPathString());
								}
								System.out.println("  -  PRMs: "
										+ Utils.JoinIntArray(
												specNode.GetPRMSpectrum().PRMs,
												" "));

							}

							int[] shiftedPRMs = null;
							if (specNode.getModdedDBAnnotation() != null) {
								int[] alignedPeaks = AntibodyUtils
										.GetSeqPRMScaled(Utils
												.GetUnModded(specNode
														.getModdedDBAnnotation()));

								int massShift = currSeed.GetMassShift(
										alignedPeaks, this.PRMTolerance);
								shiftedPRMs = new int[specNode.GetPRMSpectrum().PRMs.length];
								for (int p = 0; p < shiftedPRMs.length; ++p) {
									shiftedPRMs[p] = specNode.GetPRMSpectrum().PRMs[p]
											+ massShift;
								}
							} else {
								// shiftedPRMs =
								// specNode.HMMAnnotation.GetShiftedPRMs();
								int[] alignedPeaks = specNode.HMMAnnotation
										.getAlignedPeaks();
								int massShiftA = specNode.HMMAnnotation
										.GetMassShift();
								for (int p = 0; p < alignedPeaks.length; ++p)
									alignedPeaks[p] = alignedPeaks[p]
											- massShiftA;

								// if(LocalDebug)
								// System.out.println("AlignedPeaks: " +
								// Utils.JoinIntArray(alignedPeaks, " "));
								int massShift = currSeed.GetMassShift(
										alignedPeaks, this.PRMTolerance);
								shiftedPRMs = new int[specNode.GetPRMSpectrum().PRMs.length];
								for (int p = 0; p < shiftedPRMs.length; ++p) {
									shiftedPRMs[p] = specNode.GetPRMSpectrum().PRMs[p]
											+ massShift;
								}
							}

							if (LocalDebug) {
								System.out.println("  -  Shifted?: "
										+ Utils.JoinIntArray(shiftedPRMs, " "));
							}

							int PRMIndex = AntibodyUtils.MassIndexClosest(
									shiftedPRMs, currSeed.SeedPRMs[j]);
							// The alignment doesn't contain this mass, it's
							// possible that this is a 2+ round alignment
							// so the shifted amount is not accurate.
							// if(PRMIndex < 0 && specNode.HMMAnnotation !=
							// null)
							// {
							// shiftedPRMs =
							// currSeed.GetShiftedPRMSpectrum(specNode,this.PRMTolerance);
							// if(LocalDebug)
							// {
							// System.out.println("  -  Shifted(2)?: " +
							// Utils.JoinIntArray(shiftedPRMs, " "));
							// }
							// PRMIndex =
							// AntibodyUtils.MassIndexClosest(shiftedPRMs,
							// currSeed.SeedPRMs[j], this.PRMTolerance);
							// }

							if (PRMIndex < 0) {
								// System.err.println("ERROR: Unable to find the allegedly supporting peak of "
								// +
								// currSeed.SeedPRMs[j] + " in spectrum:" +
								// specNode.GetScanNumber());
								// currSeed.DebugPrint();
								// System.err.println("ShiftedPRMs: " +
								// Utils.JoinIntArray(shiftedPRMs, " "));
								// specNode.DebugPrint();
								ErrorThrower.ThrowErrorCustum(
										ErrorThrower.CUSTOM_ERROR,
										"Unable to find the allegedly supporting peak of "
												+ currSeed.SeedPRMs[j]
												+ " in spectrum:"
												+ specNode.GetSpecIndex());

							}
							float PRM = ((float) (specNode.GetPRMSpectrum().PRMs[PRMIndex]))
									/ (float) AntibodyUtils.MASS_SCALE;
							if (LocalDebug)
								System.out.println("Spec Index: " + PRMIndex
										+ ", mass: " + PRM + ", is matched to "
										+ currSeed.SeedPRMs[j]);
							peakMasses.add(new Float(PRM));

							// if(ScanNumber.intValue() == 1231 && LocalDebug)
							// Utils.WaitForEnter();
							// else
							// LocalDebug = false;
						}

						// if (LocalDebug)
						// Utils.WaitForEnter();

					}
					// if (LocalDebug)
					// Utils.WaitForEnter();
					currSeed = currSeed.NextSeed;
					contigCount += 1;
				}
				currTemplate = currTemplate.NextInList;

			}
		}
		spectraCount = allSpectra.size();
		try {
			f = new DataOutputStream(new FileOutputStream(aBFileName));

			// numUsedSpectra [PROJ:int] (a)
			f.writeInt(Integer.reverseBytes(spectraCount));
			if (LocalDebug)
				System.out.println("spectraCount: " + spectraCount);

			// numUsedSpectra-by-1 unsigned int array of all used spectrum
			// indices (b)
			for (int i = 0; i < spectraCount; ++i) {
				f.writeInt(Integer.reverseBytes(((Integer) (allSpectra.get(i)))
						.intValue()));
				if (LocalDebug)
					System.out.println("spectrum[" + i + "]: "
							+ ((Integer) (allSpectra.get(i))).intValue());
			}

			// numUsedSpectra-by-2 unsigned short array of [PROJ:component
			// index,specFlipped] (c)
			for (int i = 0; i < spectraCount; ++i) {
				int[] specInfo = (int[]) (spectraInfo.get(i));
				f.writeShort(Short.reverseBytes((short) (specInfo[0])));
				f.writeShort(Short.reverseBytes((short) (specInfo[1])));
				if (LocalDebug)
					System.out.println("spectrum[" + i + "]: component: "
							+ specInfo[0] + ", flipped: " + specInfo[1]);
			}

			// numComponents [PROJ:int] (d)
			f.writeInt(Integer.reverseBytes(contigCount));
			if (LocalDebug)
				System.out.println("contigCount: " + contigCount);

			// numComponents-by-1 unsigned short array of number of ABruijn
			// vertices per component (e)
			for (int i = 0; i < peakCounts.size(); ++i) {
				f.writeShort(Short.reverseBytes(((Short) (peakCounts.get(i)))
						.shortValue()));
				if (LocalDebug)
					System.out.println("Component[" + i + "]: "
							+ ((Short) (peakCounts.get(i))).shortValue()
							+ " vertices");
			}

			// numABVertices-by-1 unsigned short array of number of spectrum
			// peaks per ABruijn vertex (f)
			for (int i = 0; i < ABVertexCounts.size(); ++i) {
				f.writeShort(Short.reverseBytes(((Short) (ABVertexCounts.get(i)))
						.shortValue()));
				if (LocalDebug)
					System.out.println("Vertex[" + i + "] has "
							+ (((Short) (ABVertexCounts.get(i))).shortValue())
							+ " spectra");
			}

			// totNumSpecPeaks-by-1 unsigned int array of spectrum index per
			// peak (g)
			// ArrayList specIndex = new ArrayList();
			// ArrayList peakMasses = new ArrayList();
			for (int i = 0; i < specIndex.size(); ++i) {
				f.writeInt(Integer.reverseBytes(((Integer) (specIndex.get(i)))
						.intValue()));
				if (LocalDebug)
					System.out.println("Peak[" + i + "] belongs to spectrum "
							+ ((Integer) (specIndex.get(i))).intValue()
							+ " and has mass "
							+ ((Float) (peakMasses.get(i))).floatValue());
			}

			// totNumSpecPeaks-by-1 float array of peak masses per spectrum peak
			// (h)
			for (int i = 0; i < peakMasses.size(); ++i) {
				// f.writeFloat(((Float)(peakMasses.get(i))).floatValue());
				f.writeInt(Integer.reverseBytes(Float
						.floatToIntBits(((Float) (peakMasses.get(i)))
								.floatValue())));
			}

			f.close();

		} catch (IOException E) {
			ErrorThrower.ThrowError(5, aBFileName);
		}

		return true;
	}

	/**
	 * Writes the 'contig' spectra from each seed. This is just the consensus
	 * prm spectrum that has been created from merging all of the overlapping
	 * spectra.
	 * 
	 * @param contigSpectrumFile
	 */
	private boolean writeSeedSpectraAsPKLBIN(String contigSpectrumFile) {

		return Template.writeAllSeedSpectraToPKLBIN(contigSpectrumFile,
				this.SelectedTemplates);

	}

	/**
	 * Writes the 'contig' spectra from each seed. This is just the consensus
	 * prm spectrum that has been created from merging all of the overlapping
	 * spectra.
	 * 
	 * @param contigSpectrumFile
	 */
	private boolean writeSeedSpectraAsMGF(String contigSpectrumFile) {

		return Template.writeAllSeedSpectraToMGF(contigSpectrumFile,
				this.SelectedTemplates);

	}

	/**
	 * Writes the contents of AllSpectra (PRM spectra only) to a pklbin format
	 * file
	 * 
	 * @param starSpectrumFile
	 *            The pklbin file to write to
	 * @return Returns true if write to file is successful, false otherwise
	 */
	private boolean writePRMSpectraAsPKLBIN(String starSpectrumFile) {
		return SpectrumNode.writeAllPRMSpectraToPKLBIN(starSpectrumFile,
				this.allPRMSpectra);
	}

	/**
	 * Writes the contents of AllSpectra (PRM spectra only) to an mgf format
	 * file
	 * 
	 * @param starSpectrumFile
	 *            The pklbin file to write to
	 * @return Returns true if write to file is successful, false otherwise
	 */
	private boolean writePRMSpectraAsMGF(String starSpectrumFile) {
		return SpectrumNode.writeAllPRMSpectraToMGF(starSpectrumFile,
				this.allPRMSpectra);
	}

	/*
	 * private void WriteSpreadsheetFile(String fileName) { FileWriter Writer =
	 * null; try { Writer = new FileWriter(fileName, false);
	 * 
	 * for (int i = 0; i < this.SelectedTemplates.length; ++i) { Template
	 * CurrTemplate = this.SelectedTemplates[i];
	 * Writer.write("\n# ---- DATA FOR SEQUENCE " + (i + 1) + "/" +
	 * this.SelectedTemplates.length + "----\n"); ArrayList seedStrings = new
	 * ArrayList(); ArrayList supportStrings = new ArrayList(); ArrayList
	 * overlapStrings = new ArrayList();
	 * 
	 * while (CurrTemplate != null) { Seed CurrSeed = CurrTemplate.anchors;
	 * 
	 * if (CurrTemplate.Sequence == null) { CurrTemplate =
	 * CurrTemplate.NextInList; continue; } String[] AAs = AntibodyUtils
	 * .SplitStringtoAA(CurrTemplate.Sequence); String TemplateString =
	 * Utils.JoinStringArray(AAs, "\t");
	 * 
	 * seedStrings.clear(); supportStrings.clear(); overlapStrings.clear(); int
	 * absoluteStartOfTemplate = 0;
	 * 
	 * String SeedString = ""; String SupportString = ""; String OverlapString =
	 * ""; // int CurrTemplatePos = 0;
	 * 
	 * // Iterate over each anchor // 1. Align it to the template // 2. Create
	 * the seed string while (CurrSeed != null) { SeedString = ""; SupportString
	 * = ""; OverlapString = ""; // Find the position of the anchor in the
	 * template int ExtensionInTemplate = CurrSeed
	 * .GetExtendedSeedPositionInTemplate(this.PRMTolerance);
	 * 
	 * // System.out.println("Considering seed : " + // CurrSeed.SeedSequence);
	 * // System.out.println("Alignment position " + // ExtensionInTemplate); if
	 * (ExtensionInTemplate == Integer.MIN_VALUE) { System.out
	 * .println("DARN! We couldn't find the base in the extended guy!!");
	 * CurrSeed.DebugPrint(); Utils.WaitForEnter(); CurrSeed =
	 * CurrSeed.NextSeed; continue; } // The extension goes beyond the template
	 * so add some // buffer to template string else if (ExtensionInTemplate <
	 * absoluteStartOfTemplate) { //
	 * System.out.println("Extension is < than start of Template " // +
	 * absoluteStartOfTemplate); if (ExtensionInTemplate <
	 * absoluteStartOfTemplate) absoluteStartOfTemplate = ExtensionInTemplate;
	 * while (ExtensionInTemplate < absoluteStartOfTemplate) { // Add buffer to
	 * template TemplateString = " \t" + TemplateString; // Add buffer to any
	 * previous seed that needs it for (int k = 0; k < seedStrings.size(); ++k)
	 * { String tempString = (String) (seedStrings .remove(k)); tempString =
	 * "*\t" + tempString; seedStrings.add(k, tempString);
	 * 
	 * tempString = (String) (supportStrings .remove(k)); tempString = "*\t" +
	 * tempString; supportStrings.add(k, tempString);
	 * 
	 * tempString = (String) (overlapStrings .remove(k)); tempString = "*\t" +
	 * tempString; overlapStrings.add(k, tempString); } ExtensionInTemplate +=
	 * 1; } } // This is deprecated since each anchor is now given its // own
	 * line (like a contig) // else if(ExtensionInTemplate - CurrTemplatePos <
	 * 0) // //The seeds overlap so we need to put it in a new // line // { //
	 * SeedString += "\nPredicted\t"; // SupportString += "\nSupport\t"; //
	 * OverlapString += "\nOverlap\t"; // CurrTemplatePos =
	 * absoluteStartOfTemplate; // }
	 * 
	 * while (ExtensionInTemplate - absoluteStartOfTemplate > 0) { //
	 * System.out.println("Extension " + // ExtensionInTemplate + " is beyond "
	 * + // absoluteStartOfTemplate); SeedString += "*\t"; SupportString +=
	 * "*\t"; OverlapString += "*\t"; ExtensionInTemplate -= 1; } String[]
	 * SeedAAs = AntibodyUtils .SplitStringtoAA(CurrSeed.SeedSequence); int[]
	 * OverlapArray = CurrSeed .GetSiteOverlap(this.PRMTolerance); int[]
	 * SupportArray = CurrSeed .GetSiteSupport(this.PRMTolerance); SeedString +=
	 * SeedAAs[0]; SupportString += SupportArray[0]; OverlapString +=
	 * OverlapArray[0]; for (int j = 1; j < SeedAAs.length; ++j) { SeedString +=
	 * "\t" + SeedAAs[j]; SupportString += "\t" + SupportArray[j]; OverlapString
	 * += "\t" + OverlapArray[j]; } // System.out.println("Seed: " +
	 * SeedString); // System.out.println("Support: " + SupportString); //
	 * System.out.println("Overlap: " + OverlapString); // Utils.WaitForEnter();
	 * seedStrings.add(SeedString); supportStrings.add(SupportString);
	 * overlapStrings.add(OverlapString); CurrSeed = CurrSeed.NextSeed; }
	 * 
	 * // Write template and seed info Writer.write(TemplateString + "\n");
	 * Writer.write(SeedString + "\n"); Writer.write(SupportString + "\n");
	 * Writer.write(OverlapString + "\n\n");
	 * 
	 * Writer.write(CurrTemplate.TemplateName + "\t");
	 * 
	 * int StrPos = 0; // Keeps track of the position in the // template string
	 * int SeedStrPos = 0; // Keeps track of the position in the // seed string
	 * int SupportStrPos = 0; int OverlapStrPos = 0;
	 * 
	 * int WrittenCols = 1; while (StrPos < TemplateString.length()) { if
	 * (TemplateString.charAt(StrPos) != '\t' && TemplateString.charAt(StrPos)
	 * != '\n') WrittenCols += 1; Writer.write(TemplateString.charAt(StrPos));
	 * StrPos += 1; if (WrittenCols == this.MaxXLSWidth || StrPos ==
	 * TemplateString.length()) // We // wrote // a // full // Template // line
	 * // so // write // a // seed // line { for (int k = 0; k <
	 * seedStrings.size(); ++k) { SeedString = (String) (seedStrings.get(k));
	 * SupportString = (String) (supportStrings.get(k)); OverlapString =
	 * (String) (overlapStrings.get(k)); Writer.write("\nPredicted (" + k +
	 * ")"); WrittenCols = 1; if (SeedStrPos < SeedString.length() &&
	 * SeedString.charAt(SeedStrPos) != '\t') Writer.write("\t");
	 * 
	 * while (WrittenCols < this.MaxXLSWidth && SeedStrPos <
	 * SeedString.length()) { if (SeedString.charAt(SeedStrPos) != '\t')
	 * WrittenCols += 1;
	 * 
	 * if (SeedString.charAt(SeedStrPos) == '[') // This // is // the //
	 * beginning // of // a // mass // gap { while
	 * (SeedString.charAt(SeedStrPos) != ']') { Writer.write(SeedString
	 * .charAt(SeedStrPos)); SeedStrPos += 1; } }
	 * Writer.write(SeedString.charAt(SeedStrPos)); SeedStrPos += 1; }
	 * Writer.write("\nSupport (" + k + ")"); WrittenCols = 1; if (SupportStrPos
	 * < SupportString.length() && SupportString.charAt(SupportStrPos) != '\t')
	 * Writer.write("\t"); while (WrittenCols < this.MaxXLSWidth &&
	 * SupportStrPos < SupportString .length()) { if
	 * (SupportString.charAt(SupportStrPos) != '\t') WrittenCols += 1;
	 * Writer.write(SupportString .charAt(SupportStrPos)); SupportStrPos += 1; }
	 * 
	 * Writer.write("\nOverlap (" + k + ")"); WrittenCols = 1; if (OverlapStrPos
	 * < OverlapString.length() && OverlapString.charAt(OverlapStrPos) != '\t')
	 * Writer.write("\t"); while (WrittenCols < this.MaxXLSWidth &&
	 * OverlapStrPos < OverlapString .length()) { if
	 * (OverlapString.charAt(OverlapStrPos) != '\t') WrittenCols += 1;
	 * Writer.write(OverlapString .charAt(OverlapStrPos)); OverlapStrPos += 1; }
	 * }
	 * 
	 * WrittenCols = 1; Writer.write("\n"); if (StrPos <
	 * TemplateString.length()) Writer.write("\n" + CurrTemplate.TemplateName);
	 * } } Writer.write("\n");
	 * 
	 * CurrTemplate = CurrTemplate.NextInList; }
	 * 
	 * } Writer.close(); } catch (Exception E) { ErrorThrower.ThrowError(5,
	 * fileName); } }
	 */
	/**
	 * Reconstruct the path for a given set of selected templates
	 * 
	 * @param i
	 * @return
	 */
	private String[] ReconstructPath(int TemplateSet) {
		// int Count = 0;
		String[] Ret = new String[3];
		Template Ptr = this.SelectedTemplates[TemplateSet];
		while (Ptr != null) {
			Seed Head = Ptr.anchors;
			while (Head != null) {
				// Count += 1;
				Head = Head.NextSeed;

			}
			// if(this.DivergenceTestFlag)
			// break;
			Ptr = Ptr.NextInList;
		}

		String FinalPath = "";
		String Overlap = "";
		String Support = "";
		// int ClassNum = 0;
		Seed Head = this.SelectedTemplates[TemplateSet].anchors;

		while (Head != null) {
			// Head.DebugPrint();
			if (Head.HasSequenceOverlapWithNext(this.SelectedTemplates,
					this.PRMTolerance)) {
				if (!Head.MergeSequenceWithNext(this.SelectedTemplates,
						this.PRMTolerance)) {
					ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
							"ERROR: Error merging seed with next");
				}

			} else {
				if (Head == null)
					System.out.println("WTF!!");
				FinalPath += " ... " + Head.SeedSequence;

				Overlap += " ... "
						+ Utils.IntArrayToString(Head
								.GetSiteOverlap(this.PRMTolerance));
				Support += " ... "
						+ Utils.IntArrayToString(Head
								.GetSiteSupport(this.PRMTolerance));
				Head = Head.GetNextSeed(this.SelectedTemplates);
			}
			// AntibodyUtils.WaitForEnter();
		}

		Ret[0] = FinalPath.substring(5);
		Ret[1] = Support.substring(5);
		Ret[2] = Overlap.substring(5);
		return Ret;
	}

	private void WriteOutputHeaderInfo(String OutputFile) {

		FileWriter Writer = null;
		String FinalOutputFile = OutputFile;
		if (OutputFile == null)
			FinalOutputFile = this.detailsOutputFileName;

		// Attempt to Write out the parameters
		try {
			Writer = new FileWriter(FinalOutputFile, false);
			Writer.write("#***GenoMS Version 2.0***\n");
			Writer.write("#Config File: " + this.configFileName + "\n");
			Writer.write("#Spectrum Files:\n");
			for (int i = 0; i < SpectraFileNames.length; ++i)
				Writer.write("#" + SpectraFileNames[i] + "\n");

			Writer.write("#MSGFDB Data (" + this.annCount + ")\n");

			if (this.DatabaseRootName != null)
				Writer.write("#DatabaseRootName: " + this.DatabaseRootName
						+ "\n");
			if (this.CombinedShuffledDBName != null)
				Writer.write("#ShuffledDBName: " + this.CombinedShuffledDBName
						+ "\n");

			Writer.write("#PRMFileNames (" + this.allPRMSpectra.length + "):\n");
			for (int i = 0; i < this.PRMFileNames.length; ++i)
				Writer.write("#" + this.PRMFileNames[i] + "\n");
			Writer.write("#Parameters:\n");
			Writer.write("# Min Template Likelihood: "
					+ AntibodyUtils.MIN_TEMPLATE_LIKELIHOOD + "\n");
			Writer.write("# SpecAlignment SCORE Cutoff: "
					+ AntibodyUtils.MSALIGN_SCORE_CUTOFF + "\n");
			Writer.write("# Min peaks spec aligned: "
					+ AntibodyUtils.MIN_OVERLAPPING_PEAKS + "\n");
			Writer.write("# MaxAlignedSpectra kept per Seed: "
					+ AntibodyUtils.MAX_ALIGNED_SPECTRA + "\n");
			Writer.write("# Min HMM Match State Score: "
					+ AntibodyUtils.MIN_PEAK_GROUP_SCORE + "\n");
			Writer.write("# Min Votes to consider a new match: "
					+ AntibodyUtils.MIN_PEAK_GROUP_SIZE + "\n");
			Writer.write("# HMM Score worsening multiplier: "
					+ AntibodyUtils.HMM_SCORE_SCALE + "\n");
			if (this.ClusterFlag)
				Writer.write("# All2All Spectra Similarity in HMM: "
						+ AntibodyUtils.MSALIGN_MIN_SIMILARITY + "\n");
			// Writer.write("# Min Spectral Similarity: " +
			// Utils.MSALIGN_MIN_SIMILARITY + "\n");
			if (this.applyMatchPenalty)
				Writer.write("# HMM Match Miss Penalty: "
						+ this.matchPenalty + "\n");

			for (int i = 0; this.FixedMods != null && i < this.FixedMods.length; ++i)
				if (this.FixedMods[i] != 0)
					Writer.write("# Fixed Mod on " + (char) (i + 65) + ":"
							+ this.FixedMods[i] + "\n");
			Writer.close();
		} catch (Exception E) {
			ErrorThrower.ThrowError(5, FinalOutputFile);
		}

	}

	/*
	 * private Seed ExtendSeedLeft(String OutputFile, Seed Head, int NumRounds)
	 * { String FinalOutputFile = OutputFile; if (OutputFile == null)
	 * FinalOutputFile = this.detailsOutputFileName; MSAlignmentType[]
	 * Candidates = null; // if(Debug) // Head.DebugPrint(); FileWriter Writer =
	 * null; // EXTEND THE SEED LEFT if (Debug)
	 * System.out.println("\n--- Extending Seed " + Head.GetSeedSequence() +
	 * " left ---"); int Round = 0; do { if (Debug)
	 * System.out.println("[Round: " + Round + "]"); //
	 * System.out.println("CAN WE MERGE WITH PREVSEED?"); if
	 * (Head.HasSequenceOverlapWithPrev(this.SelectedTemplates,
	 * this.PRMTolerance)) { // System.out.println("Yes!!!"); return Head; }
	 * 
	 * try { Writer = new FileWriter(FinalOutputFile, true);
	 * Writer.write("\n[LEFT EXTENSION] Round: " + Round + " for seed " +
	 * Head.GetSeedSequence() + "\n"); Writer.close(); } catch (Exception E) {
	 * ErrorThrower.ThrowError(5, FinalOutputFile); } Round++;
	 * 
	 * Candidates = this.GatherSpectraLeft(Head, this.allPRMSpectra);
	 * 
	 * int[] ScaledConsensusPeaks = null; double[] PeakScores = null;
	 * ConsensusPRMSpectrum Consensus = null; if (Candidates != null) { try {
	 * Writer = new FileWriter(FinalOutputFile, true);
	 * 
	 * for (int c = 0; c < Candidates.length; ++c) {
	 * 
	 * Writer.write("[" + c + "]:" + Candidates[c].toStringDebug() + "\n"); }
	 * Writer.close(); } catch (Exception E) { ErrorThrower.ThrowError(5,
	 * FinalOutputFile); } MSAlignmentType[] ReducedCandidates = Candidates;
	 * 
	 * Consensus = this.PerformAlignmentWithScoresLeft( ReducedCandidates,
	 * Head); ScaledConsensusPeaks = Consensus.ScaledPeaks; PeakScores =
	 * Consensus.PeakScores; try { Writer = new FileWriter(FinalOutputFile,
	 * true); Writer.write("Consensus:"); for (int k = 0; k <
	 * ScaledConsensusPeaks.length; ++k) Writer.write(" " +
	 * ScaledConsensusPeaks[k]); Writer.write("\n"); } catch (Exception E) {
	 * ErrorThrower.ThrowError(5, FinalOutputFile); }
	 * 
	 * // if (Debug) // System.out.println("Finished gathering ConsensusPeaks");
	 * 
	 * String Peptide = DeNovo.Reconstruct(ScaledConsensusPeaks, PeakScores,
	 * this.PRMTolerance);
	 * 
	 * // // System.out.println("De Novo: " + Peptide);
	 * 
	 * ConsensusPRMSpectrum MergedPRMs = ConsensusPRMSpectrum
	 * .MergePRMs(ScaledConsensusPeaks, PeakScores,
	 * Consensus.OverlappingSpectra, Consensus.SupportingSpectra,
	 * Head.GetSeedPRMs(), Head.GetSeedPRMScores(), Head.PeakOverlapping,
	 * Head.PeakSupporting, this.PRMTolerance); String Merged = ""; System.out
	 * .println("* Step 4: Merge the consensus with the current seed sequence");
	 * if (MergedPRMs != null) { Merged =
	 * DeNovo.Reconstruct(MergedPRMs.ScaledPeaks, MergedPRMs.PeakScores,
	 * this.PRMTolerance);
	 * 
	 * System.out.println("Extended Seed Sequence: " + Merged); }
	 * 
	 * if (Merged.length() == 0) {
	 * 
	 * if (Debug) System.out
	 * .println("Extension is not valid, not merging with larger seed!");
	 * 
	 * try { Writer.write("De Novo: " + Peptide + ", NOT MERGEABLE!!\n\n");
	 * Writer.close(); } catch (Exception E) { ErrorThrower.ThrowError(7,
	 * FinalOutputFile); } break; } else { String MergedPRMStr = ""; String
	 * MergedScoreStr = ""; for (int r = 0; r < MergedPRMs.ScaledPeaks.length;
	 * ++r) { MergedPRMStr += MergedPRMs.ScaledPeaks[r] + " "; MergedScoreStr +=
	 * MergedPRMs.PeakScores[r] + " "; } try { Writer.write("De Novo: " +
	 * Peptide + ", Merged: " + Merged + "\n"); Writer.write(MergedPRMStr +
	 * "\n"); Writer.write(MergedScoreStr + "\n"); Writer.close(); } catch
	 * (Exception E) { ErrorThrower.ThrowError(7, FinalOutputFile); } } if
	 * (AntibodyUtils.CompareToPRMs(Merged, Head.GetSeedSequence(),
	 * this.PRMTolerance) == 0) { if (Debug) System.out
	 * .println("Merged sequence is same as seed sequence, finishing seeds");
	 * break; } // Head.SetSeedSequence(Merged); //
	 * Head.SetSeedPRMs(MergedPRMs); // Head.SetSeedPRMs(MergedPRMs.ScaledPeaks,
	 * // MergedPRMs.PeakScores); double[][] UsedPRMs =
	 * DeNovo.ReconstructPRMsAndScores( MergedPRMs.ScaledPeaks,
	 * MergedPRMs.PeakScores, this.PRMTolerance); double[] Scores = new
	 * double[UsedPRMs.length]; int[] PRMs = new int[UsedPRMs.length]; for (int
	 * p = 0; p < Scores.length; ++p) { Scores[p] = UsedPRMs[p][1]; PRMs[p] =
	 * (int) (UsedPRMs[p][0]); } if (Scores.length != PRMs.length) { System.out
	 * .println("ERROR: Length of Score array != PRM array!! " + Scores.length +
	 * " != " + PRMs.length); Utils.WaitForEnter(); } Object[] Overlapping =
	 * MergedPRMs.GetOverlappingSpectra(PRMs); Object[] Supporting =
	 * MergedPRMs.GetSupportingSpectra(PRMs); //
	 * System.out.println("----Finished Round of Extension----"); //
	 * System.out.println(DeNovo.Reconstruct(PRMs, Scores)); //
	 * System.out.println(AntibodyUtils.DoubleArrayToString(IntervalScores)); //
	 * AntibodyUtils.WaitForEnter(); Head.SetSeedPRMs(PRMs, Scores, Overlapping,
	 * Supporting, this.PRMTolerance);
	 * 
	 * } } while (Candidates != null && Round < NumRounds); return Head; }
	 */
	/*
	 * private Seed ExtendSeedRight(String OutputFile, Seed Head, int NumRounds)
	 * { String FinalOutputFile = OutputFile; if (OutputFile == null)
	 * FinalOutputFile = this.detailsOutputFileName; MSAlignmentType[]
	 * Candidates = null; FileWriter Writer = null; // EXTEND THE SEED Right if
	 * (Debug) System.out.println("\n--- Extending Seed " +
	 * Head.GetSeedSequence() + " right ---"); int Round = 0; do { if (Debug)
	 * System.out.println("[Round: " + Round + "]");
	 * 
	 * try { Writer = new FileWriter(FinalOutputFile, true);
	 * Writer.write("\n[RIGHT EXTENSION] Round: " + Round + " for seed " +
	 * Head.GetSeedSequence() + "\n"); Writer.close(); } catch (Exception E) {
	 * ErrorThrower.ThrowError(5, FinalOutputFile); } Round++;
	 * 
	 * Candidates = this.GatherSpectraRight(Head, this.allPRMSpectra);
	 * 
	 * int[] ScaledConsensusPeaks = null; double[] PeakScores = null;
	 * ConsensusPRMSpectrum Consensus = null; if (Candidates != null) {
	 * 
	 * try { Writer = new FileWriter(FinalOutputFile, true);
	 * 
	 * for (int c = 0; c < Candidates.length; ++c) { Writer.write("[" + c + "]:"
	 * + Candidates[c].toStringDebug() + "\n"); } Writer.close(); } catch
	 * (Exception E) { ErrorThrower.ThrowError(7, FinalOutputFile); }
	 * 
	 * // MSAlignmentType [] ReducedCandidates = //
	 * MSAlignmentType.ClusterRankAlignments(Candidates); MSAlignmentType[]
	 * ReducedCandidates = Candidates;
	 * 
	 * Consensus = this.PerformAlignmentWithScoresRight( ReducedCandidates,
	 * Head); ScaledConsensusPeaks = Consensus.ScaledPeaks; PeakScores =
	 * Consensus.PeakScores; try { Writer = new FileWriter(this.outputFileName,
	 * true); Writer.write("Consensus:"); for (int k = 0; k <
	 * ScaledConsensusPeaks.length; ++k) Writer.write(" " +
	 * ScaledConsensusPeaks[k]); Writer.write("\n"); } catch (Exception E) {
	 * ErrorThrower.ThrowError(5, this.outputFileName); }
	 * 
	 * String Peptide = DeNovo.Reconstruct(ScaledConsensusPeaks, PeakScores,
	 * this.PRMTolerance); // int[] UsedPeaks = //
	 * DeNovo.ReconstructPRMs(ScaledConsensusPeaks,PeakScores); // double[]
	 * PeptideScores = // Consensus.Model.ComputeProbabilityTable(UsedPeaks);
	 * 
	 * // if (Debug) // System.out.println("De Novo: " + Peptide);
	 * 
	 * 
	 * ConsensusPRMSpectrum MergedPRMs = ConsensusPRMSpectrum
	 * .MergePRMs(Head.GetSeedPRMs(), Head.GetSeedPRMScores(),
	 * Head.PeakOverlapping, Head.PeakSupporting, ScaledConsensusPeaks,
	 * PeakScores, Consensus.OverlappingSpectra, Consensus.SupportingSpectra,
	 * this.PRMTolerance); String Merged = ""; System.out
	 * .println("* Step 4: Merge the consensus with the current seed sequence");
	 * if (MergedPRMs != null) { Merged =
	 * DeNovo.Reconstruct(MergedPRMs.ScaledPeaks, MergedPRMs.PeakScores,
	 * this.PRMTolerance);
	 * 
	 * System.out.println("Extended Seed Sequence: " + Merged);
	 * 
	 * } // }
	 * 
	 * if (Merged.length() == 0) {
	 * 
	 * if (Debug) System.out
	 * .println("Extension is not valid, not merging with larger seed!");
	 * 
	 * try { Writer.write("De Novo: " + Peptide + ", NOT MERGEABLE!!\n\n");
	 * Writer.close(); } catch (Exception E) { ErrorThrower.ThrowError(7,
	 * this.outputFileName); } break; } else { String MergedPRMStr = ""; String
	 * MergedScoreStr = ""; for (int r = 0; r < MergedPRMs.ScaledPeaks.length;
	 * ++r) { MergedPRMStr += MergedPRMs.ScaledPeaks[r] + " "; MergedScoreStr +=
	 * MergedPRMs.PeakScores[r] + " "; } try { Writer.write("De Novo: " +
	 * Peptide + ", Merged: " + Merged + "\n"); Writer.write(MergedPRMStr +
	 * "\n"); Writer.write(MergedScoreStr + "\n"); Writer.close(); } catch
	 * (Exception E) { ErrorThrower.ThrowError(7, this.outputFileName); } } if
	 * (AntibodyUtils.CompareToPRMs(Merged, Head.GetSeedSequence(),
	 * this.PRMTolerance) == 0) {
	 * 
	 * if (Debug) System.out
	 * .println("Merged sequence is same as seed sequence, finishing seeds");
	 * break; } double[][] UsedPRMs = DeNovo.ReconstructPRMsAndScores(
	 * MergedPRMs.ScaledPeaks, MergedPRMs.PeakScores, this.PRMTolerance);
	 * double[] Scores = new double[UsedPRMs.length]; int[] PRMs = new
	 * int[UsedPRMs.length]; for (int p = 0; p < Scores.length; ++p) { Scores[p]
	 * = UsedPRMs[p][1]; PRMs[p] = (int) (UsedPRMs[p][0]); } if (Scores.length
	 * != PRMs.length) { System.out
	 * .println("ERROR: Length of Score array != PRM array!! " + Scores.length +
	 * " != " + PRMs.length); Utils.WaitForEnter(); } Object[] Overlapping =
	 * MergedPRMs.GetOverlappingSpectra(PRMs); Object[] Supporting =
	 * MergedPRMs.GetSupportingSpectra(PRMs); //
	 * System.out.println("----Finished Round of Extension----"); //
	 * System.out.println(DeNovo.Reconstruct(PRMs, Scores)); //
	 * System.out.println(AntibodyUtils.DoubleArrayToString(IntervalScores)); //
	 * AntibodyUtils.WaitForEnter(); Head.SetSeedPRMs(PRMs, Scores, Overlapping,
	 * Supporting, this.PRMTolerance);
	 * 
	 * // See if we can merge this anchor with the next one
	 * 
	 * if (Head.HasSequenceOverlapWithNext(this.SelectedTemplates,
	 * this.PRMTolerance)) { return Head; } } } while (Candidates != null &&
	 * Round < NumRounds); return Head; }
	 */
	private void BuildAndExtendSeeds(String OutputFile, int NumRounds)
			throws IOException {
		// AddDBSearchResultsToSpectrumNodes();

		System.out.println("\n[GenoMS Part 7]: Attempt seed extension");
		for (int i = 0; i < this.SelectedTemplates.length; ++i) {
			Template Ptr = this.SelectedTemplates[i];
			while (Ptr != null) {
				if (Ptr.anchors == null) {
					Ptr = Ptr.NextInList;
					continue;
				}
				// System.out.println("**Performing extensions for template: "
				// + Ptr.ProteinID);
				Seed Head = Ptr.anchors;
				while (Head != null) {

					// Head.ReconcileSequence();
					// this.ExtendSeedLeft(OutputFile, Head, NumRounds);

					// this.ExtendSeedRight(OutputFile, Head, NumRounds);

					if (Debug) {
						Head.DebugPrint();
						// AntibodyUtils.GetFullSequenceCorrectness(Head.SeedSequence);
						// AntibodyUtils.WaitForEnter();
					}
					Head = Head.NextSeed;
				}
				Ptr = Ptr.NextInList;
			}
		}
	}

	private void WriteSeedInfoToFile(String FileName) throws IOException {

		FileWriter Writer = null;
		try {
			Writer = new FileWriter(FileName, true);
			Writer.write("# ---- FINAL RECONSTRUCTION INFO ----\n");

			for (int i = 0; i < this.SelectedTemplates.length; ++i) {
				int TemplateCount = 0;
				Template Ptr = this.SelectedTemplates[i];
				Writer.write("\n# ---- DATA FOR SEQUENCE " + i + "/"
						+ this.SelectedTemplates.length + "----\n");
				while (Ptr != null) {
					Seed CurrSeed = Ptr.anchors;

					int SeedCount = 0;
					if (CurrSeed != null) {
						Writer.write("# Template " + TemplateCount + ":"
								+ Ptr.TemplateName + " seeds: \n");
						Writer.write("# Template: " + Ptr.Sequence + "\n");
					}
					while (CurrSeed != null) {
						Writer.write("# Seed[" + TemplateCount + "]["
								+ SeedCount + "]: "
								+ CurrSeed.AnchoredSeedSequence + " -> "
								+ CurrSeed.SeedSequence + "\n");
						Hashtable<String, Integer> Peps = new Hashtable<String, Integer>();
						for (int k = 0; k < CurrSeed.SupportingSpectra.size(); ++k) {
							SpectrumNode CurrAnn = (SpectrumNode) (CurrSeed.SupportingSpectra
									.get(k));
							Writer.write("  [" + k + "]: "
									+ CurrAnn.getModdedDBAnnotation() + "\n");

							String CurrUnModdedAnn = Utils.GetUnModded(CurrAnn
									.getModdedDBAnnotation());
							Peps.put(CurrUnModdedAnn, new Integer(1));
						}
						Writer.write("  Peptides: "
								+ AntibodyUtils.GetHashtableStringKeys(Peps,
										",") + "\n");
						Writer.write("  Total Peptides: " + Peps.size() + "\n");
						CurrSeed = CurrSeed.NextSeed;
						SeedCount++;
					}
					Ptr = Ptr.NextInList;
					TemplateCount += 1;
				}

				Hashtable<String, Integer> Peptides = new Hashtable<String, Integer>();
				for (int j = 0; j < this.allPRMSpectra.length; ++j) {
					String CurrUnModdedAnn = Utils
							.GetUnModded(this.allPRMSpectra[j]
									.getModdedDBAnnotation());
					if (CurrUnModdedAnn == null)
						continue;
					Peptides.put(CurrUnModdedAnn, new Integer(1));
				}
				Writer.write("# Unique Peptides: " + Peptides.size() + "\n");
				String[] Final = this.ReconstructPath(i);

				Writer.write("# ----FINAL SEQUENCE INFO----\n");
				Writer.write("# Final Seq: " + Final[0] + "\n");
				Writer.write("# Sitewise Support: " + Final[1] + "\n");
				Writer.write("# Sitewise Overlap: " + Final[2] + "\n");
				Final[0] = Final[0].replace("-", "");
				String[] AAs = AntibodyUtils.SplitStringtoAA(Final[0]);

				Writer.write("# FinalLen: " + AAs.length + "\n");

			}
			Writer.close();
		} catch (Exception E) {
			ErrorThrower.ThrowError(5, FileName);
		}
	}

	public void redirectSystemOut(String logFileName) {

		try {

			System.setOut(new PrintStream(new FileOutputStream(logFileName)));

		} catch (FileNotFoundException ex) {
			ErrorThrower.ThrowError(5, logFileName);
		}
	}

	public static void main(String[] args) throws IOException {
		// Parse Command Line Stuff
		// System.exit(26);
		Date Start = new Date();
		long StartTime = Start.getTime();

		String[] Commands = { "-i", "-o", "-r", "-x", "-p", "-s", "-a", "-w",
				"-f", "-q", "-l", "-k", "-e" };
		boolean[] Values = { true, true, true, false, false, false, false,
				true, false, false, true, true, true };
		Hashtable<String,String> Options = Utils.ParseCommandLine(args, Commands, Values);

		// Some checking for all parameters
		if (!Options.containsKey("-i") || !Options.containsKey("-o")) {
			System.err.println(UsageInfo);
			ErrorThrower.ThrowError(2,
					"Must specify both an input and an output file");
		}
		
		// Initialize Driver
				AntibodyDriver Driver = new AntibodyDriver();
		if (Options.containsKey("-p"))
			Driver.applyMatchPenalty = true;

		
		// Driver.isPC = Utils.IsRunningWindows();
		Driver.SetConfigFileName((String) (Options.get("-i")));
		Driver.SetOutputFileName((String) (Options.get("-o")));
		Driver.detailsOutputFileName = Utils
				.GetFileNameNoExtension(Driver.outputFileName) + ".Details.txt";
		if (Options.containsKey("-s"))
			Driver.ClusterFlag = true;

		File Temp;
		if (Options.containsKey("-r")) {
			Driver.SetResourceDir((String) (Options.get("-r")));
			Temp = new File(Driver.ResourceDir);
			Driver.ResourceDir = Temp.getCanonicalPath();
			// System.out.println("Old resource dir: " +
			// System.getProperty("user.dir"));
			// System.setProperty("user.dir", Driver.ResourceDir);
			// System.out.println("New resource dir: " +
			// System.getProperty("user.dir"));
		} else {
			Driver.SetResourceDir(System.getProperty("user.dir"));
			if (!Utils.IsDir(Driver.ResourceDir))
				Utils.MakeDir(Driver.ResourceDir);
		}

		if (Options.containsKey("-k")) {
			Driver.projectDir = ((String) (Options.get("-k")));
			Temp = new File(Driver.projectDir);
			Driver.projectDir = Temp.getCanonicalPath();
			// System.out.println("Old resource dir: " +
			// System.getProperty("user.dir"));
			// System.setProperty("user.dir", Driver.ResourceDir);
			// System.out.println("New resource dir: " +
			// System.getProperty("user.dir"));
		} else {
			Driver.projectDir = (System.getProperty("user.dir"));
		}

		if (Options.containsKey("-e")) {
			Driver.exeDir = ((String) (Options.get("-e")));
			Temp = new File(Driver.exeDir);
			Driver.exeDir = Temp.getCanonicalPath();
			// System.out.println("Old resource dir: " +
			// System.getProperty("user.dir"));
			// System.setProperty("user.dir", Driver.ResourceDir);
			// System.out.println("New resource dir: " +
			// System.getProperty("user.dir"));
		} else {
			Driver.exeDir = (System.getProperty("user.dir"));
		}

		if (Options.containsKey("-q")) {
			Driver.reportFileDir = Driver.projectDir + File.separator;

			if (!Utils.IsDir(Driver.reportFileDir))
				Utils.MakeDir(Driver.reportFileDir);
		}

		if (Options.containsKey("-x"))
			Driver.runDBSearchFlag = true;
		if (Options.containsKey("-a"))
			Driver.AddMatchStatesAtEndFlag = true;

		if (Options.containsKey("-w")) {
			double D = Double.parseDouble((String) (Options.get("-w")));
			if (D <= 0 || D > 1.0) {
				System.err.println("WARNING: PValue cutoff " + D
						+ " is invalid, using default " + Driver.fdrCutoff);
			} else
				Driver.fdrCutoff = D;

		}

		if (Options.containsKey("-l")) {
			Driver.redirectSystemOut(AntibodyUtils.makePathAbsolute(
					(String) (Options.get("-l")), Driver.projectDir));

		}
		// Make sure we can find the config file
		Driver.configFileName = AntibodyUtils.makePathAbsolute(
				Driver.configFileName, Driver.projectDir);
		Driver.outputFileName = AntibodyUtils.makePathAbsolute(
				Driver.outputFileName, Driver.projectDir);
		Driver.detailsOutputFileName = AntibodyUtils.makePathAbsolute(
				Driver.detailsOutputFileName, Driver.projectDir);

		// if(Options.containsKey("-c"))
		// Driver.aBTLACoverage =
		// Double.parseDouble((String)(Options.get("-c")));
		System.out.println("---- GenoMS Version " + AntibodyUtils.VersionInfo
				+ " ----\n");
		System.out.println("PARAMETERS:");
		System.out.println("Config File Name: " + Driver.GetConfigFileName());
		System.out.println("Output FileName: " + Driver.GetOutputFileName());
		System.out.println("Resource Dir: " + Driver.GetResourceDir());
		System.out.flush();

		// if(Driver.FreeModsFlag && Driver.PathtoMSGF != null)
		// {
		// System.err.println("ERROR: Searching for mutations is currently incompatible with MSGF");
		// return;
		// }

		Driver.RunDriver();

		Date End = new Date();
		long EndTime = End.getTime();

		System.out.println("Elapsed Time: " + (EndTime - StartTime)
				+ " ms, aka " + ((EndTime - StartTime) / 1000) + " s");
		System.out.flush();
		System.exit(0);
	}

}
