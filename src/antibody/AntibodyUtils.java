package antibody;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import trieUtils.TrieDB;
import basicUtils.Utils;
import enosi.db.SixFrameBuilder;
import errorUtils.ErrorThrower;
//import basicUtils.Utils;

/**
 * 
 * @author natalie
 * @version 12.16.2008
 */
public class AntibodyUtils {

	public static final String VersionInfo = "2015.11.05.1";

	/**
	 * These are the sequences for the antibody aBTLA, the OracleSequence is
	 * what we are trying to figure out. It's included here just for reference
	 * and should never be used except in 'oracle' modes for some methods.
	 */
	public static final String aBTLA_OracleSequence = "QVQLKESGPGLVAPSQSLSITCTVSGFSLTSYGVSWVRQPPGKGLEWLGVIWGDGSTNYHSALISRLSISKDNSKSQVFLKLNSLQTDDTATYYCAKGGYRFYYAMDYWGQGTSVTVSSAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSQTVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVLTITLTPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTQPREEQFNSTFRSVSELPIMHQDWLNGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTIPPPKEQMAKDKVSLTCMITDFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKLNVQKSNWEAGNTFTCSVLHEGLHNHHTEKSLSHSPGK";
	public static final String aBTLA_OracleSequence_OLD = "QVQLKESGPGLVAPSQSLSITCTVSGFSLTSYGVSWVRQPPGKGLEWLGVIWGDGSTNYHSALISRLSISKDNSKSQVFLKLNSLQTDDTATYYCAKGGYRFYYAMDYWGQGTSVTVSSAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVMSSPRPSETVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVLTITLTPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTQPREEQFNSTFRSVSELPIMHQDWLNGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTIPPPKEQMAKDKVSLTCMITDFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKLNVQKSNWEAGNTFTCSVLHEGLHNHHTEKSLSHSPGK";

	public static final String aBTLA_VRegion = "QVQLKESGPGLVAPSQSLSITCTVSGFSLTSYGVSWVRQPPGKGLEWLGVIWGDGSTNYHSALISRLSISKDNSKSQVFLKLNSLQTDDTATYYCAK";
	public static final String aBTLA_DRegion = "GGYRF";
	public static final String aBTLA_JRegion = "YYAMDYWGQGTSVTVSS";
	public static final String aBTLA_CRegion = "AKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVMSSPRPSETVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVLTITLTPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTQPREEQFNSTFRSVSELPIMHQDWLNGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTIPPPKEQMAKDKVSLTCMITDFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKLNVQKSNWEAGNTFTCSVLHEGLHNHHTEKSLSHSPGK";

	public static final String aBTLA_LC = "DIVMSQSPSSLAVSAGEKVTMSCKSSQSLLNSRTRKNYLAWYQQKPGQSPKLLIYWASTRESGVPDRFTGSGSGTDFTLTISSVQAEDLAVYYCKQSYNLPTFGSGTKIEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC";

	public static final String BSA = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA";

	/**
	 * C is +57 by default Currently C is +57 TO DO: We should allow users to
	 * specify fixed modifications like C+57 instead of using it by default.
	 */
	public static double[] AAMasses = { 71.03711, Double.NEGATIVE_INFINITY,
			103.08406, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891,
			113.08406, Double.NEGATIVE_INFINITY, 128.09496, 113.08406,
			131.04049, 114.04293, Double.NEGATIVE_INFINITY, 97.05276,
			128.05858, 156.10111, 87.03203, 101.04768,
			Double.NEGATIVE_INFINITY, 99.06841, 186.07931,
			Double.NEGATIVE_INFINITY, 163.06333 };

	public static int[] AAMassesScaled = { 71037, Integer.MIN_VALUE, 103084,
			115026, 129042, 147068, 57021, 137058, 113084, Integer.MIN_VALUE,
			128094, 113084, 131040, 114042, Integer.MIN_VALUE, 97052, 128058,
			156101, 87032, 101047, Integer.MIN_VALUE, 99068, 186079,
			Integer.MIN_VALUE, 163063 };

	/*
	 * 
	 * A 71037 C 160084 D 115026 E 129042 F 147068 G 57021 H 137058 I 113084 K
	 * 128094 L 113084 M 131040 N 114042 P 97052 Q 128058 R 156101 S 87032 T
	 * 101047 U 150800 V 99068 W 186079 Y 163063
	 */

	/**
	 * In order to use integer operations instead of float operations, most
	 * things are scaled by MASS_SCALE (allowing up to 3 decimal precision when
	 * MASS_SCALE = 1000. It's important to specify when masses are expected to
	 * be scaled and when the output is scaled
	 */
	public static final int MASS_SCALE = 1000;

	/**
	 * Different mass spectrometers can achieve different mass accuracies. These
	 * Tolerances represent the mass tolerance (in Daltons, and not scaled by
	 * MASS_SCALE). TO DO: Currently the flow uses the LTQFTIonTolerance by
	 * default, we should change this to be an option associated with either a
	 * PRMSpectrum object or an MSType object.
	 */
	// public static final double OrbiIonTolerance = 0.05;
	// public static final double LTQFTIonTolerance = 0.5;
	// public static final double LTQIonTolerance = 0.3;
	public static final double DefaultFragTolerance = 0.5;
	public static final double DefaultPMTolerance = 2.0;



	/**
	 * Constants for template selection
	 */
	public static final int MAX_SHARED_PEPTIDES = 2;
	public static final char CLASS_DELIM = '.';
	// public static final String CLASS_DELIM_REGEX = "\\.";

	/**
	 * Constants for Initial Inspect Searches, etc.
	 */
	// The minimum number of spectra a sequence in the search database must have
	// to be considered the
	// winner of the majority vote
	public static final int MIN_SUPPORTING_SPECTRA = 50;
	public static final double MIN_TEMPLATE_LIKELIHOOD = 1.0;
	public static final double MIN_TEMPLATE_LIKELIHOOD_GENOMIC = 3.00;

	// The min number of peaks a spectrum must have to be considered.
	public static final int MIN_PRM_PEAKS = 10;

	/**
	 * Constants for GatherConsensus
	 */
	// Initial match state score
	public static final double INIT_MATCH_STATE_SCORE = 15;

	// Initial transition probabilities
	public static final double INIT_MATCH_MATCH = 0.8;
	public static final double INIT_MATCH_INSERT = 0.15;
	public static final double INIT_MATCH_DELETE = 0.05;

	public static final double INIT_INSERT_MATCH = 0.80;
	public static final double INIT_INSERT_INSERT = 0.15;
	public static final double INIT_INSERT_DELETE = 0.05;

	public static final double INIT_DELETE_MATCH = 0.85;
	public static final double INIT_DELETE_INSERT = 0.05;
	public static final double INIT_DELETE_DELETE = 0.1;

	// Minimum number of peaks that must be aligned for a spectrum to be
	// included
	// in the model
	public static int MIN_PEAKS_ALIGNED = 4;

	// Pseudocounts when updating transition probabilities
	public static final int PSEUDOCOUNT_MATCH = 7; // Must be non-zero for
	// ProducePRMConsensus
	public static final int PSEUDOCOUNT_DELETE = 1; // Must be non-zero for
	// ProducePRMConsensus
	public static final int PSEUDOCOUNT_INSERT = 1;

	// The fraction of the mean match score that a new state must meet
	public static final double HMM_MIN_SCORE_FRAC = 3.0;
	//public static double MISSING_PEAK_PENALTY = -1.0;
	//public static boolean APPLY_MATCH_PENALTY = false;

	// Learning rate for the transitions of the model
	public static final double TRANSITION_LEARNING_RATE = 0.7;

	public static final double HMM_SCORE_SCALE = 1.2;

	// Min number of peaks that are close together in a insert state to be
	// considered as
	// a candidate for a match state
	public static int MIN_PEAK_GROUP_SIZE = 2;
	public static double MIN_PEAK_GROUP_SCORE = 15;

	// Random model parameters, determined by empirical evidence from
	// DetermePA.py
	// public static final double pA = 0.07; //probability of observing a real
	// peak
	// public static final double pNA = 0.04; //probability of observing a fake
	// peak
	public static final double[] pA = { 0.31096, 0.24751, 0.007 };
	public static final double NotpA = (1 - 0.31096 - 0.24751 - 0.007);
	public static final double[] pNA = { 0.03878, 0.0019, 0.00003 };
	public static final double NotpNA = (1 - 0.03878 - 0.0019 - 0.00003);
	public static final double[] ScoreRanges = { 12.856, 25.712,
			Double.MAX_VALUE };
	public static final double MAX_PRM_SCORE = 38.569;

	/*
	 * Constants for SpecAlign
	 */
	// Max number of spectra to align
	// public static final double MAX_PRM_SCORE = 15;
	public static int MAX_ALIGNED_SPECTRA = 5;
	public static int MIN_ALIGNED_SPECTRA = MIN_PEAK_GROUP_SIZE;
	public static final int NUM_DECOY_ALIGNMENTS = 1000;
	public static final double MSALIGN_PVALUE_CUTOFF = 1.00;
	public static double MSALIGN_SCORE_CUTOFF = 70.00;
	public static final int DECOY_EXPLOSION = 1000;
	public static int MIN_OVERLAPPING_PEAKS = MIN_PEAKS_ALIGNED;
	public static final int MSALIGN_MIN_SIMILARITY = MIN_PEAKS_ALIGNED;

	// The min number of exactly matching AA's for two adjacent seeds to be
	// merged
	public static final int MIN_OVERLAPPING_AA = 3;
	public static final int MIN_OVERLAPPING_AA_INTERNAL = 6;
	public static final int MAX_SKIPPED_PRMs = 2;

	// Amount of a seed sequence used for MSAlignment
	public static final int MAX_ALIGNMENT_SEQ = 10;

	/*
	 * Constants for path reconstruction
	 */
	public static final int MIN_ALIGNMENT_LENGTH = 4;
	public static final int MAX_BEST_PATHS = 3;

	/*
	 * Constants for Template TEST
	 */
	// The number of seeds to randomly generate for the Template Test
	public static final int NUM_TEST_SEEDS = AntibodyUtils.aBTLA_OracleSequence
			.length() - 3 * MAX_ALIGNMENT_SEQ;

	// Histogram parameters
	public static final double SCORE_BIN_SIZE = 1.0;
	public static final double PVALUE_BIN_SIZE = 0.01;

	public static final double DEFAULT_PRM_SCORE = 10.0;
	public static final String CLASS_DELIM_REGEX = "[.]";

	// Sequence alignment parameters
	public static final int MATCH_SCORE = 1;
	public static final int GAP_SCORE = 0;
	public static final int MISMATCH_SCORE = 0;

	public static double GetModMassFromName(String ModStr) {
		if (ModStr.charAt(0) == '+' || ModStr.charAt(0) == '-')
			return Double.parseDouble(ModStr);
		else if (ModStr.indexOf("->") >= 0) {
			double OldAAMass = AntibodyUtils.AAMasses[ModStr.charAt(0) - 65];
			double NewAAMass = AntibodyUtils.AAMasses[ModStr.charAt(ModStr
					.length() - 1) - 65];
			return NewAAMass - OldAAMass;
		} else {
			System.err.println("ERROR: Unable to determine delta mass of '"
					+ ModStr + "'");
			Utils.WaitForEnter();
			return 0;
		}
	}

	public static ModSite[] GetModsFromPeptide(String OrigPeptide) {
		String Peptide = OrigPeptide;
		String UnModded = Utils.GetUnModded(OrigPeptide);
		if (OrigPeptide.indexOf(".") >= 0) {
			// System.out.println("Peptide: " + Peptide);
			// String[] Parts = Peptide.split(".");
			// System.out.println("Parts: " + Parts.length);
			Peptide = OrigPeptide.substring(2, OrigPeptide.length() - 2);
		}
		System.out.println("Parsing mods from : " + OrigPeptide);
		ArrayList TempMods = new ArrayList();
		for (int i = 0; i < Peptide.length(); ++i) {
			if (Peptide.charAt(i) == '+'
					|| (Peptide.charAt(i) == '-' && Peptide.charAt(i + 1) != '>')) {
				int j = i + 1;
				String ModStr = "" + Peptide.charAt(i);
				System.out.println("Found mod at pos " + i);
				while (Character.isDigit(Peptide.charAt(j))
						|| Peptide.charAt(j) == '.') {
					ModStr += Peptide.charAt(j);
					j += 1;
				}
				System.out.println("modStr: " + ModStr);
				ModSite M = new ModSite(ModStr, i - 1,
						Double.parseDouble(ModStr), UnModded);
				TempMods.add(M);
				i = j - 1;
			} else if (Peptide.charAt(i) == '>') {
				String ModStr = Peptide.substring(i - 2, i + 2).toUpperCase();
				System.out.println("Found mod at pos " + i);
				System.out.println("modStr: " + ModStr);
				double OldAAMass = AntibodyUtils.AAMasses[ModStr.charAt(0) - 65];
				double NewAAMass = AntibodyUtils.AAMasses[ModStr.charAt(ModStr
						.length() - 1) - 65];
				ModSite M = new ModSite(ModStr, i - 3, NewAAMass - OldAAMass,
						UnModded);
				TempMods.add(M);
				i = i + 1;
			}
		}
		if (TempMods.size() == 0)
			return null;

		ModSite[] Ret = new ModSite[TempMods.size()];
		for (int j = 0; j < Ret.length; ++j)
			Ret[j] = (ModSite) (TempMods.get(j));

		return Ret;
	}

	public static int Round(double Val) {
		int IntPart = (int) (Val);

		if (Val >= IntPart + 0.5)
			return IntPart + 1;
		return IntPart;
	}

	public static boolean IdenticalPRMs(String A, String B) {
		int[] PRMsA = AntibodyUtils.GetSeqPRMScaled(A);
		int[] PRMsB = AntibodyUtils.GetSeqPRMScaled(B);

		if (PRMsA.length != PRMsB.length)
			return false;

		for (int i = 0; i < PRMsA.length; ++i) {
			if (PRMsA[i] != PRMsB[i])
				return false;
		}
		return true;

	}

	public static int DeterminePRMScoreZone(double CurrScore) {
		for (int i = 0; i < AntibodyUtils.ScoreRanges.length; ++i) {
			if (CurrScore < AntibodyUtils.ScoreRanges[i])
				return i;

		}
		return -1;// The ranges should never miss a score!!
	}

	/**
	 * Performs a deep copy of an Arraylist
	 * 
	 * @param oldList
	 *            The ArrayList of Integers to copy
	 * @return A new ArrayList, with deep copies of the Integer objects in
	 *         oldList.
	 */
	public static ArrayList CopyIntArrayList(ArrayList oldList) {
		ArrayList newList = new ArrayList();
		for (int i = 0; i < oldList.size(); i++) {
			int value = ((Integer) (oldList.get(i))).intValue();
			newList.add(new Integer(value));
		}
		return newList;
	}

	public static int CountChar(String S, char c) {
		int count = 0;
		for (int i = 0; i < S.length(); ++i) {
			if (S.charAt(i) == c)
				count++;
		}
		return count;
	}

	public static boolean IntArrayListContainedIn(ArrayList ListofLists,
			ArrayList NewList) {

		for (int i = 0; i < ListofLists.size(); ++i) {
			ArrayList CurrList = (ArrayList) (ListofLists.get(i));

			if (CurrList.size() == NewList.size()) {
				boolean equal = true;
				for (int j = 0; j < CurrList.size(); ++j) {
					int A = ((Integer) (CurrList.get(j))).intValue();
					int B = ((Integer) (NewList.get(j))).intValue();
					if (A != B) {
						equal = false;
						break;
					}
				}
				if (equal) {
					return true;
				}
			}
		}
		return false;

	}

	/**
	 * Computes teh log base 10 of the values of a list. Uses Utils.Log10
	 * 
	 * @param Values
	 *            list of values to take the log of
	 * @return Returns the original list, but with the log value of each element
	 */
	public static double[] Log10Array(double[] Values) {
		double[] logArray = new double[Values.length];
		for (int i = 0; i < Values.length; ++i) {
			logArray[i] = Log10(Values[i]);
		}
		return logArray;
	}

	/**
	 * Computes the log base 10 of the value, or returns
	 * Double.NEGATIVE_INFINITY if Value = 0
	 * 
	 * @param Value
	 *            The number of compute the log of
	 * @return Returns NEGATIVE_INFINITY if Value = 0.0, otherwise returns the
	 *         log of the Value.
	 */
	public static double Log10(double Value) {
		if (Value == 0.0)
			return Double.NEGATIVE_INFINITY;
		else
			return Math.log10(Value);
	}

	public static double GetVariance(double Mean, double Tolerance) {
		// return 0.5/2*0.5/2;
		return Math.pow((int) ((MASS_SCALE * Tolerance) / 2), 2);
		// return 300*300;
	}

	/**
	 * Converts an ArrayList of Integer objects to an array of ints
	 * 
	 * @param List
	 *            the ArrayList of Integers
	 * @return a standard array of the ints
	 */
	public static int[] ConvertIntegerArrayList(ArrayList List) {
		int[] newList = new int[List.size()];
		for (int i = 0; i < List.size(); ++i) {

			Integer Temp = (Integer) List.get(i);
			newList[i] = Temp.intValue();
		}
		return newList;
	}

	public static double[] ConvertDoubleArrayList(ArrayList List) {
		double[] newList = new double[List.size()];
		for (int i = 0; i < List.size(); ++i) {

			Double Temp = (Double) (List.get(i));
			newList[i] = Temp.doubleValue();
		}
		return newList;
	}

	/**
	 * Converts a list of doubles to a scaled list of ints by multiplying by
	 * MASS_SCALE
	 * 
	 * @param oldList
	 *            The list of doubles to convert
	 * @return Returns the scaled array of ints
	 */
	public static int[] ScaleArray(double[] oldList) {
		int[] newList = new int[oldList.length];
		for (int i = 0; i < oldList.length; ++i) {
			int ScaledMass = (int) (oldList[i] * MASS_SCALE);
			newList[i] = ScaledMass;
		}
		return newList;
	}

	/**
	 * Converts a list of ints to an unscaled list of doubles by dividing by
	 * MASS_SCALE Note that the precision is limited by MASS_SCALE
	 * 
	 * @param oldList
	 *            The list of ints to convert
	 * @return Returns a list of doubles
	 */
	public static double[] UnScaleArray(int[] oldList) {
		double[] newList = new double[oldList.length];
		for (int i = 0; i < oldList.length; ++i) {
			double UnScaledMass = (double) (oldList[i]) / MASS_SCALE;
			newList[i] = UnScaledMass;
		}
		return newList;
	}

	/**
	 * Returns the string version of the mass of two amino acids, which are
	 * within Utils.IonTolerance of Mass. The first pair of amino acids
	 * (alphabetically) which falls within IonTolerance of Mass is returned, NOT
	 * necessarily the closest pair.
	 * 
	 * @param Mass
	 *            The mass that we are trying to explain with 2 amino acids. The
	 *            Mass is scaled.
	 * @return The original mass (Scaled). Returns "0" if no pair is found.
	 */
	public static String FindClosestDoubleAAScaled(int Mass, double Tolerance) {
		for (int i = 0; i < AAMasses.length; ++i) {
			if ((int) (AAMassesScaled[i]) == Integer.MIN_VALUE)
				continue;
			for (int j = 0; j < AAMasses.length; ++j) {
				if ((int) (AAMassesScaled[j]) == Integer.MIN_VALUE)
					continue;
				int ScaledMass = AAMassesScaled[i] + (AAMassesScaled[j]);
				if (Math.abs(ScaledMass - Mass) <= (int) (Tolerance * MASS_SCALE))
					return Mass + "";
			}
		}
		return "0";
	}

	/**
	 * Returns the two amino acids, which are within Utils.IonTolerance of and
	 * closest to Mass.
	 * 
	 * @param Mass
	 *            The mass that we are trying to explain with 2 amino acids. The
	 *            Mass is scaled.
	 * @return The original mass (Scaled). Returns "0" if no pair is found.
	 */
	public static String FindClosestDoubleAAScaledSeqs(int Mass,
			double Tolerance) {
		int distance = Integer.MAX_VALUE;
		String currAAs = null;

		for (int i = 0; i < AAMasses.length; ++i) {
			if ((int) (AAMassesScaled[i]) == Integer.MIN_VALUE)
				continue;
			for (int j = 0; j < AAMasses.length; ++j) {
				if ((int) (AAMassesScaled[j]) == Integer.MIN_VALUE)
					continue;
				int ScaledMass = AAMassesScaled[i] + (AAMassesScaled[j]);
				int diff = Math.abs(ScaledMass - Mass);
				if (diff <= (int) (Tolerance * MASS_SCALE) && diff < distance) {
					distance = diff;
					currAAs = (char) (i + 65) + "" + (char) (j + 65);
				}
			}
		}
		return currAAs;
	}

	/**
	 * Returns the amino acid which falls within Utils.IonTolerance of Mass. If
	 * multiple amino acids match the Mass, only the first one alphabetically is
	 * reported.
	 * 
	 * @param Mass
	 *            The scaled mass in question
	 * @param Tolerance
	 *            The unscaled tolerance of the mass in Da
	 * @return The amino acid which is close to Mass. Returns '0' if no amino
	 *         acid is found.
	 */
	public static char FindClosestAAScaled(int Mass, double Tolerance) {
		int SmallestDelta = Integer.MAX_VALUE;
		char Ret = '0';
		for (int i = 0; i < AAMasses.length; ++i) {
			int ScaledMass = (int) (AAMassesScaled[i]);
			if (AAMassesScaled[i] == Integer.MIN_VALUE)
				continue;
			if (Math.abs(ScaledMass - Mass) <= (int) (Tolerance * MASS_SCALE)) {
				if (Math.abs(ScaledMass - Mass) < SmallestDelta) {
					Ret = (char) (i + 65);
					SmallestDelta = Math.abs(ScaledMass - Mass);
				}
			}
		}
		return Ret;
	}

	/**
	 * Converts the PepNovo reported M+H to the sum of the amino acid masses.
	 * Subtracts the hydrogen, the water lost, and scaled by MASS_SCALE.
	 * 
	 * @param Value
	 *            The M+H reported by PepNovo(unscaled)
	 * @return The scaled mass of the amino acid sequence.
	 */
	public static int ComputeMassFromMHScaled(double Value) {
		return (int) ((Value - Utils.HydrogenMass - Utils.WaterMass) * MASS_SCALE);
	}

	/**
	 * Determines the mass of the given sequence, not scaled by MASS_SCALE
	 * 
	 * @param Sequence
	 *            The sequence to be converted
	 * @return The sum of the masses of the amino acids, not scaled.
	 */
	public static int GetSeqMassScaled(String Sequence) {

		double TotalMass = getMassFromSeq(Sequence);
		return (int) (TotalMass*MASS_SCALE);
	}

	/**
	 * Produces the prefix residue masses for a given sequence, including 0 as
	 * the first PRM.
	 * 
	 * @param Sequence
	 *            the sequence to calculate the PRM
	 * @return a list of the PRM masses, scaled by MASS_SCALE
	 */
	public static int[] GetSeqPRMScaled(String Sequence) {

		if (Sequence == null)
			ErrorThrower.ThrowError(20,
					"AntibodyUtils.GetSeqPRMScaled.Sequence");
		// int [] PRMs = new int[Sequence.length()+1];
		ArrayList PRMs = new ArrayList();
		PRMs.add(new Integer(0));
		// PRMs[0] = 0;
		int CurrSum = 0;

		for (int i = 0; i < Sequence.length(); ++i) {
			if (Sequence.charAt(i) == '[') {
				int End = Sequence.indexOf(']', i);
				if (End < 0) {
					System.err
							.println("ERROR: Illformed sequence with numerical mass '"
									+ Sequence + "'!");
					return null;
				}
				int PRM = (int) (Double.parseDouble(Sequence.substring(i + 1,
						End)) * AntibodyUtils.MASS_SCALE);
				CurrSum += PRM;
				PRMs.add(new Integer(CurrSum));
				i = End;
			} else {
				int Index = Sequence.charAt(i) - 65;
				if (Index < 0 || Index > AAMasses.length)
					System.out.println("**ERROR:Inappropriate Char: "
							+ Sequence.charAt(i) + " in **" + Sequence + "**");

				CurrSum += (int) (AAMassesScaled[Sequence.charAt(i) - 65]);
				// System.out.println(" " + Sequence.charAt(i) + "-" +
				// AAMassesScaled[Sequence.charAt(i)-65]);
				PRMs.add(new Integer(CurrSum));
			}
		}

		return AntibodyUtils.ConvertIntegerArrayList(PRMs);
	}

	public static boolean EquivalentScaled(int A, int B, double Tolerance) {
		if (Math.abs(A - B) <= (int) (Tolerance * AntibodyUtils.MASS_SCALE))
			return true;
		return false;
	}

	public static int[] MergePRMs(int[] PRMA, int[] PRMB, double Tolerance) {
		int[][] Alignment = new int[PRMA.length][PRMB.length];
		int[][][] PrevSteps = new int[PRMA.length][PRMB.length][2];

		for (int i = 0; i < PRMA.length; ++i) {
			Alignment[i][0] = 1;
			PrevSteps[i][0][0] = -1;
			PrevSteps[i][0][1] = -1;
		}

		int TotalMaxScore = 0;
		int[] BestCells = new int[2];
		for (int i = 1; i < PRMA.length; ++i) {
			for (int j = 1; j < PRMB.length; ++j) {
				int MaxScore = 0;
				PrevSteps[i][j][0] = -1;
				PrevSteps[i][j][1] = -1;

				// Consider matching PRMA[i] with PRMB[j], look at all previous
				// matchings
				for (int Previ = 0; Previ < i; Previ++) {
					for (int Prevj = 0; Prevj < j; Prevj++) {
						int Diff1 = PRMA[i] - PRMA[Previ];
						int Diff2 = PRMB[j] - PRMB[Prevj];

						if (AntibodyUtils.EquivalentScaled(Diff1, Diff2,
								Tolerance)) {
							if (Alignment[Previ][Prevj] + 1 > MaxScore) {
								MaxScore = Alignment[Previ][Prevj] + 1;
								PrevSteps[i][j][0] = Previ;
								PrevSteps[i][j][1] = Prevj;
							}
						}
					}
				}
				Alignment[i][j] = MaxScore;
				if ((j == PRMB.length - 1 || i == PRMA.length - 1)
						&& MaxScore > TotalMaxScore) {
					TotalMaxScore = MaxScore;
					BestCells[0] = i;
					BestCells[1] = j;
				}
			}
		}

		if (TotalMaxScore <= AntibodyUtils.MIN_OVERLAPPING_AA)
			return null;

		// Sanity check
		if (BestCells[0] != PRMA.length - 1 && BestCells[1] != PRMA.length)
			return null;

		int MassShift = PRMA[BestCells[0]] - PRMB[BestCells[1]];

		// Sanity Check, B should start somewhere after A
		if (MassShift < 0
				&& Math.abs(MassShift) > Tolerance * AntibodyUtils.MASS_SCALE)
			return null;

		ArrayList TempPRMs = new ArrayList();
		int AIndex = 0;
		int BIndex = 0;

		while (AIndex < PRMA.length && BIndex < PRMB.length) {
			if (AntibodyUtils.EquivalentScaled(PRMA[AIndex], PRMB[BIndex]
					+ MassShift, Tolerance)) {
				TempPRMs.add(new Integer(PRMA[AIndex]));
				AIndex++;
				BIndex++;
			} else if (PRMA[AIndex] < PRMB[BIndex] + MassShift) {
				TempPRMs.add(new Integer(PRMA[AIndex]));
				AIndex++;
			} else {
				TempPRMs.add(new Integer(PRMB[BIndex] + MassShift));
				BIndex++;
			}
		}
		while (AIndex < PRMA.length) {
			TempPRMs.add(new Integer(PRMA[AIndex]));
			AIndex++;
		}
		while (BIndex < PRMB.length) {
			TempPRMs.add(new Integer(PRMB[BIndex] + MassShift));
			BIndex++;
		}
		return AntibodyUtils.ConvertIntegerArrayList(TempPRMs);

	}

	/**
	 * Sometimes we may compare two sequences that are wrong at the edges. We
	 * will allow merging where the strongest alignment path does not end at the
	 * end of one of the sequences.
	 * 
	 * @param PRMA
	 * @param ScoresA
	 * @param PRMB
	 * @param ScoresB
	 * @return
	 */
	public static ConsensusPRMSpectrum MergePRMsSingleSkipB(int[] PRMA,
			double[] ScoresA, int[] PRMB, double[] ScoresB, double Tolerance) {
		int[][] Alignment = new int[PRMA.length][PRMB.length];
		int[][][] PrevSteps = new int[PRMA.length][PRMB.length][2];

		for (int i = 0; i < PRMA.length; ++i)
			for (int j = 0; j < PRMB.length; ++j)
				Alignment[i][j] = Integer.MIN_VALUE;

		for (int i = 0; i < PRMA.length; ++i) {
			Alignment[i][0] = 1;
			PrevSteps[i][0][0] = -1;
			PrevSteps[i][0][1] = -1;
			Alignment[i][1] = 1; // Can allow skip of first AA
			PrevSteps[i][1][0] = -1;
			PrevSteps[i][1][1] = -1;
		}

		int TotalMaxScore = 0;
		int[] BestCells = new int[2];
		for (int i = 1; i < PRMA.length; ++i) {
			for (int j = 1; j < PRMB.length; ++j) {
				int MaxScore = 0;
				PrevSteps[i][j][0] = -1;
				PrevSteps[i][j][1] = -1;

				// Consider matching PRMA[i] with PRMB[j], look at previous
				// matchings, allowing 1 skipped PRM
				for (int Previ = Math.max(0, i - 2); Previ < i; Previ++) {
					for (int Prevj = Math.max(0, j - 2); Prevj < j; Prevj++) {
						int Diff1 = PRMA[i] - PRMA[Previ];
						int Diff2 = PRMB[j] - PRMB[Prevj];

						if (AntibodyUtils.EquivalentScaled(Diff1, Diff2,
								Tolerance)) {
							if ((Alignment[Previ][Prevj] > Integer.MIN_VALUE)
									&& Alignment[Previ][Prevj] + 1 > MaxScore) {
								MaxScore = Alignment[Previ][Prevj] + 1;
								PrevSteps[i][j][0] = Previ;
								PrevSteps[i][j][1] = Prevj;
							}
						}
					}
				}
				Alignment[i][j] = MaxScore;
				if (MaxScore > TotalMaxScore) {
					TotalMaxScore = MaxScore;
					BestCells[0] = i;
					BestCells[1] = j;
				}
			}
		}

		// If we end at the end of one of our Strings, then we only need
		// MIN_OVERLAPPING_AA
		if ((BestCells[0] == PRMA.length - 1 || BestCells[1] == PRMA.length)
				&& TotalMaxScore < AntibodyUtils.MIN_OVERLAPPING_AA) {
			System.out.println("MergePRMs: TotalMaxScore " + TotalMaxScore
					+ " < " + AntibodyUtils.MIN_OVERLAPPING_AA);
			return null;
		}
		// If we end in the middle of one of our Strings, then we only need
		// MIN_OVERLAPPING_AA_INTERNAL
		if ((BestCells[0] != PRMA.length - 1 && BestCells[1] != PRMA.length)
				&& TotalMaxScore < AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL) {
			System.out.println("MergePRMs: TotalMaxScore " + TotalMaxScore
					+ " < " + AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL);
			return null;
		}

		int MassShift = PRMA[BestCells[0]] - PRMB[BestCells[1]];

		// Sanity Check, B should start somewhere after A
		// if(MassShift < 0 && Math.abs(MassShift) >
		// Utils.LTQFTIonTolerance*Utils.MASS_SCALE)
		// return null;

		System.out.println("Seed should be shifted by " + MassShift);
		int[] CurrCells = { BestCells[0], BestCells[1] };
		while (CurrCells[0] != -1) {
			System.out.println("PRMA[" + CurrCells[0] + "]="
					+ PRMA[CurrCells[0]]);
			System.out.println("PRMB[" + CurrCells[1] + "]="
					+ (PRMB[CurrCells[1]] + MassShift) + "("
					+ PRMB[CurrCells[1]] + ")");
			int i = CurrCells[0];
			int j = CurrCells[1];
			CurrCells[0] = PrevSteps[i][j][0];
			CurrCells[1] = PrevSteps[i][j][1];

		}
		ArrayList TempPRMs = new ArrayList();
		ArrayList PRMScores = new ArrayList();
		int AIndex = 0;
		int BIndex = 0;

		while (AIndex < PRMA.length && BIndex < PRMB.length) {
			if (AntibodyUtils.EquivalentScaled(PRMA[AIndex], PRMB[BIndex]
					+ MassShift, Tolerance)) {
				TempPRMs.add(new Integer(PRMA[AIndex]));
				PRMScores.add(new Double(
						(ScoresA[AIndex] + ScoresB[BIndex]) / 2));
				AIndex++;
				BIndex++;
			} else if (PRMA[AIndex] < PRMB[BIndex] + MassShift) {
				TempPRMs.add(new Integer(PRMA[AIndex]));
				PRMScores.add(new Double(ScoresA[AIndex]));
				AIndex++;
			} else {
				TempPRMs.add(new Integer(PRMB[BIndex] + MassShift));
				PRMScores.add(new Double(ScoresB[BIndex]));
				BIndex++;
			}
		}
		while (AIndex < PRMA.length) {
			TempPRMs.add(new Integer(PRMA[AIndex]));
			PRMScores.add(new Double(ScoresA[AIndex]));
			AIndex++;
		}
		while (BIndex < PRMB.length) {
			TempPRMs.add(new Integer(PRMB[BIndex] + MassShift));
			PRMScores.add(new Double(ScoresB[BIndex]));
			BIndex++;
		}
		ConsensusPRMSpectrum Ret = new ConsensusPRMSpectrum();
		Ret.ScaledPeaks = AntibodyUtils.ConvertIntegerArrayList(TempPRMs);
		Ret.PeakScores = AntibodyUtils.ConvertDoubleArrayList(PRMScores);
		return Ret;

	}

	public static String[] GetIntersection(String[] A, String[] B) {
		ArrayList Temp = new ArrayList();
		for (int i = 0; i < A.length; ++i) {
			if (A[i] == null)
				continue;
			for (int j = 0; j < B.length; ++j) {
				if (B[j] == null)
					continue;
				if (A[i].compareTo(B[j]) == 0) {
					Temp.add(A[i]);
					break;
				}
			}
		}
		return AntibodyUtils.ConvertStringArrayList(Temp);
	}

	public static String[] GetUnion(String[] TempA, String[] TempB) {
		Hashtable Sps = new Hashtable();

		for (int k = 0; k < TempA.length; ++k) {
			if (TempA[k] != null)
				Sps.put(TempA[k], new Integer(1));
		}

		for (int k = 0; k < TempB.length; ++k) {
			if (TempB[k] != null)
				Sps.put(TempB[k], new Integer(1));
		}
		String[] ret = AntibodyUtils.HashKeysToStringArray(Sps);
		for (int k = 0; k < ret.length; ++k) {
			if (ret[k] == null)
				return null;
		}
		return ret;
	}

	public static String[] HashKeysToStringArray(Hashtable H) {
		Enumeration Keys = H.keys();
		String[] Ret = new String[H.size()];
		int i = 0;
		while (Keys.hasMoreElements()) {
			Ret[i] = (String) (Keys.nextElement());
			i++;
		}
		return Ret;
	}

	/**
	 * We assume that the alignment contains the end of SequenceA, and the
	 * beginning of SequenceB or the beginning of SequenceA A: ----------- B:
	 * -----------
	 * 
	 * A: ---------- B: -----------------
	 * 
	 * @param PRMA
	 * @param ScoresA
	 * @param PRMB
	 * @param ScoresB
	 * @return
	 */

	public static int[] GetBestAlignmentRecursive(String SeqA, String SeqB,
			double Tolerance) {
		int[] AnchorPRMs = AntibodyUtils.GetSeqPRMScaled(SeqA);
		int[] OraclePRMs = AntibodyUtils.GetSeqPRMScaled(SeqB);
		int[][] Alignment = new int[AnchorPRMs.length][OraclePRMs.length];
		int[][][] PrevSteps = new int[AnchorPRMs.length][OraclePRMs.length][2];

		// System.out.println("SeqA length: " + AnchorPRMs.length);
		// System.out.println("SeqB length: " + OraclePRMs.length);
		boolean LocalDebug = false;

		PrevSteps[0][0][0] = -1;
		PrevSteps[0][0][1] = -1;
		if (LocalDebug)
			System.out.println("Aligning: " + SeqA + " to " + SeqB);
		int MaxScore = 0;
		int[] BestCells = { 0, 0 };
		for (int i = 1; i < AnchorPRMs.length; ++i) {
			for (int j = 1; j < OraclePRMs.length; ++j) {
				PrevSteps[i][j][0] = -1;
				PrevSteps[i][j][1] = -1;

				for (int Previ = i - 1; Previ >= Math.max(0, i - 3); --Previ) {
					for (int Prevj = j - 1; Prevj >= Math.max(0, j - 3); --Prevj) {
						int Diff1 = AnchorPRMs[Previ] - AnchorPRMs[i];
						int Diff2 = OraclePRMs[Prevj] - OraclePRMs[j];
						if (AntibodyUtils.EquivalentScaled(Diff1, Diff2,
								Tolerance)) {
							if (Alignment[Previ][Prevj] + 1 > Alignment[i][j]) {
								Alignment[i][j] = Alignment[Previ][Prevj] + 1;
								PrevSteps[i][j][0] = Previ;
								PrevSteps[i][j][1] = Prevj;
							}
						}
					}
				}

				if (Alignment[i][j] > MaxScore) {
					MaxScore = Alignment[i][j];
					BestCells[0] = i;
					BestCells[1] = j;

				}
			}
		}

		if (MaxScore < AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL)
			return null;

		int SkippedA = 0;
		int SkippedB = 0;
		int AStart = 0;
		int AEnd = 0;
		int BStart = 0;
		int BEnd = 0;

		if (LocalDebug)
			System.out.println("Best end is at position " + BestCells[1]
					+ " with score " + MaxScore);

		int Curri = BestCells[0];
		int Currj = BestCells[1];

		int Previ = PrevSteps[Curri][Currj][0];
		int Prevj = PrevSteps[Curri][Currj][1];

		// SkippedA += (AnchorPRMs.length - Curri - 1);
		// SkippedB += (OraclePRMs.length - Currj - 1);

		AEnd = BestCells[0];
		BEnd = BestCells[1];

		while (Previ != -1 && Prevj != -1) {

			SkippedA += Math.max(0, (Curri - Previ - 1));
			SkippedB += Math.max(0, (Currj - Prevj - 1));

			Curri = Previ;
			Currj = Prevj;
			if (Previ == 0 || Prevj == 0)
				break;

			Previ = PrevSteps[Curri][Currj][0];
			Prevj = PrevSteps[Curri][Currj][1];
		}
		AStart = Curri;
		BStart = Currj;

		if (LocalDebug) {
			System.out.println("Alignment of " + SeqA + " to " + SeqB + ":"
					+ MaxScore);
			System.out.println("A: " + AStart + "-" + AEnd + " skipped "
					+ SkippedA);
			System.out.println("B: " + BStart + "-" + BEnd + " skipped "
					+ SkippedB);
		}
		String[] SeqA_AAs = AntibodyUtils.SplitStringtoAA(SeqA);
		int[] PrefRet = null;
		int[] SuffixRet = null;
		// See if we can get a prefix alignment!!!)
		if (AStart >= AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL
				&& BStart >= AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL) {

			String PrefixString = "";
			for (int i = 0; i < AStart; ++i)
				PrefixString += SeqA_AAs[i];
			if (LocalDebug)
				System.out.println("Finding prefix alignment of "
						+ PrefixString + " and " + SeqB.substring(0, BStart));
			PrefRet = AntibodyUtils.GetBestAlignmentRecursive(PrefixString,
					SeqB.substring(0, BStart), Tolerance);
		}
		if (SeqA_AAs.length - AEnd >= AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL
				&& SeqB.length() - BEnd >= AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL) {
			String SuffixString = "";
			for (int i = AEnd; i < SeqA_AAs.length; ++i)
				SuffixString += SeqA_AAs[i];
			if (LocalDebug)
				System.out.println("Finding suffix alignment of "
						+ SuffixString + " and " + SeqB.substring(BEnd));
			SuffixRet = AntibodyUtils.GetBestAlignmentRecursive(SuffixString,
					SeqB.substring(BEnd), Tolerance);
		}
		if (LocalDebug) {
			System.out.println("FInished alignment of " + SeqA + " to " + SeqB);
			System.out.println("Main alignment: " + MaxScore + ", " + AStart
					+ " to " + AEnd + " skipping " + SkippedA + ", and "
					+ BStart + " to " + BEnd + " skipping " + SkippedB);
			System.out.println("SeqAAs: " + SeqA_AAs.length);
		}
		int[] Ret = { MaxScore, AStart, AEnd, BStart, BEnd, SkippedA, SkippedB };
		if (PrefRet != null) {
			if (LocalDebug) {
				System.out.println("Prefix alignment: " + PrefRet[0] + ", "
						+ PrefRet[1] + " to " + PrefRet[2] + " skipping "
						+ PrefRet[5] + ", and " + PrefRet[3] + " to "
						+ PrefRet[4] + " skipping " + PrefRet[6]);
			}
			Ret[0] += PrefRet[0]; // Update score of alignment
			Ret[1] = PrefRet[1]; // AStart is start of prefix alignment
			Ret[3] = PrefRet[3]; // BStart is start of prefix alignment

			Ret[5] += PrefRet[5] + (AStart - PrefRet[2]); // Skipped = skipped
			// in this alignment
			// + skipped in
			// prefix alignment
			// + preceding
			// prefix alignment
			// + succeeding
			// prefix alignment
			Ret[6] += PrefRet[6] + (BStart - PrefRet[4]);
			if (LocalDebug) {
				System.out.println("Summed alignment: " + Ret[0] + ", "
						+ Ret[1] + " to " + Ret[2] + " skipping " + Ret[5]
						+ ", and " + Ret[3] + " to " + Ret[4] + " skipping "
						+ Ret[6]);
			}
		} else {
			if (LocalDebug)
				System.out.println("No prefix alignment");
			// Ret[5] += AStart;
			// Ret[6] += BStart;
		}
		if (SuffixRet != null) {
			if (LocalDebug) {
				System.out.println("Suffix alignment: " + SuffixRet[0] + ", "
						+ (SuffixRet[1] + AEnd) + " to "
						+ (SuffixRet[2] + AEnd) + " skipping " + SuffixRet[5]
						+ ", and " + (SuffixRet[3] + BEnd) + " to "
						+ (SuffixRet[4] + BEnd) + " skipping " + SuffixRet[6]);
			}
			Ret[0] += SuffixRet[0];
			Ret[2] = SuffixRet[2] + AEnd;
			Ret[4] = SuffixRet[4] + BEnd;

			Ret[5] += SuffixRet[5] + (SuffixRet[1]);
			Ret[6] += SuffixRet[6] + (SuffixRet[3]);
		} else {
			if (LocalDebug)
				System.out.println("No suffix alignment");
			// Ret[5] += (SeqA_AAs.length - (AEnd));
			// Ret[6] += (SeqB.length() - (BEnd)+1);
		}

		if (LocalDebug) {
			System.out.println("Summed alignment: " + Ret[0] + ", " + Ret[1]
					+ " to " + Ret[2] + " skipping " + Ret[5] + ", and "
					+ Ret[3] + " to " + Ret[4] + " skipping " + Ret[6]);
			Utils.WaitForEnter();
		}
		return Ret;

	}

	/**
	 * Determines the best alignment, no gaps, between SeqA and SeqB
	 * 
	 * @param SeqA
	 * @param SeqB
	 * @return Returns an integer array [Score,
	 *         SeqAStart,SeqAEnd,SeqBStart,SeqBEnd, SkippedSeqA,SkippedSeqB]
	 */
	public static int[] GetBestAlignment(String SeqA, String SeqB,
			double Tolerance) {
		int[] AnchorPRMs = AntibodyUtils.GetSeqPRMScaled(SeqA);
		int[] OraclePRMs = AntibodyUtils.GetSeqPRMScaled(SeqB);
		int[][] Alignment = new int[AnchorPRMs.length][OraclePRMs.length];
		int[][][] PrevSteps = new int[AnchorPRMs.length][OraclePRMs.length][2];

		// System.out.println("SeqA length: " + AnchorPRMs.length);
		// System.out.println("SeqB length: " + OraclePRMs.length);
		boolean LocalDebug = false;

		PrevSteps[0][0][0] = -1;
		PrevSteps[0][0][1] = -1;
		if (LocalDebug)
			System.out.println("Aligning: " + SeqA + " to " + SeqB);
		int MaxScore = 0;
		int[] BestCells = { 0, 0 };
		for (int i = 1; i < AnchorPRMs.length; ++i) {
			for (int j = 1; j < OraclePRMs.length; ++j) {
				PrevSteps[i][j][0] = -1;
				PrevSteps[i][j][1] = -1;

				for (int Previ = i - 1; Previ >= Math.max(0, i - 3); --Previ) {
					for (int Prevj = j - 1; Prevj >= Math.max(0, j - 3); --Prevj) {
						int Diff1 = AnchorPRMs[Previ] - AnchorPRMs[i];
						int Diff2 = OraclePRMs[Prevj] - OraclePRMs[j];
						if (AntibodyUtils.EquivalentScaled(Diff1, Diff2,
								Tolerance)) {
							if (Alignment[Previ][Prevj] + 1 > Alignment[i][j]) {
								Alignment[i][j] = Alignment[Previ][Prevj] + 1;
								PrevSteps[i][j][0] = Previ;
								PrevSteps[i][j][1] = Prevj;
							}
						}
					}
				}

				if (Alignment[i][j] > MaxScore) {
					MaxScore = Alignment[i][j];
					BestCells[0] = i;
					BestCells[1] = j;

				}
			}
		}

		int SkippedA = 0;
		int SkippedB = 0;
		int AStart = 0;
		int AEnd = 0;
		int BStart = 0;
		int BEnd = 0;

		if (LocalDebug)
			System.out.println("Best end is at position " + BestCells[1]
					+ " with score " + MaxScore);

		int Curri = BestCells[0];
		int Currj = BestCells[1];

		int Previ = PrevSteps[Curri][Currj][0];
		int Prevj = PrevSteps[Curri][Currj][1];

		// SkippedA += (AnchorPRMs.length - Curri - 1);
		// SkippedB += (OraclePRMs.length - Currj - 1);

		AEnd = BestCells[0];
		BEnd = BestCells[1];

		while (Previ != -1 && Prevj != -1) {

			SkippedA += Math.max(0, (Curri - Previ - 1));
			SkippedB += Math.max(0, (Currj - Prevj - 1));

			Curri = Previ;
			Currj = Prevj;
			if (Previ == 0 || Prevj == 0)
				break;

			Previ = PrevSteps[Curri][Currj][0];
			Prevj = PrevSteps[Curri][Currj][1];
		}
		AStart = Curri;
		BStart = Currj;

		if (LocalDebug) {
			System.out.println("Alignment of " + SeqA + " to " + SeqB + ":"
					+ MaxScore);
			System.out.println("A: " + AStart + "-" + AEnd + " skipped "
					+ SkippedA);
			System.out.println("B: " + BStart + "-" + BEnd + " skipped "
					+ SkippedB);
		}
		int[] Ret = { MaxScore, AStart, AEnd, BStart, BEnd, SkippedA, SkippedB };
		return Ret;

	}

	public static int GetFullSequenceCorrectness(String FinalSeq,
			double Tolerance) {
		String[] Anchors = FinalSeq.split("-");
		// boolean[] Covered = new
		// boolean[AntibodyUtils.aBTLA_OracleSequence.length() + 1];
		int IncorrectPRMCount = 0;
		int CorrectPRMCount = 0;
		boolean LocalDebug = false;
		// int[] OraclePRMs =
		// AntibodyUtils.GetSeqPRMScaled(AntibodyUtils.aBTLA_OracleSequence);
		int SkipCount = 0;
		for (int a = 0; a < Anchors.length; ++a) {
			// while(SkipCount < AntibodyUtils.MAX_BEST_PATHS)
			// {
			if (LocalDebug)
				System.out.println("Examining anchor " + Anchors[a]);
			String[] AnchorPRMs = AntibodyUtils.SplitStringtoAA(Anchors[a]);
			// int[] PathInfo =
			// AntibodyUtils.GetBestAlignmentRecursive(Anchors[a],
			// AntibodyUtils.BSA);
			int[] PathInfo = AntibodyUtils.GetBestAlignmentRecursive(
					Anchors[a], AntibodyUtils.aBTLA_OracleSequence, Tolerance);

			if (PathInfo == null
					|| PathInfo[0] < Math.min(Anchors[a].length(),
							AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL)) {
				System.out.println("WARNING: Anchor " + Anchors[a]
						+ " not found on Oracle!!");
				// AntibodyUtils.WaitForEnter();
				IncorrectPRMCount += AnchorPRMs.length;
				continue;
			}
			CorrectPRMCount += PathInfo[2] - PathInfo[1] - PathInfo[5];
			IncorrectPRMCount += PathInfo[5] + PathInfo[1]
					+ (AnchorPRMs.length - PathInfo[2]);

			/*
			 * //If there is significant anchor un-aligned preceding the alinged
			 * part if(PathInfo[1] > AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL
			 * && PathInfo[3] > AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL) {
			 * 
			 * String PrefixAnchor = ""; for(int i = 0; i < PathInfo[1]; ++i) {
			 * PrefixAnchor += AnchorPRMs[i]; } if(LocalDebug) {
			 * System.out.println("WOrking with the prefix of anchor: " +
			 * PrefixAnchor);
			 * 
			 * System.out.println("WOrking with the prefix of Oracle: " +
			 * AntibodyUtils.aBTLA_OracleSequence.substring(0, PathInfo[3])); }
			 * SkipCount += 1; int [] PrefixPathInfo =
			 * AntibodyUtils.GetBestAlignment
			 * (PrefixAnchor,AntibodyUtils.aBTLA_OracleSequence.substring(0,
			 * PathInfo[3])); if(PrefixPathInfo[0] <
			 * Math.min(PathInfo[1],AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL))
			 * { System.out.println("WARNING: Anchor " + PrefixAnchor +
			 * " not found on Oracle!!"); //AntibodyUtils.WaitForEnter();
			 * IncorrectPRMCount += PathInfo[1]; if(LocalDebug)
			 * System.out.println("No good prefix path found, adding " +
			 * PathInfo[1] + " incorrect"); } else { if(LocalDebug)
			 * System.out.println("Good prefix path found, adding " +
			 * (PrefixPathInfo[2]-PrefixPathInfo[1] - PrefixPathInfo[5]) +
			 * " correct and " + (PrefixPathInfo[5] + PrefixPathInfo[1] +
			 * (PathInfo[1] - PrefixPathInfo[2])) + " incorrect");
			 * CorrectPRMCount += PrefixPathInfo[2]-PrefixPathInfo[1] -
			 * PrefixPathInfo[5]; //PathEndinAnchor - PathStartinAnchor -
			 * SKippedinAnchor IncorrectPRMCount += PrefixPathInfo[5] +
			 * PrefixPathInfo[1] + (PathInfo[1] - PrefixPathInfo[2]);
			 * //SkippedinAnchor + PathStartinAnchor + (Gap) }
			 * 
			 * } if(AnchorPRMs.length - PathInfo[2]>
			 * AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL &&
			 * AntibodyUtils.aBTLA_OracleSequence.length() - PathInfo[4] >
			 * AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL) { String PrefixAnchor
			 * = ""; int Length = 0; for(int i = PathInfo[2]; i <
			 * AnchorPRMs.length; ++i) { PrefixAnchor += AnchorPRMs[i]; Length
			 * += 1;
			 * 
			 * } if(LocalDebug) {
			 * System.out.println("WOrking with the suffix of anchor: " +
			 * PrefixAnchor);
			 * System.out.println("WOrking with the suffix of Oracle: " +
			 * AntibodyUtils.aBTLA_OracleSequence.substring(PathInfo[4])); }
			 * SkipCount += 1; int[] SuffixPathInfo =
			 * AntibodyUtils.GetBestAlignment(PrefixAnchor,
			 * AntibodyUtils.aBTLA_OracleSequence.substring(PathInfo[4]));
			 * if(SuffixPathInfo[0] <
			 * Math.min(Length,AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL)) {
			 * System.out.println("WARNING: Anchor " +
			 * Anchors[a].substring(PathInfo[2]) + " not found on Oracle!!");
			 * AntibodyUtils.WaitForEnter(); IncorrectPRMCount +=
			 * AnchorPRMs.length - PathInfo[2] - 1; if(LocalDebug)
			 * System.out.println("No good suffix path found, adding " +
			 * (AnchorPRMs.length - PathInfo[2] - 1) + " incorrect"); } else {
			 * if(LocalDebug)
			 * System.out.println("Good suffix path found, adding " +
			 * (SuffixPathInfo[2]-SuffixPathInfo[1] - SuffixPathInfo[5]) +
			 * " correct and " + (SuffixPathInfo[5] + (Length -
			 * SuffixPathInfo[2]) + SuffixPathInfo[1]) + " incorrect");
			 * CorrectPRMCount += SuffixPathInfo[2]-SuffixPathInfo[1] -
			 * SuffixPathInfo[5]; //PathEndAnchor - PathStartAnchor -
			 * SkippedAnchor IncorrectPRMCount += SuffixPathInfo[5] + (Length -
			 * SuffixPathInfo[2]) + SuffixPathInfo[1]; //SkippedAnchor +
			 * (Remaining Suffix) + Gap } }
			 */
			System.out.println("Corrects: " + CorrectPRMCount + ", Incorrect: "
					+ IncorrectPRMCount + ", Total: "
					+ (CorrectPRMCount + IncorrectPRMCount));
			if (LocalDebug) {
				System.out.println("Corrects: " + CorrectPRMCount
						+ ", Incorrect: " + IncorrectPRMCount + ", Total: "
						+ (CorrectPRMCount + IncorrectPRMCount));
				Utils.WaitForEnter();
			}
		}
		return CorrectPRMCount;

	}

	public static int GetFullSequenceCorrectnessLC(String FinalSeq,
			double Tolerance) {
		String[] Anchors = FinalSeq.split("-");
		// boolean[] Covered = new boolean[AntibodyUtils.aBTLA_LC.length() + 1];
		int IncorrectPRMCount = 0;
		int CorrectPRMCount = 0;
		boolean LocalDebug = false;
		// int[] OraclePRMs =
		// AntibodyUtils.GetSeqPRMScaled(AntibodyUtils.aBTLA_LC);

		int SkipCount = 0;
		for (int a = 0; a < Anchors.length; ++a) {
			// while(SkipCount < AntibodyUtils.MAX_BEST_PATHS)
			// {
			if (LocalDebug)
				System.out.println("Examining anchor " + Anchors[a]);
			String[] AnchorPRMs = AntibodyUtils.SplitStringtoAA(Anchors[a]);
			int[] PathInfo = AntibodyUtils.GetBestAlignmentRecursive(
					Anchors[a], AntibodyUtils.aBTLA_LC, Tolerance);
			if (PathInfo[0] < Math.min(Anchors[a].length(),
					AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL)) {
				System.out.println("WARNING: Anchor " + Anchors[a]
						+ " not found on Oracle!!");
				// AntibodyUtils.WaitForEnter();
				continue;
			}
			CorrectPRMCount += PathInfo[2] - PathInfo[1] - PathInfo[5];
			IncorrectPRMCount += PathInfo[5] + PathInfo[1]
					+ (AnchorPRMs.length - PathInfo[2]);

			if (LocalDebug) {
				System.out.println("Corrects: " + CorrectPRMCount
						+ ", Incorrect: " + IncorrectPRMCount + ", Total: "
						+ (CorrectPRMCount + IncorrectPRMCount));
				Utils.WaitForEnter();
			}
		}
		return CorrectPRMCount;

	}

	/**
	 * We treat A as teh base sequence, and string B starts at some position
	 * after A.
	 * 
	 * @param A
	 *            the base sequence
	 * @param B
	 *            the sequence to merge with the base
	 * @return Returns teh updated base sequence
	 */
	public static String MergeSequences(String A, String B, double Tolerance) {
		// int ShiftIndex = 0;
		//
		// System.out.println("Merging: " + A + " and " + B);
		boolean LocalDebug = false;

		int[] PRMA = AntibodyUtils.GetSeqPRMScaled(A);
		int[] PRMB = AntibodyUtils.GetSeqPRMScaled(B);

		String[] StringAParts = AntibodyUtils.SplitStringtoAA(A);
		String[] StringBParts = AntibodyUtils.SplitStringtoAA(B);

		int[][] Alignment = new int[PRMA.length][PRMB.length];
		int[][][] PrevSteps = new int[PRMA.length][PRMB.length][2];

		String StringPRMA = "";
		String StringPRMB = "";
		for (int i = 0; i < PRMA.length; ++i) {
			StringPRMA += PRMA[i] + " ";
			Alignment[i][0] = 1;
			PrevSteps[i][0][0] = -1;
			PrevSteps[i][0][1] = -1;
		}
		System.out.println("PRMA: " + StringPRMA);

		for (int i = 0; i < PRMB.length; ++i) {
			StringPRMB += PRMB[i] + " ";
			Alignment[0][i] = 1;
			PrevSteps[0][i][0] = -1;
			PrevSteps[0][i][1] = -1;
		}

		// System.out.println("PRMB: " + StringPRMB);
		for (int i = 1; i < PRMA.length; ++i) {
			for (int j = 1; j < PRMB.length; ++j) {
				int MaxScore = 0;
				PrevSteps[i][j][0] = -1;
				PrevSteps[i][j][1] = -1;

				if (LocalDebug && !LocalDebug) {
					System.out.println("Considering " + PRMA[i] + " to "
							+ PRMB[j]);
				}

				// Consider matching PRMA[i] with PRMB[j], look at all previous
				// matchings
				for (int Previ = 0; Previ < i; Previ++) {
					for (int Prevj = 0; Prevj < j; Prevj++) {
						int Diff1 = PRMA[i] - PRMA[Previ];
						int Diff2 = PRMB[j] - PRMB[Prevj];

						if (AntibodyUtils.EquivalentScaled(Diff1, Diff2,
								Tolerance)) {
							if (Alignment[Previ][Prevj] + 1 > MaxScore) {
								MaxScore = Alignment[Previ][Prevj] + 1;
								PrevSteps[i][j][0] = Previ;
								PrevSteps[i][j][1] = Prevj;
							}
						}
					}
				}
				if (LocalDebug && !LocalDebug) {
					if (MaxScore > 1)
						System.out.println("MaxScore: " + MaxScore + " from ["
								+ PRMA[PrevSteps[i][j][0]] + ","
								+ PRMB[PrevSteps[i][j][1]] + "]");
					else
						System.out.println("MaxScore: " + MaxScore
								+ ",but no prev");

				}
				Alignment[i][j] = MaxScore;
			}
		}

		// Find the best path
		int MaxPath = -1;
		int[] MaxEntry = new int[2];
		String FinalString = "";

		for (int i = 0; i < PRMA.length; ++i) {
			if (Alignment[i][PRMB.length - 1] > MaxPath) {
				MaxPath = Alignment[i][PRMB.length - 1];
				MaxEntry[0] = i;
				MaxEntry[1] = PRMB.length - 1;
				FinalString = "";
				for (int j = i; j < StringAParts.length; ++j)
					FinalString += StringAParts[j];
			}
		}

		for (int i = 0; i < PRMB.length; ++i) {
			if (Alignment[PRMA.length - 1][i] >= MaxPath) {
				MaxPath = Alignment[PRMA.length - 1][i];
				MaxEntry[0] = PRMA.length - 1;
				MaxEntry[1] = i;
				FinalString = "";
				for (int j = i; j < StringBParts.length; ++j)
					FinalString += StringBParts[j];
			}
		}
		if (MaxPath <= AntibodyUtils.MIN_OVERLAPPING_AA)
			return A + B;
		if (LocalDebug && !LocalDebug) {
			System.out.println("Best path length: " + MaxPath + " ends at ["
					+ MaxEntry[0] + "," + MaxEntry[1] + "]");
			System.out.println("Suffix: " + FinalString);
		}

		int Curri = MaxEntry[0];
		int Currj = MaxEntry[1];
		int Previ = 0;
		int Prevj = 0;
		for (int i = 1; i < MaxPath; ++i) {
			Previ = PrevSteps[Curri][Currj][0];
			Prevj = PrevSteps[Curri][Currj][1];
			int Stepsi = Curri - Previ;
			int Stepsj = Currj - Prevj;

			if (LocalDebug)
				System.out.println("To [" + Curri + "," + Currj + "] from ["
						+ Previ + "," + Prevj + "]");

			if (Stepsi == Stepsj || Stepsi > Stepsj) {
				if (LocalDebug)
					System.out
							.println("Took an even number of steps, or more A than B");
				for (int j = Curri - 1; j >= Previ; --j) {
					FinalString = StringAParts[j] + FinalString;
				}
				if (LocalDebug)
					System.out.println(FinalString);
			}

			else {
				if (LocalDebug)
					System.out.println("Took more B steps");
				for (int j = Currj - 1; j >= Prevj; --j) {
					FinalString = StringBParts[j] + FinalString;
				}
				if (LocalDebug)
					System.out.println(FinalString);
			}
			if (LocalDebug)
				Utils.WaitForEnter();
			Curri = Previ;
			Currj = Prevj;
		}
		while (Curri > 0) {
			FinalString = StringAParts[Curri - 1] + FinalString;
			Curri--;
		}
		while (Currj > 0) {
			FinalString = StringBParts[Currj - 1] + FinalString;
			Currj--;
		}
		return FinalString;
	}

	public static String MergeSequences2(String A, String B, double Tolerance) {

		int ShiftIndex = 0;

		int[] PRMA = AntibodyUtils.GetSeqPRMScaled(A);
		int[] PRMB = AntibodyUtils.GetSeqPRMScaled(B);

		int Watch = 5;
		int MaxMatching = 0;
		int MaxShift = 0;
		String MaxString = "";

		String[] StringAParts = AntibodyUtils.SplitStringtoAA(A);
		String[] StringBParts = AntibodyUtils.SplitStringtoAA(B);

		while (ShiftIndex < PRMA.length - 1) {
			int IndexA = ShiftIndex + 1;
			int IndexB = 1;
			int CurrMatching = 0;
			String CurrString = A.substring(0, ShiftIndex + 1);

			while (IndexB < PRMB.length) {
				int IntStartA = -1;
				int IntEndA = -1;
				if (IndexA < PRMA.length) {
					IntStartA = PRMA[IndexA - 1];
					IntEndA = PRMA[IndexA];
				}

				int IntStartB = PRMB[IndexB - 1] + PRMA[ShiftIndex];
				int IntEndB = PRMB[IndexB] + PRMA[ShiftIndex];

				if (IndexA < PRMA.length && IntStartA == IntStartB
						&& IntEndA == IntEndB) {
					int Diff = IntEndB - IntStartB;

					// char AA = Utils.FindClosestAAScaled(Diff);
					// if(AA == '0')
					// CurrString += "[" + (((double)Diff)/Utils.MASS_SCALE) +
					// "]";
					// else
					// CurrString += AA;
					CurrString += StringAParts[IndexA - 1];
					IndexA++;
					IndexB++;
					CurrMatching += 1;
				} else if (IndexA < PRMA.length && IntStartA == IntStartB) {
					int Diff = Math.min(IntEndA, IntEndB) - IntStartA;
					char AA = AntibodyUtils
							.FindClosestAAScaled(Diff, Tolerance);
					if (AA == '0')
						CurrString += "["
								+ (((double) Diff) / AntibodyUtils.MASS_SCALE)
								+ "]";
					else
						CurrString += AA;
					if (IntEndA < IntEndB)

						IndexA++;
					else
						IndexB++;
				} else if (IndexA < PRMA.length && IntEndA == IntEndB) {
					int Diff = IntEndA - Math.max(IntStartA, IntStartB);
					char AA = AntibodyUtils
							.FindClosestAAScaled(Diff, Tolerance);
					if (AA == '0')
						CurrString += "["
								+ (((double) Diff) / AntibodyUtils.MASS_SCALE)
								+ "]";
					else
						CurrString += AA;
					IndexA++;
					IndexB++;
					CurrMatching += 1;
				} else if (IndexA < PRMA.length && IntStartA != IntStartB
						&& IntEndA != IntEndB) {

					int Diff = Math.max(IntStartA, IntStartB)
							- Math.min(IntStartA, IntStartB);
					char AA = AntibodyUtils
							.FindClosestAAScaled(Diff, Tolerance);
					if (AA == '0')
						CurrString += "["
								+ (((double) Diff) / AntibodyUtils.MASS_SCALE)
								+ "]";
					else
						CurrString += AA;
					if (IntStartA < IntStartB)
						IndexA++;
					else
						IndexB++;

				} else {
					int Diff = IntEndB - IntStartB;
					char AA = AntibodyUtils
							.FindClosestAAScaled(Diff, Tolerance);
					if (AA == '0')
						CurrString += "["
								+ (((double) Diff) / AntibodyUtils.MASS_SCALE)
								+ "]";
					else
						CurrString += AA;
					IndexB++;
				}

			}
			if (CurrMatching > MaxMatching) {
				MaxMatching = CurrMatching;
				MaxShift = ShiftIndex;
				MaxString = CurrString;

			}
			ShiftIndex++;
		}

		return MaxString;

	}

	public static String[] SplitStringtoAA(String StringA) {

		if (StringA == null) {
			System.err
					.println("ERROR: SplitStringToAA: Unable to split a null string!");
			return null;
		}
		ArrayList Temp = new ArrayList();

		int i = 0;
		while (i < StringA.length()) {
			if (StringA.charAt(i) == '[') {
				int End = StringA.indexOf(']', i);
				String s = StringA.substring(i, End + 1);
				i = End + 1;
				Temp.add(s);
			} else {
				Temp.add((new Character(StringA.charAt(i))).toString());
				i++;
			}
		}
		String[] Parts = new String[Temp.size()];
		for (i = 0; i < Parts.length; ++i)
			Parts[i] = (String) (Temp.get(i));
		return Parts;
	}

	public static String[] LoadLinesFromFile(String inputFile) {
		BufferedReader Buf = null;
		String Line = null;
		try {
			Buf = new BufferedReader(new FileReader(inputFile));
			Line = Buf.readLine();
		} catch (IOException E) {
			System.err.println(E);
			inputFile = null;
			return null;
		}

		ArrayList InputFiles = new ArrayList();
		while (Line != null) {
			Line = Line.trim();
			if (Line.length() > 0) {
				InputFiles.add(Line);
			}
			try {
				Line = Buf.readLine();
			} catch (IOException E) {
				System.err.println(E);
			}

		}

		try {
			Buf.close();
		} catch (IOException E) {
			System.err.println(E);
		}
		return Utils.ConvertArraylistToStringArray(InputFiles);
	}

	public static void FindAllSubsetsHelper(int[] Els, int k, int[] Prefix,
			ArrayList Results) {
		if (k == 0 && Results.size() < 10000)
			Results.add(Prefix);
		else if (Results.size() < 10000) {
			for (int i = 0; i < Els.length; i++) {
				int[] ReducedEls = new int[Els.length - i - 1];
				for (int j = i + 1; j < Els.length; ++j)
					ReducedEls[j - i - 1] = Els[j];
				int[] ExtendedPrefix = new int[Prefix.length + 1];
				for (int j = 0; j < Prefix.length; ++j)
					ExtendedPrefix[j] = Prefix[j];
				ExtendedPrefix[Prefix.length] = Els[i];
				FindAllSubsetsHelper(ReducedEls, k - 1, ExtendedPrefix, Results);
			}
		}
	}

	public static ArrayList FindAllSubsets(int n, int k) {

		// System.out.println("Find All Subsets " + n + "," + k);
		ArrayList Results = new ArrayList();
		int[] Els = new int[n];
		for (int i = 0; i < n; ++i)
			Els[i] = i;
		int[] Prefix = {};
		FindAllSubsetsHelper(Els, k, Prefix, Results);

		return Results;
	}

	public static String StringArrayToString(String[] input) {
		String Ret = "";
		for (int i = 0; i < input.length; ++i)
			Ret += input[i] + " ";
		return Ret;
	}

	public static int InsertInDecreasingOrder(double Value, double[] List) {
		int Index = 0;
		for (int i = 0; i < List.length; ++i) {
			if (Value > List[i]) {
				Index = i;
				for (int j = List.length - 1; j > i; --j) {
					List[j] = List[j - 1];
				}
				List[i] = Value;
				return Index;
			}
		}
		return Index;
	}

	public static boolean ContainsMass(int[] List, int Mass, double Tolerance) {
		for (int i = 0; i < List.length; ++i) {
			if (AntibodyUtils.EquivalentScaled(List[i], Mass, Tolerance))
				return true;
		}
		return false;
	}

	public static int MassIndex(int[] List, int Mass, double Tolerance) {
		for (int i = 0; i < List.length; ++i) {
			if (AntibodyUtils.EquivalentScaled(List[i], Mass, Tolerance))
				return i;
		}
		return -1;
	}

	/**
	 * Finds the index of the value in the List that is closest to the query
	 * Mass, within the defined Tolerance
	 * 
	 * @param List
	 * @param Mass
	 * @param Tolerance
	 * @return
	 */
	public static int MassIndexClosest(int[] List, int Mass, double Tolerance) {
		int currBestDistance = Integer.MAX_VALUE;
		int currBestIndex = -1;
		for (int i = 0; i < List.length; ++i) {
			if (AntibodyUtils.EquivalentScaled(List[i], Mass, Tolerance)
					&& Math.abs(List[i] - Mass) < currBestDistance) {
				currBestDistance = Math.abs(List[i] - Mass);
				currBestIndex = i;
			}
		}
		return currBestIndex;
	}

	/**
	 * Finds the index of the value in the List that is closest to the query
	 * Mass
	 * 
	 * @param List
	 * @param Mass
	 * @param Tolerance
	 * @return
	 */
	public static int MassIndexClosest(int[] List, int Mass) {
		return AntibodyUtils.MassIndexClosest(List, Mass, Integer.MAX_VALUE);
	}

	public static int MassIndexExact(int[] List, int Mass) {
		for (int i = 0; i < List.length; ++i) {
			if (List[i] == Mass)
				return i;
		}
		return -1;
	}

	public static boolean ContainsString(ArrayList List, String Item) {
		for (int i = 0; i < List.size(); ++i) {
			String CurrItem = (String) (List.get(i));
			if (CurrItem.compareTo(Item) == 0)
				return true;
		}
		return false;
	}

	public static int FindORFPosition(ArrayList ORFList, String ORFName) {
		if (ORFList.size() == 0)
			return 0;

		String[] Parts = SixFrameBuilder.ParseSixFrameHeader(ORFName);

		if (Parts[3].compareTo("1") >= 0) {
			int Start = Integer.parseInt(Parts[1]);
			int End = Integer.parseInt(Parts[2]);

			for (int i = 0; i < ORFList.size(); ++i) {
				String CurrORF = (String) (ORFList.get(i));
				String[] CurrParts = SixFrameBuilder
						.ParseSixFrameHeader(CurrORF);

				int CurrStart = Integer.parseInt(CurrParts[1]);
				int CurrEnd = Integer.parseInt(CurrParts[2]);

				if (Start < CurrStart)
					return i;
				if (Start == CurrStart && End < CurrEnd)
					return i;
			}
		} else {
			int Start = Integer.parseInt(Parts[1]);
			int End = Integer.parseInt(Parts[2]);

			for (int i = 0; i < ORFList.size(); ++i) {
				String CurrORF = (String) (ORFList.get(i));
				String[] CurrParts = SixFrameBuilder
						.ParseSixFrameHeader(CurrORF);

				int CurrStart = Integer.parseInt(CurrParts[1]);
				int CurrEnd = Integer.parseInt(CurrParts[2]);

				if (End > CurrEnd)
					return i;
				if (End == CurrEnd && Start > CurrStart)
					return i;
			}

		}
		return ORFList.size();
	}

	public static boolean MassWithinPRMRange(int[] List, int Mass,
			double Tolerance) {
		if (Mass < List[List.length - 1] && Mass > List[0])
			return true;
		if (AntibodyUtils.EquivalentScaled(Mass, List[0], Tolerance))
			return true;
		return false;
	}

	public static void Remove(int Index, double[] List, double BlankValue) {
		for (int i = Index; i < List.length - 1; ++i) {
			List[i] = List[i + 1];
		}
		List[List.length - 1] = BlankValue;
	}

	public static void Remove(int Index, int[] List, int BlankValue) {
		for (int i = Index; i < List.length - 1; ++i) {
			List[i] = List[i + 1];
		}
		List[List.length - 1] = BlankValue;
	}

	public static void Insert(int Index, double Value, double[] List) {
		for (int j = List.length - 1; j > Index; --j) {
			List[j] = List[j - 1];
		}
		List[Index] = Value;
	}

	public static void Insert(int Index, int Value, int[] List) {
		for (int j = List.length - 1; j > Index; --j) {
			List[j] = List[j - 1];
		}
		List[Index] = Value;
	}

	public static void DebugPrintList(double[] List, double BlankValue) {
		String S = "";
		for (int i = 0; i < List.length; ++i) {
			if (List[i] == BlankValue) {
				System.out.println("Len: " + i + " : " + S);
				return;
			}
			S += List[i] + " ";
		}
		System.out.println(S);
	}

	public static void DebugPrintList(int[] List, int BlankValue) {
		String S = "";
		for (int i = 0; i < List.length; ++i) {
			if (List[i] == BlankValue) {
				System.out.println("Len: " + i + " : " + S);
				return;
			}
			S += List[i] + " ";
		}
		System.out.println(S);
	}

	public static int GetNumExtended(int[] SeedPRM, int[] ExtendedPRM,
			double Tolerance) {
		for (int i = SeedPRM.length - 1; i >= 0; --i) {
			int SeedMass = SeedPRM[i];
			for (int j = 0; j < ExtendedPRM.length; ++j) {
				int ExtendMass = ExtendedPRM[j];
				if (AntibodyUtils.EquivalentScaled(SeedMass, ExtendMass,
						Tolerance))
					return (ExtendedPRM.length - j - 1);
			}
		}
		return 0;
	}

	public static int GetNumExtendedLeft(int[] SeedPRM, int[] ExtendedPRM,
			double Tolerance) {
		for (int i = 0; i < SeedPRM.length; ++i) {
			int SeedMass = SeedPRM[i];
			for (int j = 0; j < ExtendedPRM.length; ++j) {
				int ExtendMass = ExtendedPRM[j];
				if (AntibodyUtils.EquivalentScaled(SeedMass, ExtendMass,
						Tolerance))
					return (j);
			}
		}
		return 0;
	}

	public static int GetNumCorrect(int[] SeedPRM, int[] ExtendedPRM,
			double Tolerance) {
		int Count = 0;
		for (int i = 0; i < SeedPRM.length; ++i) {
			int SeedMass = SeedPRM[i];
			for (int j = 0; j < ExtendedPRM.length; ++j) {
				int ExtendMass = ExtendedPRM[j];
				if (AntibodyUtils.EquivalentScaled(SeedMass, ExtendMass,
						Tolerance))
					Count += 1;
			}
		}
		return Count;
	}

	public static String GetNumCorrectString(int[] SeedPRM, int[] ExtendedPRM,
			double Tolerance) {
		boolean LocalDebug = false;
		if (LocalDebug) {
			String A = "";
			String B = "";
			for (int i = 0; i < SeedPRM.length; ++i) {
				A += SeedPRM[i] + " ";
			}
			for (int i = 0; i < ExtendedPRM.length; ++i) {
				B += ExtendedPRM[i] + " ";
			}
			System.out.println("Seed: " + A);
			System.out.println("Extension: " + B);
			Utils.WaitForEnter();
		}

		String Ret = "";
		for (int i = 0; i < ExtendedPRM.length; ++i) {
			boolean found = false;
			int ExtendMass = ExtendedPRM[i];
			for (int j = 0; j < SeedPRM.length; ++j) {
				int SeedMass = SeedPRM[j];
				if (AntibodyUtils.EquivalentScaled(SeedMass, ExtendMass,
						Tolerance)) {
					found = true;
					break;
				}
			}
			if (found)
				Ret += "1";
			else
				Ret += "0";
		}
		if (LocalDebug) {
			System.out.println(Ret);
			Utils.WaitForEnter();
		}
		return Ret;
	}

	public static String IntegerListToString(int[] List) {
		String Ret = "";
		for (int i = 0; i < List.length; ++i)
			Ret += List[i] + " ";
		return Ret;
	}

	public static int GetNumMissing(int[] SeedPRM, int[] ExtendedPRM,
			double Tolerance) {
		int Count = 0;
		for (int i = 0; i < SeedPRM.length; ++i) {
			int SeedMass = SeedPRM[i];
			boolean Found = false;
			for (int j = 0; j < ExtendedPRM.length; ++j) {
				int ExtendMass = ExtendedPRM[j];
				if (AntibodyUtils.EquivalentScaled(SeedMass, ExtendMass,
						Tolerance)) {
					Found = true;
					break;
				}
			}
			if (!Found)
				Count += 1;
		}
		return Count;
	}

	public static double AlignmentScoreLengthWithMassShiftRight(
			MSAlignmentType A, MSAlignmentType B, double Tolerance) {
		double Score = 0.0;
		int[] NewPRMsA = new int[A.GetSpectrum().PRMs.length];
		int[] NewPRMsB = new int[B.GetSpectrum().PRMs.length];

		// int MassShiftAPos = A.GetAnchorSequence().indexOf(A.GetSequence());
		int MassShiftA = A.GetScaledMassShift();

		// int MassShiftBPos = B.GetAnchorSequence().indexOf(B.GetSequence());
		int MassShiftB = B.GetScaledMassShift();

		// System.out.println(A.GetSequence() + ":" + A.GetScore() + "-" +
		// MassShiftA);
		// System.out.println(B.GetSequence() + ":" + B.GetScore() + "-" +
		// MassShiftB);

		for (int i = 0; i < NewPRMsA.length; ++i)
			NewPRMsA[i] = A.GetSpectrum().PRMs[i] + MassShiftA;

		for (int i = 0; i < NewPRMsB.length; ++i)
			NewPRMsB[i] = B.GetSpectrum().PRMs[i] + MassShiftB;

		for (int i = 0; i < NewPRMsA.length; ++i) {
			for (int j = 0; j < NewPRMsB.length; ++j) {
				if (AntibodyUtils.EquivalentScaled(NewPRMsA[i], NewPRMsB[j],
						Tolerance))
					Score += 1;

			}
		}
		return Score;
	}

	public static double AlignmentScoreLengthWithMassShiftLeft(
			MSAlignmentType A, MSAlignmentType B, double Tolerance) {
		double Score = 0.0;
		boolean LocalDebug = false;
		int[] NewPRMsA = new int[A.GetSpectrum().PRMs.length];
		int[] NewPRMsB = new int[B.GetSpectrum().PRMs.length];

		int[] AlignedPeaksA = A.GetAlignedPeaksScaled();
		int[] AlignedPeaksB = B.GetAlignedPeaksScaled();

		// int MassShiftAPos = A.GetAnchorSequence().indexOf(A.GetSequence()) +
		// A.GetSequence().length();
		int MassShiftA = A.GetScaledMassShift();
		// if(MassShiftAPos < A.GetAnchorSequence().length())
		// MassShiftA +=
		// AntibodyUtils.GetSeqMassScaled(A.GetAnchorSequence().substring(MassShiftAPos));
		// int MassShiftBPos = B.GetAnchorSequence().indexOf(B.GetSequence()) +
		// B.GetSequence().length();
		int MassShiftB = B.GetScaledMassShift();
		// if(MassShiftBPos < B.GetAnchorSequence().length())
		// MassShiftB +=
		// AntibodyUtils.GetSeqMassScaled(B.GetAnchorSequence().substring(MassShiftBPos));
		// int MassShiftB =
		// Utils.GetSeqMassScaled(B.GetAnchorSequence().substring(0,MassShiftBPos))-B.GetAlignedPeaksScaled()[0];

		// System.out.println(A.GetSequence() + ":" + A.GetScore() + "-" +
		// MassShiftA);
		// System.out.println(B.GetSequence() + ":" + B.GetScore() + "-" +
		// MassShiftB);

		for (int i = 0; i < NewPRMsA.length; ++i)
			NewPRMsA[i] = A.GetSpectrum().PRMs[i] - MassShiftA;

		for (int i = 0; i < NewPRMsB.length; ++i)
			NewPRMsB[i] = B.GetSpectrum().PRMs[i] - MassShiftB;

		for (int i = 0; i < NewPRMsA.length; ++i) {
			for (int j = 0; j < NewPRMsB.length; ++j) {
				if (AntibodyUtils.EquivalentScaled(NewPRMsA[i], NewPRMsB[j],
						Tolerance))
					Score += 1;

			}
		}
		return Score;
	}

	public static int CompareToPRMs(String A, String B, double Tolerance) {
		boolean LocalDebug = false;
		if (A.length() != B.length()) {
			if (LocalDebug)
				System.out.println("CompareToPRMs: A length " + A.length()
						+ " != B length " + B.length());
			return -1;
		}
		int[] PRMsA = AntibodyUtils.GetSeqPRMScaled(A);
		int[] PRMsB = AntibodyUtils.GetSeqPRMScaled(B);

		if (PRMsA.length != PRMsB.length) {
			if (LocalDebug)
				System.out.println("CompareToPRMs: PRMA length " + PRMsA.length
						+ " != PRMB length " + PRMsB.length);

			return 1;
		}
		for (int i = 0; i < PRMsA.length; ++i) {
			if (!AntibodyUtils.EquivalentScaled(PRMsA[i], PRMsB[i], Tolerance)) {
				if (LocalDebug)
					System.out.println("CompareToPRMs: A[" + i + "] = "
							+ PRMsA[i] + " != B[" + i + "] = " + PRMsB[i]);
				return 1;
			}
		}
		return 0;

	}

	/**
	 * Returns 2 2-D arrays containing the DP tables of alignment scores and
	 * backtracking.
	 * 
	 * @param A
	 * @param B
	 * @param Tolerance
	 * @return
	 */
	public static Object[] GetLongestAlignmentStrictDetailedSeqs(String A,
			String B, double Tolerance) {
		int[] PRMsA = AntibodyUtils.GetSeqPRMScaled(A);
		int[] PRMsB = AntibodyUtils.GetSeqPRMScaled(B);

		return AntibodyUtils.GetLongestAlignmentStrictDetailedPRMs(PRMsA,
				PRMsB, Tolerance);
	}

	public static Object[] GetLongestAlignmentStrictDetailedPRMs(int[] PRMsA,
			int[] PRMsB, double Tolerance) {
		int[][] Alignment = new int[PRMsA.length][PRMsB.length];
		int[][] prevCell = new int[PRMsA.length][PRMsB.length];
		boolean LocalDebug = false;
		for (int a = 0; a < PRMsA.length; ++a) {
			Alignment[a][0] = 1;
			prevCell[a][0] = -1;
		}
		for (int b = 0; b < PRMsB.length; ++b) {
			Alignment[0][b] = 1;
			prevCell[0][b] = -1;
		}
		for (int a = 1; a < PRMsA.length; ++a) {

			int AMass = PRMsA[a];
			for (int b = 1; b < PRMsB.length; ++b) {
				int BMass = PRMsB[b];
				int BestScore = 1;

				for (int prevA = Math.max(0, a - 2); prevA < a; ++prevA) {
					int prevAMass = PRMsA[prevA];
					for (int prevB = Math.max(0, b - 2); prevB < b; ++prevB) {
						int prevBMass = PRMsB[prevB];

						if (AntibodyUtils.EquivalentScaled(BMass - prevBMass,
								AMass - prevAMass, Tolerance)) {
							if (LocalDebug) {
								System.out.println("Prevs [" + prevAMass + ","
										+ prevBMass + "] -> [" + AMass + ","
										+ BMass + "]");
							}
							if (Alignment[prevA][prevB] + 1 > BestScore) {
								BestScore = Alignment[prevA][prevB] + 1;
								prevCell[a][b] = prevA * PRMsB.length + prevB;
							}

						}
					}

				}
				Alignment[a][b] = BestScore;

			}
		}
		Object[] ret = new Object[2];
		ret[0] = prevCell;
		ret[1] = Alignment;
		return ret;
	}

	/**
	 * Determines the best (longest) alignment of 2 PRM spectra. It is
	 * considered relaxed because the alignment can start from anywhere, and any
	 * number of peaks can be skipped.
	 * 
	 * @param PRMsA
	 * @param PRMsB
	 * @param Tolerance
	 * @return
	 */
	public static Object[] GetLongestAlignmentRelaxedDetailedPRMs(int[] PRMsA,
			int[] PRMsB, double Tolerance) {
		return AntibodyUtils.GetLongestAlignmentRelaxedDetailedPRMs(PRMsA,
				PRMsB, Tolerance, false);
	}

	/**
	 * Determines the best (longest) alignment of 2 PRM spectra. It is
	 * considered relaxed because the alignment can start from anywhere, and any
	 * number of peaks can be skipped.
	 * 
	 * @param PRMsA
	 * @param PRMsB
	 * @param Tolerance
	 * @return
	 */
	public static Object[] GetLongestAlignmentRelaxedDetailedPRMs(int[] PRMsA,
			int[] PRMsB, double Tolerance, boolean LocalDebug) {
		return AntibodyUtils.GetLongestAlignmentDetailedPRMs(PRMsA, PRMsB,
				Tolerance, Integer.MAX_VALUE, Integer.MAX_VALUE);
	}

	/**
	 * Computes the length of the best alignment between the two strings. Single
	 * PRM skips are allowed The length of the best alignment (the number of
	 * aligned PRMs) is returned
	 * 
	 * @param A
	 * @param B
	 * @return Returns an 3-tuple containing the length of the alignment, the
	 *         start of the alignment in A, the start of the alignment in B
	 */
	public static int[] GetLongestAlignmentStrictSeqs(String A, String B,
			double Tolerance) {
		int[] PRMsA = AntibodyUtils.GetSeqPRMScaled(A);
		int[] PRMsB = AntibodyUtils.GetSeqPRMScaled(B);
		return AntibodyUtils.GetLongestAlignmentStrictPRMs(PRMsA, PRMsB,
				Tolerance);
	}

	/**
	 * Computes the length of the best (longest) alignment between the two PRM
	 * lsits. Single PRM skips are allowed The length of the best alignment (the
	 * number of aligned PRMs) is returned
	 * 
	 * @param A
	 * @param B
	 * @return Returns an 3-tuple containing the length of the alignment, the
	 *         start of the alignment in A, the start of the alignment in B
	 */
	public static int[] GetLongestAlignmentStrictPRMs(int[] PRMsA, int[] PRMsB,
			double Tolerance) {
		int[][] Alignment = null;
		int[][] prevCell = null;
		Object[] details = AntibodyUtils.GetLongestAlignmentStrictDetailedPRMs(
				PRMsA, PRMsB, Tolerance);
		Alignment = (int[][]) (details[1]);
		prevCell = (int[][]) (details[0]);
		int OverallBestScore = 0;
		int BestA = 0, BestB = 0;
		for (int i = 0; i < Alignment.length; ++i) {
			for (int j = 0; j < Alignment[0].length; ++j) {
				if (Alignment[i][j] > OverallBestScore) {
					OverallBestScore = Alignment[i][j];
					BestA = i;
					BestB = j;
				}
			}
		}
		int[] ret = new int[3];
		ret[0] = OverallBestScore;
		int prevA = BestA;
		int prevB = BestB;
		while (prevCell[prevA][prevB] >= 0) {
			int currCell = prevCell[prevA][prevB];
			prevA = currCell / PRMsB.length;
			prevB = currCell % PRMsB.length;
		}
		ret[1] = prevA;
		ret[2] = prevB;
		return ret;
	}

	public static double AlignmentScore(PRMSpectrum A, PRMSpectrum B,
			double Tolerance) {
		double Score = 0.0;

		double[][] Alignment = new double[A.PRMs.length][B.PRMs.length];
		for (int a = 0; a < A.PRMs.length; ++a) {
			double AScore = A.PRMScores[a];
			int AMass = A.PRMs[a];
			for (int b = 0; b < B.PRMs.length; ++b) {
				double BScore = B.PRMScores[b];
				int BMass = B.PRMs[b];

				double MaxScore = AScore + BScore;
				for (int prevA = 0; prevA < a; ++prevA) {
					int prevAMass = A.PRMs[prevA];
					for (int prevB = 0; prevB < b; ++prevB) {
						int prevBMass = B.PRMs[prevB];

						if (AntibodyUtils.EquivalentScaled(BMass - prevBMass,
								AMass - prevAMass, Tolerance)) {
							double NewScore = AScore + BScore
									+ Alignment[prevA][prevB];
							if (NewScore > MaxScore)
								MaxScore = NewScore;

						}
					}

				}
				Alignment[a][b] = MaxScore;
			}
		}

		for (int a = 0; a < A.PRMs.length; ++a) {
			if (Alignment[a][B.PRMs.length - 1] > Score)
				Score = Alignment[a][B.PRMs.length - 1];
		}
		for (int b = 0; b < B.PRMs.length; ++b) {
			if (Alignment[A.PRMs.length - 1][b] > Score)
				Score = Alignment[A.PRMs.length - 1][b];
		}

		return Score;

	}

	public static double AlignmentScore(String A, PRMSpectrum B,
			double Tolerance) {
		double Score = 0.0;
		int[] APRMs = AntibodyUtils.GetSeqPRMScaled(A);
		double[][] Alignment = new double[APRMs.length][B.PRMs.length];

		for (int a = 0; a < APRMs.length; ++a) {
			// double AScore = A.PRMScores[a];
			int AMass = APRMs[a];
			for (int b = 0; b < B.PRMs.length; ++b) {
				double BScore = B.PRMScores[b];
				int BMass = B.PRMs[b];

				double MaxScore = BScore;
				for (int prevA = 0; prevA < a; ++prevA) {
					int prevAMass = APRMs[prevA];
					for (int prevB = 0; prevB < b; ++prevB) {
						int prevBMass = B.PRMs[prevB];

						if (AntibodyUtils.EquivalentScaled(BMass - prevBMass,
								AMass - prevAMass, Tolerance)) {
							double NewScore = BScore + Alignment[prevA][prevB];
							if (NewScore > MaxScore)
								MaxScore = NewScore;
							// System.out.println("Extending [" + prevA + "," +
							// prevB + "] -> [" + a + "," + b + "] = " +
							// NewScore);
						}
					}

				}
				Alignment[a][b] = MaxScore;
				if (MaxScore > Score)
					Score = MaxScore;
			}
		}
		return Score;

	}

	/*
	 * Attempt at the SPS scoring of spectrum alignments
	 */
	public static double alignmentScoreSpecAlign(String seq,
			PRMSpectrum spectrum, double tol) {

		boolean LocalDebug = false;
		int[] currSequencePRMs = AntibodyUtils.GetSeqPRMScaled(seq);

		if (LocalDebug) {
			System.out.println("alignmentScoreSpecAlign");
			System.out.println("Seq PRMs: "
					+ Utils.JoinIntArray(currSequencePRMs, " "));
			System.out.println("Spec PRMs: "
					+ Utils.JoinIntArray(spectrum.PRMs, " "));
			System.out.println("Spec PRMScores: "
					+ Utils.JoinDoubleArray(spectrum.PRMScores, " "));
		}

		// The sum of the PRM peaks scores overlapping
		double[][] oScoreA = new double[currSequencePRMs.length][spectrum.PRMs.length];
		double[][] oScoreB = new double[currSequencePRMs.length][spectrum.PRMs.length];

		// Precompute the sums of peaks in different ranges.
		for (int i = 0; i < currSequencePRMs.length; ++i) {
			for (int j = 0; j < spectrum.PRMs.length; ++j) {

				if (LocalDebug)
					System.out.println("Aligning Seq " + i + ": "
							+ currSequencePRMs[i] + " to Spec " + j + ": "
							+ spectrum.PRMs[j]);

				// Compute for the sequence first

				int startingIndex = Utils.getNextGreatestPeak(currSequencePRMs,
						currSequencePRMs[i]
								- (spectrum.PRMs[j] - spectrum.PRMs[0]));
				int endingIndex = Utils
						.getNextSmallestPeak(
								currSequencePRMs,
								currSequencePRMs[i]
										+ (spectrum.PRMs[spectrum.PRMs.length - 1] - spectrum.PRMs[j]));
				oScoreA[i][j] = endingIndex - startingIndex + 1;

				if (LocalDebug) {
					System.out
							.println("Looking for starting mass in seq prms: "
									+ (currSequencePRMs[i] - (spectrum.PRMs[j] - spectrum.PRMs[0])));
					System.out.println("Starting peak is at " + startingIndex
							+ ": " + currSequencePRMs[startingIndex]);
					System.out
							.println("Looking for ending mass in seq prms: "
									+ (currSequencePRMs[i] + (spectrum.PRMs[spectrum.PRMs.length - 1] - spectrum.PRMs[j])));
					System.out.println("Ending peak is at " + endingIndex
							+ ": " + currSequencePRMs[endingIndex]);
				}

				// Compute for the spectrum second
				startingIndex = Utils.getNextGreatestPeak(spectrum.PRMs,
						spectrum.PRMs[j] - (currSequencePRMs[i]));
				endingIndex = Utils
						.getNextSmallestPeak(
								spectrum.PRMs,
								spectrum.PRMs[j]
										+ (currSequencePRMs[currSequencePRMs.length - 1] - currSequencePRMs[i]));
				oScoreB[i][j] = 0;
				for (int k = startingIndex; k <= endingIndex; ++k)
					oScoreB[i][j] += spectrum.PRMScores[k];

				if (LocalDebug) {
					System.out
							.println("Looking for starting mass in spec prms: "
									+ (spectrum.PRMs[j] - (currSequencePRMs[i])));
					System.out.println("Starting peak is at " + startingIndex
							+ ": " + spectrum.PRMs[startingIndex]);
					System.out
							.println("Looking for ending mass in spec prms: "
									+ (spectrum.PRMs[j] + (currSequencePRMs[currSequencePRMs.length - 1] - currSequencePRMs[i])));
					System.out.println("Ending peak is at " + endingIndex
							+ ": " + spectrum.PRMs[endingIndex]);
					System.out.println("oScoreB[" + i + "][" + j + "]: "
							+ oScoreB[i][j]);
					// Utils.WaitForEnter();
				}
			}
		}

		if (LocalDebug) {
			System.out.println("Overlap matrices:");
			System.out.println("A: ");
			String sA = "";
			String sB = "";
			for (int i = 0; i < currSequencePRMs.length; ++i) {
				for (int j = 0; j < spectrum.PRMs.length; ++j) {
					sA += oScoreA[i][j] + " ";
					sB += oScoreB[i][j] + " ";
				}
				sA += "\n";
				sB += "\n";
			}
			System.out.println(sA);
			System.out.println("\nB: ");
			System.out.println(sB);
			Utils.WaitForEnter();
		}

		// matchScoreA[i][j] is the sum of the PRMs in seq that are matched in
		// the
		// best alignment ending with seq peak i matching spectrum peak j
		double[][] matchScoreA = new double[currSequencePRMs.length][spectrum.PRMs.length];
		double[][] matchScoreB = new double[currSequencePRMs.length][spectrum.PRMs.length];
		double[][] ratioScoresA = new double[currSequencePRMs.length][spectrum.PRMs.length];
		double[][] ratioScoresB = new double[currSequencePRMs.length][spectrum.PRMs.length];
		int[][] prev = new int[currSequencePRMs.length][spectrum.PRMs.length];

		for (int i = 0; i < currSequencePRMs.length; ++i) {

			matchScoreA[i][0] = 1;
			matchScoreB[i][0] = spectrum.PRMScores[0];
			ratioScoresA[i][0] = matchScoreA[i][0] / oScoreA[i][0];
			ratioScoresB[i][0] = matchScoreB[i][0] / oScoreB[i][0];
			prev[i][0] = -1;
		}
		for (int i = 0; i < spectrum.PRMs.length; ++i) {
			matchScoreA[0][i] = 1;
			matchScoreB[0][i] = spectrum.PRMScores[i];
			ratioScoresA[0][i] = matchScoreA[0][i] / oScoreA[0][i];
			ratioScoresB[0][i] = matchScoreB[0][i] / oScoreB[0][i];
			prev[0][i] = -1;
		}

		for (int i = 1; i < currSequencePRMs.length; ++i) {
			for (int j = 1; j < spectrum.PRMs.length; ++j) {

				if (LocalDebug)
					System.out.println("Considering seq PRM " + i + ": "
							+ currSequencePRMs[i] + " and spec PRM " + j + ": "
							+ spectrum.PRMs[j]);
				// Find the best prevScore
				double bestScore = -1;
				int[] bestPrev = { -1, -1 };
				for (int prevSeqIndex = Math.max(0, i - 2); prevSeqIndex < i; ++prevSeqIndex) {

					int seqDelta = currSequencePRMs[i]
							- currSequencePRMs[prevSeqIndex];

					for (int prevSpecIndex = j - 1; prevSpecIndex >= 0; --prevSpecIndex) {

						int specDelta = spectrum.PRMs[j]
								- spectrum.PRMs[prevSpecIndex];
						if (LocalDebug) {
							System.out.println("Consdering prev seq "
									+ prevSeqIndex + " and prev spec "
									+ prevSpecIndex);
							System.out.println("seqDelta: " + seqDelta);
							System.out.println("specDelta: " + specDelta);
						}
						if (Math.abs(specDelta - seqDelta) <= tol
								* AntibodyUtils.MASS_SCALE) {
							if (LocalDebug)
								System.out.println("FOUND A PREV: Seq PRM "
										+ prevSeqIndex + ": "
										+ currSequencePRMs[prevSeqIndex]
										+ " and spec PRM " + prevSpecIndex
										+ ": " + spectrum.PRMs[prevSpecIndex]);
							if (Math.min(
									ratioScoresA[prevSeqIndex][prevSpecIndex],
									ratioScoresB[prevSeqIndex][prevSpecIndex]) > bestScore) {
								bestScore = Math
										.min(ratioScoresA[prevSeqIndex][prevSpecIndex],
												ratioScoresB[prevSeqIndex][prevSpecIndex]);
								bestPrev[0] = prevSeqIndex;
								bestPrev[1] = prevSpecIndex;
								if (LocalDebug)
									System.out
											.println("Found a new best score: "
													+ bestScore);
							}

						}
					}
				}

				if (bestPrev[0] == -1) // We didn't find a reasonable previous,
										// so we just count this peak as aligned
				{
					if (LocalDebug)
						System.out.println("No prev found for this matching");
					matchScoreA[i][j] = 1;
					matchScoreB[i][j] = spectrum.PRMScores[j];
					ratioScoresA[i][j] = matchScoreA[i][j] / oScoreA[i][j];
					ratioScoresB[i][j] = matchScoreB[i][j] / oScoreB[i][j];
					prev[i][j] = -1;
				}
				// We found a good prev!
				else {
					if (LocalDebug)
						System.out.println("A prev was found!");
					matchScoreA[i][j] = 1 + matchScoreA[bestPrev[0]][bestPrev[1]];
					matchScoreB[i][j] = spectrum.PRMScores[j]
							+ matchScoreB[bestPrev[0]][bestPrev[1]];
					ratioScoresA[i][j] = matchScoreA[i][j] / oScoreA[i][j];
					ratioScoresB[i][j] = matchScoreB[i][j] / oScoreB[i][j];
				}
				// if(LocalDebug)
				// Utils.WaitForEnter();
				prev[i][j] = bestPrev[0] * spectrum.PRMs.length + bestPrev[1];
			}
		}
		if (LocalDebug) {
			for (int i = 0; i < currSequencePRMs.length; ++i) {
				String currLine = "";
				for (int j = 0; j < spectrum.PRMs.length; ++j) {
					currLine += Math
							.min(ratioScoresA[i][j], ratioScoresB[i][j]) + " ";
				}
				System.out.println(currLine);
			}
			Utils.WaitForEnter();
		}

		double bestScore = -1;
		int[] bestIndx = { -1, -1 };
		for (int i = MIN_OVERLAPPING_PEAKS - 1; i < currSequencePRMs.length; ++i) {
			for (int j = MIN_OVERLAPPING_PEAKS - 1; j < spectrum.PRMs.length; ++j) {
				if (Math.min(ratioScoresA[i][j], ratioScoresB[i][j]) >= bestScore) {
					bestScore = Math
							.min(ratioScoresA[i][j], ratioScoresB[i][j]);
					bestIndx[0] = i;
					bestIndx[1] = j;
				}
			}
		}

		if (LocalDebug)
			System.out.println("Best align is with " + bestIndx[0] + " to "
					+ bestIndx[1]);
		String matchStr = "";
		String prmStr = "";
		int length = 0;
		while (bestIndx[0] != -1) {
			matchStr = "(" + bestIndx[0] + "," + bestIndx[1] + ") " + matchStr;
			prmStr = "(" + currSequencePRMs[bestIndx[0]] + ","
					+ spectrum.PRMs[bestIndx[1]] + ") " + prmStr;
			length += 1;

			int prevVal = prev[bestIndx[0]][bestIndx[1]];
			if (prevVal == -1)
				break;

			bestIndx[0] = prevVal / spectrum.PRMs.length;
			bestIndx[1] = prevVal % spectrum.PRMs.length;
			if (LocalDebug)
				System.out.println("prev: " + bestIndx[0] + "," + bestIndx[1]);
		}

		if (length < AntibodyUtils.MIN_OVERLAPPING_PEAKS) {
			if (LocalDebug)
				System.out.println("Alignment is short: " + length + " < "
						+ (AntibodyUtils.MIN_OVERLAPPING_PEAKS));
			return -1;
		}
		return bestScore;
	}

	/*
	 * public static double AlignmentScoreSpecAlign_old(String Sequence,
	 * PRMSpectrum Spectrum, double Tolerance) { boolean LocalDebug = false;
	 * int[] CurrSequencePRMs = AntibodyUtils.GetSeqPRMScaled(Sequence); //
	 * Inititalize the alignment matrix double[][] H = new
	 * double[CurrSequencePRMs.length][Spectrum.PRMs.length]; int[][] PrevH =
	 * new int[CurrSequencePRMs.length][Spectrum.PRMs.length];
	 * 
	 * // MassShift represents how much we have to add to the SpectrumPRM to //
	 * get to the Sequence PRM int[][] MassShift = new
	 * int[CurrSequencePRMs.length][Spectrum.PRMs.length];
	 * 
	 * H[0][0] = Spectrum.PRMScores[0]; PrevH[0][0] = -1;
	 * 
	 * for (int i = 0; i < CurrSequencePRMs.length; ++i) { H[i][0] = H[0][0];
	 * PrevH[i][0] = -1; MassShift[i][0] = CurrSequencePRMs[i]; } for (int i =
	 * 0; i < Spectrum.PRMs.length; ++i) { H[0][i] = Spectrum.PRMScores[i];
	 * 
	 * PrevH[0][i] = -1; MassShift[0][i] = -1 * Spectrum.PRMs[i]; //
	 * MassShift[0][i] = Spectrum.PRMs[i]; }
	 * 
	 * double BestValue = Double.MIN_VALUE; int BestCell = -1;
	 * 
	 * // Fill in the table for (int SeqIndex = 1; SeqIndex <
	 * CurrSequencePRMs.length; ++SeqIndex) { int CurrSeqMass =
	 * CurrSequencePRMs[SeqIndex]; for (int SpecIndex = 1; SpecIndex <
	 * Spectrum.PRMs.length; ++SpecIndex) { int CurrSpecMass =
	 * Spectrum.PRMs[SpecIndex];
	 * 
	 * double Value = Spectrum.PRMScores[SpecIndex]; double MaxValue = Value;
	 * 
	 * int CurrAlignmentShift = CurrSeqMass - CurrSpecMass; int
	 * MaxAlignmentShift = CurrAlignmentShift;
	 * 
	 * int PrevCell = -1; // Index of previous match int MaxPrevCell = PrevCell;
	 * 
	 * // Look back no more more than 1 missed AA for (int PrevSeqIndex =
	 * SeqIndex - 1; PrevSeqIndex > Math.max( SeqIndex - 4, -1); PrevSeqIndex--)
	 * { int PrevSeqMass = CurrSequencePRMs[PrevSeqIndex];
	 * 
	 * for (int PrevSpecIndex = SpecIndex - 1; PrevSpecIndex >= 0;
	 * PrevSpecIndex--) { int PrevSpecMass = Spectrum.PRMs[PrevSpecIndex];
	 * 
	 * // Check that the delta in Spectrum mass is about the // same as teh
	 * delta in Sequence masses // Check that the CurrentAlignment Shift is
	 * about the // same as the previous alignment shift if
	 * (Math.abs((CurrSeqMass - PrevSeqMass) - (CurrSpecMass - PrevSpecMass)) <
	 * (int) (Tolerance * AntibodyUtils.MASS_SCALE) && Math
	 * .abs(CurrAlignmentShift - MassShift[PrevSeqIndex][PrevSpecIndex]) < (int)
	 * (Tolerance * AntibodyUtils.MASS_SCALE)) { if
	 * (H[PrevSeqIndex][PrevSpecIndex] > 0) Value =
	 * H[PrevSeqIndex][PrevSpecIndex] + Spectrum.PRMScores[SpecIndex];
	 * 
	 * // if(LocalDebug) // { //
	 * System.out.println("Considering Previous SeqPRM[" // + PrevSeqIndex +
	 * "]="+ PrevSeqMass + // " to SpecPRM[" + PrevSpecIndex + "]=" + //
	 * PrevSpecMass); // System.out.println("Successful extension!!"); // } } if
	 * (Value > MaxValue) { MaxValue = Value; MaxPrevCell = PrevSeqIndex *
	 * Spectrum.PRMs.length + PrevSpecIndex; MaxAlignmentShift =
	 * MassShift[PrevSeqIndex][PrevSpecIndex]; } }
	 * 
	 * }
	 * 
	 * H[SeqIndex][SpecIndex] = MaxValue; PrevH[SeqIndex][SpecIndex] =
	 * MaxPrevCell; MassShift[SeqIndex][SpecIndex] = MaxAlignmentShift;
	 * 
	 * if (MaxValue > BestValue) { BestValue = MaxValue; BestCell = SeqIndex *
	 * Spectrum.PRMs.length + SpecIndex; if (LocalDebug) {
	 * 
	 * System.out.println("New Max at [" + SeqIndex + "][" + SpecIndex + "]=" +
	 * BestValue); System.out.println("Prev [" + MaxPrevCell /
	 * Spectrum.PRMs.length + "," + MaxPrevCell % Spectrum.PRMs.length + "]");
	 * System.out.println("MassShift: " + MaxAlignmentShift);
	 * Utils.WaitForEnter(); } }
	 * 
	 * }
	 * 
	 * } int Previ = 0;
	 * 
	 * int BestSeqIndex = BestCell / Spectrum.PRMs.length; int BestSpecIndex =
	 * BestCell % Spectrum.PRMs.length;
	 * 
	 * if (BestSeqIndex < CurrSequencePRMs.length - 3) { if (LocalDebug)
	 * System.out.println("Alignment doesn't end near the end of Seq " +
	 * BestSeqIndex + " < " + (CurrSequencePRMs.length - 3)); return -1; }
	 * ArrayList SeqMasses = new ArrayList(); ArrayList SpecMasses = new
	 * ArrayList();
	 * 
	 * int TempSeqIndex = BestSeqIndex; int TempSpecIndex = BestSpecIndex; int
	 * Temp = BestCell;
	 * 
	 * while (Temp >= 0) { TempSeqIndex = Temp / Spectrum.PRMs.length;
	 * TempSpecIndex = Temp % Spectrum.PRMs.length; SeqMasses.add(0, new
	 * Integer(CurrSequencePRMs[TempSeqIndex])); SpecMasses.add(0, new
	 * Integer(Spectrum.PRMs[TempSpecIndex])); Previ = TempSeqIndex; Temp =
	 * PrevH[TempSeqIndex][TempSpecIndex];
	 * 
	 * }
	 * 
	 * // String[] SequenceSegments = AntibodyUtils.SplitStringtoAA(Sequence);
	 * // String Prefix = ""; // for(int i = Previ; i < BestSeqIndex; ++i) // {
	 * // Prefix += SequenceSegments[i]; // }
	 * 
	 * // String Prefix = this.Sequence.substring(Previ,BestSeqIndex); if
	 * (LocalDebug) { System.out.println("AlignedPeaks: " +
	 * AntibodyUtils.IntegerListToString(AntibodyUtils
	 * .ConvertIntegerArrayList(SeqMasses)));
	 * System.out.println("SeqMasses.size() = " + SeqMasses.size()); }
	 * 
	 * // Check for a valid prefix or suffix mass to match int[] SpecMassesFinal
	 * = AntibodyUtils .ConvertIntegerArrayList(SpecMasses); int[]
	 * SeqMassesFinal = AntibodyUtils.ConvertIntegerArrayList(SeqMasses);
	 * 
	 * boolean HasPrefix = SpecAlign.IsValidPrefix(SpecMassesFinal[0],
	 * SeqMassesFinal[0], CurrSequencePRMs, Tolerance);
	 * 
	 * int SpectrumPrefixMass = SeqMassesFinal[0] - SpecMassesFinal[0]; // int
	 * SpectrumSuffixMass = //
	 * this.GetSuffixMass(SpecMassesFinal[SpecMassesFinal.length-1], //
	 * SeqMassesFinal
	 * [SeqMassesFinal.length-1],Spectrum.PrecursorMass,SequenceMass);
	 * 
	 * if (HasPrefix && SpectrumPrefixMass >= -1 * Tolerance
	 * AntibodyUtils.MASS_SCALE) { if (LocalDebug)
	 * System.out.println("HAS PREFIX!"); BestValue += 15; }
	 * 
	 * if (SpectrumPrefixMass < -1 * Tolerance * AntibodyUtils.MASS_SCALE) { if
	 * (SeqMasses.size() < AntibodyUtils.MIN_OVERLAPPING_PEAKS + 1) { if
	 * (LocalDebug) System.out
	 * .println("Spectrum totally overlaps anchor, but alignment is short: " +
	 * SeqMasses.size() + " < " + (AntibodyUtils.MIN_OVERLAPPING_PEAKS + 1));
	 * return -1; } } if (SeqMasses.size() <
	 * AntibodyUtils.MIN_OVERLAPPING_PEAKS) { if (LocalDebug)
	 * System.out.println("Alignment is short: " + SeqMasses.size() + " < " +
	 * (AntibodyUtils.MIN_OVERLAPPING_PEAKS)); return -1; } return BestValue;
	 * 
	 * }
	 */
	public static double AlignmentScoreLength(PRMSpectrum A, PRMSpectrum B,
			double Tolerance) {
		double Score = 0.0;

		double[][] Alignment = new double[A.PRMs.length][B.PRMs.length];
		for (int a = 0; a < A.PRMs.length; ++a) {
			double AScore = A.PRMScores[a];
			int AMass = A.PRMs[a];
			for (int b = 0; b < B.PRMs.length; ++b) {
				double BScore = B.PRMScores[b];
				int BMass = B.PRMs[b];

				double MaxScore = 1;
				for (int prevA = 0; prevA < a; ++prevA) {
					int prevAMass = A.PRMs[prevA];
					for (int prevB = 0; prevB < b; ++prevB) {
						int prevBMass = B.PRMs[prevB];

						if (AntibodyUtils.EquivalentScaled(BMass - prevBMass,
								AMass - prevAMass, Tolerance)) {
							double NewScore = 1 + Alignment[prevA][prevB];
							if (NewScore > MaxScore)
								MaxScore = NewScore;

						}
					}

				}
				Alignment[a][b] = MaxScore;
			}
		}

		for (int a = 0; a < A.PRMs.length; ++a) {
			if (Alignment[a][B.PRMs.length - 1] > Score)
				Score = Alignment[a][B.PRMs.length - 1];
		}
		for (int b = 0; b < B.PRMs.length; ++b) {
			if (Alignment[A.PRMs.length - 1][b] > Score)
				Score = Alignment[A.PRMs.length - 1][b];
		}

		return Score;

	}

	public static int[] SortIntArrayDecreasing(int[] OldArray) {
		int[] Ret = new int[OldArray.length];
		for (int i = 0; i < OldArray.length - 1; ++i) {
			int CurrMax = OldArray[i];
			for (int j = i + 1; j < OldArray.length; ++j) {
				if (OldArray[j] > CurrMax) {
					CurrMax = OldArray[j];
					int temp = OldArray[i];
					OldArray[i] = OldArray[j];
					OldArray[j] = temp;
				}
			}

			Ret[i] = CurrMax;
		}
		return Ret;
	}

	public static int[] RankIntArrayDecreasing(int[] OldArray) {

		if (OldArray.length == 1) {
			int[] Ret = { 0 };
			return Ret;
		}
		boolean LocalDebug = false;

		if (LocalDebug)
			System.out.println("Ranking: " + Utils.IntArrayToString(OldArray));
		int[] Ret = new int[OldArray.length];
		for (int i = 0; i < Ret.length; ++i)
			Ret[i] = -1;
		for (int i = 0; i < OldArray.length; ++i) {
			int Count = 0;
			if (LocalDebug)
				System.out.println("Considering OldArray[" + i + "]="
						+ OldArray[i]);
			for (int j = 0; j < OldArray.length; ++j) {
				if (i == j)
					continue;
				if (OldArray[j] > OldArray[i])
					Count++;
				if (LocalDebug) {
					System.out.println("   Comparing to OldArray[" + j + "]="
							+ OldArray[j]);
					System.out.println("Count: " + Count);
				}
			}

			while (Ret[Count] != -1)
				Count++;
			Ret[Count] = i;
			if (LocalDebug)
				System.out.println("Final resting spot: " + Count);
		}

		if (LocalDebug) {
			for (int i = 0; i < Ret.length; ++i) {
				System.out.println("Rank " + i + " = " + OldArray[Ret[i]]);
			}
			Utils.WaitForEnter();
		}
		return Ret;
	}

	/**
	 * Creates a unique key for this spetrum file and scan number File_Scan,
	 * where File contains only the base name without extension
	 * 
	 * @param FileName
	 * @param ScanNumber
	 * @return
	 */
	public static String CreateSpectrumKey(String FileName, int[] ScanNumber) {
		String key = Utils.GetBaseNameNoExtension(FileName);
		for (int i = 0; i < ScanNumber.length; ++i)
			key += "_" + ScanNumber[i];
		return key;
	}

	/**
	 * Returns a 'delim' delimitted list of the hashtable keys (if they are
	 * Strings)
	 * 
	 * @param Hash
	 * @param delim
	 * @return
	 */
	public static String GetHashtableStringKeys(Hashtable Hash, String delim) {
		String Ret = "";
		Enumeration E = Hash.keys();

		while (E.hasMoreElements()) {
			Ret += delim + (String) (E.nextElement());
		}

		return Ret.substring(1);

	}

	public static String DetermineDigestTypePepNovo(String FileName) {
		if (FileName.toLowerCase().indexOf("tryp") >= 0
				&& FileName.toLowerCase().indexOf("chymo") < 0)
			return "TRYPSIN";

		return "NON_SPECIFIC";

	}

	public static String DetermineDigestTypeInspect(String FileName) {
		if (FileName.toLowerCase().indexOf("tryp") >= 0
				&& FileName.toLowerCase().indexOf("chymo") < 0)
			return "Trypsin";
		else if (FileName.toLowerCase().indexOf("tryp") >= 0
				&& FileName.toLowerCase().indexOf("chymo") >= 0)
			return "Chymotrypsin";

		return "None";

	}

	/**
	 * Looks up a file into a list of files, comparing only the Base names with
	 * no extension
	 * 
	 * @param QueryFileName
	 *            The query file name
	 * @param SubjectFiles
	 *            A list of file names
	 * @return Returns the version of the file contained in the file list
	 */
	public static String GetStringMatchNoExtension(String QueryFileName,
			String[] SubjectFiles) {

		if (SubjectFiles == null || SubjectFiles.length == 0
				|| QueryFileName == null) {
			System.err
					.println("ERROR: AntibodyUtils.GetStringMatchNoExtension requires non-null parameters");
			return null;
		}
		String Query = Utils.GetBaseNameNoExtension(QueryFileName);
		for (int i = 0; SubjectFiles != null && SubjectFiles[i] != null
				&& i < SubjectFiles.length; ++i) {
			String Subject = Utils.GetBaseNameNoExtension(SubjectFiles[i]);
			if (Subject.compareTo(Query) == 0)
				return SubjectFiles[i];
		}
		return null;
	}

	/**
	 * Loads a file into a string array
	 * 
	 * @param FileName
	 * @return
	 */
	public static String[] LoadEachFileLine(String FileName) {

		BufferedReader Buf = null;
		String Line = "";
		try {
			Buf = new BufferedReader(new FileReader(FileName));
			Line = Buf.readLine();
		} catch (Exception E) {
			E.printStackTrace();
			return null;
		}

		ArrayList Lines = new ArrayList();
		while (Line != null) {
			Line = Line.trim();
			if (Line.length() > 0)
				Lines.add(Line);

			try {
				Line = Buf.readLine();
			} catch (Exception E) {
				E.printStackTrace();
				return null;
			}

		}

		try {
			Buf.close();
		} catch (Exception E) {
			E.printStackTrace();
			return null;
		}
		String[] Ret = new String[Lines.size()];
		for (int i = 0; i < Ret.length; ++i)
			Ret[i] = (String) (Lines.get(i));

		return Ret;

	}

	/**
	 * Given a PRMList and a minimum mass, returns the PRMList with only masses
	 * >= MinMass, and shifted so that the first PRM is 0.
	 * 
	 * @param OrigList
	 *            * @return
	 */
	public static Object[] ReturnAllFollowing(int[] OrigList, int MinMass,
			double[] OrigPRMScores, double Tolerance) {
		boolean Debug = false;
		int MassShift = Integer.MIN_VALUE;
		ArrayList Temp = new ArrayList();
		ArrayList ScoreTemp = new ArrayList();
		for (int i = 0; i < OrigList.length; ++i) {
			if (OrigList[i] >= MinMass - Tolerance * AntibodyUtils.MASS_SCALE) {

				if (MassShift == Integer.MIN_VALUE)
					MassShift = OrigList[i];
				if (Debug) {
					System.out.println("PRM[" + i + "] = " + OrigList[i] + "->"
							+ (OrigList[i] - MassShift));
					System.out.println("Score[" + i + "]=" + OrigPRMScores[i]);
				}
				Temp.add(new Integer(OrigList[i] - MassShift));
				ScoreTemp.add(new Double(OrigPRMScores[i]));
			}

		}
		Object[] Ret = new Object[2];
		Ret[0] = AntibodyUtils.ConvertIntegerArrayList(Temp);
		Ret[1] = AntibodyUtils.ConvertDoubleArrayList(ScoreTemp);
		if (Debug)
			System.out.println("PRMLength: " + Temp.size());
		return Ret;
	}

	/**
	 * Given a PRMList and a minimum mass, returns the PRMList with only masses
	 * <= MinMass, and shifted so that the first PRM is 0.
	 * 
	 * @param OrigList
	 * @param MinMass
	 * @param OrigPRMScores
	 * @param Tolerance
	 * @return
	 */
	public static Object[] ReturnAllPreceeding(int[] OrigList, int MinMass,
			double[] OrigPRMScores, double Tolerance) {
		boolean Debug = false;
		int MassShift = -1;
		ArrayList Temp = new ArrayList();
		ArrayList ScoreTemp = new ArrayList();
		for (int i = 0; i < OrigList.length; ++i) {
			if (OrigList[i] <= MinMass + Tolerance * AntibodyUtils.MASS_SCALE) {

				if (i == 0)
					MassShift = OrigList[i];
				if (Debug) {
					System.out.println("PRM[" + i + "] = " + OrigList[i] + "->"
							+ (OrigList[i] - MassShift));
					System.out.println("Score[" + i + "]=" + OrigPRMScores[i]);
				}
				Temp.add(new Integer(OrigList[i] - MassShift));
				ScoreTemp.add(new Double(OrigPRMScores[i]));
			}

		}
		Object[] Ret = new Object[2];
		Ret[0] = AntibodyUtils.ConvertIntegerArrayList(Temp);
		Ret[1] = AntibodyUtils.ConvertDoubleArrayList(ScoreTemp);
		if (Debug)
			System.out.println("PRMLength: " + Temp.size());
		return Ret;
	}

	/**
	 * Sorts a list of peaks. Each peak in the arraylist is a double array,
	 * first element is the mass, second element is the intensity
	 * 
	 * @param Peaks
	 * @return
	 */
	public static ArrayList SortPeaks(ArrayList Peaks) {
		for (int i = 0; i < Peaks.size(); ++i) {
			int MinMass = (int) ((double[]) (Peaks.get(i)))[0];
			int MinIndex = i;
			for (int j = i + 1; j < Peaks.size(); ++j) {
				int CurrMass = (int) ((double[]) (Peaks.get(j)))[0];
				if (CurrMass < MinMass) {
					MinMass = CurrMass;
					MinIndex = j;
				}
			}
			if (MinIndex != i) {
				double[] Temp = (double[]) (Peaks.get(MinIndex));
				Peaks.set(MinIndex, (double[]) (Peaks.get(i)));
				Peaks.set(i, Temp);
			}

		}
		return Peaks;

	}

	public static int[] GetStrongestPath(int[][] AdjacencyScores) {
		boolean Debug = true;
		ArrayList[] BestPaths = null;
		int BestScore = 0;
		int BestStart = 0;
		for (int StartSeed = 0; StartSeed < AdjacencyScores.length; ++StartSeed) {
			int[] scores = new int[AdjacencyScores.length];
			ArrayList[] paths = new ArrayList[scores.length];
			ArrayList Temp;
			for (int i = 0; i < scores.length; ++i) {
				scores[i] = 0;
				paths[i] = new ArrayList();
				Temp = new ArrayList();
				Temp.add(new Integer(i));
				paths[i].add(Temp);
			}

			for (int CurrSeed = StartSeed + 1; CurrSeed < scores.length; ++CurrSeed) {
				for (int PrevSeed = StartSeed; PrevSeed < CurrSeed; ++PrevSeed) {
					if (AdjacencyScores[PrevSeed][CurrSeed] > 0) {
						if (scores[PrevSeed]
								+ AdjacencyScores[PrevSeed][CurrSeed] > scores[CurrSeed]) {
							scores[CurrSeed] = scores[PrevSeed]
									+ AdjacencyScores[PrevSeed][CurrSeed];
							paths[CurrSeed] = new ArrayList();
							// distances[CurrNode] = distances[PrevNode] + 1;
							for (int i = 0; i < paths[PrevSeed].size(); ++i) {
								Temp = AntibodyUtils
										.CopyIntArrayList((ArrayList) (paths[PrevSeed]
												.get(i)));
								Temp.add(new Integer(CurrSeed));
								paths[CurrSeed].add(Temp);
							}
						}

					}

				}

			}

			/**
			 * Once we have computed the longest paths to each node, determine
			 * the score of the strongest paths overall. maxNode is only one of
			 * several nodes that could have a strongest path from node 0.
			 */
			int maxScore = Integer.MIN_VALUE;
			int maxNode = 0;
			// int maxDist = 0;
			for (int i = 0; i < scores.length; ++i) {
				if (scores[i] >= maxScore) {
					// maxDist = distances[i];
					maxScore = scores[i];
					maxNode = i;
				}
			}

			if (Debug) {
				System.out.println("Distance to " + (scores.length - 1) + " = "
						+ scores[scores.length - 1]);

				System.out.println("Strongest path starting at " + StartSeed
						+ " ends at node " + maxNode + " at score " + maxScore);

			}
			if (maxScore > BestScore) {

				BestScore = maxScore;
				BestStart = StartSeed;
				BestPaths = new ArrayList[paths[maxNode].size()];
				/**
				 * Determine the final paths, by converting the path[maxNode].
				 * Keep in mind that this is only 1 of several possible longest
				 * paths from node 0 to maxNode, and only 1 of several nodes
				 * which may have a path from node 0 of maxDist.
				 */
				for (int p = 0; p < paths[maxNode].size(); ++p) {
					BestPaths[p] = (ArrayList) (paths[maxNode].get(p));

					// finalPaths[p].add(new Integer(maxNode));
				}
			}

		}

		if (Debug)
			System.out.println("Best scoring = " + BestScore + " starts at "
					+ BestStart);

		int Longest = 0;
		int LongestIndex = 0;
		for (int i = 0; i < BestPaths.length; ++i) {
			if (BestPaths[i].size() >= Longest) {
				Longest = BestPaths[i].size();
				LongestIndex = i;
			}
		}
		if (Debug)
			System.out.println("Longest Path: " + Longest);

		return AntibodyUtils.ConvertIntegerArrayList(BestPaths[LongestIndex]);

	}

	/**
	 * Converst a string array list into a string array
	 * 
	 * @param Input
	 * @return
	 */
	public static String[] ConvertStringArrayList(ArrayList Input) {
		String[] Ret = new String[Input.size()];
		for (int i = 0; i < Ret.length; ++i)
			Ret[i] = (String) (Input.get(i));
		return Ret;
	}

	/**
	 * Sorts a list of Comparables
	 * 
	 * @param Input
	 * @return
	 */
	public static ArrayList SortArrayList(ArrayList Input) {
		for (int i = 0; i < Input.size() - 1; ++i) {
			Comparable CurrMin = (Comparable) (Input.get(i));
			int MinIndex = i;
			for (int j = i + 1; j < Input.size(); ++j) {
				Comparable CurrSeed = (Comparable) (Input.get(j));
				if (CurrMin.compareTo((Object) CurrSeed) > 0) {
					CurrMin = CurrSeed;
					MinIndex = j;
				}
			}
			if (MinIndex != i) {
				Object Temp = Input.get(MinIndex);
				Input.set(MinIndex, Input.get(i));
				Input.set(i, Temp);
			}
		}
		return Input;
	}

	public static void main(String[] args) throws IOException {

		// Check file extension stuff
		String test = "./Examples/Test.txt.out";
		String testWin = ".\\Examples\\Test.txt";
		String test2 = Utils.GetFileNameNoExtension(test);
		String test2Win = Utils.GetFileNameNoExtension(testWin);
		System.out.println(test);
		System.out.println(test2);
		System.out.println(Utils.GetBaseName(test));
		System.out.println(Utils.GetFileExtension(test));
		System.out.println(Utils.GetFilePath(test));
		System.out.println(testWin);
		System.out.println(test2Win);
		System.out.println(Utils.GetBaseName(testWin));
		System.out.println(Utils.GetFileExtension(testWin));
		System.out.println(Utils.GetFilePath(testWin));

		// Check path reconstruction
		/*
		 * String A =
		 * "IKESGPGIVAPSQSISITCTVSGFSITSYGVSWVRQPPGKGIEWIGVIWGDGSTNYHSAIISRISISKDNSKSQVFIKINSIQTDDTATYYCAKG[220.05][174.424]E-FYYAMDYWGQGTSVTVS[157.941]K-IFIVAKTTPPSVYPIAPGSAAQTNSMVTIGCIVKGYFPEPVTVTWNSGSISSGVHTFPAVIQSDIYTISSSVTVPSSTWPSQTVTCNVAHPASSTKVDKKIV[253.234]DC"
		 * +
		 * "-DCGCKPCICTV[226.12]V[273.135]FI[342.126][352.361]-SSSIVPEVSSVFIFPPKPKDVITITITPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTKPREEQF-TFRSVSEIPIMHQDWINGKEFKCRVNSAAFPAPIEKTISKTK[213.066][225.191]AP[227.024]-[351.217][224.965]APKVYTIPPPKEQMAKDKVSITCMITNFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKINVQKSNWEAGNTFTCSVIHEGIHNHHTEKSISHSPGK"
		 * ;
		 * 
		 * String Path = AntibodyUtils.ReconstructPath(A);
		 * System.out.println(Path); int Res =
		 * AntibodyUtils.GetFullSequenceCorrectness(Path);
		 * System.out.println("Heavy Chain Genomic DB:");
		 * System.out.println("  Correct: " + Res); Path = Path.replace("-","");
		 * System.out.println("  Length: " +
		 * AntibodyUtils.SplitStringtoAA(Path).length);
		 * AntibodyUtils.WaitForEnter();
		 * 
		 * A =
		 * "IKESGPGIVAPSQSISITCTVSGFSITSYGVSWVRQPPGKGIEWIGVIWGDGSTNYHSAIISRISISKDNSKSQVFIKINSIQTDDTATYYCAKG[220.05][174.424]E-"
		 * + "FYYAMDYWGQGTSVTVS[157.941]K-" +
		 * "DYWGQGTWTWSAKTTPPSVYPIAPGSAAQTNSMVTIGCIVKGYFPEPVTVTWNSGSISSGVHTFPAVIQSDIYTISSSVTVP[275.894]-"
		 * +
		 * "WPSETVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVITITITPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTQPREEQFNSTFRSVSEIPIMHQDWINGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTIPPPKEQMAKDKVSITCMITDFFPEDITVEWQWNGQPAENYKNTQPIMNTNGSYFVYSKINVQKSNWEAGNTFTCSVIHEGIHNHHTEKSISHSPGK"
		 * ; Path = AntibodyUtils.ReconstructPath(A); System.out.println(Path);
		 * Res = AntibodyUtils.GetFullSequenceCorrectness(Path);
		 * System.out.println("Heavy Chain Protein DB:");
		 * System.out.println("  Correct: " + Res); Path = Path.replace("-","");
		 * System.out.println("  Length: " +
		 * AntibodyUtils.SplitStringtoAA(Path).length);
		 * AntibodyUtils.WaitForEnter();
		 * 
		 * A = "NNVMSQS[184.636]SIAVSAGEKVTMSCKSSQSIINSR-" +
		 * "ENYIAWYQKKPGQSPKIIIYWASTRESGVPDRFTGSGSGTDFTITISSVQAEDIAVYYCK-" +
		 * "ITFGAGTKIEIKSA[158.077][314.211]-" +
		 * "RADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
		 * ; Path = AntibodyUtils.ReconstructPath(A); System.out.println(Path);
		 * Res = AntibodyUtils.GetFullSequenceCorrectnessLC(Path);
		 * System.out.println("Light Chain Protein DB:");
		 * System.out.println("  Correct: " + Res); Path = Path.replace("-","");
		 * System.out.println("  Length: " +
		 * AntibodyUtils.SplitStringtoAA(Path).length);
		 * AntibodyUtils.WaitForEnter();
		 */
		/*
		 * AntibodyUtils.AAMasses[2] += 57.00; AntibodyUtils.AAMassesScaled[2]
		 * += 57000;
		 * 
		 * String A =
		 * "TFRSVSEIPIMHQDWINGKEFKCRVNSAAFPAPIEKTISKTK[213.066][225.191]AP[227.024]"
		 * ; String B =
		 * "[351.226][224.955]APKVYTIPPPKEQMAKDKVSITCMITNFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKINVQKSNWEAGNTFTCSVIHEGIHNHHTEKSISHSPGK"
		 * ; int[] PRMA = AntibodyUtils.GetSeqPRMScaled(A); double [] ScoresA =
		 * new double[PRMA.length]; Object[] OverlappingSpectraA = new
		 * Object[PRMA.length]; Object[] SupportingSpectraA = new
		 * Object[PRMA.length];
		 * 
		 * for(int i = 0; i < ScoresA.length; ++i) { String[] Temp = new
		 * String[1]; Temp[0] = "Test"; OverlappingSpectraA[i] = Temp;
		 * SupportingSpectraA[i] = Temp; }
		 * 
		 * int[] PRMB = AntibodyUtils.GetSeqPRMScaled(B); double [] ScoresB =
		 * new double[PRMA.length]; Object[] OverlappingSpectraB = new
		 * Object[PRMA.length]; Object[] SupportingSpectraB = new
		 * Object[PRMA.length]; for(int i = 0; i < ScoresB.length; ++i) {
		 * String[] Temp = new String[1]; Temp[0] = "Test";
		 * OverlappingSpectraB[i] = Temp; SupportingSpectraB[i] = Temp; }
		 * 
		 * ConsensusPRMSpectrum S = ConsensusPRMSpectrum.MergeAnchors(PRMA,
		 * ScoresA, OverlappingSpectraA, SupportingSpectraA, PRMB, ScoresB,
		 * OverlappingSpectraB, SupportingSpectraB,
		 * AntibodyUtils.DefaultFragTolerance);
		 *//*
			 * FileWriter Writer = new FileWriter(
			 * "/home/natalie/Projects/MySVNProjects/ImmunoSeq/AntibodyData/aBTLA/Reconstructions/DIVERGENCE_SUMMARY_10162009.txt"
			 * ); for (int Cov = 65; Cov <= 95; Cov += 5) { for (int i = 0; i <
			 * 20; ++i) { String DivergenceFile =
			 * "/home/natalie/Projects/MySVNProjects/ImmunoSeq/AntibodyData/aBTLA/Reconstructions/aBTLA_HC_DIVERGENCE_10072009_Cov"
			 * + Cov + "_" + i + ".txt"; String ConcatedAnchors = "";
			 * System.out.println("Reading " + DivergenceFile); BufferedReader
			 * buf = new BufferedReader(new FileReader( DivergenceFile)); String
			 * Line = buf.readLine(); // Skip all the header stuff while (Line
			 * != null && Line.indexOf("# ----SEED INFO----") < 0) Line =
			 * buf.readLine();
			 * 
			 * while (Line != null) { if (Line.indexOf("Seed[") < 0) { Line =
			 * buf.readLine(); continue; } Line = Line.trim();
			 * System.out.println("Line: " + Line); int Index =
			 * Line.indexOf(">"); Line = Line.substring(Index + 1).trim();
			 * System.out.println("Got seed: " + Line); ConcatedAnchors += "-" +
			 * Line; Line = buf.readLine(); } buf.close();
			 * System.out.println("All: " + ConcatedAnchors.substring(1));
			 * String Path = AntibodyUtils.ReconstructPath(ConcatedAnchors
			 * .substring(1)); int Res =
			 * AntibodyUtils.GetFullSequenceCorrectness(Path);
			 * System.out.println("Heavy Chain Protein DB:");
			 * System.out.println("  Correct: " + Res); Path = Path.replace("-",
			 * ""); int L = AntibodyUtils.SplitStringtoAA(Path).length;
			 * System.out.println("  Length: " + L); Writer.write(DivergenceFile
			 * + "\t" + Cov + "\t" + i + "\t" + Res + "\t" + L + "\t" +
			 * ((double) (Res) / L) + "\t" + ConcatedAnchors + "\t" + Path +
			 * "\n"); Writer.flush(); // AntibodyUtils.WaitForEnter(); } }
			 */

		/*
		 * String HC_DB =
		 * "IKESGPGIVAPSQSISITCTVSGFSITSYGVSWVRQPPGKGIEWIGVIWGDGSTNYHSAIISRISISKDNSKSQVFIKINSIQTDDTATYYCAKG[220.05][174.424]E-FYYAMDYWGQGTSVTVSSAKTTPPSVYPIAPGSAAQTNSMVTIGCIVKGYFPEPVTVTWNSGSISSGVHTFPAVIQSDIYTISSSVTVP[275.894]-[188.169]WPSETVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVITITITPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTQPREEQFNSTFRSVSEIPIMHQDWINGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTIPPPKEQMAKDKVSITCMITDFFPEDITVEWQWNGQPAENYKNTQPIMNTNGSYFVYSKINVQKSNWEAGNTFTCSVIHEGIHNHHTEKSISHSPGK"
		 * ; String HC_Genomic =
		 * "IKESGPGIVAPSQSISITCTVSGFSITSYGVSWVRQPPGKGIEWIGVIWGDGSTNYHSAIISRISISKDNSKSQVFIKINSIQTDDTATYYCAKG[220.05][174.424]E-FYYAMDYWGQGTSVTVS[157.941]K-IFIVAKTTPPSVYPIAPGSAAQTNSMVTIGCIVKGYFPEPVTVTWNSGSISSGVHTFPAVIQSDIYTISSSVTVPSSTWPSQTVTCNVAHPASSTKVDKKIV[253.234]DC-DCGCKPCICTVPEVSSVFIFPPKPKDVITITITPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTKPREEQF-TFRSVSEIPIMHQDWINGKEFKCRVNSAAFPAPIEKTISKTK[213.066][225.191]APKVYTIPPPKEQMAKDKVSITCMITNFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKINVQKSNWEAGNTFTCSVIHEGIHNHHTEKSISHSPGK"
		 * ; //String A =
		 * "IKESGPGIVAPSKSISITCTVSGFSITSYGVSWVRKPPGKGIEWIGVIWGDGSTNYHSAIISRISISKDNSKSKVFIKINSIKTDDTATYYCAKG[220.05][174.492]E-[310.138]YAMDYWGKGTSVTVSSAKTTPPSVYPIAPGSAAKTNSMVTIGCIVKGYFPEPVTVTWNSGSISSGVHTFPAVIKSDIYTISSSVTVPSSTWPSKTVTCNVAHPASSTKVDKKIV[253.224]D-DCGCKPCIC[320.451]-SSSIVPEVSSVFIFPPKPKDVITITITPKVTCVVVDISKDDPEVKFSWFVDDVEVHTAKTKPREEKINSTFRSVSEIPIMHKDWINGKEFKCRVNSAAFPAPIEKTISKTK[213.066][225.191]APKVYTIPPPKEKMAKDKVSITCMITNFFPEDITVEWKWNGKPAENYKNTKPIMDTDGSYFVYSKINVKKSNWEAGNTFTCSVIHEGIHNHHTEKSISHSPGK"
		 * ; String A =
		 * "NNVMSQS[184.636]SIAVSAGEKVTMSCKSSQSIINSR-KN[276.443]AWYQKKPGQSPKIIIYWASTRESGVPDRFTGSGSGTDFTITISSVQAEDIAVYYCK-ITFGAGTKIEIKS[229.021]-RADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
		 * ; //String A =
		 * "[276.743]NNVMSKS[184.653]SIAVSAGEKVTMSCKSSKSIINSR-RNYIAWYKKKPGKSPKII-[198.045]IIIYWASTRESGVPDRFTGSGSGTDFTITISSVKAEDIAVYYCK-KSYNI[198.071]FGSGTKIEIK-RADAAPTVSIFPPSSEKITSGGASVVCFINNFYPKDINVKWKIDGSERKNGVINSWTDKDSKDSTYSMSSTITITKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
		 * ; //String A =
		 * "[276.743]NNVMSKS[184.653]SIAVSAGEKVTMS[310.407]W[310.205]-NYIAWYKKKPGKSPKII-[198.045]IIIYWASTRESGVPDRFTGSGSGTDFTITISSVKAEDIAVYYCK-KSYNI[198.071]FGSGTKIEIK-RADAAPTVSIFPPSSEKITSGGASVVCFINNFYPKDINVKWKIDGSERKNGVINSWTDKDSKDSTYSMSSTITITKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC"
		 * ; //String A =
		 * "IKESGPGIVAPSKSISITCTVSGFSITSYGVSWVRKPPGKGIEWIGVIWGDGSTNYHSAIISRISISKDNSKSKVFIKINSIKTDDTATYYCAKGGYRFYYAMDYWGKGTSVTVSSAKTTPPSVYPIAPGSAAKTNSMVTIGCIVKGYFPEPVTVTWNSGSISSGVHTFPAVIKSDIYTISSSVTVMSSPRPSETVTCNVAHPASSTKVDIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVITITITPKVTCVVVDISKDDKFSWFVDDVEVHTAKTK[253.149]E-EEQFNSTFRSRSVSELPIMHQDWLNGKEFKCRVNSAAFKTISKTKGRPKAPQVVYTIPPPKEQMAKDKVSLTCMITDFFPEDITVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKLNVQKSNWEAGNTFTCSVLHEGLHNHHTEKSLSHSPGK"
		 * ; //String A =
		 * "WWWWAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSAAAAGVHTFPAVLQSDLYTLSSSVTVPSSTWPSQTVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVLTITLTPKVTCVVVDISKDDPEVQFSWFVDDVEVHTAQTQPREEQFNSTFRSVSELPIMHQDWLNGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPQVYTIPPPKEQMAKDKVSLTCMITDFFPEDITAAAAAVEWQWNGQPAENYKNTQPIMDTDGSYFVYSKLNVQKSNWEAGNTFTCSVLHEGLHNHHTEKSLSHSPGKWWWWW"
		 * ; int Res = AntibodyUtils.GetFullSequenceCorrectness(HC_DB);
		 * System.out.println("Heavy Chain Protein DB:");
		 * System.out.println("  Correct: " + Res); HC_DB =
		 * HC_DB.replace("-",""); System.out.println("  Length: " +
		 * AntibodyUtils.SplitStringtoAA(HC_DB).length);
		 * 
		 * Res = AntibodyUtils.GetFullSequenceCorrectness(HC_Genomic);
		 * System.out.println("Heavy Chain Genomic DB:");
		 * System.out.println("  Correct: " + Res); HC_Genomic =
		 * HC_Genomic.replace("-",""); System.out.println("  Length: " +
		 * AntibodyUtils.SplitStringtoAA(HC_Genomic).length);
		 * 
		 * 
		 * 
		 * Res = AntibodyUtils.GetFullSequenceCorrectnessLC(A);
		 * System.out.println("Light Chain Protein DB:");
		 * System.out.println("  Correct: " + Res); A = A.replace("-","");
		 * System.out.println("  Length: " +
		 * AntibodyUtils.SplitStringtoAA(A).length);
		 * //System.out.println(((double)Res)/AntibodyUtils.aBTLA_LC.length());
		 */

		/*
		 * String A = "FYYAMDYWGKGTSVTVS[158.185]K[298.862]"; String B =
		 * "VSSAKTTPPSVYPIAPGSAAKTNSMVTIGCIVKGYFPEPVTVTWNSGSISSGVHTFPAVIKSDIYTISSSVTVPSSTWPSKTVTCNVAHPASSTKVDKKIVPRDCGCKPCICTVPEVSSVFIFPPKPKDVITITITPKVTCVVVDISKDDPEVKFSWFVDDVEVHTAKTKPREEKINSTFRSVSEIPIMHKDWINGKEFKCRVNSAAFPAPIEKTISKTKGRPKAPKVYTIPPPKEKMAKDKVSITCMITN"
		 * ;
		 * 
		 * int[] PRMsA = AntibodyUtils.GetSeqPRMScaled(A); int[] PRMsB =
		 * AntibodyUtils.GetSeqPRMScaled(B); double[] ScoresA = new
		 * double[PRMsA.length]; double[] ScoresB = new double[PRMsB.length];
		 * for(int i = 0; i < ScoresA.length; ++i) ScoresA[i] = 1; for(int i =
		 * 0; i < ScoresB.length; ++i) ScoresB[i] = 1; System.out.println("A: "
		 * + AntibodyUtils.IntArrayToString(PRMsA)); System.out.println("A: " +
		 * AntibodyUtils.StringArrayToString(AntibodyUtils.SplitStringtoAA(A)));
		 * System.out.println("B: " + AntibodyUtils.IntArrayToString(PRMsB));
		 * System.out.println("B: " +
		 * AntibodyUtils.StringArrayToString(AntibodyUtils.SplitStringtoAA(B)));
		 * 
		 * ConsensusPRMSpectrum Output =
		 * AntibodyUtils.MergePRMsDoubleSkip(PRMsA, ScoresA, PRMsB, ScoresB);
		 * System.out.println("Merged: " +
		 * AntibodyUtils.IntArrayToString(Output.ScaledPeaks));
		 * System.out.println("Seq: " + DeNovo.Reconstruct(Output.ScaledPeaks,
		 * Output.PeakScores));
		 */
		/*
		 * String A = "FFFFAAAA[142.074]DD"; String B = "FAAAAAAD"; //String B =
		 * "[" + Utils.AAMasses['A'-65] + "]AAA[" + (2*Utils.AAMasses['A'-65])
		 * +"]A[" + Utils.AAMasses['D'-65] + "]DDDAAAAA[" +
		 * Utils.AAMasses['A'-65]+"]"; //String B = "AAAA[156.090]G[" +
		 * Utils.AAMasses['D'-65] + "]DDDAAAAA[" +
		 * ((Utils.AAMassesScaled['A'-65])/MASS_SCALE)+"]";
		 * 
		 * System.out.println(A); System.out.println(B); String C =
		 * AntibodyUtils.MergeSequences(A, B); System.out.println(C);
		 * 
		 * String D = "YHSAIISRI[0.995]"; int[] PRMs =
		 * AntibodyUtils.GetSeqPRMScaled(D); System.out.println(D); for(int i =
		 * 0; i < PRMs.length; ++i) System.out.println("Prefix " + i + ":" +
		 * PRMs[i]);
		 * 
		 * D = "YHSALISRLD"; PRMs = AntibodyUtils.GetSeqPRMScaled(D);
		 * System.out.println(D); for(int i = 0; i < PRMs.length; ++i)
		 * System.out.println("Prefix " + i + ":" + PRMs[i]);
		 */
	}

	public static String[] ListSpectrumDir(String inFile) {

		boolean LocalDebug = false;
		String[] fList = Utils.ListDir(inFile);
		ArrayList temp = new ArrayList();
		for (int i = 0; i < fList.length; ++i) {
			String ext = Utils.GetFileExtension(fList[i]).toLowerCase();
			if (LocalDebug)
				System.out.println("Considering file " + fList[i]
						+ " with ext " + ext);
			if (ext.compareTo(".mgf") != 0 && ext.compareTo(".mzxml") != 0
					&& ext.compareTo(".dta") != 0 && ext.compareTo(".pkl") != 0)
				continue;

			String newFile = fList[i];
			temp.add(newFile);
		}
		return Utils.ConvertArraylistToStringArray(temp);
	}

	/*
	 * public static boolean RunCommand2StdOut(String command, String[] env,
	 * String Dir) throws Exception{ System.out.println(">> " + command);
	 * Runtime CurrRuntime = Runtime.getRuntime(); Process P; P =
	 * CurrRuntime.exec(command, env, new File(Dir));
	 * 
	 * InputStream stdin = P.getInputStream(); InputStreamReader isr = new
	 * InputStreamReader(stdin); BufferedReader br = new BufferedReader(isr);
	 * String Line = br.readLine();
	 * 
	 * InputStream stderr = P.getErrorStream(); InputStreamReader esr = new
	 * InputStreamReader(stderr); BufferedReader ebr = new BufferedReader(esr);
	 * String LineE = ebr.readLine();
	 * 
	 * while(Line != null) { System.out.println(Line); Line = br.readLine();
	 * 
	 * }
	 * 
	 * while(LineE != null) { System.err.println(LineE); LineE = ebr.readLine();
	 * }
	 * 
	 * P.waitFor(); P.destroy(); br.close(); isr.close(); esr.close();
	 * ebr.close(); stderr.close(); stdin.close();
	 * 
	 * System.out.println("Command finished"); return true; }
	 */

	/**
	 * Find the position of the substring in the String array aAs. THe match
	 * must be exact allowing Ks to match Qs and Is to match Ls
	 */
	public static int MatchStringToStringArray(String[] aAs, String substring) {

		boolean LocalDebug = false;
		for (int AAIndex = 0; AAIndex < aAs.length - substring.length() + 1; ++AAIndex) {
			for (int strIndex = 0; strIndex < substring.length(); ++strIndex) {
				if (aAs[AAIndex + strIndex].compareTo(substring
						.charAt(strIndex) + "") == 0
						|| (aAs[AAIndex + strIndex].compareTo("K") == 0 && substring
								.charAt(strIndex) == 'Q')
						|| (aAs[AAIndex + strIndex].compareTo("Q") == 0 && substring
								.charAt(strIndex) == 'K')
						|| (aAs[AAIndex + strIndex].compareTo("I") == 0 && substring
								.charAt(strIndex) == 'L')
						|| (aAs[AAIndex + strIndex].compareTo("L") == 0 && substring
								.charAt(strIndex) == 'I')) {
					// If we match, and it's the last one!
					if (strIndex == substring.length() - 1
							&& aAs[AAIndex + strIndex].compareTo(substring
									.charAt(strIndex) + "") == 0) {
						if (LocalDebug)
							System.out.println("Found " + substring
									+ " at position " + AAIndex + "!");
						return AAIndex;
					}
				} else
					break;

			}
		}
		return -1;
	}

	public static int[] GetLongestAlignmentRelaxedPRMs(int[] prms,
			int[] seedPRMs, double tol) {
		int[][] Alignment = null;// GetLongestAlignmentStrictRelaxedPRMs
		int[][] prevCell = null;
		Object[] details = AntibodyUtils
				.GetLongestAlignmentRelaxedDetailedPRMs(prms, seedPRMs, tol);
		Alignment = (int[][]) (details[1]);
		prevCell = (int[][]) (details[0]);
		int OverallBestScore = 0;
		int BestA = 0, BestB = 0;
		for (int i = 0; i < Alignment.length; ++i) {
			for (int j = 0; j < Alignment[0].length; ++j) {
				if (Alignment[i][j] > OverallBestScore) {
					OverallBestScore = Alignment[i][j];
					BestA = i;
					BestB = j;
				}
			}
		}
		int[] ret = new int[3];
		ret[0] = OverallBestScore;
		int prevA = BestA;
		int prevB = BestB;

		while (prevCell[prevA][prevB] >= 0) {
			// System.out.println("prevs: [" + prevA + "," + prevB + "]");
			int currCell = prevCell[prevA][prevB];
			prevA = currCell / seedPRMs.length;
			prevB = currCell % seedPRMs.length;
		}
		ret[1] = prevA;
		ret[2] = prevB;
		return ret;
	}

	public static Object[] GetLongestAlignmentRelaxedDetailedSeqs(
			String seedSequence, String updatedTemplateSeq, double prmTolerance) {
		int[] PRMsA = AntibodyUtils.GetSeqPRMScaled(seedSequence);
		int[] PRMsB = AntibodyUtils.GetSeqPRMScaled(updatedTemplateSeq);

		return AntibodyUtils.GetLongestAlignmentDetailedPRMs(PRMsA, PRMsB,
				prmTolerance, Integer.MAX_VALUE, Integer.MAX_VALUE);
	}

	public static Object[] GetLongestAlignmentRelaxedDetailedSeqs(
			String seedSequence, String updatedTemplateSeq,
			double prmTolerance, boolean LocalDebug) {
		int[] PRMsA = AntibodyUtils.GetSeqPRMScaled(seedSequence);
		int[] PRMsB = AntibodyUtils.GetSeqPRMScaled(updatedTemplateSeq);

		return AntibodyUtils.GetLongestAlignmentDetailedPRMs(PRMsA, PRMsB,
				prmTolerance, Integer.MAX_VALUE, Integer.MAX_VALUE);
	}

	public static Object[] GetLongestAlignmentDetailedPRMs(int[] PRMsA,
			int[] PRMsB, double Tolerance, int numSkipA, int numSkipB) {
		int[][] Alignment = new int[PRMsA.length][PRMsB.length];
		int[][] prevCell = new int[PRMsA.length][PRMsB.length];
		boolean LocalDebug = false;
		for (int a = 0; a < PRMsA.length; ++a) {
			for (int b = 0; b < PRMsB.length; ++b) {
				Alignment[a][b] = 1;
				prevCell[a][b] = -1;
			}
		}
		for (int a = 1; a < PRMsA.length; ++a) {

			int AMass = PRMsA[a];
			for (int b = 1; b < PRMsB.length; ++b) {
				int BMass = PRMsB[b];
				int BestScore = 1;

				for (int prevA = Math.max(0, a - numSkipA); prevA < a; ++prevA) {
					int prevAMass = PRMsA[prevA];
					for (int prevB = Math.max(0, b - numSkipB); prevB < b; ++prevB) {
						int prevBMass = PRMsB[prevB];

						if (AntibodyUtils.EquivalentScaled(BMass - prevBMass,
								AMass - prevAMass, Tolerance)) {
							if (LocalDebug) {
								System.out.println("Prevs [" + prevAMass + ","
										+ prevBMass + "] -> [" + AMass + ","
										+ BMass + "]");
							}
							if (Alignment[prevA][prevB] + 1 > BestScore) {
								BestScore = Alignment[prevA][prevB] + 1;
								prevCell[a][b] = prevA * PRMsB.length + prevB;
							}

						}
					}

				}
				Alignment[a][b] = BestScore;

			}
		}
		Object[] ret = new Object[2];
		ret[0] = prevCell;
		ret[1] = Alignment;
		return ret;
	}

	public static Object[] GetLongestAlignmentDetailedSeqs(String seedSequence,
			String updatedTemplateSeq, double prmTolerance, int i, int j) {
		int[] PRMsA = AntibodyUtils.GetSeqPRMScaled(seedSequence);
		int[] PRMsB = AntibodyUtils.GetSeqPRMScaled(updatedTemplateSeq);
		return AntibodyUtils.GetLongestAlignmentDetailedPRMs(PRMsA, PRMsB,
				prmTolerance, i, j);
	}

	/**
	 * Creates a constraint file from a concatenated DB for GenoMS. FASTA
	 * entries must contain the DBIndex and the RootName such as
	 * >4.DBRootName.XXXXXXXXXX
	 * 
	 * @param combinedDBName
	 *            The FASTA file created with the appropriate FASTA entries
	 * @param constraintFileName
	 *            The constraint file to be created
	 * @return Returns true if the constraint file creation is successful, false
	 *         otherwise
	 */
	public static boolean CreateConstraintsFromDB(String combinedDBName,
			String constraintFileName) {
		if (constraintFileName == null || combinedDBName == null) {
			return false;
		}
		String trieDB = Utils.GetFileNameNoExtension(combinedDBName) + ".trie";
		TrieDB.prepDB(combinedDBName, trieDB);
		TrieDB TrieDatabase = new TrieDB(trieDB);

		Template[] AllTemplates = Template.CreateTemplatesFromDB(trieDB);
		if (AllTemplates == null) {
			return false;
		}
		FragmentDB.WriteConstraintsToFileBinary(AllTemplates, constraintFileName);

		return true;
	}

	static String makePathAbsolute(String fileName, String dir) {
		File Test = new File(fileName);
		if (!Test.isAbsolute()) {
			Test = new File(dir + File.separator + fileName);
			try {
				return Test.getCanonicalPath();
			} catch (Exception E) {
				ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
						E.toString());
			}
		}
		return fileName;
	}

	public static String CreateDatabase(String[] DatabaseRootName,
			String ResourceDir, boolean Debug) {

		BufferedReader Reader = null;
		String CombinedDBFasta = Utils.GetFilePath(DatabaseRootName[0])
				+ File.separator;

		for (int i = 0; i < DatabaseRootName.length; ++i) {
			makePathAbsolute(DatabaseRootName[i], ResourceDir);
			CombinedDBFasta += Utils
					.GetBaseNameNoExtension(DatabaseRootName[i]) + "_";
		}
		CombinedDBFasta += "combined.fa";
		BufferedWriter Writer = null;
		try {
			Writer = new BufferedWriter(new FileWriter(CombinedDBFasta));
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, CombinedDBFasta);
		}
		if (Debug)
			System.out.println("Attempting to write: " + CombinedDBFasta);

		int seqClassIndex = 0;
		for (int j = 0; j < DatabaseRootName.length; ++j) {
			int DBIndex = 1; // Database index

			while (true) {
				String CurrFile = DatabaseRootName[j]
						+ AntibodyUtils.CLASS_DELIM + DBIndex;

				try {
					Reader = new BufferedReader(new FileReader(CurrFile));
				} catch (IOException E) {
					// System.err.println(E.getMessage());
					break;
				}
				seqClassIndex += 1;
				if (Debug)
					System.out.println("Loading DB " + CurrFile);
				String Line = null;
				try {
					Line = Reader.readLine();
				} catch (IOException E) {
					ErrorThrower.ThrowError(6, CurrFile);
				}
				while (Line != null) {
					// Write the line
					if (Line.length() > 0 && Line.charAt(0) == '>')
						Line = ">"
								+ seqClassIndex
								+ AntibodyUtils.CLASS_DELIM
								+ Utils.GetBaseNameNoExtension(DatabaseRootName[j])
								+ AntibodyUtils.CLASS_DELIM
								+ Line.substring(1, Line.length());
					if (Line.length() > 0) {
						try {
							Writer.write(Line);
							Writer.newLine();
						} catch (IOException E) {
							ErrorThrower.ThrowError(7, CombinedDBFasta);
						}
					}
					// Read another line
					try {
						Line = Reader.readLine();
					} catch (IOException E) {
						ErrorThrower.ThrowError(6, CurrFile);
					}

				}
				if (Line == null) {
					try {
						Reader.close();
					} catch (IOException E) {
						ErrorThrower.ThrowError(8, CurrFile);
					}
				}
				DBIndex += 1;

			}
		}
		try {
			Writer.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(8, CombinedDBFasta);
		}

		return CombinedDBFasta;
	}

	/**
	 * Concatenates an array of files
	 */
	public static boolean concatenateFiles(String[] inputFiles,
			String outputFileName) {
		BufferedReader Reader = null;
		FileWriter Writer = Utils.openFileWriter(outputFileName);
		boolean Debug = false;

		for (int j = 0; j < inputFiles.length; ++j) {

			String CurrFile = inputFiles[j];
			int DBIndex = Integer.parseInt(Utils.GetFileExtension(CurrFile)
					.substring(1));
			System.out.println("Looking for file " + CurrFile);
			Reader = Utils.openBufferedReader(CurrFile);

			if (Debug)
				System.out.println("Loading DB " + CurrFile);
			String Line = null;
			Line = Utils.readNextLine(Reader, CurrFile);
			while (Line != null) {
				// Write the line
				if (Line.length() > 0 && Line.charAt(0) == '>')
					Line = ">" + DBIndex + AntibodyUtils.CLASS_DELIM
							+ Utils.GetBaseNameNoExtension(inputFiles[j])
							+ AntibodyUtils.CLASS_DELIM
							+ Line.substring(1, Line.length());
				else // This is a sequence line
				{
					// Make sure all letters are amino acids
					Line = Utils.AminoAcidify(Line);
				}
				if (Line.length() > 0) {
					Utils.writeLine(Writer, outputFileName, Line + "\n");
				}
				// Read another line
				Line = Utils.readNextLine(Reader, CurrFile);
			}
			Utils.closeBufferedReader(Reader, CurrFile);

		}

		Utils.closeFileWriter(Writer, outputFileName);
		return true;
	}

	/**
	 * Performs the alignment between seq and dbSeq, returns the score and the
	 * start and end indexes
	 * 
	 * @param seq
	 * @param dbSeq
	 * @return
	 */
	public static int[] getBestSeqAlignment(String seq, String dbSeq) {
		int[][] scoreTable = new int[seq.length() + 1][dbSeq.length() + 1];
		int[][] backtrackTable = new int[seq.length() + 1][dbSeq.length() + 1];

		for (int i = 0; i <= seq.length(); ++i) {
			scoreTable[i][0] = i * GAP_SCORE;
			backtrackTable[i][0] = (i - 1) * dbSeq.length();
		}
		for (int j = 0; j <= dbSeq.length(); ++j) {
			scoreTable[0][j] = j * GAP_SCORE;
			backtrackTable[0][j] = j - 1;
		}
		for (int i = 1; i <= seq.length(); ++i) {
			for (int j = 1; j <= dbSeq.length(); ++j) {

				int skipAScore = 0;
				int skipAPrev = 0;

				int skipBScore = 0;
				int skipBPrev = 0;

				int matchScore = 0;
				int matchPrev = 0;

				// Try a skip in dbSeq
				skipAScore = scoreTable[i - 1][j] + GAP_SCORE;
				skipAPrev = (i - 1) * (dbSeq.length() + 1) + j;

				// Try a skip in seq
				skipBScore = scoreTable[i][j - 1] + GAP_SCORE;
				skipBPrev = i * (dbSeq.length() + 1) + (j - 1);

				// Try a match
				matchScore = scoreTable[i - 1][j - 1];
				matchPrev = (i - 1) * (dbSeq.length() + 1) + (j - 1);

				if (seq.charAt(i - 1) == dbSeq.charAt(j - 1))
					matchScore += MATCH_SCORE;
				else
					matchScore += MISMATCH_SCORE;

				// Put the score in the table
				if (skipAScore > skipBScore && skipAScore > matchScore) {
					scoreTable[i][j] = skipAScore;
					backtrackTable[i][j] = skipAPrev;
					// System.out.println("Best score is skipping dbSeq");
				} else if (skipBScore > skipAScore && skipBScore > matchScore) {
					scoreTable[i][j] = skipBScore;
					backtrackTable[i][j] = skipBPrev;
					// System.out.println("Best score is skipping seq");
				} else if (matchScore >= skipAScore && matchScore >= skipBScore) {
					scoreTable[i][j] = matchScore;
					backtrackTable[i][j] = matchPrev;
					// System.out.println("Best score is match/mismatch");
				}
				/*
				 * if(i == j) { System.out.println("Comparing " +
				 * seq.charAt(i-1) + " to " + dbSeq.charAt(j-1));
				 * System.out.println("INSERT: " + skipAScore + ", " +
				 * skipAPrev); System.out.println("DELETE: " + skipBScore + ", "
				 * + skipBPrev); System.out.println("MATCH: " + matchScore +
				 * ", " + matchPrev); System.out.println(" For [" + i + "][" + j
				 * + "] has score " + scoreTable[i][j] + " with previous " +
				 * backtrackTable[i][j] + " [" +
				 * backtrackTable[i][j]/(dbSeq.length()+1) + "][" +
				 * backtrackTable[i][j]%(dbSeq.length()+1) + "]");
				 * Utils.WaitForEnter(); }
				 */
				// Utils.WaitForEnter();
			}
			// Utils.WaitForEnter();
		}

		// Get the results
		int currR = seq.length();
		int currC = dbSeq.length();

		int[] ret = new int[5];
		ret[0] = scoreTable[currR][currC];
		ret[1] = seq.length() - 1; // start of alignment in sequence 1
		ret[2] = 0; // end of alignment in sequence 1
		ret[3] = dbSeq.length() - 1; // start of alignment in sequence 2
		ret[4] = 0; // end of alignment in sequence 2

		int prevIndex = backtrackTable[currR][currC];

		while (true) {
			int r = prevIndex / (dbSeq.length() + 1);
			int c = prevIndex % (dbSeq.length() + 1);
			// System.out.println("[" + r + "][" + c + "] = " +
			// scoreTable[r][c]);

			if (scoreTable[currR][currC] > scoreTable[r][c]) {
				// See if we need to update the start or end
				if (currR > ret[2])
					ret[2] = currR;
				if (currR < ret[1])
					ret[1] = currR;
				if (currC > ret[4])
					ret[4] = currC;
				if (currC < ret[3])
					ret[3] = currC;
			}

			if (prevIndex == 0)
				break;

			currR = r;
			currC = c;

			prevIndex = backtrackTable[currR][currC];

		}
		ret[1] -= 1;
		ret[3] -= 1;
		return ret;
	}

	

	public static Object[] createSeedsFromDBSearchResults(
			ArrayList<Peptide> dbPepObj, TrieDB trieDatabase,
			boolean isShuffled, boolean localDebug) {
		Hashtable ProteinIDs = new Hashtable(); // ProteinID -> ArrayList of
												// Peptide sequences
		Hashtable Peptides = new Hashtable(); // Peptide -> ArrayList of Peptide
												// indexes
		Hashtable Peptide2ProteinIDs = new Hashtable(); // Peptide -> ArrayList
														// of locations of
														// peptide

		for (int i = 0; i < dbPepObj.size(); ++i) {
			Peptide CurrAnnotation = dbPepObj.get(i);

			String CurrUnModdedAnn = CurrAnnotation.getUnModdedPeptideSeq();
			if (Peptides.containsKey(CurrUnModdedAnn)) {
				ArrayList Temp = (ArrayList) (Peptides.get(CurrUnModdedAnn));
				Temp.add(new Integer(i));
				Peptides.put(CurrUnModdedAnn, Temp);

			} else {
				ArrayList Temp = new ArrayList();
				Temp.add(new Integer(i));

				ArrayList Locations = null;

				Locations = trieDatabase.GetAllLocations(CurrUnModdedAnn);

				if (Locations.size() == 0)
					continue;
				// System.out.println("[" + i + "]: " +
				// CurrAnnotation.getPeptideSeq() + "-" + Locations.size());
				Peptides.put(CurrUnModdedAnn, Temp);

				Peptide2ProteinIDs.put(CurrUnModdedAnn, Locations);
				for (int LocIndex = 0; LocIndex < Locations.size(); ++LocIndex) {
					Object[] Loc = (Object[]) (Locations.get(LocIndex));
					int ProteinID = ((Integer) (Loc[TrieDB.TrieLocationColumns.ProteinID]))
							.intValue();

					// Now that we are using MSGF+, there will be no shuffled
					// proteins
					// String ProteinName = (String)
					// (Loc[TrieDB.TrieLocationColumns.ProteinName]);
					// ASSUMPTION: Shuffled file is of the form
					// Shuffled1*True1*Shuffled2*True2...
					/*
					 * if (isShuffled) ProteinID = ProteinID / 2; if
					 * (ProteinName.substring(0, 3).compareTo("XXX") == 0) {
					 * continue; }
					 */

					ArrayList PepList;

					if (!ProteinIDs.containsKey(new Integer(ProteinID))) {
						PepList = new ArrayList();
					} else {

						PepList = (ArrayList) (ProteinIDs.get(new Integer(
								ProteinID)));
					}
					PepList.add(CurrUnModdedAnn);
					ProteinIDs.put(new Integer(ProteinID), PepList);
				}
			}
			// if (LocalDebug)
			// AntibodyUtils.WaitForEnter();
		}
		// System.out.println("Loaded " + ProteinIDs.size() + " templates and "
		// + Peptides.size() + " peptides");

		Object[] ret = new Object[3];
		ret[0] = ProteinIDs;
		ret[1] = Peptides;
		ret[2] = Peptide2ProteinIDs;
		return ret;
	}

	public static Object[] createSeedsFromDBSearchResults(
			SpectrumNode[] allSpectra, TrieDB TrieDatabase, boolean isShuffled,
			boolean LocalDebug) {

		Hashtable ProteinIDs = new Hashtable(); // ProteinID -> ArrayList of
												// Peptide sequences
		Hashtable Peptides = new Hashtable(); // Peptide -> ArrayList of
												// SpectrumNode indexes
		Hashtable Peptide2ProteinIDs = new Hashtable(); // Peptide -> ArrayList
														// of locations of
														// peptide

		for (int i = 0; i < allSpectra.length; ++i) {
			SpectrumNode CurrAnnotation = allSpectra[i];
			if (CurrAnnotation.getModdedDBAnnotation() == null)
				continue;
			String CurrUnModdedAnn = CurrAnnotation.getUnModdedDBAnnotation();
			if (Peptides.containsKey(CurrUnModdedAnn)) {
				ArrayList Temp = (ArrayList) (Peptides.get(CurrUnModdedAnn));
				Temp.add(new Integer(i));
				Peptides.put(CurrUnModdedAnn, Temp);
				if (LocalDebug) {
					System.out.println("Already saw " + CurrUnModdedAnn);
					System.out.flush();
				}

			} else {
				ArrayList Temp = new ArrayList();
				Temp.add(new Integer(i));

				ArrayList Locations = null;

				Locations = TrieDatabase.GetAllLocations(CurrUnModdedAnn);
				// if(this.RunMSAlignmentFlag)
				// Locations.addAll(this.TrieDatabase.GetAllLocationsSingleMod(CurrAnnotation.UnModdedAnn));
				if (LocalDebug) {
					System.out.println("New peptide: " + CurrUnModdedAnn);
					System.out.flush();
				}

				// if (Locations.size() == 1) {
				// CurrAnnotation.IsUnique = true;
				// }
				if (LocalDebug) {
					System.out.println("Found locations: " + Locations.size());
					System.out.flush();
				}
				if (Locations.size() == 0)
					continue;
				Peptides.put(CurrUnModdedAnn, Temp);

				Peptide2ProteinIDs.put(CurrUnModdedAnn, Locations);
				for (int LocIndex = 0; LocIndex < Locations.size(); ++LocIndex) {
					Object[] Loc = (Object[]) (Locations.get(LocIndex));
					int ProteinID = ((Integer) (Loc[TrieDB.TrieLocationColumns.ProteinID]))
							.intValue();
					String ProteinName = (String) (Loc[TrieDB.TrieLocationColumns.ProteinName]);
					// ASSUMPTION: Shuffled file is of the form
					// Shuffled1*True1*Shuffled2*True2...
					if (isShuffled)
						ProteinID = ProteinID / 2;
					if (ProteinName.substring(0, 3).compareTo("XXX") == 0) {
						if (LocalDebug) {
							System.out.println("Protein " + Loc[0] + " - "
									+ ProteinName);
							System.out.flush();
						}
						continue;
					}

					ArrayList PepList;

					if (!ProteinIDs.containsKey(new Integer(ProteinID))) {
						if (LocalDebug && !LocalDebug)
							System.out
									.println("This is the first peptide for this protein!");
						PepList = new ArrayList();
					} else {

						PepList = (ArrayList) (ProteinIDs.get(new Integer(
								ProteinID)));
						if (LocalDebug && !LocalDebug)
							System.out.println("This protein already has "
									+ PepList.size() + " peptides");
					}
					PepList.add(CurrUnModdedAnn);
					ProteinIDs.put(new Integer(ProteinID), PepList);
				}
			}
			// if (LocalDebug)
			// AntibodyUtils.WaitForEnter();
		}
		// System.out.println("Loaded " + ProteinIDs.size() + " templates and "
		// + Peptides.size() + " peptides");

		Object[] ret = new Object[3];
		ret[0] = ProteinIDs;
		ret[1] = Peptides;
		ret[2] = Peptide2ProteinIDs;
		return ret;
	}

	public static String getMostPopular(String[] currPeps) {
		if (currPeps == null || currPeps.length == 0)
			return null;
		if (currPeps.length == 1)
			return currPeps[0];

		Hashtable cnt = new Hashtable();
		int max = 0;
		String best = null;
		for (int i = 0; i < currPeps.length; ++i) {
			int currCount = 0;
			if (currPeps[i] == null)
				continue;
			if (cnt.containsKey(currPeps[i]))
				currCount = ((Integer) (cnt.get(currPeps[i]))).intValue();
			currCount += 1;
			cnt.put(currPeps[i], new Integer(currCount));

			if (currCount > max) {
				max = currCount;
				best = currPeps[i];
			}

		}
		return best;

	}

	public static String getMostPopular(ArrayList<String> currPeps) {
		if (currPeps == null || currPeps.size() == 0)
			return null;
		if (currPeps.size() == 1)
			return currPeps.get(0);

		Hashtable cnt = new Hashtable();
		int max = 0;
		String best = null;
		for (int i = 0; i < currPeps.size(); ++i) {
			int currCount = 0;
			String currPep = currPeps.get(i);
			if (currPep == null || currPep.equals(""))
				continue;
			if (cnt.containsKey(currPep))
				currCount = ((Integer) (cnt.get(currPep))).intValue();
			currCount += 1;
			cnt.put(currPep, new Integer(currCount));

			if (currCount > max) {
				max = currCount;
				best = currPep;
			}

		}
		return best;

	}

	public static double getMassFromSeq(String sequence) {
		double TotalMass = 0;
		double CurrMass;
		for (int i = 0; i < sequence.length(); ++i) {
			if (sequence.charAt(i) == '[') {
				int End = sequence.indexOf(']', i);

				String AAStr = sequence.substring(i + 1, End);
				// System.out.println("GetSeqMass: [ at " + i + ", ] at " + End
				// + " : " + AAStr);
				TotalMass += Double.parseDouble(AAStr);
				i = End;
				continue;
			}
			int AAIndex = sequence.charAt(i) - 65;
			if (AAIndex < 0 || AAIndex > AAMasses.length) {
				System.err.println("ERROR: Unrecognized char in AA Sequence: "
						+ sequence.charAt(i) + " in " + sequence);
				return 0;
			}
			CurrMass = (AAMasses[AAIndex]);
			if (CurrMass > 0)
				TotalMass += CurrMass;
			else {
				// System.out.println("WARNING: " + Sequence.charAt(i) + " = " +
				// CurrMass);
				ErrorThrower.ThrowWarning(11, "" + sequence.charAt(i) + " = "
						+ CurrMass);
			}
		}
		return TotalMass;
	}
	
	
	/**
	 * Returns the MZ value for this sequence with one charge (so the sum of the masses plus water plus Hydrogen).
	 * @param seedSequence
	 * @return
	 */
	public static double getPepMHFromSeq(String seedSequence) {
		double mass = AntibodyUtils.getMassFromSeq(seedSequence);
		mass += Utils.MHMass;
		return mass;
	}

	
	/**
	 * Takes a peptide with modifications of the form AC+57D and converts them to A(C,57)D
	 * @param modifiedPeptide
	 * @return Returns the updated string
	 */
	public static String ConvertMSGFPlusPeptide2SPS(String modifiedPeptide) {
		String[] splitPep = Utils.splitSequenceWithModsToArray(modifiedPeptide);
		String ret = "";
		for(int i = 0; i < splitPep.length; ++i) {
			if(splitPep[i].length() == 1) //If it's just a single amino acid, just add it to the returning peptide
				ret += splitPep[i];
			else {
				ret += "(" + splitPep[i].charAt(0) + ",";
				if(splitPep[i].charAt(1) == '+') 
					ret += splitPep[i].substring(2) + ")";//If it's a positive mod, drop the '+'
				else
					ret += splitPep[i].substring(1) + ")";//If it's a negative mod, keep the '-'
			}
		}
		return ret;
	}

}
