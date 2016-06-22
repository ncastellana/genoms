package antibody;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

import basicUtils.InspectAnnotation;
import basicUtils.Utils;
import errorUtils.ErrorThrower;

public class FDRCalculator {

	public static final int MSGF = 0;
	public static final int FScore = 1;
	public static final int MinLength = 7;

	public static void performFDRCutoffWithDecoys(String InputFile,
			String OutputFile, double PValueCutOff, int ScoreMode) {
		// Load Inspect results
		System.out.println("Loading results from " + InputFile
				+ " for FDR calculation");
		InspectAnnotation[] RawAnns = null; // InspectAnnotation.LoadInspectResultsFile(InputFile);
		// String[] colHeaders = RawAnns[0].columnHeaders;

		if (ScoreMode == MSGF) {
			RawAnns = InspectAnnotation
					.LoadInspectResultsFileIncreasingOrderTwoVotes(InputFile,
							"SpecProb");
		} else if (ScoreMode == FScore) {
			RawAnns = InspectAnnotation
					.LoadInspectResultsFileDecreasingOrderTwoVotes(InputFile,
							"F-Score");
		} else {
			ErrorThrower
					.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
							"Sorry score modes besides MSGF and FScore are not currently supported!");
		}
		System.out.println("Loaded " + RawAnns.length + " Inspect results");

		if (RawAnns.length < 100) {
			System.err
					.println("WARNING: Fewer than 100 PSMs to compute FDR, calculations may be erroneous!");
		}

		InspectAnnotation[] FilteredAnns = null;
		FilteredAnns = FDRCalculator.ApplyQValueCutOff(RawAnns, PValueCutOff,
				ScoreMode, false);

		// Write Results
		FileWriter Writer = null;
		boolean wroteHeader = false;
		try {
			Writer = new FileWriter(OutputFile);
			// Writer.write(InspectAnnotation.SplicedHeaderLine);
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, OutputFile);
		}
		for (int i = 0; i < FilteredAnns.length; ++i) {

			try {

				if (!wroteHeader) {
					Writer.write("#"
							+ Utils.JoinStringArray(
									FilteredAnns[i].columnHeaders, "\t") + "\n");
					wroteHeader = true;
				}
				Writer.write(FilteredAnns[i].toString() + "\n");
			} catch (IOException E) {
				ErrorThrower.ThrowError(7, OutputFile);
			}
		}
		try {
			Writer.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(8, OutputFile);
		}
	}

	private static InspectAnnotation[] ApplyQValueCutOff(
			InspectAnnotation[] SortedAnns, double fdrCutOff, int scoreMode,
			boolean reportFalse) {

		Hashtable<Double,Double> fdr = new Hashtable<Double,Double>();

		ArrayList<InspectAnnotation> finalList = new ArrayList<InspectAnnotation>();
		boolean debug = false;
		if (debug)
			System.out.println("Applying QValue cutoff!");
		// For each charge state, we consider the distributions differently
		for (int Charge = 1; Charge < 4; ++Charge) {
			System.out.println("Considering charge state: " + Charge);

			int cumTrue = 0;
			int cumFalse = 0;
			int Count = 0;

			/**
			 * First determine the pdf + cdf (the number of true and false
			 * matches at each score)
			 */
			for (int i = 0; i < SortedAnns.length; ++i) {
				if (SortedAnns[i].Charge != Charge)
					continue;

				if (debug)
					System.out.println(SortedAnns[i].OriginalLine);
				double Score = FDRCalculator.GetScore(SortedAnns[i], scoreMode);
				if (SortedAnns[i].ProteinName.length() < 3)
					System.out.println(SortedAnns[i].toString());
				if (Utils.IsDecoyProtein(SortedAnns[i].ProteinName)) {
					cumFalse += 1;
				} else {
					cumTrue += 1;
				}
				double currFDR = ((double) cumFalse) / cumTrue;
				fdr.put(new Double(Score), new Double(currFDR));
				Count += 1;
			}
			System.out.println("At charge state " + Charge + " we found "
					+ Count + " PSMs");

			/**
			 * Now we find the qvalue for each guy
			 */
			double minFDR = 1.0;
			boolean found = false;
			for (int i = SortedAnns.length - 1; i >= 0; --i) {
				if (SortedAnns[i].Charge != Charge)
					continue;

				double Score = FDRCalculator.GetScore(SortedAnns[i], scoreMode);
				double currFDR = ((Double) (fdr.get(new Double(Score))))
						.doubleValue();
				if (currFDR < minFDR)
					minFDR = currFDR;
				if (minFDR == 0)
					minFDR = 0.0001;
				SortedAnns[i].FDR = minFDR;
				if (minFDR < fdrCutOff) {
					if (!found) {
						System.out.println("Score cutoff: " + Score);
						found = true;
					}
					finalList.add(0, SortedAnns[i]);
				}
			}
		}
		System.out.println("Accepted " + finalList.size() + " PSMs");
		return InspectAnnotation.ArrayListToList(finalList);
	}

	private static InspectAnnotation[] ApplyCutOff(
			InspectAnnotation[] SortedAnns, double PValueCutOff, int ScoreMode) {
		ArrayList<InspectAnnotation> FilteredResults = new ArrayList<InspectAnnotation>();
		Hashtable<Double,int[]> ScoreCounts = new Hashtable<Double,int[]>();

		boolean LocalDebug = false;
		if (LocalDebug)
			System.out.println("BLAH!!: " + SortedAnns.length);
		// For each charge state, we consider the distributions differently
		// for(int Charge = 1; Charge < 4; ++Charge)
		// {
		// System.out.println("Considering charge state: " + Charge);
		ScoreCounts.clear();
		int TrueCount = 0;
		int FalseCount = 0;

		// Populate the hashtable
		for (int i = 0; i < SortedAnns.length; ++i) {
			// if(SortedAnns[i].Charge != Charge)
			// continue;

			double Score = FDRCalculator.GetScore(SortedAnns[i], ScoreMode);
			if (LocalDebug) {
				System.out.println("Considering New Ann: "
						+ SortedAnns[i].toString());
				System.out.println("Score: " + Score);
			}
			int[] Counts = null;
			if (ScoreCounts.containsKey(new Double(Score))) {
				if (LocalDebug) {
					System.out.println("We've seen this score before!!");
					// AntibodyUtils.WaitForEnter();
				}
				Counts = (int[]) (ScoreCounts.get(new Double(Score)));
			} else {
				if (LocalDebug) {
					System.out
							.println("This is the first time we've seen this score!");
					// AntibodyUtils.WaitForEnter();
				}
				Counts = new int[2];
				Counts[0] = TrueCount;
				Counts[1] = FalseCount;
			}
			if (SortedAnns[i].ProteinName.substring(0, 3).compareTo("XXX") == 0) {
				if (LocalDebug)
					System.out.println("FALSE PROTEIN!!");
				Counts[1] += 1;
				FalseCount += 1;
			} else {
				if (LocalDebug)
					System.out.println("TRUE PROTEIN!");
				Counts[0] += 1;
				TrueCount += 1;
			}
			ScoreCounts.put(new Double(Score), Counts);
			// if(LocalDebug)
			// AntibodyUtils.WaitForEnter();
		}

		Hashtable<String,Integer> ScansSeen = new Hashtable<String,Integer>();
		double WorstScore = -1;
		if (ScoreMode == FDRCalculator.MSGF)
			WorstScore = Double.MIN_VALUE;
		else if (ScoreMode == FDRCalculator.FScore)
			WorstScore = Double.MAX_VALUE;
		for (int i = 0; i < SortedAnns.length; ++i) {
			// if(SortedAnns[i].Charge != Charge)
			// continue;
			// String[] CurrScan =
			// {AntibodyUtils.GetBaseName(SortedAnns[i].SpectrumFile),"" +
			// SortedAnns[i].ScanNumber};
			String CurrScan = Utils.GetBaseName(SortedAnns[i].SpectrumFile)
					+ "_" + SortedAnns[i].ScanNumber;
			if (ScansSeen.containsKey(CurrScan))
				continue;
			ScansSeen.put(CurrScan, new Integer(1));
			if (SortedAnns[i].ProteinName.substring(0, 3).compareTo("XXX") == 0
					|| SortedAnns[i].Length < FDRCalculator.MinLength)
				continue;
			double Score = FDRCalculator.GetScore(SortedAnns[i], ScoreMode);
			int[] Counts = (int[]) (ScoreCounts.get(new Double(Score)));
			SortedAnns[i].FDR = Math.max(0.0001, ((double) (Counts[1]))
					/ ((double) (Counts[0] + Counts[1])));

			if (SortedAnns[i].FDR <= PValueCutOff) {
				if (LocalDebug) {
					System.out.println("For Score: " + Score + ", FDR: "
							+ SortedAnns[i].FDR);
					SortedAnns[i].DebugPrint();
					System.out.println("KEEPING IT!");
					Utils.WaitForEnter();
				}
				if (ScoreMode == FDRCalculator.MSGF
						&& SortedAnns[i].SpecProb > WorstScore)
					WorstScore = SortedAnns[i].SpecProb;
				if (ScoreMode == FDRCalculator.FScore
						&& SortedAnns[i].FScore < WorstScore)
					WorstScore = SortedAnns[i].FScore;
				FilteredResults.add(SortedAnns[i]);
			}
		}
		// System.out.println("Total for ChargeState: " + (TrueCount +
		// FalseCount));
		System.out.println("Score Cutoff: " + WorstScore);
		// AntibodyUtils.WaitForEnter();
		// }
		return InspectAnnotation.ArrayListToList(FilteredResults);
	}

	public static double GetScore(InspectAnnotation Ann, int ScoreMode) {
		if (ScoreMode == FDRCalculator.MSGF)
			return Ann.SpecProb;
		if (ScoreMode == FDRCalculator.FScore)
			return Ann.FScore;
		System.err
				.println("ERROR: No other form of scoring is available than MSGF or FScore:"
						+ FDRCalculator.MSGF + " or " + FDRCalculator.FScore);
		return 0;
	}

	public static void main(String[] args) {
		String in = "/Users/Natalie/Documents/workspace/ImmunoSeq/TestData//Test1Output_CombinedInspectOutput.out";
		String out = "/Users/Natalie/Documents/workspace/ImmunoSeq/TestData//Test1Output_CombinedInspectOutput.p01.test.out";
		FDRCalculator.performFDRCutoffWithDecoys(in, out, 0.01, MSGF);
	}
}
