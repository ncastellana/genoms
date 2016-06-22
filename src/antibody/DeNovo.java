package antibody;

import java.util.ArrayList;

import basicUtils.Utils;

import errorUtils.ErrorThrower;

/**
 * This class is responsible for reconstructing the amino acid sequence from a
 * list of Prefix Residue Masses (PRMs).
 * 
 * @author natalie
 * 
 */
public class DeNovo {

	private static boolean Debug = false;

	/*
	 * public static int[] GetTagSignificance(int TagLength, double TagScore,
	 * int[] ScaledPeaks, double[] PeakScores, double MinTagScore) { int
	 * BetterScoringTags = 0; int TotalTags = 0;
	 */
	/**
	 * This variable represents the graph structure. There is an edge between
	 * two nodes, i and j, in the graph if adjacencyMatrix[i][j] == true
	 */
	/*
	 * boolean[][] adjacencyMatrix =
	 * BuildAdjacencyMatrixWithDoubles(ScaledPeaks);
	 * 
	 * int[] Ret = { BetterScoringTags, TotalTags }; return Ret; }
	 */
	/*
	 * public static int[] GetTagSignificance_OLD(int TagLength, double
	 * TagScore, int[] ScaledPeaks, double[] PeakScores, double MinTagScore) {
	 * int BetterScoringTags = 0; int TotalTags = 0;
	 */
	/**
	 * This variable represents the graph structure. There is an edge between
	 * two nodes, i and j, in the graph if adjacencyMatrix[i][j] == true
	 */
	/*
	 * boolean[][] adjacencyMatrix =
	 * BuildAdjacencyMatrixWithDoubles(ScaledPeaks, Tolerance); int[] distances
	 * = new int[adjacencyMatrix.length];
	 * 
	 * ArrayList[] paths = new ArrayList[distances.length]; for (int i = 0; i <
	 * paths.length; ++i) paths[i] = new ArrayList();
	 * 
	 * for (int StartingNode = 0; StartingNode < ScaledPeaks.length;
	 * ++StartingNode) { // Initialize all the datastructures
	 * distances[StartingNode] = 0; paths[StartingNode].clear(); //
	 * scores[StartingNode] = PeakScores[StartingNode]; ArrayList Temp = new
	 * ArrayList(); Temp.add(new Integer(StartingNode));
	 * paths[StartingNode].add(Temp); for (int i = 0; i < distances.length; ++i)
	 * { if (i != StartingNode) { distances[i] = Integer.MIN_VALUE;
	 * paths[i].clear(); } }
	 * 
	 * // We know we can't reach any Node with index < Starting node //
	 * (property of DAG) for (int CurrNode = StartingNode + 1; CurrNode <
	 * distances.length; ++CurrNode) { for (int PrevNode = StartingNode;
	 * PrevNode < CurrNode; ++PrevNode) { if
	 * (adjacencyMatrix[PrevNode][CurrNode]) { for (int i = 0; i <
	 * paths[PrevNode].size(); ++i) { ArrayList CurrList = (ArrayList)
	 * (paths[PrevNode] .get(i)); Temp =
	 * AntibodyUtils.CopyIntArrayList(CurrList); Temp.add(new
	 * Integer(CurrNode)); if (Temp.size() < TagLength)
	 * paths[CurrNode].add(Temp); else { double CurrScore = 0.0; for (int Node =
	 * 0; Node < TagLength; ++Node) { int NodeIndex = ((Integer)
	 * (Temp.get(Node))) .intValue(); CurrScore += PeakScores[NodeIndex]; } if
	 * (CurrScore > TagScore) { BetterScoringTags += 1; } if (CurrScore >=
	 * MinTagScore) TotalTags += 1; if (TotalTags > 1000000) { int[] Ret = {
	 * BetterScoringTags, TotalTags }; return Ret; }
	 * 
	 * } } } } } } int[] Ret = { BetterScoringTags, TotalTags }; return Ret; }
	 */

	public static double[] FindAllPathScores(int PathLength, int[] ScaledPeaks,
			double[] PeakScores, int MaxNum, double Tolerance) {
		/**
		 * This variable represents the graph structure. There is an edge
		 * between two nodes, i and j, in the graph if adjacencyMatrix[i][j] ==
		 * true
		 */
		boolean[][] adjacencyMatrix = BuildAdjacencyMatrixWithDoubles(
				ScaledPeaks, Tolerance);
		int[] distances = new int[adjacencyMatrix.length];
		// double[] scores = new double[adjacencyMatrix.length];

		ArrayList[] paths = new ArrayList[distances.length];
		for (int i = 0; i < paths.length; ++i)
			paths[i] = new ArrayList();

		ArrayList AllPaths = new ArrayList();

		for (int StartingNode = 0; StartingNode < ScaledPeaks.length; ++StartingNode) {
			// Initialize all the datastructures
			distances[StartingNode] = 0;
			paths[StartingNode].clear();
			// scores[StartingNode] = PeakScores[StartingNode];
			ArrayList Temp = new ArrayList();
			Temp.add(new Integer(StartingNode));
			paths[StartingNode].add(Temp);
			for (int i = 0; i < distances.length; ++i) {
				if (i != StartingNode) {
					distances[i] = Integer.MIN_VALUE;
					paths[i].clear();
				}
			}

			// We know we can't reach any Node with index < Starting node
			// (property of DAG)
			for (int CurrNode = StartingNode + 1; CurrNode < distances.length
					&& AllPaths.size() < MaxNum; ++CurrNode) {
				for (int PrevNode = StartingNode; PrevNode < CurrNode
						&& AllPaths.size() < MaxNum; ++PrevNode) {
					if (adjacencyMatrix[PrevNode][CurrNode]) {
						for (int i = 0; i < paths[PrevNode].size()
								&& AllPaths.size() < MaxNum; ++i) {
							ArrayList CurrList = (ArrayList) (paths[PrevNode]
									.get(i));
							if (CurrList.size() < PathLength) {
								Temp = AntibodyUtils.CopyIntArrayList(CurrList);
								Temp.add(new Integer(CurrNode));
								if (Temp.size() < PathLength)
									paths[CurrNode].add(Temp);
								else
									AllPaths.add(Temp);
							}
							if (AllPaths.size() >= MaxNum) {
								break;
							}
						}
					}
				}
			}
		}

		double[] Scores = new double[AllPaths.size()];
		for (int i = 0; i < Scores.length; ++i) {
			ArrayList CurrList = (ArrayList) (AllPaths.get(i));
			double CurrScore = 0.0;
			if (CurrList.size() != PathLength) {
				System.out.println("Path " + i + " of length "
						+ CurrList.size() + "!=" + PathLength);
				continue;
			}
			for (int j = 0; j < CurrList.size(); ++j) {
				CurrScore += PeakScores[((Integer) (CurrList.get(j)))
						.intValue()];
			}
			Scores[i] = CurrScore;
		}

		return Scores;

	}

	public static String ReconstructEnd2End(int[] Peaks, double[] PeakScores,
			double Tolerance) {
		String PeakString = "";
		String ScoreString = "";
		for (int i = 0; i < Peaks.length; ++i) {
			PeakString += Peaks[i] + " ";
			ScoreString += PeakScores[i] + " ";

		}
		// System.out.println("Peaks: " + PeakString);
		// System.out.println("Scores: " + ScoreString);
		/**
		 * This variable represents the graph structure. There is an edge
		 * between two nodes, i and j, in the graph if adjacencyMatrix[i][j] ==
		 * true
		 */
		boolean[][] adjacencyMatrix = BuildAdjacencyMatrixWithDoubles(Peaks,
				Tolerance);

		/**
		 * This variable represents the sequence of nodes in the graph for the
		 * longest path.
		 */
		ArrayList[] longestPaths = FindStrongestPath(adjacencyMatrix,
				PeakScores);

		String finalSeq = "";
		/**
		 * This loop constructs the amino acid sequence from the peaks reported
		 * in the longest path variable.
		 */
		if (Debug)
			System.out.println("Found " + longestPaths.length
					+ " different highest scoring paths");
		double Score = 0.0;
		for (int p = 0; p < longestPaths.length; ++p) {
			finalSeq = "";
			ArrayList CurrList = longestPaths[p];

			for (int i = 1; i < longestPaths[p].size(); ++i) {
				int NodeNum = ((Integer) (CurrList.get(i))).intValue();
				int PrevNode = ((Integer) (CurrList.get(i - 1))).intValue();
				Score += PeakScores[PrevNode];
				if (i == longestPaths[p].size() - 1)
					Score += PeakScores[NodeNum];
				int CurrPeak = Peaks[NodeNum];
				int PrevPeak = Peaks[PrevNode];
				int Diff = CurrPeak - PrevPeak;

				// The Sequence didn't start at 0!!
				if (i == 1 && PrevPeak != 0) {
					finalSeq = "[" + (PrevPeak - 0) / AntibodyUtils.MASS_SCALE
							+ "]";

				}

				// The amino acid added is determined by the difference between
				// the current peak and
				// the previous peak.
				char AA = AntibodyUtils.FindClosestAAScaled(Diff, Tolerance);
				String AAStr = AA + "";
				double AAMass = 0.0;
				// If no single amino acid explains the mass difference, see if
				// pair of amino acids does.
				if (AA == '0') {
					AAMass = Double.parseDouble(AntibodyUtils
							.FindClosestDoubleAAScaled(Diff, Tolerance));
					AAStr = "[" + (AAMass / AntibodyUtils.MASS_SCALE) + "]";
				}
				if (AAStr == "0") {
					ErrorThrower
							.ThrowErrorCustum(
									ErrorThrower.CUSTOM_ERROR,
									"ERROR: Mass "
											+ Diff
											+ " is not explained by either 1 or 2 amino acids!");
				}

				if (Debug) {
					System.out.println("Peak " + PrevPeak + " to " + CurrPeak);
					System.out.println("For Diff: " + Diff + " AA: " + AA);
				}

				finalSeq += AAStr;

				// The sequence didn't end at the PM
				if (i == longestPaths[p].size() - 1
						&& CurrPeak != Peaks[Peaks.length - 1]) {
					finalSeq += "[" + (Peaks[Peaks.length - 1] - CurrPeak)
							/ AntibodyUtils.MASS_SCALE + "]";
				}
			}

			// if(Debug || !Debug)
			// System.out.println("Final Seq: " + finalSeq);
		}
		// System.out.println("DENOVO Score: " + Score);
		return finalSeq;

	}

	/**
	 * This is the gateway method for finding the highest scoring reconstruction
	 * of the peaks,
	 * 
	 * @param Peaks
	 *            The list of scaled peak masses (PRMs)
	 * @param PeakScores
	 *            The scores for the PRMs in Peaks
	 * @return Returns the string reconstruction. A jump of 2 amino acids
	 *         appears as the mass
	 */
	public static String Reconstruct(int[] Peaks, double[] PeakScores,
			double Tolerance) {

		String PeakString = "";
		String ScoreString = "";
		for (int i = 0; i < Peaks.length; ++i) {
			PeakString += Peaks[i] + " ";
			ScoreString += PeakScores[i] + " ";

		}
		// System.out.println("Peaks: " + PeakString);
		// System.out.println("Scores: " + ScoreString);
		/**
		 * This variable represents the graph structure. There is an edge
		 * between two nodes, i and j, in the graph if adjacencyMatrix[i][j] ==
		 * true
		 */
		boolean[][] adjacencyMatrix = BuildAdjacencyMatrixWithDoubles(Peaks,
				Tolerance);

		/**
		 * This variable represents the sequence of nodes in the graph for the
		 * longest path.
		 */
		ArrayList[] longestPaths = FindStrongestPath(adjacencyMatrix,
				PeakScores);

		String finalSeq = "";
		/**
		 * This loop constructs the amino acid sequence from the peaks reported
		 * in the longest path variable.
		 */
		if (Debug) {
			System.out.println("Found " + longestPaths.length
					+ " different highest scoring paths");
			System.out.println("Len: " + longestPaths[0].size());
		}
		double Score = 0.0;
		for (int p = 0; p < longestPaths.length; ++p) {
			finalSeq = "";
			ArrayList CurrList = longestPaths[p];

			for (int i = 1; i < longestPaths[p].size(); ++i) {
				int NodeNum = ((Integer) (CurrList.get(i))).intValue();
				int PrevNode = ((Integer) (CurrList.get(i - 1))).intValue();
				Score += PeakScores[PrevNode];
				if (i == longestPaths[p].size() - 1)
					Score += PeakScores[NodeNum];
				int CurrPeak = Peaks[NodeNum];
				int PrevPeak = Peaks[PrevNode];
				int Diff = CurrPeak - PrevPeak;

				// The amino acid added is determined by the difference between
				// the current peak and
				// the previous peak.
				char AA = AntibodyUtils.FindClosestAAScaled(Diff, Tolerance);
				String AAStr = AA + "";
				double AAMass = 0.0;
				// If no single amino acid explains the mass difference, see if
				// pair of amino acids does.
				if (AA == '0') {
					AAMass = Double.parseDouble(AntibodyUtils
							.FindClosestDoubleAAScaled(Diff, Tolerance));
					AAStr = "[" + (AAMass / AntibodyUtils.MASS_SCALE) + "]";
				}
				if (AAStr == "0") {
					ErrorThrower
							.ThrowErrorCustum(
									ErrorThrower.CUSTOM_ERROR,
									"ERROR: Mass "
											+ Diff
											+ " is not explained by either 1 or 2 amino acids!");
				}

				if (Debug) {
					System.out.println("Peak " + PrevPeak + " to " + CurrPeak);
					System.out.println("For Diff: " + Diff + " AA: " + AA);
				}

				finalSeq += AAStr;
			}

			if (Debug)
				System.out.println("Final Seq: " + finalSeq);
		}
		// System.out.println("DENOVO Score: " + Score);
		return finalSeq;
	}

	/**
	 * This is the gateway method for finding the highest scoring reconstruction
	 * of the peaks,
	 * 
	 * @param Peaks
	 *            The list of scaled peak masses (PRMs)
	 * @param PeakScores
	 *            The scores for the PRMs in Peaks
	 * @return Returns the string reconstruction. A jump of 2 amino acids
	 *         appears as the mass
	 */
	public static int[] ReconstructPRMs(int[] Peaks, double[] PeakScores,
			double Tolerance) {

		String PeakString = "";
		String ScoreString = "";
		for (int i = 0; i < Peaks.length; ++i) {
			PeakString += Peaks[i] + " ";
			ScoreString += PeakScores[i] + " ";

		}
		// System.out.println("Peaks: " + PeakString);
		// System.out.println("Scores: " + ScoreString);
		/**
		 * This variable represents the graph structure. There is an edge
		 * between two nodes, i and j, in the graph if adjacencyMatrix[i][j] ==
		 * true
		 */
		boolean[][] adjacencyMatrix = BuildAdjacencyMatrixWithDoubles(Peaks,
				Tolerance);

		/**
		 * This variable represents the sequence of nodes in the graph for the
		 * longest path.
		 */
		ArrayList[] longestPaths = FindStrongestPath(adjacencyMatrix,
				PeakScores);

		String finalSeq = "";
		int[] Ret = null;
		/**
		 * This loop constructs the amino acid sequence from the peaks reported
		 * in the longest path variable.
		 */
		if (Debug)
			System.out.println("Found " + longestPaths.length
					+ " different highest scoring paths");
		double Score = 0.0;
		for (int p = 0; p < longestPaths.length; ++p) {
			finalSeq = "";
			ArrayList CurrList = longestPaths[p];
			Ret = new int[CurrList.size()];
			for (int i = 1; i < longestPaths[p].size(); ++i) {
				int NodeNum = ((Integer) (CurrList.get(i))).intValue();
				int PrevNode = ((Integer) (CurrList.get(i - 1))).intValue();

				Score += PeakScores[PrevNode];
				if (i == longestPaths[p].size() - 1)
					Score += PeakScores[NodeNum];
				int CurrPeak = Peaks[NodeNum];
				int PrevPeak = Peaks[PrevNode];
				Ret[i - 1] = PrevPeak;
				Ret[i] = CurrPeak;
				int Diff = CurrPeak - PrevPeak;

				// The amino acid added is determined by the difference between
				// the current peak and
				// the previous peak.
				char AA = AntibodyUtils.FindClosestAAScaled(Diff, Tolerance);
				String AAStr = AA + "";
				double AAMass = 0.0;
				// If no single amino acid explains the mass difference, see if
				// pair of amino acids does.
				if (AA == '0') {
					AAMass = Double.parseDouble(AntibodyUtils
							.FindClosestDoubleAAScaled(Diff, Tolerance));
					AAStr = "[" + (AAMass / AntibodyUtils.MASS_SCALE) + "]";
				}
				if (AAStr == "0") {
					ErrorThrower
							.ThrowErrorCustum(
									ErrorThrower.CUSTOM_ERROR,
									"ERROR: Mass "
											+ Diff
											+ " is not explained by either 1 or 2 amino acids!");
				}

				if (Debug) {
					System.out.println("Peak " + PrevPeak + " to " + CurrPeak);
					System.out.println("For Diff: " + Diff + " AA: " + AA);
				}

				finalSeq += AAStr;
			}

			// if(Debug || !Debug)
			// System.out.println("Final Seq: " + finalSeq);
		}
		// System.out.println("DENOVO Score: " + Score);
		return Ret;
	}

	public static double[][] ReconstructPRMsAndScores(int[] Peaks,
			double[] PeakScores, double Tolerance) {

		String PeakString = "";
		String ScoreString = "";
		for (int i = 0; i < Peaks.length; ++i) {
			PeakString += Peaks[i] + " ";
			ScoreString += PeakScores[i] + " ";

		}
		// System.out.println("Peaks: " + PeakString);
		// System.out.println("Scores: " + ScoreString);
		/**
		 * This variable represents the graph structure. There is an edge
		 * between two nodes, i and j, in the graph if adjacencyMatrix[i][j] ==
		 * true
		 */
		boolean[][] adjacencyMatrix = BuildAdjacencyMatrixWithDoubles(Peaks,
				Tolerance);

		/**
		 * This variable represents the sequence of nodes in the graph for the
		 * longest path.
		 */
		ArrayList[] longestPaths = FindStrongestPath(adjacencyMatrix,
				PeakScores);

		String finalSeq = "";
		double[][] Ret = null;
		/**
		 * This loop constructs the amino acid sequence from the peaks reported
		 * in the longest path variable.
		 */
		if (Debug)
			System.out.println("Found " + longestPaths.length
					+ " different highest scoring paths");
		double Score = 0.0;
		for (int p = 0; p < longestPaths.length; ++p) {
			finalSeq = "";
			ArrayList CurrList = longestPaths[p];
			Ret = new double[CurrList.size()][2];
			for (int i = 1; i < longestPaths[p].size(); ++i) {
				int NodeNum = ((Integer) (CurrList.get(i))).intValue();
				int PrevNode = ((Integer) (CurrList.get(i - 1))).intValue();

				Score += PeakScores[PrevNode];
				if (i == longestPaths[p].size() - 1)
					Score += PeakScores[NodeNum];
				int CurrPeak = Peaks[NodeNum];
				int PrevPeak = Peaks[PrevNode];
				Ret[i - 1][0] = PrevPeak;
				Ret[i][0] = CurrPeak;

				Ret[i - 1][1] = PeakScores[PrevNode];
				Ret[i][1] = PeakScores[NodeNum];

				int Diff = CurrPeak - PrevPeak;

				// The amino acid added is determined by the difference between
				// the current peak and
				// the previous peak.
				char AA = AntibodyUtils.FindClosestAAScaled(Diff, Tolerance);
				String AAStr = AA + "";
				double AAMass = 0.0;
				// If no single amino acid explains the mass difference, see if
				// pair of amino acids does.
				if (AA == '0') {
					AAMass = Double.parseDouble(AntibodyUtils
							.FindClosestDoubleAAScaled(Diff, Tolerance));
					AAStr = "[" + (AAMass / AntibodyUtils.MASS_SCALE) + "]";
				}
				if (AAStr == "0") {
					ErrorThrower
							.ThrowErrorCustum(
									ErrorThrower.CUSTOM_ERROR,
									"ERROR: Mass "
											+ Diff
											+ " is not explained by either 1 or 2 amino acids!");
				}

				if (Debug) {
					System.out.println("Peak " + PrevPeak + " to " + CurrPeak);
					System.out.println("For Diff: " + Diff + " AA: " + AA);
				}

				finalSeq += AAStr;
			}

			// if(Debug || !Debug)
			// System.out.println("Final Seq: " + finalSeq);
		}
		// System.out.println("DENOVO Score: " + Score);
		return Ret;
	}

	/**
	 * This is the gateway method to be called by anyone who wants a sequence
	 * for the list of PRM peaks that they are interested in. This interface
	 * should not change without talking to Natalie about it.
	 * 
	 * @param Peaks
	 *            The list of peaks (all scalled by Utils.MASS_SCALE)
	 * @return The best amino acid sequence constructed from these peaks. If a
	 *         jump of 2 amino acids occurs in the graph, then a mass is output
	 *         instead of an amino acid
	 */
	public static String Reconstruct(int[] Peaks, double Tolerance) {
		/**
		 * This variable represents the graph structure. There is an edge
		 * between two nodes, i and j, in the graph if adjacencyMatrix[i][j] ==
		 * true
		 */
		boolean[][] adjacencyMatrix = BuildAdjacencyMatrixWithDoubles(Peaks,
				Tolerance);
		// boolean[][] adjacencyMatrix = BuildAdjacencyMatrix(Peaks);

		/**
		 * This variable represents the sequence of nodes in the graph for the
		 * longest path.
		 */
		int[][] longestPaths = FindLongestPath(adjacencyMatrix);

		String finalSeq = "";
		/**
		 * This loop constructs the amino acid sequence from the peaks reported
		 * in the longest path variable.
		 */
		if (Debug)
			System.out.println("Found " + longestPaths.length
					+ " different longest paths");
		for (int p = 0; p < longestPaths.length; ++p) {
			finalSeq = "";
			for (int i = 1; i < longestPaths[p].length; ++i) {
				int CurrPeak = Peaks[longestPaths[p][i]];
				int PrevPeak = Peaks[longestPaths[p][i - 1]];
				int Diff = CurrPeak - PrevPeak;

				// The amino acid added is determined by the difference between
				// the current peak and
				// the previous peak.
				char AA = AntibodyUtils.FindClosestAAScaled(Diff, Tolerance);
				String AAStr = AA + "";

				// If no single amino acid explains the mass difference, see if
				// pair of amino acids does.
				if (AA == '0') {
					double AAMass = Double.parseDouble(AntibodyUtils
							.FindClosestDoubleAAScaled(Diff, Tolerance));
					AAStr = "[" + (AAMass / AntibodyUtils.MASS_SCALE) + "]";
				}
				if (AAStr == "0") {
					ErrorThrower
							.ThrowErrorCustum(
									ErrorThrower.CUSTOM_ERROR,
									"ERROR: Mass "
											+ Diff
											+ " is not explained by either 1 or 2 amino acids!");
				}

				if (Debug) {
					System.out.println("Peak " + PrevPeak + " to " + CurrPeak);
					System.out.println("For Diff: " + Diff + " AA: " + AA);
				}

				finalSeq += AAStr;
			}

			if (Debug)
				System.out.println("Final Seq: " + finalSeq);
		}
		System.out.println("DENOVO Length: " + longestPaths[0].length);
		return finalSeq;
	}

	public static int[] ReconstructPRMs(int[] Peaks, double Tolerance) {
		/**
		 * This variable represents the graph structure. There is an edge
		 * between two nodes, i and j, in the graph if adjacencyMatrix[i][j] ==
		 * true
		 */
		boolean[][] adjacencyMatrix = BuildAdjacencyMatrixWithDoubles(Peaks,
				Tolerance);
		// boolean[][] adjacencyMatrix = BuildAdjacencyMatrix(Peaks);

		/**
		 * This variable represents the sequence of nodes in the graph for the
		 * longest path.
		 */
		int[][] longestPaths = FindLongestPath(adjacencyMatrix);

		String finalSeq = "";
		int[] Ret = null;
		/**
		 * This loop constructs the amino acid sequence from the peaks reported
		 * in the longest path variable.
		 */
		if (Debug)
			System.out.println("Found " + longestPaths.length
					+ " different longest paths");
		for (int p = 0; p < longestPaths.length; ++p) {
			finalSeq = "";
			Ret = new int[longestPaths[p].length];
			for (int i = 1; i < longestPaths[p].length; ++i) {
				int CurrPeak = Peaks[longestPaths[p][i]];
				int PrevPeak = Peaks[longestPaths[p][i - 1]];
				Ret[i - 1] = PrevPeak;
				Ret[i] = CurrPeak;
				int Diff = CurrPeak - PrevPeak;

				// The amino acid added is determined by the difference between
				// the current peak and
				// the previous peak.
				char AA = AntibodyUtils.FindClosestAAScaled(Diff, Tolerance);
				String AAStr = AA + "";

				// If no single amino acid explains the mass difference, see if
				// pair of amino acids does.
				if (AA == '0') {
					double AAMass = Double.parseDouble(AntibodyUtils
							.FindClosestDoubleAAScaled(Diff, Tolerance));
					AAStr = "[" + (AAMass / AntibodyUtils.MASS_SCALE) + "]";
				}
				if (AAStr == "0") {
					ErrorThrower
							.ThrowErrorCustum(
									ErrorThrower.CUSTOM_ERROR,
									"ERROR: Mass "
											+ Diff
											+ " is not explained by either 1 or 2 amino acids!");
				}

				if (Debug) {
					System.out.println("Peak " + PrevPeak + " to " + CurrPeak);
					System.out.println("For Diff: " + Diff + " AA: " + AA);
				}

				finalSeq += AAStr;
			}

			if (Debug)
				System.out.println("Final Seq: " + finalSeq);
		}
		System.out.println("DENOVO Length: " + longestPaths[0].length);
		return Ret;
	}

	/**
	 * Builds the adjacencyMatrix. If the difference between two peaks is
	 * explained by a single amino acid, then the nodes are connected. If
	 * adjacencyMatrix[i][j] == true, then peaks i and j are connected.
	 * 
	 * @param Peaks
	 *            The list of unscaled peaks to construct the spectrum graph
	 *            from
	 * @return Returns the adjacency matrix for the spectrum graph constructed
	 *         from Peaks.
	 */
	public static boolean[][] BuildAdjacencyMatrix(int[] Peaks, double Tolerance) {
		boolean[][] adjacencyMatrix = new boolean[Peaks.length][Peaks.length];

		for (int i = 0; i < Peaks.length; ++i) {
			for (int j = i + 1; j < Peaks.length; ++j) {
				int Diff = Peaks[j] - Peaks[i];
				if (AntibodyUtils.FindClosestAAScaled(Diff, Tolerance) != '0') {
					if (Debug) {
						System.out.println("Peak " + Peaks[i] + "(" + i
								+ ") has child " + Peaks[j] + "(" + j + ")");
						System.out.println("AA Diff = "
								+ AntibodyUtils.FindClosestAAScaled(Diff,
										Tolerance));
					}

					adjacencyMatrix[i][j] = true;
				}
			}
		}
		return adjacencyMatrix;
	}

	/**
	 * Builds the adjacency matrix. If the difference between two peaks is
	 * explained by a single amino acid or a pair of amino acids, then the nodes
	 * are connected. If adjacencyMatrix[i][j] == true, then peaks i and j are
	 * connected.
	 * 
	 * @param Peaks
	 *            The list of unscaled peaks from which to construct the
	 *            spectrum graph.
	 * @return Returns the adjacency matrix for the spectrum graph constructed
	 *         from Peaks.
	 */
	public static boolean[][] BuildAdjacencyMatrixWithDoubles(int[] Peaks,
			double Tolerance) {
		boolean[][] adjacencyMatrix = new boolean[Peaks.length][Peaks.length];

		for (int i = 0; i < Peaks.length; ++i) {
			for (int j = i + 1; j < Peaks.length; ++j) {
				int Diff = Peaks[j] - Peaks[i];
				if (AntibodyUtils.FindClosestAAScaled(Diff, Tolerance) != '0'
						|| AntibodyUtils.FindClosestDoubleAAScaled(Diff,
								Tolerance).compareTo("0") != 0) {

					if (Debug) {
						System.out.println("Peak " + Peaks[i] + "(" + i
								+ ") has child " + Peaks[j] + "(" + j + ")");
						System.out.println("AA Diff = "
								+ AntibodyUtils.FindClosestAAScaled(Diff,
										Tolerance)
								+ " Double AA Diff = "
								+ AntibodyUtils.FindClosestDoubleAAScaled(Diff,
										Tolerance));
					}

					adjacencyMatrix[i][j] = true;
				}
			}
		}

		return adjacencyMatrix;
	}

	public static int[][] FindLongestPath(boolean[][] adjacencyMatrix) {

		// distances[i] is the length of the longest path so far from node 0 to
		// i.
		int[] distances = new int[adjacencyMatrix.length];
		ArrayList[] paths = new ArrayList[distances.length];
		distances[0] = 0;
		paths[0] = new ArrayList();
		ArrayList Temp = new ArrayList();
		Temp.add(new Integer(0));
		paths[0].add(Temp);
		for (int i = 1; i < distances.length; ++i) {
			distances[i] = Integer.MIN_VALUE;
			paths[i] = new ArrayList();
		}

		for (int CurrNode = 0; CurrNode < distances.length; ++CurrNode) {
			for (int PrevNode = 0; PrevNode < CurrNode; ++PrevNode) {
				if (adjacencyMatrix[PrevNode][CurrNode]) {
					if (distances[PrevNode] + 1 > distances[CurrNode]) {
						distances[CurrNode] = distances[PrevNode] + 1;
						paths[CurrNode] = new ArrayList();
						for (int i = 0; i < paths[PrevNode].size(); ++i) {
							Temp = AntibodyUtils
									.CopyIntArrayList((ArrayList) (paths[PrevNode]
											.get(i)));
							Temp.add(new Integer(CurrNode));
							paths[CurrNode].add(Temp);
						}
					}

				}

			}

		}

		/**
		 * Once we have computed the longest paths to each node, determine the
		 * distance of the longest paths overall. maxNode is only one of several
		 * nodes that could have a longest path from node 0.
		 */
		int maxDist = 0;
		int maxNode = 0;
		for (int i = 0; i < distances.length; ++i) {
			if (distances[i] >= maxDist) {
				maxDist = distances[i];
				maxNode = i;
			}
		}

		if (Debug && !Debug) {
			System.out.println("Distance to " + (distances.length - 1) + " = "
					+ distances[distances.length - 1]);

			System.out.println("Longest path ends at node " + maxNode
					+ " at distance " + maxDist);

		}
		int[][] finalPaths = new int[paths[maxNode].size()][maxDist + 1];
		/**
		 * Determine the final paths, by converting the path[maxNode]. Keep in
		 * mind that this is only 1 of several possible longest paths from node
		 * 0 to maxNode, and only 1 of several nodes which may have a path from
		 * node 0 of maxDist.
		 */
		for (int p = 0; p < paths[maxNode].size(); ++p) {
			ArrayList CurrPath = (ArrayList) (paths[maxNode].get(p));
			for (int i = 0; i < maxDist; ++i) {
				finalPaths[p][i] = ((Integer) (CurrPath.get(i))).intValue();

			}
			finalPaths[p][maxDist] = maxNode;
		}
		return finalPaths;

	}

	public static ArrayList[] FindStrongestPath(boolean[][] adjacencyMatrix,
			double[] NodeScores) {

		// scores[i] is the score of the strongest path so far from node 0 to i.
		double[] scores = new double[adjacencyMatrix.length];
		ArrayList[] paths = new ArrayList[scores.length];
		scores[0] = NodeScores[0];
		paths[0] = new ArrayList();
		ArrayList Temp = new ArrayList();
		Temp.add(new Integer(0));
		paths[0].add(Temp);
		for (int i = 1; i < scores.length; ++i) {
			scores[i] = NodeScores[i];
			paths[i] = new ArrayList();
			Temp = new ArrayList();
			Temp.add(new Integer(i));
			paths[i].add(Temp);
		}

		for (int CurrNode = 1; CurrNode < scores.length; ++CurrNode) {
			for (int PrevNode = 0; PrevNode < CurrNode; ++PrevNode) {
				if (adjacencyMatrix[PrevNode][CurrNode]) {
					if (scores[PrevNode] + NodeScores[CurrNode] > scores[CurrNode]) {
						scores[CurrNode] = scores[PrevNode]
								+ NodeScores[CurrNode];
						paths[CurrNode] = new ArrayList();
						// distances[CurrNode] = distances[PrevNode] + 1;
						for (int i = 0; i < paths[PrevNode].size(); ++i) {
							Temp = AntibodyUtils
									.CopyIntArrayList((ArrayList) (paths[PrevNode]
											.get(i)));
							Temp.add(new Integer(CurrNode));
							paths[CurrNode].add(Temp);
						}
					}

				}

			}

		}

		/**
		 * Once we have computed the longest paths to each node, determine the
		 * score of the strongest paths overall. maxNode is only one of several
		 * nodes that could have a strongest path from node 0.
		 */
		double maxScore = Double.NEGATIVE_INFINITY;
		int maxNode = 0;
		// int maxDist = 0;
		for (int i = 0; i < scores.length; ++i) {
			if (scores[i] >= maxScore) {
				// maxDist = distances[i];
				maxScore = scores[i];
				maxNode = i;
			}
		}

		if (Debug && !Debug) {
			System.out.println("Distance to " + (scores.length - 1) + " = "
					+ scores[scores.length - 1]);

			System.out.println("Strongest path ends at node " + maxNode
					+ " at score " + maxScore);

		}
		ArrayList[] finalPaths = new ArrayList[paths[maxNode].size()];
		/**
		 * Determine the final paths, by converting the path[maxNode]. Keep in
		 * mind that this is only 1 of several possible longest paths from node
		 * 0 to maxNode, and only 1 of several nodes which may have a path from
		 * node 0 of maxDist.
		 */
		for (int p = 0; p < paths[maxNode].size(); ++p) {
			finalPaths[p] = (ArrayList) (paths[maxNode].get(p));

			// finalPaths[p].add(new Integer(maxNode));
		}
		return finalPaths;

	}

	/**
	 * Determines the path of nodes representing a longest path from Node 0.
	 * Currently it only records 1 longest path for each end node, but in
	 * reality there may be multiple paths. The method works much like
	 * Dijkstra's for finding the longest path. A depth first search is used to
	 * continuously update a table of path lengths to each node. TO DO: 1.
	 * Reconsider scoring function beyond simply the longest path 2. Record
	 * multiple paths for each end node
	 * 
	 * @param adjacencyMatrix
	 *            The matrix representing the spectrum graph
	 * @return Returns the sequence of nodes visited by following a longest path
	 *         in the graph.
	 */
	public static int[][] FindLongestPathOld(boolean[][] adjacencyMatrix) {
		// Keeps track of the distances to each node, i.e. distances[i] is the
		// length of the longest path so far from node 0 to node i.
		int[] distances = new int[adjacencyMatrix.length];

		// Keeps track of the longest paths, i.e. paths[i] is all longest paths
		// (of distance[i]) from node 0 to node i.
		ArrayList[] paths = new ArrayList[adjacencyMatrix.length];

		// Initialize the distances to 0 and the paths.
		for (int i = 0; i < adjacencyMatrix.length; ++i) {
			distances[i] = 0;
			paths[i] = new ArrayList();
		}

		/**
		 * Starting from node 0 we perform a depth first search, and the queue
		 * maintains a list of nodes in the current depth. So initially, the
		 * queue only holds node 0, then we remove node 0 and add all children
		 * of node 0 (all nodes at depth 1 from node 0).
		 */
		ArrayList queue = new ArrayList();
		queue.add(new Integer(0));

		/**
		 * Continue, considering all nodes at the current depth until we reach
		 * the leaf nodes (leaving the queue empty). We are guaranteed to reach
		 * an end because the graph is directed, acyclic.
		 */
		while (queue.size() > 0) {
			// Remove this guy
			int Node = ((Integer) (queue.remove(0))).intValue();
			if (Debug) {
				System.out.println("Considering Node: " + Node);
			}
			for (int i = Node + 1; i < adjacencyMatrix.length; ++i) {
				// Add the children to the queue, and update the distance to
				// them through the current Node.
				if (adjacencyMatrix[Node][i]) {
					if (Debug)
						System.out.println("Node " + i + " is a child of "
								+ Node);
					queue.add(new Integer(i));
					if (distances[i] < distances[Node] + 1) {
						distances[i] = distances[Node] + 1;
						paths[i].clear();
						for (int p = 0; p < paths[Node].size(); ++p) {
							ArrayList NewPath = (AntibodyUtils
									.CopyIntArrayList((ArrayList) (paths[Node]
											.get(p))));
							NewPath.add(new Integer(Node));
							if (Debug)
								System.out.println("New path: " + NewPath);
							if (!AntibodyUtils.IntArrayListContainedIn(
									paths[i], NewPath))
								paths[i].add(NewPath);
						}
						if (paths[Node].size() == 0) {
							ArrayList NewPath = new ArrayList();
							NewPath.add(new Integer(Node));
							if (!AntibodyUtils.IntArrayListContainedIn(
									paths[i], NewPath))
								paths[i].add(NewPath);
						}
					} else if (distances[i] == distances[Node] + 1) {
						// distances[i] = distances[Node] + 1;
						for (int p = 0; p < paths[Node].size(); ++p) {
							ArrayList NewPath = (AntibodyUtils
									.CopyIntArrayList((ArrayList) (paths[Node]
											.get(p))));
							NewPath.add(new Integer(Node));
							// paths[i] = new ArrayList();
							if (!AntibodyUtils.IntArrayListContainedIn(
									paths[i], NewPath))
								paths[i].add(NewPath);
						}

					}
					if (Debug) {
						for (int k = 0; k < distances.length; ++k) {
							System.out.println(k + ":" + distances[k]);
							System.out.println(paths[k]);
						}
					}
				}
			}
		}

		/**
		 * Once we have computed the longest paths to each node, determine the
		 * distance of the longest paths overall. maxNode is only one of several
		 * nodes that could have a longest path from node 0.
		 */
		int maxDist = 0;
		int maxNode = 0;
		for (int i = 0; i < distances.length; ++i) {
			if (distances[i] >= maxDist) {
				maxDist = distances[i];
				maxNode = i;
			}
		}

		if (Debug || !Debug) {
			System.out.println("Distance to " + (distances.length - 1) + " = "
					+ distances[distances.length - 1]);

			System.out.println("Longest path ends at node " + maxNode
					+ " at distance " + maxDist);

		}
		int[][] finalPaths = new int[paths[maxNode].size()][maxDist + 1];
		/**
		 * Determine the final paths, by converting the path[maxNode]. Keep in
		 * mind that this is only 1 of several possible longest paths from node
		 * 0 to maxNode, and only 1 of several nodes which may have a path from
		 * node 0 of maxDist.
		 */
		for (int p = 0; p < paths[maxNode].size(); ++p) {
			ArrayList CurrPath = (ArrayList) (paths[maxNode].get(p));
			for (int i = 0; i < maxDist; ++i) {
				finalPaths[p][i] = ((Integer) (CurrPath.get(i))).intValue();

			}
			finalPaths[p][maxDist] = maxNode;
		}
		return finalPaths;

	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Each example has the peaks for the consensus after 10,20,30,... runs
		// and we predict the sequence
		// for each intermediate consensus
		// double[] Peaks1 = {0.0,57.03,114.06};

		// EXAMPLE 1
		// FSLTSYG - VSWVRQP - PGKGLE
		// 10 iterations
		/*
		 * double[] Peaks1 = {0.079, 147.11, 234.137, 347.2, 448.329, 535.325,
		 * 698.396, 755.367}; int [] Peaks = Utils.ScaleArray(Peaks1);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 * 
		 * //20 iterations double[] Peaks2 = {0.079, 147.11, 347.192, 448.293,
		 * 535.303, 698.378, 755.318, 782.674, 791.8705, 849.539, 854.422,
		 * 873.783666667, 941.438, 1012.643, 1099.64116667, 1127.53533333,
		 * 1226.5725, 1382.74383333}; Peaks = Utils.ScaleArray(Peaks2);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 * 
		 * //30 iterations double[] Peaks3 = {0.089, 147.106, 347.189, 448.302,
		 * 535.303, 698.374, 755.33, 782.672, 791.892, 849.513, 854.418,
		 * 873.791, 885.517, 941.452, 1012.657, 1030.527, 1099.62, 1127.507,
		 * 1226.538, 1382.744}; Peaks = Utils.ScaleArray(Peaks3);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 */
		// EXAMPLE 2
		// GSTNYHSAL - ISRLSISK
		/*
		 * //10 iterations double[] Peaks1 = {0.17, 57.021, 144.087, 245.108,
		 * 359.179, 522.242, 659.31, 746.324, 817.408, 930.415, 1025.4635,
		 * 1027.608, 1043.548, 1130.664, 1286.616}; int [] Peaks =
		 * Utils.ScaleArray(Peaks1);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 * 
		 * //20 iterations double[] Peaks2 = {0.17, 56.716, 144.087, 245.108,
		 * 359.187, 522.231, 659.303, 746.347, 817.403, 930.434, 1004.63133333,
		 * 1025.492, 1027.608, 1043.593, 1060.496, 1130.652, 1286.656}; Peaks =
		 * Utils.ScaleArray(Peaks2);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 * 
		 * //30 iterations double[] Peaks3 = {0.17, 56.716, 144.087, 245.108,
		 * 359.183, 522.244, 659.288, 746.331, 817.394, 930.471, 934.326333333,
		 * 948.271, 987.544, 1004.602, 1009.52, 1015.733, 1017.48733333,
		 * 1025.466, 1027.532, 1043.571, 1057.44866667, 1060.496, 1061.56133333,
		 * 1112.626, 1130.665, 1286.668, 1399.762}; Peaks =
		 * Utils.ScaleArray(Peaks3);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 */
		// EXAMPLE 3
		// SQVFLKL -
		// 10 iterations
		/*
		 * double[] Peaks1 = {0.017, 87.032, 215.09, 314.175, 461.249, 574.34,
		 * 702.416, 815.488}; int [] Peaks = Utils.ScaleArray(Peaks1);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 * 
		 * //20 iterations double[] Peaks2 = {0.024, 87.393, 215.09, 314.174,
		 * 461.215, 574.345, 702.368, 815.38}; Peaks = Utils.ScaleArray(Peaks2);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 * 
		 * //29 iterations double[] Peaks3 = {0.024, 87.393, 215.09, 314.174,
		 * 461.209, 574.335, 702.373, 815.43, 819.55, 825.623, 835.8, 929.542,
		 * 938.5225, 1016.533, 1042.8695, 1051.6315, 1055.272, 1083.218,
		 * 1100.0385, 1100.625, 1109.328, 1129.767, 1138.584, 1144.5205,
		 * 1147.7065, 1155.631, 1183.2925, 1240.792, 1245.541, 1247.7065,
		 * 1252.563, 1257.77, 1262.737, 1275.711, 1300.5995, 1318.288,
		 * 1347.7725, 1358.909, 1363.705, 1474.158, 1588.88, 1632.358}; Peaks =
		 * Utils.ScaleArray(Peaks3);
		 * System.out.println(DeNovo.Reconstruct(Peaks));
		 */
		int [] Peaks = {137058,200160,243100,261096,362157,390152,561251,720268,741302,776331,791305,833352,905400,987416};
		System.out.println(DeNovo.Reconstruct(Peaks, 0.3));

		
		//PRMSpectrum[] prms = PRMSpectrum.LoadAllPRMSpectrum("/Users/natalie/Documents/DigitalProteomics/Customers/RinatPfizer/140003/RMP1_chymo_MSMS_678_high_res.prms");
		/*PRMSpectrum[] prms = PRMSpectrum.LoadAllPRMSpectrum("/Users/natalie/Documents/DigitalProteomics/Customers/RinatPfizer/140003/CDR3/CIDETD_specs_scored.mgf");
		int matchMass = (int)(Utils.GetMFromMZ(678.79, 2)*AntibodyUtils.MASS_SCALE);
		System.out.println("Mz: 678.79, M: " + ((double)matchMass)/1000);
		for(int i = 0; i < prms.length; ++i) {
			if(Utils.isSamePM(prms[i].PrecursorMass,matchMass,40))
					System.out.println("[" + i + "]: (" + prms[i].PrecursorMass + "), " + DeNovo.Reconstruct(prms[i].PRMs, prms[i].PRMScores, 0.5));
		}
		*/
		//double[] peaks = {0,262.12,363.167,521.237,554.107,695.276,760.239,766.276,774.243,778.321,891.332,992.328,1093.418,1337.528};
		//int[] peaksScaled = AntibodyUtils.ScaleArray(peaks);
		//System.out.println(DeNovo.Reconstruct(peaksScaled, 0.1));
		
		// EXAMPLE 4
		// TATYYCAK - GGYRF (VD region)
		// 30 iterations
		//double[] Peaks1 = {0, 186.505, 189.139, 189.546, 195.414, 198.165, 202.507, 207.163, 209.632, 212.386, 214.044, 217.133, 218.072, 219.066, 224.447, 228.43, 233.165, 234.103, 234.816, 241.193, 244.211, 244.605, 245.689, 252.798, 254.662, 256.82, 258.6, 259.055, 263.139, 264.142, 270.396, 273.471, 297.11, 364.186, 365.19, 388.175, 395.971, 455.821, 479.679, 486.235, 488.192, 488.694, 497.199, 504.246, 513.213, 514.22, 516.707, 521.197, 521.697, 522.252, 524.721, 525.219, 529.71, 530.204, 530.704, 533.723, 538.716, 539.216, 539.709, 547.721, 548.222, 565.229, 574.235, 574.716, 578.237, 578.736, 587.242, 587.748, 590.235, 596.298, 642.536, 647.769, 648.27, 652.265, 652.765, 653.3, 654.3, 656.774, 657.275, 660.776, 661.277, 661.77, 665.259, 669.781, 670.292, 695.042, 704.278, 705.28, 751.073, 751.408, 751.741, 752.076, 762.156, 765.471, 767.339, 777.646, 777.896, 778.145, 779.404, 779.655, 779.906, 800.767, 802.162, 802.414, 802.664, 802.918, 808.423, 809.429, 835.318, 836.321, 853.224, 892.338, 893.342, 976.377, 993.387, 994.39, 1051.408, 1076.424, 1077.418, 1078.41, 1094.434, 1095.438, 1267.578, 1267.945, 1268.256};
		/*double[] Peaks1 = { 0.020, 100.90, 172.060, 273.150, 436.160, 599.260,
				759.250, 830.270, 958.360, 977.490, 988.16, 997.21, 999.59,
				1015.39, 1041.61, 1070.25, 1072.42, 1077.69, 1078.68, 1086.08,
				1086.58, 1090.52, 1122.72, 1123.37, 1151.31, 1161.49, 1162.57,
				1193.21, 1196.62, 1197.45, 1213.68, 1230.61, 1233.61, 1235.52,
				1250.55, 1253.52, 1325.51, 1375.9, 1391.6, 1409.9, 1538.72 };*/
		
		
		

	}

}
