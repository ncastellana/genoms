package antibody;

import java.util.ArrayList;
import java.util.Hashtable;

import basicUtils.Utils;
import errorUtils.ErrorThrower;

public abstract class SpectrumHMM {

	public boolean Debug = false;
	// Seed sequence for the model
	protected String Seed = "";

	// The number of spectra overlapping each match state, OverlapCount[i] =
	// [#of spectra supporting match state i,# of spectra overlapping match
	// state i]
	// This is volatile, and changes at each alignment.
	protected int[][] OverlapCount;
	protected ArrayList[] MatchedPeaks;

	protected int[] InsertOverlapCount;

	// Alignments of spectra that are used, the MassShift parameter should be
	// set for these spectra
	protected ArrayList AlignedSpectra; // Should be of type HMMAlignment

	protected boolean OracleFlag = false;

	// The scaled prefix residue masses for the Seed
	protected int[] SequencePRM;
	protected double[] SequencePRMScores;

	// Maintains the observed masses aligned to each insert state
	// This is volatile, and changes when the model changes
	protected ArrayList[] InsertMasses; // InsertMasses[i] =
	// {[SpectrumKey,Mass,Score],...}

	// HMM MODEL PARAMETERS

	// Transition probabilities
	// This is volatile and changes at each alignment
	protected double[][] Transitions;
	protected double[][] LogTransitions; // Same as above, but logs

	// Prior probabilities of starting at each state
	// This is equal for all match states and delete states, but is 0 for insert
	// states
	protected double[] Priors;
	protected double[] LogPriors; // Same as above, but logs

	// Emission parameters for each state
	// Match states are the only states with emission params
	protected DistributionType[] EmissionParams;

	// Model bookkeeping
	// These are volatile, and change with each model change
	protected int NumStates = 0;
	protected int NumMatchStates = 0; // Currently, there are an equal numbers
	// of match,insert, and delete states

	// The average of the peak scores aligned to this match state
	// Just bookkeepting
	// This is volatile, and changes at each alignment
	// protected double[] MatchStateScores;

	// Pretty descriptive name
	protected double MeanMatchStateScore;

	// The score of the current alignment
	// Currently this is a sum of the scores of the states
	protected double CurrAlignmentScore;
	
	protected boolean applyMatchPenalty = false;
	protected double matchPenalty = -1.0;

	/**
	 * Produces the pretty state string for a path in the current model.
	 * Depedendent on the layout of the match, insert, delete states in the
	 * state space. Assumes match states are numbered 0 to NumMatchStates,
	 * insert states are numbered NumMatchStates to 2*NumMatchStates, and delete
	 * states are numbered 2*NumMatchStates to 3*NumMatchStates.
	 * 
	 * @param Path
	 * @return
	 */
	public String GetStateString(int[] Path) {
		String StateString = "";
		for (int i = 0; i < Path.length; ++i) {
			if (IsMatchState(Path[i]))
				StateString += "s" + Path[i] + " ";
			else if (IsInsertState(Path[i]))
				StateString += "i" + (Path[i] - this.NumMatchStates) + " ";
			else if (IsDeleteState(Path[i]))
				StateString += "d" + (Path[i] - 2 * this.NumMatchStates) + " ";
			else {
				ErrorThrower
						.ThrowErrorCustum(
								ErrorThrower.CUSTOM_ERROR,
								"State "
										+ Path[i]
										+ " is not a recognized state in GetStateString");
			}
		}
		return StateString;
	}

	/**
	 * Returns the group with the highest 'score', with atleast the minimum
	 * number of peaks, doesn't require a central PRM
	 * 
	 * @param MatchState
	 *            The match state directly preceeding teh Insert state of
	 *            interest
	 * @return Returns the group of masses.
	 */
	protected ArrayList GetBestGroup(int MatchState, double Tolerance) {
		ArrayList CurrMasses = this.InsertMasses[MatchState];
		boolean LocalDebug = false;
		CurrMasses = SortInsertMasses(CurrMasses);
		ArrayList CurrBestGroup = null;
		double CurrBestGroupScore = 0.0;
		if (LocalDebug) {
			System.out.println("Considering insert after state " + MatchState
					+ ", " + CurrMasses.size() + " peaks");
			Utils.WaitForEnter();
		}
		for (int LeftMassIndex = 0; LeftMassIndex < CurrMasses.size(); ++LeftMassIndex) {
			Object[] LeftPeak = (Object[]) (CurrMasses.get(LeftMassIndex));
			int LeftMass = ((Integer) (LeftPeak[1])).intValue();
			ArrayList CurrGroup = new ArrayList();
			Hashtable CurrSpectraKeys = new Hashtable();
			CurrGroup.add(LeftPeak);
			CurrSpectraKeys.put((String) (LeftPeak[0]), LeftPeak);

			double CurrGroupScore = ((Double) (LeftPeak[2])).doubleValue();
			if (LocalDebug)
				System.out.println("Considering Left mass " + LeftMass);
			for (int RightMassIndex = LeftMassIndex + 1; RightMassIndex < CurrMasses
					.size(); ++RightMassIndex) {
				Object[] RightPeak = (Object[]) (CurrMasses.get(RightMassIndex));
				int RightMass = ((Integer) (RightPeak[1])).intValue();
				String RightSpectrumKey = (String) (RightPeak[0]);
				if (LocalDebug)
					System.out.println("Considering Right mass: " + RightMass);
				// See if we already have a mass from this spectrum!!
				if (CurrSpectraKeys.containsKey(RightSpectrumKey)) {
					if (LocalDebug)
						System.out.println("We already have a mass from "
								+ RightSpectrumKey);
					continue;
				}

				// Would this be a good addition?
				int PotentialMassSum = RightMass + GetInsertMassSum(CurrGroup);

				int PotentialMean = PotentialMassSum / (CurrGroup.size() + 1);
				double StdDev = Math.pow(PotentialMean - RightMass, 2)
						+ GetInsertMassVar(PotentialMean, CurrGroup);

				// StdDev = StdDev/(CurrGroup.size()+1);
				StdDev = Math.sqrt(StdDev);
				if (LocalDebug)
					System.out.println("PotentialMean: " + PotentialMean
							+ " PotentialStdDev: " + StdDev);
				if (2 * StdDev <= Tolerance * AntibodyUtils.MASS_SCALE) {
					if (LocalDebug)
						System.out
								.println("CLose enough to be in same grouP!!");
					CurrGroup.add(RightPeak);
					CurrSpectraKeys.put(RightSpectrumKey, RightPeak);
					CurrGroupScore += ((Double) (RightPeak[2])).doubleValue();
					if (LocalDebug) {
						System.out
								.println("CUrrGroupSize: " + CurrGroup.size());
						System.out.println("Curr Group Score: "
								+ CurrGroupScore);
					}
				} else
					break;
			}
			if (CurrGroupScore > CurrBestGroupScore
					&& CurrGroup.size() >= AntibodyUtils.MIN_PEAK_GROUP_SIZE) {
				CurrBestGroupScore = CurrGroupScore;
				CurrBestGroup = CurrGroup;
				if (LocalDebug)
					System.out.println("Current Best Group!!");
			}

		}
		if (LocalDebug) {
			Utils.WaitForEnter();
		}
		return CurrBestGroup;
	}

	public static ArrayList SortInsertMasses(ArrayList Peaks) {
		for (int i = 0; i < Peaks.size(); ++i) {
			Object[] CurrPeak = (Object[]) (Peaks.get(i));
			int MinMass = ((Integer) (CurrPeak[1])).intValue();
			int MinIndex = i;
			for (int j = i + 1; j < Peaks.size(); ++j) {
				Object[] OtherPeak = (Object[]) (Peaks.get(j));
				int CurrMass = ((Integer) (OtherPeak[1])).intValue();
				if (CurrMass < MinMass) {
					MinMass = CurrMass;
					MinIndex = j;
				}
			}
			if (MinIndex != i) {
				Object[] Temp = (Object[]) (Peaks.get(MinIndex));
				Peaks.set(MinIndex, (Peaks.get(i)));
				Peaks.set(i, Temp);
			}

		}
		return Peaks;
	}

	public static int GetInsertMassSum(ArrayList Peaks) {
		int Sum = 0;
		for (int i = 0; i < Peaks.size(); ++i) {
			Object[] CurrPeak = (Object[]) (Peaks.get(i));
			Sum += ((Integer) (CurrPeak[1])).intValue();

		}
		return Sum;

	}

	public static double GetInsertMassVar(int Mean, ArrayList Peaks) {
		double Var = 0.0;
		for (int k = 0; k < Peaks.size(); ++k) {
			Object[] CurrPeak = (Object[]) (Peaks.get(k));
			Var += Math.pow(Mean - ((Integer) (CurrPeak[1])).intValue(), 2);
		}
		return Var / (Peaks.size());
	}

	public int GetNumInserts(int State) {
		if (this.IsMatchState(State))
			return this.InsertMasses[State].size();
		else
			return -1;
	}

	public int GetNumMatchStates() {
		return this.NumMatchStates;

	}

	public double GetMatchMass(int State) {
		if (this.IsMatchState(State))
			return this.EmissionParams[State].GetMean();
		return Double.NEGATIVE_INFINITY;
	}

	public ArrayList GetAlignedSpectra() {
		return this.AlignedSpectra;
	}

	/**
	 * Create the data structure LogTransitions, identical to Transitions except
	 * with log probabilities instead of probabilities
	 * 
	 * @return Returns true if construction succeeds.
	 */

	protected boolean ConstructLogTransitions() {
		if (this.Transitions == null)
			return false;

		this.LogTransitions = new double[this.Transitions.length][this.Transitions[0].length];

		for (int startState = 0; startState < this.NumStates; startState++) {
			for (int endState = 0; endState < this.NumStates; endState++) {
				this.LogTransitions[startState][endState] = AntibodyUtils
						.Log10(this.Transitions[startState][endState]);

			}
		}
		return true;
	}

	/**
	 * Returns the log probability of emittion the observed peak from the given
	 * state. If the state is an insert state the peak must be between the means
	 * of the two adjacent match states, and the probability is exp(-PeakScore).
	 * This penalizes strong peaks that do not match a match state. For match
	 * states, the probability is sampled from the mass distribution for that
	 * state.
	 * 
	 * @param ObservedPeak
	 *            The observed peak mass (SCALED)
	 * @param MassScore
	 *            The PRM score of the observed peak
	 * @param State
	 *            The state the peak is observed at
	 * @return Returns the log probability of the emission of observed peak at
	 *         State.
	 */
	public abstract double GetEmissionLogProb(int ObservedPeak,
			double MassScore, int State);

	/**
	 * Basically the same as AlignSpectrum but updates the current HMMAlignment
	 * object, instead of creatining a new one.
	 * 
	 * @param Alignment
	 *            the HMMAlignment to update
	 * @return Returns true if update is successful
	 */
	public boolean ReAlignSpectrum(HMMAlignment Alignment, double Tolerance) {
		MSAlignmentType SpectrumAlignment = Alignment.GetMSAlignmentType();
		HMMAlignment NewAlignment = AlignSpectrum(SpectrumAlignment, Tolerance);
		Alignment.UpdateAlignment(NewAlignment.GetMLPath(),
				NewAlignment.GetPathScore(), this, Tolerance);

		if (Debug)
			System.out.println("Score: " + Alignment.GetPathScore() + " Path: "
					+ Alignment.GetPathString());
		return true;

	}

	/**
	 * Computes the mean score of the match states where score of the match
	 * state = Pi = pA^k(1-pA)^(n-k) where n is the number of overlapping
	 * spectra, k is the number supporting the match state i Qi =
	 * pNA^k(1-pNA)^(n-k) where n is the nubmer of overlapping spectra, k is the
	 * number supporting the match state i Score_i = log(Pi/Qi)
	 */
	public void ComputeMeanMatchScore() {
		double CumScore = 0.0; // Cumulative Score (sum of all match scores)
		for (int i = 0; i < this.NumMatchStates; ++i) {

			// CumScore += ComputeMatchStateLogOddsScore(i);
			if (!this.applyMatchPenalty)
				CumScore += ComputeMatchStateSumScore(i);
			else
				CumScore += ComputeMatchStateSumScoreWithPenalty(i);
			// CumScore += NumSupporting/NumOverlapping;
		}

		this.MeanMatchStateScore = CumScore / this.NumMatchStates;
	}

	/**
	 * Returns the mean score of the match states.
	 * 
	 * @see ComputeMeanMatchScore
	 * @return Returns the mean score of the match states
	 */
	public double GetMeanMatchScore() {
		this.ComputeMeanMatchScore();
		return this.MeanMatchStateScore;
	}

	public double ComputeMatchStateLogOddsScore(int MatchState) {
		int NumSupporting = this.OverlapCount[MatchState][0];
		int NumOverlapping = this.OverlapCount[MatchState][1];
		ArrayList AlignedPeaks = this.MatchedPeaks[MatchState];
		double P = 1;
		double Q = 1;
		for (int j = 0; j < AlignedPeaks.size(); ++j) {
			double[] CurrPeak = (double[]) (AlignedPeaks.get(j));
			int ScoreZone = AntibodyUtils.DeterminePRMScoreZone(CurrPeak[1]);
			P *= AntibodyUtils.pA[ScoreZone];
			Q *= AntibodyUtils.pNA[ScoreZone];
		}

		P *= Math.pow(AntibodyUtils.NotpA, NumOverlapping - NumSupporting);
		Q *= Math.pow(AntibodyUtils.NotpNA, NumOverlapping - NumSupporting);

		return AntibodyUtils.Log10(P / Q);
	}

	public double ComputeMatchStateSumScore(int MatchState) {

		ArrayList AlignedPeaks = this.MatchedPeaks[MatchState];

		double Total = 0.0;

		for (int j = 0; j < AlignedPeaks.size(); ++j) {
			double[] CurrPeak = (double[]) (AlignedPeaks.get(j));
			Total += CurrPeak[1];
		}
		return Total;

	}

	public double ComputeMatchStateSumScoreWithPenalty(int MatchState) {
		int NumSupporting = this.OverlapCount[MatchState][0];
		int NumOverlapping = this.OverlapCount[MatchState][1];
		ArrayList AlignedPeaks = this.MatchedPeaks[MatchState];

		if (AlignedPeaks == null) {
			System.out.println("WTF, Aligned peaks are null!!");
			System.out.println("MatchState: " + MatchState);
			System.out.println("NumMatchStates: " + NumMatchStates);
			System.out.println("MatchStateMean: "
					+ this.EmissionParams[MatchState].GetMean());
			DebugPrintSimple();
		}

		double Total = 0.0;

		for (int j = 0; j < AlignedPeaks.size(); ++j) {
			double[] CurrPeak = (double[]) (AlignedPeaks.get(j));
			Total += CurrPeak[1];
		}
		return Total -= (NumOverlapping - NumSupporting)
				* this.matchPenalty;

	}

	public abstract boolean IsMatchState(int stateNum);

	public abstract boolean IsInsertState(int stateNum);

	public abstract boolean IsDeleteState(int stateNum);

	public abstract int GetNextMatch(int stateNum);

	public abstract int GetNextInsert(int stateNum);

	public abstract int GetNextDelete(int stateNum);

	public abstract int GetPrevMatch(int stateNum);

	public abstract int GetPrevInsert(int stateNum);

	public abstract int GetPrevDelete(int stateNum);

	public abstract void DebugPrintSimple();

	public abstract void DebugPrint();

	public abstract void DebugPrintTransitions();

	public abstract int GetNumOverlapping(int State);

	public abstract int GetNumSupporting(int State);

	protected abstract boolean BuildBaseTransitions();

	public abstract boolean CreateBaseModel(double Tolerance);

	/**
	 * Performs Viterbi to align this spectrum to the model
	 * 
	 * @param NewAlignment
	 *            HMMAlignment that we are updating by aligning it's spectrum to
	 *            the model
	 * @return Returns true if the path finding succeeds, false if the path is
	 *         of poor quality or fails
	 */
	protected abstract boolean GetLikeliestPath(HMMAlignment NewAlignment,
			double Tolerance);

	/**
	 * This method should change the model topology by adding a new match state.
	 * It will change the model parameters, resetting everyting to initial
	 * values and will clear the HMMAlignments. ALL SPECTRA NEED TO BE
	 * REALIGNED. The realignment does not occur automatically in case the
	 * spectrum is not desirably matched after realignment.
	 * 
	 * @param MatchMass
	 *            The mean of the new match state emission distribution.
	 * @param MassDist
	 *            The list of masses that contribute to this matched mass.
	 * @return Returns true if the added state match succeeds, false otherwise.
	 */
	public abstract boolean AddMatchState(int[] MassDist, double Tolerance);

	/**
	 * This method should be called when a new alignment is added to the model.
	 * It updates the following parameters: Transitions
	 * 
	 * @return Returns true if there is success
	 */
	public abstract boolean UpdateModelTransitions();

	/**
	 * It shouldn't change the model topology, but it has been decided that it
	 * can be added.
	 * 
	 * @param A
	 *            an alignment of a specturm to the HMM
	 * @return Returns true if success, false otherwise
	 */
	public abstract boolean AddAlignmentToModel(HMMAlignment A);

	public abstract boolean AddAlignmentsToModel(HMMAlignment[] A);

	/**
	 * Computes the ML Path and returns the Consensus masses (SCALED)
	 * 
	 * @return Returns an array of the masses in the ML Consensus
	 */
	public abstract int[] ProduceConsensus();

	public abstract ConsensusPRMSpectrum ProduceConsensusWithScores();

	/**
	 * Populates the probability table P, such that P[i][j] is the probability
	 * of the residue represented by the PRM peaks i and j is correct. Populates
	 * the following intermediate tables alpha(i,t) = the probability of a path
	 * ending at peak i at time t (includes emission at state i) alpha(i,t) =
	 * Sum_k T_k,i * e_i(sym_t) * alpha(k,t-1) beta(j,t) = the probability of a
	 * path continuing from peak j at time t (does not include emission at state
	 * j) beta(j,t) = Sum_k T_j,k * e_k(sym_t+1) * beta(k,t+1) tau(i,j) = the
	 * probability of transitioning from state i to state j in one time step
	 * (does not include emission at state i, includes emission at state j)
	 * tau(i,j) = Sum_k T_i,k * e_k(sym_t1+1) * tau(t1+1,t2,k,j)
	 * 
	 * @param MatchedMasses
	 */
	public double[] ComputeProbabilityTable(int[] MatchedMasses) {
		// DetermineStates of the masses

		boolean LocalDebug = true;
		int[] StateNums = new int[MatchedMasses.length];

		for (int i = 0; i < MatchedMasses.length; ++i) {
			boolean Found = false;
			for (int j = 0; j < this.NumMatchStates; ++j) {
				if ((int) MatchedMasses[i] == (int) this.EmissionParams[j]
						.GetMean()) {
					StateNums[i] = j;
					Found = true;
					break;
				}
			}
			if (!Found) {
				System.err.println("ERROR: Could not find mass "
						+ MatchedMasses[i] + " in our model!!");
				for (int j = 0; j < this.NumMatchStates; ++j) {
					System.out.println((int) this.EmissionParams[j].GetMean());
				}
				Utils.WaitForEnter();
			}
		}

		if (LocalDebug) {
			this.DebugPrint();
			this.DebugPrintTransitions();
			Utils.WaitForEnter();
		}

		double[][] alpha = this.ComputeAlphaUsingInserts();
		double[][] beta = this.ComputeBetaUsingInserts();

		// Do some preprocessing
		// See how many times a single spectrum iterated on each loop
		int[] InsertMassCounts = this.GetInsertMassCounts();
		double[] InternalPathProbs = this.GetInsertStatePaths(InsertMassCounts);

		if (LocalDebug) {
			System.out.println("Insert MassCounts: "
					+ AntibodyUtils.IntegerListToString(InsertMassCounts));
			System.out.println("PathProbs: "
					+ Utils.DoubleArrayToString(InternalPathProbs));
		}
		int RunTime = 2 * this.NumMatchStates; // We could traverse every match
		// and every insert

		double ForwardProb = 0.0;
		for (int t = AntibodyUtils.MIN_PEAKS_ALIGNED - 1; t < RunTime; ++t) {
			for (int EndState = 0; EndState < this.NumMatchStates; ++EndState) {
				ForwardProb += alpha[EndState][t];
				ForwardProb += alpha[EndState + this.NumMatchStates][t];
				if (LocalDebug) {
					System.out.println(" += alpha[" + EndState + "][" + t
							+ "] = " + alpha[EndState][t]);
					System.out
							.println(" += alpha["
									+ (EndState + this.NumMatchStates) + "]["
									+ t + "] = "
									+ alpha[EndState + this.NumMatchStates][t]);
				}
			}
		}
		double ReverseProb = 0.0;
		for (int t = RunTime - 3 - 1; t >= 0; --t) {
			for (int StartState = 0; StartState < this.NumMatchStates; ++StartState) {
				ReverseProb += beta[StartState][t]
						* this.Priors[StartState]
						* this.EmissionParams[StartState]
								.ApproximatePDF(this.EmissionParams[StartState]
										.GetMean());
				ReverseProb += beta[StartState + this.NumMatchStates][t]
						* this.Priors[StartState + this.NumMatchStates]
						* InternalPathProbs[StartState];
				if (LocalDebug) {
					System.out
							.println(" += beta["
									+ StartState
									+ "]["
									+ t
									+ "]("
									+ beta[StartState][t]
									+ ") * "
									+ this.Priors[StartState]
									+ " * "
									+ this.EmissionParams[StartState]
											.ApproximatePDF(this.EmissionParams[StartState]
													.GetMean())
									+ " = "
									+ (beta[StartState][t]
											* this.Priors[StartState] * this.EmissionParams[StartState]
												.ApproximatePDF(this.EmissionParams[StartState]
														.GetMean())));
					System.out
							.println(" += beta["
									+ (StartState + this.NumMatchStates)
									+ "]["
									+ t
									+ "] * "
									+ this.Priors[StartState
											+ this.NumMatchStates]
									+ " * "
									+ InternalPathProbs[StartState]
									+ " = "
									+ (beta[StartState][t]
											* this.Priors[StartState] * InternalPathProbs[StartState]));

				}
			}
		}

		if (LocalDebug) {
			System.out.println("ForwardProb: " + ForwardProb);
			System.out.println("ReverseProb: " + ReverseProb);
			Utils.WaitForEnter();
		}

		double[][][] tau = this.ComputeTauUsingInserts();

		if (LocalDebug) {
			System.out.println("ALPHA: ");
			for (int i = 0; i < alpha.length; ++i) {
				String S = "";
				for (int j = 0; j < RunTime; ++j) {
					S += alpha[i][j] + " ";
				}
				System.out.println(S);
			}

			System.out.println("\nBETA: ");
			for (int i = 0; i < beta.length; ++i) {
				String S = "";
				for (int j = 0; j < RunTime; ++j) {
					S += beta[i][j] + " ";
				}
				System.out.println(S);
			}
			Utils.WaitForEnter();

			System.out.println("\nTau: ");
			for (int t = 0; t < RunTime; ++t) {
				System.out.println("Distance = " + t);
				for (int i = 0; i < tau.length; ++i) {
					String S = "";
					for (int j = 0; j < tau[0].length; ++j) {
						S += tau[i][j][t] + " ";
					}
					System.out.println(S);
				}
				Utils.WaitForEnter();
			}
		}

		double[] Scores = new double[MatchedMasses.length - 1];
		for (int i = 0; i < Scores.length; ++i) {
			Scores[i] = 0;
			for (int SpecLength = AntibodyUtils.MIN_PEAKS_ALIGNED; SpecLength < RunTime; SpecLength++) {
				int StartState = StateNums[i];
				int EndState = StateNums[i + 1];
				for (int t = 0; t < SpecLength; ++t) {
					for (int t2 = t + 1; t2 < SpecLength; ++t2) {

						int ActualTime = RunTime - (SpecLength - t2);
						if (LocalDebug) {
							System.out.println("t2: " + t2 + ", ActualTime: "
									+ ActualTime);
							System.out.println("alpha[" + StartState + "][" + t
									+ "] = " + alpha[StartState][t]);
							System.out.println("tau[" + StartState + "]["
									+ EndState + "][" + (t2 - t) + "] = "
									+ tau[StartState][EndState][t2 - t]);
							System.out.println("beta[" + EndState + "]["
									+ ActualTime + "] = "
									+ beta[EndState][ActualTime]);
						}
						Scores[i] += (alpha[StartState][t]
								* tau[StartState][EndState][t2 - t] * beta[EndState][ActualTime]);
						if (LocalDebug)
							System.out.println("Scores[" + i + "] = "
									+ Scores[i]);
					}
				}
				if (LocalDebug)
					Utils.WaitForEnter();
			}
			Scores[i] /= ForwardProb;
			System.out.println("Scores[" + i + "]:" + Scores[i]);
			Utils.WaitForEnter();
		}
		System.out.println("DONE");
		Utils.WaitForEnter();
		return Scores;
	}

	protected int[] GetInsertMassCounts() {
		int[] Ret = new int[this.NumMatchStates];
		for (int i = 0; i < this.NumMatchStates; ++i) {
			ArrayList Masses = this.InsertMasses[i];
			Hashtable Spectra = new Hashtable();
			int MaxCount = 0;
			for (int j = 0; j < Masses.size(); ++j) {
				Object[] Peak = (Object[]) (Masses.get(j));
				String Spectrum = (String) (Peak[0]);
				if (!Spectra.containsKey(Spectrum)) {
					Spectra.put(Spectrum, new Integer(1));
					if (1 > MaxCount)
						MaxCount = 1;
				} else {
					int Count = ((Integer) (Spectra.get(Spectrum))).intValue();
					Count += 1;
					Spectra.put(Spectrum, new Integer(Count));
					if (Count > MaxCount)
						MaxCount = Count;
				}
			}
			Ret[i] = MaxCount;
		}
		return Ret;
	}

	protected double[] GetInsertStatePaths(int[] InsertMassCounts) {
		double[] Ret = new double[InsertMassCounts.length];
		for (int i = 0; i < Ret.length; ++i) {
			for (int j = 0; j <= InsertMassCounts[i]; ++j) {
				Ret[i] += Math.pow(this.Transitions[i + this.NumMatchStates][i
						+ this.NumMatchStates], j);
			}
		}
		return Ret;
	}

	/**
	 * We count inserts as emitters, and take up one time step regardless of how
	 * many times we iterate on them
	 * 
	 * @return
	 */
	protected double[][] ComputeAlphaUsingInserts() {

		boolean LocalDebug = false;
		int RunTime = 2 * this.NumMatchStates; // We could traverse every match
		// and every insert
		double[][] Ret = new double[this.NumMatchStates * 3][RunTime];

		// Do some preprocessing
		// See how many times a single spectrum iterated on each loop
		int[] InsertMassCounts = this.GetInsertMassCounts();
		double[] InternalPathProbs = this.GetInsertStatePaths(InsertMassCounts);

		for (int i = 0; i < this.NumMatchStates; ++i) {
			// Match state = prior*emission
			Ret[i][0] = this.Priors[i]
					* this.EmissionParams[i]
							.ApproximatePDF(this.EmissionParams[i].GetMean());
			if (LocalDebug) {
				System.out.println("Ret["
						+ i
						+ "][0] = "
						+ this.Priors[i]
						+ "*"
						+ this.EmissionParams[i]
								.ApproximatePDF(this.EmissionParams[i]
										.GetMean()) + " = " + Ret[i][0]);
			}

			// Insert state = prior*zero or more transitions
			Ret[i + this.NumMatchStates][0] = this.Priors[i
					+ this.NumMatchStates]
					* InternalPathProbs[i];
			if (LocalDebug) {
				System.out.println("Ret[" + (i + this.NumMatchStates)
						+ "][0] = " + this.Priors[this.NumMatchStates] + "*"
						+ InternalPathProbs[i] + " = "
						+ Ret[i + this.NumMatchStates][0]);
			}
			// Delete state = sum of paths to previous states
			Ret[i + 2 * this.NumMatchStates][0] = 0;
			int PrevMatchState = this.GetPrevMatch(i + 2 * this.NumMatchStates);
			int PrevInsertState = this.GetPrevInsert(i + 2
					* this.NumMatchStates);
			int PrevDeleteState = this.GetPrevDelete(i + 2
					* this.NumMatchStates);
			if (PrevMatchState >= 0) {
				Ret[i + 2 * this.NumMatchStates][0] += Ret[PrevMatchState][0]
						* this.Transitions[PrevMatchState][i + 2
								* this.NumMatchStates];
				if (LocalDebug) {
					System.out.println("For delete state "
							+ (i + 2 * this.NumMatchStates)
							+ ", prev match is " + PrevMatchState);
					System.out.println("Ret["
							+ (i + 2 * this.NumMatchStates)
							+ "][0] += "
							+ Ret[PrevMatchState][0]
							+ "*"
							+ this.Transitions[PrevMatchState][i + 2
									* this.NumMatchStates] + " = "
							+ Ret[i + 2 * this.NumMatchStates][0]);
				}
			}
			if (PrevInsertState >= 0) {
				Ret[i + 2 * this.NumMatchStates][0] += Ret[PrevInsertState][0]
						* this.Transitions[PrevInsertState][i + 2
								* this.NumMatchStates];
				if (LocalDebug) {
					System.out.println("For delete state "
							+ (i + 2 * this.NumMatchStates)
							+ ", prev insert is " + PrevInsertState);
					System.out.println("Ret["
							+ (i + 2 * this.NumMatchStates)
							+ "][0] += "
							+ Ret[PrevInsertState][0]
							+ "*"
							+ this.Transitions[PrevInsertState][i + 2
									* this.NumMatchStates] + " = "
							+ Ret[i + 2 * this.NumMatchStates][0]);
				}
			}
			if (PrevDeleteState >= 0) {
				Ret[i + 2 * this.NumMatchStates][0] += Ret[PrevDeleteState][0]
						* this.Transitions[PrevDeleteState][i + 2
								* this.NumMatchStates];
				if (LocalDebug) {
					System.out.println("For delete state "
							+ (i + 2 * this.NumMatchStates)
							+ ", prev delete is " + PrevDeleteState);
					System.out.println("Ret["
							+ (i + 2 * this.NumMatchStates)
							+ "][0] += "
							+ Ret[PrevDeleteState][0]
							+ "*"
							+ this.Transitions[PrevDeleteState][i + 2
									* this.NumMatchStates] + " = "
							+ Ret[i + 2 * this.NumMatchStates][0]);
				}
			}
			if (LocalDebug)
				Utils.WaitForEnter();
		}

		for (int t = 1; t < RunTime; ++t) {
			for (int State = 0; State < this.NumMatchStates; ++State) {
				int CurrInsertState = State + this.NumMatchStates;
				int CurrDeleteState = State + 2 * this.NumMatchStates;

				/*
				 * Consider arrival at the match state State
				 */
				int PrevInsertState = this.GetPrevInsert(State);
				int PrevDeleteState = this.GetPrevDelete(State);
				int PrevMatchState = this.GetPrevMatch(State);

				if (PrevMatchState >= 0) {
					// Match State - Match State - traverse from previous time
					// point and emit
					Ret[State][t] += Ret[PrevMatchState][t - 1]
							* this.Transitions[PrevMatchState][State]
							* this.EmissionParams[State]
									.ApproximatePDF(this.EmissionParams[State]
											.GetMean());
					if (LocalDebug) {
						System.out.println("For match state " + (State)
								+ ", prev match is " + PrevMatchState);
						System.out
								.println("Ret["
										+ State
										+ "]["
										+ t
										+ "] += "
										+ Ret[PrevMatchState][t - 1]
										+ "*"
										+ this.Transitions[PrevMatchState][State]
										+ "*"
										+ this.EmissionParams[State]
												.ApproximatePDF(this.EmissionParams[State]
														.GetMean()) + " = "
										+ Ret[State][t]);
					}
				}

				if (PrevInsertState >= 0) {
					// Insert State - Match State - traverse from previous time
					// point and emit
					Ret[State][t] += Ret[PrevInsertState][t - 1]
							* this.Transitions[PrevInsertState][State]
							* this.EmissionParams[State]
									.ApproximatePDF(this.EmissionParams[State]
											.GetMean());
					if (LocalDebug) {
						System.out.println("For match state " + (State)
								+ ", prev insert is " + PrevInsertState);
						System.out
								.println("Ret["
										+ State
										+ "]["
										+ t
										+ "] += "
										+ Ret[PrevInsertState][t - 1]
										+ "*"
										+ this.Transitions[PrevInsertState][State]
										+ "*"
										+ this.EmissionParams[State]
												.ApproximatePDF(this.EmissionParams[State]
														.GetMean()) + " = "
										+ Ret[State][t]);
					}
				}
				if (PrevDeleteState >= 0) {
					// Delete State - Match State - traverse from previous time
					// point and emit
					Ret[State][t] += Ret[PrevDeleteState][t - 1]
							* this.Transitions[PrevDeleteState][State]
							* this.EmissionParams[State]
									.ApproximatePDF(this.EmissionParams[State]
											.GetMean());
					if (LocalDebug) {
						System.out.println("For match state " + (State)
								+ ", prev delete is " + PrevDeleteState);
						System.out
								.println("Ret["
										+ State
										+ "]["
										+ t
										+ "] += "
										+ Ret[PrevDeleteState][t - 1]
										+ "*"
										+ this.Transitions[PrevDeleteState][State]
										+ "*"
										+ this.EmissionParams[State]
												.ApproximatePDF(this.EmissionParams[State]
														.GetMean()) + " = "
										+ Ret[State][t]);
					}
				}

				/*
				 * Consider arrival at the Insert state CurrInsertState
				 */
				PrevDeleteState = this.GetPrevDelete(CurrInsertState);
				PrevMatchState = this.GetPrevMatch(CurrInsertState);

				if (PrevMatchState >= 0) {
					// Match State - Insert State - traverse from previous time
					// point and cycle up to k times
					Ret[CurrInsertState][t] += Ret[PrevMatchState][t - 1]
							* this.Transitions[PrevMatchState][CurrInsertState]
							* InternalPathProbs[State];
					if (LocalDebug) {
						System.out.println("For insert state "
								+ (CurrInsertState) + ", prev match is "
								+ PrevMatchState);
						System.out
								.println("Ret["
										+ CurrInsertState
										+ "]["
										+ t
										+ "] += "
										+ Ret[PrevMatchState][t - 1]
										+ "*"
										+ this.Transitions[PrevMatchState][CurrInsertState]
										+ "*" + InternalPathProbs[State]
										+ " = " + Ret[CurrInsertState][t]);
					}
				}
				if (PrevDeleteState >= 0) {
					// Delete State - Insert State = traverse from previous time
					// point and cycle up to k times
					Ret[CurrInsertState][t] += Ret[PrevDeleteState][t - 1]
							* this.Transitions[PrevDeleteState][CurrInsertState]
							* InternalPathProbs[State];
					if (LocalDebug) {
						System.out.println("For insert state "
								+ (CurrInsertState) + ", prev delete is "
								+ PrevDeleteState);
						System.out
								.println("Ret["
										+ CurrInsertState
										+ "]["
										+ t
										+ "] += "
										+ Ret[PrevDeleteState][t - 1]
										+ "*"
										+ this.Transitions[PrevDeleteState][CurrInsertState]
										+ "*" + InternalPathProbs[State]
										+ " = " + Ret[CurrInsertState][t]);
					}

				}

				/*
				 * Consider arrival at the delete state CurrDeleteState
				 */
				PrevDeleteState = this.GetPrevDelete(CurrDeleteState);
				PrevMatchState = this.GetPrevMatch(CurrDeleteState);
				PrevInsertState = this.GetPrevInsert(CurrDeleteState);

				if (PrevMatchState >= 0) {
					// Match State - Delete State - traverse from same time
					// point
					Ret[CurrDeleteState][t] += Ret[PrevMatchState][t]
							* this.Transitions[PrevMatchState][CurrDeleteState];
					if (LocalDebug) {
						System.out.println("For delete state "
								+ (CurrDeleteState) + ", prev match is "
								+ PrevMatchState);
						System.out
								.println("Ret["
										+ CurrDeleteState
										+ "]["
										+ t
										+ "] += "
										+ Ret[PrevMatchState][t]
										+ "*"
										+ this.Transitions[PrevMatchState][CurrDeleteState]
										+ " = " + Ret[CurrDeleteState][t]);
					}
				}
				if (PrevInsertState >= 0) {
					// Insert State - Delete State - traverse from same time
					// point
					Ret[CurrDeleteState][t] += Ret[PrevInsertState][t]
							* this.Transitions[PrevInsertState][CurrDeleteState];
					if (LocalDebug) {
						System.out.println("For delete state "
								+ (CurrDeleteState) + ", prev insert is "
								+ PrevInsertState);
						System.out
								.println("Ret["
										+ CurrDeleteState
										+ "]["
										+ t
										+ "] += "
										+ Ret[PrevInsertState][t]
										+ "*"
										+ this.Transitions[PrevInsertState][CurrDeleteState]
										+ " = " + Ret[CurrDeleteState][t]);
					}
				}
				if (PrevDeleteState >= 0) {
					// Delete State- Delete State - traverse from same time
					// point
					Ret[CurrDeleteState][t] += Ret[PrevDeleteState][t]
							* this.Transitions[PrevDeleteState][CurrDeleteState];
					if (LocalDebug) {
						System.out.println("For delete state "
								+ (CurrDeleteState) + ", prev delete is "
								+ PrevDeleteState);
						System.out
								.println("Ret["
										+ CurrDeleteState
										+ "]["
										+ t
										+ "] += "
										+ Ret[PrevDeleteState][t]
										+ "*"
										+ this.Transitions[PrevDeleteState][CurrDeleteState]
										+ " = " + Ret[CurrDeleteState][t]);
					}
				}
				if (LocalDebug)
					Utils.WaitForEnter();
			}
		}

		return Ret;
	}

	private double[][] ComputeBetaUsingInserts() {
		int RunTime = 2 * this.NumMatchStates; // We could traverse every match

		boolean LocalDebug = false;
		// and every insert
		double[][] Ret = new double[this.NumMatchStates * 3][RunTime];

		// Do some preprocessing
		// See how many times a single spectrum iterated on each loop
		int[] InsertMassCounts = this.GetInsertMassCounts();
		double[] InternalPathProbs = this.GetInsertStatePaths(InsertMassCounts);

		for (int i = this.NumMatchStates - 1; i >= 0; --i) {
			Ret[i][RunTime - 1] = 1.0;
			Ret[i + this.NumMatchStates][RunTime - 1] = 1.0;

			int CurrDeleteState = this.GetCorrespondingDelete(i);
			int NextInsertState = this.GetNextInsert(CurrDeleteState);
			int NextMatchState = this.GetNextMatch(CurrDeleteState);
			int NextDeleteState = this.GetNextDelete(CurrDeleteState);
			if (LocalDebug)

				// System.out.println("Delete: " + CurrDeleteState + " goes to "
				// + NextDeleteState + "," + NextInsertState + "," +
				// NextMatchState);
				if (NextMatchState >= 0) {
					Ret[CurrDeleteState][RunTime - 1] += Ret[NextMatchState][RunTime - 1]
							* this.Transitions[CurrDeleteState][NextMatchState]
							* this.EmissionParams[NextMatchState]
									.ApproximatePDF(this.EmissionParams[NextMatchState]
											.GetMean());
					if (LocalDebug) {
						System.out.println("For delete state "
								+ CurrDeleteState + ", next Match is "
								+ NextMatchState);
						System.out
								.println("Beta["
										+ CurrDeleteState
										+ "]["
										+ (RunTime - 1)
										+ "] += "
										+ Ret[NextMatchState][RunTime - 1]
										+ "*"
										+ this.Transitions[CurrDeleteState][NextMatchState]
										+ "*"
										+ this.EmissionParams[NextMatchState]
												.ApproximatePDF(this.EmissionParams[NextMatchState]
														.GetMean()) + " = "
										+ Ret[CurrDeleteState][RunTime - 1]);
					}
				}
			if (NextInsertState >= 0) {
				Ret[CurrDeleteState][RunTime - 1] += Ret[NextInsertState][RunTime - 1]
						* this.Transitions[CurrDeleteState][NextInsertState]
						* InternalPathProbs[NextInsertState
								- this.NumMatchStates];
				if (LocalDebug) {
					System.out.println("For delete state " + CurrDeleteState
							+ ", next Insert is " + NextInsertState);
					System.out
							.println("Beta["
									+ CurrDeleteState
									+ "]["
									+ (RunTime - 1)
									+ "] += "
									+ Ret[NextInsertState][RunTime - 1]
									+ "*"
									+ this.Transitions[CurrDeleteState][NextInsertState]
									+ "*"
									+ InternalPathProbs[NextInsertState
											- this.NumMatchStates] + " = "
									+ Ret[CurrDeleteState][RunTime - 1]);
				}
			}
			if (NextDeleteState >= 0) {
				Ret[CurrDeleteState][RunTime - 1] += Ret[NextDeleteState][RunTime - 1]
						* this.Transitions[CurrDeleteState][NextDeleteState];
				if (LocalDebug) {
					System.out.println("For delete state " + CurrDeleteState
							+ ", next delete is " + NextDeleteState);
					System.out
							.println("Beta["
									+ CurrDeleteState
									+ "]["
									+ (RunTime - 1)
									+ "] += "
									+ Ret[NextDeleteState][RunTime - 1]
									+ "*"
									+ this.Transitions[CurrDeleteState][NextDeleteState]
									+ " = " + Ret[CurrDeleteState][RunTime - 1]);
				}
			}
			if (LocalDebug)
				Utils.WaitForEnter();
		}

		for (int t = RunTime - 2; t >= 0; --t) {
			for (int State = this.NumMatchStates - 1; State >= 0; --State) {
				int CurrInsertState = State + this.NumMatchStates;
				int CurrDeleteState = State + 2 * this.NumMatchStates;

				int NextMatchState = this.GetNextMatch(State);
				int NextInsertState = this.GetNextInsert(State);
				int NextDeleteState = this.GetNextDelete(State);

				if (NextMatchState >= 0) {
					// Match State - Match State
					Ret[State][t] += Ret[NextMatchState][t + 1]
							* this.Transitions[State][NextMatchState]
							* this.EmissionParams[NextMatchState]
									.ApproximatePDF(this.EmissionParams[NextMatchState]
											.GetMean());

				}
				if (NextDeleteState >= 0) {
					// Match State - Delete State
					Ret[State][t] += Ret[NextDeleteState][t + 1]
							* this.Transitions[State][NextDeleteState];
				}
				if (NextInsertState >= 0) {
					// Match State - Insert State
					Ret[State][t] += Ret[NextInsertState][t + 1]
							* this.Transitions[State][NextInsertState]
							* InternalPathProbs[NextInsertState
									- this.NumMatchStates];
				}
				NextMatchState = this.GetNextMatch(CurrInsertState);
				NextDeleteState = this.GetNextDelete(CurrInsertState);
				if (NextMatchState >= 0) {
					// Insert State - Match State
					Ret[CurrInsertState][t] += Ret[NextMatchState][t + 1]
							* this.Transitions[CurrInsertState][NextMatchState]
							* this.EmissionParams[NextMatchState]
									.ApproximatePDF(this.EmissionParams[NextMatchState]
											.GetMean());
				}

				if (NextDeleteState >= 0) {
					// Insert State - Delete State
					Ret[CurrInsertState][t] += Ret[NextDeleteState][t + 1]
							* this.Transitions[CurrInsertState][NextDeleteState];
				}

				NextMatchState = this.GetNextMatch(CurrDeleteState);
				NextInsertState = this.GetNextInsert(CurrDeleteState);
				NextDeleteState = this.GetNextDelete(CurrDeleteState);

				if (NextMatchState >= 0) {
					// Delete State - Match State
					Ret[CurrDeleteState][t] += Ret[NextMatchState][t]
							* this.Transitions[CurrDeleteState][NextMatchState]
							* this.EmissionParams[NextMatchState]
									.ApproximatePDF(this.EmissionParams[NextMatchState]
											.GetMean());
				}
				if (NextInsertState >= 0) {
					// Delete State - Insert State
					Ret[CurrDeleteState][t] += Ret[NextInsertState][t]
							* this.Transitions[CurrDeleteState][NextInsertState]
							* InternalPathProbs[NextInsertState
									- this.NumMatchStates];
				}
				if (NextDeleteState >= 0) {
					// Delete State - Delete State
					Ret[CurrDeleteState][t] += Ret[NextDeleteState][t]
							* this.Transitions[CurrDeleteState][NextDeleteState];
				}

			}
		}

		return Ret;
	}

	/**
	 * Computed the probability of starting in a state i and ending in state j
	 * with no internal emissions from paths which traverse Delete and Match
	 * states only, and emit ideal masses.
	 * 
	 * @return
	 */
	private double[][][] ComputeTauUsingInserts() {
		// Try only considering match and delete states
		int RunTime = 2 * this.NumMatchStates; // We could traverse every match
		// and every insert
		double[][][] Ret = new double[this.NumMatchStates][3 * this.NumMatchStates][RunTime];

		int[] InsertMassCounts = this.GetInsertMassCounts();
		double[] InternalPathProbs = this.GetInsertStatePaths(InsertMassCounts);

		boolean LocalDebug = false;
		for (int State = 0; State < this.NumMatchStates; ++State) {

			Ret[State][State][0] = 1.0;
			int NextDeleteState = this.GetNextDelete(State);
			if (LocalDebug) {
				System.out.println("Ret[" + State + "][" + State + "][0] = "
						+ Ret[State][State][0]);
			}

			while (NextDeleteState >= 0) {
				int PrevDeleteState = this.GetPrevDelete(NextDeleteState);
				Ret[State][NextDeleteState][0] += Ret[State][State][0]
						* this.Transitions[State][NextDeleteState];
				if (LocalDebug) {
					System.out.println("After state " + State
							+ ", the next delete is " + NextDeleteState);
					System.out.println("Ret[" + State + "][" + NextDeleteState
							+ "][0] += " + Ret[State][State][0] + "*"
							+ this.Transitions[State][NextDeleteState] + " = "
							+ Ret[State][NextDeleteState][0]);
				}
				if (PrevDeleteState >= 0) {
					Ret[State][NextDeleteState][0] += Ret[State][PrevDeleteState][0]
							* this.Transitions[PrevDeleteState][NextDeleteState];
					if (LocalDebug)
						System.out
								.println("Ret["
										+ State
										+ "]["
										+ NextDeleteState
										+ "][0] += "
										+ Ret[State][PrevDeleteState][0]
										+ "*"
										+ this.Transitions[PrevDeleteState][NextDeleteState]
										+ " = "
										+ Ret[State][NextDeleteState][0]);
				}
				NextDeleteState = this.GetNextDelete(NextDeleteState);
			}
			if (LocalDebug)
				Utils.WaitForEnter();

		}
		// How many intervening emiting states
		for (int State = 0; State < this.NumMatchStates; ++State) // The
		// starting
		// match
		// state
		{
			for (int Interval = 1; Interval < RunTime; ++Interval) {
				for (int EndState = State; EndState < this.NumMatchStates; ++EndState) // The
																						// ending
																						// state
																						// (could
				// be match or insert or delete)
				{
					int MatchStateNum = EndState;
					int InsertStateNum = EndState + this.NumMatchStates;
					int DeleteStateNum = EndState + 2 * this.NumMatchStates;

					/*
					 * Handle transitions to match state, EndState first
					 */
					int PrevMatchState = this.GetPrevMatch(MatchStateNum);
					int PrevInsertState = this.GetPrevInsert(MatchStateNum);
					int PrevDeleteState = this.GetPrevDelete(MatchStateNum);

					if (PrevMatchState >= 0) {
						// Consider PrevState Match - CurrState Match, but only
						// if it's the starting match
						if (PrevMatchState == State)
							Ret[State][MatchStateNum][Interval] += Ret[State][PrevMatchState][Interval - 1]
									* this.Transitions[PrevMatchState][MatchStateNum]
									* this.EmissionParams[MatchStateNum]
											.ApproximatePDF(this.EmissionParams[MatchStateNum]
													.GetMean());
						if (LocalDebug) {
							System.out.println("For match " + MatchStateNum
									+ ", prev match is " + PrevMatchState);
							System.out
									.println("Ret["
											+ State
											+ "]["
											+ MatchStateNum
											+ "]["
											+ Interval
											+ "] += "
											+ Ret[State][PrevMatchState][Interval - 1]
											+ "*"
											+ this.Transitions[PrevMatchState][MatchStateNum]
											+ "*"
											+ this.EmissionParams[MatchStateNum]
													.ApproximatePDF(this.EmissionParams[MatchStateNum]
															.GetMean())
											+ " = "
											+ Ret[State][MatchStateNum][Interval]);
						}

					}
					if (PrevDeleteState >= 0) {
						// Consider PrevState Delete - CurrState Match
						Ret[State][MatchStateNum][Interval] += Ret[State][PrevDeleteState][Interval - 1]
								* this.Transitions[PrevDeleteState][MatchStateNum]
								* this.EmissionParams[MatchStateNum]
										.ApproximatePDF(this.EmissionParams[MatchStateNum]
												.GetMean());
						if (LocalDebug) {
							System.out.println("For match " + MatchStateNum
									+ ", prev delete is " + PrevDeleteState);
							System.out
									.println("Ret["
											+ State
											+ "]["
											+ MatchStateNum
											+ "]["
											+ Interval
											+ "] += "
											+ Ret[State][PrevDeleteState][Interval - 1]
											+ "*"
											+ this.Transitions[PrevDeleteState][MatchStateNum]
											+ "*"
											+ this.EmissionParams[MatchStateNum]
													.ApproximatePDF(this.EmissionParams[MatchStateNum]
															.GetMean())
											+ " = "
											+ Ret[State][MatchStateNum][Interval]);
						}
					}
					if (PrevInsertState >= 0) {
						// Consider PrevState Insert - CurrState Match
						Ret[State][MatchStateNum][Interval] += Ret[State][PrevInsertState][Interval - 1]
								* this.Transitions[PrevInsertState][MatchStateNum]
								* this.EmissionParams[MatchStateNum]
										.ApproximatePDF(this.EmissionParams[MatchStateNum]
												.GetMean());
						if (LocalDebug) {
							System.out.println("For match " + MatchStateNum
									+ ", prev insert is " + PrevInsertState);
							System.out
									.println("Ret["
											+ State
											+ "]["
											+ MatchStateNum
											+ "]["
											+ Interval
											+ "] += "
											+ Ret[State][PrevInsertState][Interval - 1]
											+ "*"
											+ this.Transitions[PrevInsertState][MatchStateNum]
											+ "*"
											+ this.EmissionParams[MatchStateNum]
													.ApproximatePDF(this.EmissionParams[MatchStateNum]
															.GetMean())
											+ " = "
											+ Ret[State][MatchStateNum][Interval]);
						}
					}

					/*
					 * Handle current insert state
					 */
					PrevMatchState = this.GetPrevMatch(InsertStateNum);
					PrevDeleteState = this.GetPrevDelete(InsertStateNum);
					if (PrevMatchState >= 0) {
						// Consider PrevState Match - CurrState Insert
						if (PrevMatchState == State)
							Ret[State][InsertStateNum][Interval] += Ret[State][PrevMatchState][Interval - 1]
									* this.Transitions[PrevMatchState][InsertStateNum]
									* InternalPathProbs[EndState];
						if (LocalDebug) {
							System.out.println("For insert " + InsertStateNum
									+ ", prev match is " + PrevMatchState);
							System.out
									.println("Ret["
											+ State
											+ "]["
											+ InsertStateNum
											+ "]["
											+ Interval
											+ "] += "
											+ Ret[State][PrevMatchState][Interval - 1]
											+ "*"
											+ this.Transitions[PrevMatchState][InsertStateNum]
											+ "*"
											+ InternalPathProbs[EndState]
											+ " = "
											+ Ret[State][InsertStateNum][Interval]);
						}
					}
					if (PrevDeleteState >= 0) {
						// Consider PrevState Delete - CurrState Insert
						Ret[State][InsertStateNum][Interval] += Ret[State][PrevDeleteState][Interval - 1]
								* this.Transitions[PrevDeleteState][InsertStateNum]
								* InternalPathProbs[EndState];
						if (LocalDebug) {
							System.out.println("For insert " + InsertStateNum
									+ ", prev delete is " + PrevDeleteState);
							System.out
									.println("Ret["
											+ State
											+ "]["
											+ InsertStateNum
											+ "]["
											+ Interval
											+ "] += "
											+ Ret[State][PrevDeleteState][Interval - 1]
											+ "*"
											+ this.Transitions[PrevDeleteState][InsertStateNum]
											+ "*"
											+ InternalPathProbs[EndState]
											+ " = "
											+ Ret[State][InsertStateNum][Interval]);
						}

					}
					// Since the delete probabilities at time t depend on the
					// match probabilities
					// at time t, we compute these last
					PrevInsertState = this.GetPrevInsert(DeleteStateNum);
					PrevDeleteState = this.GetPrevDelete(DeleteStateNum);

					if (PrevInsertState >= 0) {
						// Consider PrevState Insert - CurrState Delete
						Ret[State][DeleteStateNum][Interval] += Ret[State][PrevInsertState][Interval]
								* this.Transitions[PrevInsertState][DeleteStateNum];
						if (LocalDebug) {
							System.out.println("For delete " + DeleteStateNum
									+ ", prev insert is " + PrevInsertState);
							System.out
									.println("Ret["
											+ State
											+ "]["
											+ DeleteStateNum
											+ "]["
											+ Interval
											+ "] += "
											+ Ret[State][PrevInsertState][Interval]
											+ "*"
											+ this.Transitions[PrevInsertState][DeleteStateNum]
											+ " = "
											+ Ret[State][DeleteStateNum][Interval]);
						}
					}
					if (PrevDeleteState >= 0) {
						// Consider PrevState Delete - CurrState Delete
						Ret[State][DeleteStateNum][Interval] += Ret[State][PrevDeleteState][Interval]
								* this.Transitions[PrevDeleteState][DeleteStateNum];
						if (LocalDebug) {
							System.out.println("For delete " + DeleteStateNum
									+ ", prev delete is " + PrevDeleteState);
							System.out
									.println("Ret["
											+ State
											+ "]["
											+ DeleteStateNum
											+ "]["
											+ Interval
											+ "] += "
											+ Ret[State][PrevDeleteState][Interval]
											+ "*"
											+ this.Transitions[PrevDeleteState][DeleteStateNum]
											+ " = "
											+ Ret[State][DeleteStateNum][Interval]);
						}
					}
					if (LocalDebug)
						Utils.WaitForEnter();
				}
			}
		}

		return Ret;
	}

	/**
	 * Considers each insert state, and returns possible candidates for adding
	 * as a match state, only 1 per insert state. Currently the set of masses
	 * from the insert state that are no more than 2*IonTolerance apart with the
	 * highest ODDS score above ThresholdScore are considered for candidacy. A
	 * candidate (i.e. result[i]) is simply the mean of the masses within a
	 * 2*IonTolerance window that are the highest scoring (and above
	 * ThresholdScore) at Insert state i.
	 * 
	 * @param ThresholdScore
	 *            the cutoff for considering candidates
	 * @return Returns
	 */
	public abstract ArrayList[] GetMatchCandidatesOddsScore(
			double ThresholdScore, double Tolerance);

	public abstract ArrayList[] GetMatchCandidatesSumScore(
			double ThresholdScore, double Tolerance);

	/**
	 * Creates the HMMAlignment of this MSAlignment to the model. THe model is
	 * unchanged. MSAlignmentType may be unalignable if 1) the spectrum does not
	 * overlap the seed sequence, 2) the spectrum has a negative mass shift
	 * relative to the seed, 3) the seed sequence is not a part of the anchor
	 * sequence from SpecAlign for this MSAlignmentTYpe object.
	 * 
	 * @param SpectrumAlignment
	 *            an aligned spectrum result of the SpecAlign step
	 * @return Returns the HMMAlignment for this MSAlignmentType object, or null
	 *         if it is unalignable.
	 */
	public abstract HMMAlignment AlignSpectrum(
			MSAlignmentType SpectrumAlignment, double Tolerance);

	/**
	 * Considers each insert state, and returns possible candidates for adding
	 * as a match state, only 1 per insert state. Currently the set of masses
	 * from the insert state that are no more than 2*IonTolerance apart with the
	 * highest score above ThresholdScore are considered for candidacy. A
	 * candidate (i.e. result[i]) is simply the mean of the masses within a
	 * 2*IonTolerance window that are the highest scoring (and above
	 * ThresholdScore) at Insert state i.
	 * 
	 * @param ThresholdScore
	 *            the cutoff for considering candidates
	 * @return Returns
	 */
	public abstract ArrayList[] GetMatchCandidates(double ThresholdScore,
			double Tolerance);

	public int GetCorrespondingDelete(int MatchState) {
		if (!this.IsMatchState(MatchState)) {
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Cannot get corresponding delete for non-match state "
							+ MatchState);
		}

		return MatchState + 2 * this.NumMatchStates;
	}
}
