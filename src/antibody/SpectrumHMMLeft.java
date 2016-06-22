package antibody;

import java.util.ArrayList;
import java.util.Hashtable;

import basicUtils.Utils;
import errorUtils.ErrorThrower;

public class SpectrumHMMLeft extends SpectrumHMM {

	public SpectrumHMMLeft(String Seed, double Tolerance, boolean applyMatchPenalty, double matchPenalty) {
		this.Seed = Seed;
		this.applyMatchPenalty = applyMatchPenalty;
		this.matchPenalty = matchPenalty;
		String[] Els = AntibodyUtils.SplitStringtoAA(Seed);
		this.Seed = "";
		for (int i = 0; i < Math.min(10, Els.length); ++i)
			this.Seed += Els[i];

		if (Debug)
			System.out.println("Profile HMM Seq: " + this.Seed);

		// Create the SequencePRM for this seed (these are the initial
		// masses for the model)
		this.SequencePRM = AntibodyUtils.GetSeqPRMScaled(this.Seed);
		if (Debug) {
			String CurrString = "";
			for (int i = 0; i < this.SequencePRM.length; ++i)
				CurrString += this.SequencePRM[i] + " ";
			System.out.println("PRM: " + CurrString);
			Utils.WaitForEnter();
		}

		// Create the Base Model
		if (!this.CreateBaseModel(Tolerance)) {
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Unable to create the base model for the SpectrumHMM");
		}
	}

	public SpectrumHMMLeft(Seed SeqSeed, double Tolerance, boolean applyMatchPenalty, double matchPenalty) {
		this.Seed = SeqSeed.GetSeedSequence();
		this.applyMatchPenalty = applyMatchPenalty;
		this.matchPenalty = matchPenalty;
		String[] Els = AntibodyUtils.SplitStringtoAA(Seed);
		this.Seed = "";
		for (int i = 0; i < Math.min(10, Els.length); ++i)
			this.Seed += Els[i];

		Object[] StartSpectrum = AntibodyUtils.ReturnAllPreceeding(
				SeqSeed.GetSeedPRMs(),
				SeqSeed.UsedPRMs[Math.min(10, SeqSeed.UsedPRMs.length - 1)],
				SeqSeed.GetSeedPRMScores(), Tolerance);

		this.SequencePRM = (int[]) StartSpectrum[0];
		this.SequencePRMScores = (double[]) StartSpectrum[1];

		System.out
				.println("* Step 3: Construct a consensus spectrum from the recruited spectra");
		// System.out.println("**Creating a new SpectrumHMMLeft:");
		// System.out.println("Sequence: " + this.Seed);
		// System.out.println("PRMs: "
		// + AntibodyUtils.IntArrayToString(this.SequencePRM));
		// System.out.println("Scores: "
		// + AntibodyUtils.DoubleArrayToString(this.SequencePRMScores));
		// Create the Base Model
		if (!this.CreateBaseModelFromPRMs(Tolerance)) {
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Unable to create the base model for the SpectrumHMM");
		}
	}

	public boolean AddAlignmentToModel(HMMAlignment A) {
		// Add this spectrum to the list of AlignedSpectra
		this.AlignedSpectra.add(A);
		int[] CurrPath = A.GetMLPath();
		// int[] PRMPeaks = A.GetPRMSpectrum().PRMs;
		// double[] PRMScores = A.GetPRMSpectrum().PRMScores;

		int[] PRMPeaks = A.GetShiftedPRMs();
		double[] PRMScores = A.GetShiftedPRMScores();
		for (int i = 0; i < this.NumMatchStates; ++i)
			if (this.MatchedPeaks[i] == null)
				this.MatchedPeaks[i] = new ArrayList();

		// Update overlap count, emiissionparams and insertmasses
		int PRMPeakIndex = 0;
		int PrevState = -1;
		for (int i = 0; i < CurrPath.length; ++i) {
			if (IsMatchState(CurrPath[i])) {

				double[] CurrPeak = new double[2];
				CurrPeak[0] = PRMPeaks[PRMPeakIndex];
				CurrPeak[1] = PRMScores[PRMPeakIndex];

				OverlapCount[CurrPath[i]][0] += 1;// Update the supporting
				// number
				OverlapCount[CurrPath[i]][1] += 1;// Update the total number
				if (this.MatchedPeaks[CurrPath[i]] == null)
					this.MatchedPeaks[CurrPath[i]] = new ArrayList();
				this.MatchedPeaks[CurrPath[i]].add(CurrPeak);

				this.EmissionParams[CurrPath[i]]
						.AddNewSample(PRMPeaks[PRMPeakIndex]);
				PRMPeakIndex += 1;

				if (i > 0)
					this.InsertOverlapCount[CurrPath[i]] += 1;
			} else if (IsInsertState(CurrPath[i])) {
				// [] NewPeak =
				// {PRMPeaks[PRMPeakIndex],PRMScores[PRMPeakIndex]};
				Object[] NewPeak = {
						A.GetPRMSpectrum().FileName + ":"
								+ A.GetPRMSpectrum().SpecIndex,
						new Integer(PRMPeaks[PRMPeakIndex]),
						new Double(PRMScores[PRMPeakIndex]) };

				// this.InsertMasses[this.GetNextMatch(CurrPath[i]) -
				// 1].add(NewPeak);
				this.InsertMasses[CurrPath[i] - this.NumMatchStates]
						.add(NewPeak);
				if (i == CurrPath.length - 1)
					this.InsertOverlapCount[CurrPath[i] - this.NumMatchStates] += 1;
				PRMPeakIndex += 1;
			} else if (IsDeleteState(CurrPath[i])) {
				if (i > 0)
					this.InsertOverlapCount[CurrPath[i] - 2
							* this.NumMatchStates] += 1;
				OverlapCount[CurrPath[i] - 2 * this.NumMatchStates][1] += 1; // Update
				// the
				// total
				// number

			}
			PrevState = CurrPath[i];
		}
		this.ComputeMeanMatchScore();
		this.UpdateModelTransitions();

		return true;
	}

	public boolean AddAlignmentsToModel(HMMAlignment[] A) {
		// Add this spectrum to the list of AlignedSpectra

		for (int i = 0; i < this.NumMatchStates; ++i)
			if (this.MatchedPeaks[i] == null)
				this.MatchedPeaks[i] = new ArrayList();
		for (int k = 0; k < A.length; ++k) {
			if (A[k] == null)
				continue;
			this.AlignedSpectra.add(A[k]);
			int[] CurrPath = A[k].GetMLPath();
			// int[] PRMPeaks = A[k].GetPRMSpectrum().PRMs;
			// double[] PRMScores = A[k].GetPRMSpectrum().PRMScores;

			int[] PRMPeaks = A[k].GetShiftedPRMs();
			double[] PRMScores = A[k].GetShiftedPRMScores();

			for (int i = 0; i < this.NumMatchStates; ++i)
				if (this.MatchedPeaks[i] == null) {
					System.out.println("** A:MatchedPeaks[" + i
							+ "] is inexplicably null!!");
					Utils.WaitForEnter();
				}
			// Update overlap count, emissionparams and insertmasses
			int PRMPeakIndex = 0;
			for (int i = 0; i < CurrPath.length; ++i) {
				for (int r = 0; r < this.NumMatchStates; ++r)
					if (this.MatchedPeaks[r] == null) {
						System.out.println("** B" + i + ":MatchedPeaks[" + r
								+ "] is inexplicably null!!");
						Utils.WaitForEnter();
					}
				if (IsMatchState(CurrPath[i])) {
					double[] CurrPeak = new double[2];
					CurrPeak[0] = PRMPeaks[PRMPeakIndex];
					CurrPeak[1] = PRMScores[PRMPeakIndex];

					OverlapCount[CurrPath[i]][0] += 1;// Update the supporting
					// number
					OverlapCount[CurrPath[i]][1] += 1;// Update the total number
					if (this.MatchedPeaks[CurrPath[i]] == null)
						this.MatchedPeaks[CurrPath[i]] = new ArrayList();
					this.MatchedPeaks[CurrPath[i]].add(CurrPeak);

					this.EmissionParams[CurrPath[i]]
							.AddNewSample(PRMPeaks[PRMPeakIndex]);

					PRMPeakIndex += 1;
					if (i > 0)
						this.InsertOverlapCount[CurrPath[i]] += 1;
				} else if (IsInsertState(CurrPath[i])) {
					// double[] NewPeak =
					// {PRMPeaks[PRMPeakIndex],PRMScores[PRMPeakIndex]};
					Object[] NewPeak = {
							A[k].GetPRMSpectrum().FileName + ":"
									+ A[k].GetPRMSpectrum().SpecIndex,
							new Integer(PRMPeaks[PRMPeakIndex]),
							new Double(PRMScores[PRMPeakIndex]) };

					this.InsertMasses[CurrPath[i] - this.NumMatchStates]
							.add(NewPeak);
					PRMPeakIndex += 1;
					if (i == CurrPath.length - 1)
						this.InsertOverlapCount[CurrPath[i]
								- this.NumMatchStates] += 1;
				} else if (IsDeleteState(CurrPath[i])) {
					if (i > 0)
						this.InsertOverlapCount[CurrPath[i] - 2
								* this.NumMatchStates] += 1;
					OverlapCount[CurrPath[i] - 2 * this.NumMatchStates][1] += 1; // Update
					// the
					// total
					// number

				}

				for (int r = 0; r < this.NumMatchStates; ++r)
					if (this.MatchedPeaks[r] == null) {
						System.out.println("** B" + i + "End :MatchedPeaks["
								+ r + "] is inexplicably null!!");
						Utils.WaitForEnter();
					}

			}
			for (int r = 0; r < this.NumMatchStates; ++r)
				if (this.MatchedPeaks[r] == null) {
					System.out.println("** C" + k + "End :MatchedPeaks[" + r
							+ "] is inexplicably null!!");
					Utils.WaitForEnter();
				}
		}
		for (int i = 0; i < this.NumMatchStates; ++i)
			if (this.MatchedPeaks[i] == null) {
				System.out.println("** C:MatchedPeaks[" + i
						+ "] is inexplicably null!!");
				Utils.WaitForEnter();
			}

		// this.DebugPrintSimple();
		this.ComputeMeanMatchScore();
		this.UpdateModelTransitions();
		return true;
	}

	public boolean AddMatchState(int[] MassDist, double Tolerance) {
		boolean LocalDebug = false;
		int NewMatchMean = 0; // The mean of the new match state
		double[] Samples = new double[MassDist.length];
		for (int i = 0; i < MassDist.length; ++i) {
			if (LocalDebug)
				System.out
						.println("Contributor to match state: " + MassDist[i]);
			NewMatchMean += MassDist[i];
			Samples[i] = (double) (MassDist[i]);
		}
		NewMatchMean = NewMatchMean / MassDist.length;

		if (Debug || LocalDebug)
			System.out.println("New MatchMean: " + NewMatchMean);
		int NewMatchIndex = 0; // The index of the new match state
		while (NewMatchIndex < this.NumMatchStates
				&& NewMatchMean > this.EmissionParams[NewMatchIndex].GetMean())
			NewMatchIndex++;

		// The mapping of new states numbers to old state numbers
		int[] New2OldStates = new int[this.NumStates + 3];
		for (int i = 0; i < New2OldStates.length; ++i)
			New2OldStates[i] = -1;

		if (Debug || LocalDebug)
			System.out.println("Index for the next match: " + NewMatchIndex);
		// Compute the mapping
		for (int i = 0; i < this.NumMatchStates + 1; ++i) {
			if (i < NewMatchIndex) {
				New2OldStates[i] = i;
				New2OldStates[i + this.NumMatchStates + 1] = i
						+ this.NumMatchStates;
				New2OldStates[i + 2 * this.NumMatchStates + 2] = i + 2
						* this.NumMatchStates;
			}
			if (i > NewMatchIndex) {
				New2OldStates[i] = i - 1;
				New2OldStates[i + this.NumMatchStates + 1] = i
						+ this.NumMatchStates - 1;
				New2OldStates[i + 2 * this.NumMatchStates + 2] = i + 2
						* this.NumMatchStates - 1;
			}
			if (Debug && !Debug) {
				System.out.println("New: " + i + " = " + New2OldStates[i]);
				System.out.println("New: " + (i + this.NumMatchStates + 1)
						+ " = " + New2OldStates[i + this.NumMatchStates + 1]);
				System.out.println("New: " + (i + 2 * this.NumMatchStates + 2)
						+ " = "
						+ New2OldStates[i + 2 * this.NumMatchStates + 2]);
				Utils.WaitForEnter();
			}
		}

		if (Debug || LocalDebug) {
			System.out.println("New to old state mappings:");
			for (int k = 0; k < New2OldStates.length; ++k) {
				System.out.println(k + " - " + New2OldStates[k]);
			}
			Utils.WaitForEnter();
		}

		// Recreate model parameters (need to add 3 more states)
		this.NumMatchStates += 1;
		this.NumStates += 3;
		this.Priors = new double[3 * this.NumMatchStates];
		this.InsertMasses = new ArrayList[this.NumMatchStates];
		this.OverlapCount = new int[this.NumMatchStates][2];
		this.MatchedPeaks = new ArrayList[this.NumMatchStates];
		this.InsertOverlapCount = new int[this.NumMatchStates];
		this.AlignedSpectra = new ArrayList();

		// Add the new emission param element
		DistributionType[] NewEmissionParams = new DistributionType[this.NumMatchStates];
		for (int i = 0; i < this.NumMatchStates; ++i) {
			if (New2OldStates[i] == -1)
				NewEmissionParams[i] = new DistributionType(Samples,
						AntibodyUtils.GetVariance(NewMatchMean, Tolerance));
			else
				NewEmissionParams[i] = this.EmissionParams[New2OldStates[i]];
			this.Priors[i] = (1.0 / (2 * this.NumMatchStates));
			this.Priors[i + this.NumMatchStates] = (1.0 / (2 * this.NumMatchStates));
			this.Priors[2 * this.NumMatchStates + i] = 0.0;
			this.InsertMasses[i] = new ArrayList();
			this.MatchedPeaks[i] = null;

			this.InsertOverlapCount[i] = 0;

			// If this state is from the original sequence, set the initial
			// pseudocounts
			int OrigIndex = AntibodyUtils.MassIndex(this.SequencePRM,
					(int) (NewEmissionParams[i].GetMean()), Tolerance);
			if (OrigIndex >= 0) {
				this.OverlapCount[i][0] = AntibodyUtils.MIN_ALIGNED_SPECTRA;
				this.OverlapCount[i][1] = AntibodyUtils.MIN_ALIGNED_SPECTRA;
				double[] CurrPeak = new double[2];
				CurrPeak[0] = NewEmissionParams[i].GetMean();
				CurrPeak[1] = this.SequencePRMScores[OrigIndex]
						/ AntibodyUtils.MIN_ALIGNED_SPECTRA;
				if (this.MatchedPeaks[i] == null)
					this.MatchedPeaks[i] = new ArrayList();
				for (int k = 0; k < AntibodyUtils.MIN_ALIGNED_SPECTRA; ++k)
					this.MatchedPeaks[i].add(CurrPeak);

				// Assumes no internal insertions
				if (i > this.NumMatchStates - this.SequencePRM.length)
					this.InsertOverlapCount[i] = AntibodyUtils.MIN_ALIGNED_SPECTRA;

			} else if (AntibodyUtils.MassWithinPRMRange(this.SequencePRM,
					(int) (NewEmissionParams[i].GetMean()), Tolerance)) {
				this.OverlapCount[i][1] = AntibodyUtils.MIN_ALIGNED_SPECTRA;
				// this.InsertOverlapCount[i] = 0;
				// Assumes no internal insertions
				if (i > this.NumMatchStates - this.SequencePRM.length)
					this.InsertOverlapCount[i] = AntibodyUtils.MIN_ALIGNED_SPECTRA;
			}
		}
		this.EmissionParams = NewEmissionParams;

		double[][] NewTransitions = new double[this.NumStates][this.NumStates];
		for (int i = 0; i < this.NumStates; ++i) {
			int OldStateStart = New2OldStates[i];
			for (int j = 0; j < this.NumStates; ++j) {
				int OldStateEnd = New2OldStates[j];

				if (Debug && !Debug) {
					System.out.println("Updating " + i + " to " + j
							+ " which is old " + OldStateStart + " to "
							+ OldStateEnd);
				}
				NewTransitions[i][j] = 0.0;

				boolean DirectlyPreceeds = (NewMatchIndex != 0)
						&& ((i == NewMatchIndex - 1)
								|| (i - this.NumMatchStates == NewMatchIndex) || i
								- 2 * this.NumMatchStates == NewMatchIndex - 1);

				// This transition involves a new state or the old source
				// insert, so take the base transition prob
				if (DirectlyPreceeds || OldStateStart == -1
						|| OldStateEnd == -1) {
					if (this.IsMatchState(i)) {
						if (this.IsMatchState(j) && j == this.GetNextMatch(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_MATCH_MATCH;
						else if (this.IsInsertState(j)
								&& j == this.GetNextInsert(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_MATCH_INSERT;
						else if (this.IsDeleteState(j)
								&& j == this.GetNextDelete(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_MATCH_DELETE;
					} else if (this.IsInsertState(i)) {
						if (this.IsMatchState(j) && j == this.GetNextMatch(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_INSERT_MATCH;
						else if (this.IsInsertState(j)
								&& j == this.GetNextInsert(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_INSERT_INSERT;
						else if (this.IsDeleteState(j)
								&& j == this.GetNextDelete(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_INSERT_DELETE;
					} else if (this.IsDeleteState(i)) {
						if (this.IsMatchState(j) && j == this.GetNextMatch(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_DELETE_MATCH;
						else if (this.IsInsertState(j)
								&& j == this.GetNextInsert(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_DELETE_INSERT;
						else if (this.IsDeleteState(j)
								&& j == this.GetNextDelete(i))
							NewTransitions[i][j] = AntibodyUtils.INIT_DELETE_DELETE;
					}

				} else {
					if (this.IsMatchState(i)) {
						if (this.IsMatchState(j) && j == this.GetNextMatch(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
						else if (this.IsInsertState(j)
								&& j == this.GetNextInsert(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
						else if (this.IsDeleteState(j)
								&& j == this.GetNextDelete(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
					} else if (this.IsInsertState(i)) {
						if (this.IsMatchState(j) && j == this.GetNextMatch(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
						else if (this.IsInsertState(j)
								&& j == this.GetNextInsert(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
						else if (this.IsDeleteState(j)
								&& j == this.GetNextDelete(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
					} else if (this.IsDeleteState(i)) {
						if (this.IsMatchState(j) && j == this.GetNextMatch(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
						else if (this.IsInsertState(j)
								&& j == this.GetNextInsert(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
						else if (this.IsDeleteState(j)
								&& j == this.GetNextDelete(i))
							NewTransitions[i][j] = this.Transitions[OldStateStart][OldStateEnd];
					}
				}

			}
		}
		this.Transitions = NewTransitions;
		// this.Transitions[this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;
		// this.Transitions[2*this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;
		// this.Transitions[3*this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;

		if (Debug || LocalDebug) {
			this.DebugPrintTransitions();
			Utils.WaitForEnter();
		}
		if (!ConstructLogTransitions()) {
			System.err.println("ERROR: Unable to create base transition table");
			return false;
		}

		this.LogPriors = AntibodyUtils.Log10Array(this.Priors);
		return true;
	}

	protected boolean BuildBaseTransitions() {
		if (this.NumStates == 0)
			return false;

		this.Transitions = new double[this.NumStates][this.NumStates];

		for (int startState = 0; startState < this.NumStates; startState++) {
			for (int endState = 0; endState < this.NumStates; endState++) {
				this.Transitions[startState][endState] = 0.0;
				if (IsMatchState(startState)) {
					if (IsMatchState(endState)
							&& endState == GetNextMatch(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_MATCH_MATCH;
					else if (IsInsertState(endState)
							&& endState == GetNextInsert(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_MATCH_INSERT;
					else if (IsDeleteState(endState)
							&& endState == GetNextDelete(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_MATCH_DELETE;

				} else if (IsInsertState(startState)) {
					if (IsMatchState(endState)
							&& endState == GetNextMatch(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_INSERT_MATCH;
					else if (IsInsertState(endState)
							&& endState == GetNextInsert(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_INSERT_INSERT;
					else if (IsDeleteState(endState)
							&& endState == GetNextDelete(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_INSERT_DELETE;

				} else if (IsDeleteState(startState)) {
					if (IsMatchState(endState)
							&& endState == GetNextMatch(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_DELETE_MATCH;
					else if (IsInsertState(endState)
							&& endState == GetNextInsert(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_DELETE_INSERT;
					else if (IsDeleteState(endState)
							&& endState == GetNextDelete(startState))
						this.Transitions[startState][endState] = AntibodyUtils.INIT_DELETE_DELETE;

				}

			}
		}

		// this.Transitions[this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;
		// this.Transitions[2*this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;
		// this.Transitions[3*this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;

		return true;
	}

	public boolean CreateBaseModel(double Tolerance) {

		this.NumStates = 0;
		this.NumMatchStates = 0;

		this.EmissionParams = new DistributionType[this.SequencePRM.length];
		// this.MatchStateScores = new double[this.SequencePRM.length];
		this.Priors = new double[3 * this.SequencePRM.length];
		this.InsertMasses = new ArrayList[this.SequencePRM.length];

		this.NumStates = this.SequencePRM.length;
		this.NumMatchStates = this.SequencePRM.length;

		this.AlignedSpectra = new ArrayList();

		this.OverlapCount = new int[this.NumMatchStates][2];
		this.MatchedPeaks = new ArrayList[this.NumMatchStates];
		this.InsertOverlapCount = new int[this.NumMatchStates];
		// Initialize start states (states 0 - SequencePRM.length)
		for (int i = 0; i < this.NumMatchStates; ++i) {
			this.OverlapCount[i][0] = AntibodyUtils.MIN_ALIGNED_SPECTRA;
			this.OverlapCount[i][1] = AntibodyUtils.MIN_ALIGNED_SPECTRA;

			if (i > 0)
				this.InsertOverlapCount[i] = AntibodyUtils.MIN_ALIGNED_SPECTRA;
			this.EmissionParams[i] = new DistributionType(this.SequencePRM[i],
					AntibodyUtils.GetVariance(this.SequencePRM[i], Tolerance));

			this.MatchedPeaks[i] = new ArrayList();
			double[] CurrPeak = new double[2];
			CurrPeak[0] = this.EmissionParams[i].GetMean();
			CurrPeak[1] = AntibodyUtils.MAX_PRM_SCORE;
			for (int k = 0; k < AntibodyUtils.MIN_ALIGNED_SPECTRA; ++k)
				this.MatchedPeaks[i].add(CurrPeak);

			// this.MatchStateScores[i] = Utils.INIT_MATCH_STATE_SCORE;
			this.Priors[i] = (1.0 / (2 * this.SequencePRM.length));
		}

		// this.MeanMatchStateScore =
		// GetLogOdds(Utils.MIN_ALIGNED_SPECTRA,Utils.MIN_ALIGNED_SPECTRA);
		// this.MeanMatchStateScore = Utils.INIT_MATCH_STATE_SCORE;
		// this.MeanMatchStateScore = 1.0;
		this.ComputeMeanMatchScore();

		// Initialize insert states (states SequencePRM.length -
		// 2*SequencePRM.length)
		for (int i = 0; i < this.NumMatchStates; ++i) {
			this.NumStates += 1;
			this.Priors[this.NumMatchStates + i] = (1.0 / (2 * this.SequencePRM.length));
			this.InsertMasses[i] = new ArrayList();
		}

		// Initialize Delete States
		for (int i = 0; i < this.NumMatchStates; ++i) {
			this.NumStates += 1;
			this.Priors[2 * this.NumMatchStates + i] = 0.0;
		}

		if (!BuildBaseTransitions() || !ConstructLogTransitions()) {
			System.err.println("ERROR: Unable to create base transition table");
			return false;
		}

		this.LogPriors = AntibodyUtils.Log10Array(this.Priors);

		if (Debug && !Debug) {
			this.DebugPrintSimple();
			System.out.println("Sequence: " + this.Seed);
			this.DebugPrintTransitions();
		}
		return true;
	}

	public boolean CreateBaseModelFromPRMs(double Tolerance) {

		this.NumStates = 0;
		this.NumMatchStates = 0;

		this.EmissionParams = new DistributionType[this.SequencePRM.length];
		// this.MatchStateScores = new double[this.SequencePRM.length];
		this.Priors = new double[3 * this.SequencePRM.length];
		this.InsertMasses = new ArrayList[this.SequencePRM.length];

		this.NumStates = this.SequencePRM.length;
		this.NumMatchStates = this.SequencePRM.length;

		this.AlignedSpectra = new ArrayList();

		this.OverlapCount = new int[this.NumMatchStates][2];
		this.MatchedPeaks = new ArrayList[this.NumMatchStates];
		this.InsertOverlapCount = new int[this.NumMatchStates];
		// Initialize start states (states 0 - SequencePRM.length)
		for (int i = 0; i < this.NumMatchStates; ++i) {
			this.OverlapCount[i][0] = AntibodyUtils.MIN_ALIGNED_SPECTRA;
			this.OverlapCount[i][1] = AntibodyUtils.MIN_ALIGNED_SPECTRA;

			if (i > 0)
				this.InsertOverlapCount[i] = AntibodyUtils.MIN_ALIGNED_SPECTRA;
			this.EmissionParams[i] = new DistributionType(this.SequencePRM[i],
					AntibodyUtils.GetVariance(this.SequencePRM[i], Tolerance));

			this.MatchedPeaks[i] = new ArrayList();
			double[] CurrPeak = new double[2];
			CurrPeak[0] = this.EmissionParams[i].GetMean();
			CurrPeak[1] = this.SequencePRMScores[i]
					/ AntibodyUtils.MIN_ALIGNED_SPECTRA;
			for (int k = 0; k < AntibodyUtils.MIN_ALIGNED_SPECTRA; ++k)
				this.MatchedPeaks[i].add(CurrPeak);

			// this.MatchStateScores[i] = Utils.INIT_MATCH_STATE_SCORE;
			this.Priors[i] = (1.0 / (2 * this.SequencePRM.length));
		}

		// this.MeanMatchStateScore =
		// GetLogOdds(Utils.MIN_ALIGNED_SPECTRA,Utils.MIN_ALIGNED_SPECTRA);
		// this.MeanMatchStateScore = Utils.INIT_MATCH_STATE_SCORE;
		// this.MeanMatchStateScore = 1.0;
		this.ComputeMeanMatchScore();

		// Initialize insert states (states SequencePRM.length -
		// 2*SequencePRM.length)
		for (int i = 0; i < this.NumMatchStates; ++i) {
			this.NumStates += 1;
			this.Priors[this.NumMatchStates + i] = (1.0 / (2 * this.SequencePRM.length));
			this.InsertMasses[i] = new ArrayList();
		}

		// Initialize Delete States
		for (int i = 0; i < this.NumMatchStates; ++i) {
			this.NumStates += 1;
			this.Priors[2 * this.NumMatchStates + i] = 0.0;
		}

		if (!BuildBaseTransitions() || !ConstructLogTransitions()) {
			System.err.println("ERROR: Unable to create base transition table");
			return false;
		}

		this.LogPriors = AntibodyUtils.Log10Array(this.Priors);

		if (Debug && !Debug) {
			this.DebugPrintSimple();
			System.out.println("Sequence: " + this.Seed);
			this.DebugPrintTransitions();
		}
		return true;
	}

	public void DebugPrint() {
		System.out.println("----------------SpectrumHMM---------------");
		System.out.println("Num States: " + this.NumStates
				+ " Num Match States: " + this.NumMatchStates);
		for (int i = 0; i < this.NumMatchStates; ++i) {
			System.out.println("Match State[" + i + "] : "
					+ this.EmissionParams[i].GetMean() + " "
					+ this.OverlapCount[i][0] + "/" + this.OverlapCount[i][1]
					+ " prior: " + this.Priors[i]);
			System.out.println("Insert State[" + i + "] prior: "
					+ this.Priors[i + this.NumMatchStates]);
			for (int j = 0; j < this.InsertMasses[i].size(); ++j) {
				Object[] CurrMass = (Object[]) (this.InsertMasses[i].get(j));
				System.out.println("  " + ((Integer) (CurrMass[1])).intValue()
						+ " - " + ((Double) (CurrMass[2])).doubleValue());
			}
		}

	}

	public void DebugPrintSimple() {
		System.out.println("----------------SpectrumHMM---------------");
		System.out.println("Num States: " + this.NumStates
				+ " Num Match States: " + this.NumMatchStates);
		for (int i = 0; i < this.NumMatchStates; ++i) {
			System.out.println("Match State[" + i + "] : "
					+ this.EmissionParams[i].GetMean() + " "
					+ this.OverlapCount[i][0] + "/" + this.OverlapCount[i][1]
					+ " Score: " + this.ComputeMatchStateSumScore(i));
		}
	}

	public void DebugPrintTransitions() {
		for (int i = 0; i < this.NumStates; ++i) {
			int[] temp = { i };
			String CurrString = this.GetStateString(temp) + " to ";
			for (int j = 0; j < this.NumStates; ++j) {
				int[] temp2 = { j };
				String NewString = CurrString + this.GetStateString(temp2)
						+ " : " + this.Transitions[i][j];
				if (this.Transitions[i][j] != 0.0)
					System.out.println(NewString);
			}
			System.out.println("");
		}
	}

	/*
	 * protected ArrayList GetBestGroupOld(int MatchState) { ArrayList
	 * CurrMasses = this.InsertMasses[MatchState]; ArrayList CurrBestGroup =
	 * null; double CurrBestGroupScore = 0.0;
	 * 
	 * for(int MainMassIndex = 0; MainMassIndex < CurrMasses.size();
	 * ++MainMassIndex) { int MainMass =
	 * (int)((double[])(CurrMasses.get(MainMassIndex)))[0]; ArrayList CurrGroup
	 * = new ArrayList(); //double CurrGroupScore = 0.0;
	 * 
	 * //double P = 1.0; //double Q = 1.0;
	 * 
	 * if(Debug) System.out.println("Considering centeral mass " + MainMass);
	 * for(int SideMassIndex = 0; SideMassIndex < CurrMasses.size();
	 * ++SideMassIndex) { int SideMass =
	 * (int)((double[])(CurrMasses.get(SideMassIndex)))[0]; if(Debug && ! Debug)
	 * System.out.println("Considering side mass " + SideMass);
	 * if(Math.abs(MainMass - SideMass) <=
	 * AntibodyUtils.LTQFTIonTolerance*AntibodyUtils.MASS_SCALE/2) { if(Debug)
	 * System.out.println("CLose enough to be in same grouP!!");
	 * CurrGroup.add(CurrMasses.get(SideMassIndex)); //CurrGroupScore +=
	 * ((double[])(CurrMasses.get(SideMassIndex)))[1]; } }
	 * 
	 * int NumOverlapping = this.InsertOverlapCount[MatchState]; int
	 * NumSupporting = CurrGroup.size();
	 * 
	 * double P = 1; double Q = 1; for(int j = 0; j < CurrGroup.size(); ++j) {
	 * double[] CurrPeak = (double[])(CurrGroup.get(j)); int ScoreZone =
	 * AntibodyUtils.DeterminePRMScoreZone(CurrPeak[1]); P *=
	 * AntibodyUtils.pA[ScoreZone]; Q *= AntibodyUtils.pNA[ScoreZone]; }
	 * 
	 * P *= Math.pow(AntibodyUtils.NotpA, NumOverlapping - NumSupporting); Q *=
	 * Math.pow(AntibodyUtils.NotpNA, NumOverlapping - NumSupporting);
	 * 
	 * double CurrGroupScore = AntibodyUtils.Log10(P/Q); if(CurrGroupScore >
	 * CurrBestGroupScore && CurrGroup.size() >=
	 * AntibodyUtils.MIN_PEAK_GROUP_SIZE) { CurrBestGroupScore = CurrGroupScore;
	 * CurrBestGroup = CurrGroup; //System.out.println("Current Best Group: " +
	 * MainMass + "!!"); }
	 * 
	 * } return CurrBestGroup; }
	 */

	protected boolean GetLikeliestPath(HMMAlignment NewAlignment,
			double Tolerance) {

		boolean LocalDebug = false;
		int[] PRMPeaks = NewAlignment.GetPRMSpectrum().PRMs;

		// TO DO: create a map to add the MassShift to all peaks
		int[] NewPeaks = NewAlignment.GetShiftedPRMs();
		double[] PRMScores = NewAlignment.GetShiftedPRMScores();
		if (Debug || LocalDebug) {
			String CurrString = "";
			String CurrString2 = "";
			for (int i = 0; i < NewPeaks.length; ++i) {
				CurrString += NewPeaks[i] + " ";
				CurrString2 += PRMPeaks[i] + " ";
			}
			System.out.println("Shifted PRMs: " + CurrString);
			System.out.println("Orig PRMs: " + CurrString2);
			System.out.println("Last match state: "
					+ this.EmissionParams[this.NumMatchStates - 1].GetMean());
			Utils.WaitForEnter();
			// this.DebugPrintTransitions();
			// Utils.WaitForEnter();
		}

		double MaxProb;
		ArrayList MaxPath = new ArrayList();
		int StartPeak = 0;

		// CurrProbs[i][j] is the probability of a path ending with PRMPeak[i]
		// and state j
		double[][] CurrProbs = new double[NewPeaks.length][this.NumStates];

		// CurrPaths[i][j] is a most likely path ending at PRMPeak[i] and state
		// j
		ArrayList[][] CurrPaths = new ArrayList[NewPeaks.length][this.NumStates];

		// Initialize the arrays
		for (int i = 0; i < NewPeaks.length; ++i) {
			for (int j = 0; j < this.NumStates; ++j) {
				CurrProbs[i][j] = Double.NEGATIVE_INFINITY;
				CurrPaths[i][j] = new ArrayList();
			}
		}

		// Set up initial paths, probabilities of paths of length 0, or through
		// delete states
		// Since delete states are silent, an alignment could begin with a bunch
		// of deletes
		for (int State = 0; State < this.NumStates; ++State) {
			if (IsDeleteState(State)) {
				MaxProb = Double.NEGATIVE_INFINITY;

				for (int PrevState = 0; PrevState < State; ++PrevState) {
					// if(Debug || LocalDebug)
					// {
					// int[] temp1 = {PrevState};
					// int[] temp2 = {State};
					// System.out.println("LogPrior of " +
					// this.GetStateString(temp1) + " to " +
					// this.GetStateString(temp2) + ":" +
					// this.LogTransitions[PrevState][State]);
					// }
					double P = CurrProbs[StartPeak][PrevState]
							+ this.LogTransitions[PrevState][State];
					if (P > MaxProb) {
						MaxProb = P;
						MaxPath = AntibodyUtils
								.CopyIntArrayList(CurrPaths[StartPeak][PrevState]);
						MaxPath.add(new Integer(State));
					}
				}
				CurrPaths[StartPeak][State] = MaxPath;
				CurrProbs[StartPeak][State] = MaxProb;
				// if(Debug || LocalDebug)
				// {
				// int[] temp = {State};
				// System.out.println("Max path to peak " + StartPeak +
				// " and State " + this.GetStateString(temp) + ":" +
				// CurrProbs[StartPeak][State]);
				// System.out.println(this.GetStateString(Utils.ConvertIntegerArrayList(CurrPaths[StartPeak][State])));
				// }
			} else if (IsMatchState(State) || IsInsertState(State)) {
				CurrProbs[StartPeak][State] = this.LogPriors[State]
						+ GetEmissionLogProb(NewPeaks[StartPeak],
								PRMScores[StartPeak], State);
				CurrPaths[StartPeak][State].clear();
				CurrPaths[StartPeak][State].add(new Integer(State));
				// if(Debug | LocalDebug)
				// {
				// System.out.println("LogPrior: " + this.LogPriors[State] +
				// " emission: " +
				// GetEmissionLogProb(NewPeaks[StartPeak],PRMScores[StartPeak],State));
				// int[] temp = {State};
				// System.out.println("Max path to peak " + StartPeak +
				// " and State " + this.GetStateString(temp) + ":" +
				// CurrProbs[StartPeak][State]);
				// System.out.println(this.GetStateString(Utils.ConvertIntegerArrayList(CurrPaths[StartPeak][State])));
				// }

			}
		}

		// if(LocalDebug)
		// Utils.WaitForEnter();
		// Now compute the rest of the table
		int LastPeak = 0;
		for (int PeakIndex = StartPeak + 1; PeakIndex < NewPeaks.length; ++PeakIndex) {
			LastPeak = PeakIndex;
			int ObservedPeak = NewPeaks[PeakIndex];
			if (ObservedPeak > this.EmissionParams[this.NumMatchStates - 1]
					.GetMean() + Tolerance * AntibodyUtils.MASS_SCALE) {
				// System.err
				// .println("WARNING: Skipping Mass shifted PRMPeaks greather than model!! "
				// + ObservedPeak);
				LastPeak = PeakIndex - 1;
				break;
			}
			// if(ObservedPeak < 0)
			// {
			// System.err.println("ERROR: Mass shifted PRMPeaks are less than 0!! "
			// + ObservedPeak);
			// return false;
			// }

			double MaxPeakProb = Double.NEGATIVE_INFINITY;
			ArrayList MaxPeakPath = null;// = new ArrayList();
			for (int State = 0; State < this.NumStates; ++State) {
				MaxPath = new ArrayList();
				MaxProb = Double.NEGATIVE_INFINITY;

				if (IsMatchState(State) || IsInsertState(State)) {
					for (int PrevState = 0; PrevState < this.NumStates; ++PrevState) {
						double TempProb = CurrProbs[PeakIndex - 1][PrevState];
						ArrayList TempPath = CurrPaths[PeakIndex - 1][PrevState];

						TempProb += GetEmissionLogProb(ObservedPeak,
								PRMScores[PeakIndex], State)
								+ this.LogTransitions[PrevState][State];
						if (TempProb > MaxProb) {
							MaxProb = TempProb;
							MaxPath = AntibodyUtils.CopyIntArrayList(TempPath);
							MaxPath.add(new Integer(State));
						}
					}
				} else if (IsDeleteState(State)) {
					for (int PrevState = 0; PrevState < State; PrevState++) {
						double TempProb = CurrProbs[PeakIndex][PrevState];
						ArrayList TempPath = CurrPaths[PeakIndex][PrevState];
						TempProb += this.LogTransitions[PrevState][State];
						if (TempProb > MaxProb) {
							MaxProb = TempProb;
							MaxPath = AntibodyUtils.CopyIntArrayList(TempPath);
							MaxPath.add(new Integer(State));
						}
					}
				}
				/*
				 * if(LocalDebug) { System.out.println("Best probl for peak[" +
				 * PeakIndex + "] and state[" + State + "]: " + MaxProb); String
				 * T = ""; for(int p = 0; p < MaxPath.size(); ++p) T +=
				 * ((Integer)(MaxPath.get(p))).intValue() + " ";
				 * System.out.println("Best Path: " + T);
				 * //Utils.WaitForEnter(); }
				 */

				CurrProbs[PeakIndex][State] = MaxProb;
				CurrPaths[PeakIndex][State] = MaxPath;
				if (MaxProb > MaxPeakProb) {
					MaxPeakProb = MaxProb;
					MaxPeakPath = MaxPath;
				}

			}
			if (LocalDebug) {
				System.out.println("Best probl for peak[" + PeakIndex + "] = "
						+ NewPeaks[PeakIndex] + " : " + MaxPeakProb);
				String T = "";
				for (int p = 0; p < MaxPeakPath.size(); ++p)
					T += ((Integer) (MaxPeakPath.get(p))).intValue() + " ";
				System.out.println("Best Path: " + T);
				Utils.WaitForEnter();
			}
		}
		MaxProb = Double.NEGATIVE_INFINITY;
		for (int State = 0; State < this.NumStates; ++State) {
			double TempProb = CurrProbs[LastPeak][State];
			if (TempProb > MaxProb) {
				MaxProb = TempProb;
				MaxPath = CurrPaths[LastPeak][State];
			}
		}

		// String PathString =
		// this.GetStateString(Utils.ConvertIntegerArrayList(MaxPath));
		// int NumMatches = Utils.CountChar(PathString,'s');
		// int NumInserts = Utils.CountChar(PathString,'i');
		// int NumDeletes = Utils.CountChar(PathString,'d');
		// double NullProb =
		// Utils.Log10(Math.pow(Utils.pNA,NumMatches+NumInserts)*Math.pow(1-Utils.pNA,NumDeletes));
		NewAlignment.UpdateAlignment(
				AntibodyUtils.ConvertIntegerArrayList(MaxPath), MaxProb, this,
				Tolerance);
		return true;
	}

	public ArrayList[] GetMatchCandidates(double ThresholdScore,
			double Tolerance) {
		Debug = true;
		ArrayList[] Candidates = new ArrayList[this.NumMatchStates];
		for (int MatchState = 0; MatchState < this.NumMatchStates; ++MatchState) {

			if (MatchState > 0) {
				// If these are internal nodes
				double PrevMass = this.EmissionParams[MatchState - 1].GetMean();
				double NextMass = this.EmissionParams[MatchState].GetMean();

				if (AntibodyUtils.ContainsMass(this.SequencePRM,
						(int) PrevMass, Tolerance)
						&& AntibodyUtils.ContainsMass(this.SequencePRM,
								(int) NextMass, Tolerance)) {
					double MassDiff = NextMass - PrevMass;
					if (AntibodyUtils.FindClosestAAScaled((int) MassDiff,
							Tolerance) != '0') // If
						// the
						// different
						// is
						// an
						// AA
						// then
						// no
						// inserts
						continue;
					if (Debug)
						System.out.println("!!!Can insert between " + PrevMass
								+ " and " + NextMass + "!!");

				}
			}

			if (Debug)
				System.out.println("Looking for best insert after state s"
						+ MatchState);
			ArrayList BestGroup = this.GetBestGroup(MatchState, Tolerance);
			Candidates[MatchState] = null;
			if (BestGroup != null) {
				// int SumMass = 0;
				double SumScore = 0.0;
				for (int i = 0; i < BestGroup.size(); ++i) {
					// SumMass += (int)(((double[])(BestGroup.get(i)))[0]);
					SumScore += (((double[]) (BestGroup.get(i)))[1]);
				}
				if (SumScore >= ThresholdScore)
					Candidates[MatchState] = BestGroup;
				else if (Debug)
					System.out.println("BestGroup is not good enough =(");

			}
		}
		Debug = false;
		return Candidates;
	}

	public ArrayList[] GetMatchCandidatesOddsScore(double ThresholdScore,
			double Tolerance) {

		boolean LocalDebug = false;
		ArrayList[] Candidates = new ArrayList[this.NumMatchStates];
		for (int MatchState = 0; MatchState < this.NumMatchStates; ++MatchState) {
			// For now, we don't allow insertions in the sequence
			// Do not allow inserts into the internals of the sequence.
			/*
			 * if(MatchState > this.NumMatchStates - this.SequencePRM.length) {
			 * Candidates[MatchState] = null; continue; }
			 */

			// Allow inserts into the internals of the sequence, IF the abutting
			// matches are a distance of two amino acids away.
			// i.e. ABCD[144.1]GHI, can have internal states inserted between
			// match states 5 and 6
			if (MatchState > 0) {
				// If these are internal nodes
				double PrevMass = this.EmissionParams[MatchState - 1].GetMean();
				double NextMass = this.EmissionParams[MatchState].GetMean();

				if (AntibodyUtils.ContainsMass(this.SequencePRM,
						(int) PrevMass, Tolerance)
						&& AntibodyUtils.ContainsMass(this.SequencePRM,
								(int) NextMass, Tolerance)) {
					double MassDiff = NextMass - PrevMass;
					if (AntibodyUtils.FindClosestAAScaled((int) MassDiff,
							Tolerance) != '0') // If
					// the
					// different
					// is
					// an
					// AA
					// then
					// no
					// inserts
					{

						if (LocalDebug)
							System.out.println("!!!Cannot insert between "
									+ PrevMass + " and " + NextMass + "!!!");
						continue;
					}
					if (LocalDebug)
						System.out.println("!!!Can insert between " + PrevMass
								+ " and " + NextMass + "!!");

				}
			}

			if (LocalDebug) {
				System.out.println("Looking for best insert before state s"
						+ MatchState);
				this.DebugPrintInsertMasses(MatchState);
			}
			// this.DebugPrintInsertMasses(MatchState);
			ArrayList BestGroup = this.GetBestGroup(MatchState, Tolerance);
			Candidates[MatchState] = null;
			if (BestGroup != null) {
				int NumOverlapping = this.InsertOverlapCount[MatchState];
				int NumSupporting = BestGroup.size();

				double P = 1;
				double Q = 1;
				for (int j = 0; j < BestGroup.size(); ++j) {
					double[] CurrPeak = (double[]) (BestGroup.get(j));
					int ScoreZone = AntibodyUtils
							.DeterminePRMScoreZone(CurrPeak[1]);
					P *= AntibodyUtils.pA[ScoreZone];
					Q *= AntibodyUtils.pNA[ScoreZone];
				}

				P *= Math.pow(AntibodyUtils.NotpA, NumOverlapping
						- NumSupporting);
				Q *= Math.pow(AntibodyUtils.NotpNA, NumOverlapping
						- NumSupporting);

				double SumScore = AntibodyUtils.Log10(P / Q);

				// double SumScore = ((double)NumSupporting)/NumOverlapping;
				if (LocalDebug) {
					System.out.println("Considering Match State: "
							+ ((double[]) (BestGroup.get(0)))[0] + " size: "
							+ BestGroup.size() + " overlap: " + NumOverlapping);

					System.out.println("Score: " + SumScore);
					System.out.println("Threshold: " + ThresholdScore);
				}
				if (SumScore >= ThresholdScore)
					Candidates[MatchState] = BestGroup;
				else if (Debug || LocalDebug)
					System.out.println("BestGroup is not good enough =(");
			}
		}
		return Candidates;
	}

	/**
	 * Given a particular mass, determine how many spectra at this insert state
	 * overlap the mass. It is used to determine the penalty given for a new
	 * match state.
	 * 
	 * @param StateNum
	 * @param Mass
	 * @return
	 */
	public int GetNumOverlappingInsert(int StateNum, int Mass, double Tolerance) {
		ArrayList CurrMasses = this.InsertMasses[StateNum];
		int Count = 0;
		Hashtable SpectrumHash = new Hashtable();
		for (int i = 0; i < CurrMasses.size(); ++i) {
			Object[] Peak = (Object[]) (CurrMasses.get(i));
			int PeakMass = ((Integer) (Peak[1])).intValue();
			String Hashkey = (String) (Peak[0]);
			// Since we are doing left extension, we know we should have atleast
			// one match on the right
			// ASSUMPTION: No insertions internal to the anchor
			if (PeakMass <= Mass - Tolerance * AntibodyUtils.MASS_SCALE) {
				if (!SpectrumHash.containsKey(Hashkey))
					Count += 1;
				SpectrumHash.put(Hashkey, new Integer(1));
			}
		}
		return Count;

	}

	public ArrayList[] GetMatchCandidatesSumScore(double ThresholdScore,
			double Tolerance) {

		boolean LocalDebug = false;
		ArrayList[] Candidates = new ArrayList[this.NumMatchStates];
		for (int MatchState = 0; MatchState < this.NumMatchStates; ++MatchState) {
			// For now, we don't allow insertions in the sequence
			// Do not allow inserts into the internals of the sequence.
			/*
			 * if(MatchState > this.NumMatchStates - this.SequencePRM.length) {
			 * Candidates[MatchState] = null; continue; }
			 */

			// Allow inserts into the internals of the sequence, IF the abutting
			// matches are a distance of two amino acids away.
			// i.e. ABCD[144.1]GHI, can have internal states inserted between
			// match states 5 and 6

			if (MatchState > 0) {
				// If these are internal nodes
				double PrevMass = this.EmissionParams[MatchState - 1].GetMean();
				double NextMass = this.EmissionParams[MatchState].GetMean();

				// if(AntibodyUtils.ContainsMass(this.SequencePRM,
				// (int)PrevMass) &&
				// AntibodyUtils.ContainsMass(this.SequencePRM, (int)NextMass))
				if (MatchState >= this.NumMatchStates - this.SequencePRM.length) {
					double MassDiff = NextMass - PrevMass;
					if (AntibodyUtils.FindClosestAAScaled((int) MassDiff,
							Tolerance) != '0') // If
					// the
					// different
					// is
					// an
					// AA
					// then
					// no
					// inserts
					{

						if (LocalDebug)
							System.out.println("!!!Cannot insert between "
									+ PrevMass + " and " + NextMass + "!!!");
						continue;
					}
					if (LocalDebug)
						System.out.println("!!!Can insert between " + PrevMass
								+ " and " + NextMass + "!!");

				}
			}
			if (LocalDebug) {
				System.out.println("Looking for best insert before state s"
						+ MatchState);
				this.DebugPrintInsertMasses(MatchState);
			}
			// this.DebugPrintInsertMasses(MatchState);
			ArrayList BestGroup = this.GetBestGroup(MatchState, Tolerance);
			Candidates[MatchState] = null;
			if (BestGroup != null) {
				// int NumOverlapping = this.InsertOverlapCount[MatchState];
				int NumSupporting = BestGroup.size();

				double SumScore = 0.0;
				int MeanMass = 0;
				for (int j = 0; j < BestGroup.size(); ++j) {
					Object[] CurrPeak = (Object[]) (BestGroup.get(j));
					SumScore += ((Double) (CurrPeak[2])).doubleValue();
					MeanMass += ((Integer) (CurrPeak[1])).intValue();
				}
				MeanMass /= BestGroup.size();
				int NumOverlapping = this.GetNumOverlappingInsert(MatchState,
						MeanMass, Tolerance);
				if (this.applyMatchPenalty)
					SumScore -= (NumOverlapping - NumSupporting)
							* this.matchPenalty;
				// double SumScore = ((double)NumSupporting)/NumOverlapping;
				if (LocalDebug) {
					Object[] CurrPeak = (Object[]) (BestGroup.get(0));
					System.out.println("Considering Match State: "
							+ ((Integer) (CurrPeak[1])).intValue() + " size: "
							+ BestGroup.size());
					System.out.println("Overlapping: " + NumOverlapping
							+ ", Supporting: " + NumSupporting);
					System.out.println("Score: " + SumScore);
					System.out.println("Threshold: " + ThresholdScore);
					Utils.WaitForEnter();
				}
				if (BestGroup.size() >= AntibodyUtils.MIN_PEAK_GROUP_SIZE
						&& SumScore >= ThresholdScore)
					Candidates[MatchState] = BestGroup;
				else if (Debug || LocalDebug)
					System.out.println("BestGroup is not good enough =(");
			}
		}
		return Candidates;
	}

	public int GetNextDelete(int stateNum) {
		if (IsMatchState(stateNum) && stateNum != this.NumMatchStates - 1)
			return stateNum + 2 * this.NumMatchStates + 1;
		if (IsInsertState(stateNum))
			return stateNum + this.NumMatchStates;
		if (IsDeleteState(stateNum) && stateNum != this.NumStates - 1)
			return stateNum + 1;
		return -1;
	}

	public int GetNextInsert(int stateNum) {
		if (IsMatchState(stateNum) && stateNum < this.NumMatchStates - 1)
			return stateNum + this.NumMatchStates + 1;
		if (IsInsertState(stateNum))
			return stateNum;
		if (IsDeleteState(stateNum) && stateNum < this.NumStates - 1)
			return stateNum - this.NumMatchStates + 1;
		return -1;
	}

	public int GetNextMatch(int stateNum) {
		if (IsMatchState(stateNum) && stateNum != this.NumMatchStates - 1)
			return stateNum + 1;
		if (IsInsertState(stateNum))
			return stateNum - this.NumMatchStates;
		if (IsDeleteState(stateNum) && stateNum != this.NumStates - 1)
			return stateNum - 2 * this.NumMatchStates + 1;
		return -1;
	}

	public int GetNumOverlapping(int State) {
		if (this.IsMatchState(State))
			return this.OverlapCount[State][1];
		else if (this.IsDeleteState(State))
			return this.OverlapCount[State - 2 * this.NumMatchStates][1];
		else if (this.IsInsertState(State))
			return this.InsertOverlapCount[State];
		return -1;
	}

	public int GetNumSupporting(int State) {
		if (this.IsMatchState(State))
			return this.OverlapCount[State][0];
		else if (this.IsDeleteState(State))
			return this.OverlapCount[State - 2 * this.NumMatchStates][1]
					- this.OverlapCount[State - 2 * this.NumMatchStates][0];
		return -1;
	}

	public boolean IsDeleteState(int stateNum) {
		if (stateNum >= 2 * this.NumMatchStates && stateNum < this.NumStates)
			return true;
		return false;
	}

	public boolean IsInsertState(int stateNum) {
		if (stateNum >= this.NumMatchStates
				&& stateNum < 2 * this.NumMatchStates)
			return true;
		return false;

	}

	public boolean IsMatchState(int stateNum) {
		if (stateNum >= 0 && stateNum < this.NumMatchStates)
			return true;
		return false;
	}

	public int[] ProduceConsensus() {
		boolean LocalDebug = Debug || false;
		System.out.println("Finding Consensus...");
		double[] BestScores = new double[this.NumStates];
		ArrayList[] BestPaths = new ArrayList[this.NumStates];

		BestScores[0] = AntibodyUtils
				.Log10(((double) AntibodyUtils.PSEUDOCOUNT_MATCH)
						/ (AntibodyUtils.PSEUDOCOUNT_MATCH + AntibodyUtils.PSEUDOCOUNT_DELETE));
		// BestScores[0] = this.ComputeMatchStateLogOddsScore(0);
		BestPaths[0] = new ArrayList();
		BestPaths[0].add(new Integer(0));

		BestScores[2 * this.NumMatchStates] = AntibodyUtils
				.Log10(((double) AntibodyUtils.PSEUDOCOUNT_DELETE)
						/ (AntibodyUtils.PSEUDOCOUNT_MATCH + AntibodyUtils.PSEUDOCOUNT_DELETE));
		// BestScores[2*this.NumMatchStates] =
		// this.ComputeDeleteStateLogOddsScore(0);
		BestPaths[2 * this.NumMatchStates] = new ArrayList();
		BestPaths[2 * this.NumMatchStates].add(new Integer(
				2 * this.NumMatchStates));

		for (int i = 1; i < this.NumMatchStates; ++i) {
			BestScores[i] = Double.NEGATIVE_INFINITY;
			BestPaths[i] = new ArrayList();
			BestPaths[i].add(new Integer(i));

			BestScores[i + 2 * this.NumMatchStates] = Double.NEGATIVE_INFINITY;
			BestPaths[i + 2 * this.NumMatchStates] = new ArrayList();
			BestPaths[i + 2 * this.NumMatchStates].add(new Integer(i + 2
					* this.NumMatchStates));
		}

		for (int t = 1; t < this.NumMatchStates; ++t) {
			if (LocalDebug)
				System.out.println("Time: " + t);
			double[] NewScores = new double[this.NumStates];
			ArrayList[] NewPaths = new ArrayList[this.NumStates];

			for (int EndState = 0; EndState < this.NumMatchStates; ++EndState) {
				if (LocalDebug)
					System.out.println("Considering end point " + EndState);
				double CurrMatchScore = Double.NEGATIVE_INFINITY; // best path
				// score
				// ending in
				// match
				// state E
				double CurrDeleteScore = Double.NEGATIVE_INFINITY; // best path
				// score
				// ending in
				// delete
				// state
				// 2*self.NumMatchStates
				// + E

				ArrayList CurrMatchPath = new ArrayList();
				ArrayList CurrDeletePath = new ArrayList();

				for (int StartState = 0; StartState < this.NumMatchStates; ++StartState) {
					if (LocalDebug && !LocalDebug) {
						System.out.println("CurrBestMatchScore: "
								+ CurrMatchScore);
						System.out.println("MatchPath: " + CurrMatchPath);
						System.out.println("CurrBestDeleteScore: "
								+ CurrDeleteScore);
						System.out.println("DeletePath: " + CurrDeletePath);
						Utils.WaitForEnter();

						System.out.println("\nNow considering prev State "
								+ StartState);

					}
					// If the last transition is a match to a match
					double TempScoreM2M = BestScores[StartState]
							+ this.LogTransitions[StartState][EndState];
					if (TempScoreM2M > CurrMatchScore) {
						CurrMatchScore = TempScoreM2M;
						CurrMatchPath = AntibodyUtils
								.CopyIntArrayList(BestPaths[StartState]);
						CurrMatchPath.add(new Integer(EndState));
					}
					// If the last transition is a delete to a match
					double TempScoreD2M = BestScores[2 * this.NumMatchStates
							+ StartState]
							+ this.LogTransitions[2 * this.NumMatchStates
									+ StartState][EndState];
					if (TempScoreD2M > CurrMatchScore) {
						CurrMatchScore = TempScoreD2M;
						CurrMatchPath = AntibodyUtils
								.CopyIntArrayList(BestPaths[2
										* this.NumMatchStates + StartState]);
						CurrMatchPath.add(new Integer(EndState));
					}
					// If the last transitions is a match to a delete
					double TempScoreM2D = BestScores[StartState]
							+ this.LogTransitions[StartState][2
									* this.NumMatchStates + EndState];
					if (TempScoreM2D > CurrDeleteScore) {
						CurrDeleteScore = TempScoreM2D;
						CurrDeletePath = AntibodyUtils
								.CopyIntArrayList(BestPaths[StartState]);
						CurrDeletePath.add(new Integer(2 * this.NumMatchStates
								+ EndState));
					}
					// If the last transitions is a delete to a delete
					double TempScoreD2D = BestScores[2 * this.NumMatchStates
							+ StartState]
							+ this.LogTransitions[2 * this.NumMatchStates
									+ StartState][2 * this.NumMatchStates
									+ EndState];
					if (TempScoreD2D > CurrDeleteScore) {
						CurrDeleteScore = TempScoreD2D;
						CurrDeletePath = AntibodyUtils
								.CopyIntArrayList(BestPaths[2
										* this.NumMatchStates + StartState]);
						CurrDeletePath.add(new Integer(2 * this.NumMatchStates
								+ EndState));
					}

					if (LocalDebug & !LocalDebug) {
						System.out.println("Score of PrevState M: "
								+ BestScores[StartState]);
						System.out.println("Score of PrevState D: "
								+ BestScores[2 * this.NumMatchStates
										+ StartState]);

						System.out.println("Transition M2M: "
								+ this.LogTransitions[StartState][EndState]);
						System.out.println("Transition M2D: "
								+ this.LogTransitions[StartState][2
										* this.NumMatchStates + EndState]);
						System.out.println("Transition D2M: "
								+ this.LogTransitions[2 * this.NumMatchStates
										+ StartState][EndState]);
						System.out.println("Transition D2D: "
								+ this.LogTransitions[2 * this.NumMatchStates
										+ StartState][2 * this.NumMatchStates
										+ EndState]);

						System.out.println("M2M: " + TempScoreM2M);
						System.out.println("M2D: " + TempScoreM2D);
						System.out.println("D2M: " + TempScoreD2M);
						System.out.println("D2D: " + TempScoreD2D);
						Utils.WaitForEnter();
					}
				}
				NewScores[EndState] = CurrMatchScore;
				NewPaths[EndState] = CurrMatchPath;
				NewScores[2 * this.NumMatchStates + EndState] = CurrDeleteScore;
				NewPaths[2 * this.NumMatchStates + EndState] = CurrDeletePath;
			}

			BestScores = NewScores;
			BestPaths = NewPaths;
			if (LocalDebug && !LocalDebug) {
				System.out.println("Current BestScores:");
				for (int k = 0; k < this.NumMatchStates; ++k) {
					System.out.println("s" + k + " : " + BestScores[k]);
					System.out.println("d" + k + " : "
							+ BestScores[2 * this.NumMatchStates + k]);
				}
				Utils.WaitForEnter();
			}
		}

		ArrayList PRMMasses = new ArrayList();
		double ConsensusScore;
		ArrayList FinalPath;
		if (BestScores[this.NumMatchStates - 1] > BestScores[3 * this.NumMatchStates - 1]) {
			ConsensusScore = BestScores[this.NumMatchStates - 1];
			FinalPath = BestPaths[this.NumMatchStates - 1];
		} else {
			ConsensusScore = BestScores[3 * this.NumMatchStates - 1];
			FinalPath = BestPaths[3 * this.NumMatchStates - 1];
		}

		if (LocalDebug) {
			System.out.println("Consensus Score: " + ConsensusScore);
		}
		for (int i = 0; i < FinalPath.size(); ++i) {

			int State = ((Integer) (FinalPath.get(i))).intValue();
			if (LocalDebug) {
				System.out.println("FInalPath[" + i + "] = " + State);
			}
			if (this.IsMatchState(State))
				PRMMasses.add(new Integer((int) (this.EmissionParams[State]
						.GetMean())));
		}
		int[] FinalMasses = AntibodyUtils.ConvertIntegerArrayList(PRMMasses);

		return FinalMasses;
	}

	public ConsensusPRMSpectrum ProduceConsensusWithScores() {

		boolean LocalDebug = Debug || false;
		if (Debug)
			System.out.println("Finding Consensus...");
		double[] BestScores = new double[this.NumStates];
		ArrayList[] BestPaths = new ArrayList[this.NumStates];

		BestScores[0] = AntibodyUtils
				.Log10(((double) AntibodyUtils.PSEUDOCOUNT_MATCH)
						/ (AntibodyUtils.PSEUDOCOUNT_MATCH + AntibodyUtils.PSEUDOCOUNT_DELETE));
		// BestScores[0] = this.ComputeMatchStateLogOddsScore(0);
		BestPaths[0] = new ArrayList();
		BestPaths[0].add(new Integer(0));

		BestScores[2 * this.NumMatchStates] = AntibodyUtils
				.Log10(((double) AntibodyUtils.PSEUDOCOUNT_DELETE)
						/ (AntibodyUtils.PSEUDOCOUNT_MATCH + AntibodyUtils.PSEUDOCOUNT_DELETE));
		// BestScores[2*this.NumMatchStates] =
		// this.ComputeDeleteStateLogOddsScore(0);
		BestPaths[2 * this.NumMatchStates] = new ArrayList();
		BestPaths[2 * this.NumMatchStates].add(new Integer(
				2 * this.NumMatchStates));

		for (int i = 1; i < this.NumMatchStates; ++i) {
			BestScores[i] = Double.NEGATIVE_INFINITY;
			BestPaths[i] = new ArrayList();
			BestPaths[i].add(new Integer(i));

			BestScores[i + 2 * this.NumMatchStates] = Double.NEGATIVE_INFINITY;
			BestPaths[i + 2 * this.NumMatchStates] = new ArrayList();
			BestPaths[i + 2 * this.NumMatchStates].add(new Integer(i + 2
					* this.NumMatchStates));
		}

		for (int t = 1; t < this.NumMatchStates; ++t) {
			if (LocalDebug)
				System.out.println("Time: " + t);
			double[] NewScores = new double[this.NumStates];
			ArrayList[] NewPaths = new ArrayList[this.NumStates];

			for (int EndState = 0; EndState < this.NumMatchStates; ++EndState) {
				if (LocalDebug)
					System.out.println("Considering end point " + EndState);
				double CurrMatchScore = Double.NEGATIVE_INFINITY; // best path
				// score
				// ending in
				// match
				// state E
				double CurrDeleteScore = Double.NEGATIVE_INFINITY; // best path
				// score
				// ending in
				// delete
				// state
				// 2*self.NumMatchStates
				// + E

				ArrayList CurrMatchPath = new ArrayList();
				ArrayList CurrDeletePath = new ArrayList();

				for (int StartState = 0; StartState < this.NumMatchStates; ++StartState) {
					if (LocalDebug & !LocalDebug) {
						System.out.println("CurrBestMatchScore: "
								+ CurrMatchScore);
						System.out.println("MatchPath: " + CurrMatchPath);
						System.out.println("CurrBestDeleteScore: "
								+ CurrDeleteScore);
						System.out.println("DeletePath: " + CurrDeletePath);
						Utils.WaitForEnter();

						System.out.println("\nNow considering prev State "
								+ StartState);

					}
					// If the last transition is a match to a match
					double TempScoreM2M = BestScores[StartState]
							+ this.LogTransitions[StartState][EndState];
					if (TempScoreM2M > CurrMatchScore) {
						CurrMatchScore = TempScoreM2M;
						CurrMatchPath = AntibodyUtils
								.CopyIntArrayList(BestPaths[StartState]);
						CurrMatchPath.add(new Integer(EndState));
					}
					// If the last transition is a delete to a match
					double TempScoreD2M = BestScores[2 * this.NumMatchStates
							+ StartState]
							+ this.LogTransitions[2 * this.NumMatchStates
									+ StartState][EndState];
					if (TempScoreD2M > CurrMatchScore) {
						CurrMatchScore = TempScoreD2M;
						CurrMatchPath = AntibodyUtils
								.CopyIntArrayList(BestPaths[2
										* this.NumMatchStates + StartState]);
						CurrMatchPath.add(new Integer(EndState));
					}
					// If the last transitions is a match to a delete
					double TempScoreM2D = BestScores[StartState]
							+ this.LogTransitions[StartState][2
									* this.NumMatchStates + EndState];
					if (TempScoreM2D > CurrDeleteScore) {
						CurrDeleteScore = TempScoreM2D;
						CurrDeletePath = AntibodyUtils
								.CopyIntArrayList(BestPaths[StartState]);
						CurrDeletePath.add(new Integer(2 * this.NumMatchStates
								+ EndState));
					}
					// If the last transitions is a delete to a delete
					double TempScoreD2D = BestScores[2 * this.NumMatchStates
							+ StartState]
							+ this.LogTransitions[2 * this.NumMatchStates
									+ StartState][2 * this.NumMatchStates
									+ EndState];
					if (TempScoreD2D > CurrDeleteScore) {
						CurrDeleteScore = TempScoreD2D;
						CurrDeletePath = AntibodyUtils
								.CopyIntArrayList(BestPaths[2
										* this.NumMatchStates + StartState]);
						CurrDeletePath.add(new Integer(2 * this.NumMatchStates
								+ EndState));
					}

					if (LocalDebug && !LocalDebug) {
						System.out.println("Score of PrevState M: "
								+ BestScores[StartState]);
						System.out.println("Score of PrevState D: "
								+ BestScores[2 * this.NumMatchStates
										+ StartState]);

						System.out.println("Transition M2M: "
								+ this.LogTransitions[StartState][EndState]);
						System.out.println("Transition M2D: "
								+ this.LogTransitions[StartState][2
										* this.NumMatchStates + EndState]);
						System.out.println("Transition D2M: "
								+ this.LogTransitions[2 * this.NumMatchStates
										+ StartState][EndState]);
						System.out.println("Transition D2D: "
								+ this.LogTransitions[2 * this.NumMatchStates
										+ StartState][2 * this.NumMatchStates
										+ EndState]);

						System.out.println("M2M: " + TempScoreM2M);
						System.out.println("M2D: " + TempScoreM2D);
						System.out.println("D2M: " + TempScoreD2M);
						System.out.println("D2D: " + TempScoreD2D);
						Utils.WaitForEnter();
					}
				}
				NewScores[EndState] = CurrMatchScore;
				NewPaths[EndState] = CurrMatchPath;
				NewScores[2 * this.NumMatchStates + EndState] = CurrDeleteScore;
				NewPaths[2 * this.NumMatchStates + EndState] = CurrDeletePath;
			}

			BestScores = NewScores;
			BestPaths = NewPaths;
			if (LocalDebug && !LocalDebug) {
				System.out.println("Current BestScores:");
				for (int k = 0; k < this.NumMatchStates; ++k) {
					System.out.println("s" + k + " : " + BestScores[k]);
					System.out.println("d" + k + " : "
							+ BestScores[2 * this.NumMatchStates + k]);
				}
				Utils.WaitForEnter();
			}
		}

		ArrayList PRMMasses = new ArrayList();
		double ConsensusScore;
		ArrayList FinalPath;
		if (BestScores[this.NumMatchStates - 1] > BestScores[3 * this.NumMatchStates - 1]) {
			ConsensusScore = BestScores[this.NumMatchStates - 1];
			FinalPath = BestPaths[this.NumMatchStates - 1];
		} else {
			ConsensusScore = BestScores[3 * this.NumMatchStates - 1];
			FinalPath = BestPaths[3 * this.NumMatchStates - 1];
		}

		if (LocalDebug) {
			System.out.println("Consensus Score: " + ConsensusScore);
		}

		ArrayList PRMScores = new ArrayList();
		ArrayList StateNums = new ArrayList();
		for (int i = 0; i < FinalPath.size(); ++i) {

			int State = ((Integer) (FinalPath.get(i))).intValue();

			if (this.IsMatchState(State)) {
				// int NumSupporting = this.OverlapCount[State][0];
				// int NumOverlapping = this.OverlapCount[State][1];

				// double P = Math.pow(Utils.pA,
				// NumSupporting)*Math.pow(Utils.pA, NumOverlapping -
				// NumSupporting);
				// double Q = Math.pow(Utils.pNA,
				// NumSupporting)*Math.pow(Utils.pNA, NumOverlapping -
				// NumSupporting);
				// PRMScores.add(new Double(Utils.Log10(P/Q)));
				// PRMScores.add(new Double(((double)NumSupporting)/(Math.max(1,
				// NumOverlapping))));
				// PRMScores.add(new
				// Double(this.ComputeMatchStateLogOddsScore(State)));
				PRMScores
						.add(new Double(this.ComputeMatchStateSumScore(State)));
			}
			if (LocalDebug) {
				System.out.println("FInalPath[" + i + "] = " + State + ","
						+ (Double) (PRMScores.get(i)));
			}
			if (this.IsMatchState(State)) {
				PRMMasses.add(new Integer((int) (this.EmissionParams[State]
						.GetMean())));
				StateNums.add(new Integer(State));
			}
		}
		ConsensusPRMSpectrum Consensus = new ConsensusPRMSpectrum();
		Consensus.PeakScores = AntibodyUtils.ConvertDoubleArrayList(PRMScores);

		// LocalDebug = true;

		Consensus.ScaledPeaks = AntibodyUtils
				.ConvertIntegerArrayList(PRMMasses);
		if (LocalDebug)
			System.out.println("ConsensusPeaks: "
					+ Utils.IntArrayToString(Consensus.ScaledPeaks));
		int[] ConsensusStateNums = AntibodyUtils
				.ConvertIntegerArrayList(StateNums);
		if (LocalDebug)
			System.out.println("States: "
					+ Utils.IntArrayToString(ConsensusStateNums));
		Consensus.OverlappingSpectra = new Object[Consensus.PeakScores.length];
		Consensus.SupportingSpectra = new Object[Consensus.PeakScores.length];
		for (int i = 0; i < Consensus.OverlappingSpectra.length; ++i) {
			String[] Overlapping = new String[this.OverlapCount[ConsensusStateNums[i]][1]];
			String[] Supporting = new String[this.OverlapCount[ConsensusStateNums[i]][0]];
			int oIndex = 0;
			int sIndex = 0;
			if (LocalDebug)
				System.out.println("Considering state " + i
						+ " of the consensus which is model state "
						+ ConsensusStateNums[i]);
			for (int j = 0; j < this.AlignedSpectra.size(); ++j) {
				HMMAlignment H = (HMMAlignment) (this.AlignedSpectra.get(j));
				String SpecKey = AntibodyUtils.CreateSpectrumKey(
						H.GetPRMSpectrum().FileName,
						H.GetPRMSpectrum().SpecIndex);
				if (H.MLPathSupports(ConsensusStateNums[i])) {
					if (LocalDebug) {
						System.out.println("Spectrum " + SpecKey + ": "
								+ H.GetPathString());
						System.out.println(" supports state "
								+ ConsensusStateNums[i]);
						Utils.WaitForEnter();
					}
					Supporting[sIndex] = SpecKey;
					Overlapping[oIndex] = Supporting[sIndex];
					oIndex++;
					sIndex++;
				} else {
					int[] MLPath = H.GetMLPath();
					int StartState = MLPath[0];
					int EndState = MLPath[MLPath.length - 1];
					if (this.IsInsertState(StartState)) {
						StartState = this.GetNextMatch(StartState);
					} else if (this.IsDeleteState(StartState)) {
						System.out
								.println("WARNING: Why does our ML path start at a Delete!!?");
						Utils.WaitForEnter();
						StartState -= 2 * this.NumMatchStates;
					}
					if (this.IsInsertState(EndState)) {
						EndState = this.GetPrevMatch(EndState);
					} else if (this.IsDeleteState(EndState)) {
						System.out
								.println("WARNING: Why does our ML path end at a Delete!!?");
						Utils.WaitForEnter();
						EndState -= 2 * this.NumMatchStates;
					}
					if (StartState <= ConsensusStateNums[i]
							&& EndState >= ConsensusStateNums[i]) {
						if (LocalDebug) {
							System.out.println("Path for spectrum (" + SpecKey
									+ "): " + H.GetPathString()
									+ "\n overlaps but does not support state "
									+ ConsensusStateNums[i]);
							Utils.WaitForEnter();
						}
						Overlapping[oIndex] = SpecKey;
						oIndex++;
					} else if (LocalDebug) {
						if (LocalDebug) {
							System.out.println("Path for spectrum (" + SpecKey
									+ "): " + H.GetPathString()
									+ "\n does not overlap state "
									+ ConsensusStateNums[i]);
							Utils.WaitForEnter();
						}
					}
				}

			}

			if (sIndex != Supporting.length) {
				// System.err.println("WARNING: We only found " + sIndex +
				// " of " + Supporting.length +
				// " expected supporting spectra!!!");
				// System.err.println("MatchState = " + ConsensusStateNums[i]);
				// Utils.WaitForEnter();
				String[] newSup = new String[sIndex];
				for (int p = 0; p < sIndex; ++p)
					newSup[p] = Supporting[p];
				Supporting = newSup;
			}
			Consensus.OverlappingSpectra[i] = Overlapping;
			Consensus.SupportingSpectra[i] = Supporting;
		}

		return Consensus;
		// return FinalMasses;
	}

	public boolean UpdateModelTransitions() {

		int[][] TransitionCounts = new int[this.NumStates][this.NumStates];
		int[] TransitionStarts = new int[this.NumStates];

		for (int s = 0; s < this.AlignedSpectra.size(); ++s) {
			int[] CurrPath = ((HMMAlignment) (this.AlignedSpectra.get(s)))
					.GetMLPath();

			for (int i = 0; i < CurrPath.length; ++i) {
				if (i < CurrPath.length - 1) {
					TransitionCounts[CurrPath[i]][CurrPath[i + 1]] += 1;
					TransitionStarts[CurrPath[i]] += 1;
				}

			}
		}

		this.BuildBaseTransitions();
		for (int Start = 0; Start < this.NumStates; ++Start) {
			for (int End = 0; End < this.NumStates; ++End) {
				if (this.Transitions[Start][End] == 0.0)
					continue;
				double Scale = 1.0 / (TransitionStarts[Start] + 1.0);

				if (this.IsMatchState(End))
					this.Transitions[Start][End] = TransitionCounts[Start][End]
							+ AntibodyUtils.PSEUDOCOUNT_MATCH + Scale
							* this.Transitions[Start][End];
				else if (this.IsInsertState(End))
					this.Transitions[Start][End] = TransitionCounts[Start][End]
							+ AntibodyUtils.PSEUDOCOUNT_INSERT + Scale
							* this.Transitions[Start][End];
				else if (this.IsDeleteState(End))
					this.Transitions[Start][End] = TransitionCounts[Start][End]
							+ AntibodyUtils.PSEUDOCOUNT_DELETE + Scale
							* this.Transitions[Start][End];

				this.Transitions[Start][End] = this.Transitions[Start][End]
						/ (TransitionStarts[Start]
								+ AntibodyUtils.PSEUDOCOUNT_MATCH
								+ AntibodyUtils.PSEUDOCOUNT_INSERT
								+ AntibodyUtils.PSEUDOCOUNT_DELETE + Scale);
			}
		}

		// this.Transitions[this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;
		// this.Transitions[2*this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;
		// this.Transitions[3*this.NumMatchStates-1][2*this.NumMatchStates-1] =
		// 1.0;

		this.ConstructLogTransitions();
		return true;
	}

	public double GetEmissionLogProb(int ObservedPeak, double MassScore,
			int State) {

		if (IsInsertState(State)) {
			if (State > this.NumMatchStates) {
				DistributionType DPrev = this.EmissionParams[State
						- this.NumMatchStates - 1];

				if (ObservedPeak < DPrev.GetMean()) {
					return Double.NEGATIVE_INFINITY;
				}
			}
			if (GetNextMatch(State) != -1) {
				DistributionType DNext = this.EmissionParams[GetNextMatch(State)];
				if (ObservedPeak > DNext.GetMean()) {

					return Double.NEGATIVE_INFINITY;
				}
			}

			return -1 * MassScore;
		}
		if (!IsMatchState(State)) {
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Cannot produce emission params for a delete state!");
		}
		DistributionType D = this.EmissionParams[State];
		double P = D.ApproximatePDF(ObservedPeak);

		return AntibodyUtils.Log10(P);
	}

	public HMMAlignment AlignSpectrum(MSAlignmentType SpectrumAlignment,
			double Tolerance) {
		// Determine MassShift
		if (Debug) {
			// System.out.println("Aligning: " +
			// SpectrumAlignment.GetSequence());
			SpectrumAlignment.Print();
			String P = "";
			int[] Peaks = SpectrumAlignment.GetAlignedPeaksScaled();
			for (int i = 0; i < Peaks.length; ++i)
				P += Peaks[i] + " ";
			System.out.println("AlignedPeaks: " + P);
			System.out.println("SeqPeaks: "
					+ Utils.IntArrayToString(SpectrumAlignment
							.GetAlignedSequencePeaksScaled()));
			Utils.WaitForEnter();
		}
		Debug = false;
		String FullAlignmentSequence = SpectrumAlignment.GetAnchorSequence();
		// String AlignedSeq = SpectrumAlignment.GetSequence();

		// Check that the seed sequence is in the Anchor for the MSAlignment
		// int SeedIndex = FullAlignmentSequence.indexOf(this.Seed);
		// if(SeedIndex < 0)
		// //{
		// System.err.println("ERROR: Seed sequence " + this.Seed +
		// " is not present in MSAlignment anchor " + FullAlignmentSequence);
		// return null;
		// }

		// Check that the aligned sequence is in the anchor for the MSAligment
		// If this fails, then our MSAlignment is crap
		// int AlignedIndex = FullAlignmentSequence.indexOf(AlignedSeq);
		// if(AlignedIndex < 0)
		// {
		// System.err.println("ERROR: Aligned sequence " + AlignedSeq +
		// " is not present in MSAlignment anchor " + FullAlignmentSequence);
		// return null;
		// }

		// Check that our seed starts before the aligned seq starts in the
		// anchor sequence
		// Also check that the seed and aligned seq overlap
		// These are not errors, we just don't align this spectrum
		// if (SeedIndex > AlignedIndex || SeedIndex + Seed.length() <=
		// AlignedIndex)
		// {
		// return null;
		// }

		// The mass shift is the mass of the amino acids between the start of
		// the seed and the start of the
		// AlignedSeq
		int MassShift = SpectrumAlignment.GetScaledMassShift();
		// MassShift =
		// (int)(AntibodyUtils.MASS_SCALE*AntibodyUtils.GetSeqMass(FullAlignmentSequence.substring(SeedIndex,AlignedIndex)));

		// for(int i = SeedIndex; i < AlignedIndex; ++i)
		// {
		// MassShift += (int)(Utils.AAMasses[FullAlignmentSequence.charAt(i) -
		// 65]*Utils.MASS_SCALE);
		// }
		// MassShift -= SpectrumAlignment.GetAlignedPeaksScaled()[0];

		HMMAlignment newAlignment = new HMMAlignment(SpectrumAlignment,
				MassShift, HMMAlignment.HMMAlignmentType.LEFT_ALIGN, Tolerance);

		Debug = false;
		if (Debug) {
			// System.out.println("SeedIndex: " + SeedIndex + ",AlignedIndex: "
			// + AlignedIndex);
			System.out.println("MassShift: " + MassShift);
			Utils.WaitForEnter();
		}
		Debug = false;
		if (!GetLikeliestPath(newAlignment, Tolerance)) {
			return null;
		}

		return newAlignment;
	}

	public void DebugPrintInsertMasses(int State) {
		if (!IsMatchState(State))
			System.out.println("State " + State + " is not a match state!!");
		else {
			ArrayList Masses = this.InsertMasses[State];
			System.out.println("**InsertMasses preceeding state " + State
					+ ":(" + Masses.size() + ")");
			for (int i = 0; i < Masses.size(); ++i) {
				Object[] Mass = (Object[]) (Masses.get(i));
				System.out.println("Mass: " + ((Integer) (Mass[1])).intValue()
						+ " Score: " + ((Double) (Mass[2])).doubleValue()
						+ " Spectrum: " + (String) Mass[0]);
			}

		}
	}

	public int GetPrevMatch(int stateNum) {
		if (IsMatchState(stateNum) && stateNum > 0)
			return stateNum - 1;
		if (IsInsertState(stateNum) && (stateNum > this.NumMatchStates))
			return stateNum - this.NumMatchStates - 1;
		if (IsDeleteState(stateNum) && (stateNum - 2 * this.NumMatchStates) > 0)
			return stateNum - 2 * this.NumMatchStates - 1;
		return -1;
	}

	public int GetPrevInsert(int stateNum) {
		if (IsMatchState(stateNum))
			return stateNum + this.NumMatchStates;
		if (IsInsertState(stateNum))
			return stateNum;
		if (IsDeleteState(stateNum))
			return stateNum - this.NumMatchStates;
		return -1;
	}

	public int GetPrevDelete(int stateNum) {
		if (IsMatchState(stateNum) && stateNum > 0)
			return stateNum - 1 + 2 * this.NumMatchStates;
		if (IsInsertState(stateNum) && stateNum > this.NumMatchStates)
			return stateNum + this.NumMatchStates - 1;
		if (IsDeleteState(stateNum) && (stateNum - 2 * this.NumMatchStates) > 0)
			return stateNum - 1;
		return -1;
	}

}
