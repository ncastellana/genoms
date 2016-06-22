package antibody;

import java.util.ArrayList;
import java.util.Hashtable;

import basicUtils.Utils;

public class GatherConsensus {

	public boolean Debug = false;

	// File to write the alignment information to
	private String OutputFileName;

	// List of alignments being considered for this model building
	// Alignments here are to the anchor sequence from SpecAlign, not alignment
	// to the HMM.
	private MSAlignmentType[] AllAlignments;

	// List of alignments to the HMM, this needs to get updated
	// every time the model changes
	private ArrayList HMMAlignments;

	private SpectrumHMM Model;

	private boolean applyMatchPenalty;

	private double matchPenalty;

	public GatherConsensus(boolean applyMatchPenalty, double matchPenalty) {
		this.applyMatchPenalty = applyMatchPenalty;
		this.matchPenalty = matchPenalty;
	}

	/**
	 * Point of access to the Gather method
	 * 
	 * @param Alignments
	 * @param OutputFileName
	 * @param SeedSequence
	 * @return
	 */
	public int[] RunGatherAlign(MSAlignmentType[] Alignments,
			String OutputFileName, Seed SeqSeed, double Tolerance) {
		this.AllAlignments = Alignments;
		//String SeedSequence = SeqSeed.GetSeedSequence();
		// Model = new SpectrumHMM(SeedSequence);
		Model = new SpectrumHMMRight(SeqSeed, Tolerance, this.applyMatchPenalty, this.matchPenalty);
		if (Debug && !Debug) {
			Model.DebugPrintSimple();
			Utils.WaitForEnter();
		}
		/*
		 * For each Alignment in AllAlignments: Attempt to align to the model if
		 * score of alignment > threshold: AddAlignmentToModel If Model should
		 * be changed: 1)Change Model 2) Realign spectra 3) Readd good spectra
		 * Create a consensus
		 */

		int SkippedAtRealignment = 0;
		for (int i = 0; i < this.AllAlignments.length; ++i) {
			// If this Spectrum already is used in an HMM then skip it
			if (this.AllAlignments[i].GetSpectrumNode().GetHMMAnnotation() != null)
				continue;

			HMMAlignment NewAlignment = Model.AlignSpectrum(
					this.AllAlignments[i], Tolerance);

			if (Debug && NewAlignment != null) {
				System.out.println("Aligned "
						+ NewAlignment.GetPRMSpectrum().FileName + " "
						+ NewAlignment.GetPRMSpectrum().SpecIndex);
				this.AllAlignments[i].Print();
				System.out.println(NewAlignment.GetPathString() + "\nScore:"
						+ NewAlignment.GetPathScore());
				Utils.WaitForEnter();
			}
			if (NewAlignment != null
					&& NewAlignment.GetNumMatches() >= AntibodyUtils.MIN_PEAKS_ALIGNED) {
				// SpectrumNode ParentTemp =
				// this.AllAlignments[i].GetSpectrumNode();
				this.AllAlignments[i].GetSpectrumNode().SetHMMAnnotation(
						NewAlignment);
				((SpectrumHMMRight) Model).AddAlignmentToModel(NewAlignment);

			} else {
				// If this spectrum doesn't align well, then reset it so it can
				// be used again later
				SpectrumNode ParentTemp = this.AllAlignments[i]
						.GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				this.AllAlignments[i].SetSpectrumNode(null);
				ParentTemp.SetHMMAnnotation(null);
				if (Debug) {
					System.out.println("EXCLUDED:This alignment has less than "
							+ AntibodyUtils.MIN_PEAKS_ALIGNED + " matches");
					System.out.println(" ");
					Utils.WaitForEnter();
				}
				continue;
			}
			if (Debug) {
				Model.DebugPrintSimple();
				Utils.WaitForEnter();
			}

			double MeanMatchScore = Model.GetMeanMatchScore();
			// ArrayList[] InsertCandidates =
			// Model.GetMatchCandidatesOddsScore(MeanMatchScore);
			ArrayList[] InsertCandidates = Model.GetMatchCandidates(
					MeanMatchScore, Tolerance);

			boolean AddedState = false;
			ArrayList Spectra = Model.GetAlignedSpectra();
			if (Debug)
				System.out.println("Num spectra in alignment: "
						+ Spectra.size());
			for (int k = 0; k < InsertCandidates.length; ++k) {
				if (InsertCandidates[k] != null) {
					ArrayList CurrMasses = InsertCandidates[k];
					int[] Masses = new int[CurrMasses.size()];
					for (int j = 0; j < Masses.length; ++j) {
						Masses[j] = (int) ((double[]) (CurrMasses.get(j)))[0];
					}
					if (Debug) {
						System.out.println("Adding an insert mass near "
								+ Masses[0]);
						Utils.WaitForEnter();
					}
					Model.AddMatchState(Masses, Tolerance);
					if (Debug && !Debug) {
						// Model.DebugPrint();
						System.out.println("Added the new state!!");
						Utils.WaitForEnter();
					}
					AddedState = true;
				}
			}
			// If we changed the topology, realign the spectra
			if (AddedState) {
				if (Debug) {
					System.out.println("Holy Shit we changed teh topology!!");

					System.out.println("Num spectra in alignment: "
							+ Spectra.size());
					Utils.WaitForEnter();
				}
				HMMAlignment[] NewAlignments = new HMMAlignment[Spectra.size()];
				for (int j = 0; j < Spectra.size(); ++j) {
					HMMAlignment NextAlignment = (HMMAlignment) (Spectra.get(j));
					if (Debug)
						System.out.println("Realigning "
								+ NextAlignment.GetPRMSpectrum().FileName + ":"
								+ NextAlignment.GetPRMSpectrum().SpecIndex);
					double OldScore = NextAlignment.GetPathScore();
					Model.ReAlignSpectrum(NextAlignment, Tolerance);
					if (NextAlignment.GetNumMatches() >= AntibodyUtils.MIN_PEAKS_ALIGNED
							&& NextAlignment.GetPathScore() > OldScore
									* AntibodyUtils.HMM_SCORE_SCALE) {
						NewAlignments[j] = NextAlignment;
						// System.out.println("**New Alignment Score: " +
						// NewAlignments[j].GetPathScore());
					} else // If the realignment is poor, clear this
					// spectrumNode for later use
					{
						SkippedAtRealignment += 1;
						if (Debug) {
							System.out
									.println("**This alignment is bad, throwing it out!!");
						}
						SpectrumNode ParentTemp = NextAlignment
								.GetMSAlignmentType().GetSpectrumNode();
						ParentTemp.SetMSAlignment(null);
						NextAlignment.GetMSAlignmentType().GetSpectrumNode()
								.SetHMMAnnotation(null);
						NewAlignments[j] = null;
					}
				}
				Model.AddAlignmentsToModel(NewAlignments);
				if (Debug) {
					Model.DebugPrintSimple();
					Utils.WaitForEnter();
				}
			}
			if (Debug) {
				int[] Consensus = Model.ProduceConsensus();
				String S = "";
				for (int k = 0; k < Consensus.length; ++k) {
					S += Consensus[k] + " ";
				}
				System.out.println("Consensus: " + S);
				Utils.WaitForEnter();
			}
			// System.out.println("Skipped At Realignment: " +
			// SkippedAtRealignment);

		}
		return Model.ProduceConsensus();

	}

	public ConsensusPRMSpectrum RunGatherAlignWithScoresRight(
			MSAlignmentType[] Alignments, String OutputFileName, Seed SeqSeed,
			double Tolerance) {
		return RunGatherAlignWithScoresRight(Alignments, OutputFileName,
				SeqSeed, null, Tolerance);

	}

	public ConsensusPRMSpectrum RunGatherAlignWithScoresRight(
			MSAlignmentType[] Alignments, String OutputFileName, Seed SeqSeed,
			Hashtable Oracle, double Tolerance) {
		return RunGatherAlignWithScoresRight(Alignments, OutputFileName,
				SeqSeed, Oracle, false, false, Tolerance);
	}

	public ConsensusPRMSpectrum RunGatherAlignWithScoresRight(
			MSAlignmentType[] Alignments, String OutputFileName, Seed SeqSeed,
			Hashtable Oracle, boolean ClusterFlag, boolean AddAtEndFlag,
			double Tolerance) {
		this.AllAlignments = Alignments;
		Model = new SpectrumHMMRight(SeqSeed, Tolerance, this.applyMatchPenalty, this.matchPenalty);
		// Model = new SpectrumHMM(SeedSequence);

		boolean TestDebug = false;
		double MeanMatchScore = Model.GetMeanMatchScore();

		double[] InitialScores = new double[Alignments.length];
		if (Debug)
			System.out.println("MeanMatchScore: " + MeanMatchScore);
		// if(Model.Seed.indexOf("KGACLLPK") >= 0)
		// TestDebug = true;
		if (TestDebug) {
			Model.DebugPrintSimple();
			Utils.WaitForEnter();
		}
		/*
		 * For each Alignment in AllAlignments: Attempt to align to the model if
		 * score of alignment > threshold: AddAlignmentToModel If Model should
		 * be changed: 1)Change Model 2) Realign spectra 3) Read good spectra
		 * Create a consensus
		 */
		int SkippedAtRealignment = 0;
		for (int i = 0; i < this.AllAlignments.length; ++i) {
			if (this.AllAlignments[i] == null)
				System.out.println("SHIT!! " + i);
			// If this Spectrum already is used in an HMM then skip it
			// if (this.AllAlignments[i].GetSpectrumNode().GetHMMAnnotation() !=
			// null) {
			// System.out.println("We used to skip this one, but now its ok!");
			// this.AllAlignments[i].Print();
			// Utils.WaitForEnter();
			// }

			if (ClusterFlag) {
				boolean Fails = false;

				for (int j = 0; j < Model.AlignedSpectra.size(); ++j) {
					HMMAlignment PrevAlignment = (HMMAlignment) (Model.AlignedSpectra
							.get(j));

					double Score = AntibodyUtils
							.AlignmentScoreLengthWithMassShiftRight(
									this.AllAlignments[i],
									PrevAlignment.GetMSAlignmentType(),
									Tolerance);

					if (Score < AntibodyUtils.MSALIGN_MIN_SIMILARITY) {
						System.out
								.println("Removing this spectrum because it conflicts!!");
						System.out.println("CurrSpectrum: "
								+ this.AllAlignments[i].toStringDebug(Oracle));
						System.out.println("PrevSpectrum: "
								+ PrevAlignment.GetMSAlignmentType()
										.toStringDebug(Oracle));

						Fails = true;
						break;
					}
				}
				if (Fails)
					continue;
			}

			if (Debug || TestDebug) {
				System.out.println("Aligning next spectrum...");
				System.out.println(this.AllAlignments[i].toStringDebug());
			}

			HMMAlignment NewAlignment = Model.AlignSpectrum(
					this.AllAlignments[i], Tolerance);

			if (NewAlignment == null)
				System.out.println("HMM Alignment not possible!!");

			if (Debug && NewAlignment != null) {
				System.out.println("Aligned "
						+ NewAlignment.GetPRMSpectrum().FileName + " "
						+ NewAlignment.GetPRMSpectrum().SpecIndex);
				this.AllAlignments[i].Print();
				System.out.println(Utils.IntArrayToString(NewAlignment
						.GetShiftedPRMs()) + "\n");
				System.out.println(NewAlignment.GetPathString() + "\nScore:"
						+ NewAlignment.GetPathScore());
				System.out.println("NumMatches: "
						+ NewAlignment.GetNumMatches());
				System.out.println("NumDeletes: "
						+ NewAlignment.GetNumDeletes());
				System.out.println("RunOfDeletes: "
						+ NewAlignment.GetLargestRunDeletes());
				System.out
						.println("Percent of Peaks Deletes: "
								+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
										.GetPRMSpectrum().PRMs.length));
				System.out.println(" Path: " + NewAlignment.GetPathString());
				Utils.WaitForEnter();
			}
			// NEC CHANGE
			if (Debug)
				System.out.println("Min Match States: "
						+ (4.0 / 11 * Model.NumMatchStates));
			// if(NewAlignment != null && NewAlignment.GetNumMatches() >=
			// Utils.Round(4.0/11*Model.NumMatchStates) &&
			// NewAlignment.GetNumDeletes() < NewAlignment.GetNumMatches() &&
			// NewAlignment.GetLargestRunDeletes() <=
			// Math.min(NewAlignment.GetNumMatches(), Utils.MIN_PEAKS_ALIGNED))
			if (NewAlignment != null
					&& HMMAlignment.IsGoodHMMAlignment(NewAlignment, Model)) {
				// SpectrumNode ParentTemp =
				// this.AllAlignments[i].GetSpectrumNode();
				this.AllAlignments[i].GetSpectrumNode().SetHMMAnnotation(
						NewAlignment);
				// InitialScores[i] = NewAlignment.GetPathScore();
				((SpectrumHMMRight) Model).AddAlignmentToModel(NewAlignment);
				InitialScores[Model.AlignedSpectra.size() - 1] = NewAlignment
						.GetPathScore();
				if (Debug)
					System.out.println("InitialScores["
							+ (Model.AlignedSpectra.size() - 1) + "] = "
							+ NewAlignment.GetPathScore());
				if (TestDebug) {
					System.out.println("Aligned "
							+ NewAlignment.GetPRMSpectrum().FileName + " "
							+ NewAlignment.GetPRMSpectrum().SpecIndex);
					System.out
							.println(" Score: " + NewAlignment.GetPathScore());
					System.out.println(" NumMatches: "
							+ NewAlignment.GetNumMatches());
					System.out.println(" NumDeletes: "
							+ NewAlignment.GetNumDeletes());
					System.out.println(" RunOfDeletes: "
							+ NewAlignment.GetLargestRunDeletes());
					System.out
							.println(" Percent of Peaks Deletes: "
									+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
											.GetPRMSpectrum().PRMs.length));
					System.out
							.println(" Path: " + NewAlignment.GetPathString());
				}

			} else {
				// If this spectrum doesn't align well, then reset it so it can
				// be used again later
				SpectrumNode ParentTemp = this.AllAlignments[i]
						.GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				this.AllAlignments[i].SetSpectrumNode(null);
				ParentTemp.SetHMMAnnotation(null);

				if (TestDebug && NewAlignment != null) {
					System.out.println("Aligned "
							+ NewAlignment.GetPRMSpectrum().FileName + " "
							+ NewAlignment.GetPRMSpectrum().SpecIndex);
					System.out
							.println("**This alignment is bad, throwing it out!!");
					System.out
							.println(" Score: " + NewAlignment.GetPathScore());
					System.out.println(" NumMatches: "
							+ NewAlignment.GetNumMatches());
					System.out.println(" NumDeletes: "
							+ NewAlignment.GetNumDeletes());
					System.out.println(" RunOfDeletes: "
							+ NewAlignment.GetLargestRunDeletes());
					System.out
							.println(" Percent of Peaks Deletes: "
									+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
											.GetPRMSpectrum().PRMs.length));
					System.out
							.println(" Path: " + NewAlignment.GetPathString());
					// System.out.println(" Percent of Peaks Matched: " +
					// (((double)NewAlignment.GetNumMatches())/NewAlignment.GetPRMSpectrum().PRMs.length));

					Utils.WaitForEnter();
				}
				continue;
			}
			if (Debug || TestDebug) {
				Model.DebugPrint();
				Utils.WaitForEnter();
			}

			// ArrayList[] InsertCandidates =
			// Model.GetMatchCandidates(MeanMatchScore*2/3);
			// ArrayList[] InsertCandidates =
			// Model.GetMatchCandidatesOddsScore(MeanMatchScore*Utils.HMM_MIN_SCORE_FRAC);
			// System.out.println("Calling GetMatchCanddiates SumScore for " +
			// SeqSeed.AnchoredSeedSequence + " after adding " +
			// Model.AlignedSpectra.size() + " spectra");
			ArrayList[] InsertCandidates = Model.GetMatchCandidatesSumScore(
					AntibodyUtils.MIN_PEAK_GROUP_SCORE, Tolerance);
			// System.out.println("GATHERCONSENSUS: New state must have score >= "
			// + AntibodyUtils.MIN_PEAK_GROUP_SCORE);

			boolean AddedState = false;
			ArrayList Spectra = Model.GetAlignedSpectra();
			if (Debug || TestDebug)
				System.out.println("Num spectra in alignment: "
						+ Spectra.size());
			for (int k = 0; k < InsertCandidates.length; ++k) {
				if (InsertCandidates[k] != null) {
					// Debug = true;
					ArrayList CurrMasses = InsertCandidates[k];

					int[] Masses = new int[CurrMasses.size()];
					int Index = 0;
					double Score = 0;
					for (int j = 0; j < CurrMasses.size(); ++j) {
						if (CurrMasses.get(j) != null) {
							Object[] CurrPeak = (Object[]) (CurrMasses.get(j));
							Masses[Index] = ((Integer) (CurrPeak[1]))
									.intValue();
							Score += ((Double) (CurrPeak[2])).doubleValue();
							Index += 1;
						}
					}
					if (Score < AntibodyUtils.MIN_PEAK_GROUP_SCORE) {
						if (Debug) {
							System.out
									.println("New score is not good enough!!:"
											+ Score);
							Utils.WaitForEnter();
						}
						continue;
					}

					// System.out.println("Insert Candidate: " + Masses[0] +
					// "!!");
					if (Debug || TestDebug) {
						System.out.println("Adding an insert mass near "
								+ Masses[0] + " with " + CurrMasses.size()
								+ " supporting spectra");
						Utils.WaitForEnter();
					}
					Model.AddMatchState(Masses, Tolerance);
					if (Debug && !Debug) {
						// Model.DebugPrint();
						System.out.println("Added the new state!!");
						Utils.WaitForEnter();
					}
					AddedState = true;
				}
			}
			// If we changed the topology, realign the spectra
			if (AddedState) {
				if (Debug || TestDebug) {
					System.out.println("Holy Shit we changed teh topology!!");

					System.out.println("Num spectra in alignment: "
							+ Spectra.size());
					Utils.WaitForEnter();
				}
				HMMAlignment[] NewAlignments = new HMMAlignment[Spectra.size()];
				for (int j = 0; j < Spectra.size(); ++j) {
					HMMAlignment NextAlignment = (HMMAlignment) (Spectra.get(j));
					if (Debug || TestDebug)
						System.out.println("Realigning "
								+ NextAlignment.GetPRMSpectrum().FileName + ":"
								+ NextAlignment.GetPRMSpectrum().SpecIndex);
					double OldScore = InitialScores[j];
					String OldPath = NextAlignment.GetPathString();
					Model.ReAlignSpectrum(NextAlignment, Tolerance);
					if (Debug || TestDebug) {
						System.out.println(" Old Score: " + InitialScores[j]
								+ " -> NewScore: "
								+ NextAlignment.GetPathScore());
						System.out.println(" NumMatches: "
								+ NextAlignment.GetNumMatches());
						System.out.println(" NumDeletes: "
								+ NextAlignment.GetNumDeletes());
						System.out.println(" RunOfDeletes: "
								+ NextAlignment.GetLargestRunDeletes());
						System.out
								.println(" Percent of Peaks Deletes: "
										+ (((double) NextAlignment
												.GetNumDeletes()) / NextAlignment
												.GetPRMSpectrum().PRMs.length));
						System.out.println(Utils.IntArrayToString(NextAlignment
								.GetShiftedPRMs()) + "\n");
						System.out.println(" OldPath: " + OldPath);
						System.out.println(" Path: "
								+ NextAlignment.GetPathString());
					}
					// NEC_CHANGE!!
					// if(NextAlignment.GetNumMatches() >=
					// Utils.Round(4.0/11*Model.NumMatchStates) &&
					// NextAlignment.GetPathScore() >
					// OldScore*Utils.HMM_SCORE_SCALE &&
					// NextAlignment.GetNumDeletes() <
					// NextAlignment.GetNumMatches() &&
					// NextAlignment.GetLargestRunDeletes() <=
					// Math.min(NextAlignment.GetNumMatches(),
					// Utils.MIN_PEAKS_ALIGNED))
					if (HMMAlignment.IsGoodHMMAlignment(NextAlignment, Model)) {
						NewAlignments[j] = NextAlignment;
						// System.out.println("**New Alignment Score: " +
						// NewAlignments[j].GetPathScore());
					} else // If the realignment is poor, clear this
					// spectrumNode for later use
					{
						SkippedAtRealignment += 1;
						if (TestDebug || Debug) {
							System.out.println("Realigned "
									+ NextAlignment.GetPRMSpectrum().FileName
									+ ":"
									+ NextAlignment.GetPRMSpectrum().SpecIndex);

							System.out
									.println("**This alignment is bad, throwing it out!!");
							System.out.println(" Old Score: " + OldScore
									+ " -> NewScore: "
									+ NextAlignment.GetPathScore());
							System.out.println(" NumMatches: "
									+ NextAlignment.GetNumMatches());
							System.out.println(" NumDeletes: "
									+ NextAlignment.GetNumDeletes());
							System.out.println(" RunOfDeletes: "
									+ NextAlignment.GetLargestRunDeletes());
							System.out
									.println(" Percent of Peaks Deletes: "
											+ (((double) NextAlignment
													.GetNumDeletes()) / NextAlignment
													.GetPRMSpectrum().PRMs.length));
							Utils.WaitForEnter();
						}
						SpectrumNode ParentTemp = NextAlignment
								.GetMSAlignmentType().GetSpectrumNode();
						ParentTemp.SetMSAlignment(null);
						ParentTemp.SetHMMAnnotation(null);
						NewAlignments[j] = null;
					}
				}
				Model.AddAlignmentsToModel(NewAlignments);
				if (Debug || TestDebug) {
					Model.DebugPrintSimple();
					Utils.WaitForEnter();
				}
			}
			if (Debug) {
				int[] Consensus = Model.ProduceConsensus();
				String S = "";

				for (int k = 0; k < Consensus.length; ++k) {
					S += Consensus[k] + " ";

				}
				System.out.println("Consensus: " + S);
				Utils.WaitForEnter();
			}
			// System.out.println("Skipped At Realignment: " +
			// SkippedAtRealignment);

		}

		if (AddAtEndFlag) {
			// There might be more good states to add, but no more spectra let's
			// see
			System.out
					.println("Finished adding spectra, considering more possible states");
			while (true) {
				ArrayList[] InsertCandidates = Model
						.GetMatchCandidatesSumScore(
								AntibodyUtils.MIN_PEAK_GROUP_SCORE, Tolerance);
				// System.out.println("GATHERCONSENSUS: New state must have score >= "
				// + AntibodyUtils.MIN_PEAK_GROUP_SCORE);

				boolean AddedState = false;
				ArrayList Spectra = Model.GetAlignedSpectra();
				if (Debug || TestDebug)
					System.out.println("Num spectra in alignment: "
							+ Spectra.size());
				for (int k = 0; k < InsertCandidates.length; ++k) {
					if (InsertCandidates[k] != null) {
						// Debug = true;
						ArrayList CurrMasses = InsertCandidates[k];

						int[] Masses = new int[CurrMasses.size()];
						int Index = 0;
						double Score = 0;
						for (int j = 0; j < CurrMasses.size(); ++j) {
							if (CurrMasses.get(j) != null) {
								Object[] CurrPeak = (Object[]) (CurrMasses
										.get(j));
								Masses[Index] = ((Integer) (CurrPeak[1]))
										.intValue();
								Score += ((Double) (CurrPeak[2])).doubleValue();
								Index += 1;
							}
						}
						if (Score < AntibodyUtils.MIN_PEAK_GROUP_SCORE) {
							if (Debug) {
								System.out
										.println("New score is not good enough!!:"
												+ Score);
								Utils.WaitForEnter();
							}
							continue;
						}

						// System.out.println("Insert Candidate: " + Masses[0] +
						// "!!");
						if (Debug || TestDebug) {
							System.out.println("Adding an insert mass near "
									+ Masses[0] + " with " + CurrMasses.size()
									+ " supporting spectra");
							Utils.WaitForEnter();
						}
						Model.AddMatchState(Masses, Tolerance);
						if (Debug && !Debug) {
							// Model.DebugPrint();
							System.out.println("Added the new state!!");
							Utils.WaitForEnter();
						}
						AddedState = true;
					}
				}
				// If we changed the topology, realign the spectra
				if (AddedState) {
					if (Debug || TestDebug) {
						System.out
								.println("Holy Shit we changed teh topology!!");

						System.out.println("Num spectra in alignment: "
								+ Spectra.size());
						Utils.WaitForEnter();
					}
					HMMAlignment[] NewAlignments = new HMMAlignment[Spectra
							.size()];
					for (int j = 0; j < Spectra.size(); ++j) {
						HMMAlignment NextAlignment = (HMMAlignment) (Spectra
								.get(j));
						if (Debug || TestDebug)
							System.out.println("Realigning "
									+ NextAlignment.GetPRMSpectrum().FileName
									+ ":"
									+ NextAlignment.GetPRMSpectrum().SpecIndex);
						double OldScore = InitialScores[j];
						String OldPath = NextAlignment.GetPathString();
						Model.ReAlignSpectrum(NextAlignment, Tolerance);
						if (Debug || TestDebug) {
							System.out.println(" Old Score: "
									+ InitialScores[j] + " -> NewScore: "
									+ NextAlignment.GetPathScore());
							System.out.println(" NumMatches: "
									+ NextAlignment.GetNumMatches());
							System.out.println(" NumDeletes: "
									+ NextAlignment.GetNumDeletes());
							System.out.println(" RunOfDeletes: "
									+ NextAlignment.GetLargestRunDeletes());
							System.out
									.println(" Percent of Peaks Deletes: "
											+ (((double) NextAlignment
													.GetNumDeletes()) / NextAlignment
													.GetPRMSpectrum().PRMs.length));
							System.out.println(Utils
									.IntArrayToString(NextAlignment
											.GetShiftedPRMs())
									+ "\n");
							System.out.println(" OldPath: " + OldPath);
							System.out.println(" Path: "
									+ NextAlignment.GetPathString());
						}
						// NEC_CHANGE!!
						// if(NextAlignment.GetNumMatches() >=
						// Utils.Round(4.0/11*Model.NumMatchStates) &&
						// NextAlignment.GetPathScore() >
						// OldScore*Utils.HMM_SCORE_SCALE &&
						// NextAlignment.GetNumDeletes() <
						// NextAlignment.GetNumMatches() &&
						// NextAlignment.GetLargestRunDeletes() <=
						// Math.min(NextAlignment.GetNumMatches(),
						// Utils.MIN_PEAKS_ALIGNED))
						if (HMMAlignment.IsGoodHMMAlignment(NextAlignment,
								Model)) {
							NewAlignments[j] = NextAlignment;
							// System.out.println("**New Alignment Score: " +
							// NewAlignments[j].GetPathScore());
						} else // If the realignment is poor, clear this
						// spectrumNode for later use
						{
							SkippedAtRealignment += 1;
							if (TestDebug || Debug) {
								System.out
										.println("Realigned "
												+ NextAlignment
														.GetPRMSpectrum().FileName
												+ ":"
												+ NextAlignment
														.GetPRMSpectrum().SpecIndex);

								System.out
										.println("**This alignment is bad, throwing it out!!");
								System.out.println(" Old Score: " + OldScore
										+ " -> NewScore: "
										+ NextAlignment.GetPathScore());
								System.out.println(" NumMatches: "
										+ NextAlignment.GetNumMatches());
								System.out.println(" NumDeletes: "
										+ NextAlignment.GetNumDeletes());
								System.out.println(" RunOfDeletes: "
										+ NextAlignment.GetLargestRunDeletes());
								System.out
										.println(" Percent of Peaks Deletes: "
												+ (((double) NextAlignment
														.GetNumDeletes()) / NextAlignment
														.GetPRMSpectrum().PRMs.length));
								Utils.WaitForEnter();
							}
							SpectrumNode ParentTemp = NextAlignment
									.GetMSAlignmentType().GetSpectrumNode();
							ParentTemp.SetMSAlignment(null);
							ParentTemp.SetHMMAnnotation(null);
							NewAlignments[j] = null;
						}
					}
					Model.AddAlignmentsToModel(NewAlignments);
					if (Debug || TestDebug) {
						Model.DebugPrintSimple();
						Utils.WaitForEnter();
					}
				} else {
					if (Debug || TestDebug)
						System.out.println("No more states can be added!!!");

					break;
				}
				if (Debug) {
					int[] Consensus = Model.ProduceConsensus();
					String S = "";

					for (int k = 0; k < Consensus.length; ++k) {
						S += Consensus[k] + " ";

					}
					System.out.println("Consensus: " + S);
					Utils.WaitForEnter();
				}

			}
		}

		// Model.DebugPrintSimple();
		// ConsensusPRMSpectrum Consensus =
		// Model.ProduceConsensusWithScoresWithFilter();
		ConsensusPRMSpectrum Consensus = Model.ProduceConsensusWithScores();
		Consensus.Model = Model;
		// Consensus.IntervalScores = Model
		// .ComputeProbabilityTable(Consensus.ScaledPeaks);

		return Consensus;

	}

	public SpectrumHMM RunGatherAlignWithScoresRightKeepModel(
			MSAlignmentType[] Alignments, String OutputFileName, Seed SeqSeed,
			Hashtable Oracle, boolean ClusterFlag, boolean AddAtEndFlag,
			double Tolerance) {
		return this.RunGatherAlignWithScoresRightKeepModel(Alignments,
				OutputFileName, SeqSeed, Oracle, ClusterFlag, AddAtEndFlag,
				null, Tolerance);
	}

	public SpectrumHMM RunGatherAlignWithScoresRightKeepModel(
			MSAlignmentType[] Alignments, String OutputFileName, Seed SeqSeed,
			Hashtable Oracle, boolean ClusterFlag, boolean AddAtEndFlag,
			SpectrumHMMRight OldModel, double Tolerance) {
		this.AllAlignments = Alignments;
		if (OldModel == null)
			Model = new SpectrumHMMRight(SeqSeed, Tolerance, this.applyMatchPenalty, this.matchPenalty);
		// Model = new SpectrumHMM(SeedSequence);

		boolean TestDebug = false;
		double MeanMatchScore = Model.GetMeanMatchScore();

		double[] InitialScores = new double[Alignments.length];
		if (Debug)
			System.out.println("MeanMatchScore: " + MeanMatchScore);
		if (Model.Seed.indexOf("KGACLLPK") >= 0)
			TestDebug = true;
		if (TestDebug) {
			Model.DebugPrintSimple();
			Utils.WaitForEnter();
		}
		/*
		 * For each Alignment in AllAlignments: Attempt to align to the model if
		 * score of alignment > threshold: AddAlignmentToModel If Model should
		 * be changed: 1)Change Model 2) Realign spectra 3) Read good spectra
		 * Create a consensus
		 */

		int SkippedAtRealignment = 0;
		for (int i = 0; i < this.AllAlignments.length; ++i) {
			if (this.AllAlignments[i] == null)
				System.out.println("SHIT!! " + i);
			// If this Spectrum already is used in an HMM then skip it
			if (this.AllAlignments[i].GetSpectrumNode().GetHMMAnnotation() != null)
				continue;

			// Ensure similarity to previous guys
			if (ClusterFlag) {
				boolean Fails = false;

				for (int j = 0; j < Model.AlignedSpectra.size(); ++j) {
					HMMAlignment PrevAlignment = (HMMAlignment) (Model.AlignedSpectra
							.get(j));

					double Score = AntibodyUtils
							.AlignmentScoreLengthWithMassShiftRight(
									this.AllAlignments[i],
									PrevAlignment.GetMSAlignmentType(),
									Tolerance);

					if (Score < AntibodyUtils.MSALIGN_MIN_SIMILARITY) {
						System.out
								.println("Removing this spectrum because it conflicts!!");
						System.out.println("CurrSpectrum: "
								+ this.AllAlignments[i].toStringDebug(Oracle));
						System.out.println("PrevSpectrum: "
								+ PrevAlignment.GetMSAlignmentType()
										.toStringDebug(Oracle));

						Fails = true;
						break;
					}
				}
				if (Fails)
					continue;
			}

			if (Debug || TestDebug) {
				System.out.println("Aligning next spectrum...");
				System.out.println(this.AllAlignments[i].toStringDebug());
			}

			HMMAlignment NewAlignment = Model.AlignSpectrum(
					this.AllAlignments[i], Tolerance);

			if (NewAlignment == null)
				System.out.println("HMM Alignment not possible!!");

			if (Debug && NewAlignment != null) {
				System.out.println("Aligned "
						+ NewAlignment.GetPRMSpectrum().FileName + " "
						+ NewAlignment.GetPRMSpectrum().SpecIndex);
				this.AllAlignments[i].Print();
				System.out.println(Utils.IntArrayToString(NewAlignment
						.GetShiftedPRMs()) + "\n");
				System.out.println(NewAlignment.GetPathString() + "\nScore:"
						+ NewAlignment.GetPathScore());
				System.out.println("NumMatches: "
						+ NewAlignment.GetNumMatches());
				System.out.println("NumDeletes: "
						+ NewAlignment.GetNumDeletes());
				System.out.println("RunOfDeletes: "
						+ NewAlignment.GetLargestRunDeletes());
				System.out
						.println("Percent of Peaks Deletes: "
								+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
										.GetPRMSpectrum().PRMs.length));
				System.out.println(" Path: " + NewAlignment.GetPathString());
				Utils.WaitForEnter();
			}
			// NEC CHANGE
			if (Debug)
				System.out.println("Min Match States: "
						+ (4.0 / 11 * Model.NumMatchStates));
			// if(NewAlignment != null && NewAlignment.GetNumMatches() >=
			// Utils.Round(4.0/11*Model.NumMatchStates) &&
			// NewAlignment.GetNumDeletes() < NewAlignment.GetNumMatches() &&
			// NewAlignment.GetLargestRunDeletes() <=
			// Math.min(NewAlignment.GetNumMatches(), Utils.MIN_PEAKS_ALIGNED))
			if (NewAlignment != null
					&& HMMAlignment.IsGoodHMMAlignment(NewAlignment, Model)) {
				// SpectrumNode ParentTemp =
				// this.AllAlignments[i].GetSpectrumNode();
				this.AllAlignments[i].GetSpectrumNode().SetHMMAnnotation(
						NewAlignment);
				// InitialScores[i] = NewAlignment.GetPathScore();
				((SpectrumHMMRight) Model).AddAlignmentToModel(NewAlignment);
				InitialScores[Model.AlignedSpectra.size() - 1] = NewAlignment
						.GetPathScore();
				if (Debug)
					System.out.println("InitialScores["
							+ (Model.AlignedSpectra.size() - 1) + "] = "
							+ NewAlignment.GetPathScore());
				if (TestDebug) {
					System.out.println("Aligned "
							+ NewAlignment.GetPRMSpectrum().FileName + " "
							+ NewAlignment.GetPRMSpectrum().SpecIndex);
					System.out
							.println(" Score: " + NewAlignment.GetPathScore());
					System.out.println(" NumMatches: "
							+ NewAlignment.GetNumMatches());
					System.out.println(" NumDeletes: "
							+ NewAlignment.GetNumDeletes());
					System.out.println(" RunOfDeletes: "
							+ NewAlignment.GetLargestRunDeletes());
					System.out
							.println(" Percent of Peaks Deletes: "
									+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
											.GetPRMSpectrum().PRMs.length));
					System.out
							.println(" Path: " + NewAlignment.GetPathString());
				}

			} else {
				// If this spectrum doesn't align well, then reset it so it can
				// be used again later
				SpectrumNode ParentTemp = this.AllAlignments[i]
						.GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				this.AllAlignments[i].SetSpectrumNode(null);
				ParentTemp.SetHMMAnnotation(null);

				if (TestDebug && NewAlignment != null) {
					System.out.println("Aligned "
							+ NewAlignment.GetPRMSpectrum().FileName + " "
							+ NewAlignment.GetPRMSpectrum().SpecIndex);
					System.out
							.println("**This alignment is bad, throwing it out!!");
					System.out
							.println(" Score: " + NewAlignment.GetPathScore());
					System.out.println(" NumMatches: "
							+ NewAlignment.GetNumMatches());
					System.out.println(" NumDeletes: "
							+ NewAlignment.GetNumDeletes());
					System.out.println(" RunOfDeletes: "
							+ NewAlignment.GetLargestRunDeletes());
					System.out
							.println(" Percent of Peaks Deletes: "
									+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
											.GetPRMSpectrum().PRMs.length));
					System.out
							.println(" Path: " + NewAlignment.GetPathString());
					// System.out.println(" Percent of Peaks Matched: " +
					// (((double)NewAlignment.GetNumMatches())/NewAlignment.GetPRMSpectrum().PRMs.length));

					Utils.WaitForEnter();
				}
				continue;
			}
			if (Debug || TestDebug) {
				Model.DebugPrint();
				Utils.WaitForEnter();
			}

			// ArrayList[] InsertCandidates =
			// Model.GetMatchCandidates(MeanMatchScore*2/3);
			// ArrayList[] InsertCandidates =
			// Model.GetMatchCandidatesOddsScore(MeanMatchScore*Utils.HMM_MIN_SCORE_FRAC);
			// System.out.println("Calling GetMatchCanddiates SumScore for " +
			// SeqSeed.AnchoredSeedSequence + " after adding " +
			// Model.AlignedSpectra.size() + " spectra");
			ArrayList[] InsertCandidates = Model.GetMatchCandidatesSumScore(
					AntibodyUtils.MIN_PEAK_GROUP_SCORE, Tolerance);
			// System.out.println("GATHERCONSENSUS: New state must have score >= "
			// + AntibodyUtils.MIN_PEAK_GROUP_SCORE);

			boolean AddedState = false;
			ArrayList Spectra = Model.GetAlignedSpectra();
			if (Debug || TestDebug)
				System.out.println("Num spectra in alignment: "
						+ Spectra.size());
			for (int k = 0; k < InsertCandidates.length; ++k) {
				if (InsertCandidates[k] != null) {
					// Debug = true;
					ArrayList CurrMasses = InsertCandidates[k];

					int[] Masses = new int[CurrMasses.size()];
					int Index = 0;
					double Score = 0;
					for (int j = 0; j < CurrMasses.size(); ++j) {
						if (CurrMasses.get(j) != null) {
							Object[] CurrPeak = (Object[]) (CurrMasses.get(j));
							Masses[Index] = ((Integer) (CurrPeak[1]))
									.intValue();
							Score += ((Double) (CurrPeak[2])).doubleValue();
							Index += 1;
						}
					}
					if (Score < AntibodyUtils.MIN_PEAK_GROUP_SCORE) {
						if (Debug) {
							System.out
									.println("New score is not good enough!!:"
											+ Score);
							Utils.WaitForEnter();
						}
						continue;
					}

					// System.out.println("Insert Candidate: " + Masses[0] +
					// "!!");
					if (Debug || TestDebug) {
						System.out.println("Adding an insert mass near "
								+ Masses[0] + " with " + CurrMasses.size()
								+ " supporting spectra");
						Utils.WaitForEnter();
					}
					Model.AddMatchState(Masses, Tolerance);
					if (Debug && !Debug) {
						// Model.DebugPrint();
						System.out.println("Added the new state!!");
						Utils.WaitForEnter();
					}
					AddedState = true;
				}
			}
			// If we changed the topology, realign the spectra
			if (AddedState) {
				if (Debug || TestDebug) {
					System.out.println("Holy Shit we changed teh topology!!");

					System.out.println("Num spectra in alignment: "
							+ Spectra.size());
					Utils.WaitForEnter();
				}
				HMMAlignment[] NewAlignments = new HMMAlignment[Spectra.size()];
				for (int j = 0; j < Spectra.size(); ++j) {
					HMMAlignment NextAlignment = (HMMAlignment) (Spectra.get(j));
					if (Debug || TestDebug)
						System.out.println("Realigning "
								+ NextAlignment.GetPRMSpectrum().FileName + ":"
								+ NextAlignment.GetPRMSpectrum().SpecIndex);
					double OldScore = InitialScores[j];
					String OldPath = NextAlignment.GetPathString();
					Model.ReAlignSpectrum(NextAlignment, Tolerance);
					if (Debug || TestDebug) {
						System.out.println(" Old Score: " + InitialScores[j]
								+ " -> NewScore: "
								+ NextAlignment.GetPathScore());
						System.out.println(" NumMatches: "
								+ NextAlignment.GetNumMatches());
						System.out.println(" NumDeletes: "
								+ NextAlignment.GetNumDeletes());
						System.out.println(" RunOfDeletes: "
								+ NextAlignment.GetLargestRunDeletes());
						System.out
								.println(" Percent of Peaks Deletes: "
										+ (((double) NextAlignment
												.GetNumDeletes()) / NextAlignment
												.GetPRMSpectrum().PRMs.length));
						System.out.println(Utils.IntArrayToString(NextAlignment
								.GetShiftedPRMs()) + "\n");
						System.out.println(" OldPath: " + OldPath);
						System.out.println(" Path: "
								+ NextAlignment.GetPathString());
					}
					// NEC_CHANGE!!
					// if(NextAlignment.GetNumMatches() >=
					// Utils.Round(4.0/11*Model.NumMatchStates) &&
					// NextAlignment.GetPathScore() >
					// OldScore*Utils.HMM_SCORE_SCALE &&
					// NextAlignment.GetNumDeletes() <
					// NextAlignment.GetNumMatches() &&
					// NextAlignment.GetLargestRunDeletes() <=
					// Math.min(NextAlignment.GetNumMatches(),
					// Utils.MIN_PEAKS_ALIGNED))
					if (HMMAlignment.IsGoodHMMAlignment(NextAlignment, Model)) {
						NewAlignments[j] = NextAlignment;
						// System.out.println("**New Alignment Score: " +
						// NewAlignments[j].GetPathScore());
					} else // If the realignment is poor, clear this
					// spectrumNode for later use
					{
						SkippedAtRealignment += 1;
						if (TestDebug || Debug) {
							System.out.println("Realigned "
									+ NextAlignment.GetPRMSpectrum().FileName
									+ ":"
									+ NextAlignment.GetPRMSpectrum().SpecIndex);

							System.out
									.println("**This alignment is bad, throwing it out!!");
							System.out.println(" Old Score: " + OldScore
									+ " -> NewScore: "
									+ NextAlignment.GetPathScore());
							System.out.println(" NumMatches: "
									+ NextAlignment.GetNumMatches());
							System.out.println(" NumDeletes: "
									+ NextAlignment.GetNumDeletes());
							System.out.println(" RunOfDeletes: "
									+ NextAlignment.GetLargestRunDeletes());
							System.out
									.println(" Percent of Peaks Deletes: "
											+ (((double) NextAlignment
													.GetNumDeletes()) / NextAlignment
													.GetPRMSpectrum().PRMs.length));
							Utils.WaitForEnter();
						}
						SpectrumNode ParentTemp = NextAlignment
								.GetMSAlignmentType().GetSpectrumNode();
						ParentTemp.SetMSAlignment(null);
						NextAlignment.GetMSAlignmentType().GetSpectrumNode()
								.SetHMMAnnotation(null);
						NewAlignments[j] = null;
					}
				}
				Model.AddAlignmentsToModel(NewAlignments);
				if (Debug || TestDebug) {
					Model.DebugPrintSimple();
					Utils.WaitForEnter();
				}
			}
			if (Debug) {
				int[] Consensus = Model.ProduceConsensus();
				String S = "";

				for (int k = 0; k < Consensus.length; ++k) {
					S += Consensus[k] + " ";

				}
				System.out.println("Consensus: " + S);
				Utils.WaitForEnter();
			}
			// System.out.println("Skipped At Realignment: " +
			// SkippedAtRealignment);

		}

		if (AddAtEndFlag) {
			// There might be more good states to add, but no more spectra let's
			// see
			System.out
					.println("Finished adding spectra, considering more possible states");
			while (true) {
				ArrayList[] InsertCandidates = Model
						.GetMatchCandidatesSumScore(
								AntibodyUtils.MIN_PEAK_GROUP_SCORE, Tolerance);
				// System.out.println("GATHERCONSENSUS: New state must have score >= "
				// + AntibodyUtils.MIN_PEAK_GROUP_SCORE);

				boolean AddedState = false;
				ArrayList Spectra = Model.GetAlignedSpectra();
				if (Debug || TestDebug)
					System.out.println("Num spectra in alignment: "
							+ Spectra.size());
				for (int k = 0; k < InsertCandidates.length; ++k) {
					if (InsertCandidates[k] != null) {
						// Debug = true;
						ArrayList CurrMasses = InsertCandidates[k];

						int[] Masses = new int[CurrMasses.size()];
						int Index = 0;
						double Score = 0;
						for (int j = 0; j < CurrMasses.size(); ++j) {
							if (CurrMasses.get(j) != null) {
								Object[] CurrPeak = (Object[]) (CurrMasses
										.get(j));
								Masses[Index] = ((Integer) (CurrPeak[1]))
										.intValue();
								Score += ((Double) (CurrPeak[2])).doubleValue();
								Index += 1;
							}
						}
						if (Score < AntibodyUtils.MIN_PEAK_GROUP_SCORE) {
							if (Debug) {
								System.out
										.println("New score is not good enough!!:"
												+ Score);
								Utils.WaitForEnter();
							}
							continue;
						}

						// System.out.println("Insert Candidate: " + Masses[0] +
						// "!!");
						if (Debug || TestDebug) {
							System.out.println("Adding an insert mass near "
									+ Masses[0] + " with " + CurrMasses.size()
									+ " supporting spectra");
							Utils.WaitForEnter();
						}
						Model.AddMatchState(Masses, Tolerance);
						if (Debug && !Debug) {
							// Model.DebugPrint();
							System.out.println("Added the new state!!");
							Utils.WaitForEnter();
						}
						AddedState = true;
					}
				}
				// If we changed the topology, realign the spectra
				if (AddedState) {
					if (Debug || TestDebug) {
						System.out
								.println("Holy Shit we changed teh topology!!");

						System.out.println("Num spectra in alignment: "
								+ Spectra.size());
						Utils.WaitForEnter();
					}
					HMMAlignment[] NewAlignments = new HMMAlignment[Spectra
							.size()];
					for (int j = 0; j < Spectra.size(); ++j) {
						HMMAlignment NextAlignment = (HMMAlignment) (Spectra
								.get(j));
						if (Debug || TestDebug)
							System.out.println("Realigning "
									+ NextAlignment.GetPRMSpectrum().FileName
									+ ":"
									+ NextAlignment.GetPRMSpectrum().SpecIndex);
						double OldScore = InitialScores[j];
						String OldPath = NextAlignment.GetPathString();
						Model.ReAlignSpectrum(NextAlignment, Tolerance);
						if (Debug || TestDebug) {
							System.out.println(" Old Score: "
									+ InitialScores[j] + " -> NewScore: "
									+ NextAlignment.GetPathScore());
							System.out.println(" NumMatches: "
									+ NextAlignment.GetNumMatches());
							System.out.println(" NumDeletes: "
									+ NextAlignment.GetNumDeletes());
							System.out.println(" RunOfDeletes: "
									+ NextAlignment.GetLargestRunDeletes());
							System.out
									.println(" Percent of Peaks Deletes: "
											+ (((double) NextAlignment
													.GetNumDeletes()) / NextAlignment
													.GetPRMSpectrum().PRMs.length));
							System.out.println(Utils
									.IntArrayToString(NextAlignment
											.GetShiftedPRMs())
									+ "\n");
							System.out.println(" OldPath: " + OldPath);
							System.out.println(" Path: "
									+ NextAlignment.GetPathString());
						}
						// NEC_CHANGE!!
						// if(NextAlignment.GetNumMatches() >=
						// Utils.Round(4.0/11*Model.NumMatchStates) &&
						// NextAlignment.GetPathScore() >
						// OldScore*Utils.HMM_SCORE_SCALE &&
						// NextAlignment.GetNumDeletes() <
						// NextAlignment.GetNumMatches() &&
						// NextAlignment.GetLargestRunDeletes() <=
						// Math.min(NextAlignment.GetNumMatches(),
						// Utils.MIN_PEAKS_ALIGNED))
						if (HMMAlignment.IsGoodHMMAlignment(NextAlignment,
								Model)) {
							NewAlignments[j] = NextAlignment;
							// System.out.println("**New Alignment Score: " +
							// NewAlignments[j].GetPathScore());
						} else // If the realignment is poor, clear this
						// spectrumNode for later use
						{
							SkippedAtRealignment += 1;
							if (TestDebug || Debug) {
								System.out
										.println("Realigned "
												+ NextAlignment
														.GetPRMSpectrum().FileName
												+ ":"
												+ NextAlignment
														.GetPRMSpectrum().SpecIndex);

								System.out
										.println("**This alignment is bad, throwing it out!!");
								System.out.println(" Old Score: " + OldScore
										+ " -> NewScore: "
										+ NextAlignment.GetPathScore());
								System.out.println(" NumMatches: "
										+ NextAlignment.GetNumMatches());
								System.out.println(" NumDeletes: "
										+ NextAlignment.GetNumDeletes());
								System.out.println(" RunOfDeletes: "
										+ NextAlignment.GetLargestRunDeletes());
								System.out
										.println(" Percent of Peaks Deletes: "
												+ (((double) NextAlignment
														.GetNumDeletes()) / NextAlignment
														.GetPRMSpectrum().PRMs.length));
								Utils.WaitForEnter();
							}
							SpectrumNode ParentTemp = NextAlignment
									.GetMSAlignmentType().GetSpectrumNode();
							ParentTemp.SetMSAlignment(null);
							ParentTemp.SetHMMAnnotation(null);
							NewAlignments[j] = null;
						}
					}
					Model.AddAlignmentsToModel(NewAlignments);
					if (Debug || TestDebug) {
						Model.DebugPrintSimple();
						Utils.WaitForEnter();
					}
				} else {
					if (Debug || TestDebug)
						System.out.println("No more states can be added!!!");

					break;
				}
				if (Debug) {
					int[] Consensus = Model.ProduceConsensus();
					String S = "";

					for (int k = 0; k < Consensus.length; ++k) {
						S += Consensus[k] + " ";

					}
					System.out.println("Consensus: " + S);
					Utils.WaitForEnter();
				}

			}
		}

		// Model.DebugPrintSimple();
		// ConsensusPRMSpectrum Consensus =
		// Model.ProduceConsensusWithScoresWithFilter();

		return Model;

	}

	/**
	 * We align all spectra, then update the model topology
	 * 
	 * @param Alignments
	 * @param OutputFileName
	 * @param SeqSeed
	 * @param Oracle
	 * @param ClusterFlag
	 * @return
	 */
	public ConsensusPRMSpectrum RunGatherAlignWithScoresRightAddAtEnd(
			MSAlignmentType[] Alignments, String OutputFileName, Seed SeqSeed,
			Hashtable Oracle, boolean ClusterFlag, double Tolerance) {
		this.AllAlignments = Alignments;
		Model = new SpectrumHMMRight(SeqSeed, Tolerance, this.applyMatchPenalty, this.matchPenalty);
		// Model = new SpectrumHMM(SeedSequence);

		boolean TestDebug = true;
		double MeanMatchScore = Model.GetMeanMatchScore();

		double[] InitialScores = new double[Alignments.length];
		if (Debug)
			System.out.println("MeanMatchScore: " + MeanMatchScore);
		if (TestDebug) {
			Model.DebugPrintSimple();
			Utils.WaitForEnter();
		}
		/*
		 * For each Alignment in AllAlignments: Attempt to align to the model if
		 * score of alignment > threshold: AddAlignmentToModel If Model should
		 * be changed: 1)Change Model 2) Realign spectra 3) Read good spectra
		 * Create a consensus
		 */

		int SkippedAtRealignment = 0;
		for (int i = 0; i < this.AllAlignments.length; ++i) {
			if (this.AllAlignments[i] == null)
				System.out.println("SHIT!! " + i);
			// If this Spectrum already is used in an HMM then skip it
			if (this.AllAlignments[i].GetSpectrumNode().GetHMMAnnotation() != null)
				continue;

			// Ensure similarity to previous guys
			if (ClusterFlag) {
				boolean Fails = false;

				for (int j = 0; j < Model.AlignedSpectra.size(); ++j) {
					HMMAlignment PrevAlignment = (HMMAlignment) (Model.AlignedSpectra
							.get(j));

					double Score = AntibodyUtils
							.AlignmentScoreLengthWithMassShiftRight(
									this.AllAlignments[i],
									PrevAlignment.GetMSAlignmentType(),
									Tolerance);

					if (Score < AntibodyUtils.MSALIGN_MIN_SIMILARITY) {
						System.out
								.println("Removing this spectrum because it conflicts!!");
						System.out.println("CurrSpectrum: "
								+ this.AllAlignments[i].toStringDebug(Oracle));
						System.out.println("PrevSpectrum: "
								+ PrevAlignment.GetMSAlignmentType()
										.toStringDebug(Oracle));

						Fails = true;
						break;
					}
				}
				if (Fails)
					continue;
			}

			if (Debug || TestDebug) {
				System.out.println("Aligning next spectrum...");
				System.out.println(this.AllAlignments[i].toStringDebug());
			}

			HMMAlignment NewAlignment = Model.AlignSpectrum(
					this.AllAlignments[i], Tolerance);

			if (NewAlignment == null)
				System.out.println("HMM Alignment not possible!!");

			if (Debug && NewAlignment != null) {
				System.out.println("Aligned "
						+ NewAlignment.GetPRMSpectrum().FileName + " "
						+ NewAlignment.GetPRMSpectrum().SpecIndex);
				this.AllAlignments[i].Print();
				System.out.println(Utils.IntArrayToString(NewAlignment
						.GetShiftedPRMs()) + "\n");
				System.out.println(NewAlignment.GetPathString() + "\nScore:"
						+ NewAlignment.GetPathScore());
				System.out.println("NumMatches: "
						+ NewAlignment.GetNumMatches());
				System.out.println("NumDeletes: "
						+ NewAlignment.GetNumDeletes());
				System.out.println("RunOfDeletes: "
						+ NewAlignment.GetLargestRunDeletes());
				System.out
						.println("Percent of Peaks Deletes: "
								+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
										.GetPRMSpectrum().PRMs.length));
				System.out.println(" Path: " + NewAlignment.GetPathString());
				Utils.WaitForEnter();
			}
			// NEC CHANGE
			// if(NewAlignment != null && NewAlignment.GetNumMatches() >=
			// Utils.Round(4.0/11*Model.NumMatchStates) &&
			// NewAlignment.GetNumDeletes() < NewAlignment.GetNumMatches() &&
			// NewAlignment.GetLargestRunDeletes() <=
			// Math.min(NewAlignment.GetNumMatches(), Utils.MIN_PEAKS_ALIGNED))
			if (NewAlignment != null
					&& HMMAlignment.IsGoodHMMAlignment(NewAlignment, Model)) {
				// SpectrumNode ParentTemp =
				// this.AllAlignments[i].GetSpectrumNode();
				this.AllAlignments[i].GetSpectrumNode().SetHMMAnnotation(
						NewAlignment);
				// InitialScores[i] = NewAlignment.GetPathScore();
				((SpectrumHMMRight) Model).AddAlignmentToModel(NewAlignment);
				InitialScores[Model.AlignedSpectra.size() - 1] = NewAlignment
						.GetPathScore();
				if (Debug)
					System.out.println("InitialScores["
							+ (Model.AlignedSpectra.size() - 1) + "] = "
							+ NewAlignment.GetPathScore());
				if (TestDebug) {
					System.out.println("Aligned "
							+ NewAlignment.GetPRMSpectrum().FileName + " "
							+ NewAlignment.GetPRMSpectrum().SpecIndex);
					System.out
							.println(" Score: " + NewAlignment.GetPathScore());
					System.out.println(" NumMatches: "
							+ NewAlignment.GetNumMatches());
					System.out.println(" NumDeletes: "
							+ NewAlignment.GetNumDeletes());
					System.out.println(" RunOfDeletes: "
							+ NewAlignment.GetLargestRunDeletes());
					System.out
							.println(" Percent of Peaks Deletes: "
									+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
											.GetPRMSpectrum().PRMs.length));
					System.out
							.println(" Path: " + NewAlignment.GetPathString());
				}

			} else {
				// If this spectrum doesn't align well, then reset it so it can
				// be used again later
				SpectrumNode ParentTemp = this.AllAlignments[i]
						.GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				this.AllAlignments[i].SetSpectrumNode(null);
				ParentTemp.SetHMMAnnotation(null);

				if (TestDebug && NewAlignment != null) {
					System.out.println("Aligned "
							+ NewAlignment.GetPRMSpectrum().FileName + " "
							+ NewAlignment.GetPRMSpectrum().SpecIndex);
					System.out
							.println("**This alignment is bad, throwing it out!!");
					System.out
							.println(" Score: " + NewAlignment.GetPathScore());
					System.out.println(" NumMatches: "
							+ NewAlignment.GetNumMatches());
					System.out.println(" NumDeletes: "
							+ NewAlignment.GetNumDeletes());
					System.out.println(" RunOfDeletes: "
							+ NewAlignment.GetLargestRunDeletes());
					System.out
							.println(" Percent of Peaks Deletes: "
									+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
											.GetPRMSpectrum().PRMs.length));
					System.out
							.println(" Path: " + NewAlignment.GetPathString());
					// System.out.println(" Percent of Peaks Matched: " +
					// (((double)NewAlignment.GetNumMatches())/NewAlignment.GetPRMSpectrum().PRMs.length));

					Utils.WaitForEnter();
				}
				continue;
			}
			if (Debug || TestDebug) {
				Model.DebugPrintSimple();
				Utils.WaitForEnter();
			}

		}
		// ArrayList[] InsertCandidates =
		// Model.GetMatchCandidates(MeanMatchScore*2/3);
		// ArrayList[] InsertCandidates =
		// Model.GetMatchCandidatesOddsScore(MeanMatchScore*Utils.HMM_MIN_SCORE_FRAC);
		// System.out.println("Calling GetMatchCanddiates SumScore for " +
		// SeqSeed.AnchoredSeedSequence + " after adding " +
		// Model.AlignedSpectra.size() + " spectra");
		if (TestDebug) {
			System.out.println("Finished adding all alignments");
			Utils.WaitForEnter();
		}
		int StatesAdded = 0;
		while (true) {
			ArrayList[] InsertCandidates = Model.GetMatchCandidatesSumScore(
					AntibodyUtils.MIN_PEAK_GROUP_SCORE, Tolerance);
			if (TestDebug)
				System.out.println("States added so far: " + StatesAdded);

			boolean AddedState = false;
			ArrayList Spectra = Model.GetAlignedSpectra();
			if (Debug || TestDebug)
				System.out.println("Num spectra in alignment: "
						+ Spectra.size());
			for (int k = 0; k < InsertCandidates.length; ++k) {
				if (InsertCandidates[k] != null) {
					// Debug = true;
					ArrayList CurrMasses = InsertCandidates[k];
					// int[] Masses = new int[CurrMasses.size()];
					// Check that masses do not come from same spectra

					int[] Masses = new int[CurrMasses.size()];
					int Index = 0;
					double Score = 0;
					for (int j = 0; j < CurrMasses.size(); ++j) {
						if (CurrMasses.get(j) != null) {
							Object[] CurrPeak = (Object[]) CurrMasses.get(j);
							Masses[Index] = ((Integer) (CurrPeak[1]))
									.intValue();
							Score += ((Double) (CurrPeak[2])).doubleValue();
							Index += 1;
						}
					}

					// System.out.println("Insert Candidate: " + Masses[0] +
					// "!!");
					if (Debug || TestDebug) {
						System.out.println("Adding an insert mass near "
								+ Masses[0] + " with " + CurrMasses.size()
								+ " supporting spectra");
						Utils.WaitForEnter();
					}
					Model.AddMatchState(Masses, Tolerance);
					if (Debug && !Debug) {
						// Model.DebugPrint();
						System.out.println("Added the new state!!");
						Utils.WaitForEnter();
					}
					AddedState = true;
				}
			}
			// If we changed the topology, realign the spectra
			if (AddedState) {
				if (Debug || TestDebug) {
					System.out.println("Holy Shit we changed teh topology!!");

					System.out.println("Num spectra in alignment: "
							+ Spectra.size());
					Utils.WaitForEnter();
				}
				HMMAlignment[] NewAlignments = new HMMAlignment[Spectra.size()];
				for (int j = 0; j < Spectra.size(); ++j) {
					HMMAlignment NextAlignment = (HMMAlignment) (Spectra.get(j));
					if (Debug || TestDebug)
						System.out.println("Realigning "
								+ NextAlignment.GetPRMSpectrum().FileName + ":"
								+ NextAlignment.GetPRMSpectrum().SpecIndex);
					double OldScore = InitialScores[j];
					String OldPath = NextAlignment.GetPathString();
					Model.ReAlignSpectrum(NextAlignment, Tolerance);
					if (Debug || TestDebug) {
						System.out.println(" Old Score: " + InitialScores[j]
								+ " -> NewScore: "
								+ NextAlignment.GetPathScore());
						System.out.println(" NumMatches: "
								+ NextAlignment.GetNumMatches());
						System.out.println(" NumDeletes: "
								+ NextAlignment.GetNumDeletes());
						System.out.println(" RunOfDeletes: "
								+ NextAlignment.GetLargestRunDeletes());
						System.out
								.println(" Percent of Peaks Deletes: "
										+ (((double) NextAlignment
												.GetNumDeletes()) / NextAlignment
												.GetPRMSpectrum().PRMs.length));
						System.out.println(Utils.IntArrayToString(NextAlignment
								.GetShiftedPRMs()) + "\n");
						System.out.println(" OldPath: " + OldPath);
						System.out.println(" Path: "
								+ NextAlignment.GetPathString());
					}
					// NEC_CHANGE!!
					// if(NextAlignment.GetNumMatches() >=
					// Utils.Round(4.0/11*Model.NumMatchStates) &&
					// NextAlignment.GetPathScore() >
					// OldScore*Utils.HMM_SCORE_SCALE &&
					// NextAlignment.GetNumDeletes() <
					// NextAlignment.GetNumMatches() &&
					// NextAlignment.GetLargestRunDeletes() <=
					// Math.min(NextAlignment.GetNumMatches(),
					// Utils.MIN_PEAKS_ALIGNED))
					if (HMMAlignment.IsGoodHMMAlignment(NextAlignment, Model)) {
						NewAlignments[j] = NextAlignment;
						// System.out.println("**New Alignment Score: " +
						// NewAlignments[j].GetPathScore());
					} else // If the realignment is poor, clear this
					// spectrumNode for later use
					{
						SkippedAtRealignment += 1;
						if (TestDebug || Debug) {
							System.out.println("Realigned "
									+ NextAlignment.GetPRMSpectrum().FileName
									+ ":"
									+ NextAlignment.GetPRMSpectrum().SpecIndex);

							System.out
									.println("**This alignment is bad, throwing it out!!");
							System.out.println(" Old Score: " + OldScore
									+ " -> NewScore: "
									+ NextAlignment.GetPathScore());
							System.out.println(" NumMatches: "
									+ NextAlignment.GetNumMatches());
							System.out.println(" NumDeletes: "
									+ NextAlignment.GetNumDeletes());
							System.out.println(" RunOfDeletes: "
									+ NextAlignment.GetLargestRunDeletes());
							System.out
									.println(" Percent of Peaks Deletes: "
											+ (((double) NextAlignment
													.GetNumDeletes()) / NextAlignment
													.GetPRMSpectrum().PRMs.length));
							Utils.WaitForEnter();
						}
						SpectrumNode ParentTemp = NextAlignment
								.GetMSAlignmentType().GetSpectrumNode();
						ParentTemp.SetMSAlignment(null);
						ParentTemp.SetHMMAnnotation(null);
						NewAlignments[j] = null;
					}
				}
				Model.AddAlignmentsToModel(NewAlignments);
				if (Debug || TestDebug) {
					Model.DebugPrintSimple();
					Utils.WaitForEnter();
				}
			} else
				// we added no new states, we must have no states to add!!
				break;
			if (Debug) {
				int[] Consensus = Model.ProduceConsensus();
				String S = "";

				for (int k = 0; k < Consensus.length; ++k) {
					S += Consensus[k] + " ";

				}
				System.out.println("Consensus: " + S);
				Utils.WaitForEnter();
			}
			// System.out.println("Skipped At Realignment: " +
			// SkippedAtRealignment);
		}
		// Model.DebugPrintSimple();
		// ConsensusPRMSpectrum Consensus =
		// Model.ProduceConsensusWithScoresWithFilter();
		ConsensusPRMSpectrum Consensus = Model.ProduceConsensusWithScores();
		// String S = "";
		// String P = "";
		// for(int i = 0; i < Consensus.ScaledPeaks.length; ++i)
		// {
		// S += Consensus.ScaledPeaks[i] + " ";
		// P += Consensus.PeakScores[i] + " ";
		// }
		// System.out.println("DONE GATHER CONSENSUS:");
		// System.out.println("  " +S);
		// System.out.println("  " + P);
		// ConsensusPRMSpectrum Consensus =
		// Model.ProduceConsensusWithScoresSimple();

		// Determine mean and standard deviation of match and insert states
		/*
		 * for(int i = 0; i < Model.NumMatchStates; ++i) { ArrayList Peaks =
		 * Model.MatchedPeaks[i]; double [] PeakMasses = new
		 * double[Peaks.size()]; double Mean = 0; for(int j = 0; j <
		 * PeakMasses.length; ++j) { double[] CurrPeak =
		 * (double[])(Peaks.get(j)); PeakMasses[j] = CurrPeak[0]; Mean +=
		 * CurrPeak[0]; } Mean = Mean/PeakMasses.length; double StdDev = 0.0;
		 * for(int j = 0; j < PeakMasses.length; ++j) { StdDev +=
		 * Math.pow(PeakMasses[j]-Mean,2); } StdDev = StdDev/PeakMasses.length;
		 * StdDev = Math.sqrt(StdDev); System.out.println("Match State " + i +
		 * ", " + PeakMasses.length + " peaks, mean=" + Mean + ",stddev=" +
		 * StdDev);
		 * 
		 * Peaks = Model.InsertMasses[i]; PeakMasses = new double[Peaks.size()];
		 * Mean = 0; for(int j = 0; j < PeakMasses.length; ++j) { double[]
		 * CurrPeak = (double[])(Peaks.get(j)); PeakMasses[j] = CurrPeak[0];
		 * Mean += CurrPeak[0]; } Mean = Mean/PeakMasses.length; StdDev = 0.0;
		 * for(int j = 0; j < PeakMasses.length; ++j) { StdDev +=
		 * Math.pow(PeakMasses[j]-Mean,2); } StdDev = StdDev/PeakMasses.length;
		 * StdDev = Math.sqrt(StdDev); System.out.println("Insert State " + i +
		 * ", " + PeakMasses.length + " peaks, mean=" + Mean + ",stddev=" +
		 * StdDev);
		 * 
		 * } Utils.WaitForEnter();
		 */

		return Consensus;
		// Model.ProduceConsensus(ScaledPeaks,PeakScores);
		// return Model.ProduceConsensus();

	}

	public ConsensusPRMSpectrum RunGatherAlignWithScoresLeft(
			MSAlignmentType[] Alignments, String OutputFileName, Seed SeqSeed,
			Hashtable Oracle, boolean ClusterFlag, boolean AddAtEndFlag,
			double Tolerance) {
		// this.AllAlignments =
		// MSAlignmentType.ClusterRankAlignmentsLeft(Alignments);
		// if(Debug || !Debug)
		// System.out.println("\nClustering reduced spectra to " +
		// this.AllAlignments.length);
		this.AllAlignments = Alignments;
		// String SeedSequence = SeqSeed.GetSeedSequence();
		// Model = new SpectrumHMM(SeedSequence);
		Model = new SpectrumHMMLeft(SeqSeed, Tolerance, this.applyMatchPenalty, this.matchPenalty);

		if (Debug) {
			Model.DebugPrintSimple();
			Utils.WaitForEnter();
		}
		double MeanMatchScore = Model.GetMeanMatchScore();
		/*
		 * For each Alignment in AllAlignments: Attempt to align to the model if
		 * score of alignment > threshold: AddAlignmentToModel If Model should
		 * be changed: 1)Change Model 2) Realign spectra 3) Readd good spectra
		 * Create a consensus
		 */

		int SkippedAtRealignment = 0;
		int Rank = 1;
		for (int i = 0; i < this.AllAlignments.length; ++i) {
			// If this Spectrum already is used in an HMM then skip it
			// if (this.AllAlignments[i].GetSpectrumNode().GetHMMAnnotation() !=
			// null) {
			// System.out.println("We used to skip this one, but now its ok!");
			// this.AllAlignments[i].Print();
			// Utils.WaitForEnter();
			// }

			// Ensure similarity to previous guys
			if (ClusterFlag) {
				boolean Fails = false;

				for (int j = 0; j < Model.AlignedSpectra.size(); ++j) {
					HMMAlignment PrevAlignment = (HMMAlignment) (Model.AlignedSpectra
							.get(j));

					double Score = AntibodyUtils
							.AlignmentScoreLengthWithMassShiftRight(
									this.AllAlignments[i],
									PrevAlignment.GetMSAlignmentType(),
									Tolerance);

					if (Score < AntibodyUtils.MSALIGN_MIN_SIMILARITY) {
						System.out
								.println("Removing this spectrum because it conflicts!!");
						System.out.println("CurrSpectrum: "
								+ this.AllAlignments[i].toStringDebug(Oracle));
						System.out.println("PrevSpectrum: "
								+ PrevAlignment.GetMSAlignmentType()
										.toStringDebug(Oracle));
						// Utils.WaitForEnter();
						Fails = true;
						break;
					}

				}
				if (Fails)
					continue;
			}

			HMMAlignment NewAlignment = Model.AlignSpectrum(
					this.AllAlignments[i], Tolerance);

			if (Debug && NewAlignment != null) {
				System.out.println("Aligned "
						+ NewAlignment.GetPRMSpectrum().FileName + " "
						+ NewAlignment.GetPRMSpectrum().SpecIndex);
				this.AllAlignments[i].Print();
				System.out.println(NewAlignment.GetPathString() + "\nScore:"
						+ NewAlignment.GetPathScore());
				Utils.WaitForEnter();
			}
			if (NewAlignment != null
					&& HMMAlignment.IsGoodHMMAlignment(NewAlignment, Model)) {
				// SpectrumNode ParentTemp =
				// this.AllAlignments[i].GetSpectrumNode();
				this.AllAlignments[i].GetSpectrumNode().SetHMMAnnotation(
						NewAlignment);
				Model.AddAlignmentToModel(NewAlignment);
				Rank++;
				if (Debug) {
					System.out.println(" NumMatches: "
							+ NewAlignment.GetNumMatches());
					System.out.println(" NumDeletes: "
							+ NewAlignment.GetNumDeletes());
					System.out.println(" RunOfDeletes: "
							+ NewAlignment.GetLargestRunDeletes());
					System.out
							.println(" Percent of Peaks Deletes: "
									+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
											.GetPRMSpectrum().PRMs.length));
					System.out.println(Utils.IntArrayToString(NewAlignment
							.GetShiftedPRMs()) + "\n");
					// System.out.println(" OldPath: " + OldPath);
					System.out
							.println(" Path: " + NewAlignment.GetPathString());
				}

			} else {
				// If this spectrum doesn't align well, then reset it so it can
				// be used again later
				SpectrumNode ParentTemp = this.AllAlignments[i]
						.GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				this.AllAlignments[i].SetSpectrumNode(null);
				ParentTemp.SetHMMAnnotation(null);
				if (Debug) {
					System.out.println("EXCLUDED:This alignment has less than "
							+ AntibodyUtils.MIN_PEAKS_ALIGNED + " matches");
					System.out.println(" NumMatches: "
							+ NewAlignment.GetNumMatches());
					System.out.println(" NumDeletes: "
							+ NewAlignment.GetNumDeletes());
					System.out.println(" RunOfDeletes: "
							+ NewAlignment.GetLargestRunDeletes());
					System.out
							.println(" Percent of Peaks Deletes: "
									+ (((double) NewAlignment.GetNumDeletes()) / NewAlignment
											.GetPRMSpectrum().PRMs.length));
					System.out.println(Utils.IntArrayToString(NewAlignment
							.GetShiftedPRMs()) + "\n");
					// System.out.println(" OldPath: " + OldPath);
					System.out
							.println(" Path: " + NewAlignment.GetPathString());
					System.out.println(" ");
					Utils.WaitForEnter();
				}
				continue;
			}
			if (Debug) {
				Model.DebugPrintSimple();
				Utils.WaitForEnter();
			}

			// ArrayList[] InsertCandidates =
			// Model.GetMatchCandidates(MeanMatchScore*2/3);
			// System.out.println("GATHERCONSENSUS: PeakGroup Threshold Score: "
			// + AntibodyUtils.MIN_PEAK_GROUP_SCORE);
			ArrayList[] InsertCandidates = Model.GetMatchCandidatesSumScore(
					AntibodyUtils.MIN_PEAK_GROUP_SCORE, Tolerance);
			// System.out.println("GATHERCONSENSUS: New state must have score >= "
			// + (MeanMatchScore*Utils.HMM_MIN_SCORE_FRAC));

			boolean AddedState = false;
			ArrayList Spectra = Model.GetAlignedSpectra();
			if (Debug)
				System.out.println("Num spectra in alignment: "
						+ Spectra.size());
			for (int k = 0; k < InsertCandidates.length; ++k) {
				if (InsertCandidates[k] != null) {
					// Debug = true;
					ArrayList CurrMasses = InsertCandidates[k];

					int[] Masses = new int[CurrMasses.size()];
					int Index = 0;
					double Score = 0;
					for (int j = 0; j < CurrMasses.size(); ++j) {
						if (CurrMasses.get(j) != null) {
							Object[] CurrPeak = (Object[]) (CurrMasses.get(j));
							Masses[Index] = ((Integer) (CurrPeak[1]))
									.intValue();
							Score += ((Double) (CurrPeak[2])).doubleValue();
							if (Debug)
								System.out.println("Adding mass "
										+ Masses[Index]);
							Index += 1;
						} else if (Debug)
							System.out.println("Original Mass at index " + j
									+ " was removed!!");
					}
					if (Score < AntibodyUtils.MIN_PEAK_GROUP_SCORE) {
						if (Debug) {
							System.out
									.println("New score is not good enough!!:"
											+ Score);
							Utils.WaitForEnter();
						}
						continue;
					}
					// System.out.println("Insert Candidate: " + Masses[0] +
					// "!!");
					if (Debug) {
						System.out.println("Adding an insert mass near "
								+ Masses[0] + " with " + CurrMasses.size()
								+ " supporting spectra");
						Utils.WaitForEnter();
					}
					// Debug = false;
					Model.AddMatchState(Masses, Tolerance);
					if (Debug && !Debug) {
						// Model.DebugPrint();
						System.out.println("Added the new state!!");
						Utils.WaitForEnter();
					}
					AddedState = true;
				}
			}
			// If we changed the topology, realign the spectra
			if (AddedState) {
				if (Debug) {
					System.out.println("Holy Shit we changed teh topology!!");

					System.out.println("Num spectra in alignment: "
							+ Spectra.size());
					Utils.WaitForEnter();
				}
				HMMAlignment[] NewAlignments = new HMMAlignment[Spectra.size()];
				for (int j = 0; j < Spectra.size(); ++j) {
					HMMAlignment NextAlignment = (HMMAlignment) (Spectra.get(j));
					if (Debug)
						System.out.println("Realigning "
								+ NextAlignment.GetPRMSpectrum().FileName + ":"
								+ NextAlignment.GetPRMSpectrum().SpecIndex);
					double OldScore = NextAlignment.GetPathScore();
					String OldPath = NextAlignment.GetPathString();
					Model.ReAlignSpectrum(NextAlignment, Tolerance);
					if (HMMAlignment.IsGoodHMMAlignment(NextAlignment, Model)) {
						NewAlignments[j] = NextAlignment;
						if (Debug) {
							System.out.println(" NumMatches: "
									+ NextAlignment.GetNumMatches());
							System.out.println(" NumDeletes: "
									+ NextAlignment.GetNumDeletes());
							System.out.println(" RunOfDeletes: "
									+ NextAlignment.GetLargestRunDeletes());
							System.out
									.println(" Percent of Peaks Deletes: "
											+ (((double) NextAlignment
													.GetNumDeletes()) / NextAlignment
													.GetPRMSpectrum().PRMs.length));
							System.out.println(Utils
									.IntArrayToString(NextAlignment
											.GetShiftedPRMs())
									+ "\n");
							System.out.println(" OldPath: " + OldPath);
							System.out.println(" Path: "
									+ NextAlignment.GetPathString());
							// System.out.println("**New Alignment Score: " +
							// NewAlignments[j].GetPathScore());
						}
					} else // If the realignment is poor, clear this
					// spectrumNode for later use
					{
						SkippedAtRealignment += 1;
						if (Debug) {
							System.out
									.println("**This alignment is bad, throwing it out!!");
							System.out.println(" NumMatches: "
									+ NextAlignment.GetNumMatches());
							System.out.println(" NumDeletes: "
									+ NextAlignment.GetNumDeletes());
							System.out.println(" RunOfDeletes: "
									+ NextAlignment.GetLargestRunDeletes());
							System.out
									.println(" Percent of Peaks Deletes: "
											+ (((double) NextAlignment
													.GetNumDeletes()) / NextAlignment
													.GetPRMSpectrum().PRMs.length));
							System.out.println(Utils
									.IntArrayToString(NextAlignment
											.GetShiftedPRMs())
									+ "\n");
							System.out.println(" OldPath: " + OldPath);
							System.out.println(" Path: "
									+ NextAlignment.GetPathString());
							Utils.WaitForEnter();
						}
						SpectrumNode ParentTemp = NextAlignment
								.GetMSAlignmentType().GetSpectrumNode();
						ParentTemp.SetMSAlignment(null);
						ParentTemp.SetHMMAnnotation(null);
						NewAlignments[j] = null;
					}
				}
				Model.AddAlignmentsToModel(NewAlignments);
				if (Debug) {
					Model.DebugPrintSimple();
					Utils.WaitForEnter();
				}
			}
			if (Debug) {
				int[] Consensus = Model.ProduceConsensus();
				String S = "";

				for (int k = 0; k < Consensus.length; ++k) {
					S += Consensus[k] + " ";

				}
				System.out.println("Consensus: " + S);
				Utils.WaitForEnter();
			}
			// System.out.println("Skipped At Realignment: " +
			// SkippedAtRealignment);

		}
		if (AddAtEndFlag) {
			System.out
					.println("Finished adding spectra, considering further match states");
			while (true) {
				ArrayList[] InsertCandidates = Model
						.GetMatchCandidatesSumScore(
								AntibodyUtils.MIN_PEAK_GROUP_SCORE, Tolerance);
				// System.out.println("GATHERCONSENSUS: New state must have score >= "
				// + (MeanMatchScore*Utils.HMM_MIN_SCORE_FRAC));

				boolean AddedState = false;
				ArrayList Spectra = Model.GetAlignedSpectra();
				if (Debug)
					System.out.println("Num spectra in alignment: "
							+ Spectra.size());
				for (int k = 0; k < InsertCandidates.length; ++k) {
					if (InsertCandidates[k] != null) {
						// Debug = true;
						ArrayList CurrMasses = InsertCandidates[k];

						int[] Masses = new int[CurrMasses.size()];
						int Index = 0;
						double Score = 0;
						for (int j = 0; j < CurrMasses.size(); ++j) {
							if (CurrMasses.get(j) != null) {
								Object[] CurrPeak = (Object[]) (CurrMasses
										.get(j));
								Masses[Index] = ((Integer) (CurrPeak[1]))
										.intValue();
								Score += ((Double) (CurrPeak[2])).doubleValue();
								if (Debug)
									System.out.println("Adding mass "
											+ Masses[Index]);
								Index += 1;
							} else if (Debug)
								System.out.println("Original Mass at index "
										+ j + " was removed!!");
						}
						if (Score < AntibodyUtils.MIN_PEAK_GROUP_SCORE) {
							if (Debug) {
								System.out
										.println("New score is not good enough!!:"
												+ Score);
								Utils.WaitForEnter();
							}
							continue;
						}
						// System.out.println("Insert Candidate: " + Masses[0] +
						// "!!");
						if (Debug) {
							System.out.println("Adding an insert mass near "
									+ Masses[0] + " with " + CurrMasses.size()
									+ " supporting spectra");
							Utils.WaitForEnter();
						}
						// Debug = false;
						Model.AddMatchState(Masses, Tolerance);
						if (Debug && !Debug) {
							// Model.DebugPrint();
							System.out.println("Added the new state!!");
							Utils.WaitForEnter();
						}
						AddedState = true;
					}
				}
				// If we changed the topology, realign the spectra
				if (AddedState) {
					if (Debug) {
						System.out
								.println("Holy Shit we changed teh topology!!");

						System.out.println("Num spectra in alignment: "
								+ Spectra.size());
						Utils.WaitForEnter();
					}
					HMMAlignment[] NewAlignments = new HMMAlignment[Spectra
							.size()];
					for (int j = 0; j < Spectra.size(); ++j) {
						HMMAlignment NextAlignment = (HMMAlignment) (Spectra
								.get(j));
						if (Debug)
							System.out.println("Realigning "
									+ NextAlignment.GetPRMSpectrum().FileName
									+ ":"
									+ NextAlignment.GetPRMSpectrum().SpecIndex);
						double OldScore = NextAlignment.GetPathScore();
						String OldPath = NextAlignment.GetPathString();
						Model.ReAlignSpectrum(NextAlignment, Tolerance);
						if (HMMAlignment.IsGoodHMMAlignment(NextAlignment,
								Model)) {
							NewAlignments[j] = NextAlignment;
							if (Debug) {
								System.out.println(" NumMatches: "
										+ NextAlignment.GetNumMatches());
								System.out.println(" NumDeletes: "
										+ NextAlignment.GetNumDeletes());
								System.out.println(" RunOfDeletes: "
										+ NextAlignment.GetLargestRunDeletes());
								System.out
										.println(" Percent of Peaks Deletes: "
												+ (((double) NextAlignment
														.GetNumDeletes()) / NextAlignment
														.GetPRMSpectrum().PRMs.length));
								System.out.println(Utils
										.IntArrayToString(NextAlignment
												.GetShiftedPRMs())
										+ "\n");
								System.out.println(" OldPath: " + OldPath);
								System.out.println(" Path: "
										+ NextAlignment.GetPathString());
								// System.out.println("**New Alignment Score: "
								// + NewAlignments[j].GetPathScore());
							}
						} else // If the realignment is poor, clear this
						// spectrumNode for later use
						{
							SkippedAtRealignment += 1;
							if (Debug) {
								System.out
										.println("**This alignment is bad, throwing it out!!");
								System.out.println(" NumMatches: "
										+ NextAlignment.GetNumMatches());
								System.out.println(" NumDeletes: "
										+ NextAlignment.GetNumDeletes());
								System.out.println(" RunOfDeletes: "
										+ NextAlignment.GetLargestRunDeletes());
								System.out
										.println(" Percent of Peaks Deletes: "
												+ (((double) NextAlignment
														.GetNumDeletes()) / NextAlignment
														.GetPRMSpectrum().PRMs.length));
								System.out.println(Utils
										.IntArrayToString(NextAlignment
												.GetShiftedPRMs())
										+ "\n");
								System.out.println(" OldPath: " + OldPath);
								System.out.println(" Path: "
										+ NextAlignment.GetPathString());
								Utils.WaitForEnter();
							}
							SpectrumNode ParentTemp = NextAlignment
									.GetMSAlignmentType().GetSpectrumNode();
							ParentTemp.SetMSAlignment(null);
							ParentTemp.SetHMMAnnotation(null);
							NewAlignments[j] = null;
						}
					}
					Model.AddAlignmentsToModel(NewAlignments);
					if (Debug) {
						Model.DebugPrintSimple();
						Utils.WaitForEnter();
					}
				} else {
					System.out.println("No more states to add!!");
					break;
				}
				if (Debug) {
					int[] Consensus = Model.ProduceConsensus();
					String S = "";

					for (int k = 0; k < Consensus.length; ++k) {
						S += Consensus[k] + " ";

					}
					System.out.println("Consensus: " + S);
					Utils.WaitForEnter();
				}

			}
		}
		// Model.DebugPrintSimple();
		ConsensusPRMSpectrum Consensus = Model.ProduceConsensusWithScores();
		Consensus.Model = Model;
		// Consensus.IntervalScores = Model
		// .ComputeProbabilityTable(Consensus.ScaledPeaks);
		// String S = "";
		// String P = "";
		// for(int i = 0; i < Consensus.ScaledPeaks.length; ++i)
		// {
		// S += Consensus.ScaledPeaks[i] + " ";
		// P += Consensus.PeakScores[i] + " ";
		// }
		// System.out.println("DONE GATHER CONSENSUS:");
		// System.out.println("  " +S);
		// System.out.println("  " + P);
		// ConsensusPRMSpectrum Consensus =
		// Model.ProduceConsensusWithScoresSimple();

		return Consensus;
		// Model.ProduceConsensus(ScaledPeaks,PeakScores);
		// return Model.ProduceConsensus();

	}

	public String GetModelMatchString() {
		String Ret = "SpectrumHMM Model:";
		for (int i = 0; i < Model.GetNumMatchStates(); ++i) {
			Ret += " * M[" + i + "] - Mass:" + Model.GetMatchMass(i)
					+ " Overlap:" + Model.GetNumSupporting(i) + "/"
					+ Model.GetNumOverlapping(i) + " LogOddsScore:"
					+ Model.ComputeMatchStateLogOddsScore(i) + ", SumScore:"
					+ Model.ComputeMatchStateSumScore(i) + ", InsertsAfter: "
					+ Model.GetNumInserts(i) + "\n";

		}
		return Ret;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
