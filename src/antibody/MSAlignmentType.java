package antibody;

import java.util.ArrayList;
import java.util.Hashtable;

import basicUtils.InspectAnnotation;
import basicUtils.Utils;

public class MSAlignmentType {

	private double Score;
	private double OddsScore;
	// private String Sequence; //The matched portion of the sequence for this
	// alignment
	private String AnchorSequence; // Full sequence searched
	private int[] AnchorPRMs;
	private PRMSpectrum Spectrum;
	private int[] AlignedPeaks; // Masses of peaks from the spectrum that are
	// aligned
	private int[] AlignedSequence; // Masses of theoretical peaks of the
	// sequence that are aligned

	private boolean Debug = false;

	private SpectrumNode ContainingNode;
	private Seed AnchorSeed;

	// private HMMAlignment HMMAnnotation;

	// private DistributionType NullDistribution;

	/**
	 * Constructor which creates an empty alignment, to signify that no good
	 * alignment could be found between this spectrum and this anchor sequence.
	 * The score is Double.MIN_VALUE and the Sequence is empty.
	 * 
	 * @param Spectrum
	 *            the spectrum for this alignment
	 * @param AnchorSequence
	 *            the sequence for this alignment
	 */
	public MSAlignmentType(PRMSpectrum Spectrum, String AnchorSequence) {
		this.AnchorSequence = AnchorSequence;
		this.Spectrum = Spectrum;

		// this.Sequence = "";
		this.Score = Double.MIN_VALUE;
	}

	/**
	 * Constructor which creates an object representing a valid alignment
	 * 
	 * @param Spectrum
	 *            the aligned spectrum object
	 * @param AnchorSequence
	 *            the sequence to which the spectrum is aligned
	 * @param Score
	 *            the score of the best partial alignment
	 * @param Sequence
	 *            the partial sequence that is aligned
	 * @param AlignedPeaks
	 *            the specific peaks of the spectrum which are aligned
	 * @param AlignedSequence
	 *            the specific theoretical peaks of the sequence that are
	 *            aligned
	 */
	public MSAlignmentType(PRMSpectrum Spectrum, String AnchorSequence,
			double Score, String Sequence, int[] AlignedPeaks,
			int[] AlignedSequence, SpectrumNode Container, Seed CurrSeed) {
		this.ContainingNode = Container;
		this.Spectrum = Spectrum;
		this.AnchorSequence = AnchorSequence;
		this.AnchorPRMs = AntibodyUtils.GetSeqPRMScaled(this.AnchorSequence);
		// this.Sequence = Sequence;
		this.Score = Score;
		this.AlignedPeaks = AlignedPeaks;
		this.AlignedSequence = AlignedSequence;
		this.AnchorSeed = CurrSeed;
		// DistributionType NullDistribution = this.ContainingNode.GetNull();
		// **Permute Subsets of PRMScores, generates many unfeasible subsets
		// DistributionType TempDistribution = new
		// DistributionType(this.AlignedPeaks.length, this.Spectrum.PRMScores);

		// **Attempt to approximate permutation with average mean and variance
		// this.NullDistribution = new
		// DistributionType(this.AlignedPeaks.length*Spectrum.MeanScore,this.AlignedPeaks.length*Spectrum.ScoreVariance);

		// **Align spectrum against decoy database, generates too few decoy
		// alignments
		// if(NullDistribution == null)
		// {

		// String ShuffledProtein =
		// DecoyDB.GetShuffledDB(AnchorSequence.length()*Utils.DECOY_EXPLOSION);
		// double[] NullScores = SpecAlign.AlignSimple(Spectrum,
		// ShuffledProtein, Utils.NUM_DECOY_ALIGNMENTS, AlignedPeaks.length);
		// System.out.println("Found " + NullScores.length +
		// " alignments to the decoy");
		// NullDistribution = new DistributionType(NullScores);

		// **Attempt to generate all equal length paths in the DeNovo graph
		// NullDistribution = new DistributionType(this.AlignedPeaks.length,
		// this.Spectrum.PRMs, this.Spectrum.PRMScores);
		// this.ContainingNode.SetNull(NullDistribution);

		// }
		// else
		// System.out.println("**Reusing the Null Distribution already computed!!!");
		// **Attempt to generate all equal length paths in the DeNovo graph
		// this.NullDistribution = new
		// DistributionType(this.AlignedPeaks.length, this.Spectrum.PRMs,
		// this.Spectrum.PRMScores);

		// this.OddsScore = this.ContainingNode.GetNull().PValue(this.Score);
		if (this.ContainingNode != null)
			this.ContainingNode.SetMSAlignment(this);
		// this.Print();
		// Utils.WaitForEnter();
	}

	public MSAlignmentType(PRMSpectrum Spectrum, String AnchorSequence,
			int[] AnchorPRMs, double Score, int[] AlignedPeaks,
			int[] AlignedSequence, SpectrumNode Container, Seed CurrSeed) {
		this.ContainingNode = Container;
		this.Spectrum = Spectrum;
		// Reset the Spectrum's stuff
		// this.ContainingNode.SetMSAlignment(this);
		// this.ContainingNode.SetHMMAnnotation(null);

		this.AnchorSequence = AnchorSequence;
		this.AnchorPRMs = AnchorPRMs;
		// this.Sequence = Sequence;
		this.Score = Score;
		this.AlignedPeaks = AlignedPeaks;
		this.AlignedSequence = AlignedSequence;
		this.AnchorSeed = CurrSeed;
		// DistributionType NullDistribution = this.ContainingNode.GetNull();
		// **Permute Subsets of PRMScores, generates many unfeasible subsets
		// DistributionType TempDistribution = new
		// DistributionType(this.AlignedPeaks.length, this.Spectrum.PRMScores);

		// **Attempt to approximate permutation with average mean and variance
		// this.NullDistribution = new
		// DistributionType(this.AlignedPeaks.length*Spectrum.MeanScore,this.AlignedPeaks.length*Spectrum.ScoreVariance);

		// **Align spectrum against decoy database, generates too few decoy
		// alignments
		// if(NullDistribution == null)
		// {

		// String ShuffledProtein =
		// DecoyDB.GetShuffledDB(AnchorSequence.length()*Utils.DECOY_EXPLOSION);
		// double[] NullScores = SpecAlign.AlignSimple(Spectrum,
		// ShuffledProtein, Utils.NUM_DECOY_ALIGNMENTS, AlignedPeaks.length);
		// System.out.println("Found " + NullScores.length +
		// " alignments to the decoy");
		// NullDistribution = new DistributionType(NullScores);

		// **Attempt to generate all equal length paths in the DeNovo graph
		// NullDistribution = new DistributionType(this.AlignedPeaks.length,
		// this.Spectrum.PRMs, this.Spectrum.PRMScores);
		// this.ContainingNode.SetNull(NullDistribution);

		// }
		// else
		// System.out.println("**Reusing the Null Distribution already computed!!!");
		// **Attempt to generate all equal length paths in the DeNovo graph
		// this.NullDistribution = new
		// DistributionType(this.AlignedPeaks.length, this.Spectrum.PRMs,
		// this.Spectrum.PRMScores);

		// this.OddsScore = this.ContainingNode.GetNull().PValue(this.Score);
		if (this.ContainingNode != null)
			this.ContainingNode.SetMSAlignment(this);
		// this.Print();
		// Utils.WaitForEnter();
	}

	/**
	 * Constructor which creates an object representing a valid alignment
	 * 
	 * @param Spectrum
	 *            the aligned spectrum object
	 * @param AnchorSequence
	 *            the sequence to which the spectrum is aligned
	 * @param Score
	 *            the score of the best partial alignment
	 * @param Sequence
	 *            the partial sequence that is aligned
	 * @param AlignedPeaks
	 *            the specific peaks of the spectrum which are aligned
	 * @param AlignedSequence
	 *            the specific theoretical peaks of the sequence that are
	 *            aligned
	 */
	/*
	 * public MSAlignmentType(PRMSpectrum Spectrum, String AnchorSequence,double
	 * Score,String Sequence, int[] AlignedPeaks,int[] AlignedSequence) {
	 * this.Spectrum = Spectrum; this.AnchorSequence = AnchorSequence;
	 * this.Sequence = Sequence; this.Score = Score; this.AlignedPeaks =
	 * AlignedPeaks; this.AlignedSequence = AlignedSequence;
	 * 
	 * 
	 * 
	 * this.NullDistribution = new DistributionType(this.AlignedPeaks.length,
	 * this.Spectrum.PRMScores); //this.NullDistribution = new
	 * DistributionType(this
	 * .AlignedPeaks.length*Spectrum.MeanScore,this.AlignedPeaks
	 * .length*Spectrum.ScoreVariance);
	 * 
	 * this.OddsScore = this.NullDistribution.PValue(this.Score);
	 * 
	 * 
	 * //System.out.println("PValue: " + this.OddsScore + " PDF: " +
	 * this.NullDistribution.ApproximatePDF(this.Score)); //DistributionType
	 * Test = new
	 * DistributionType(this.AlignedPeaks.length*Spectrum.MeanScore,this
	 * .AlignedPeaks.length*Spectrum.ScoreVariance);
	 * 
	 * //System.out.println("Exact OddsScore: " + this.OddsScore +
	 * ", Approximate OddsScore: " + Test.PValue(this.Score));
	 * //Utils.WaitForEnter(); }
	 */
	public void SetScore(double Value) {
		this.Score = Value;
		// this.OddsScore = this.NullDistribution.ApproximatePDF(this.Score);
		// this.OddsScore = this.ContainingNode.GetNull().PValue(this.Score);
		// System.out.println("SetScore: PValue: " + this.OddsScore + " PDF: " +
		// this.NullDistribution.ApproximatePDF(this.Score));

	}

	public double GetScore() {
		return this.Score;
	}

	/*
	 * public double GetOddsScore() { DistributionType NullDistribution =
	 * this.ContainingNode.GetNull(); if (NullDistribution == null) {
	 * NullDistribution = new DistributionType(this.AlignedPeaks.length,
	 * this.Spectrum.PRMs, this.Spectrum.PRMScores);
	 * this.ContainingNode.SetNull(NullDistribution); } this.OddsScore =
	 * this.ContainingNode.GetNull().PValue(this.Score); return this.OddsScore;
	 * }
	 */

	// public void SetSequence(String Seq)
	// {
	// this.Sequence = Seq;
	// }

	// public String GetSequence()
	// {
	// return this.Sequence;
	// }
	public void SetAnchorSequence(String Seq) {
		this.AnchorSequence = Seq;
		this.AnchorPRMs = AntibodyUtils.GetSeqPRMScaled(this.AnchorSequence);
	}

	public void SetAnchorPRMs(int[] PRMs) {
		this.AnchorPRMs = PRMs;
		// for(int i = this.AnchorPRMs.length - 1; i >= 0; --i)
		// this.AnchorPRMs[i] -= this.AnchorPRMs[0];
	}

	public String GetAnchorSequence() {
		return this.AnchorSequence;
	}

	public Seed GetAnchorSeed() {
		return this.AnchorSeed;
	}

	public int[] GetAnchorPRMs() {
		return this.AnchorPRMs;
	}

	public void SetAlignedPeaksScaled(int[] Peaks) {
		this.AlignedPeaks = Peaks;
	}

	public void SetAlignedPeaksUnScaled(double[] Peaks) {
		this.AlignedPeaks = AntibodyUtils.ScaleArray(Peaks);
	}

	public int[] GetAlignedPeaksScaled() {
		return this.AlignedPeaks;
	}

	public double[] GetAlignedPeaksUnScaled() {
		return AntibodyUtils.UnScaleArray(this.AlignedPeaks);
	}

	public void SetAlignedSequencePeaksScaled(int[] Peaks) {
		this.AlignedSequence = Peaks;
	}

	public void SetAlignedSequencePeaksUnScaled(double[] Peaks) {
		this.AlignedSequence = AntibodyUtils.ScaleArray(Peaks);
	}

	public int[] GetAlignedSequencePeaksScaled() {
		return this.AlignedSequence;
	}

	public double[] GetAlignedSequencePeaksUnScaled() {
		return AntibodyUtils.UnScaleArray(this.AlignedPeaks);
	}

	public void SetSpectrum(PRMSpectrum Spectrum) {
		this.Spectrum = Spectrum;
	}

	public PRMSpectrum GetSpectrum() {
		return this.Spectrum;
	}

	public void SetSpectrumNode(SpectrumNode Node) {
		this.ContainingNode = Node;
	}

	public SpectrumNode GetSpectrumNode() {
		return this.ContainingNode;
	}

	// public void SetHMMAnnotation(HMMAlignment Annotation) {
	// this.HMMAnnotation = Annotation;
	// }

	// public HMMAlignment GetHMMAnnotation() {
	// return this.HMMAnnotation;
	// }

	public int GetScaledMassShift() {
		// System.out.println("MassShift = " + this.AlignedPeaks[0] + " - (" +
		// this.AlignedPeaks[0] + " - " + this.AlignedSequence[0] + ") = " +
		// (this.AlignedPeaks[0] - (this.AlignedPeaks[0] -
		// this.AlignedSequence[0])));
		// return this.GetAlignedPeaksScaled()[0] -
		// (this.GetAlignedPeaksScaled()[0] -
		// this.GetAlignedSequencePeaksScaled()[0]);
		return this.AlignedSequence[0] - this.AlignedPeaks[0];
	}

	public static MSAlignmentType[] GetBestNonRedundant(int NumBest,
			MSAlignmentType[] Alignments) {
		ArrayList TempList = new ArrayList();
		double CurrWorstScore = Double.MAX_VALUE;
		int CurrWorstIndex = -1;

		for (int i = 0; i < Alignments.length; ++i) {
			if (TempList.size() < NumBest) // If we haven't filled the list yet
			{
				TempList.add(Alignments[i]);
				if (Alignments[i].Score < CurrWorstScore) {
					CurrWorstScore = Alignments[i].Score;
					CurrWorstIndex = TempList.size() - 1;
				}
			} else if (Alignments[i].Score > CurrWorstScore) // If this is
			// better than
			// the current
			// worst
			{
				MSAlignmentType Temp = (MSAlignmentType) (TempList
						.get(CurrWorstIndex));
				SpectrumNode ParentTemp = Temp.GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				Temp.SetSpectrumNode(null);

				TempList.remove(CurrWorstIndex);

				TempList.add(Alignments[i]);
				CurrWorstIndex = MSAlignmentType.GetWorstIndex(TempList);
				CurrWorstScore = ((MSAlignmentType) (TempList
						.get(CurrWorstIndex))).Score;
			} else {
				SpectrumNode ParentTemp = Alignments[i].GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				Alignments[i].SetSpectrumNode(null);
			}

		}

		MSAlignmentType[] FinalList = new MSAlignmentType[TempList.size()];
		for (int i = 0; i < TempList.size(); ++i) {
			FinalList[i] = (MSAlignmentType) (TempList.get(i));
		}

		return FinalList;
	}

	public static MSAlignmentType[] GetBestOddsCutOff(double CutOff,
			MSAlignmentType[] Alignments, int MaxNum) {
		ArrayList TempList = new ArrayList();
		if (MaxNum == -1)
			MaxNum = Integer.MAX_VALUE;

		for (int i = 0; i < Alignments.length; ++i) {
			if (Alignments[i].OddsScore < CutOff && TempList.size() < MaxNum) {
				boolean Added = false;
				if (TempList.size() == 0) {
					TempList.add(Alignments[i]);
					Added = true;
				}
				for (int j = 0; j < TempList.size(); ++j) {
					MSAlignmentType CurrAlignment = (MSAlignmentType) (TempList
							.get(j));
					if (Alignments[i].OddsScore < CurrAlignment.OddsScore) {
						TempList.add(j, Alignments[i]);
						Added = true;
						break;
					}
				}
				if (!Added)
					TempList.add(Alignments[i]);

			} else if (Alignments[i].OddsScore < CutOff) {
				MSAlignmentType WorstAlignment = (MSAlignmentType) (TempList
						.get(TempList.size() - 1));
				if (Alignments[i].OddsScore < WorstAlignment.OddsScore) {
					TempList.remove(TempList.size() - 1);
					SpectrumNode ParentTemp = WorstAlignment.GetSpectrumNode();
					ParentTemp.SetMSAlignment(null);
					WorstAlignment.SetSpectrumNode(null);

					boolean Added = false;
					if (TempList.size() == 0) {
						TempList.add(Alignments[i]);
						Added = true;
					}
					for (int j = 0; j < TempList.size(); ++j) {
						MSAlignmentType CurrAlignment = (MSAlignmentType) (TempList
								.get(j));
						if (Alignments[i].OddsScore < CurrAlignment.OddsScore) {
							TempList.add(j, Alignments[i]);
							Added = true;
							break;
						}
					}
					if (!Added)
						TempList.add(Alignments[i]);
				}
			} else {
				SpectrumNode ParentTemp = Alignments[i].GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				Alignments[i].SetSpectrumNode(null);
			}
		}
		MSAlignmentType[] FinalList = new MSAlignmentType[TempList.size()];
		for (int i = 0; i < TempList.size(); ++i) {
			FinalList[i] = (MSAlignmentType) (TempList.get(i));
		}

		return FinalList;
	}

	/**
	 * Reduces a list of MSAlignmentType alignments to a maximum number above
	 * the score cutoff.
	 * 
	 * @param CutOff
	 *            The minimum score of an alignment kept
	 * @param Alignments
	 *            The array of MSAlignmentTypes
	 * @param MaxNum
	 *            The maximum number of alignments to keep
	 * @return Returns a sorted list of MSAlignmentType alignments.
	 */
	public static MSAlignmentType[] GetBestScoreCutOff(double CutOff,
			MSAlignmentType[] Alignments, int MaxNum) {

		ArrayList TempList = new ArrayList();
		if (MaxNum == -1)
			MaxNum = Integer.MAX_VALUE;

		for (int i = 0; i < Alignments.length; ++i) {
			if (Alignments[i] == null)
				continue;

			if (Alignments[i].GetSpectrumNode() == null) {
				System.out
						.println("WTF: msalignment has no parent node - Score: "
								+ Alignments[i].Score + "!!!");
				Utils.WaitForEnter();
			}
			if (Alignments[i].Score > CutOff && TempList.size() < MaxNum) {
				boolean Added = false;
				if (TempList.size() == 0) {
					TempList.add(Alignments[i]);
					Added = true;
				}
				for (int j = 0; j < TempList.size(); ++j) {
					MSAlignmentType CurrAlignment = (MSAlignmentType) (TempList
							.get(j));
					if (Alignments[i].Score > CurrAlignment.Score) {
						TempList.add(j, Alignments[i]);
						Added = true;
						break;
					}
				}
				if (!Added)
					TempList.add(Alignments[i]);

			} else if (Alignments[i].Score > CutOff) {
				MSAlignmentType WorstAlignment = (MSAlignmentType) (TempList
						.get(TempList.size() - 1));
				if (Alignments[i].Score > WorstAlignment.Score) {
					TempList.remove(TempList.size() - 1);
					SpectrumNode ParentTemp = WorstAlignment.GetSpectrumNode();
					if (ParentTemp != null) {
						if (ParentTemp.HMMAnnotation != null) {
							System.out
									.println("We are removing this spectrum,but it has an HMM annotation, so we should not reset the parent!!");
							WorstAlignment.Print();
							ParentTemp.DebugPrint();
							Utils.WaitForEnter();
						} else {
							ParentTemp.SetMSAlignment(null);
							ParentTemp.SetHMMAnnotation(null);
							WorstAlignment.SetSpectrumNode(null);
						}
					} else {
						System.out
								.println("Removing a spectrum with a good alignment ("
										+ WorstAlignment.Score
										+ ") but not better than the new guy ("
										+ Alignments[i].Score + ")");
					}

					boolean Added = false;
					if (TempList.size() == 0) {
						TempList.add(Alignments[i]);
						Added = true;
					}
					for (int j = 0; j < TempList.size(); ++j) {
						MSAlignmentType CurrAlignment = (MSAlignmentType) (TempList
								.get(j));
						if (Alignments[i].Score > CurrAlignment.Score) {
							TempList.add(j, Alignments[i]);
							Added = true;
							break;
						}
					}
					if (!Added)
						TempList.add(Alignments[i]);
				} else {
					SpectrumNode ParentTemp = Alignments[i].GetSpectrumNode();
					if (ParentTemp != null) {
						if (ParentTemp.HMMAnnotation != null) {
							System.out
									.println("We are removing ths spectrum,but it has an HMM annotation, so we should not reset the parent!!");
							Alignments[i].Print();
							ParentTemp.DebugPrint();
							Utils.WaitForEnter();
						} else {
							ParentTemp.SetHMMAnnotation(null);
							ParentTemp.SetMSAlignment(null);
							Alignments[i].SetSpectrumNode(null);
						}
					} else {
						System.out
								.println("Removing a spectrum with a good alignment ("
										+ WorstAlignment.Score
										+ ") but not better than the new guy ("
										+ Alignments[i].Score + ")");
					}
					Alignments[i].SetSpectrumNode(null);

				}

			} else {

				SpectrumNode ParentTemp = Alignments[i].GetSpectrumNode();
				if (ParentTemp != null) {
					if (ParentTemp.HMMAnnotation != null) {
						System.out
								.println("We are removing this spectrum,but it has an HMM annotation, so we should not reset the parent!!");
						Alignments[i].Print();
						ParentTemp.DebugPrint();
						Utils.WaitForEnter();
					} else {
						ParentTemp.SetMSAlignment(null);
						ParentTemp.SetHMMAnnotation(null);
						Alignments[i].SetSpectrum(null);
					}

				}
				Alignments[i].SetSpectrumNode(null);
			}
		}
		MSAlignmentType[] FinalList = new MSAlignmentType[TempList.size()];
		for (int i = 0; i < TempList.size(); ++i) {
			FinalList[i] = (MSAlignmentType) (TempList.get(i));
		}

		// return MSAlignmentType.SortedByNumPeaks(FinalList);
		return FinalList;
	}

	public static MSAlignmentType[] SortedByNumPeaks(MSAlignmentType[] OrigList) {
		MSAlignmentType[] Ret = new MSAlignmentType[OrigList.length];
		for (int i = 0; i < Ret.length; ++i)
			Ret[i] = OrigList[i];
		for (int i = 0; i < Ret.length; ++i) {
			int MaxIndex = i;
			int MaxValue = Ret[i].AlignedPeaks.length;
			for (int j = i + 1; j < Ret.length; ++j) {
				if (Ret[j].AlignedPeaks.length > MaxValue) {
					MaxIndex = j;
					MaxValue = Ret[j].AlignedPeaks.length;
				}
			}
			MSAlignmentType Temp = Ret[i];
			Ret[i] = Ret[MaxIndex];
			Ret[MaxIndex] = Temp;
		}

		return Ret;

	}

	public static MSAlignmentType[] GetBestNoCutOff(
			MSAlignmentType[] Alignments, int MaxNum) {
		ArrayList TempList = new ArrayList();
		double WorstScore = Double.MAX_VALUE;
		if (MaxNum == -1)
			MaxNum = Integer.MAX_VALUE;

		for (int i = 0; i < Alignments.length; ++i) {
			if (TempList.size() < MaxNum) {
				boolean Added = false;
				if (TempList.size() == 0) {
					TempList.add(Alignments[i]);
					Added = true;
				}
				for (int j = 0; j < TempList.size(); ++j) {
					MSAlignmentType CurrAlignment = (MSAlignmentType) (TempList
							.get(j));
					if (Alignments[i].Score > CurrAlignment.Score) {
						TempList.add(j, Alignments[i]);
						Added = true;
						break;
					}
				}
				if (!Added)
					TempList.add(Alignments[i]);
				if (Alignments[i].Score < WorstScore)
					WorstScore = Alignments[i].Score;

			} else if (Alignments[i].Score > WorstScore) {
				MSAlignmentType WorstAlignment = (MSAlignmentType) (TempList
						.get(TempList.size() - 1));
				TempList.remove(TempList.size() - 1);
				SpectrumNode ParentTemp = WorstAlignment.GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				WorstAlignment.SetSpectrumNode(null);

				boolean Added = false;
				if (TempList.size() == 0) {
					TempList.add(Alignments[i]);
					Added = true;
				}
				for (int j = 0; j < TempList.size(); ++j) {
					MSAlignmentType CurrAlignment = (MSAlignmentType) (TempList
							.get(j));
					if (Alignments[i].Score > CurrAlignment.Score) {
						TempList.add(j, Alignments[i]);
						Added = true;
						break;
					}
				}
				if (!Added)
					TempList.add(Alignments[i]);
				WorstAlignment = (MSAlignmentType) (TempList.get(TempList
						.size() - 1));
				WorstScore = WorstAlignment.Score;
			} else {
				SpectrumNode ParentTemp = Alignments[i].GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				Alignments[i].SetSpectrumNode(null);
			}
		}
		MSAlignmentType[] FinalList = new MSAlignmentType[TempList.size()];
		for (int i = 0; i < TempList.size(); ++i) {
			FinalList[i] = (MSAlignmentType) (TempList.get(i));
		}

		return FinalList;
	}

	/**
	 * Returns a list of the top Odds scoring MSAlignments, in decreasing order
	 * of odds score.
	 * 
	 * @param NumBest
	 *            The number of alignments to keep
	 * @param Alignments
	 *            the list of alignments to filter
	 * @return Returns the filtered, sorted list of Alignments.
	 */
	public static MSAlignmentType[] GetBestNonRedundantOdds(int NumBest,
			MSAlignmentType[] Alignments) {
		ArrayList TempList = new ArrayList();
		double CurrWorstScore = 0.0;
		int CurrWorstIndex = -1;

		for (int i = 0; i < Alignments.length; ++i) {
			if (TempList.size() < NumBest) // If we haven't filled the list yet
			{
				boolean Added = false;
				if (TempList.size() == 0) {
					TempList.add(Alignments[i]);
					Added = true;
				}
				for (int j = 0; j < TempList.size(); ++j) {
					MSAlignmentType CurrAlignment = (MSAlignmentType) (TempList
							.get(j));
					if (Alignments[i].OddsScore < CurrAlignment.OddsScore) {
						TempList.add(j, Alignments[i]);
						Added = true;
						break;
					}
				}
				if (!Added)
					TempList.add(Alignments[i]);
				if (Alignments[i].OddsScore > CurrWorstScore) {
					CurrWorstScore = Alignments[i].OddsScore;

				}
				CurrWorstIndex = TempList.size() - 1;
			} else if (Alignments[i].OddsScore < CurrWorstScore) // If this is
			// better
			// than the
			// current
			// worst
			{
				MSAlignmentType Temp = (MSAlignmentType) (TempList
						.get(CurrWorstIndex));
				SpectrumNode ParentTemp = Temp.GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				Temp.SetSpectrumNode(null);

				TempList.remove(CurrWorstIndex);
				boolean Added = false;
				for (int j = 0; j < TempList.size(); ++j) {
					MSAlignmentType CurrAlignment = (MSAlignmentType) (TempList
							.get(j));
					if (Alignments[i].OddsScore < CurrAlignment.OddsScore) {
						TempList.add(j, Alignments[i]);
						Added = true;
						break;
					}
				}
				if (!Added)
					TempList.add(Alignments[i]);

				// TempList.add(Alignments[i]);
				CurrWorstIndex = TempList.size() - 1;
				CurrWorstScore = ((MSAlignmentType) (TempList
						.get(CurrWorstIndex))).OddsScore;
			} else {
				SpectrumNode ParentTemp = Alignments[i].GetSpectrumNode();
				ParentTemp.SetMSAlignment(null);
				Alignments[i].SetSpectrumNode(null);
			}

		}

		MSAlignmentType[] FinalList = new MSAlignmentType[TempList.size()];
		for (int i = 0; i < TempList.size(); ++i) {
			FinalList[i] = (MSAlignmentType) (TempList.get(i));
		}

		return FinalList;
	}

	public static int GetWorstIndex(ArrayList AlignmentList) {
		int WorstIndex = -1;
		double WorstScore = Double.MAX_VALUE;
		for (int i = 0; i < AlignmentList.size(); ++i) {
			MSAlignmentType Temp = (MSAlignmentType) (AlignmentList.get(i));
			if (Temp.Score < WorstScore) {
				WorstScore = Temp.Score;
				WorstIndex = i;
			}
		}
		return WorstIndex;
	}

	public static int GetWorstIndexOdds(ArrayList AlignmentList) {
		int WorstIndex = -1;
		double WorstScore = 0.0;
		for (int i = 0; i < AlignmentList.size(); ++i) {
			MSAlignmentType Temp = (MSAlignmentType) (AlignmentList.get(i));
			if (Temp.OddsScore > WorstScore) {
				WorstScore = Temp.OddsScore;
				WorstIndex = i;
			}
		}
		return WorstIndex;
	}

	public boolean ExtendsSeed(double Tolerance) {
		boolean LocalDebug = false;
		int SuffixMass = this.Spectrum.PrecursorMass
				- this.AlignedPeaks[this.AlignedPeaks.length - 1];

		boolean HasGap = (this.AnchorSequence.indexOf("[") >= 0);
		int RemainingMass = this.AnchorPRMs[this.AnchorPRMs.length - 1]
				- this.AlignedSequence[this.AlignedSequence.length - 1];
		if (Debug || LocalDebug) {
			System.out.println("Seeing if this Alignment extends teh Seed:");
			this.PrintPeaks();
			System.out.println("Suffix mass (" + this.Spectrum.PrecursorMass
					+ " - " + this.AlignedPeaks[this.AlignedPeaks.length - 1]
					+ ") = " + SuffixMass);
			System.out.println("SeqSuffix Mass ("
					+ this.AnchorPRMs[this.AnchorPRMs.length - 1] + " - "
					+ this.AlignedSequence[this.AlignedSequence.length - 1]
					+ " = " + RemainingMass);
		}

		if (SuffixMass <= 0 && !HasGap) {
			if (LocalDebug)
				System.out.println("NOPE!");
			return false;
		}
		if (Debug || LocalDebug) {
			System.out.println("Remaining Mass: " + RemainingMass);
			// Utils.WaitForEnter();
		}
		if (AntibodyUtils
				.EquivalentScaled(SuffixMass, RemainingMass, Tolerance)
				&& HasGap) {
			if (LocalDebug) {
				System.out
						.println("**We have a spectrum which ends at the same place as the seed, but there is a gap in the Seed");
				System.out.println("Seed: " + this.AnchorSequence);
			}
			return true;

		}
		if (SuffixMass - RemainingMass >= -1 * Tolerance
				* AntibodyUtils.MASS_SCALE) {
			if (LocalDebug)
				System.out.println("YEP!");
			return true;
		}

		if (LocalDebug)
			System.out.println("NOPE!");
		return false;

	}

	public boolean ExtendsSeedLeft(double Tolerance) {
		boolean LocalDebug = false;
		int PrefixMass = this.AlignedPeaks[0];

		boolean HasGap = (this.AnchorSequence.indexOf("[") >= 0);
		int RemainingMass = this.AlignedSequence[0] - this.AnchorPRMs[0];
		if (Debug || LocalDebug) {
			System.out.println("Seeing if this Alignment extends teh Seed:");
			this.PrintPeaks();
			System.out.println("Prefix mass (" + this.AlignedPeaks[0] + ") = "
					+ PrefixMass);
			System.out.println("RemainingMass (" + this.AlignedSequence[0]
					+ " - " + this.AnchorPRMs[0] + " = " + RemainingMass);
		}
		if (PrefixMass <= 0 && !HasGap) {
			if (LocalDebug)
				System.out.println("NOPE!");
			return false;
		}

		// int RemainingMass =
		// (int)(AntibodyUtils.GetSeqMass(RemainingSequence)*AntibodyUtils.MASS_SCALE);
		if (Debug || LocalDebug) {
			System.out.println("Remaining Mass: " + RemainingMass);
			// Utils.WaitForEnter();
		}
		if (AntibodyUtils
				.EquivalentScaled(PrefixMass, RemainingMass, Tolerance)
				&& HasGap) {
			if (Debug || LocalDebug) {
				System.out
						.println("**We have a spectrum which ends at the same place as the seed, but there is a gap in the Seed");
				System.out.println("Seed: " + this.AnchorSequence);
			}
			return true;

		}
		if (PrefixMass - RemainingMass >= -1 * Tolerance
				* AntibodyUtils.MASS_SCALE) {
			if (LocalDebug)
				System.out.println("YEP!");
			return true;
		}

		if (LocalDebug)
			System.out.println("NOPE!");
		return false;

	}

	public void Print() {
		System.out.println("PRMSpectrum: " + this.Spectrum.FileName + ":"
				+ this.Spectrum.SpecIndex);
		System.out.println("Anchored Seq: " + this.AnchorSequence);
		System.out.println("Score: " + Score);

	}

	public void PrintPeaks() {
		this.Print();
		System.out.println("Peak     SeqPeak");
		for (int i = 0; i < this.AlignedPeaks.length; ++i) {
			System.out.println(this.AlignedPeaks[i] + " "
					+ this.AlignedSequence[i]);
		}
	}

	public String toString() {
		String AlignedPeaksString = "";
		String SeqPeaksString = Utils.IntArrayToString(this.AlignedSequence);
		for (int i = 0; i < this.AlignedPeaks.length; ++i)
			AlignedPeaksString += this.AlignedPeaks[i] + " ";
		String Line = this.Spectrum.FileName + ":" + this.Spectrum.SpecIndex
				+ "\n";
		Line += " Score: " + this.Score + ", AlignedPeaks: "
				+ AlignedPeaksString + "\n";
		return Line;
	}

	public String toStringDebug() {
		String Line = this.toString();
		String PeakString = "";
		String ScoreString = "";
		for (int i = 0; i < this.Spectrum.PRMs.length; ++i) {
			PeakString += this.Spectrum.PRMs[i] + " ";
			ScoreString += this.Spectrum.PRMScores[i] + " ";
		}
		Line += PeakString + "\n";
		Line += ScoreString + "\n";
		return Line;
	}

	public static MSAlignmentType[] ConvertToArray(ArrayList Alignments) {
		MSAlignmentType[] Ret = new MSAlignmentType[Alignments.size()];
		for (int i = 0; i < Ret.length; ++i) {
			Ret[i] = (MSAlignmentType) (Alignments.get(i));
		}
		return Ret;
	}

	public String toStringDebug(Hashtable Oracle) {
		if (Oracle == null)
			return this.toStringDebug();
		String Line = this.toString();
		String PeakString = "";
		String ScoreString = "";
		for (int i = 0; i < this.Spectrum.PRMs.length; ++i) {
			PeakString += this.Spectrum.PRMs[i] + " ";
			ScoreString += this.Spectrum.PRMScores[i] + " ";
		}
		Line += PeakString + "\n";
		Line += ScoreString + "\n";
		String HashKey = Utils.GetBaseNameNoExtension(this.Spectrum.FileName)
				+ "_" + this.Spectrum.SpecIndex;

		if (Oracle != null && Oracle.containsKey(HashKey)) {
			InspectAnnotation Ann = (InspectAnnotation) (Oracle.get(HashKey));
			Line += "Oracle: " + Ann.Annotation + "\n";
		} else
			Line += "Oracle: Unknown\n";
		return Line;
	}

	public static MSAlignmentType[] ClusterRankAlignments(
			MSAlignmentType[] OldAlignments, double Tolerance) {
		if (OldAlignments.length <= 1)
			return OldAlignments;
		boolean LocalDebug = false;

		int[][] All2All = new int[OldAlignments.length][OldAlignments.length];
		// int MaxClusterSize = 0;

		for (int i = 0; i < OldAlignments.length; ++i) {
			for (int j = i; j < OldAlignments.length; ++j) {
				All2All[i][j] = (int) (AntibodyUtils
						.AlignmentScoreLengthWithMassShiftRight(
								OldAlignments[i], OldAlignments[j], Tolerance));
				All2All[j][i] = All2All[i][j];

			}
		}

		double MaxClusterScore = 0;
		int MaxClusterCenter = 1;
		int MaxClusterSize = 0;
		boolean[] BestKept = null;
		for (int i = 0; i < OldAlignments.length; ++i) {
			int[] CurrScores = All2All[i];
			double CurrClusterScore = OldAlignments[i].Score;
			int CurrClusterSize = 1;

			int[] RankIndexes = AntibodyUtils
					.RankIntArrayDecreasing(CurrScores);
			if (LocalDebug) {
				System.out.println("Considring cluster centered at [" + i
						+ "]:" + OldAlignments[i].toString());
				System.out.println(Utils.IntArrayToString(CurrScores));
				System.out.println("Ranks: "
						+ Utils.IntArrayToString(RankIndexes));
			}
			boolean[] Kept = new boolean[CurrScores.length];
			Kept[i] = true;
			for (int j = 0; j < RankIndexes.length; ++j) {
				if (i == RankIndexes[j]
						|| All2All[i][RankIndexes[j]] < AntibodyUtils.MSALIGN_MIN_SIMILARITY)
					continue;
				if (LocalDebug)
					System.out.println("Considering new element ["
							+ RankIndexes[j] + "] at rank " + j
							+ " with distance " + All2All[i][RankIndexes[j]]);
				Kept[RankIndexes[j]] = true;
				for (int k = 0; k < Kept.length; ++k) {
					if (Kept[k]) {
						if (All2All[k][RankIndexes[j]] < AntibodyUtils.MSALIGN_MIN_SIMILARITY) {
							Kept[RankIndexes[j]] = false;
							break;
						}
					}
				}
				if (Kept[RankIndexes[j]]) {
					if (LocalDebug)
						System.out.println("CLOSE ENOUGH, Adding!!");
					CurrClusterScore += OldAlignments[RankIndexes[j]].Score;
					CurrClusterSize++;
				} else if (LocalDebug)
					System.out.println("NOT CLOSE ENOUGH!!");
			}
			if (CurrClusterScore > MaxClusterScore) {
				MaxClusterScore = CurrClusterScore;
				MaxClusterCenter = i;
				MaxClusterSize = CurrClusterSize;
				BestKept = Kept;
				if (LocalDebug) {
					System.out.println("***Best Cluster So Far: "
							+ MaxClusterCenter);
					System.out.println("Size: " + CurrClusterSize + ", Score: "
							+ MaxClusterScore);
					for (int k = 0; k < BestKept.length; ++k)
						if (BestKept[k])
							System.out.println(OldAlignments[k].toString());
					Utils.WaitForEnter();
				}
			}
		}

		MSAlignmentType[] Ret = new MSAlignmentType[MaxClusterSize];
		int Index = 0;
		for (int i = 0; i < BestKept.length; ++i) {
			if (BestKept[i]) {
				Ret[Index] = OldAlignments[i];
				Index++;
			}
		}

		if (LocalDebug) {
			System.out.println("Final cluster: " + Ret.length);
			System.out.println(Utils.BooleanArrayToString(BestKept));
			Utils.WaitForEnter();
		}
		return Ret;
	}

	public static MSAlignmentType[] ClusterRankAlignmentsLeft(
			MSAlignmentType[] OldAlignments, double Tolerance) {
		if (OldAlignments.length <= 1)
			return OldAlignments;
		boolean LocalDebug = false;

		int[][] All2All = new int[OldAlignments.length][OldAlignments.length];
		// int MaxClusterSize = 0;

		for (int i = 0; i < OldAlignments.length; ++i) {
			for (int j = i; j < OldAlignments.length; ++j) {
				All2All[i][j] = (int) (AntibodyUtils
						.AlignmentScoreLengthWithMassShiftLeft(
								OldAlignments[i], OldAlignments[j], Tolerance));
				All2All[j][i] = All2All[i][j];

			}
		}

		double MaxClusterScore = 0;
		int MaxClusterCenter = 1;
		int MaxClusterSize = 0;
		boolean[] BestKept = null;
		// boolean LocalDebug = true;
		for (int i = 0; i < OldAlignments.length; ++i) {
			int[] CurrScores = All2All[i];
			double CurrClusterScore = OldAlignments[i].Score;
			int CurrClusterSize = 1;
			int[] RankIndexes = AntibodyUtils
					.RankIntArrayDecreasing(CurrScores);
			boolean[] Kept = new boolean[CurrScores.length];
			Kept[i] = true;
			if (LocalDebug)
				System.out.println("Considering Center at " + i + "="
						+ OldAlignments[i].Score);
			for (int j = 0; j < RankIndexes.length; ++j) {
				if (i == RankIndexes[j]
						|| All2All[i][RankIndexes[j]] < AntibodyUtils.MSALIGN_MIN_SIMILARITY)
					continue;
				if (LocalDebug)
					System.out.println("Comparing to rank " + j + "="
							+ OldAlignments[RankIndexes[j]].Score);
				Kept[RankIndexes[j]] = true;
				for (int k = 0; k < Kept.length; ++k) {
					if (LocalDebug)
						System.out.println("Is item " + k + "="
								+ OldAlignments[k] + " in the cluster?");
					if (Kept[k]) {
						if (LocalDebug)
							System.out.println("Yup!");
						if (All2All[k][RankIndexes[j]] < AntibodyUtils.MSALIGN_MIN_SIMILARITY) {
							if (LocalDebug)
								System.out.println("Too far from this guy!! ("
										+ All2All[k][RankIndexes[j]] + ")");
							Kept[RankIndexes[j]] = false;
							break;
						} else if (LocalDebug)
							System.out.println("Close enough! ("
									+ All2All[k][RankIndexes[j]] + ")");
					}
				}
				if (Kept[RankIndexes[j]]) {
					CurrClusterScore += OldAlignments[RankIndexes[j]].Score;
					CurrClusterSize++;
				}
			}
			if (CurrClusterScore > MaxClusterScore) {
				MaxClusterScore = CurrClusterScore;
				MaxClusterCenter = i;
				MaxClusterSize = CurrClusterSize;
				BestKept = Kept;
				if (LocalDebug) {
					System.out.println("***Best Cluster So Far: "
							+ MaxClusterCenter);
					System.out.println("Size: " + CurrClusterSize + ", Score: "
							+ MaxClusterScore);
					for (int k = 0; k < BestKept.length; ++k)
						if (BestKept[k])
							System.out.println(OldAlignments[k].toString());
					Utils.WaitForEnter();
				}
			}
		}

		MSAlignmentType[] Ret = new MSAlignmentType[MaxClusterSize];
		int Index = 0;
		for (int i = 0; i < BestKept.length; ++i) {
			if (BestKept[i]) {
				Ret[Index] = OldAlignments[i];
				Index++;
			}
		}
		return Ret;
	}

}
