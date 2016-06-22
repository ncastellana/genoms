package antibody;

import java.util.ArrayList;

import basicUtils.Utils;

public class ConsensusPRMSpectrum {

	public int[] ScaledPeaks;
	public double[] PeakScores;

	public Object[] SupportingSpectra; // SupportingSpectra[i] = {list of
										// SpectrumNode indexes};
	public Object[] OverlappingSpectra;

	// These are not used currently, only for scoring based on
	public double[] IntervalScores;

	public SpectrumHMM Model;

	/**
	 * Often when we reconstruct a sequence from a consensus spectrum, we want
	 * scores only for the peaks that we used. Each PRM in UsedPRMs should
	 * appear in ScaledPeaks. Though not all peaks in ScaledPeaks will be in
	 * UsedPRMs. If an interval in UsedPRMs is covered by more than 1 interval
	 * in ScaledPeaks, the maximum probability is used for the probability of
	 * the gap
	 * 
	 * @param UsedPRMs
	 */
	public double[] GetIntervalScores(int[] UsedPRMs) {
		int[] StateNums = new int[UsedPRMs.length];
		for (int i = 0; i < StateNums.length; ++i) {
			int CurrMass = UsedPRMs[i];
			boolean Found = false;
			for (int j = 0; j < this.ScaledPeaks.length; ++j) {
				if (this.ScaledPeaks[j] == CurrMass) {
					Found = true;
					StateNums[i] = j;
					break;
				}
			}
			if (!Found) {
				System.err.println("ERROR: Could not find mass " + CurrMass
						+ " in our consensus!!");
				for (int j = 0; j < this.ScaledPeaks.length; ++j) {
					System.out.println(this.ScaledPeaks[j]);
				}
				Utils.WaitForEnter();
			}
		}
		double[] Ret = new double[UsedPRMs.length - 1];
		for (int i = 0; i < Ret.length; ++i) {
			Ret[i] = 0.0;
			int LeftPRMIndex = StateNums[i];
			int RightPRMIndex = StateNums[i + 1];
			for (int j = LeftPRMIndex; j < RightPRMIndex; ++j) {
				Ret[i] = Math.max(Ret[i], this.IntervalScores[j]);
			}
		}
		return Ret;
	}

	/**
	 * Forces the alignment like: A -------- B --------
	 * 
	 * Must begin at BStart and end at AEnd. Single PRMSkips are allowed
	 * 
	 * @param PRMA
	 * @param ScoresA
	 * @param PRMB
	 * @param ScoresB
	 * @return
	 */
	public static ConsensusPRMSpectrum MergePRMsForceAEndBStart(int[] PRMA,
			double[] ScoresA, Object[] OverlappingSpectraA,
			Object[] SupportingSpectraA, int[] PRMB, double[] ScoresB,
			Object[] OverlappingSpectraB, Object[] SupportingSpectraB,
			double Tolerance) {
		int[][] Alignment = new int[PRMA.length][PRMB.length];
		int[][][] PrevSteps = new int[PRMA.length][PRMB.length][2];
		boolean LocalDebug = false;
		if (LocalDebug) {
			System.out.println("ALIGNING, FORCING AStart BEnd");
			Utils.WaitForEnter();
		}
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
				if (j == 1)
					MaxScore = 1;
				PrevSteps[i][j][0] = -1;
				PrevSteps[i][j][1] = -1;

				// Consider matching PRMA[i] with PRMB[j], look at previous
				// matchings, allowing 1 skipped PRM
				for (int Previ = Math.max(0, i - 2); Previ < i; Previ++) {
					for (int Prevj = Math.max(0, j - 2); Prevj < j; Prevj++) {
						// if(LocalDebug)
						// System.out.println("   prevs" + PRMA[Previ] + "(" +
						// Previ + ")" + " and " + PRMB[Prevj] + "(" + Prevj +
						// ")");
						if (Alignment[Previ][Prevj] <= 0) {
							// System.out.println("   -no alignment ends here");
							continue;
						}
						int Diff1 = PRMA[i] - PRMA[Previ];
						int Diff2 = PRMB[j] - PRMB[Prevj];
						// int Penalty = Math.max(i - Previ - 1, j - Prevj - 1);

						if (AntibodyUtils.EquivalentScaled(Diff1, Diff2,
								Tolerance)) {
							if (Alignment[Previ][Prevj] + 1 > MaxScore) {
								MaxScore = Alignment[Previ][Prevj] + 1;
								PrevSteps[i][j][0] = Previ;
								PrevSteps[i][j][1] = Prevj;
								// if(LocalDebug)
								// {
								// System.out.println("BestScore to [" + PRMA[i]
								// + "(" + i + ")," + PRMB[j] + "(" + j +
								// ")] is from [" + PRMA[Previ] + "(" + Previ +
								// ")," + PRMB[Prevj] + "(" + Prevj + ")]");
								// }
							}
						}
						// else if(LocalDebug)
						// System.out.println("   -not close enough");
					}
				}
				Alignment[i][j] = MaxScore;
				// if(LocalDebug)
				// {
				// System.out.println("Best score: " + MaxScore);
				// Utils.WaitForEnter();
				// }
			}
		}
		for (int i = 0; i < PRMB.length; ++i) {
			if (Alignment[PRMA.length - 1][i] > TotalMaxScore) {
				TotalMaxScore = Alignment[PRMA.length - 1][i];
				BestCells[0] = PRMA.length - 1;
				BestCells[1] = i;
			}
			if (Alignment[PRMA.length - 2][i] > TotalMaxScore) {
				TotalMaxScore = Alignment[PRMA.length - 2][i];
				BestCells[0] = PRMA.length - 2;
				BestCells[1] = i;

			}
		}

		if (LocalDebug) {
			System.out.println("Max Score: " + TotalMaxScore + ", EndA: "
					+ BestCells[0] + " EndB: " + BestCells[1]);
			Utils.WaitForEnter();
		}
		if (BestCells[0] < PRMA.length - 2) {
			if (LocalDebug)
				System.out.println(" We don't end at end of A!!");
			return null;
		}

		if (TotalMaxScore <= 2) {
			if (LocalDebug)
				System.out.println("Alignment is only 1 prm!!");
			return null;
		}
		int MassShift = PRMA[BestCells[0]] - PRMB[BestCells[1]];
		// Sanity Check, B should start somewhere after A
		// if(MassShift < 0 && Math.abs(MassShift) >
		// Utils.LTQFTIonTolerance*Utils.MASS_SCALE)
		// return null;

		if (LocalDebug) {
			System.out.println("Seed should be shifted by " + MassShift);
			Utils.WaitForEnter();
		}

		int[] CurrCells = { BestCells[0], BestCells[1] };
		int Length = 1;
		int ASkips = PRMA.length - 1 - BestCells[0];
		int BSkips = 0;
		while (true) {
			if (LocalDebug) {
				System.out.println("PRMA[" + CurrCells[0] + "]="
						+ PRMA[CurrCells[0]]);
				System.out.println("PRMB[" + CurrCells[1] + "]="
						+ (PRMB[CurrCells[1]] + MassShift) + "("
						+ PRMB[CurrCells[1]] + ")");
			}
			int i = CurrCells[0];
			int j = CurrCells[1];
			if (PrevSteps[i][j][0] == -1 || PrevSteps[i][j][1] == -1)
				break;
			CurrCells[0] = PrevSteps[i][j][0];
			CurrCells[1] = PrevSteps[i][j][1];

			Length += 1;
			ASkips += (i - CurrCells[0] - 1);
			BSkips += (j - CurrCells[1] - 1);
		}
		BSkips += CurrCells[1];

		if (LocalDebug) {
			System.out.println("LenA: " + PRMA.length + ", LenB: "
					+ PRMB.length);
			System.out.println("ASkips: " + ASkips);
			System.out.println("BSkips: " + BSkips);
			System.out.println("Length: " + Length);
		}

		if (ASkips >= 2 || BSkips >= 2) {
			if (LocalDebug)
				System.out.println("NOOO!!!");
			return null;
		}

		return ConsensusPRMSpectrum.ConstructSpectrumFromAlignment(MassShift,
				PRMA, PRMB, ScoresA, ScoresB, SupportingSpectraA,
				SupportingSpectraB, OverlappingSpectraA, OverlappingSpectraB,
				Tolerance);

	}

	public static ConsensusPRMSpectrum MergeAnchors(int[] PRMA,
			double[] ScoresA, Object[] OverlappingSpectraA,
			Object[] SupportingSpectraA, int[] PRMB, double[] ScoresB,
			Object[] OverlappingSpectraB, Object[] SupportingSpectraB,
			double Tolerance) {
		boolean LocalDebug = false;
		// Try to merge according to the strongest alignment.
		ConsensusPRMSpectrum Ret = ConsensusPRMSpectrum.MergePRMsSkipPenalty(
				PRMA, ScoresA, OverlappingSpectraA, SupportingSpectraA, PRMB,
				ScoresB, OverlappingSpectraB, SupportingSpectraB, Tolerance);
		if (Ret != null) {
			if (LocalDebug)
				System.out.println("MERGEABLE BY Scores!!");
			return Ret;
		}
		if (LocalDebug) {
			System.out.println("Not Mergeable by scores!!");
			Utils.WaitForEnter();
		}
		// Try the easy case
		// A --------
		// ---------
		Ret = ConsensusPRMSpectrum.MergePRMsForceAEndBStart(PRMA, ScoresA,
				OverlappingSpectraA, SupportingSpectraA, PRMB, ScoresB,
				OverlappingSpectraB, SupportingSpectraB, Tolerance);
		if (Ret == null && LocalDebug) {
			System.out.println("Not Mergeable by simple strict!!");
			Utils.WaitForEnter();
		}
		return Ret;
	}

	public static ConsensusPRMSpectrum MergePRMs(int[] PRMA, double[] ScoresA,
			Object[] OverlappingSpectraA, Object[] SupportingSpectraA,
			int[] PRMB, double[] ScoresB, Object[] OverlappingSpectraB,
			Object[] SupportingSpectraB, double Tolerance) {
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

		if (TotalMaxScore < AntibodyUtils.MIN_OVERLAPPING_AA) {
			System.out.println("MergePRMs: TotalMaxScore " + TotalMaxScore
					+ " < " + AntibodyUtils.MIN_OVERLAPPING_AA);
			return null;
		}

		// Sanity check
		if (BestCells[0] != PRMA.length - 1 && BestCells[1] != PRMA.length) {
			System.out.println("BestCells[0] = " + BestCells[0]
					+ ", BestCells[1] = " + BestCells[1]);
			return null;
		}
		int MassShift = PRMA[BestCells[0]] - PRMB[BestCells[1]];

		// Sanity Check, B should start somewhere after A
		// if(MassShift < 0 && Math.abs(MassShift) >
		// Utils.LTQFTIonTolerance*Utils.MASS_SCALE)
		// return null;

		return ConsensusPRMSpectrum.ConstructSpectrumFromAlignment(MassShift,
				PRMA, PRMB, ScoresA, ScoresB, SupportingSpectraA,
				SupportingSpectraB, OverlappingSpectraA, OverlappingSpectraB,
				Tolerance);
	}

	public static ConsensusPRMSpectrum MergePRMsSkipPenalty(int[] PRMA,
			double[] ScoresA, Object[] OverlappingSpectraA,
			Object[] SupportingSpectraA, int[] PRMB, double[] ScoresB,
			Object[] OverlappingSpectraB, Object[] SupportingSpectraB,
			double Tolerance) {
		double[][] Alignment = new double[PRMA.length][PRMB.length]; // Score
		// it!
		int[][][] PrevSteps = new int[PRMA.length][PRMB.length][2];

		boolean LocalDebug = false;
		if (LocalDebug) {
			System.out.println("ALIGNING, Skip Penalty");
			Utils.WaitForEnter();
		}
		// Initialize the alignment. We allow starting anywhere.
		// Initialize guys to have the average PRM score
		for (int i = 0; i < PRMA.length; ++i) {
			// Alignment[i][0] = (ScoresA[i] + ScoresB[0])/2;
			Alignment[i][0] = 1;
			PrevSteps[i][0][0] = -1;
			PrevSteps[i][0][1] = -1;

		}
		// Initialize the alignment. We allow starting anywhere
		for (int i = 0; i < PRMB.length; ++i) {
			// Alignment[0][i] = (ScoresA[0] + ScoresB[i])/2;
			Alignment[0][i] = 1;
			PrevSteps[0][i][0] = -1;
			PrevSteps[0][i][1] = -1;
		}
		double TotalMaxScore = 0; // Best score so far
		int[] BestCells = new int[2]; // Best Cell so far

		for (int i = 1; i < PRMA.length; ++i) {
			for (int j = 1; j < PRMB.length; ++j) {
				// double MaxScore = (ScoresA[i] + ScoresB[j])/2; //The best
				// score is just these peaks aligned
				double MaxScore = 1; // The best score is just these peaks
				// aligned
				PrevSteps[i][j][0] = -1;
				PrevSteps[i][j][1] = -1;

				// Consider matching PRMA[i] with PRMB[j], look at previous
				// matchings, allowing any number of Skipped
				for (int Previ = Math.max(i - AntibodyUtils.MAX_SKIPPED_PRMs
						- 1, 0); Previ < i; Previ++) {
					for (int Prevj = Math.max(j
							- AntibodyUtils.MAX_SKIPPED_PRMs - 1, 0); Prevj < j; Prevj++) {
						int Diff1 = PRMA[i] - PRMA[Previ];
						int Diff2 = PRMB[j] - PRMB[Prevj];

						// The penalty is the likelihood of the peaks in the gap
						// missing. The peak
						// penalty is just the MISSING PEAK PENALTY from the
						// HMM.
						// double Penalty = Math.min(i-Previ-1,
						// j-Prevj-1)*AntibodyUtils.MISSING_PEAK_PENALTY;
						if (AntibodyUtils.EquivalentScaled(Diff1, Diff2,
								Tolerance)) {

							// double NewScore = Alignment[Previ][Prevj] +
							// (ScoresA[i] + ScoresB[j])/2 - Penalty;
							double NewScore = Alignment[Previ][Prevj] + 1;
							if (NewScore > MaxScore) {

								MaxScore = NewScore;
								PrevSteps[i][j][0] = Previ;
								PrevSteps[i][j][1] = Prevj;
							}
						}
					}
				}
				Alignment[i][j] = MaxScore;
				// if(MaxScore -
				// Math.min(PRMA.length-1-i,PRMB.length-1-j)*AntibodyUtils.MISSING_PEAK_PENALTY
				// > TotalMaxScore)
				if (MaxScore > TotalMaxScore) {
					// TotalMaxScore = MaxScore -
					// Math.min(PRMA.length-1-i,PRMB.length-1-j)*AntibodyUtils.MISSING_PEAK_PENALTY;
					TotalMaxScore = MaxScore;
					BestCells[0] = i;
					BestCells[1] = j;
				}
			}
		}
		if (LocalDebug) {
			System.out.println("Max Score: " + TotalMaxScore + ", EndA: "
					+ BestCells[0] + " EndB: " + BestCells[1]);
		}

		int MassShift = PRMA[BestCells[0]] - PRMB[BestCells[1]];

		if (LocalDebug) {
			System.out.println("Seed should be shifted by " + MassShift);
			Utils.WaitForEnter();
		}

		int[] CurrCells = { BestCells[0], BestCells[1] };
		int Length = 1;
		int ASkips = 0;
		int BSkips = 0;
		while (true) {
			if (LocalDebug) {
				System.out.println("PRMA[" + CurrCells[0] + "]="
						+ PRMA[CurrCells[0]]);
				System.out.println("PRMB[" + CurrCells[1] + "]="
						+ (PRMB[CurrCells[1]] + MassShift) + "("
						+ PRMB[CurrCells[1]] + ")");
			}
			int i = CurrCells[0];
			int j = CurrCells[1];
			if (PrevSteps[i][j][0] == -1 || PrevSteps[i][j][1] == -1)
				break;

			Length += 1;
			CurrCells[0] = PrevSteps[i][j][0];
			CurrCells[1] = PrevSteps[i][j][1];
			ASkips += (i - CurrCells[0] - 1);
			BSkips += (j - CurrCells[1] - 1);
		}
		if (LocalDebug) {
			System.out.println("LenA: " + PRMA.length + ", LenB: "
					+ PRMB.length);
			System.out.println("ASkips: " + ASkips);
			System.out.println("BSkips: " + BSkips);
			System.out.println("Length: " + Length);
		}

		if (ASkips > 2 || BSkips > 2) {
			if (LocalDebug)
				System.out.println("NOOO!!!");
			return null;
		}// If we have one of the easy cases, then we need a shorter length
			// A ------------
			// B ------------
		if (BestCells[0] >= PRMA.length - 2 && CurrCells[1] <= 1) {
			if (Length < AntibodyUtils.MIN_OVERLAPPING_AA + 1) {
				if (LocalDebug)
					System.out
							.println("We end at the end of A, and begin at the beginning of B, but length is too short: "
									+ Length);
				return null;
			}
		}
		// A --------------
		// B -------
		else if (BestCells[1] == PRMB.length - 1 && CurrCells[1] == 0) {
			if (Length < AntibodyUtils.MIN_OVERLAPPING_AA + 1) {
				if (LocalDebug)
					System.out.println("A contains B but length is too short: "
							+ Length);
				return null;
			}
		}
		// A --------
		// B --------------
		else if (BestCells[0] == PRMA.length - 1 && CurrCells[0] == 0) {
			if (Length < AntibodyUtils.MIN_OVERLAPPING_AA + 1) {
				if (LocalDebug)
					System.out.println("B contains A but length is too short: "
							+ Length);
				return null;
			}

		}

		// We don't have an easy alignment, so we need to have at least
		// MIN_OVERLAPPING_AA_INTERNAL
		else {
			if (Length < AntibodyUtils.MIN_OVERLAPPING_AA_INTERNAL + 1) {
				if (LocalDebug)
					System.out.println("Alignment is internal and too short: "
							+ Length);
				return null;
			}

		}
		return ConsensusPRMSpectrum.ConstructSpectrumFromAlignment(MassShift,
				PRMA, PRMB, ScoresA, ScoresB, SupportingSpectraA,
				SupportingSpectraB, OverlappingSpectraA, OverlappingSpectraB,
				Tolerance);
	}

	private static ConsensusPRMSpectrum ConstructSpectrumFromAlignment(
			int MassShift, int[] PRMA, int[] PRMB, double[] ScoresA,
			double[] ScoresB, Object[] SupportingSpectraA,
			Object[] SupportingSpectraB, Object[] OverlappingSpectraA,
			Object[] OverlappingSpectraB, double Tolerance) {
		ArrayList TempPRMs = new ArrayList();
		ArrayList PRMScores = new ArrayList();
		ArrayList SupportingSpectra = new ArrayList();
		ArrayList Attributor = new ArrayList();
		ArrayList AttIndex = new ArrayList();
		int AIndex = 0;
		int BIndex = 0;

		while (AIndex < PRMA.length && BIndex < PRMB.length) {
			if (AntibodyUtils.EquivalentScaled(PRMA[AIndex], PRMB[BIndex]
					+ MassShift, Tolerance)) {
				TempPRMs.add(new Integer(PRMA[AIndex]));
				PRMScores.add(new Double(
						(ScoresA[AIndex] + ScoresB[BIndex]) / 2));

				// Compute the total number of overlapping for this peak

				String[] TempA = (String[]) (SupportingSpectraA[AIndex]);
				String[] TempB = (String[]) (SupportingSpectraB[BIndex]);
				Attributor.add("AB");
				int[] Pos = { AIndex, BIndex };
				AttIndex.add(Pos);
				String[] CombinedSupporting = AntibodyUtils.GetUnion(TempA,
						TempB);
				SupportingSpectra.add(CombinedSupporting);

				AIndex++;
				BIndex++;

			} else if (PRMA[AIndex] < PRMB[BIndex] + MassShift) {
				TempPRMs.add(new Integer(PRMA[AIndex]));
				PRMScores.add(new Double(ScoresA[AIndex]));

				Attributor.add("A");
				int[] Pos = { AIndex, -1 };
				AttIndex.add(Pos);
				String[] CombinedSupporting = (String[]) (SupportingSpectraA[AIndex]);
				SupportingSpectra.add(CombinedSupporting);

				AIndex++;
			} else {
				TempPRMs.add(new Integer(PRMB[BIndex] + MassShift));
				PRMScores.add(new Double(ScoresB[BIndex]));

				Attributor.add("B");
				int[] Pos = { -1, BIndex };
				AttIndex.add(Pos);
				String[] CombinedSupporting = (String[]) (SupportingSpectraB[BIndex]);
				SupportingSpectra.add(CombinedSupporting);
				BIndex++;
			}
		}
		while (AIndex < PRMA.length) {
			TempPRMs.add(new Integer(PRMA[AIndex]));
			PRMScores.add(new Double(ScoresA[AIndex]));

			Attributor.add("A");
			int[] Pos = { AIndex, -1 };
			AttIndex.add(Pos);
			String[] CombinedSupporting = (String[]) (SupportingSpectraA[AIndex]);
			SupportingSpectra.add(CombinedSupporting);

			AIndex++;
		}
		while (BIndex < PRMB.length) {
			TempPRMs.add(new Integer(PRMB[BIndex] + MassShift));
			PRMScores.add(new Double(ScoresB[BIndex]));

			Attributor.add("B");
			int[] Pos = { -1, BIndex };
			AttIndex.add(Pos);
			String[] CombinedSupporting = (String[]) (SupportingSpectraB[BIndex]);
			SupportingSpectra.add(CombinedSupporting);

			BIndex++;
		}

		ConsensusPRMSpectrum Ret = new ConsensusPRMSpectrum();
		Ret.ScaledPeaks = AntibodyUtils.ConvertIntegerArrayList(TempPRMs);
		Ret.PeakScores = AntibodyUtils.ConvertDoubleArrayList(PRMScores);

		boolean LocalDebug = false;
		if (LocalDebug) {
			System.out.println("Attributions: "
					+ AntibodyUtils.StringArrayToString(AntibodyUtils
							.ConvertStringArrayList(Attributor)));
			String S = "";
			for (int i = 0; i < Attributor.size(); ++i) {
				int[] Pos = (int[]) (AttIndex.get(i));
				S += "[" + Pos[0] + "," + Pos[1] + "] ";
			}
			System.out.println("Indexes: " + S);

		}
		// On the first pass we computed the number of supporting spectra for
		// each peak, but now we need to figure out overlapping spectra
		Ret.SupportingSpectra = new Object[SupportingSpectra.size()];
		Ret.OverlappingSpectra = new Object[SupportingSpectra.size()];
		int PrevA = -1;
		int PrevB = -1;
		for (int i = 0; i < Ret.SupportingSpectra.length; ++i) {

			Ret.SupportingSpectra[i] = (String[]) (SupportingSpectra.get(i));
			if (LocalDebug) {
				System.out.println("Considering position: " + i);
				System.out.println("SupportingSpectra: "
						+ ((String[]) (Ret.SupportingSpectra[i])).length);
			}
			String Attribution = (String) (Attributor.get(i));
			int[] Pos = (int[]) (AttIndex.get(i));

			if (LocalDebug) {
				System.out.println("Attribution: " + Attribution);
				System.out.println("Pos: [" + Pos[0] + "," + Pos[1] + "]");
			}
			// If the peak was from both guys, then take the union of the two
			// overlapping spectra sets
			if (Attribution.compareTo("AB") == 0) {
				String[] TempA = (String[]) (OverlappingSpectraA[Pos[0]]);
				String[] TempB = (String[]) (OverlappingSpectraB[Pos[1]]);

				String[] CombinedOverlapping = AntibodyUtils.GetUnion(TempA,
						TempB);
				Ret.OverlappingSpectra[i] = (CombinedOverlapping);
				if (LocalDebug) {
					System.out.println("Ended at a shared peak");
					System.out.println("A overlapping: "
							+ AntibodyUtils.StringArrayToString(TempA));
					System.out.println("B overlapping: "
							+ AntibodyUtils.StringArrayToString(TempB));
					System.out.println("Combined: "
							+ AntibodyUtils
									.StringArrayToString(CombinedOverlapping));
				}

				PrevA = Pos[0];
				PrevB = Pos[1];
			}

			// If the peak was only from A, then find the previous B and next B
			else if (Attribution.compareTo("A") == 0) {

				if (PrevB == -1) // We have had a b yet, so the overlapping are
									// only the overlaping of A
				{
					if (LocalDebug)
						System.out.println("No previous Bs...");
					Ret.OverlappingSpectra[i] = (String[]) (OverlappingSpectraA[Pos[0]]);
				} else

				{
					int NextB = i + 1;
					String Attribution2 = "";
					if (NextB < Ret.SupportingSpectra.length)
						Attribution2 = (String) (Attributor.get(NextB));
					while (NextB < Ret.SupportingSpectra.length
							&& Attribution2.indexOf('B') < 0) {
						NextB += 1;
						if (NextB < Ret.SupportingSpectra.length)
							Attribution2 = (String) (Attributor.get(NextB));
					}
					if (NextB == Ret.SupportingSpectra.length) // There is no
																// next B peak,
																// so the
																// overalapping
																// are only the
																// overlapping
																// of A
					{
						if (LocalDebug)
							System.out.println("No previous Bs...");
						Ret.OverlappingSpectra[i] = (String[]) (OverlappingSpectraA[Pos[0]]);
					} else {
						int[] BPos = (int[]) AttIndex.get(NextB);
						String[] PrevOverlapping = (String[]) (OverlappingSpectraB[PrevB]);
						String[] NextOverlapping = (String[]) (OverlappingSpectraB[BPos[1]]);
						String[] TempB = AntibodyUtils.GetIntersection(
								PrevOverlapping, NextOverlapping);
						String[] TempA = (String[]) (OverlappingSpectraA[Pos[0]]);
						Ret.OverlappingSpectra[i] = AntibodyUtils.GetUnion(
								TempA, TempB);
					}
				}
				PrevA = Pos[0];
			}
			// If the peak was only from A, then find the previous B and next B
			else if (Attribution.compareTo("B") == 0) {
				if (PrevA == -1) // We have had a b yet, so the overlapping are
									// only the overlaping of A
				{
					if (LocalDebug)
						System.out.println("No previous As...");
					Ret.OverlappingSpectra[i] = (String[]) (OverlappingSpectraB[Pos[1]]);
				} else {
					int NextA = i + 1;
					String Attribution2 = "";
					if (NextA < Ret.SupportingSpectra.length)
						Attribution2 = (String) (Attributor.get(NextA));
					while (NextA < Ret.SupportingSpectra.length
							&& Attribution2.indexOf('A') < 0) {
						NextA += 1;
						if (NextA < Ret.SupportingSpectra.length)
							Attribution2 = (String) (Attributor.get(NextA));
					}
					if (NextA == Ret.SupportingSpectra.length) // There is no
																// next B peak,
																// so the
																// overalapping
																// are only the
																// overlapping
																// of A
					{
						if (LocalDebug)
							System.out.println("No following As...");
						Ret.OverlappingSpectra[i] = (String[]) (OverlappingSpectraB[Pos[1]]);
					} else {
						int[] APos = (int[]) AttIndex.get(NextA);
						String[] PrevOverlapping = (String[]) (OverlappingSpectraA[PrevA]);
						String[] NextOverlapping = (String[]) (OverlappingSpectraA[APos[0]]);
						String[] TempA = AntibodyUtils.GetIntersection(
								PrevOverlapping, NextOverlapping);
						String[] TempB = (String[]) (OverlappingSpectraB[Pos[1]]);
						Ret.OverlappingSpectra[i] = AntibodyUtils.GetUnion(
								TempA, TempB);
						if (LocalDebug) {
							System.out.println("Foudn prev A: " + PrevA
									+ " and next A: " + APos[0]);
							System.out
									.println("Prev A: "
											+ AntibodyUtils
													.StringArrayToString(PrevOverlapping));
							System.out
									.println("Next A: "
											+ AntibodyUtils
													.StringArrayToString(NextOverlapping));
							System.out.println("Intersect: "
									+ AntibodyUtils.StringArrayToString(TempA));
							System.out.println("From B: "
									+ AntibodyUtils.StringArrayToString(TempB));
							System.out
									.println("Combined: "
											+ AntibodyUtils
													.StringArrayToString((String[]) (Ret.OverlappingSpectra[i])));
						}
					}
				}
				PrevB = Pos[1];
			}
			if (LocalDebug)
				Utils.WaitForEnter();

		}
		return Ret;
	}

	public Object[] GetSupportingSpectra(int[] ThesePRMs) {
		Object[] Ret = new Object[ThesePRMs.length];
		for (int i = 0; i < ThesePRMs.length; ++i) {
			for (int j = 0; j < this.ScaledPeaks.length; ++j) {
				if (ThesePRMs[i] == this.ScaledPeaks[j]) {
					Ret[i] = this.SupportingSpectra[j];
					break;
				}
			}
		}
		return Ret;
	}

	public Object[] GetOverlappingSpectra(int[] ThesePRMs) {
		Object[] Ret = new Object[ThesePRMs.length];
		for (int i = 0; i < ThesePRMs.length; ++i) {
			for (int j = 0; j < this.ScaledPeaks.length; ++j) {
				if (ThesePRMs[i] == this.ScaledPeaks[j]) {
					Ret[i] = this.OverlappingSpectra[j];
					break;
				}
			}
		}
		return Ret;
	}
}
