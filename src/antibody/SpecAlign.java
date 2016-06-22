package antibody;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

import trieUtils.TrieDB;
import basicUtils.Utils;
import errorUtils.ErrorThrower;
//import basicUtils.InspectAnnotation;

/**
 * Gathers spectra for use in progressive multiple alignment. 
 */

/**
 * The entry point for this method is through Run. The result of Run is a list
 * of MSAlignmentType objects.
 * 
 * @author natalie
 * @version 12.16.2008
 */
public class SpecAlign {

	public final String UsageInfo = "Test";

	// public static final double AveragePRMScore = 0.0;

	private String Sequence;
	private boolean OracleFlag;
	private boolean WriteFlag;
	private Integer ScanNum;
	private Seed Seed;

	private PRMSpectrum[] PRMScans;

	private String PRMScoreFileName;
	private String PeptideOracleFileName;
	private String OutputFileName;

	private int[] SequencePRM;
	private SpectrumNode[] AllSpectra;

	private boolean Debug = false;
	private boolean TestMode = false;
	private TrieDB DecoyDB;

	private String[] Competitors = null;

	public SpecAlign(Seed SeedSequence, SpectrumNode[] Spectra,
			boolean TestMode, TrieDB DecoyDB) {
		this.Sequence = SeedSequence.GetSeedSequence();
		this.TestMode = TestMode;
		this.DecoyDB = DecoyDB;
		this.Seed = SeedSequence;

		// this.Sequence = this.Sequence.substring(Math.max(0,
		// this.Sequence.length() - Utils.MAX_ALIGNMENT_SEQ));

		this.OracleFlag = false;
		this.WriteFlag = false;
		this.ScanNum = new Integer(-1);
		this.AllSpectra = Spectra;

		// this.SequencePRM = SeedSequence.GetSeedPRMs();
		if (Debug) {
			System.out.println("Attempting Alignment to Seed: ");
			SeedSequence.DebugPrint();
		}
	}

	public MSAlignmentType[] RunAlignmentRight(double Tolerance) {
		return RunAlignmentRight(null, Tolerance);
	}

	public MSAlignmentType[] RunAlignmentRight(Template[] Templates,
			double Tolerance) {
		return this.RunAlignmentRight(Templates, null, Tolerance);
	}

	public MSAlignmentType[] RunAlignmentRight(Template[] Templates,
			String[] OtherCompetitors, double Tolerance) {
		if (!TestMode) {
			/*
			 * String[] Els = AntibodyUtils.SplitStringtoAA(this.Sequence);
			 * this.Sequence = ""; for(int i = Math.max(0,Els.length-10); i <
			 * Els.length; ++i) this.Sequence += Els[i]; this.SequencePRM =
			 * AntibodyUtils.GetSeqPRMScaled(this.Sequence);
			 */
			// We only want to use a small fraction of the sequence. So let's
			// lookat the first 10
			// 'used' peaks and take all peaks in that range.

			this.SequencePRM = this.Seed.GetSeedPRMsRight(10);//

		}
		boolean LocalDebug = false;

		// Determining if competitors are needed
		if (Templates != null) {
			ArrayList TempSeqs = new ArrayList();
			System.out.println("* Step 1: Look for competitors");
			for (int i = 0; i < Templates.length; ++i) {
				Template Ptr = Templates[i];
				while (Ptr != null) {
					if (Ptr.Sequence == null) {
						Ptr = Ptr.NextInList;
						continue;
					}
					// System.out.println("Looking at template " +
					// Templates[i].TemplateName);
					Seed Head = Ptr.anchors;
					while (Head != null) {
						if (Head.ClassNum == this.Seed.ClassNum
								&& Head.SeedStart == this.Seed.SeedStart
								&& Head.SeedEnd == this.Seed.SeedEnd) {
							Head = Head.NextSeed;
							continue;
						}
						String AnchorSeq = Head.AnchoredSeedSequence;
						// System.out.println("Comparing to anchor " +
						// AnchorSeq);
						for (int j = 0; j < AnchorSeq.length() - 10; ++j) {
							String TestSeq = AnchorSeq.substring(j, j + 10);
							PRMSpectrum TestSpec = new PRMSpectrum(TestSeq);
							MSAlignmentType TestResult = this.alignSimpleRight(
									TestSpec, null, this.SequencePRM,
									this.Sequence, this.Seed, Tolerance);
							// if(TestResult != null)
							// System.out.println(this.Sequence + " to " +
							// TestSeq + " = " + TestResult.GetScore());
							if (TestResult != null
									&& TestResult.GetScore() >= AntibodyUtils.MSALIGN_SCORE_CUTOFF) {
								System.out.println("Competitor: " + TestSeq);
								TempSeqs.add(TestSeq);
							}

						}
						Head = Head.NextSeed;
					}
					Ptr = Ptr.NextInList;
				}
			}
			// Utils.WaitForEnter();
			if (OtherCompetitors != null) {
				// System.out.println("Have added competitors to consider!!!");
				for (int k = 0; k < OtherCompetitors.length; ++k) {
					String TestSeq = OtherCompetitors[k];
					PRMSpectrum TestSpec = new PRMSpectrum(TestSeq);
					MSAlignmentType TestResult = this.alignSimpleRight(
							TestSpec, null, this.SequencePRM, this.Sequence,
							this.Seed, Tolerance);
					if (TestResult != null
							&& TestResult.GetScore() >= AntibodyUtils.MSALIGN_SCORE_CUTOFF) {
						System.out.println("Competitor: " + TestSeq);
						TempSeqs.add(TestSeq);
					}
				}
			}
			if (TempSeqs.size() > 0)
				this.Competitors = AntibodyUtils
						.ConvertStringArrayList(TempSeqs);
		} else if (OtherCompetitors != null) {
			// System.out.println("Have added competitors to consider!!!");
			ArrayList TempSeqs = new ArrayList();
			for (int k = 0; k < OtherCompetitors.length; ++k) {
				String TestSeq = OtherCompetitors[k];
				PRMSpectrum TestSpec = new PRMSpectrum(TestSeq);
				MSAlignmentType TestResult = this.alignSimpleRight(TestSpec,
						null, this.SequencePRM, this.Sequence, this.Seed,
						Tolerance);
				if (TestResult != null
						&& TestResult.GetScore() >= AntibodyUtils.MSALIGN_SCORE_CUTOFF) {
					System.out.println("Competitor: " + TestSeq);
					TempSeqs.add(TestSeq);
				}
			}
			if (TempSeqs.size() > 0)
				this.Competitors = AntibodyUtils
						.ConvertStringArrayList(TempSeqs);
		}

		// System.out.println("Spec aligning to " + this.Sequence);
		// System.out.println("PRMS: "
		// + AntibodyUtils.IntArrayToString(this.SequencePRM));
		// Utils.WaitForEnter();

		ArrayList Alignments = new ArrayList();
		int TotalSpectra = 0;
		int GoodCount = 0;

		for (int i = 0; i < this.AllSpectra.length; ++i) {

			MSAlignmentType PreviousAlignment = null;// Keeps track of the
			// previous MSAlignment
			// for this node, in
			// case we need to
			// revert
			HMMAlignment prevHMMAlignment = null;

			// The current spectrum to align to the seed (and it's competitors)
			SpectrumNode CurrSpectrum = this.AllSpectra[i];

			PreviousAlignment = this.AllSpectra[i].GetMSAnnotation();
			prevHMMAlignment = this.AllSpectra[i].GetHMMAnnotation();

			// if((CurrSpectrum.GetScanNumber() == 46 &&
			// CurrSpectrum.GetSourceFileName().indexOf("BSA_Other_0_0_550") >=
			// 0 && this.Sequence.indexOf("CADDR") >= 0) ||
			// (CurrSpectrum.GetScanNumber() == 16 &&
			// CurrSpectrum.GetSourceFileName().indexOf("BSA_Other_0_0_850") >=
			// 0 && this.Sequence.indexOf("CADDR") >= 0) ||
			// (CurrSpectrum.GetScanNumber() == 73 &&
			// CurrSpectrum.GetSourceFileName().indexOf("BSA_Other_0_0_825") >=
			// 0 && this.Sequence.indexOf("CADDR") >= 0))
			// LocalDebug = true;
			// else
			// LocalDebug = false;

			if (Debug || LocalDebug) {
				System.out.println("Considering "
						+ CurrSpectrum.GetPRMSpectrum().FileName + " : "
						+ CurrSpectrum.GetPRMSpectrum().SpecIndex);
			}
			// Skip this spectrum if it has an inspect annotation
			if (CurrSpectrum.getModdedDBAnnotation() != null) {

				if (Debug || LocalDebug) {
					System.out
							.println("Skipping this Spectrum because it has an InspectAnnotation:");
					CurrSpectrum.DebugPrint();
					Utils.WaitForEnter();
				}
				continue;
			}

			// Skip this spectrum if it has two few PRM peaks
			if (CurrSpectrum.GetPRMSpectrum().PRMs.length < AntibodyUtils.MIN_PRM_PEAKS)
				continue;

			// It's ok if it has a previous MSAlignment, as long as it is from
			// the same seed and alignment side (right or left extension)

			// For the purpose of making this work with ExecModules I had to
			// take this part out. Previously, it was possible
			// For a spectrum to be successfully recruited and aligned to the
			// HMM, and then in a subsequent round it would be recruited again.
			// The second recruitment would wipe the previous HMM alignment.
			if (CurrSpectrum.GetMSAnnotation() != null) {
				/*
				 * if (Debug || LocalDebug) { System.out
				 * .println("This Spectrum already has an MSAnnotation");
				 * CurrSpectrum.GetMSAnnotation().Print(); Utils.WaitForEnter();
				 * } if (CurrSpectrum.GetMSAnnotation().GetAnchorSeed().Equals(
				 * this.Seed) && CurrSpectrum.GetHMMAnnotation() != null &&
				 * CurrSpectrum.GetHMMAnnotation() .GetAlignmentType() ==
				 * HMMAlignment.HMMAlignmentType.RIGHT_ALIGN) { if (LocalDebug)
				 * {
				 * 
				 * System.out.println(
				 * "But it has already been used for this seed, so let's see if we can use it again!!"
				 * ); CurrSpectrum.GetMSAnnotation().Print();
				 * Utils.WaitForEnter(); } } else continue;
				 */
				continue;
			}

			// Find the best alignment of the spectrum to the seed.
			MSAlignmentType NewAlignment = this.AlignSimple(CurrSpectrum,
					this.SequencePRM, this.Sequence, this.Seed, Tolerance);

			// Check that the NewAlignment is at least as good as the alignemnt
			// of the spectrum to the competitor sequences
			if (NewAlignment != null && this.Competitors != null) {
				for (int c = 0; c < this.Competitors.length; ++c) {

					MSAlignmentType TestAlignment = this.AlignSimple(
							CurrSpectrum,
							AntibodyUtils.GetSeqPRMScaled(this.Competitors[c]),
							this.Competitors[c], null, Tolerance);

					// The alignment to the competitor sequence is better!, so
					// we toss the new alignment
					if (TestAlignment != null
							&& TestAlignment.GetScore() > NewAlignment
									.GetScore()) {
						if (LocalDebug) {
							System.out.println("ELIMINATING MSALIGNMENT for "
									+ CurrSpectrum.GetSourceFileName() + ":"
									+ CurrSpectrum.GetSpecIndex());
							System.out.println("Alignment to " + this.Sequence
									+ " = " + NewAlignment.GetScore());
							System.out.println("Alignment to "
									+ this.Competitors[c] + " = "
									+ TestAlignment.GetScore());
							Utils.WaitForEnter();
						}
						// We created an alignment but we need to clear it out
						CurrSpectrum.SetMSAlignment(PreviousAlignment);
						CurrSpectrum.SetHMMAnnotation(prevHMMAlignment);
						TestAlignment.SetSpectrumNode(null);
						NewAlignment = null;
						break;
					}
				}
				if (NewAlignment != null) {
					CurrSpectrum.SetMSAlignment(NewAlignment);
					CurrSpectrum.SetHMMAnnotation(null);
				}
			}
			TotalSpectra += 1;

			if (NewAlignment != null
					&& (NewAlignment.ExtendsSeed(Tolerance) || TestMode)) {
				if (Debug || LocalDebug) {
					System.out
							.println("Keeping this Spectrum because it has a good alignment");
					CurrSpectrum.DebugPrint();
					NewAlignment.Print();
					Utils.WaitForEnter();
				}
				GoodCount += 1;
				Alignments.add(NewAlignment);

				// NewAlignment.SetSpectrumNode(this.AllSpectra[i]);
				if (TestMode && GoodCount % 100 == 0)
					System.out.println("Good Alignments: " + GoodCount);
			} else {
				if (Debug || LocalDebug) {
					System.out
							.println("Skipping this Spectrum because it has a weak alignment");
					CurrSpectrum.DebugPrint();
					Utils.WaitForEnter();
				}
				if (NewAlignment != null) {
					NewAlignment.SetSpectrumNode(null);

					if (Debug || LocalDebug)

						System.out.println("Doesn't Extend Seed!!" + " for "
								+ this.Sequence);
				}
				this.AllSpectra[i].SetMSAlignment(PreviousAlignment);
				this.AllSpectra[i].SetHMMAnnotation(prevHMMAlignment);
				// CurrSpectrum.DebugPrint();
				// if(Debug || LocalDebug)
				// Utils.WaitForEnter();
			}

		}

		System.out.println("* Step 2: Recruit overlapping spectra");
		System.out.println("Total spectra considered: " + TotalSpectra);
		System.out.println("Total kept: " + Alignments.size());
		MSAlignmentType[] Ret = new MSAlignmentType[Alignments.size()];
		for (int i = 0; i < Alignments.size(); ++i) {
			Ret[i] = ((MSAlignmentType) (Alignments.get(i)));
			Ret[i].GetSpectrumNode().SetMSAlignment(Ret[i]);
			Ret[i].GetSpectrumNode().SetHMMAnnotation(null);
		}
		return Ret;

	}

	public MSAlignmentType[] RunAlignmentLeft(double Tolerance) {
		return RunAlignmentLeft(null, null, Tolerance);
	}

	public MSAlignmentType[] RunAlignmentLeft(Hashtable Oracle, double Tolerance) {
		return this.RunAlignmentLeft(Oracle, null, Tolerance);
	}

	public MSAlignmentType[] RunAlignmentLeft(Hashtable Oracle,
			Template[] Templates, double Tolerance) {

		if (!TestMode) {
			/*
			 * String[] Els = AntibodyUtils.SplitStringtoAA(this.Sequence);
			 * this.Sequence = ""; for(int i = 0; i < Math.min(10,Els.length);
			 * ++i) this.Sequence += Els[i]; this.SequencePRM =
			 * AntibodyUtils.GetSeqPRMScaled(this.Sequence);
			 */
			this.SequencePRM = this.Seed.GetSeedPRMsLeft(10);
		}
		boolean LocalDebug = false;
		// if(this.Sequence.compareTo("LNNFYPKDI") == 0)
		// LocalDebug = true;

		if (LocalDebug) {
			System.out.println("We've got our seed: " + this.Sequence);
			Utils.WaitForEnter();
		}
		// Determining if competitors are needed
		if (Templates != null) {
			ArrayList<String> TempSeqs = new ArrayList<String>();
			System.out.println("* Step 1: Look for competitors");
			for (int i = 0; i < Templates.length; ++i) {
				Template Ptr = Templates[i];
				while (Ptr != null) {
					if (Ptr.Sequence == null) {
						Ptr = Ptr.NextInList;
						continue;
					}
					// System.out.println("Looking at template " +
					// Ptr.TemplateName);
					Seed Head = Ptr.anchors;
					while (Head != null) {
						if (Head.ClassNum == this.Seed.ClassNum
								&& Head.SeedStart == this.Seed.SeedStart
								&& Head.SeedEnd == this.Seed.SeedEnd) {
							Head = Head.NextSeed;
							continue;
						}
						String AnchorSeq = Head.AnchoredSeedSequence;
						// System.out.println("Comparing to anchor " +
						// AnchorSeq);
						for (int j = AnchorSeq.length(); j >= 10; --j) {
							String TestSeq = AnchorSeq.substring(j - 10, j);
							PRMSpectrum TestSpec = new PRMSpectrum(TestSeq);
							MSAlignmentType TestResult = SpecAlign
									.alignSimpleLeft(TestSpec, null,
											this.SequencePRM, this.Sequence,
											this.Seed, Tolerance);
							// if(TestResult != null)
							// System.out.println(this.Sequence + " to " +
							// TestSeq + " = " + TestResult.GetScore());
							if (TestResult != null
									&& TestResult.GetScore() >= AntibodyUtils.MSALIGN_SCORE_CUTOFF) {
								System.out.println("Competitor: " + TestSeq);
								TempSeqs.add(TestSeq);
							}

						}
						Head = Head.NextSeed;
					}
					Ptr = Ptr.NextInList;
				}
			}
			// Utils.WaitForEnter();
			if (TempSeqs.size() > 0)
				this.Competitors = AntibodyUtils
						.ConvertStringArrayList(TempSeqs);
		}
		// System.out.println("Done considering competitors...");
		ArrayList Alignments = new ArrayList();
		int TotalSpectra = 0;
		int GoodCount = 0;

		for (int i = 0; i < this.AllSpectra.length; ++i) {
			SpectrumNode CurrSpectrum = this.AllSpectra[i];
			// LocalDebug = false;
			MSAlignmentType PreviousAlignment = CurrSpectrum.GetMSAnnotation();
			HMMAlignment prevHMMAlignment = CurrSpectrum.GetHMMAnnotation();
			/*
			 * if((CurrSpectrum.GetScanNumber() == 24 &&
			 * CurrSpectrum.GetSourceFileName().indexOf("BSA_Other_0_0_1225") >=
			 * 0) || (CurrSpectrum.GetScanNumber() == 7 &&
			 * CurrSpectrum.GetSourceFileName().indexOf("BSA_Other_0_0_1250") >=
			 * 0) || (CurrSpectrum.GetScanNumber() == 9 &&
			 * CurrSpectrum.GetSourceFileName().indexOf("BSA_Trypsin_0_0_1250")
			 * >= 0) || (CurrSpectrum.GetScanNumber() == 1 &&
			 * CurrSpectrum.GetSourceFileName().indexOf("BSA_Trypsin_0_0_1300")
			 * >= 0)) LocalDebug = true; else LocalDebug = false;
			 * 
			 * if(i == 1110 || i == 1663 || i == 1740 || i == 2296) LocalDebug =
			 * true;
			 */
			// LocalDebug = true;
			/*
			 * if((CurrSpectrum.GetScanNumber() == 66 &&
			 * CurrSpectrum.GetSourceFileName().indexOf("BSA_Other_0_0_500") >=
			 * 0 && this.Sequence.indexOf("LLKH") >= 0) ||
			 * (CurrSpectrum.GetScanNumber() == 79 &&
			 * CurrSpectrum.GetSourceFileName().indexOf("BSA_Other_0_0_550") >=
			 * 0 && this.Sequence.indexOf("LLKH") >= 0) ||
			 * (CurrSpectrum.GetScanNumber() == 31 &&
			 * CurrSpectrum.GetSourceFileName().indexOf("BSA_Trypsin_0_0_550")
			 * >= 0 && this.Sequence.indexOf("LLKH") >= 0)) LocalDebug = true;
			 * else LocalDebug = false;
			 */

			if (Debug || LocalDebug) {
				System.out.println("Considering "
						+ CurrSpectrum.GetPRMSpectrum().FileName + " : "
						+ CurrSpectrum.GetPRMSpectrum().SpecIndex);

			}
			// Skip this spectrum if it has an inspect annotation
			if (CurrSpectrum.getModdedDBAnnotation() != null) {

				if (Debug || LocalDebug) {
					System.out
							.println("Skipping this Spectrum because it has an InspectAnnotation:");
					CurrSpectrum.DebugPrint();
					// Utils.WaitForEnter();
				}
				continue;
			}

			if (CurrSpectrum.GetPRMSpectrum().PRMs.length < AntibodyUtils.MIN_PRM_PEAKS)
				continue;

			if (CurrSpectrum.GetMSAnnotation() != null) {
				continue;
				/*
				 * if (Debug || LocalDebug) { System.out
				 * .println("This SPectrumalready has an MSAnnotation");
				 * CurrSpectrum.GetMSAnnotation().Print(); //
				 * Utils.WaitForEnter(); } if
				 * (CurrSpectrum.GetMSAnnotation().GetAnchorSeed().Equals(
				 * this.Seed) && CurrSpectrum.GetHMMAnnotation() != null &&
				 * CurrSpectrum.GetHMMAnnotation().GetAlignmentType() ==
				 * HMMAlignment.HMMAlignmentType.LEFT_ALIGN) { if (LocalDebug) {
				 * System.out
				 * .println("But it's from the same seed, so we try again...");
				 * 
				 * CurrSpectrum.GetMSAnnotation().Print(); Utils.WaitForEnter();
				 * } } else if
				 * (this.Seed.Succeeds(CurrSpectrum.GetMSAnnotation()
				 * .GetAnchorSeed(), Templates) &&
				 * CurrSpectrum.GetHMMAnnotation() != null &&
				 * CurrSpectrum.GetHMMAnnotation() .GetAlignmentType() ==
				 * HMMAlignment.HMMAlignmentType.RIGHT_ALIGN) { if (LocalDebug)
				 * { System.out
				 * .println("But it's from the previous seed, so we try again..."
				 * ); CurrSpectrum.GetMSAnnotation().Print();
				 * 
				 * Utils.WaitForEnter(); } } else continue;
				 */

			}
			MSAlignmentType NewAlignment = SpecAlign.AlignSimpleLeft(
					CurrSpectrum, this.SequencePRM, this.Sequence, this.Seed,
					Tolerance);
			if (NewAlignment != null && this.Competitors != null) {
				for (int c = 0; c < this.Competitors.length; ++c) {
					MSAlignmentType TestAlignment = SpecAlign.AlignSimple(
							CurrSpectrum,
							AntibodyUtils.GetSeqPRMScaled(this.Competitors[c]),
							this.Competitors[c], this.Seed, Tolerance);
					if (TestAlignment != null
							&& TestAlignment.GetScore() > NewAlignment
									.GetScore()) {
						if (LocalDebug) {
							System.out.println("ELIMINATING MSALIGNMENT for "
									+ CurrSpectrum.GetSourceFileName() + ":"
									+ CurrSpectrum.GetSpecIndex());
							System.out.println("Alignment to " + this.Sequence
									+ " = " + NewAlignment.GetScore());
							System.out.println("Alignment to "
									+ this.Competitors[c] + " = "
									+ TestAlignment.GetScore());
							Utils.WaitForEnter();
						}
						// We created an alignment but we need to clear it out
						CurrSpectrum.SetMSAlignment(PreviousAlignment);
						CurrSpectrum.SetHMMAnnotation(prevHMMAlignment);
						TestAlignment.SetSpectrumNode(null);
						NewAlignment = null;
						break;
					}
				}
				if (NewAlignment != null) {
					CurrSpectrum.SetMSAlignment(NewAlignment);
					CurrSpectrum.SetHMMAnnotation(null);
				}
				CurrSpectrum.SetMSAlignment(NewAlignment);

			}
			TotalSpectra += 1;
			// if (TotalSpectra % 1000 == 0)
			// System.out.println("//: " + TotalSpectra
			// + ", Good Spectra: " + GoodCount);
			if (NewAlignment != null
					&& (NewAlignment.ExtendsSeedLeft(Tolerance) || TestMode)) {
				if (Debug || LocalDebug) {
					System.out
							.println("Keeping this Spectrum because it has a good alignment");
					CurrSpectrum.DebugPrint();
					// NewAlignment.Print();
					System.out.println(NewAlignment.toStringDebug(Oracle));
					Utils.WaitForEnter();
				}
				GoodCount += 1;
				Alignments.add(NewAlignment);

				// NewAlignment.SetSpectrumNode(this.AllSpectra[i]);
				if (TestMode && GoodCount % 100 == 0)
					System.out.println("Good Alignments: " + GoodCount);
			} else {
				if (Debug || LocalDebug)
					System.out
							.println("Skipping this Spectrum because it has a weak alignment");
				if (NewAlignment != null) {
					NewAlignment.SetSpectrumNode(null);

					if (Debug) {
						System.out.println("Doesn't Extend Seed!!" + " for "
								+ this.Sequence);
						System.out.println(NewAlignment.toStringDebug(Oracle));
						// Utils.WaitForEnter();
					}
				}
				this.AllSpectra[i].SetMSAlignment(PreviousAlignment);
				this.AllSpectra[i].SetHMMAnnotation(prevHMMAlignment);
				// CurrSpectrum.DebugPrint();
				// if(Debug || LocalDebug)
				// Utils.WaitForEnter();
				// Utils.WaitForEnter();
			}

		}

		// if(Debug || LocalDebug)
		// {
		System.out.println("* Step 2: Recruit overlapping spectra");
		System.out.println("Total spectra considered: " + TotalSpectra);
		System.out.println("Total kept: " + Alignments.size());
		// }
		MSAlignmentType[] Ret = new MSAlignmentType[Alignments.size()];
		for (int i = 0; i < Alignments.size(); ++i) {
			Ret[i] = ((MSAlignmentType) (Alignments.get(i)));
			Ret[i].GetSpectrumNode().SetMSAlignment(Ret[i]);
			Ret[i].GetSpectrumNode().SetHMMAnnotation(null);
		}
		return Ret;

	}

	public SpecAlign(String[] args) {
		// System.out.println("New SpecAlign...");
		// Set default values
		Sequence = AntibodyUtils.aBTLA_VRegion;
		OracleFlag = false;
		WriteFlag = false;
		ScanNum = new Integer(-1);

		SequencePRM = AntibodyUtils.GetSeqPRMScaled(Sequence);

		ParseCommandLine(args);
	}

	public static MSAlignmentType AlignSimple(SpectrumNode ContainingNode,
			int[] SequencePRMs, String Sequence, Seed CurrSeed, double Tolerance) {
		return alignSimpleRight(ContainingNode.GetPRMSpectrum(),
				ContainingNode, SequencePRMs, Sequence, CurrSeed, Tolerance);
	}

	/**
	 * The main alignment routine for aligning a spectrum to a sequence. Assumes
	 * the spectrum extends the sequence to the right. TODO: Implelment!
	 * 
	 * @param Spectrum
	 * @param ContainingNode
	 * @param CurrSequencePRMs
	 * @param Sequence
	 * @param CurrSeed
	 * @return
	 */
	public static MSAlignmentType alignSimpleRight(PRMSpectrum spectrum,
			SpectrumNode ContainingNode, int[] currSequencePRMs,
			String Sequence, Seed CurrSeed, double tol) {

		boolean LocalDebug = false;

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
		for (int i = AntibodyUtils.MIN_OVERLAPPING_PEAKS - 1; i < currSequencePRMs.length; ++i) {
			for (int j = AntibodyUtils.MIN_OVERLAPPING_PEAKS - 1; j < spectrum.PRMs.length; ++j) {
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
		// String matchStr = "";
		// String prmStr = "";
		int length = 0;

		/*
		 * This is the only difference between AlignSimpleLeft and
		 * AlignSimpleRight For AlignSimpleRight, the alignment needs to end
		 * close to the end of the sequence.
		 */
		if (bestIndx[0] < currSequencePRMs.length - 3) {
			if (LocalDebug)
				System.out.println("Alignment doesn't end near the end of Seq "
						+ bestIndx[0] + " < " + (currSequencePRMs.length - 3));
			return null;
		}

		ArrayList<Integer> SpecMasses = new ArrayList<Integer>();
		ArrayList<Integer> SeqMasses = new ArrayList<Integer>();

		while (bestIndx[0] != -1) {
			// matchStr = "(" + bestIndx[0] + "," + bestIndx[1] + ") " +
			// matchStr;
			// prmStr = "(" + currSequencePRMs[bestIndx[0]] + "," +
			// spectrum.PRMs[bestIndx[1]] + ") " + prmStr;
			length += 1;
			SpecMasses.add(0, new Integer(spectrum.PRMs[bestIndx[1]]));
			SeqMasses.add(0, new Integer(currSequencePRMs[bestIndx[0]]));
			int prevVal = prev[bestIndx[0]][bestIndx[1]];
			if (prevVal == -1)
				break;

			bestIndx[0] = prevVal / spectrum.PRMs.length;
			bestIndx[1] = prevVal % spectrum.PRMs.length;
			if (LocalDebug)
				System.out.println("prev: " + bestIndx[0] + "," + bestIndx[1]);
		}

		if (length < AntibodyUtils.MIN_OVERLAPPING_PEAKS) {
			// System.out.println("Alignment is short: " + length
			// + " < " + (AntibodyUtils.MIN_OVERLAPPING_PEAKS));
			return null;
		}

		return new MSAlignmentType(spectrum, Sequence, currSequencePRMs,
				bestScore, AntibodyUtils.ConvertIntegerArrayList(SpecMasses),
				AntibodyUtils.ConvertIntegerArrayList(SeqMasses),
				ContainingNode, CurrSeed);

	}

	public static MSAlignmentType AlignSimpleRight_old(PRMSpectrum Spectrum,
			SpectrumNode ContainingNode, int[] CurrSequencePRMs,
			String Sequence, Seed CurrSeed, double Tolerance) {
		boolean LocalDebug = false;

		/*
		 * if(ContainingNode != null && ((ContainingNode.GetScanNumber() == 46
		 * && ContainingNode.GetSourceFileName().indexOf("BSA_Other_0_0_550") >=
		 * 0 && Sequence.indexOf("CADDR") >= 0) ||
		 * (ContainingNode.GetScanNumber() == 16 &&
		 * ContainingNode.GetSourceFileName().indexOf("BSA_Other_0_0_850") >= 0
		 * && Sequence.indexOf("CADDR") >= 0) || (ContainingNode.GetScanNumber()
		 * == 73 &&
		 * ContainingNode.GetSourceFileName().indexOf("BSA_Other_0_0_825") >= 0
		 * && Sequence.indexOf("CADDR") >= 0))) LocalDebug = true; else
		 * LocalDebug = false;
		 */
		if (LocalDebug) {
			System.out
					.println(AntibodyUtils.IntegerListToString(Spectrum.PRMs));
			System.out.println(Utils.DoubleArrayToString(Spectrum.PRMScores));
			Utils.WaitForEnter();
		}

		// Inititalize the alignment matrix
		double[][] H = new double[CurrSequencePRMs.length][Spectrum.PRMs.length];
		int[][] PrevH = new int[CurrSequencePRMs.length][Spectrum.PRMs.length];

		// MassShift represents how much we have to add to the SpectrumPRM to
		// get to the Sequence PRM
		int[][] MassShift = new int[CurrSequencePRMs.length][Spectrum.PRMs.length];

		H[0][0] = Spectrum.PRMScores[0];
		PrevH[0][0] = -1;

		for (int i = 0; i < CurrSequencePRMs.length; ++i) {
			H[i][0] = H[0][0];
			PrevH[i][0] = -1;
			MassShift[i][0] = CurrSequencePRMs[i];
		}
		for (int i = 0; i < Spectrum.PRMs.length; ++i) {
			H[0][i] = Spectrum.PRMScores[i];

			PrevH[0][i] = -1;
			MassShift[0][i] = -1 * Spectrum.PRMs[i];
		}

		double BestValue = Double.MIN_VALUE;
		int BestCell = -1;
		int SequenceMass = CurrSequencePRMs[CurrSequencePRMs.length - 1];

		// Fill in the table
		for (int SeqIndex = 1; SeqIndex < CurrSequencePRMs.length; ++SeqIndex) {
			int CurrSeqMass = CurrSequencePRMs[SeqIndex];
			for (int SpecIndex = 1; SpecIndex < Spectrum.PRMs.length; ++SpecIndex) {
				int CurrSpecMass = Spectrum.PRMs[SpecIndex];

				// double Value = H[0][0] + Spectrum.PRMScores[SpecIndex];
				// //Value of this alignment if we align with 0,0
				// double Value = Spectrum.PRMScores[SpecIndex]; //Value of just
				// aligning this spectrum
				double Value = Spectrum.PRMScores[SpecIndex];
				double MaxValue = Value;

				int CurrAlignmentShift = CurrSeqMass - CurrSpecMass;
				int MaxAlignmentShift = CurrAlignmentShift;

				int PrevCell = -1; // Index of previous match
				int MaxPrevCell = PrevCell;

				// Look back no more more than 1 missed AA
				for (int PrevSeqIndex = SeqIndex - 1; PrevSeqIndex > Math.max(
						SeqIndex - 4, -1); PrevSeqIndex--) {
					int PrevSeqMass = CurrSequencePRMs[PrevSeqIndex];

					for (int PrevSpecIndex = SpecIndex - 1; PrevSpecIndex >= 0; PrevSpecIndex--) {
						int PrevSpecMass = Spectrum.PRMs[PrevSpecIndex];

						// Check that the delta in Spectrum mass is about the
						// same as teh delta in Sequence masses
						// Check that the CurrentAlignment Shift is about the
						// same as the previous alignment shift
						if (Math.abs((CurrSeqMass - PrevSeqMass)
								- (CurrSpecMass - PrevSpecMass)) < (int) (Tolerance * AntibodyUtils.MASS_SCALE)
								&& Math.abs(CurrAlignmentShift
										- MassShift[PrevSeqIndex][PrevSpecIndex]) < (int) (Tolerance * AntibodyUtils.MASS_SCALE)) {
							if (H[PrevSeqIndex][PrevSpecIndex] > 0)
								Value = H[PrevSeqIndex][PrevSpecIndex]
										+ Spectrum.PRMScores[SpecIndex];
						}
						if (Value > MaxValue) {
							MaxValue = Value;
							MaxPrevCell = PrevSeqIndex * Spectrum.PRMs.length
									+ PrevSpecIndex;
							MaxAlignmentShift = MassShift[PrevSeqIndex][PrevSpecIndex];
						}
					}

				}

				H[SeqIndex][SpecIndex] = MaxValue;
				PrevH[SeqIndex][SpecIndex] = MaxPrevCell;
				MassShift[SeqIndex][SpecIndex] = MaxAlignmentShift;

				if (MaxValue > BestValue) {
					BestValue = MaxValue;
					BestCell = SeqIndex * Spectrum.PRMs.length + SpecIndex;
					if (LocalDebug) {

						System.out.println("New Max at [" + SeqIndex + "]["
								+ SpecIndex + "]=" + BestValue);
						System.out.println("Prev [" + MaxPrevCell
								/ Spectrum.PRMs.length + "," + MaxPrevCell
								% Spectrum.PRMs.length + "]");
						System.out.println("MassShift: " + MaxAlignmentShift);
						Utils.WaitForEnter();
					}
				}

			}

		}
		int Previ = 0;

		int BestSeqIndex = BestCell / Spectrum.PRMs.length;
		int BestSpecIndex = BestCell % Spectrum.PRMs.length;

		if (BestSeqIndex < CurrSequencePRMs.length - 3) {
			if (LocalDebug)
				System.out.println("Alignment doesn't end near the end of Seq "
						+ BestSeqIndex + " < " + (CurrSequencePRMs.length - 3));
			return null;
		}
		ArrayList SeqMasses = new ArrayList();
		ArrayList SpecMasses = new ArrayList();

		int TempSeqIndex = BestSeqIndex;
		int TempSpecIndex = BestSpecIndex;
		int Temp = BestCell;

		while (Temp >= 0) {
			TempSeqIndex = Temp / Spectrum.PRMs.length;
			TempSpecIndex = Temp % Spectrum.PRMs.length;
			SeqMasses.add(0, new Integer(CurrSequencePRMs[TempSeqIndex]));
			SpecMasses.add(0, new Integer(Spectrum.PRMs[TempSpecIndex]));
			Previ = TempSeqIndex;
			Temp = PrevH[TempSeqIndex][TempSpecIndex];

		}

		if (LocalDebug) {
			System.out.println("AlignedPeaks: "
					+ AntibodyUtils.IntegerListToString(AntibodyUtils
							.ConvertIntegerArrayList(SeqMasses)));
			System.out.println("SeqMasses.size() = " + SeqMasses.size());
		}

		// Check for a valid prefix or suffix mass to match
		int[] SpecMassesFinal = AntibodyUtils
				.ConvertIntegerArrayList(SpecMasses);
		int[] SeqMassesFinal = AntibodyUtils.ConvertIntegerArrayList(SeqMasses);

		boolean HasPrefix = SpecAlign.IsValidPrefix(SpecMassesFinal[0],
				SeqMassesFinal[0], CurrSequencePRMs, Tolerance);
		boolean HasSuffix = SpecAlign.IsValidSuffix(
				SpecMassesFinal[SpecMassesFinal.length - 1],
				SeqMassesFinal[SeqMassesFinal.length - 1],
				Spectrum.PrecursorMass, CurrSequencePRMs, Tolerance);

		int SpectrumPrefixMass = SpecAlign.GetPrefixMass(SpecMassesFinal[0],
				SeqMassesFinal[0]);
		int SpectrumSuffixMass = SpecAlign.GetSuffixMass(
				SpecMassesFinal[SpecMassesFinal.length - 1],
				SeqMassesFinal[SeqMassesFinal.length - 1],
				Spectrum.PrecursorMass, SequenceMass);

		if (HasPrefix
				&& SpectrumPrefixMass >= -1 * Tolerance
						* AntibodyUtils.MASS_SCALE) {
			if (LocalDebug)
				System.out.println("HAS PREFIX!");
			BestValue += 15;
		}
		if (HasSuffix
				&& SpectrumSuffixMass <= SequenceMass + Tolerance
						* AntibodyUtils.MASS_SCALE) {
			if (LocalDebug)
				System.out.println("HAS SUFFIX!");
			BestValue += 15;
		}
		if (SpectrumPrefixMass < -1 * Tolerance * AntibodyUtils.MASS_SCALE
				&& SpectrumSuffixMass > SequenceMass + Tolerance
						* AntibodyUtils.MASS_SCALE) {
			if (SeqMasses.size() < AntibodyUtils.MIN_OVERLAPPING_PEAKS + 1) {
				if (LocalDebug)
					System.out
							.println("Spectrum totally overlaps anchor, but alignment is short: "
									+ SeqMasses.size()
									+ " < "
									+ (AntibodyUtils.MIN_OVERLAPPING_PEAKS + 1));
				return null;
			}
		}
		if (SeqMasses.size() < AntibodyUtils.MIN_OVERLAPPING_PEAKS) {
			if (LocalDebug)
				System.out.println("Alignment is short: " + SeqMasses.size()
						+ " < " + (AntibodyUtils.MIN_OVERLAPPING_PEAKS));
			return null;
		}

		return new MSAlignmentType(Spectrum, Sequence, CurrSequencePRMs,
				BestValue, AntibodyUtils.ConvertIntegerArrayList(SpecMasses),
				AntibodyUtils.ConvertIntegerArrayList(SeqMasses),
				ContainingNode, CurrSeed);

	}

	public static MSAlignmentType AlignSimpleLeft(SpectrumNode ContainingNode,
			int[] SequencePRMs, String Sequence, Seed CurrSeed, double Tolerance) {
		return alignSimpleLeft(ContainingNode.GetPRMSpectrum(), ContainingNode,
				SequencePRMs, Sequence, CurrSeed, Tolerance);
	}

	/**
	 * Constructs the alignment table, allowing no modifications and returns the
	 * highest scoring subpath. TODO: Verify that alignment is for the left
	 * correctly
	 * 
	 * @param Spectrum
	 *            the PRMSpectrum to be aligned to the anchor sequence
	 * @return the Alignment structure for this spectrum, or null if no good
	 *         alignment exists
	 */
	public static MSAlignmentType alignSimpleLeft(PRMSpectrum spectrum,
			SpectrumNode ContainingNode, int[] currSequencePRMs,
			String Sequence, Seed CurrSeed, double tol) {

		boolean LocalDebug = false;

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
		for (int i = AntibodyUtils.MIN_OVERLAPPING_PEAKS - 1; i < currSequencePRMs.length; ++i) {
			for (int j = AntibodyUtils.MIN_OVERLAPPING_PEAKS - 1; j < spectrum.PRMs.length; ++j) {
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
		// String matchStr = "";
		// String prmStr = "";
		int length = 0;

		ArrayList<Integer> SpecMasses = new ArrayList<Integer>();
		ArrayList<Integer> SeqMasses = new ArrayList<Integer>();

		while (bestIndx[0] != -1) {
			// matchStr = "(" + bestIndx[0] + "," + bestIndx[1] + ") " +
			// matchStr;
			// prmStr = "(" + currSequencePRMs[bestIndx[0]] + "," +
			// spectrum.PRMs[bestIndx[1]] + ") " + prmStr;
			length += 1;
			SpecMasses.add(0, new Integer(spectrum.PRMs[bestIndx[1]]));
			SeqMasses.add(0, new Integer(currSequencePRMs[bestIndx[0]]));
			int prevVal = prev[bestIndx[0]][bestIndx[1]];
			if (prevVal == -1)
				break;

			bestIndx[0] = prevVal / spectrum.PRMs.length;
			bestIndx[1] = prevVal % spectrum.PRMs.length;
			if (LocalDebug)
				System.out.println("prev: " + bestIndx[0] + "," + bestIndx[1]);
		}

		if (length < AntibodyUtils.MIN_OVERLAPPING_PEAKS) {
			// System.out.println("Alignment is short: " + length
			// + " < " + (AntibodyUtils.MIN_OVERLAPPING_PEAKS));
			return null;
		}

		/**
		 * This is the only difference between AlignSimpleLeft and
		 * AlignSimpleRight In AlignSimpleLeft we want the start of the
		 * alignment to be close to the start of the sequence.
		 */
		if (bestIndx[0] > 2) {
			return null;
		}

		return new MSAlignmentType(spectrum, Sequence, currSequencePRMs,
				bestScore, AntibodyUtils.ConvertIntegerArrayList(SpecMasses),
				AntibodyUtils.ConvertIntegerArrayList(SeqMasses),
				ContainingNode, CurrSeed);

	}

	public static MSAlignmentType AlignSimpleLeft_old(PRMSpectrum Spectrum,
			SpectrumNode ContainingNode, int[] CurrSequencePRMs,
			String Sequence, Seed CurrSeed, double Tolerance) {

		// Inititalize the alignment matrix
		double[][] H = new double[CurrSequencePRMs.length][Spectrum.PRMs.length];
		int[][] PrevH = new int[CurrSequencePRMs.length][Spectrum.PRMs.length];

		double OverallBestValue = Double.MIN_VALUE;
		int OverallBestCell = -1;

		int SequenceMass = CurrSequencePRMs[CurrSequencePRMs.length - 1];

		// MassShift represents how much we have to add to the SpectrumPRM to
		// get to the Sequence PRM
		int[][] MassShift = new int[CurrSequencePRMs.length][Spectrum.PRMs.length];

		H[0][0] = Spectrum.PRMScores[0];
		PrevH[0][0] = -1;

		for (int i = 0; i < CurrSequencePRMs.length; ++i) {
			H[i][0] = H[0][0];
			PrevH[i][0] = -1;
			MassShift[i][0] = CurrSequencePRMs[i];
		}
		for (int i = 0; i < Spectrum.PRMs.length; ++i) {
			H[0][i] = Spectrum.PRMScores[i];
			PrevH[0][i] = -1;
			MassShift[0][i] = -1 * Spectrum.PRMs[i];
			// MassShift[0][i] = Spectrum.PRMs[i];
		}
		boolean LocalDebug = false;
		/*
		 * if(ContainingNode != null && ((ContainingNode.GetScanNumber() == 66
		 * && ContainingNode.GetSourceFileName().indexOf("BSA_Other_0_0_500") >=
		 * 0 && Sequence.indexOf("LLKH") >= 0) ||
		 * (ContainingNode.GetScanNumber() == 79 &&
		 * ContainingNode.GetSourceFileName().indexOf("BSA_Other_0_0_550") >= 0
		 * && Sequence.indexOf("LLKH") >= 0) || (ContainingNode.GetScanNumber()
		 * == 31 &&
		 * ContainingNode.GetSourceFileName().indexOf("BSA_Trypsin_0_0_550") >=
		 * 0 && Sequence.indexOf("LLKH") >= 0))) LocalDebug = true; else
		 * LocalDebug = false;
		 */
		if (LocalDebug) {
			System.out.println(Utils.IntArrayToString(Spectrum.PRMs));
			System.out.println(Utils.DoubleArrayToString(Spectrum.PRMScores));
		}
		// Fill in the table
		for (int SeqIndex = 1; SeqIndex < CurrSequencePRMs.length; ++SeqIndex) {
			int CurrSeqMass = CurrSequencePRMs[SeqIndex];
			for (int SpecIndex = 1; SpecIndex < Spectrum.PRMs.length; ++SpecIndex) {
				int CurrSpecMass = Spectrum.PRMs[SpecIndex];

				double Value = Spectrum.PRMScores[SpecIndex]; // Value of just
				// aligning this
				// spectrum
				double MaxValue = Value;

				int CurrAlignmentShift = CurrSeqMass - CurrSpecMass;
				int MaxAlignmentShift = CurrAlignmentShift;

				int PrevCell = -1; // Index of previous match
				int MaxPrevCell = PrevCell;
				// Look back no more more than 1 missed AA
				for (int PrevSeqIndex = SeqIndex - 1; PrevSeqIndex > Math.max(
						SeqIndex - 3, -1); PrevSeqIndex--) {
					int PrevSeqMass = CurrSequencePRMs[PrevSeqIndex];

					for (int PrevSpecIndex = SpecIndex - 1; PrevSpecIndex >= 0; PrevSpecIndex--) {
						int PrevSpecMass = Spectrum.PRMs[PrevSpecIndex];
						// Check that the delta in Spectrum mass is about the
						// same as teh delta in Sequence masses
						// Check that the CurrentAlignment Shift is about the
						// same as the previous alignment shift
						if (Math.abs((CurrSeqMass - PrevSeqMass)
								- (CurrSpecMass - PrevSpecMass)) < (int) (Tolerance * AntibodyUtils.MASS_SCALE)
								&& Math.abs(CurrAlignmentShift
										- MassShift[PrevSeqIndex][PrevSpecIndex]) < (int) (Tolerance * AntibodyUtils.MASS_SCALE)) {
							Value = H[PrevSeqIndex][PrevSpecIndex]
									+ Spectrum.PRMScores[SpecIndex];
						}

						if (Value > MaxValue) {
							MaxValue = Value;
							MaxPrevCell = PrevSeqIndex * Spectrum.PRMs.length
									+ PrevSpecIndex;
							MaxAlignmentShift = MassShift[PrevSeqIndex][PrevSpecIndex];
						}
					}

				}
				H[SeqIndex][SpecIndex] = MaxValue;
				PrevH[SeqIndex][SpecIndex] = MaxPrevCell;
				MassShift[SeqIndex][SpecIndex] = MaxAlignmentShift;

				if (MaxValue > OverallBestValue) {
					OverallBestValue = MaxValue;
					OverallBestCell = SeqIndex * Spectrum.PRMs.length
							+ SpecIndex;
					if (LocalDebug) {
						System.out.println("New Max Value: " + MaxValue);
						System.out.println("Ends at [" + SeqIndex + ","
								+ SpecIndex + "]");
						System.out.println("Prev [" + MaxPrevCell
								/ Spectrum.PRMs.length + "," + MaxPrevCell
								% Spectrum.PRMs.length + "]");
						Utils.WaitForEnter();
					}
				}
			}

		}

		// Determine best partial alignment

		String PrevStr = "";
		String Str = "";
		int CurrCell = OverallBestCell;

		int Previ = 0;

		int BestSeqIndex = OverallBestCell / Spectrum.PRMs.length;
		int BestSpecIndex = OverallBestCell % Spectrum.PRMs.length;
		if (LocalDebug) {
			System.out.println("Best End: [" + BestSeqIndex + ","
					+ BestSpecIndex + "]");
			System.out.println("Mass Shift: "
					+ MassShift[BestSeqIndex][BestSpecIndex]);
		}
		ArrayList SeqMasses = new ArrayList();
		ArrayList SpecMasses = new ArrayList();

		int TempSeqIndex = BestSeqIndex;
		int TempSpecIndex = BestSpecIndex;
		int Temp = OverallBestCell;
		if (Temp < 0) {
			return null;
		}
		while (Temp >= 0) {
			SeqMasses.add(0, new Integer(CurrSequencePRMs[TempSeqIndex]));
			SpecMasses.add(0, new Integer(Spectrum.PRMs[TempSpecIndex]));
			Previ = TempSeqIndex;
			Temp = PrevH[TempSeqIndex][TempSpecIndex];
			if (Temp == -1 && TempSeqIndex > 2) {
				if (LocalDebug)
					System.out.println("FAIL: start in sequence is "
							+ TempSeqIndex + " > 2");
				return null;
			}
			TempSeqIndex = Temp / Spectrum.PRMs.length;
			TempSpecIndex = Temp % Spectrum.PRMs.length;

		}

		// String[] SequenceSegments = AntibodyUtils.SplitStringtoAA(Sequence);
		// String Prefix = "";
		// for(int i = Previ; i < BestSeqIndex; ++i)
		// {
		// Prefix += SequenceSegments[i];
		// }
		// String Prefix = this.Sequence.substring(Previ,BestSeqIndex);
		int[] SpecMassesFinal = AntibodyUtils
				.ConvertIntegerArrayList(SpecMasses);
		int[] SeqMassesFinal = AntibodyUtils.ConvertIntegerArrayList(SeqMasses);

		boolean HasPrefix = IsValidPrefix(SpecMassesFinal[0],
				SeqMassesFinal[0], CurrSequencePRMs, Tolerance);
		boolean HasSuffix = IsValidSuffix(
				SpecMassesFinal[SpecMassesFinal.length - 1],
				SeqMassesFinal[SeqMassesFinal.length - 1],
				Spectrum.PrecursorMass, CurrSequencePRMs, Tolerance);

		int SpectrumPrefixMass = GetPrefixMass(SpecMassesFinal[0],
				SeqMassesFinal[0]);
		int SpectrumSuffixMass = GetSuffixMass(
				SpecMassesFinal[SpecMassesFinal.length - 1],
				SeqMassesFinal[SeqMassesFinal.length - 1],
				Spectrum.PrecursorMass, SequenceMass);

		if (HasPrefix
				&& SpectrumPrefixMass >= -1 * Tolerance
						* AntibodyUtils.MASS_SCALE) {
			if (LocalDebug)
				System.out.println(" We have a prefix!!!");
			OverallBestValue += 15;
		}
		if (HasSuffix
				&& SpectrumSuffixMass <= SequenceMass + Tolerance
						* AntibodyUtils.MASS_SCALE) {
			if (LocalDebug)
				System.out.println(" We have a suffix!!!");
			OverallBestValue += 15;
		}
		if (LocalDebug)
			System.out.println("Best value: " + OverallBestValue + "!!");
		if (SpectrumPrefixMass < -1 * Tolerance * AntibodyUtils.MASS_SCALE
				&& SpectrumSuffixMass > SequenceMass + Tolerance
						* AntibodyUtils.MASS_SCALE) {
			if (SeqMasses.size() < AntibodyUtils.MIN_OVERLAPPING_PEAKS + 1) {
				if (LocalDebug)
					System.out
							.println("Spectrum totally overlaps the sequence, but we have a weak alignment of length "
									+ SeqMasses.size()
									+ " < "
									+ (AntibodyUtils.MIN_OVERLAPPING_PEAKS + 1));

				return null;
			}
		}
		if (SeqMasses.size() < AntibodyUtils.MIN_OVERLAPPING_PEAKS) {
			if (LocalDebug)
				System.out.println("Short alignment: " + SeqMasses.size()
						+ " < " + AntibodyUtils.MIN_OVERLAPPING_PEAKS);
			return null;
		}

		return new MSAlignmentType(Spectrum, Sequence, CurrSequencePRMs,
				OverallBestValue,
				AntibodyUtils.ConvertIntegerArrayList(SpecMasses),
				AntibodyUtils.ConvertIntegerArrayList(SeqMasses),
				ContainingNode, CurrSeed);

	}

	public static boolean IsValidPrefix(int CurrSpecMass, int CurrSeqMass,
			int[] SequencePRM, double Tolerance) {
		int MassDiff = CurrSeqMass - CurrSpecMass;

		for (int i = 0; i < SequencePRM.length; ++i)
			if (Math.abs(MassDiff - SequencePRM[i]) <= Tolerance
					* AntibodyUtils.MASS_SCALE)
				return true;
		return false;
	}

	/*
	 * 
	 * public boolean IsValidPrefix(int CurrSpecMass, int CurrSeqMass, boolean
	 * Debug) { int MassDiff = CurrSeqMass - CurrSpecMass;
	 * 
	 * if(Debug) System.out.println("Looking for prefix mass at " +CurrSeqMass +
	 * "-" + CurrSpecMass + "=" + MassDiff);
	 * 
	 * for(int i = 0; i < SequencePRM.length; ++i) if(Math.abs(MassDiff -
	 * SequencePRM[i]) <=
	 * AntibodyUtils.LTQFTIonTolerance*AntibodyUtils.MASS_SCALE) return true;
	 * return false; }
	 */

	public static int GetPrefixMass(int CurrSpecMass, int CurrSeqMass) {
		return CurrSeqMass - CurrSpecMass;
	}

	public static boolean IsValidSuffix(int CurrSpecMass, int CurrSeqMass,
			int PrecursorMass, int[] SequencePRM, double Tolerance) {
		int MassDiff = CurrSeqMass + (PrecursorMass - CurrSpecMass);

		for (int i = 0; i < SequencePRM.length; ++i)
			if (Math.abs(MassDiff - SequencePRM[i]) <= Tolerance
					* AntibodyUtils.MASS_SCALE)
				return true;
		return false;
	}

	public static int GetSuffixMass(int CurrSpecMass, int CurrSeqMass,
			int PrecursorMass, int SequenceMass) {
		return CurrSeqMass + (PrecursorMass - CurrSpecMass);

	}

	private void ParseCommandLine(String[] args) {
		if (Debug)
			System.out.println("SpecALign: ParseCommandLine...");
		// Parse command line
		String[] options = { "-s", "-p", "-r", "-o", "-w" };
		boolean[] values = { true, true, true, true, true };
		Hashtable Args = Utils.ParseCommandLine(args, options, values);

		if (!Args.containsKey("-r") || !Args.containsKey("-w")) {
			ErrorThrower.ThrowError(2,
					"Did not specify an MZXML file or PRM file");
		}

		Enumeration Keys = Args.keys();
		while (Keys.hasMoreElements()) {
			String Arg = (String) Keys.nextElement();
			// System.out.println("Arg: " + Arg);
			if (Arg == "-s")
				Sequence = (String) Args.get(Arg);
			else if (Arg == "-p")
				ScanNum = new Integer((String) Args.get(Arg));
			else if (Arg == "-r")
				PRMScoreFileName = (String) Args.get(Arg);
			else if (Arg == "-o") {
				PeptideOracleFileName = (String) Args.get(Arg);
				OracleFlag = true;
			} else if (Arg == "-w") {
				OutputFileName = (String) Args.get(Arg);
				WriteFlag = true;
			}
		}
	}

	/*
	 * private static boolean TestRun(String[] PRMFileNames, String[]
	 * OracleFileNames, String Seed, String OutputFileName, double Tolerance) {
	 * 
	 * // Load the PRMS PRMSpectrum[] AllPRMs; ArrayList TempList = new
	 * ArrayList(); for (int i = 0; i < PRMFileNames.length; ++i) { try {
	 * AllPRMs = PRMSpectrum.LoadAllPRMSpectrum(PRMFileNames[i]); } catch
	 * (Exception E) { System.err.println(E.getMessage()); return false; }
	 * System.out.println("Loaded " + AllPRMs.length + " PRMS from " +
	 * PRMFileNames[i] + "..."); for (int j = 0; j < AllPRMs.length; ++j)
	 * TempList.add(AllPRMs[j]); } SpectrumNode[] AllSpectra = new
	 * SpectrumNode[TempList.size()]; for (int i = 0; i < TempList.size(); ++i)
	 * { AllSpectra[i] = new SpectrumNode((PRMSpectrum) (TempList.get(i))); }
	 * 
	 * // Load the oracle InspectAnnotation[] Oracle; TempList.clear(); for (int
	 * i = 0; i < OracleFileNames.length; ++i) { Oracle = InspectAnnotation
	 * .LoadInspectResultsFile(OracleFileNames[i]); System.out.println("Loaded "
	 * + Oracle.length + " oracle seqs from " + OracleFileNames[i] + "..."); for
	 * (int j = 0; j < Oracle.length; ++j) TempList.add(Oracle[j]); } Oracle =
	 * new InspectAnnotation[TempList.size()]; for (int i = 0; i <
	 * Oracle.length; ++i) { Oracle[i] = (InspectAnnotation) (TempList.get(i));
	 * }
	 * 
	 * // Create the test seed if (Seed == null) Seed =
	 * AntibodyUtils.aBTLA_OracleSequence; ProteinSeed TestSeed = new
	 * ProteinSeed(Seed);
	 * 
	 * // Run the alignment SpecAlign Tester = new SpecAlign(TestSeed,
	 * AllSpectra, true, null); MSAlignmentType[] TestResults =
	 * Tester.RunAlignmentRight(Tolerance);
	 * 
	 * // TestResults = //
	 * MSAlignmentType.GetBestNonRedundantOdds(TestResults.length, //
	 * TestResults); TestResults = MSAlignmentType.GetBestOddsCutOff(
	 * AntibodyUtils.MSALIGN_PVALUE_CUTOFF, TestResults, -1);
	 * WriteTestOutput(TestResults, Oracle, OutputFileName, PRMFileNames,
	 * OracleFileNames, Seed);
	 * 
	 * return true; }
	 * 
	 * public static void WriteTestOutput(MSAlignmentType[] Alignments,
	 * InspectAnnotation[] Oracle, String OutputFileName, String[] PRMFileNames,
	 * String[] OracleFileNames, String SeedSeq) { BufferedWriter Writer = null;
	 * String Line =
	 * "#SpectrumFile\tScanNumber\tMSAlignment\tOracle\tCorrect\tMSA_Score\tMSA_Pos\tOraclePos\n"
	 * ; try { Writer = new BufferedWriter(new FileWriter(OutputFileName)); for
	 * (int i = 0; i < PRMFileNames.length; ++i) { Writer.write("#PRM: " +
	 * PRMFileNames[i]); Writer.newLine(); } for (int i = 0; i <
	 * OracleFileNames.length; ++i) { Writer.write("#Oracle: " +
	 * OracleFileNames[i]); Writer.newLine(); } Writer.write("#AlignmentSeq: " +
	 * SeedSeq); Writer.newLine(); Writer.newLine(); Writer.write(Line);
	 * 
	 * } catch (Exception E) { System.err.println(E.getLocalizedMessage());
	 * return; }
	 * 
	 * for (int i = 0; i < Alignments.length; ++i) { String SpectrumFile =
	 * Alignments[i].GetSpectrum().FileName; String RootName =
	 * Utils.GetFileNameNoExtension(Utils .GetBaseName(SpectrumFile)); int
	 * ScanNumber = Alignments[i].GetSpectrum().ScanNumber;
	 * 
	 * String OracleString = "*"; for (int j = 0; j < Oracle.length; ++j) {
	 * String OracleFile = Utils.GetFileNameNoExtension(Utils
	 * .GetBaseName(Oracle[j].SpectrumFile)); if (RootName.compareTo(OracleFile)
	 * != 0) continue; if (Oracle[j].ScanNumber != ScanNumber) continue;
	 * OracleString = Oracle[j].Annotation; }
	 * 
	 * int MSAStart = -1; int MSAEnd = MSAStart +
	 * Alignments[i].GetAlignedPeaksScaled().length;
	 * 
	 * int OracleStart = AntibodyUtils.aBTLA_OracleSequence
	 * .indexOf(OracleString); int OracleEnd = OracleStart +
	 * OracleString.length(); if (OracleStart == -1) OracleEnd = -1;
	 * 
	 * int Correct = 0; if (Utils.HasOverlap(MSAStart, MSAEnd, OracleStart,
	 * OracleEnd)) Correct = 1; if (OracleStart == -1) Correct = -1;
	 * 
	 * Line = SpectrumFile + "\t" + ScanNumber + "\t" + OracleString + "\t" +
	 * Correct + "\t" + "\t" + MSAStart + "\t" + OracleStart + "\n"; try {
	 * Writer.write(Line); } catch (Exception E) {
	 * System.err.println(E.getLocalizedMessage()); return; } }
	 * 
	 * // Close the file try { Writer.close();
	 * 
	 * } catch (Exception E) { System.err.println(E.getLocalizedMessage());
	 * return; }
	 * 
	 * }
	 */

	/**
	 * @param args
	 */
	/*
	 * public static void main(String[] args) { /*
	 * System.out.println("Running..."); String[] newArgs ={"-r",
	 * "/home/natalie/Projects/Antibody/Mouse/PRMScores/aBTLA_HC_pepsin_30min_042707.prm.txt"
	 * ,"-w","Temp.txt"}; SpecAlign Test = new SpecAlign(newArgs); Test.Run();
	 */
	/*
	 * String[] PRMFileNames = {
	 * "/home/natalie/Projects/MySVNProjects/ImmunoSeq/TestData/PRMs/SpecAlignTest/aBTLA_HC_pepsin_3h_042707.prm"
	 * ,
	 * "/home/natalie/Projects/MySVNProjects/ImmunoSeq/TestData/PRMs/SpecAlignTest/aBTLA_HC_pepsin_30min_042707.prm"
	 * }; String[] OracleFileNames = {
	 * "/home/natalie/Projects/MySVNProjects/ImmunoSeq/TestData/Oracle/aBTLA_HC_pepsin_3h_042707.txt"
	 * ,
	 * "/home/natalie/Projects/MySVNProjects/ImmunoSeq/TestData/Oracle/aBTLA_HC_pepsin_30min_042707.txt"
	 * };
	 * 
	 * String OutputFileName =
	 * "/home/natalie/Projects/MySVNProjects/ImmunoSeq/TestData/SpecAlign_Test4.out"
	 * ;
	 * 
	 * String TestSequence = AntibodyUtils.aBTLA_OracleSequence.substring(256,
	 * 266); // String TestSequence = Utils.aBTLA_OracleSequence;
	 * SpecAlign.TestRun(PRMFileNames, OracleFileNames, TestSequence,
	 * OutputFileName, AntibodyUtils.DefaultFragTolerance); }
	 */

}
