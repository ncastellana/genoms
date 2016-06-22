package antibody;

import java.util.ArrayList;

import basicUtils.Utils;
//import basicUtils.InspectAnnotation;

public class ProteinSeed extends Seed {

	private boolean Debug = false;

	public ProteinSeed(String Sequence) {
		this.SupportingSpectra = new ArrayList();
		this.ClassNum = -1;
		this.SeedSequence = Utils.GetUnModded(Sequence);
		this.AnchoredSeedSequence = this.SeedSequence;
		this.SeedStart = -1;
		this.SeedEnd = -1;
		this.ExtendedSeedEnd = -1;
		this.ExtendedSeedStart = -1;
		this.SeedPRMs = AntibodyUtils.GetSeqPRMScaled(Sequence);
		this.UsedPRMs = this.SeedPRMs;
		this.SeedPRMScores = this.DefaultSetSeedPRMScores();
		this.PeakOverlapping = new Object[this.SeedPRMs.length];
		this.PeakSupporting = new Object[this.SeedPRMs.length];

		this.Mods = AntibodyUtils.GetModsFromPeptide(Sequence);

		this.NextSeed = null;
		this.PrevSeed = null;
	}

	public ProteinSeed(SpectrumNode specNode, int PeptideStart, int PeptideEnd,
			int ClassNum, Template ClassTemplate) {
		this.SupportingSpectra = new ArrayList();
		this.SupportingSpectra.add(specNode);

		this.ClassNum = ClassNum;

		this.OrigTemplate = ClassTemplate;
		this.SeedSequence = specNode.getUnModdedDBAnnotation();
		this.SeedPRMs = AntibodyUtils.GetSeqPRMScaled(this.SeedSequence);
		this.UsedPRMs = this.SeedPRMs;
		this.SeedPRMScores = this.DefaultSetSeedPRMScores();
		this.AnchoredSeedSequence = this.SeedSequence;
		this.PeakOverlapping = new Object[this.SeedPRMs.length];
		this.PeakSupporting = new Object[this.SeedPRMs.length];

		String SpecKey = AntibodyUtils.CreateSpectrumKey(
				specNode.GetSourceFileName(), specNode.GetSpecIndex());
		for (int i = 0; i < this.SeedPRMs.length; ++i) {
			if (SpecKey == null) {
				System.err
						.println("ERROR: Unable to create spectrum key for this spectrum");
				specNode.DebugPrint();
				Utils.WaitForEnter();
			}
			String[] Temp = new String[1];
			Temp[0] = SpecKey;
			this.PeakOverlapping[i] = Temp;

			Temp = new String[1];
			Temp[0] = SpecKey;
			this.PeakSupporting[i] = Temp;

		}
		// this.ClassSequence = ClassTemplate.Sequence;
		// this.Mods = AntibodyUtils.GetModsFromPeptide(Annotation.Annotation);

		this.SeedStart = PeptideStart;
		this.SeedEnd = PeptideEnd;
		this.ExtendedSeedStart = PeptideStart;
		this.ExtendedSeedEnd = PeptideEnd;

		this.NextSeed = null;
		this.PrevSeed = null;

		/*
		 * System.out.println(this.SeedSequence); for(int i = 0; i <
		 * this.SeedSequence.length(); ++i) {
		 * System.out.println(this.SeedSequence.charAt(i) + ":" +
		 * AntibodyUtils.StringArrayToString
		 * ((String[])(this.PeakOverlapping[i]))); System.out.println(" - " +
		 * AntibodyUtils
		 * .StringArrayToString((String[])(this.PeakSupporting[i]))); }
		 * System.out.println(":" +
		 * AntibodyUtils.StringArrayToString((String[])(
		 * this.PeakOverlapping[this.SeedSequence.length()])));
		 * System.out.println(" - " +
		 * AntibodyUtils.StringArrayToString((String[]
		 * )(this.PeakSupporting[this.SeedSequence.length()])));
		 * Utils.WaitForEnter();
		 */

	}

	public ProteinSeed(SpectrumNode specNode, int PeptideStart, int PeptideEnd,
			int ClassNum, String ClassSequence, Seed Prev, Seed Next) {
		this.SupportingSpectra = new ArrayList();
		this.SupportingSpectra.add(specNode);

		this.ClassNum = ClassNum;

		this.SeedSequence = specNode.getUnModdedDBAnnotation();
		this.SeedPRMs = AntibodyUtils.GetSeqPRMScaled(this.SeedSequence);
		this.UsedPRMs = this.SeedPRMs;
		this.SeedPRMScores = this.DefaultSetSeedPRMScores();
		this.AnchoredSeedSequence = this.SeedSequence;
		this.PeakOverlapping = new Object[this.SeedPRMs.length];
		this.PeakSupporting = new Object[this.SeedPRMs.length];
		String SpecKey = AntibodyUtils.CreateSpectrumKey(
				specNode.GetSourceFileName(), specNode.GetSpecIndex());
		for (int i = 0; i < this.SeedPRMs.length; ++i) {
			if (SpecKey == null) {
				System.err
						.println("ERROR: Unable to create spectrum key for this spectrum");
				specNode.DebugPrint();
				Utils.WaitForEnter();
			}
			String[] Temp = new String[1];
			Temp[0] = SpecKey;
			this.PeakOverlapping[i] = Temp;

			Temp = new String[1];
			Temp[0] = SpecKey;
			this.PeakSupporting[i] = Temp;

		}
		// this.ClassSequence = ClassSequence;

		// this.Mods = AntibodyUtils.GetModsFromPeptide(Annotation.Annotation);

		this.SeedStart = PeptideStart;
		this.SeedEnd = PeptideEnd;
		this.ExtendedSeedStart = PeptideStart;
		this.ExtendedSeedEnd = PeptideEnd;

		this.NextSeed = Next;
		this.PrevSeed = Prev;

		/*
		 * System.out.println(this.SeedSequence); for(int i = 0; i <
		 * this.SeedSequence.length(); ++i) {
		 * System.out.println(this.SeedSequence.charAt(i) + ":" +
		 * AntibodyUtils.StringArrayToString
		 * ((String[])(this.PeakOverlapping[i]))); System.out.println(" - " +
		 * AntibodyUtils
		 * .StringArrayToString((String[])(this.PeakSupporting[i]))); }
		 * System.out.println(":" +
		 * AntibodyUtils.StringArrayToString((String[])(
		 * this.PeakOverlapping[this.SeedSequence.length()])));
		 * System.out.println(" - " +
		 * AntibodyUtils.StringArrayToString((String[]
		 * )(this.PeakSupporting[this.SeedSequence.length()])));
		 * Utils.WaitForEnter();
		 */

	}

	public ProteinSeed(Peptide currAnn, int start, int end, int ClassNum,
			Template template) {
		this.SupportingSpectra = new ArrayList();
		this.SupportingSpectra.addAll(currAnn.getSpectra());

		this.ClassNum = ClassNum;

		this.OrigTemplate = template;
		this.SeedSequence = currAnn.getUnModdedPeptideSeq();
		this.SeedPRMs = AntibodyUtils.GetSeqPRMScaled(this.SeedSequence);
		this.UsedPRMs = this.SeedPRMs;
		this.SeedPRMScores = this.DefaultSetSeedPRMScores();
		this.AnchoredSeedSequence = this.SeedSequence;
		this.PeakOverlapping = new Object[this.SeedPRMs.length];
		this.PeakSupporting = new Object[this.SeedPRMs.length];
		for (int s = 0; s < this.SupportingSpectra.size(); ++s) {
			SpectrumNode specNode = (SpectrumNode) (this.SupportingSpectra
					.get(s));

			String SpecKey = AntibodyUtils.CreateSpectrumKey(
					specNode.GetSourceFileName(), specNode.GetSpecIndex());
			for (int i = 0; i < this.SeedPRMs.length; ++i) {
				if (SpecKey == null) {
					System.err
							.println("ERROR: Unable to create spectrum key for this spectrum");
					specNode.DebugPrint();
					Utils.WaitForEnter();
				}
				String[] Temp = new String[1];
				Temp[0] = SpecKey;
				this.PeakOverlapping[i] = Temp;

				Temp = new String[1];
				Temp[0] = SpecKey;
				this.PeakSupporting[i] = Temp;

			}
		}
		// this.ClassSequence = ClassSequence;

		// this.Mods = AntibodyUtils.GetModsFromPeptide(Annotation.Annotation);

		this.SeedStart = start;
		this.SeedEnd = end;
		this.ExtendedSeedStart = start;
		this.ExtendedSeedEnd = end;

		this.NextSeed = null;
		this.PrevSeed = null;
	}

	public static String ReconstructPath(String Input, double Tolerance) {
		String[] Anchors = Input.split("-");
		int Count = 0;

		System.out.println("Total Seeds: " + Anchors.length);
		String FinalPath = "";
		int Index = 0;
		ProteinSeed Head = new ProteinSeed(Anchors[Index]);
		while (Index < Anchors.length) {
			int NextIndex = Index + 1;
			while (NextIndex < Anchors.length && Anchors[NextIndex] == null)
				NextIndex += 1;
			if (NextIndex == Anchors.length) {
				FinalPath += "-" + Head.SeedSequence;
				break;
			}
			ProteinSeed Next = new ProteinSeed(Anchors[NextIndex]);
			if (Head.MergeSequenceWithNext(Next, Tolerance)) {
				Anchors[NextIndex] = null;
				System.out.println("New HeadSeq: " + Head.SeedSequence);
			} else {
				FinalPath += "-" + Head.SeedSequence;
				Index += 1;
				while (Index < Anchors.length - 1 && Anchors[Index] == null)
					Index += 1;
				if (Index == Anchors.length)
					break;
				Head = new ProteinSeed(Anchors[Index]);
			}

		}
		System.out.println("Final Sequence: " + FinalPath.substring(1));
		return FinalPath.substring(1);
	}

	public boolean hasOverlap(int ClassNum, int Start, int End) {
		if (ClassNum != this.ClassNum)
			return false;

		if (Start >= this.ExtendedSeedEnd || End <= this.ExtendedSeedStart)
			return false;
		return true;
	}

	public boolean AddNewAnnotation(int ClassNum, int Start, int End,
			SpectrumNode Annotation) {
		if (!hasOverlap(ClassNum, Start, End))
			return false;

		boolean LocalDebug = false;
		// if(Annotation.UnModdedAnn.compareTo(this.AnchoredSeedSequence) != 0)
		// LocalDebug = true;

		String NewPeptide = Annotation.getUnModdedDBAnnotation();
		if (Start < this.ExtendedSeedStart) {
			this.AnchoredSeedSequence = NewPeptide.substring(0, this.SeedStart
					- Start)
					+ this.AnchoredSeedSequence;
			this.SeedStart = Start;
			this.ExtendedSeedStart = Start;

			if (Debug || LocalDebug) {
				System.out.println("Extended this sequence Forward!!");
				this.DebugPrint();
				Utils.WaitForEnter();
			}
		}
		if (End > this.ExtendedSeedEnd) {
			this.AnchoredSeedSequence = this.AnchoredSeedSequence
					+ NewPeptide.substring(this.SeedEnd - Start);
			this.SeedEnd = End;
			this.ExtendedSeedEnd = End;
			if (Debug || LocalDebug) {
				System.out.println("Extended this sequence backward!!");
				this.DebugPrint();
				Utils.WaitForEnter();
			}
		}
		this.SeedSequence = this.AnchoredSeedSequence;
		this.SeedPRMs = AntibodyUtils.GetSeqPRMScaled(this.SeedSequence);
		this.UsedPRMs = this.SeedPRMs;
		// this.SeedPRMScores = this.DefaultSetSeedPRMScores();
		this.SeedPRMScores = new double[this.SeedPRMs.length];
		this.PeakOverlapping = new Object[this.SeedPRMs.length];
		this.PeakSupporting = new Object[this.SeedPRMs.length];
		// String SpecKey =
		// AntibodyUtils.CreateSpectrumKey(Annotation.SpectrumFile,
		// Annotation.ScanNumber);
		if (LocalDebug)
			System.out
					.println("ADDING NEW PEPTIDE"
							+ Annotation.GetSourceFileName() + ","
							+ Annotation.GetSpecIndex() + ","
							+ Annotation.getUnModdedDBAnnotation() + " to "
							+ this.SeedSequence);
		this.SupportingSpectra.add(Annotation);

		for (int i = 0; i < this.SupportingSpectra.size(); ++i) {
			SpectrumNode CurrAnn = (SpectrumNode) (this.SupportingSpectra
					.get(i));
			String CurrUnModdedAnn = (CurrAnn.getUnModdedDBAnnotation());
			int Pos = this.AnchoredSeedSequence.indexOf(CurrUnModdedAnn);
			if (Pos < 0) {
				System.out.println("Oh no! We can't find " + CurrUnModdedAnn
						+ " in our anchored seequence: "
						+ this.AnchoredSeedSequence);
				System.out.println("Alleged location: " + Start + " - " + End);
				System.out.println("Anchored location: " + this.SeedStart
						+ " - " + this.SeedEnd);
			}
			String SpecKey = AntibodyUtils.CreateSpectrumKey(
					CurrAnn.GetSourceFileName(), CurrAnn.GetSpecIndex());
			if (SpecKey == null) {
				System.err
						.println("ERROR: Unable to create spectrum key for this spectrum");
				Annotation.DebugPrint();
				Utils.WaitForEnter();
			}
			if (LocalDebug) {
				System.out.println("adding to this seed: " + this.SeedSequence);
				System.out.println(SpecKey + ": " + CurrUnModdedAnn);
			}
			for (int j = Pos; j <= Pos + CurrUnModdedAnn.length(); ++j) {
				String[] Temp = new String[1];
				Temp[0] = SpecKey;
				if (this.PeakOverlapping[j] == null)
					this.PeakOverlapping[j] = Temp;
				else
					this.PeakOverlapping[j] = AntibodyUtils.GetUnion(Temp,
							(String[]) (this.PeakOverlapping[j]));
				if (SpecKey == null || this.PeakOverlapping[j] == null) {
					System.err
							.println("ERROR: Unable to create spectrum key for this spectrum");
					Annotation.DebugPrint();
					Utils.WaitForEnter();
				}
				Temp = new String[1];
				Temp[0] = SpecKey;
				if (this.PeakSupporting[j] == null)
					this.PeakSupporting[j] = Temp;
				else
					this.PeakSupporting[j] = AntibodyUtils.GetUnion(Temp,
							(String[]) (this.PeakSupporting[j]));
				if (SpecKey == null || this.PeakSupporting[j] == null) {
					System.err
							.println("ERROR: Unable to create spectrum key for this spectrum(B)");
					Annotation.DebugPrint();
					Utils.WaitForEnter();
				}

			}

		}

		for (int i = 0; i < this.SeedPRMs.length; ++i) {
			/*
			 * It is possible that we have not observed any spectra which
			 * support this peak
			 */
			this.SeedPRMScores[i] = AntibodyUtils.INIT_MATCH_STATE_SCORE
					* Math.max(AntibodyUtils.MIN_ALIGNED_SPECTRA,
							((String[]) (this.PeakSupporting[i])).length);
		}

		return true;
	}

	public boolean addNewAnnotation(int proteinID, int start, int end,
			Peptide currAnn) {
		for (int i = 0; i < currAnn.getNumSpectra(); ++i) {
			if (!this.AddNewAnnotation(proteinID, start, end,
					(SpectrumNode) (currAnn.getSpectrum(i))))
				return false;
		}
		return true;
	}

	public Seed mergeWithAdjacent() {

		ProteinSeed Ret = this;
		// Possibly merge with previous, in whcih case we return the previous
		// seed with adjustments
		if (this.PrevSeed != null
				&& this.PrevSeed.ExtendedSeedEnd > this.ExtendedSeedStart) {
			if (this.Debug)
				System.out.println("Can Merge with Previous ("
						+ this.PrevSeed.ExtendedSeedStart + "-"
						+ this.PrevSeed.ExtendedSeedEnd + ")");
			Ret = (ProteinSeed) this.PrevSeed;
			Ret.NextSeed = this.NextSeed;

			// Adjust the sequence to include the overlap
			if (this.ExtendedSeedStart < Ret.ExtendedSeedStart) {
				Ret.AnchoredSeedSequence = this.AnchoredSeedSequence.substring(
						0, Ret.ExtendedSeedStart - this.ExtendedSeedStart)
						+ Ret.AnchoredSeedSequence;
				Ret.ExtendedSeedStart = this.ExtendedSeedStart;
				Ret.SeedStart = this.SeedStart;
			}
			if (this.ExtendedSeedEnd > Ret.ExtendedSeedEnd) {
				Ret.AnchoredSeedSequence = Ret.AnchoredSeedSequence
						+ this.AnchoredSeedSequence
								.substring(Ret.ExtendedSeedEnd
										- this.ExtendedSeedStart);
				Ret.ExtendedSeedEnd = this.ExtendedSeedEnd;
				Ret.SeedEnd = this.SeedEnd;
			}

			// Transfer the supporting annotations
			for (int i = 0; i < this.SupportingSpectra.size(); ++i) {
				Ret.SupportingSpectra.add(this.SupportingSpectra.get(i));
			}

			if (Debug) {
				System.out.println("After prev merge: ");
				Ret.DebugPrint();
				Utils.WaitForEnter();
			}
		}

		// Possibly merge with teh next seed, Ret is the current seed or the
		// merged seed with the previous
		if (Ret.NextSeed != null
				&& Ret.ExtendedSeedEnd > Ret.NextSeed.ExtendedSeedStart) {

			ProteinSeed FollowingSeed = (ProteinSeed) Ret.NextSeed;
			if (this.Debug)
				System.out.println("Can Merge with Next ("
						+ FollowingSeed.ExtendedSeedStart + "-"
						+ FollowingSeed.ExtendedSeedEnd + ")");

			Ret.NextSeed = Ret.NextSeed.NextSeed;
			if (Ret.NextSeed != null)
				Ret.NextSeed.PrevSeed = Ret;

			// Adjust the sequence to include overlap
			if (FollowingSeed.ExtendedSeedStart < Ret.ExtendedSeedStart) {
				Ret.AnchoredSeedSequence = FollowingSeed.AnchoredSeedSequence
						.substring(0, Ret.ExtendedSeedStart
								- FollowingSeed.ExtendedSeedStart)
						+ Ret.AnchoredSeedSequence;
				Ret.ExtendedSeedStart = FollowingSeed.ExtendedSeedStart;
				Ret.SeedStart = FollowingSeed.SeedStart;
			}

			if (FollowingSeed.ExtendedSeedEnd > Ret.ExtendedSeedEnd) {
				Ret.AnchoredSeedSequence = Ret.AnchoredSeedSequence
						+ FollowingSeed.AnchoredSeedSequence
								.substring(Ret.ExtendedSeedEnd
										- FollowingSeed.ExtendedSeedStart);
				Ret.ExtendedSeedEnd = FollowingSeed.ExtendedSeedEnd;
				Ret.SeedEnd = FollowingSeed.SeedEnd;
			}

			// Transfer the supporting annotations
			for (int i = 0; i < FollowingSeed.SupportingSpectra.size(); ++i) {
				Ret.SupportingSpectra.add(FollowingSeed.SupportingSpectra
						.get(i));
			}
			if (Debug) {
				System.out.println("After following merge: ");
				Ret.DebugPrint();
				Utils.WaitForEnter();
			}

		}
		Ret.SeedSequence = Ret.AnchoredSeedSequence;
		Ret.SeedPRMs = AntibodyUtils.GetSeqPRMScaled(Ret.SeedSequence);
		Ret.UsedPRMs = Ret.SeedPRMs;
		// Ret.SeedPRMScores = Ret.DefaultSetSeedPRMScores();
		Ret.SeedPRMScores = new double[Ret.SeedPRMs.length];

		Ret.PeakOverlapping = new Object[Ret.SeedPRMs.length];
		Ret.PeakSupporting = new Object[Ret.SeedPRMs.length];
		for (int i = 0; i < Ret.SupportingSpectra.size(); ++i) {
			SpectrumNode CurrAnn = (SpectrumNode) (Ret.SupportingSpectra.get(i));
			String CurrUnModdedAnn = (CurrAnn.getUnModdedDBAnnotation());
			int Pos = Ret.AnchoredSeedSequence.indexOf(CurrUnModdedAnn);
			String SpecKey = AntibodyUtils.CreateSpectrumKey(
					CurrAnn.GetSourceFileName(), CurrAnn.GetSpecIndex());
			for (int j = Pos; j <= Pos + CurrUnModdedAnn.length(); ++j) {
				String[] Temp = new String[1];

				if (SpecKey == null) {
					System.err
							.println("ERROR: Unable to create spectrum key for this spectrum");
					Ret.DebugPrint();
					Utils.WaitForEnter();
				}
				Temp[0] = SpecKey;
				if (Ret.PeakOverlapping[j] == null)
					Ret.PeakOverlapping[j] = Temp;
				else
					Ret.PeakOverlapping[j] = AntibodyUtils.GetUnion(Temp,
							(String[]) (Ret.PeakOverlapping[j]));

				Temp = new String[1];
				Temp[0] = SpecKey;
				if (Ret.PeakSupporting[j] == null)
					Ret.PeakSupporting[j] = Temp;
				else
					Ret.PeakSupporting[j] = AntibodyUtils.GetUnion(Temp,
							(String[]) (Ret.PeakSupporting[j]));
				if (SpecKey == null || this.PeakSupporting[j] == null) {
					System.err
							.println("ERROR: Unable to create spectrum key for this spectrum(B)");
					Utils.WaitForEnter();
				}
			}
		}
		for (int i = 0; i < Ret.SeedPRMs.length; ++i) {
			Ret.SeedPRMScores[i] = AntibodyUtils.INIT_MATCH_STATE_SCORE
					* Math.max(AntibodyUtils.MIN_ALIGNED_SPECTRA,
							((String[]) (Ret.PeakSupporting[i])).length);
		}
		return Ret;
	}

	public boolean CanReach(int ProteinID) {
		return this.OrigTemplate.CanReach(ProteinID);
	}

	/**
	 * If the the range Start-End is greater than this anchor, then return -1,
	 * Start-End is less than this anchor's range then return 1. If they
	 * overlap, then return 0.
	 */
	public int compareTo(int ClassNum, int Start, int End) {
		if (Start >= this.ExtendedSeedEnd)
			return -1;
		if (End <= this.ExtendedSeedStart)
			return 1;
		return 0;

	}

	/**
	 * Assumes that Seed sequences do not contain mass gaps, and doesn't allow
	 * for mis predicted AAs (like I/L)
	 * 
	 * @return Returns true if there is overlap
	 */
	public boolean HasSequenceOverlapWithNext(Template[] SeedSequences,
			double Tolerance) {
		if (SeedSequences == null)
			return false;
		boolean LocalDebug = false;
		ProteinSeed TrueNext = (ProteinSeed) this.NextSeed;
		if (TrueNext == null) {
			Template NextTemplate = this.OrigTemplate.NextInList;
			while (NextTemplate != null && NextTemplate.anchors == null)
				NextTemplate = NextTemplate.NextInList;
			if (NextTemplate == null)
				return false;
			TrueNext = (ProteinSeed) NextTemplate.anchors;

		}
		int[] SeqPRMsA = this.UsedPRMs;
		double[] SeqPRMScoresA = this.GetUsedPRMScores(Tolerance);
		Object[] SupportingSpectraA = this.GetSupportingSpectraUsedPRMs();
		Object[] OverlappingSpectraA = this.GetOverlappingSpectraUsedPRMs();

		int[] SeqPRMsB = TrueNext.UsedPRMs;
		double[] SeqPRMScoresB = TrueNext.GetUsedPRMScores(Tolerance);
		Object[] SupportingSpectraB = TrueNext.GetSupportingSpectraUsedPRMs();
		Object[] OverlappingSpectraB = TrueNext.GetOverlappingSpectraUsedPRMs();
		if (LocalDebug) {
			System.out.println("For SeqA");
			System.out.println("PRMs: "
					+ AntibodyUtils.IntegerListToString(this.SeedPRMs));
			System.out.println("Used PRMs: "
					+ AntibodyUtils.IntegerListToString(this.UsedPRMs));

			System.out.println("For SeqB");
			System.out.println("PRMs: "
					+ AntibodyUtils.IntegerListToString(TrueNext.SeedPRMs));
			System.out.println("Used PRMs: "
					+ AntibodyUtils.IntegerListToString(TrueNext.UsedPRMs));
			Utils.WaitForEnter();
			// ConsensusPRMSpectrum Merged =
			// AntibodyUtils.MergePRMsDoubleSkip(this.SeedPRMs,this.SeedPRMScores,
			// TrueNext.SeedPRMs, TrueNext.SeedPRMScores);
			// ConsensusPRMSpectrum Merged =
			// AntibodyUtils.MergePRMsDoubleSkipALeads(SeqPRMsA,SeqPRMScoresA,
			// SeqPRMsB, SeqPRMScoresB);

			System.out.println("Seed 1: " + this.SeedSequence);
			System.out.println("Seed 2: " + TrueNext.SeedSequence);
			Utils.WaitForEnter();
		}
		ConsensusPRMSpectrum Merged = ConsensusPRMSpectrum.MergeAnchors(
				SeqPRMsA, SeqPRMScoresA, OverlappingSpectraA,
				SupportingSpectraA, SeqPRMsB, SeqPRMScoresB,
				OverlappingSpectraB, SupportingSpectraB, Tolerance);
		// String Merged =
		// Utils.MergeSequences(this.SeedSequence,TrueNext.SeedSequence);

		if (Merged == null) {

			// System.out.println("Not mergeable!!");
			return false;
		}
		String MergedString = DeNovo.Reconstruct(Merged.ScaledPeaks,
				Merged.PeakScores, Tolerance);
		// System.out.println("HSOWN T Merged: " + MergedString);
		if (MergedString.length() == 0)
			return false;
		return true;

	}

	/**
	 * Assumes that Seed sequences do not contain mass gaps, and doesn't allow
	 * for mis predicted AAs (like I/L)
	 * 
	 * @return Returns true if there is overlap
	 */
	public boolean HasSequenceOverlapWithPrev(Template[] SeedSequences,
			double Tolerance) {
		if (SeedSequences == null)
			return false;
		ProteinSeed TruePrev = (ProteinSeed) this.PrevSeed;
		if (TruePrev == null) {
			Template PrevTemplate = this.OrigTemplate.PrevInList;
			while (PrevTemplate != null && PrevTemplate.anchors == null)
				PrevTemplate = PrevTemplate.PrevInList;
			if (PrevTemplate == null)
				return false;
			TruePrev = (ProteinSeed) PrevTemplate.anchors;
			while (TruePrev.NextSeed != null)
				TruePrev = (ProteinSeed) TruePrev.NextSeed;

		}

		// int[] SeqPRMsB = AntibodyUtils.GetSeqPRMScaled(this.SeedSequence);
		int[] SeqPRMsB = this.UsedPRMs;
		double[] SeqPRMScoresB = this.GetUsedPRMScores(Tolerance);
		Object[] SupportingSpectraB = this.GetSupportingSpectraUsedPRMs();
		Object[] OverlappingSpectraB = this.GetOverlappingSpectraUsedPRMs();

		int[] SeqPRMsA = TruePrev.UsedPRMs;
		double[] SeqPRMScoresA = TruePrev.GetUsedPRMScores(Tolerance);
		Object[] SupportingSpectraA = TruePrev.GetSupportingSpectraUsedPRMs();
		Object[] OverlappingSpectraA = TruePrev.GetOverlappingSpectraUsedPRMs();
		/*
		 * if(LocalDebug) { System.out.println("For SeqB");
		 * 
		 * System.out.println("PRMs: " +
		 * AntibodyUtils.IntegerListToString(this.SeedPRMs));
		 * System.out.println("Intervals: " +
		 * AntibodyUtils.DoubleArrayToString(this.IntervalScores));
		 * System.out.println("Used PRMs: " +
		 * AntibodyUtils.IntegerListToString(this.UsedPRMs));
		 * System.out.println("Used Intervals: " +
		 * AntibodyUtils.DoubleArrayToString(IntervalScoresB));
		 * 
		 * System.out.println("For SeqA"); System.out.println("PRMs: " +
		 * AntibodyUtils.IntegerListToString(TruePrev.SeedPRMs));
		 * System.out.println("Intervals: " +
		 * AntibodyUtils.DoubleArrayToString(TruePrev.IntervalScores));
		 * System.out.println("Used PRMs: " +
		 * AntibodyUtils.IntegerListToString(TruePrev.UsedPRMs));
		 * System.out.println("Used Intervals: " +
		 * AntibodyUtils.DoubleArrayToString(IntervalScoresA));
		 * Utils.WaitForEnter();
		 * 
		 * }
		 */
		// ConsensusPRMSpectrum Merged =
		// AntibodyUtils.MergePRMsDoubleSkip(this.SeedPRMs,this.SeedPRMScores,
		// TrueNext.SeedPRMs, TrueNext.SeedPRMScores);
		// ConsensusPRMSpectrum Merged =
		// AntibodyUtils.MergePRMsDoubleSkipALeads(SeqPRMsA,SeqPRMScoresA,
		// SeqPRMsB, SeqPRMScoresB);
		ConsensusPRMSpectrum Merged = ConsensusPRMSpectrum.MergeAnchors(
				SeqPRMsA, SeqPRMScoresA, OverlappingSpectraA,
				SupportingSpectraA, SeqPRMsB, SeqPRMScoresB,
				OverlappingSpectraB, SupportingSpectraB, Tolerance);
		// Utils.MergeSequences(this.SeedSequence,TrueNext.SeedSequence);
		// System.out.println("Seed 1: " + TruePrev.SeedSequence);
		// System.out.println("Seed 2: " + this.SeedSequence);
		if (Merged == null) {

			// System.out.println("Not mergeable!!");
			return false;
		}
		String MergedString = DeNovo.Reconstruct(Merged.ScaledPeaks,
				Merged.PeakScores, Tolerance);
		// System.out.println("HSOWP T Merged: " + MergedString);
		if (MergedString.length() == 0)
			return false;
		return true;

	}

	public boolean MergeSequenceWithNext(Template[] SeedSequences,
			double Tolerance) {
		if (SeedSequences == null)
			return false;
		boolean LocalDebug = false;
		boolean crossTemplate = false;
		ProteinSeed TrueNext = (ProteinSeed) this.NextSeed;
		if (TrueNext == null) {
			Template NextTemplate = this.OrigTemplate.NextInList;
			while (NextTemplate != null && NextTemplate.anchors == null)
				NextTemplate = NextTemplate.NextInList;
			if (NextTemplate == null)
				return false;
			TrueNext = (ProteinSeed) NextTemplate.anchors;
			crossTemplate = true;
		}

		if (LocalDebug) {
			System.out.println("Merging SEEDS:");
			System.out.println(this.SeedSequence);
			System.out.println(TrueNext.SeedSequence);
		}
		// int[] SeqPRMsA = AntibodyUtils.GetSeqPRMScaled(this.SeedSequence);
		int[] SeqPRMsA = this.UsedPRMs;
		double[] SeqPRMScoresA = this.GetUsedPRMScores(Tolerance);
		Object[] SupportingSpectraA = this.GetSupportingSpectraUsedPRMs();
		Object[] OverlappingSpectraA = this.GetOverlappingSpectraUsedPRMs();

		int[] SeqPRMsB = TrueNext.UsedPRMs;
		double[] SeqPRMScoresB = TrueNext.GetUsedPRMScores(Tolerance);
		Object[] SupportingSpectraB = TrueNext.GetSupportingSpectraUsedPRMs();
		Object[] OverlappingSpectraB = TrueNext.GetOverlappingSpectraUsedPRMs();
		if (LocalDebug) {
			System.out.println("For SeqA");
			System.out.println("PRMs: "
					+ AntibodyUtils.IntegerListToString(this.SeedPRMs));
			System.out.println("Used PRMs: "
					+ AntibodyUtils.IntegerListToString(this.UsedPRMs));

			System.out.println("For SeqB");
			System.out.println("PRMs: "
					+ AntibodyUtils.IntegerListToString(TrueNext.SeedPRMs));
			System.out.println("Used PRMs: "
					+ AntibodyUtils.IntegerListToString(TrueNext.UsedPRMs));

			Utils.WaitForEnter();

			System.out.println("A: " + Utils.IntArrayToString(SeqPRMsA));
			System.out
					.println("A: " + Utils.DoubleArrayToString(SeqPRMScoresA));
			System.out.println("B: " + Utils.IntArrayToString(SeqPRMsB));
			System.out
					.println("B: " + Utils.DoubleArrayToString(SeqPRMScoresB));
		}
		// ConsensusPRMSpectrum Merged =
		// AntibodyUtils.MergePRMsDoubleSkip(this.SeedPRMs,this.SeedPRMScores,
		// TrueNext.SeedPRMs, TrueNext.SeedPRMScores);
		// ConsensusPRMSpectrum Merged =
		// AntibodyUtils.MergePRMsDoubleSkipALeads(SeqPRMsA,SeqPRMScoresA,
		// SeqPRMsB, SeqPRMScoresB);
		ConsensusPRMSpectrum Merged = ConsensusPRMSpectrum.MergeAnchors(
				SeqPRMsA, SeqPRMScoresA, OverlappingSpectraA,
				SupportingSpectraA, SeqPRMsB, SeqPRMScoresB,
				OverlappingSpectraB, SupportingSpectraB, Tolerance);
		// String Merged =
		// Utils.MergeSequences(this.SeedSequence,TrueNext.SeedSequence);
		if (Merged == null)
			return false;

		String MergedString = DeNovo.Reconstruct(Merged.ScaledPeaks,
				Merged.PeakScores, Tolerance);
		System.out.println("MSWN T Merged: " + MergedString);
		if (MergedString.length() == 0)
			return false;

		// this.SetSeedSequence(MergedString);
		// this.SetSeedPRMs(Merged.ScaledPeaks, Merged.PeakScores,
		// Merged.IntervalScores);
		this.SetSeedPRMs(Merged.ScaledPeaks, Merged.PeakScores,
				Merged.OverlappingSpectra, Merged.SupportingSpectra, Tolerance);

		// Only link to next one if these are on the same template
		if (this.ClassNum == TrueNext.ClassNum)
			this.NextSeed = TrueNext.NextSeed;

		if (TrueNext.PrevSeed == null)
			TrueNext.OrigTemplate.anchors = TrueNext.NextSeed;
		if (TrueNext.NextSeed != null)
			TrueNext.NextSeed.PrevSeed = TrueNext.PrevSeed;

		TrueNext.PrevSeed = null;
		TrueNext.NextSeed = null;
		this.hasMergedAcrossTemplates = crossTemplate;
		return true;
	}

	public boolean MergeSequenceWithNext(ProteinSeed TrueNext, double Tolerance) {

		boolean LocalDebug = false;

		if (LocalDebug) {
			System.out.println("Merging SEEDS:");
			System.out.println(this.SeedSequence);
			System.out.println(TrueNext.SeedSequence);
		}
		// int[] SeqPRMsA = AntibodyUtils.GetSeqPRMScaled(this.SeedSequence);
		int[] SeqPRMsA = this.UsedPRMs;
		double[] SeqPRMScoresA = this.GetUsedPRMScores(Tolerance);
		Object[] SupportingSpectraA = this.GetSupportingSpectraUsedPRMs();
		Object[] OverlappingSpectraA = this.GetOverlappingSpectraUsedPRMs();

		int[] SeqPRMsB = TrueNext.UsedPRMs;
		double[] SeqPRMScoresB = TrueNext.GetUsedPRMScores(Tolerance);
		Object[] SupportingSpectraB = TrueNext.GetSupportingSpectraUsedPRMs();
		Object[] OverlappingSpectraB = TrueNext.GetOverlappingSpectraUsedPRMs();

		if (LocalDebug) {
			System.out.println("For SeqA");
			System.out.println("PRMs: "
					+ AntibodyUtils.IntegerListToString(this.SeedPRMs));

			System.out.println("Used PRMs: "
					+ AntibodyUtils.IntegerListToString(this.UsedPRMs));

			System.out.println("For SeqB");
			System.out.println("PRMs: "
					+ AntibodyUtils.IntegerListToString(TrueNext.SeedPRMs));

			System.out.println("Used PRMs: "
					+ AntibodyUtils.IntegerListToString(TrueNext.UsedPRMs));

			Utils.WaitForEnter();

			System.out.println("A: " + Utils.IntArrayToString(SeqPRMsA));
			System.out
					.println("A: " + Utils.DoubleArrayToString(SeqPRMScoresA));
			System.out.println("B: " + Utils.IntArrayToString(SeqPRMsB));
			System.out
					.println("B: " + Utils.DoubleArrayToString(SeqPRMScoresB));
		}
		// ConsensusPRMSpectrum Merged =
		// AntibodyUtils.MergePRMsDoubleSkip(this.SeedPRMs,this.SeedPRMScores,
		// TrueNext.SeedPRMs, TrueNext.SeedPRMScores);
		// ConsensusPRMSpectrum Merged =
		// AntibodyUtils.MergePRMsDoubleSkipALeads(SeqPRMsA,SeqPRMScoresA,
		// SeqPRMsB, SeqPRMScoresB);
		ConsensusPRMSpectrum Merged = ConsensusPRMSpectrum.MergeAnchors(
				SeqPRMsA, SeqPRMScoresA, OverlappingSpectraA,
				SupportingSpectraA, SeqPRMsB, SeqPRMScoresB,
				OverlappingSpectraB, SupportingSpectraB, Tolerance);
		// String Merged =
		// Utils.MergeSequences(this.SeedSequence,TrueNext.SeedSequence);
		if (Merged == null)
			return false;

		String MergedString = DeNovo.Reconstruct(Merged.ScaledPeaks,
				Merged.PeakScores, Tolerance);
		System.out.println("MSWN Merged: " + MergedString);
		if (MergedString.length() == 0)
			return false;

		this.SetSeedSequence(MergedString);
		// this.SetSeedPRMs(Merged.ScaledPeaks, Merged.PeakScores,
		// Merged.IntervalScores);
		this.SetSeedPRMs(Merged.ScaledPeaks, Merged.PeakScores,
				Merged.OverlappingSpectra, Merged.SupportingSpectra, Tolerance);
		this.UsedPRMs = DeNovo.ReconstructPRMs(Merged.ScaledPeaks,
				Merged.PeakScores, Tolerance);
		// Only link to next one if these are on the same template
		if (this.ClassNum == TrueNext.ClassNum)
			this.NextSeed = TrueNext.NextSeed;

		if (TrueNext.PrevSeed == null)
			TrueNext.OrigTemplate.anchors = TrueNext.NextSeed;
		if (TrueNext.NextSeed != null)
			TrueNext.NextSeed.PrevSeed = TrueNext.PrevSeed;

		TrueNext.PrevSeed = null;
		TrueNext.NextSeed = null;
		return true;
	}

	public Seed GetNextSeed(Template[] SeedSequences) {
		if (SeedSequences == null)
			return null;
		ProteinSeed TrueNext = (ProteinSeed) this.NextSeed;
		if (TrueNext == null) {
			Template NextTemplate = this.OrigTemplate.NextInList;
			while (NextTemplate != null && NextTemplate.anchors == null)
				NextTemplate = NextTemplate.NextInList;
			if (NextTemplate == null)
				return null;
			TrueNext = (ProteinSeed) NextTemplate.anchors;

		}
		return TrueNext;
	}
	
	public boolean Succeeds(Seed PrevSeed, Template[] SeedSequences) {

	if (SeedSequences == null)
		return false;
	ProteinSeed TruePrev = (ProteinSeed) this.PrevSeed;
	//Seed TruePrev = this.PrevSeed;
	if (TruePrev == null) {
		Template PrevTemplate = this.OrigTemplate.PrevInList;
		while (PrevTemplate != null) {
			if (PrevTemplate.anchors == null)
				PrevTemplate = PrevTemplate.PrevInList;
			PrevTemplate = PrevTemplate.PrevInList;
		}
		if (PrevTemplate == null)
			return false;

		Seed Ptr = PrevTemplate.anchors;
		while (Ptr.NextSeed != null)
			Ptr = Ptr.NextSeed;
		return (PrevSeed.Equals(Ptr));
	}
	if (TruePrev.Equals(PrevSeed))
		return true;
	else
		return false;
	}

}
