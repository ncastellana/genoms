package antibody;

import java.util.ArrayList;
import java.util.Hashtable;

import basicUtils.Utils;
import errorUtils.ErrorThrower;

public abstract class Seed implements Comparable<Object> {

	// A list of SpectrumNode objects used to create the initial seed
	// object. This
	// does not include spectra added by extension
	protected ArrayList<SpectrumNode> SupportingSpectra;
	protected int ClassNum;
	public Template OrigTemplate;

	protected String SeedSequence;

	protected String AnchoredSeedSequence;

	// protected String ClassSequence;

	// Mods are sorted by position on the anchor sequence (not the extension)
	public ModSite[] Mods;

	public Seed NextSeed;
	public Seed PrevSeed;

	protected int SeedStart;
	protected int SeedEnd;
	protected int ExtendedSeedStart;
	protected int ExtendedSeedEnd;

	protected int[] SeedPRMs;
	protected double[] SeedPRMScores;
	protected int[] UsedPRMs; // the PRMs that are used in the sequence
	// alignment
	protected double[] IntervalScores; // deprecated
	protected Object[] PeakOverlapping;
	protected Object[] PeakSupporting;

	boolean hasMergedAcrossTemplates = false;

	public abstract boolean hasOverlap(int ClassNum, int Start, int End);

	public abstract boolean AddNewAnnotation(int ClassNum, int Start, int End,
			SpectrumNode Annotation);

	public abstract boolean addNewAnnotation(int proteinID, int start, int end,
			Peptide currAnn);

	public abstract Seed mergeWithAdjacent();

	public abstract int compareTo(int ClassNum, int Start, int End);

	public abstract boolean CanReach(int TemplateID);

	public int compareTo(Object o) {
		Seed Other = (Seed) (o);

		if (Other.ClassNum == this.ClassNum) {
			if (Other.SeedStart >= this.SeedEnd)
				return -1;
			if (Other.SeedEnd <= this.SeedStart)
				return 1;
			return 0;
		} else {
			if (this.OrigTemplate.CanReach(Other.ClassNum))
				return -1;
			if (Other.CanReach(this.ClassNum))
				return 1;
			else {
				ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
						"Seeds are not comparable!");
			}
		}
		return 0;

	}

	public boolean Equals(Seed Other) {
		if (this.ClassNum == Other.ClassNum
				&& this.SeedStart == Other.SeedStart
				&& this.SeedEnd == Other.SeedEnd)
			return true;
		return false;
	}

	/**
	 * Assumes that Seed sequences do not contain mass gaps, and doesn't allow
	 * for mis predicted AAs (like I/L)
	 * 
	 * @return Returns true if there is overlap
	 */
	public abstract boolean HasSequenceOverlapWithNext(
			Template[] SeedSequences, double Tolerance);

	/**
	 * Assumes that Seed sequences do not contain mass gaps, and doesn't allow
	 * for mis predicted AAs (like I/L)
	 * 
	 * @return Returns true if there is overlap
	 */
	public abstract boolean HasSequenceOverlapWithPrev(
			Template[] SeedSequences, double Tolerance);

	public abstract boolean MergeSequenceWithNext(Template[] SeedSequences,
			double Tolerance);

	// public abstract Seed MergeSequenceWithPrev(Template[] SeedSequences);
	/**
	 * 
	 * @param Extension
	 *            The new Seed sequence
	 */
	public void SetSeedSequence(String Extension) {
		this.SeedSequence = Extension;
		this.SeedPRMs = AntibodyUtils.GetSeqPRMScaled(this.SeedSequence);
		this.SeedPRMScores = this.DefaultSetSeedPRMScores();
		this.PeakOverlapping = new Object[this.SeedPRMs.length];
		this.PeakSupporting = new Object[this.SeedPRMs.length];
		// this.IntervalScores = this.DefaultSetSeedIntervalScores();
		this.UsedPRMs = this.SeedPRMs;
	}

	public void SetSeedPRMs(int[] PRMs, double Tolerance) {
		this.SeedPRMs = PRMs;
		this.SeedPRMScores = this.DefaultSetSeedPRMScores();
		this.SeedSequence = DeNovo.Reconstruct(this.SeedPRMs, Tolerance);
		this.UsedPRMs = DeNovo.ReconstructPRMs(this.SeedPRMs, Tolerance);
		this.PeakOverlapping = new Object[this.SeedPRMs.length];
		this.PeakSupporting = new Object[this.SeedPRMs.length];
	}

	public double[] DefaultSetSeedPRMScores() {
		double[] Scores = new double[this.SeedPRMs.length];

		double Score = AntibodyUtils.MIN_ALIGNED_SPECTRA
				* AntibodyUtils.INIT_MATCH_STATE_SCORE;
		for (int i = 0; i < Scores.length; ++i)
			Scores[i] = Score;
		return Scores;
	}

	public double[] DefaultSetSeedIntervalScores() {
		double[] Scores = new double[this.SeedPRMs.length - 1];

		double Score = ((double) AntibodyUtils.PSEUDOCOUNT_MATCH)
				/ (AntibodyUtils.PSEUDOCOUNT_MATCH + AntibodyUtils.PSEUDOCOUNT_DELETE);
		for (int i = 0; i < Scores.length; ++i)
			Scores[i] = Score;
		return Scores;
	}

	public void SetSeedPRMs(int[] PRMs, double[] Scores, Object[] Overlapping,
			Object[] Supporting, double Tolerance) {
		this.SeedPRMs = PRMs;
		this.SeedPRMScores = Scores;
		this.SeedSequence = DeNovo.Reconstruct(this.SeedPRMs,
				this.SeedPRMScores, Tolerance);
		this.UsedPRMs = DeNovo.ReconstructPRMs(this.SeedPRMs,
				this.SeedPRMScores, Tolerance);
		this.PeakOverlapping = Overlapping;
		this.PeakSupporting = Supporting;
	}

	public void SetSeedPRMs(String Sequence, int[] PRMs, double[] Scores,
			Object[] Overlapping, Object[] Supporting, double Tolerance) {
		int[] SeqPRMs = AntibodyUtils.GetSeqPRMScaled(Sequence);

		this.SeedPRMs = new int[SeqPRMs.length];
		this.SeedPRMScores = new double[SeqPRMs.length];
		for (int i = 0; i < SeqPRMs.length; ++i) {
			int SeqMass = SeqPRMs[i];
			for (int j = 0; j < PRMs.length; ++j) {
				if (AntibodyUtils.EquivalentScaled(SeqMass, PRMs[j], Tolerance)) {
					this.SeedPRMs[i] = SeqMass;
					this.SeedPRMScores[i] = Scores[j];
					break;
				}
			}
		}
		this.SeedSequence = Sequence;
		this.PeakOverlapping = Overlapping;
		this.PeakSupporting = Supporting;
	}

	public int[] GetSeedPRMs() {
		return this.SeedPRMs;
	}

	public Object[] GetSupportingSpectraUsedPRMs() {
		Object[] Ret = new Object[this.UsedPRMs.length];
		for (int i = 0; i < this.UsedPRMs.length; ++i) {
			for (int j = 0; j < this.SeedPRMs.length; ++j) {
				if (this.UsedPRMs[i] == this.SeedPRMs[j]) {
					Ret[i] = this.PeakSupporting[j];
					break;
				}
			}
		}
		return Ret;
	}

	public Object[] GetOverlappingSpectraUsedPRMs() {
		Object[] Ret = new Object[this.UsedPRMs.length];
		for (int i = 0; i < this.UsedPRMs.length; ++i) {
			for (int j = 0; j < this.SeedPRMs.length; ++j) {
				if (this.UsedPRMs[i] == this.SeedPRMs[j]) {
					Ret[i] = this.PeakOverlapping[j];
					break;
				}
			}

		}
		return Ret;
	}

	/**
	 * Retrieves the collection of spectra that have peaks supporting a given
	 * PRM spectrum
	 * 
	 * @param ThesePRMs
	 *            The list of scaled PRM masses to find support for
	 * @param Tolerance
	 *            Mass tolerance (unscaled)
	 * @return Returns an array, such that the ith element is a String array
	 *         containing the spectrum keys for spectra supporting the ith PRM
	 *         in ThesePRMs. The spectrum key is the SpecFile_ScanNum
	 */
	public Object[] GetSupportingSpectra(int[] ThesePRMs, double Tolerance) {
		Object[] Ret = new Object[ThesePRMs.length];
		for (int i = 0; i < ThesePRMs.length; ++i) {
			boolean Found = false;
			for (int j = 0; j < this.SeedPRMs.length; ++j) {
				if (AntibodyUtils.EquivalentScaled(ThesePRMs[i],
						this.SeedPRMs[j], Tolerance)) {
					Ret[i] = this.PeakSupporting[j];
					Found = true;
					break;
				}
			}
			if (!Found) {
				System.out
						.println("ERROR[GetSupportingSpectra]: Could not find "
								+ ThesePRMs[i]
								+ " in "
								+ AntibodyUtils
										.IntegerListToString(this.SeedPRMs));
				Utils.WaitForEnter();
			}
		}
		return Ret;
	}

	/**
	 * Retrieves the collection of spectra that have peaks overlapping a given
	 * PRM spectrum. Peaks may not necessarily be in support of the PRM spectrum
	 * 
	 * @param ThesePRMs
	 *            The list of scaled PRM masses to check for overlap of
	 * @param Tolerance
	 *            The mass tolerance (un-scaled)
	 * @return Returns an array, such that the ith element is a String array
	 *         containing the spectrum keys for spectra overalapping the ith PRM
	 *         in These PRMs. The spectrum key is the SpecFile_ScanNum
	 */
	public Object[] GetOverlappingSpectra(int[] ThesePRMs, double Tolerance) {
		Object[] Ret = new Object[ThesePRMs.length];
		for (int i = 0; i < ThesePRMs.length; ++i) {
			boolean Found = false;
			for (int j = 0; j < this.SeedPRMs.length; ++j) {
				if (AntibodyUtils.EquivalentScaled(ThesePRMs[i],
						this.SeedPRMs[j], Tolerance)) {
					Ret[i] = this.PeakOverlapping[j];
					Found = true;
					break;
				}
			}
			if (!Found) {
				System.out
						.println("ERROR[GetOverlappingSpectra]: Could not find "
								+ ThesePRMs[i]
								+ " in "
								+ AntibodyUtils
										.IntegerListToString(this.SeedPRMs));
				Utils.WaitForEnter();
			}
		}
		return Ret;
	}

	public double[] GetUsedPRMScores(double Tolerance) {
		double[] SeqPRMScoresA = new double[this.UsedPRMs.length];
		for (int j = 0; j < this.UsedPRMs.length; ++j) {
			for (int i = 0; i < this.SeedPRMs.length; ++i) {
				if (AntibodyUtils.EquivalentScaled(this.SeedPRMs[i],
						this.UsedPRMs[j], Tolerance)) {
					SeqPRMScoresA[j] = this.SeedPRMScores[i];
					break;
				}
			}
		}
		return SeqPRMScoresA;
	}

	/**
	 * Looks at peak number Num in the list of UsedPRMs, and retursn all PRMs
	 * preceding that peak. The peaks are shifted to begin at 0.
	 * 
	 * @param Num
	 * @return
	 */
	public int[] GetSeedPRMsLeft(int Num) {
		boolean LocalDebug = false;
		int HighestPeak = this.UsedPRMs[Math.min(Num, this.UsedPRMs.length - 1)];
		if (LocalDebug)
			System.out.println("HighestPeak: " + HighestPeak + " of "
					+ Utils.IntArrayToString(this.UsedPRMs));
		ArrayList Temp = new ArrayList();
		for (int i = 0; i < this.SeedPRMs.length; ++i)
			if (this.SeedPRMs[i] <= HighestPeak)
				Temp.add(new Integer(this.SeedPRMs[i]));

		int[] Ret = AntibodyUtils.ConvertIntegerArrayList(Temp);
		if (LocalDebug)
			System.out
					.println("Selected Peaks: " + Utils.IntArrayToString(Ret));
		for (int i = Ret.length - 1; i >= 0; --i)
			Ret[i] -= Ret[0];
		if (LocalDebug) {
			System.out.println("Shifted Peaks: " + Utils.IntArrayToString(Ret));
			Utils.WaitForEnter();
		}
		return Ret;

		/*
		 * String[] AAs =
		 * AntibodyUtils.SplitStringtoAA(this.AnchoredSeedSequence); String
		 * FinalStr = ""; for(int i = 0; i < Math.min(AAs.length,Num+1); ++i)
		 * FinalStr += AAs[i]; return AntibodyUtils.GetSeqPRMScaled(FinalStr);
		 */
	}

	/**
	 * Looks at peak number Length - Num in the list of UsedPRMs, and retursn
	 * all PRMs succeding that peak. The peaks are shifted to begin at 0.
	 * 
	 * @param Num
	 * @return
	 */
	public int[] GetSeedPRMsRight(int Num) {
		boolean LocalDebug = false;
		int LowestPeak = this.UsedPRMs[Math.max(0, this.UsedPRMs.length
				- (Num + 1))];
		if (LocalDebug)
			System.out.println("LowestPeak: " + LowestPeak + " of "
					+ Utils.IntArrayToString(this.UsedPRMs));
		ArrayList Temp = new ArrayList();
		for (int i = 0; i < this.SeedPRMs.length; ++i)
			if (this.SeedPRMs[i] >= LowestPeak)
				Temp.add(new Integer(this.SeedPRMs[i]));

		int[] Ret = AntibodyUtils.ConvertIntegerArrayList(Temp);
		if (LocalDebug)
			System.out
					.println("Selected Peaks: " + Utils.IntArrayToString(Ret));
		for (int i = Ret.length - 1; i >= 0; --i)
			Ret[i] -= Ret[0];
		if (LocalDebug) {
			System.out.println("Shifted Peaks: " + Utils.IntArrayToString(Ret));
			Utils.WaitForEnter();
		}
		return Ret;

		/*
		 * String[] AAs =
		 * AntibodyUtils.SplitStringtoAA(this.AnchoredSeedSequence); String
		 * FinalStr = ""; for(int i = Math.max(0, AAs.length - Num-1); i <
		 * AAs.length; ++i) FinalStr += AAs[i]; return
		 * AntibodyUtils.GetSeqPRMScaled(FinalStr);
		 */
	}

	public abstract Seed GetNextSeed(Template[] SeedSequences);

	/**
	 * Return only PRMs with mass greater than mass shift, and return with
	 * shifted masses.
	 * 
	 * @param MassShift
	 * @return
	 */
	public int[] GetSeedPRMs(int MassShift) {
		ArrayList Temp = new ArrayList();
		for (int i = 0; i < this.SeedPRMs.length; ++i) {
			if (this.SeedPRMs[i] >= MassShift)
				Temp.add(new Integer(this.SeedPRMs[i] - MassShift));
		}
		return AntibodyUtils.ConvertIntegerArrayList(Temp);

	}

	/*
	 * public double[] GetSeedIntervalScores() { return this.IntervalScores; }
	 */

	public Object[] GetSeedOverlappingCounts() {
		return this.PeakOverlapping;
	}

	public Object[] GetSeedSupportingCounts() {
		return this.PeakSupporting;
	}

	public double[] GetSiteProbs(double Tolerance) {
		int[] ThesePRMs = this.UsedPRMs;
		double[] Ret = new double[ThesePRMs.length - 1];
		boolean LocalDebug = false;

		Object[] Overlapping = this.GetOverlappingSpectra(ThesePRMs, Tolerance);
		Object[] Supporting = this.GetSupportingSpectra(ThesePRMs, Tolerance);
		if (LocalDebug)
			System.out.println("Scoring for seed: " + this.SeedSequence);
		for (int i = 0; i < Ret.length; ++i) {
			String[] S = AntibodyUtils.GetIntersection(
					(String[]) (Supporting[i]), (String[]) (Supporting[i + 1]));
			String[] O = AntibodyUtils.GetIntersection(
					(String[]) (Overlapping[i]),
					(String[]) (Overlapping[i + 1]));
			Ret[i] = ((double) (S.length + 1)) / (O.length + 2);
			if (LocalDebug)
				System.out.println("[" + i + "]: " + ThesePRMs[i] + ":"
						+ (S.length) + "+1/" + (O.length) + "+2");
		}

		return Ret;
	}

	public int[] GetUsedPRMs() {
		return this.UsedPRMs;
	}

	public void DebugPrint() {
		System.out.println("------------------SEED--------------------");
		System.out.println("Original Sequence (" + this.SeedStart + "-"
				+ this.SeedEnd + ") : " + this.AnchoredSeedSequence);
		System.out.println("Extended Sequence (" + this.ExtendedSeedStart + "-"
				+ this.ExtendedSeedEnd + ") : " + this.SeedSequence);
		System.out.println("Supporting Spectra from InsPecT:"
				+ this.SupportingSpectra.size());
		System.out.println("PRMs: "
				+ AntibodyUtils.IntegerListToString(this.SeedPRMs));

		System.out.println("Scores: "
				+ Utils.DoubleArrayToString(this.SeedPRMScores));
		System.out.println("UsedPRMs: "
				+ AntibodyUtils.IntegerListToString(this.UsedPRMs));
		String Str = "";
		for (int i = 0; i < this.PeakOverlapping.length; ++i) {
			String[] O = (String[]) (this.PeakOverlapping[i]);
			String[] S = (String[]) (this.PeakSupporting[i]);
			if (S != null)
				Str += "[" + S.length + "/";
			else
				Str += "[0/";
			if (O != null)
				Str += O.length + "] ";
			else
				Str += "0] ";
		}
		System.out.println(Str);
		if (this.PrevSeed != null)
			System.out.println("Prev: " + this.PrevSeed.SeedSequence);
		else
			System.out.println("Prev: null");
		if (this.NextSeed != null)
			System.out.println("Next: " + this.NextSeed.SeedSequence);
		else
			System.out.println("Next: null");
		System.out.println("--------------------------------------");
		// for(int i = 0; i < this.SupportingSpectra.size(); ++i)
		// {
		// ((InspectAnnotation)(this.SupportingSpectra.get(i))).DebugPrint();
		// }

	}

	public String GetSeedSequence() {
		return this.SeedSequence;
	}

	public int[] GetScaledSeedPRMs() {
		return this.SeedPRMs;
	}

	public double[] GetSeedPRMScores() {
		return this.SeedPRMScores;
	}

	/**
	 * Returns true if this seed succeeds PrevSeed directly. Right now, this
	 * will only happen if they are on the same template but we could update to
	 * look for previous templates
	 * 
	 * @param PrevSeed
	 * @return
	 */
	public boolean Succeeds(Seed PrevSeed) {
		if (this.PrevSeed != null && this.PrevSeed.Equals(PrevSeed))
			return true;
		return false;
	}

	public abstract boolean Succeeds(Seed PrevSeed, Template[] SeedSequences);
	/**
	 * For use after mutation search or when optional mods are allowed. Some of
	 * the peptides supporting this anchor may suggest a mutation/mod, we use a
	 * voting system to decide which mutations to keep
	 */
	public void ReconcileSequence() {
		boolean LocalDebug = true;
		int[] UnModded = new int[this.AnchoredSeedSequence.length()]; // Count
		// of
		// unmodded
		// spectra
		// at
		// that
		// position
		int[] Modded = new int[this.AnchoredSeedSequence.length()]; // Count of
		// modded
		// spectra
		// at that
		// position
		String[] ModdedSuggestion = new String[this.AnchoredSeedSequence
				.length()]; // Different mod types at this position
		if (LocalDebug)
			System.out.println("Reconciling on : " + this.AnchoredSeedSequence);
		for (int i = 0; i < this.SupportingSpectra.size(); ++i) {
			SpectrumNode CurrAnn = (SpectrumNode) (this.SupportingSpectra
					.get(i));
			String ModdedPep = CurrAnn.getModdedDBAnnotation();
			if (ModdedPep.indexOf('.') >= 0)
				ModdedPep = ModdedPep.substring(2, ModdedPep.length() - 2);
			String CurrUnModdedAnn = CurrAnn.getUnModdedDBAnnotation();
			if (ModdedPep.compareTo(CurrUnModdedAnn) == 0) // If this spectrum
															// is unmodified,
															// then ignore it
			{
				int Index = this.AnchoredSeedSequence.indexOf(CurrUnModdedAnn);
				for (int j = Index; j < Index + CurrUnModdedAnn.length(); ++j)
					UnModded[j] += 1;
				continue;
			}
			// String NewPeptide = InspectAnnotation.AdjustToMutations(CurrAnn);
			ModSite[] Sites = AntibodyUtils.GetModsFromPeptide(CurrAnn
					.getModdedDBAnnotation());
			int TrueIndex = this.AnchoredSeedSequence.indexOf(CurrUnModdedAnn);
			for (int j = 0; j < CurrUnModdedAnn.length(); ++j) {
				boolean Found = false;
				for (int k = 0; k < Sites.length; ++k) {
					// This mod appears at this position
					if (Sites[k].ModIndex == j) {

						Modded[TrueIndex + Sites[k].ModIndex] += 1;
						if (ModdedSuggestion[TrueIndex + Sites[k].ModIndex] == null)
							ModdedSuggestion[TrueIndex + Sites[k].ModIndex] = Sites[k].ModString;
						else
							ModdedSuggestion[TrueIndex + Sites[k].ModIndex] += ","
									+ Sites[k].ModString;
						Found = true;
						break;
					}
				}
				if (!Found) {
					UnModded[TrueIndex + j] += 1;
				}
			}
		}

		ArrayList Mods = new ArrayList();
		for (int i = 0; i < this.AnchoredSeedSequence.length(); ++i) {
			if (Modded[i] > 0) {
				if (LocalDebug)
					System.out.println("Mod at " + i + " modded: " + Modded[i]
							+ " unmodded: " + UnModded[i] + ", mods: "
							+ ModdedSuggestion[i]);
				if (UnModded[i] > 0 && UnModded[i] > Modded[i]) {
					if (LocalDebug) {
						System.out
								.println("Have unmodded evidence, cannot change");
						Utils.WaitForEnter();
					}
					ModdedSuggestion[i] = null;
					continue;
				}
				if (ModdedSuggestion[i] == null) {
					if (LocalDebug) {
						System.out.println("Have no mutations at this site");
						Utils.WaitForEnter();
					}
					continue;
				}
				String[] Suggestions = ModdedSuggestion[i].split(",");
				Hashtable SuggestionsTable = new Hashtable();
				for (int j = 0; j < Suggestions.length; ++j) {
					if (!SuggestionsTable.containsKey(Suggestions[j]))
						SuggestionsTable.put(Suggestions[j], new Integer(0));
					int Count = ((Integer) (SuggestionsTable
							.get(Suggestions[j]))).intValue();
					Count += 1;
					SuggestionsTable.put(Suggestions[j], new Integer(Count));
				}
				if (SuggestionsTable.size() > 1) {
					if (LocalDebug) {
						System.out
								.println("Have too many mutations at this site: "
										+ ModdedSuggestion[i]);
						Utils.WaitForEnter();
					}
					ModdedSuggestion[i] = null;
					continue;
				}
				int AnchoredIndex = this.SeedSequence
						.indexOf(this.AnchoredSeedSequence);
				if (AnchoredIndex >= 0) {

					ModSite M = new ModSite(Suggestions[0], AnchoredIndex + i,
							AntibodyUtils.GetModMassFromName(Suggestions[0]),
							this.SeedSequence);
					Mods.add(M);
					if (LocalDebug) {
						System.out.println("ADDED A NEW MOD!");
						M.DebugPrint();

						Utils.WaitForEnter();
					}
				} else {
					System.out
							.println("Anchored sequence cannot be found in extended sequence");
				}

			}
		}

		if (Mods.size() > 0) {
			StringBuffer B = new StringBuffer(this.AnchoredSeedSequence);
			this.Mods = new ModSite[Mods.size()];
			for (int i = 0; i < this.Mods.length; ++i) {
				this.Mods[i] = (ModSite) (Mods.get(i));
				char NewAA = this.Mods[i].GetMutationSub();
				if (NewAA != 0)
					B.replace(this.Mods[i].ModIndex, this.Mods[i].ModIndex + 1,
							Character.toString(NewAA));
			}
			this.SeedSequence = B.toString();
			System.out.println("AnchoredSeedSeq: " + this.AnchoredSeedSequence);
			System.out.println("NewSeedSeq: " + this.SeedSequence);
		}

	}

	public int[] GetSiteOverlap(double Tolerance) {
		// TODO Auto-generated method stub
		int[] ThesePRMs = this.UsedPRMs;
		int[] Ret = new int[ThesePRMs.length - 1];
		boolean LocalDebug = false;

		Object[] Overlapping = this.GetOverlappingSpectra(ThesePRMs, Tolerance);
		if (LocalDebug)
			System.out.println("Scoring for seed: " + this.SeedSequence);
		for (int i = 0; i < Ret.length; ++i) {
			String[] O = AntibodyUtils.GetIntersection(
					(String[]) (Overlapping[i]),
					(String[]) (Overlapping[i + 1]));
			Ret[i] = O.length;
		}

		return Ret;
	}

	public int[] GetSiteSupport(double Tolerance) {
		// TODO Auto-generated method stub
		int[] ThesePRMs = this.UsedPRMs;
		int[] Ret = new int[ThesePRMs.length - 1];
		boolean LocalDebug = false;

		Object[] Supporting = this.GetSupportingSpectra(ThesePRMs, Tolerance);
		if (LocalDebug)
			System.out.println("Scoring for seed: " + this.SeedSequence);
		for (int i = 0; i < Ret.length; ++i) {
			String[] O = AntibodyUtils.GetIntersection(
					(String[]) (Supporting[i]), (String[]) (Supporting[i + 1]));
			Ret[i] = O.length;
		}

		return Ret;
	}

	public int GetExtendedSeedPositionInTemplateStrMatch() {
		// Allow some skipped bases from the anchored guy
		String[] AAs = AntibodyUtils.SplitStringtoAA(this.SeedSequence);
		for (int i = 0; i < 3; ++i) {
			// System.out.println("ExtendedSeq: " + this.SeedSequence);
			// System.out.println("AnchoredSeq: " + this.AnchoredSeedSequence);
			System.out
					.println("DEBUG: Seed.GetExtendedSeedPositionInTemplateStrMatch()");
			System.out.println("DEBUG: " + this.SeedSequence);
			System.out.println("DEBUG: "
					+ this.AnchoredSeedSequence.substring(i,
							this.AnchoredSeedSequence.length() - i));
			System.out.println("DEBUG: Buffer: " + i);

			int index = AntibodyUtils.MatchStringToStringArray(AAs,
					this.AnchoredSeedSequence.substring(i,
							this.AnchoredSeedSequence.length() - i));
			if (index >= 0) {
				return this.SeedStart - (index - i);
			}
		}
		return Integer.MIN_VALUE;
	}

	/**
	 * Aligns the extended seed to the template.
	 * 
	 * @param Tol
	 *            PRM tolerance for alignment
	 * @return Returns a 2D array of 2 rows X N columns where N is the length of
	 *         the alignment. The ith column of the 0th row is the index of the
	 *         seed PRM aligned to the template PRM index in the 1st row of
	 *         column i.
	 */
	public int[][] GetExtendedSeedAlignmentToTemplate(double Tol) {
		/*
		 * Try the simple alignment first (but doing string matching).
		 */
		// Allow some skipped bases from the anchored guy
		String[] AAs = AntibodyUtils.SplitStringtoAA(this.SeedSequence);
		for (int i = 0; i < 3; ++i) {
			// System.out.println("ExtendedSeq: " + this.SeedSequence);
			// System.out.println("AnchoredSeq: " + this.AnchoredSeedSequence);
			int index = AntibodyUtils.MatchStringToStringArray(AAs,
					this.AnchoredSeedSequence.substring(i,
							this.AnchoredSeedSequence.length() - i));
			if (index >= 0) // We found the start of an alignment!!!
			{
				// Figure out the mass offset induced by the alignment starting
				// here
				int seedPRMAtIndex = this.UsedPRMs[index];
				int[] templatePRMs = AntibodyUtils
						.GetSeqPRMScaled(this.OrigTemplate.Sequence);
				int templatePRMAtIndex = templatePRMs[this.SeedStart
						- (index - i)];

				System.out.println(this.SeedSequence);
				System.out.println(this.AnchoredSeedSequence);
				System.out
						.println("Anchored seed pos in extended seed" + index);
				System.out.println("Extended seed pos in template: "
						+ (this.SeedStart - (index - i)));
				System.out.println(seedPRMAtIndex);
				System.out.println(templatePRMAtIndex);
				System.out.println(Utils.JoinIntArray(this.UsedPRMs, " "));
				System.out.println(Utils.JoinIntArray(this.SeedPRMs, " "));
				System.out.println(Utils.JoinIntArray(templatePRMs, " "));
				Utils.WaitForEnter();

			}

		}

		Object[] alnDetails = AntibodyUtils
				.GetLongestAlignmentStrictDetailedSeqs(
						this.AnchoredSeedSequence, this.SeedSequence, Tol);
		int[][] Alignment = (int[][]) (alnDetails[1]);
		int[][] prevCell = (int[][]) (alnDetails[0]);
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
		int[][] ret = new int[2][OverallBestScore + 1];
		int alnIndex = OverallBestScore;

		int prevA = BestA;
		int prevB = BestB;
		while (prevCell[prevA][prevB] >= 0) {
			ret[0][alnIndex] = BestA;
			ret[1][alnIndex] = BestB;
			alnIndex -= 1;
			int currCell = prevCell[prevA][prevB];
			prevA = currCell / Alignment.length;
			prevB = currCell % Alignment.length;
		}
		System.out.println(Utils.JoinIntArray(ret[0], " "));
		System.out.println(Utils.JoinIntArray(ret[1], " "));
		Utils.WaitForEnter();
		return ret;
	}

	public int GetExtendedSeedPositionInTemplate(double Tol) {

		// System.out.println("Looking for alignment between " +
		// this.AnchoredSeedSequence + " and " + this.SeedSequence);
		// Try just the string matching, maybe we'll get lucky!
		int ret = 0;
		ret = GetExtendedSeedPositionInTemplateStrMatch();
		if (ret > Integer.MIN_VALUE)
			return ret;

		// Now find the best alignment

		int[] alignment = AntibodyUtils.GetLongestAlignmentStrictSeqs(
				this.AnchoredSeedSequence, this.SeedSequence, Tol);
		if (alignment[0] > 3) {
			// System.out.println("Found alignment of length " + alignment[0] +
			// " between " + this.AnchoredSeedSequence + " and " +
			// this.SeedSequence);
			// /System.out.println("Location in anchor: " + alignment[1]);
			// System.out.println("Location in extension: " + alignment[2]);
			// System.out.println("Position of anchor: " + this.SeedStart);
			return this.SeedStart + alignment[1] - alignment[2];
		}

		/*
		 * //Allow some skipped bases from the anchored guy int[] AnchoredPRMs =
		 * AntibodyUtils.GetSeqPRMScaled(this.AnchoredSeedSequence);
		 * System.out.println("Trying to find position of " +
		 * this.AnchoredSeedSequence + " in " + this.SeedSequence);
		 * 
		 * 
		 * for(int i = 0; i < AnchoredPRMs.length; ++i) { for(int j = 0; j <
		 * this.UsedPRMs.length; ++j) { if(Math.abs(AnchoredPRMs[i] -
		 * this.UsedPRMs[j]) < ScaledTolerance) {
		 * //System.out.println("I believe the alignment starts at " +
		 * (this.SeedStart + i - j)); System.out.println("Match for masses " +
		 * AnchoredPRMs[i] + "~=" + this.UsedPRMs[j] + "<" + ScaledTolerance);
		 * System.out.println("diff: " + Math.abs(AnchoredPRMs[i] -
		 * this.UsedPRMs[j])); System.out.println("Position in anchor: " + i);
		 * System.out.println("Position in extension: " + j);
		 * System.out.println("Position of anchor: " + this.SeedStart);
		 * 
		 * return this.SeedStart + i - j; } } }
		 */
		return Integer.MIN_VALUE;
	}

	/**
	 * Takes a spectrum and aligns it to this seed and returns the PRMs of the
	 * spectrum shifted to accomodate the alignment. The PRM spectrum is
	 * expected to be noiser than the seeds
	 * 
	 * @param specNode
	 * @return
	 */
	public int[] GetShiftedPRMSpectrum(SpectrumNode specNode, double Tol) {

		int[] specPRMs = specNode.GetPRMSpectrum().PRMs;
		Object[] alignmentInfo = AntibodyUtils.GetLongestAlignmentDetailedPRMs(
				specPRMs, this.SeedPRMs, Tol, Integer.MAX_VALUE, 3);

		int[][] prevCell = (int[][]) (alignmentInfo[0]);
		int[][] alignment = (int[][]) (alignmentInfo[1]);

		int bestScore = 0;
		int bestRow = -1;
		int bestCol = -1;

		for (int i = 0; i < alignment.length; ++i) {
			for (int j = 0; j < alignment[0].length; ++j) {
				if (alignment[i][j] > bestScore) {
					bestScore = alignment[i][j];
					bestRow = i;
					bestCol = j;
				}
			}
		}

		int prevRowMass = -1;
		int prevColMass = -1;

		int currRow = bestRow;
		int currCol = bestCol;

		int qSum = 0;
		int L = 0;

		while (prevCell[currRow][currCol] >= 0) {

			// if(prevRowMass >= 0)
			// {
			// int rowDelta = prevRowMass - specPRMs[currRow];
			// int colDelta = prevColMass - this.SeedPRMs[currCol];
			// qSum += (colDelta - rowDelta);
			// L += 1;
			// }
			System.out.println("[" + currRow + "," + currCol + "] , delta = "
					+ (this.SeedPRMs[currCol] - specPRMs[currRow]));
			qSum += this.SeedPRMs[currCol] - specPRMs[currRow];
			L += 1;
			// prevRowMass = specPRMs[currRow];
			// prevColMass = specPRMs[currCol];
			int currCell = prevCell[currRow][currCol];
			currRow = currCell / prevCell[0].length;
			currCol = currCell % prevCell[0].length;
		}
		int massShift = qSum / L;

		System.out.println("Alignment length: " + L);
		// System.out.println(" start in spectrum: " + startSpecIndex);
		// System.out.println(" start in seed PRMs: " + startSeedIndex);
		System.out.println("MassShift: " + massShift);
		int[] ret = new int[specNode.GetPRMSpectrum().PRMs.length];
		for (int i = 0; i < ret.length; ++i) {
			ret[i] = specNode.GetPRMSpectrum().PRMs[i] + massShift;
		}

		return ret;
	}

	/**
	 * Takes a spectrum (this time only a filtered set of peaks), and returns
	 * the PRMs of the spectrum, shifted to accomodate the alignment. The PRM
	 * spectrum is not expected to be noisy
	 * 
	 * @param alignedPeaks
	 * @param prmTolerance
	 * @return
	 */
	public int[] GetShiftedPRMSpectrum(int[] alignedPeaks, double prmTolerance) {

		int massShift = this.GetMassShift(alignedPeaks, prmTolerance);
		System.out.println("MASS SHIFT: " + massShift);
		int[] ret = new int[alignedPeaks.length];
		for (int i = 0; i < ret.length; ++i) {
			ret[i] = alignedPeaks[i] + massShift;
		}

		return ret;
	}

	public int GetMassShift(int[] alignedPeaks, double prmTolerance) {

		Object[] alignmentInfo = AntibodyUtils.GetLongestAlignmentDetailedPRMs(
				alignedPeaks, this.SeedPRMs, prmTolerance, 2, 3);

		int[][] prevCell = (int[][]) (alignmentInfo[0]);
		int[][] alignment = (int[][]) (alignmentInfo[1]);

		int bestScore = 0;
		int bestRow = -1;
		int bestCol = -1;

		for (int i = 0; i < alignment.length; ++i) {
			for (int j = 0; j < alignment[0].length; ++j) {
				if (alignment[i][j] > bestScore) {
					bestScore = alignment[i][j];
					bestRow = i;
					bestCol = j;
				}
			}
		}

		int currRow = bestRow;
		int currCol = bestCol;

		int qSum = 0;
		int L = 0;

		while (prevCell[currRow][currCol] >= 0) {

			// if(prevRowMass >= 0)
			// {
			// int rowDelta = prevRowMass - specPRMs[currRow];
			// int colDelta = prevColMass - this.SeedPRMs[currCol];
			// qSum += (colDelta - rowDelta);
			// L += 1;
			// }
			// System.out.println("[" + currRow + "," + currCol +
			// "] , delta = (" + this.SeedPRMs[currCol] + "-" +
			// alignedPeaks[currRow] + ")=" + (this.SeedPRMs[currCol] -
			// alignedPeaks[currRow]));
			// System.out.println("shift: " + (this.SeedPRMs[currCol] -
			// alignedPeaks[currRow]));
			qSum += this.SeedPRMs[currCol] - alignedPeaks[currRow];
			L += 1;
			// prevRowMass = specPRMs[currRow];
			// prevColMass = specPRMs[currCol];
			int currCell = prevCell[currRow][currCol];
			currRow = currCell / prevCell[0].length;
			currCol = currCell % prevCell[0].length;
		}

		if (L == 0) {
			System.err.println("ERROR in Seed::GetMassShift");
			System.err.println("There were 0 aligned peaks!");
			System.err.println("Best Row: " + bestRow + ", Best Col: "
					+ bestCol);
			System.err.println("prevCell[" + bestRow + "][" + bestCol + "] = "
					+ prevCell[bestRow][bestCol]);
			System.err.println("Best Score: " + bestScore);
			System.err.println("A: " + Utils.JoinIntArray(alignedPeaks, " "));
			System.err.println("B: " + Utils.JoinIntArray(SeedPRMs, " "));
			System.exit(-1);
		}
		return qSum / L;
	}

}
