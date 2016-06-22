package antibody;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

import basicUtils.InspectAnnotation;
import basicUtils.Utils;
import errorUtils.ErrorThrower;

/*import org.jfree.chart.ChartFactory;
 import org.jfree.chart.JFreeChart;
 import org.jfree.chart.plot.PlotOrientation;
 import org.jfree.data.statistics.HistogramDataset;
 */
/**
 * The PRMSpectrum class represents a prefix residue mass spectrum
 * 
 * @author natalie
 * @version 12.16.2008
 */
public class PRMSpectrum {

	public static boolean Debug = false;

	public int[] PRMs;
	public double[] PRMScores;

	public double MeanScore;
	public double MaxScore = Double.MIN_VALUE;
	public double MinScore = Double.MAX_VALUE;
	public double ScoreVariance;

	public int PrecursorMass; //M of precursor
	public int PrecursorCharge; 

	/**
	 * Spectrum index(es) for this PRM spectrum. This should always be length 1
	 */
	public int[] SpecIndex;

	/**
	 * Scan number(s) for this PRM spectrum. Scan numbers are 1-based
	 */
	public int[] ScanNumber;
	public String FileName;

	/*
	 * public static final String UsageInfo = "java antibody.PRMSpectrum\n" +
	 * "Creates a histogram of the PRM score distribution for a given file or directory\n"
	 * + "Options:\n" +
	 * " -r [DIR/FILE] Directory or File containing PepNovo produced PRM spectra\n"
	 * + " -w [FILE] File to write histogram file (use .png as extension)\n";
	 */
	public PRMSpectrum() {
		this.PrecursorMass = -1 * AntibodyUtils.MASS_SCALE;
		this.PrecursorCharge = 0;
	}

	public PRMSpectrum(String Seq) {
		this.PRMs = AntibodyUtils.GetSeqPRMScaled(Seq);
		this.PRMScores = new double[this.PRMs.length];
		for (int i = 0; i < this.PRMScores.length; ++i)
			PRMScores[i] = AntibodyUtils.MIN_PEAK_GROUP_SCORE;
	}

	public PRMSpectrum(double Mass, int Charge, String[] Data) {
		this.PrecursorMass = (int) Mass * AntibodyUtils.MASS_SCALE;
		this.PrecursorCharge = Charge;
		this.SetData(Data);

	}

	public static PRMSpectrum MergePRMsSimple(PRMSpectrum A, PRMSpectrum B,
			boolean Debug, double Tolerance) {
		int aIndex = 0;
		int bIndex = 0;

		if (Debug) {
			System.out.println("Precursors: " + A.PrecursorMass + " and "
					+ B.PrecursorMass);
			System.out.print("A:");
			for (int i = 0; i < A.PRMs.length; ++i)
				System.out.print(" " + A.PRMs[i]);
			System.out.print("\nB:");
			for (int i = 0; i < B.PRMs.length; ++i)
				System.out.print(" " + B.PRMs[i]);
			System.out.print("\n");
		}

		ArrayList<String> TempData = new ArrayList<String>();

		while (aIndex < A.PRMs.length && bIndex < B.PRMs.length) {
			if (AntibodyUtils.EquivalentScaled(A.PRMs[aIndex], B.PRMs[bIndex],
					Tolerance)) {
				int AveMass = (A.PRMs[aIndex] + B.PRMs[bIndex]) / 2;
				double AveScore = (A.PRMScores[aIndex] + B.PRMScores[bIndex]) / 2;
				TempData.add(AveMass + " " + AveScore);
				aIndex++;
				bIndex++;
			} else if (A.PRMs[aIndex] < B.PRMs[bIndex]) {

				TempData.add(A.PRMs[aIndex] + " " + A.PRMScores[aIndex]);
				aIndex++;
			} else {
				TempData.add(B.PRMs[bIndex] + " " + B.PRMScores[bIndex]);
				bIndex++;
			}

		}
		while (aIndex < A.PRMs.length) {
			TempData.add(A.PRMs[aIndex] + " " + A.PRMScores[aIndex]);
			aIndex++;
		}
		while (bIndex < B.PRMs.length) {
			TempData.add(B.PRMs[bIndex] + " " + B.PRMScores[bIndex]);
			bIndex++;
		}

		String[] Data = new String[TempData.size()];
		for (int i = 0; i < Data.length; ++i) {
			Data[i] = (String) (TempData.get(i));
			if (Debug)
				System.out.println(Data[i]);
		}
		if (Debug)
			Utils.WaitForEnter();

		double PrecursorMass = ((double) (A.PrecursorMass + B.PrecursorMass) / 2)
				/ AntibodyUtils.MASS_SCALE;
		PRMSpectrum Ret = new PRMSpectrum(PrecursorMass, A.PrecursorCharge,
				Data);
		Ret.FileName = A.FileName;
		Ret.SpecIndex = A.SpecIndex;

		return Ret;

	}

	/**
	 * Sets the PRMs and PRMScores for this Spectrum
	 * 
	 * @param Data
	 *            An array of strings of the form 'PRMMass PRMScore'
	 */
	public boolean SetData(Object[] Data) {
		if (Data == null || Data.length == 0) {
			this.PRMs = new int[0];
			this.PRMScores = new double[0];
			return true;
		}

		this.PRMs = new int[Data.length];
		this.PRMScores = new double[Data.length];
		double ScoreSum = 0.0;

		for (int i = 0; i < Data.length; ++i) {
			String[] Bits = ((String) Data[i]).split(" ");
			try {
				this.PRMs[i] = (int) ((new Double(Bits[0])).doubleValue() * AntibodyUtils.MASS_SCALE);
				this.PRMScores[i] = (new Double(Bits[1])).doubleValue();
				if (this.PRMScores[i] > this.MaxScore)
					this.MaxScore = this.PRMScores[i];
				if (this.PRMScores[i] < this.MinScore)
					this.MinScore = this.PRMScores[i];
			} catch (Exception E) {
				// System.err.println("ERROR: Errant peak line: " + Data[i]);
				return false;
			}
			ScoreSum += this.PRMScores[i];
		}
		this.MeanScore = ScoreSum / (Data.length);

		double SampleVariance = 0.0;
		for (int i = 0; i < this.PRMScores.length; ++i) {
			SampleVariance += Math.pow(this.PRMScores[i] - this.MeanScore, 2);
		}

		// Hack tomake sure we have a 0
		if (this.PRMs[0] != 0) {
			// System.out.println("addining initial PRM of 0 before PRM of " +
			// this.PRMs[0]);
			int[] NewPRMs = new int[this.PRMs.length + 1];
			double[] NewPRMScores = new double[this.PRMs.length + 1];
			NewPRMs[0] = 0;
			NewPRMScores[0] = 15;
			for (int i = 0; i < this.PRMs.length; ++i) {
				NewPRMs[i + 1] = this.PRMs[i];
				NewPRMScores[i + 1] = this.PRMScores[i];
			}
			this.PRMs = NewPRMs;
			this.PRMScores = NewPRMScores;
		}

		this.ScoreVariance = SampleVariance / (this.PRMScores.length - 1);
		return true;
	}

	/**
	 * Expects a PRM file produced by Pepnovo, loads them into an array sorted
	 * by Parent mass. The scan number is the order of the PRMs in the file.
	 * 
	 * @param FileName
	 *            PRM file name
	 * @return Returns an array of PRMSpectrum
	 * @throws IOException
	 */
	public static PRMSpectrum[] LoadAllPRMSpectrum(String FileName) {

		boolean localDebug = false;
		// This is sort of hacky, but if the extension is .mgf, then use a
		// different loader
		if (Utils.GetFileExtension(FileName).equalsIgnoreCase(".mgf")) {
			return PRMSpectrum.LoadAllPRMSpectrumMGF(FileName);
		}

		// DataInputStream reader = new DataInputStream(new
		// FileInputStream(FileName));
		BufferedReader reader = Utils.openBufferedReader(FileName);
		String Line;
		ArrayList AllSpectra = new ArrayList();
		PRMSpectrum nextSpectrum = null;
		ArrayList currData = null;
		int NoPeakCount = 0;
		int spectrumIndex = 0;
		// int ScanNumber = 0;
		if (localDebug)
			System.out.println("Loading PRMS from " + FileName);

		while ((Line = Utils.readNextLine(reader, FileName)) != null) {
			// Line = Line.trim();
			if (Line.length() > 0 && Line.charAt(0) == '#') {
				continue;
			}

			if (Line.contains(">>")) // Line 1 of a new spectrum
			{
				if (localDebug)
					System.out.println("New Spectrum: " + Line);
				// String[] Bits = Line.split(" ");
				nextSpectrum = new PRMSpectrum();

				nextSpectrum.FileName = FileName;
				// nextSpectrum.ScanNumber = (new Integer(Bits[2])).intValue();
				// nextSpectrum.ScanNumber = Integer.parseInt(Bits[2]);
				nextSpectrum.SpecIndex = new int[1];
				nextSpectrum.SpecIndex[0] = spectrumIndex;
				spectrumIndex += 1;

				currData = new ArrayList();
			} else if (Line.contains("Charge")) // Line 2 of a new spectrum
			{
				if (Debug)
					System.out.println("Charge line: " + Line);
				String[] Bits = Line.split(" ");
				nextSpectrum.PrecursorCharge = (new Integer(Bits[1]))
						.intValue();
				nextSpectrum.PrecursorMass = AntibodyUtils
						.ComputeMassFromMHScaled((new Double(
								Bits[Bits.length - 1])).doubleValue());
			} else if (Line.length() == 0 && nextSpectrum != null) // End of a
			// spectrum
			{
				if (localDebug)
					System.out.println("Finished a spectrum!!");
				if (currData.size() > 0) {
					if (nextSpectrum.SetData(currData.toArray())) {
						AllSpectra.add(nextSpectrum);
					}
				} else {
					nextSpectrum.SetData(null);
					AllSpectra.add(nextSpectrum);
					// System.out.println("Found a spectrum with no Peaks!! " +
					// nextSpectrum.FileName + ":" + nextSpectrum.ScanNumber);
					NoPeakCount += 1;
				}
				nextSpectrum = null;
			} else if (nextSpectrum != null) // Line 3 and beyond if the
			// spectrum is initialized
			{
				// if(Debug)
				// System.out.println("This is just new data: " + Line);
				currData.add(Line);
			}
		}
		if (nextSpectrum != null) {
			if (nextSpectrum.SetData(currData.toArray())) {
				boolean Found = false;
				for (int k = 0; k < AllSpectra.size(); ++k) {
					PRMSpectrum PrevSpec = (PRMSpectrum) (AllSpectra.get(k));
					if (nextSpectrum.PrecursorMass < PrevSpec.PrecursorMass) {
						AllSpectra.add(k, nextSpectrum);
						Found = true;
						break;
					}
				}
				if (!Found)
					AllSpectra.add(nextSpectrum);
			}
			nextSpectrum = null;
		}

		PRMSpectrum[] Ret = new PRMSpectrum[AllSpectra.size()];
		for (int i = 0; i < Ret.length; ++i) {
			Ret[i] = (PRMSpectrum) (AllSpectra.get(i));
		}

		// System.out.println("File: " + FileName + " Loaded " + Ret.length +
		// ", skipped " + NoPeakCount);
		Utils.closeBufferedReader(reader, FileName);
		return Ret;

	}

	private static PRMSpectrum[] LoadAllPRMSpectrumMGF(String FileName) {

		boolean localDebug = false;
		BufferedReader reader = Utils.openBufferedReader(FileName);
		// Debug = true;
		String Line = Utils.readNextLine(reader, FileName);
		ArrayList AllSpectra = new ArrayList();
		PRMSpectrum nextSpectrum = null;
		ArrayList currData = null;
		double pepMass = 0;
		int spectrumIndex = 0;
		// int ScanNumber = 0;
		if (localDebug)
			System.out.println("Loading PRMS from " + FileName);

		while (Line != null) {
			// Line = Line.trim();
			Line = Line.toLowerCase();
			if (Line.length() > 0 && Line.charAt(0) == '#') {
				Line = Utils.readNextLine(reader, FileName);
				continue;
			}

			if (Line.contains("begin ions")) // Line 1 of a new spectrum
			{
				if (localDebug)
					System.out.println("New Spectrum: " + spectrumIndex);
				// String[] Bits = Line.split(" ");
				nextSpectrum = new PRMSpectrum();
				pepMass = 0;
				nextSpectrum.FileName = FileName;
				// nextSpectrum.ScanNumber = (new Integer(Bits[2])).intValue();
				// nextSpectrum.ScanNumber = Integer.parseInt(Bits[2]);
				nextSpectrum.SpecIndex = new int[1];
				nextSpectrum.SpecIndex[0] = spectrumIndex;
				spectrumIndex += 1;

				currData = new ArrayList();
			} else if (Line.contains("charge")) // Line 2 of a new spectrum
			{
				if (localDebug)
					System.out.println("Parse charge line: " + Line);
				String[] Bits = Line.split("=");
				nextSpectrum.PrecursorCharge = PRMSpectrum
						.ParseMGFCharge(Bits[1]);
			} else if (Line.contains("pepmass")) {
				if (localDebug)
					System.out.println("Parse pepmass line: " + Line);
				String[] Bits = Line.split("=");
				try {
					pepMass = Double.parseDouble(Bits[1]);
				} catch (Exception E) {
					ErrorThrower.ThrowWarning(20, Bits[1]);
					pepMass = 0.0;
				}
			} else if (Line.contains("scans")) {
				if (localDebug)
					System.out.println("Parse pepmass line: " + Line);
				String[] Bits = Line.split("=");

				int scanNum = nextSpectrum.SpecIndex[0] + 1;

				try {
					scanNum = Integer.parseInt(Bits[1]);
				} catch (Exception E) {

				}
				nextSpectrum.ScanNumber = new int[1];
				nextSpectrum.ScanNumber[0] = scanNum;

			} else if (Line.contains("end ions")) // End of a
			// spectrum
			{
				if (localDebug)
					System.out.println("Finished a spectrum!!");
				if (currData.size() > 0) {
					if (nextSpectrum.SetData(currData.toArray())) {
						AllSpectra.add(nextSpectrum);
						nextSpectrum.PrecursorMass = (int) (Utils.GetMFromMZ(
								pepMass, nextSpectrum.PrecursorCharge) * AntibodyUtils.MASS_SCALE);
					}
				} else {
					nextSpectrum.SetData(null);
					//nextSpectrum.PrecursorMass = 0;
					nextSpectrum.PrecursorMass = (int) (Utils.GetMFromMZ(
							pepMass, 1) * AntibodyUtils.MASS_SCALE);
					
					AllSpectra.add(nextSpectrum);
					// System.out.println("Found a spectrum with no Peaks!! " +
					// nextSpectrum.FileName + ":" + nextSpectrum.ScanNumber);
				}

				nextSpectrum = null;

			} else if (Line.contains("=")) // THere might be a lot of other
											// attributes we don't care about,
											// so just skip them.

			{
				if (localDebug)
					System.out.println("Skipping line: " + Line);
			} else if (nextSpectrum != null) // Line 3 and beyond if the
			// spectrum is initialized
			{
				if (localDebug)
					System.out.println("This is just new data: " + Line);
				currData.add(Line);
			}
			Line = Utils.readNextLine(reader, FileName);
		}
		if (nextSpectrum != null) {
			if (nextSpectrum.SetData(currData.toArray())) {
				boolean Found = false;
				for (int k = 0; k < AllSpectra.size(); ++k) {
					PRMSpectrum PrevSpec = (PRMSpectrum) (AllSpectra.get(k));
					if (nextSpectrum.PrecursorMass < PrevSpec.PrecursorMass) {
						AllSpectra.add(k, nextSpectrum);
						Found = true;
						break;
					}
				}
				if (!Found)
					AllSpectra.add(nextSpectrum);
			}
			nextSpectrum = null;
		}

		PRMSpectrum[] Ret = new PRMSpectrum[AllSpectra.size()];
		for (int i = 0; i < Ret.length; ++i) {
			Ret[i] = (PRMSpectrum) (AllSpectra.get(i));
		}

		// System.out.println("File: " + FileName + " Loaded " + Ret.length +
		// ", skipped " + NoPeakCount);
		return Ret;
	}

	/**
	 * MGF files often have stupid looking charges like 2+ or 3+. This parses
	 * that string to return an integer
	 * 
	 * @param string
	 * @return
	 */
	private static int ParseMGFCharge(String string) {
		boolean plusCharge = true;
		if (string.charAt(string.length() - 1) == '+')
			string = string.substring(0, string.length() - 1);
		else if (string.charAt(string.length() - 1) == '-') {
			string = string.substring(0, string.length() - 1);
			plusCharge = false;
		}
		int val = Integer.parseInt(string);
		if (!plusCharge)
			val = val * -1;
		return val;
	}

	/**
	 * Adds a collection of scaled, sorted PRM masses to the spectrm with the
	 * given peak score. If the mass already exists in this PRM spectrum, then
	 * nothing is added
	 * 
	 * @param pepPRMs
	 * @param mINPEAKGROUPSCORE
	 */
	public void addPRMs(int[] pepPRMs, double prmScore, double tol) {

		ArrayList newPRMs = new ArrayList();
		ArrayList newPRMScores = new ArrayList();

		int oldPRMIndex = 0;
		int newPRMIndex = 0;

		int ScoreSum = 0;

		while (oldPRMIndex < this.PRMs.length && newPRMIndex < pepPRMs.length) {
			// If these masses match, advance both pointers
			if (AntibodyUtils.EquivalentScaled(this.PRMs[oldPRMIndex],
					pepPRMs[newPRMIndex], tol)) {
				newPRMs.add(new Integer(this.PRMs[oldPRMIndex]));
				newPRMScores.add(new Double(this.PRMScores[oldPRMIndex]));
				ScoreSum += this.PRMScores[oldPRMIndex];
				oldPRMIndex += 1;
				newPRMIndex += 1;
			} else if (this.PRMs[oldPRMIndex] < pepPRMs[newPRMIndex]) {
				newPRMs.add(new Integer(this.PRMs[oldPRMIndex]));
				newPRMScores.add(new Double(this.PRMScores[oldPRMIndex]));
				ScoreSum += this.PRMScores[oldPRMIndex];
				oldPRMIndex += 1;
			} else if (this.PRMs[oldPRMIndex] > pepPRMs[newPRMIndex]) {
				newPRMs.add(new Integer(pepPRMs[newPRMIndex]));
				newPRMScores.add(new Double(prmScore));
				ScoreSum += prmScore;
				if (prmScore < this.MinScore)
					this.MinScore = prmScore;
				if (prmScore > this.MaxScore)
					this.MaxScore = prmScore;
				newPRMIndex += 1;
			}
		}
		while (oldPRMIndex < this.PRMs.length) {
			newPRMs.add(new Integer(this.PRMs[oldPRMIndex]));
			newPRMScores.add(new Double(this.PRMScores[oldPRMIndex]));
			ScoreSum += this.PRMScores[oldPRMIndex];
			oldPRMIndex += 1;
		}
		while (newPRMIndex < pepPRMs.length) {
			newPRMs.add(new Integer(pepPRMs[newPRMIndex]));
			newPRMScores.add(new Double(prmScore));
			ScoreSum += prmScore;
			if (prmScore < this.MinScore)
				this.MinScore = prmScore;
			if (prmScore > this.MaxScore)
				this.MaxScore = prmScore;
			newPRMIndex += 1;
		}

		this.MeanScore = ScoreSum / newPRMs.size();
		this.PRMs = Utils.ConvertArraylistToIntArray(newPRMs);
		this.PRMScores = Utils.ConvertArraylistToDoubleArray(newPRMScores);

		double SampleVariance = 0.0;
		for (int i = 0; i < this.PRMScores.length; ++i) {
			SampleVariance += Math.pow(this.PRMScores[i] - this.MeanScore, 2);
		}
		this.ScoreVariance = SampleVariance / (this.PRMScores.length - 1);

	}

	/**
	 * Sam as addPRMs, but removes unused guys
	 * 
	 * @param pepPRMs
	 * @param mINPEAKGROUPSCORE
	 */
	public void setPRMs(int[] pepPRMs, double prmScore, double tol) {
		ArrayList newPRMs = new ArrayList();
		ArrayList newPRMScores = new ArrayList();

		int oldPRMIndex = 0;
		int newPRMIndex = 0;

		int ScoreSum = 0;

		while (oldPRMIndex < this.PRMs.length && newPRMIndex < pepPRMs.length) {
			// If these masses match, advance both pointers
			if (AntibodyUtils.EquivalentScaled(this.PRMs[oldPRMIndex],
					pepPRMs[newPRMIndex], tol)) {
				newPRMs.add(new Integer(this.PRMs[oldPRMIndex]));
				newPRMScores.add(new Double(this.PRMScores[oldPRMIndex]));
				ScoreSum += this.PRMScores[oldPRMIndex];
				oldPRMIndex += 1;
				newPRMIndex += 1;
			} else if (this.PRMs[oldPRMIndex] < pepPRMs[newPRMIndex]) {
				// newPRMs.add(new Integer(this.PRMs[oldPRMIndex]));
				// newPRMScores.add(new Double(this.PRMScores[oldPRMIndex]));
				// ScoreSum += this.PRMScores[oldPRMIndex];
				oldPRMIndex += 1;
			} else if (this.PRMs[oldPRMIndex] > pepPRMs[newPRMIndex]) {
				newPRMs.add(new Integer(pepPRMs[newPRMIndex]));
				newPRMScores.add(new Double(prmScore));
				ScoreSum += prmScore;
				if (prmScore < this.MinScore)
					this.MinScore = prmScore;
				if (prmScore > this.MaxScore)
					this.MaxScore = prmScore;
				newPRMIndex += 1;
			}
		}
		while (oldPRMIndex < this.PRMs.length) {
			// newPRMs.add(new Integer(this.PRMs[oldPRMIndex]));
			// newPRMScores.add(new Double(this.PRMScores[oldPRMIndex]));
			// ScoreSum += this.PRMScores[oldPRMIndex];
			oldPRMIndex += 1;
		}
		while (newPRMIndex < pepPRMs.length) {
			newPRMs.add(new Integer(pepPRMs[newPRMIndex]));
			newPRMScores.add(new Double(prmScore));
			ScoreSum += prmScore;
			if (prmScore < this.MinScore)
				this.MinScore = prmScore;
			if (prmScore > this.MaxScore)
				this.MaxScore = prmScore;
			newPRMIndex += 1;
		}

		this.MeanScore = ScoreSum / newPRMs.size();
		this.PRMs = Utils.ConvertArraylistToIntArray(newPRMs);
		this.PRMScores = Utils.ConvertArraylistToDoubleArray(newPRMScores);

		double SampleVariance = 0.0;
		for (int i = 0; i < this.PRMScores.length; ++i) {
			SampleVariance += Math.pow(this.PRMScores[i] - this.MeanScore, 2);
		}
		this.ScoreVariance = SampleVariance / (this.PRMScores.length - 1);

	}

	public static PRMSpectrum[] createFakePRMSpectra(InspectAnnotation[] anns) {

		PRMSpectrum[] ret = new PRMSpectrum[anns.length];
		for (int i = 0; i < anns.length; ++i) {
			String pep = anns[i].Annotation;
			ret[i] = new PRMSpectrum(pep);
		}

		return ret;
	}

	/*
	 * public static void CreatePRMScoreHisto(String[] InputFileNames, String
	 * OutputFileName) { HistogramDataset ScoreHisto = new HistogramDataset();
	 * 
	 * ArrayList ScoresList = new ArrayList(); PRMSpectrum[] Spectra = null;
	 * double MaxScore = Double.MIN_VALUE; double MinScore = Double.MAX_VALUE;
	 * for (int i = 0; i < InputFileNames.length; ++i) {
	 * System.out.println("Loading " + InputFileNames[i] + "..."); try { Spectra
	 * = PRMSpectrum.LoadAllPRMSpectrum(InputFileNames[i]); } catch (IOException
	 * e) { e.printStackTrace(); } for (int j = 0; j < Spectra.length; ++j) {
	 * for (int k = 0; k < Spectra[j].PRMScores.length; ++k) { if
	 * (Spectra[j].PRMScores[k] > MaxScore) MaxScore = Spectra[j].PRMScores[k];
	 * if (Spectra[j].PRMScores[k] < MinScore) MinScore =
	 * Spectra[j].PRMScores[k]; ScoresList.add(new
	 * Double(Spectra[j].PRMScores[k])); } } }
	 * 
	 * int NumScoreBins = Math.max(1, (int) ((MaxScore - MinScore) / 0.1));
	 * double[] Scores = AntibodyUtils.ConvertDoubleArrayList(ScoresList);
	 * ScoreHisto.addSeries("PRM Scores", Scores, NumScoreBins); JFreeChart
	 * Chart = ChartFactory.createHistogram("Score Histogram(" +
	 * ScoresList.size() + " )", "Score", "Counts", ScoreHisto,
	 * PlotOrientation.VERTICAL, true, false, false);
	 * 
	 * BufferedImage Image = Chart.createBufferedImage(1000, 1000); try {
	 * ImageIO.write(Image, AntibodyUtils.GetFileExtension(OutputFileName), new
	 * File( OutputFileName)); } catch (Exception E) {
	 * System.err.println(E.getMessage()); return; }
	 * 
	 * }
	 */

	/*
	 * public static void main(String[] args) throws IOException { String[]
	 * options = { "-r", "-w" }; boolean[] values = { true, true }; Hashtable
	 * Options = AntibodyUtils.ParseCommandLine(args, options, values); if
	 * (!Options.containsKey("-r") || !Options.containsKey("-w")) { System.err
	 * .println("ERROR: Must specify an input and output file name");
	 * System.err.println(UsageInfo); return; }
	 * 
	 * String InputDir = (String) (Options.get("-r")); File InputFile = new
	 * File(InputDir); if (!InputFile.exists()) {
	 * System.err.println("ERROR: Input file '" + InputDir +
	 * "' does not exists!"); return; } String[] InputFileNames = null;
	 * 
	 * if (InputFile.isDirectory()) { String[] Files = InputFile.list();
	 * InputFileNames = new String[Files.length]; for (int i = 0; i <
	 * Files.length; ++i) { InputFileNames[i] = InputFile.getAbsolutePath() +
	 * Files[i]; } } else { InputFileNames = new String[1]; InputFileNames[0] =
	 * InputDir; }
	 * 
	 * String OutputFileName = (String) (Options.get("-w")); File OutputFile =
	 * new File(OutputFileName); if (!OutputFile.exists()) OutputFile.mkdir();
	 * 
	 * //PRMSpectrum.CreatePRMScoreHisto(InputFileNames, OutputFileName);
	 * 
	 * }
	 */

}
