package antibody;

/**
 * Currently only represents a Gaussian, but maybe someday will represent
 * something else, like with Isotopic distribution...
 * 
 * The type does no scaling, so make sure that the scale that you create the
 * distribution with is used consistently.
 * 
 * @author natalie
 * 
 */
public class DistributionType {

	private double Mean;
	private double Variance;
	private double[] Samples;

	/**
	 * Creates a Gaussian distribution with specified mean and variance
	 * 
	 * @param Mean
	 *            The mean of the distribution
	 * @param Variance
	 *            The variance of the distribution
	 */
	public DistributionType(double Mean, double Variance) {
		this.Mean = Mean;
		this.Variance = Variance;
		this.Samples = new double[1];
		this.Samples[0] = Mean;
	}

	public DistributionType(double[] Samples, double Variance) {

		this.Samples = Samples;
		double Sum = 0.0;
		for (int i = 0; i < this.Samples.length; ++i) {
			Sum += this.Samples[i];
		}
		this.Mean = Sum / (this.Samples.length);
		this.Variance = Variance;

	}

	public DistributionType(double[] Samples) {
		double TotalSum = 0;
		for (int i = 0; i < Samples.length; ++i) {
			TotalSum += Samples[i];
		}
		this.Mean = TotalSum / (Samples.length);
		this.Variance = 0.0;
		for (int i = 0; i < Samples.length; ++i) {
			this.Variance += Math.pow(Samples[i] - this.Mean, 2);
		}
		this.Variance = this.Variance / (Samples.length - 1);
		this.Samples = Samples;

	}

	/*
	 * public DistributionType(int PathLength, int[] ScaledPeaks, double[]
	 * PeakScores) { // Find a distributions of all DeNovoPaths of length
	 * PathLength double[] Samples = DeNovo.FindAllPathScores(PathLength,
	 * ScaledPeaks, PeakScores, 1000); // System.out.println("Found " +
	 * Samples.length + // " De novo reconstructions of length " + PathLength);
	 * double TotalSum = 0; for (int i = 0; i < Samples.length; ++i) { TotalSum
	 * += Samples[i]; } this.Mean = TotalSum / (Samples.length); this.Variance =
	 * 0.0; for (int i = 0; i < Samples.length; ++i) { this.Variance +=
	 * Math.pow(Samples[i] - this.Mean, 2); } this.Variance = this.Variance /
	 * (Samples.length - 1); this.Samples = Samples;
	 * 
	 * }
	 */

	/**
	 * We can also create a distribution for the sum of NumToSum numbers from
	 * our list of values
	 * 
	 * @param NumToSum
	 *            The number of items in Values to sum as the random variable
	 * @param Values
	 *            The items in the sample
	 */
	/*
	 * public DistributionType(int NumToSum, double[] Values) { boolean
	 * LocalDebug = false;
	 * 
	 * // boolean[] SumThese = new boolean[Values.length]; // double[]
	 * TrueValues = new double[Values.length/2]; double[] TrueValues = Values;
	 * // for(int i = 0; i < TrueValues.length; ++i) // { // TrueValues[i] =
	 * Values[2*i]; // } double TotalSums = Arithmetic.binomial((double)
	 * (TrueValues.length), (long) NumToSum); //
	 * System.out.println("We want this many subsets: " + TotalSums);
	 */
	/*
	 * for(int i = 0; i < SumThese.length; ++i) { SumThese[i] = false; if(i >=
	 * SumThese.length - NumToSum) SumThese[i] = true; if(LocalDebug)
	 * System.out.println("Bit " + i + " = " + SumThese[i]); }
	 */
	/*
	 * // ArrayList Subsets = new ArrayList(); ArrayList Subsets =
	 * AntibodyUtils.FindAllSubsets(TrueValues.length, NumToSum); //
	 * System.out.println("We only got: " + Subsets.size()); double TotalSum =
	 * 0.0; double[] Sums = new double[Subsets.size()];
	 * 
	 * for (int i = 0; i < Subsets.size(); ++i) { int[] List = (int[])
	 * (Subsets.get(i)); String CurrString = ""; double Sum = 0.0; for (int j =
	 * 0; j < List.length; ++j) { if (LocalDebug) CurrString += List[j] + " ";
	 * Sum += TrueValues[List[j]]; } Sums[i] = Sum; TotalSum += Sum; if
	 * (LocalDebug) System.out.println(CurrString); } this.Mean = TotalSum /
	 * (Subsets.size()); this.Variance = 0.0; for (int i = 0; i < Sums.length;
	 * ++i) { this.Variance += Math.pow(Sums[i] - this.Mean, 2); } this.Variance
	 * = this.Variance / (Sums.length - 1); this.Samples = Sums; }
	 */
	/**
	 * If no parameters are specified, a unit normal is produced with mean 0,
	 * variance = 1.
	 * 
	 */
	public DistributionType() {
		this.Mean = 0;
		this.Variance = 1;
	}

	public double GetMean() {
		return this.Mean;
	}

	public double GetVariance() {
		return this.Variance;
	}

	public void AddNewSample(double NewValue) {

		this.Mean = (this.Mean * Samples.length / (Samples.length + 1))
				+ NewValue / (Samples.length + 1);
		double[] NewSamples = new double[this.Samples.length + 1];
		for (int i = 0; i < this.Samples.length; ++i) {
			NewSamples[i] = this.Samples[i];
		}
		NewSamples[this.Samples.length] = NewValue;
		this.Samples = NewSamples;
	}

	/**
	 * This method assumes that the value distribution can support a 50 unit
	 * window for PDF calculation
	 * 
	 * @param Value
	 *            The observed value
	 * @return The approximate probability of observing this value with this
	 *         distribution.
	 */
	public double ApproximatePDF(double Value) {

		double Ret = (1 / (Math.sqrt(2 * Math.PI * this.Variance)))
				* Math.exp(-1 * (Math.pow(Value - this.Mean, 2))
						/ (2 * this.Variance));

		return Ret;
	}

	/**
	 * Returns the Pvalue, the probability of getting a higher score
	 * 
	 * @param Value
	 *            The Value cutoff
	 * @return Returns the PValue
	 */
	/*
	 * public double PValue(double Value) { return 1 -
	 * Probability.normal(this.Mean, this.Variance, Value);
	 * 
	 * }
	 */

	public void DebugPrint() {
		System.out.println("----------DistributionType----------");
		System.out
				.println("Mean: " + this.Mean + " Variance: " + this.Variance);
		System.out.println("");
	}

	public static void main(String[] args) {
		double[] Values = { 2000, 1900, 2000, 2200, 2600 };
		DistributionType D = new DistributionType(Values);
		D.ApproximatePDF(2.5);
		D.DebugPrint();

	}

}
