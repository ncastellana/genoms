package antibody;

import errorUtils.ErrorThrower;

public class AlignmentUtils {

	public enum alignEndType {doubleFixed, singleFixed, any};
	
	public class Alignment {
		
		private double[][] alignmentScores;
		private int[][] prevCell;
		
		private double[] alignedPeaksA;
		private double[] alignedScoresA;
		
		private double[] alignedPeaksB;
		private double[] alignedScoresB;
		
		private alignEndType startType;
		private boolean modsAllowed;
		
		public Alignment(double[][] scoreTable, int[][] prevCell, 
				double[] peaksA, double[] scoresA, double[] peaksB, 
				double[] scoresB, alignEndType start, boolean modsAllowed) {
			
			this.alignmentScores = scoreTable;
			this.prevCell = prevCell;
			this.alignedPeaksA = peaksA;
			this.alignedScoresA = scoresA;
			this.alignedPeaksB = peaksB;
			this.alignedScoresB = scoresB;
			
			this.startType = start;
			this.modsAllowed = modsAllowed;
			
			
		}
		public Alignment(double[][] scoreTable, int[][] prevCell, 
				alignEndType fixedStart, boolean modsAllowed) {
			
			this.alignmentScores = scoreTable;
			this.prevCell = prevCell;
			
			this.startType = fixedStart;
			this.modsAllowed = modsAllowed;
			
			
		}
		
		public double getBestScore(alignEndType fixedEnd) {
			if(this.alignmentScores == null || this.alignmentScores.length < 1 || this.alignmentScores[0].length < 1) {
				ErrorThrower.ThrowWarningCustom(ErrorThrower.CUSTOM_WARNING, "No alignment table set for this Alignment object!");
				return -1;
			}
			if(fixedEnd == alignEndType.doubleFixed)
				return this.alignmentScores[alignmentScores.length-1][alignmentScores[0].length-1];
			
			else if(fixedEnd == alignEndType.singleFixed) {
				double maxScore = alignmentScores[0][0];
				//Ends anywhere in the first spectrum
				for(int i = 0; i < alignmentScores.length; ++i) {
					maxScore = Math.max(maxScore, this.alignmentScores[i][this.alignmentScores[i].length-1]);
				}
				//Ends anywhere in the second spectrum
				for(int i = 0; i < alignmentScores[0].length; ++i) {
					maxScore = Math.max(maxScore, this.alignmentScores[this.alignmentScores.length-1][i]);
				}
				return maxScore;
			}
			
			double maxScore = alignmentScores[0][0];
			for(int i = 0; i < alignmentScores.length; ++i) {
				for(int j = 0; j < alignmentScores[0].length; ++j) {
					maxScore = Math.max(maxScore, this.alignmentScores[i][j]);
				}
			}
			return maxScore;
		}
		
	}
	/**
	 * Gets the score of the best (highest scoring) alignment between two spectra.  The score is
	 * the sum of the aligned peaks.  Alignment allows gaps of any size but no mods.  Alignment can start anywhere 
	 * in either spectrum.
	 * @param A PRMSpectrum object to be aligned
	 * @param B PRMSpectrum object to be aligned
	 * @param Tolerance PRM tolerance in Daltons
	 * @return
	 */
	public double getAlignmentScore(PRMSpectrum A, PRMSpectrum B,
			double Tolerance) {
		double Score = 0.0;

		
		Alignment aln = this.getAlignmentAnyStart(A,B,Tolerance);

		

		return aln.getBestScore(alignEndType.singleFixed);

	}
	/**
	 * Gets the score of the best (highest scoring) alignment between two spectra.  The score is
	 * the sum of the aligned peaks.  Alignment allows gaps of any size but no mods. Aligment must start at the beginning of
	 * each spectrum and end at the final peak of each spectrum.
	 * @param A PRMSpectrum object to be aligned
	 * @param B PRMSpectrum object to be aligned
	 * @param Tolerance PRM tolerance in Daltons
	 * @return
	 */
	public double getAlignmentScoreFixedEnds(PRMSpectrum A, PRMSpectrum B,
			double Tolerance) {
		double Score = 0.0;

		Alignment aln = this.getAlignmentScoreTableFixedStart(A,B,Tolerance);

		return aln.getBestScore(alignEndType.doubleFixed);

	}
	
	/**
	 * Returns the score table for the alignment of two PRM spectra with a given PRM tolerance.  Any number of missed PRMs is allowed
	 * but no modifications.  Alignment is expected to be between overlapping peptides with identical sequences in the overlap region.
	 * @param A A PRMSpectrum object with N peaks
	 * @param B A PRMSpectrum object with M peaks
	 * @param Tolerance
	 * @return Returns a NxM table where position (i,j) is the highest score of an alignment where peak i of spectrum A is matched to peak j of spectrum B.
	 */
	private Alignment getAlignmentAnyStart(PRMSpectrum A, PRMSpectrum B, double Tolerance) {
		double[][] Alignment = new double[A.PRMs.length][B.PRMs.length];
		int[][] prevCell = new int[A.PRMs.length][B.PRMs.length];
		
		for (int a = 0; a < A.PRMs.length; ++a) {
			double AScore = A.PRMScores[a];
			int AMass = A.PRMs[a];
			for (int b = 0; b < B.PRMs.length; ++b) {
				double BScore = B.PRMScores[b];
				int BMass = B.PRMs[b];

				double MaxScore = AScore + BScore;
				int[] bestPrev = {-1,-1};
				for (int prevA = 0; prevA < a; ++prevA) {
					int prevAMass = A.PRMs[prevA];
					for (int prevB = 0; prevB < b; ++prevB) {
						int prevBMass = B.PRMs[prevB];

						if (AntibodyUtils.EquivalentScaled(BMass - prevBMass,
								AMass - prevAMass, Tolerance)) {
							double NewScore = AScore + BScore
									+ Alignment[prevA][prevB];
							if (NewScore > MaxScore) {
								MaxScore = NewScore;
								bestPrev[0] = prevA;
								bestPrev[1] = prevB;
							}
						}

					}
				}
				Alignment[a][b] = MaxScore;
				prevCell[a][b] = bestPrev[0]*B.PRMs.length + bestPrev[1];
			}
		}
		Alignment alnObj = new Alignment(Alignment, prevCell, alignEndType.any, false);
		return alnObj;
	}
	
	/**
	 * Returns the score table for the alignment of two PRM spectra with a given PRM tolerance.  Any number of missed PRMs is allowed
	 * but no modifications.  Alignment must start at the beginning of each spectrum and end at the end.  Spectra are
	 * expected to be permutations of each other but same parent mass.
	 * @param A A PRMSpectrum object with N peaks
	 * @param B A PRMSpectrum object with M peaks
	 * @param Tolerance
	 * @return Returns a NxM table where position (i,j) is the highest score of an alignment where peak i of spectrum A is matched to peak j of spectrum B.
	 */
	private Alignment getAlignmentScoreTableFixedStart(PRMSpectrum A, PRMSpectrum B, double Tolerance) {
		
		double[][] Alignment = new double[A.PRMs.length][B.PRMs.length];
		int[][] prevCell = new int[A.PRMs.length][B.PRMs.length];
		
		Alignment[0][0] = A.PRMScores[0] + B.PRMScores[0];
		
		//The alignment must start at the beginning of both spectra, so no starting anywhere by (0,0)
		for(int a = 0; a < A.PRMs.length; ++a)
			Alignment[a][0] = Double.MIN_VALUE;
		for(int b = 0; b < B.PRMs.length; ++b)
			Alignment[0][b] = Double.MIN_VALUE;
		
		//Iterate over each row
		for (int a = 1; a < A.PRMs.length; ++a) {
			
			double AScore = A.PRMScores[a];
			int AMass = A.PRMs[a];
			
			//Iterate over each column
			for (int b = 1; b < B.PRMs.length; ++b) {
			
				double BScore = B.PRMScores[b];
				int BMass = B.PRMs[b];

				double MaxScore = AScore + BScore;
				int[] bestPrev = {-1,-1};
				//Consider every previous pair
				for (int prevA = 0; prevA < a; ++prevA) {
					int prevAMass = A.PRMs[prevA];
					
					for (int prevB = 0; prevB < b; ++prevB) {
						int prevBMass = B.PRMs[prevB];

						if (AntibodyUtils.EquivalentScaled(BMass - prevBMass,
								AMass - prevAMass, Tolerance)) {
							double NewScore = AScore + BScore
									+ Alignment[prevA][prevB];
							if (NewScore > MaxScore) {
								MaxScore = NewScore;
								bestPrev[0] = prevA;
								bestPrev[1] = prevB;
							}

						}
					}

				}
				Alignment[a][b] = MaxScore;
				prevCell[a][b] = bestPrev[0]*B.PRMs.length + bestPrev[1];
			}
		}
		Alignment alnObj = new Alignment(Alignment,prevCell,alignEndType.doubleFixed,false);
		return alnObj;
	}
}
