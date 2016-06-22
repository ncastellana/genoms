package antibody;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;

import basicUtils.PeptideSpectrumMatch;
import basicUtils.Utils;
//import basicUtils.InspectAnnotation;

/**
 * This class contains information about spectrum, as a spectrum is annotated it
 * is added to the list of SpectrumNodes. The node is part of a doubly linked
 * list
 * 
 * @author natalie
 * 
 */
public class SpectrumNode {

	private String SourceFileName;
	private int[] SpecIndex;
	private PRMSpectrum Spectrum;

	// private InspectAnnotation Annotation;
	private String dbAnnotation;

	private MSAlignmentType MSAnnotation;

	public HMMAlignment HMMAnnotation;

	private DistributionType NullDistribution = null;

	private int numSpectraMerged = 1;

	public SpectrumNode(String FileName, int ScanNumber) {
		this.SourceFileName = FileName;
		this.SpecIndex = new int[1];
		this.SpecIndex[0] = ScanNumber;
	}

	public SpectrumNode(PRMSpectrum Spectrum) {
		this.Spectrum = Spectrum;
		this.SourceFileName = Spectrum.FileName;
		this.SpecIndex = Spectrum.SpecIndex;
	}

	public void setNumSpectraMerged(int num) {
		this.numSpectraMerged = num;
	}

	public int getNumSpectraMerged() {
		return this.numSpectraMerged;
	}

	public void setDBAnnotation(PeptideSpectrumMatch Annotation) {

		this.dbAnnotation = Annotation.getModifiedPeptide();
	}

	public void setDBAnnotation(String ann) {
		this.dbAnnotation = ann;
	}

	public String GetSourceFileName() {
		return this.SourceFileName;
	}

	public int[] GetSpecIndex() {
		return this.SpecIndex;
	}

	/*
	 * public InspectAnnotation GetInspectAnnotation() { return this.Annotation;
	 * }
	 */

	public String getModdedDBAnnotation() {
		return this.dbAnnotation;
	}
	public String getUnModdedDBAnnotation() {
		return Utils.GetUnModded(this.dbAnnotation);
	}

	public void SetNull(DistributionType Dist) {
		this.NullDistribution = Dist;
	}

	public DistributionType GetNull() {
		return this.NullDistribution;
	}

	public PRMSpectrum GetPRMSpectrum() {
		return this.Spectrum;
	}

	public void SetMSAlignment(MSAlignmentType Alignment) {
		this.MSAnnotation = Alignment;

	}

	public MSAlignmentType GetMSAnnotation() {
		return this.MSAnnotation;
	}

	public static boolean WriteAllPRMSpectraToPKLBIN_OldFormat(
			String starSpectrumFile, SpectrumNode[] AllSpectra) {

		DataOutputStream f = null;

		try {
			f = new DataOutputStream(new FileOutputStream(starSpectrumFile));
			// The first value is an integer containing the number of spectra
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}
		try {
			f.writeInt(Integer.reverseBytes(AllSpectra.length));
			// f.writeShort(Short.reverseBytes(arg0))
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}

		// Next is an array of spectrum sizes (shorts)
		for (int i = 0; i < AllSpectra.length; ++i) {
			PRMSpectrum pSpec = AllSpectra[i].GetPRMSpectrum();
			try {
				f.writeShort(Short.reverseBytes((short) (pSpec.PRMs.length)));
			} catch (IOException E) {
				E.printStackTrace();
				return false;
			}
		}

		// Finally, each spectrum has a float value for the precursor mass,
		// precursor charge, and then the peak list (mass, log Odds).
		for (int i = 0; i < AllSpectra.length; ++i) {
			PRMSpectrum pSpec = AllSpectra[i].GetPRMSpectrum();

			try {
				f.writeInt(Integer.reverseBytes(Float
						.floatToIntBits(((float) (pSpec.PrecursorMass) / AntibodyUtils.MASS_SCALE) + 19)));
				if (pSpec.PrecursorCharge > 0)
					f.writeInt(Integer.reverseBytes(Float
							.floatToIntBits((float) (pSpec.PrecursorCharge))));
				else
					f.writeInt(Integer.reverseBytes(Float
							.floatToIntBits((float) (1.0))));
				// f.writeFloat((float)(pSpec.PrecursorMass));
				// f.writeFloat((float)(pSpec.PrecursorCharge));
			} catch (IOException E) {
				E.printStackTrace();
				return false;
			}
			for (int j = 0; j < pSpec.PRMs.length; ++j) {
				try {
					f.writeInt(Integer.reverseBytes(Float
							.floatToIntBits(((float) (pSpec.PRMs[j]))
									/ AntibodyUtils.MASS_SCALE)));
					f.writeInt(Integer.reverseBytes(Float
							.floatToIntBits((float) (pSpec.PRMScores[j]))));
					// f.writeFloat((float)(pSpec.PRMs[j]));
					// f.writeFloat((float)(pSpec.PRMScores[j]));
				} catch (IOException E) {
					E.printStackTrace();
					return false;
				}
			}
		}
		try {
			f.close();
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}
		return true;
	}

	/**
	 * Instead of trying to keep up with the changing pklbin formats, we just
	 * write to MGF format and convert using the SpecSet calls in main_specnets
	 * 
	 * @param starSpectrumFile
	 * @param allSpectra
	 * @return
	 */
	public static boolean writeAllPRMSpectraToMGF(String starSpectrumFile,
			SpectrumNode[] allSpectra) {

		FileWriter f = Utils.openFileWriter(starSpectrumFile);

		for (int i = 0; i < allSpectra.length; ++i) {
			PRMSpectrum pSpec = allSpectra[i].GetPRMSpectrum();
			
			
			/*
			 * BEGIN IONS PEPMASS=395.42517 CHARGE=2+
			 * FILENAME=specs_scored.pklbin ACTIVATION=PRM TITLE=Scan Number: 1
			 * SCANS=1
			 */
			Utils.writeLine(f, starSpectrumFile, "BEGIN IONS\n");
			Utils.writeLine(
					f,
					starSpectrumFile,
					"PEPMASS="
							+ Utils.GetMZFromM(((double) pSpec.PrecursorMass)
									/ AntibodyUtils.MASS_SCALE,
									pSpec.PrecursorCharge) + "\n");
			Utils.writeLine(f, starSpectrumFile, "CHARGE="
					+ pSpec.PrecursorCharge + "\n");
			Utils.writeLine(
					f,
					starSpectrumFile,
					"FILENAME="
							+ Utils.GetBaseNameNoExtension(starSpectrumFile)
							+ ".pklbin\n");
			Utils.writeLine(f, starSpectrumFile, "ACTIVATION=PRM\n");
			Utils.writeLine(f, starSpectrumFile, "TITLE=Scan Number: "
					+ pSpec.ScanNumber[0] + "\n");
			Utils.writeLine(f, starSpectrumFile, "SCANS=" + pSpec.ScanNumber[0]
					+ "\n");

			for (int j = 0; j < pSpec.PRMs.length; ++j) {
				Utils.writeLine(f, starSpectrumFile, ((double) (pSpec.PRMs[j]))
						/ AntibodyUtils.MASS_SCALE + " " + pSpec.PRMScores[j]
						+ "\n");
			}
			Utils.writeLine(f, starSpectrumFile, "END IONS\n");
		}

		Utils.closeFileWriter(f, starSpectrumFile);
		return true;

	}

	public static boolean writeAllPRMSpectraToPKLBIN(String starSpectrumFile,
			SpectrumNode[] AllSpectra) {

		DataOutputStream f = Utils.openDataOutputStream(starSpectrumFile);

		// empty 4 bytes signfies latest version?

		try {

			int size = (Integer.SIZE * 2 + Byte.SIZE * 2) / 8;
			// System.out.println("SIZE: " + size);
			ByteBuffer b = ByteBuffer.allocate(10);
			b.putInt(Integer.reverseBytes(0));
			b.put((byte) 1);
			b.put((byte) 0);
			b.putInt(Integer.reverseBytes(AllSpectra.length));
			f.write(b.array());
			// f.writeShort(Short.reverseBytes(arg0))
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}
		// Next is an array of scan numbers
		for (int i = 0; i < AllSpectra.length; ++i) {
			// PRMSpectrum pSpec = AllSpectra[i].GetPRMSpectrum();
			try {
				f.writeInt(Integer.reverseBytes(i));
			} catch (IOException E) {
				E.printStackTrace();
				return false;
			}
		}
		// Next is an array of msLevels
		for (int i = 0; i < AllSpectra.length; ++i) {
			try {
				f.writeShort(Short.reverseBytes((short) 0));
			} catch (IOException E) {
				E.printStackTrace();
				return false;
			}
		}
		// Next is an array of spectrum sizes (shorts)
		for (int i = 0; i < AllSpectra.length; ++i) {
			PRMSpectrum pSpec = AllSpectra[i].GetPRMSpectrum();
			try {
				f.writeShort(Short.reverseBytes((short) (pSpec.PRMs.length)));
			} catch (IOException E) {
				E.printStackTrace();
				return false;
			}
		}

		// Finally, each spectrum has a float value for the precursor mass,
		// precursor charge, and then the peak list (mass, log Odds).
		for (int i = 0; i < AllSpectra.length; ++i) {
			PRMSpectrum pSpec = AllSpectra[i].GetPRMSpectrum();

			try {
				f.writeInt(Integer.reverseBytes(Float
						.floatToIntBits(((float) (pSpec.PrecursorMass) / AntibodyUtils.MASS_SCALE) + 19)));
				// if(pSpec.PrecursorCharge > 0)
				// f.writeInt(Integer.reverseBytes(Float.floatToIntBits((float)(pSpec.PrecursorCharge))));
				// else
				// f.writeInt(Integer.reverseBytes(Float.floatToIntBits((float)(1.0))));
				f.writeInt(Integer.reverseBytes(Float
						.floatToIntBits((float) (0.0))));
				// f.writeFloat((float)(pSpec.PrecursorMass));
				// f.writeFloat((float)(pSpec.PrecursorCharge));
			} catch (IOException E) {
				E.printStackTrace();
				return false;
			}
			for (int j = 0; j < pSpec.PRMs.length; ++j) {
				try {
					f.writeInt(Integer.reverseBytes(Float
							.floatToIntBits(((float) (pSpec.PRMs[j]))
									/ AntibodyUtils.MASS_SCALE)));
					f.writeInt(Integer.reverseBytes(Float
							.floatToIntBits((float) (pSpec.PRMScores[j]))));
					// f.writeFloat((float)(pSpec.PRMs[j]));
					// f.writeFloat((float)(pSpec.PRMScores[j]));
				} catch (IOException E) {
					E.printStackTrace();
					return false;
				}
			}
		}
		try {
			f.close();
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}
		return true;
	}

	public void DebugPrint() {
		System.out.println("SpectrumNode: " + this.SourceFileName + ":"
				+ this.SpecIndex);
		System.out.print("  Has InspectAnn: ");
		if (this.dbAnnotation != null)
			System.out.println("true [" + this.dbAnnotation + "]");
		else
			System.out.println("false");

		System.out.print("  Has MSAlignment: ");
		if (this.MSAnnotation != null)
			System.out.println("true [" + this.MSAnnotation.GetAnchorSequence()
					+ ":" + this.MSAnnotation.GetScore() + "]");
		else
			System.out.println("false");

		System.out.print("  Has HMMAlignment: ");
		if (this.HMMAnnotation != null) {
			System.out.println("true, MassShift: "
					+ this.HMMAnnotation.GetMassShift());
			System.out.println("MLPath: " + this.HMMAnnotation.GetPathString());
			System.out.println("AlignmentType: "
					+ this.HMMAnnotation.GetAlignmentType());
		} else
			System.out.println("false");

	}

	public void SetPRMSpectrum(PRMSpectrum currSpectrum) {
		this.Spectrum = currSpectrum;
		this.SourceFileName = this.Spectrum.FileName;
		this.SpecIndex = this.Spectrum.SpecIndex;

	}

	public HMMAlignment GetHMMAnnotation() {
		return this.HMMAnnotation;
	}

	public void SetHMMAnnotation(HMMAlignment object) {
		this.HMMAnnotation = (HMMAlignment) (object);

	}

}
