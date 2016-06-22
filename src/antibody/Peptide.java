package antibody;

import java.util.ArrayList;

import basicUtils.Utils;
import errorUtils.ErrorThrower;

/**
 * The peptide class is designed to group together Spectrum Nodes from the same
 * peptide sequence to make Seed creation more efficient.
 * 
 * @author natalie
 * 
 */
public class Peptide {

	private ArrayList<SpectrumNode> spectra = null;

	private String peptideSequence = null;
	private String unModdedSequence = null;

	public Peptide(SpectrumNode s) {
		if (s == null)
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Cannot create Peptide object from null SpectrumNode");
		if (s.getModdedDBAnnotation() == null)
			ErrorThrower
					.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
							"Cannot create Peptide object from SpectrumNode with no sequence");

		this.peptideSequence = s.getModdedDBAnnotation();
		this.unModdedSequence = Utils.GetUnModded(this.peptideSequence);
		this.spectra = new ArrayList();
		this.spectra.add(s);
	}

	public Peptide(String pepSeq) {
		if (pepSeq == null || pepSeq.length() == 0)
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Cannot create Peptide object from null/empty string");
		this.peptideSequence = pepSeq;
		this.unModdedSequence = Utils.GetUnModded(this.peptideSequence);

	}

	public String getPeptideSeq() {
		return this.peptideSequence;
	}

	public String getUnModdedPeptideSeq() {
		return this.unModdedSequence;
	}

	public int getNumSpectra() {
		if (this.spectra == null)
			return 0;
		return this.spectra.size();
	}

	public ArrayList<SpectrumNode> getSpectra() {
		return this.spectra;
	}

	public void addNewSpectrum(SpectrumNode s) {
		String pep = s.getModdedDBAnnotation();
		if (!pep.equals(this.peptideSequence))
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"Cannot merge spectra with peptides '"
							+ this.peptideSequence + "' and '" + pep + "'!");
		if (this.spectra == null)
			this.spectra = new ArrayList();
		this.spectra.add(s);

	}

	public SpectrumNode getSpectrum(int i) {
		if (this.spectra == null || this.spectra.size() <= i) {
			ErrorThrower.ThrowWarningCustom(ErrorThrower.CUSTOM_WARNING,
					"Desired spectrum node (" + i
							+ ") is beyond number of spectra");
			return null;
		}
		return (SpectrumNode) (this.spectra.get(i));
	}
}
