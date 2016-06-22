package antibody;

public class ModSite {

	public double MassDelta = 0.0;
	public int ModIndex = -1;
	public String Sequence = "";
	public String ModString = "";

	public ModSite(String ModStr, int Index, double Delta, String Peptide) {
		this.ModString = Peptide;
		this.ModString = ModStr;
		this.ModIndex = Index;
		this.MassDelta = Delta;
	}

	public char GetMutationSub() {
		if (this.ModString.indexOf("->") >= 0) {
			return this.ModString.charAt(this.ModString.length() - 1);
		} else
			return 0;
	}

	public void DebugPrint() {
		System.out.println("-----ModSite-----");
		System.out.println("Delta: " + this.MassDelta);
		System.out.println("Str: " + this.ModString);
		System.out.println("At pos " + this.ModIndex + " on:");
		System.out.println(this.Sequence);
	}
}
