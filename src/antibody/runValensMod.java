/**
 * 
 */
package antibody;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import basicUtils.MSGFPlusAnnotation;
import basicUtils.PeptideSpectrumMatch;
import basicUtils.Utils;
import dbUtils.ProteinDB;
import edu.ucsd.msjava.ui.MSGFPlus;
import edu.ucsd.msjava.ui.MzIDToTsv;
import errorUtils.ErrorThrower;
import specUtils.PKLBINUtils;


/**
 * @author natalie
 * TODO Create a concatenated pklbin file?
 *
 */
public class runValensMod {

	
	
	public static double FDR_CUTOFF = 0.01;
	public static String AA_DELIM = ",";
	public static String psmFileName = "tableSpectra.txt";
	public static String coverageFileName = "tableProteinCoverage.txt";
	private static String proteinFileName = "tableProtein.txt";
	
	//borrowed from valensMod.js
	public static String TABLE_SEP_L1 = "|";
	public static String TABLE_SEP_L2 = "@";
	
	public static String usageInfo = "antibody.runValensMod version 20151120.1\n" +
			"Options:\n" +
			" -i [FILE] Parameter file containing parameters of the following form, one per line:\n" +
			"     PARAM_NAME=Value\n";
	

	
	private String[] options = {"mod_file", "input_specs_ms","proteases","frag",
								"ref_db","search_db","pm_tol","instrument",
								"project_dir","exe_dir","report_server"};
	
	private int[] numValues = {1,1,1,1,
								1,1,1,1,
								1,1,1};
	
	private Hashtable<String, ArrayList<String[]>> parameters = null;
	
	private boolean fastMode = false;
	private String modFileName;
	private String[] spectrumFileNames;
	private String[] proteases;
	private String[] frag;
	private ProteinDB refDB;
	//private ProteinDB searchDB;
	private String searchDB;
	private double pmTol;
	private String inst;
	private String projectDir;
	private String cgibinURL; 
	private String exeDir;

	private String msgfplusOutputDir;
	private String reportDataDir;
	private String reportHTMLDir;
	
	private String[] dbSearchOutputFiles;
	private String[] pklbinFileNames;
	private String statusFileName;
	private String sortedDBSearchFile;
	
	
	public runValensMod(String paramFileName) {
		parameters = Utils.ParseCSPSInputFile(paramFileName, options,
				numValues);
		
		//Check that we have what we need
		if(!parameters.containsKey("mod_file")) {
			ErrorThrower.ThrowError(3, "MOD_FILE");
		}
		this.modFileName = parameters.get("mod_file").get(0)[0];
		
		
		//Parse the spectrum information
		if(!parameters.containsKey("input_specs_ms")) {
			ErrorThrower.ThrowError(3, "INPUT_SPECS_MS");
		}
		this.spectrumFileNames = parameters.get("input_specs_ms").get(0)[0].split(";");
		
		if(!parameters.containsKey("proteases")) {
			ErrorThrower.ThrowError(3, "PROTEASES");
		}
		this.proteases = parameters.get("proteases").get(0)[0].split(";");
		
		if(!parameters.containsKey("frag")) {
			ErrorThrower.ThrowError(3, "FRAG");
		}
		this.frag = parameters.get("frag").get(0)[0].split(";");
		
		if(this.spectrumFileNames.length != this.proteases.length || 
				this.spectrumFileNames.length != this.frag.length)
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Number of spectrum files must match number of proteases and fragmentation methods in the params file");
			
		
		//Load the database information
		if(!parameters.containsKey("ref_db")) {
			ErrorThrower.ThrowError(3, "REF_DB");
		}
		this.refDB = new ProteinDB(parameters.get("ref_db").get(0)[0]);
		if(!parameters.containsKey("search_db")) {
			ErrorThrower.ThrowError(3, "SEARCH_DB");
		}
		this.searchDB = parameters.get("search_db").get(0)[0];
		
		//Load global params
		if(!parameters.containsKey("pm_tol")) {
			ErrorThrower.ThrowError(3, "PM_TOL");
		}
		this.pmTol = Double.parseDouble(parameters.get("pm_tol").get(0)[0]);
		if(!parameters.containsKey("instrument")) {
			ErrorThrower.ThrowError(3, "INSTRUMENT");
		}
		this.inst = parameters.get("instrument").get(0)[0];
		
		if(!parameters.containsKey("project_dir")) {
			ErrorThrower.ThrowError(3, "PROJECT_DIR");
		}
		this.projectDir = parameters.get("project_dir").get(0)[0];
		
		if(!parameters.containsKey("report_server")) {
			ErrorThrower.ThrowError(3, "REPORT_SERVER");
		}
		this.cgibinURL = parameters.get("report_server").get(0)[0];
		
		if(!parameters.containsKey("exe_dir")) {
			ErrorThrower.ThrowError(3, "EXE_DIR");
		}
		this.exeDir = parameters.get("exe_dir").get(0)[0];
		this.statusFileName = this.projectDir + File.separator + "status.txt";
		Utils.writeStatusToFile(this.statusFileName, "Running");
	}
	
	/**
	 * 
	 */
	public void run() {
		
		this.setupProjectDirectory();
		this.runMSGFPlus();
		this.populateTables();
		this.writeHTML();
	}

	
	private void writeHTML() {
		
		System.out.println("[3] Writing report files");
		//copy stuff to the directory
		if(!Utils.copyDirectory(this.exeDir + File.separator + "resources" + File.separator + "css" + File.separator + "css", this.reportHTMLDir + File.separator + "css")) {
			Utils.writeStatusToFile(this.statusFileName, "Error");
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to copy HTML css directory to result directory");
		}
		if(!Utils.copyDirectory(this.exeDir + File.separator + "resources" + File.separator + "css" + File.separator + "deps", this.reportHTMLDir + File.separator + "deps")) {
			Utils.writeStatusToFile(this.statusFileName, "Error");
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to copy HTML deps directory to result directory");
		}
		if(!Utils.copyDirectory(this.exeDir + File.separator + "resources" + File.separator + "css" + File.separator + "js", this.reportHTMLDir + File.separator + "js")) {
			Utils.writeStatusToFile(this.statusFileName, "Error");
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to copy HTML js directory to result directory");
		}
		if(!Utils.copyDirectory(this.exeDir + File.separator + "resources" + File.separator + "css" + File.separator + "images", this.reportHTMLDir + File.separator + "images")) {
			Utils.writeStatusToFile(this.statusFileName, "Error");
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to copy HTML images directory to result directory");
		}
		if(!Utils.copyDirectory(this.exeDir + File.separator + "resources" + File.separator + "css" + File.separator + "styles", this.reportHTMLDir + File.separator + "styles")) {
			Utils.writeStatusToFile(this.statusFileName, "Error");
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to copy HTML styles directory to result directory");
		}
		String htmlFile = this.reportHTMLDir + File.separator + "index.html";
		FileWriter writer = Utils.openFileWriter(htmlFile);
		
		Utils.writeLine(writer, htmlFile, "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n");
		Utils.writeLine(writer, htmlFile, "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"en\" lang=\"en\">\n");
		Utils.writeLine(writer, htmlFile, "<head>\n");
		Utils.writeLine(writer, htmlFile, "  <meta http-equiv=\"content-type\" content=\"text/html; charset=iso-8859-1\" />\n");
		Utils.writeLine(writer, htmlFile, "  <meta http-equiv=\"Content-Language\" content=\"en-us\" />\n");
		Utils.writeLine(writer, htmlFile, "  <title>Digital Proteomics LLC</title>\n");
		Utils.writeLine(writer, htmlFile, "  <link rel=\"shortcut icon\" href=\"images/favicon.ico\" type=\"image/icon\" />\n");
		Utils.writeLine(writer, htmlFile, "  <script type=\"text/javascript\" language=\"javascript\" src=\"js/bst.js\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <script type=\"text/javascript\" language=\"javascript\" src=\"js/valensMod.js\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <script type=\"text/javascript\" language=\"javascript\" src=\"js/XHConn.js\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <link href=\"styles/main.css\" rel=\"stylesheet\" type=\"text/css\" />\n");
		Utils.writeLine(writer, htmlFile, "  <script type=\"text/javascript\" language=\"javascript\" src=\"js/prototype.js\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <script type=\"text/javascript\" language=\"javascript\" src=\"js/scriptaculous.js?load=effects,builder\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <script type=\"text/javascript\" language=\"javascript\" src=\"js/lightbox.js\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <link rel=\"stylesheet\" href=\"css/lightbox.css\" type=\"text/css\" media=\"screen\" />\n");
		Utils.writeLine(writer, htmlFile, "  <script src=\"js/jquery-1.8.2.min.js\" type=\"text/javascript\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <script src=\"js/jquery.cookie.js\" type=\"text/javascript\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <script src=\"js/jquery.contextMenu.js\" type=\"text/javascript\"></script>\n");
		Utils.writeLine(writer, htmlFile, "  <script src='js/jquery.base64.js' type='text/javascript'></script>\n");
		Utils.writeLine(writer, htmlFile, "  <link href='js/jquery.contextMenu.css' rel='stylesheet' type='text/css' />\n");
		Utils.writeLine(writer, htmlFile, "  <script>jQuery.noConflict();</script>\n");
		Utils.writeLine(writer, htmlFile, "  <script src='js/pwd.js' type='text/javascript'></script>\n");
		Utils.writeLine(writer, htmlFile, "  <link rel='stylesheet' href='css/pwd.css' type='text/css' media='screen' />\n");
		Utils.writeLine(writer, htmlFile, "</head>\n");
		Utils.writeLine(writer, htmlFile, "<body onload='javascript:init();'>\n");
		Utils.writeLine(writer, htmlFile, "  <div class='wrapper' align='center'>\n");
		Utils.writeLine(writer, htmlFile, "    <div><h3 align=center><img src='images/logo.jpg' /></h3></div>\n");
		Utils.writeLine(writer, htmlFile, "  <br /><br />\n");
		Utils.writeLine(writer, htmlFile, "  <div id='bodyWrapper'>\n");
		Utils.writeLine(writer, htmlFile, "    <div id='textWrapper'>\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='projectDir' value='" + this.projectDir + "' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='tablesDir' value='ValensModData' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='cellsPerLine' value='21' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='serverLocation' value='" + this.cgibinURL + "' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='fastaFilename' value='" + this.refDB.getDBFileName() + "' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='noClusters' value='1' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='contigRefs' value='0' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='msFilename' value='spectra/specs_ms.pklbin' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='allowRealign' value='0' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='dynamic' value='1' />\n");
		Utils.writeLine(writer, htmlFile, "      <input type='hidden' id='aas_file' value='/data/mnt/dev/sps/bin/' />\n");
		Utils.writeLine(writer, htmlFile, "      <div id='mainDiv2'>\n");
		Utils.writeLine(writer, htmlFile, "        <table align='center'>\n");
		Utils.writeLine(writer, htmlFile, "          <tr>\n");
		Utils.writeLine(writer, htmlFile, "            <td></td>\n");
		Utils.writeLine(writer, htmlFile, "            <td width='990px'>\n");
		Utils.writeLine(writer, htmlFile, "              <table class='mainform'>\n");
		Utils.writeLine(writer, htmlFile, "                <tr><th colspan='0'>Job Status</th></tr>\n");
		Utils.writeLine(writer, htmlFile, "                <tr>\n");
		Utils.writeLine(writer, htmlFile, "                  <td>\n");
		Utils.writeLine(writer, htmlFile, "                    <table class='sched ' width='100%'>  \n");
		Utils.writeLine(writer, htmlFile, "                      <tr><th width='25%' bgcolor='#003399'><span style='color:white'>Job</span></th>\n");
		Utils.writeLine(writer, htmlFile, "                        <td>" + Utils.GetBaseName(this.projectDir) + "</td>\n");
		Utils.writeLine(writer, htmlFile, "                      </tr>\n");
		//Utils.writeLine(writer, htmlFile, "                      <tr><th width='25%' bgcolor='#003399'><span style='color:white'>User</span></th>\n");
		//Utils.writeLine(writer, htmlFile, "                        <td>anonymous</td>\n");
		//Utils.writeLine(writer, htmlFile, "                      </tr>\n");
		Utils.writeLine(writer, htmlFile, "                      <tr><th width='25%' bgcolor='#003399'><span style='color:white'>Status</span></th>\n");
		Utils.writeLine(writer, htmlFile, "                        <td style='background-color:#ffffff;' id='status'></td>\n");
		Utils.writeLine(writer, htmlFile, "                      </tr>\n");
		Utils.writeLine(writer, htmlFile, "                      <tr><th width='25%' bgcolor='#003399' rowspan='2'><span style='color:white'>Data</span></th>\n");
		//Utils.writeLine(writer, htmlFile, "                        <td><a href='#' onclick='javascript:TablesAll.loadPage(50,0);'>Group by Contig</a></td>\n");
		//Utils.writeLine(writer, htmlFile, "                      </tr>\n");
		Utils.writeLine(writer, htmlFile, "                      <tr><td><a href='#' onclick='javascript:TablesAll.loadPage(20,0);'>Group by Protein</a></td></tr>\n");
		Utils.writeLine(writer, htmlFile, "                      <tr><th width='25%' bgcolor='#003399'><span style='color:white'>Spectrum Data</span></th>\n");
		Utils.writeLine(writer, htmlFile, "                        <td>\n");
		for(int i = 0; i < this.spectrumFileNames.length; ++i) {
			Utils.writeLine(writer, htmlFile, "                          <a  href='#' onclick='javascript:TablesAll.loadPage(80," + i + ");'>Group by <i>" + Utils.GetBaseName(this.spectrumFileNames[i]) + "</i></a><br />\n");
		}
		Utils.writeLine(writer, htmlFile, "                        </td>\n");
		Utils.writeLine(writer, htmlFile, "                      </tr>\n");
		//Utils.writeLine(writer, htmlFile, "                      <tr><th width='25%' bgcolor='#003399'><span style='color:white'>pViz</span></th>\n");
		//Utils.writeLine(writer, htmlFile, "                        <td><a href='valens.pviz.html'>pViz View</a></td>\n");
		//Utils.writeLine(writer, htmlFile, "                      </tr>\n");
		//Utils.writeLine(writer, htmlFile, "                      <tr><th width='25%' bgcolor='#003399'><span style='color:white'>Tags</span></th>\n");
		//Utils.writeLine(writer, htmlFile, "                        <td><a href='tags.main.html'>PepNovo Tags View</a></td>\n");
		//Utils.writeLine(writer, htmlFile, "                      </tr>\n");
		Utils.writeLine(writer, htmlFile, "                    </table>\n"); //end class = sched
		Utils.writeLine(writer, htmlFile, "                  </td>\n");
		Utils.writeLine(writer, htmlFile, "                </tr>\n");
		Utils.writeLine(writer, htmlFile, "                <tr>\n");
		Utils.writeLine(writer, htmlFile, "                  <td colspan='0' class='bottomline'>&nbsp;</td>\n");
		Utils.writeLine(writer, htmlFile, "                </tr>\n");
		Utils.writeLine(writer, htmlFile, "              </table>\n"); //end class = mainform
		Utils.writeLine(writer, htmlFile, "            </td>\n");
		Utils.writeLine(writer, htmlFile, "            <td></td>\n");
		Utils.writeLine(writer, htmlFile, "          </tr>\n");
		Utils.writeLine(writer, htmlFile, "        </table>\n"); //end align = center
		Utils.writeLine(writer, htmlFile, "      </div>\n"); //End mainDiv2
		Utils.writeLine(writer, htmlFile, "      <div id='mainDiv' style='display: none'></div>\n");
		Utils.writeLine(writer, htmlFile, "    </div>\n"); //End textWrapper
		Utils.writeLine(writer, htmlFile, "  </div>\n"); //End bodyWrapper
		Utils.writeLine(writer, htmlFile, "  <div class='push'></div>\n");
		Utils.writeLine(writer, htmlFile, "</div>\n"); //End class = wrapper
		Utils.writeLine(writer, htmlFile, "<div class='footer'>\n");
		Utils.writeLine(writer, htmlFile, "  <table width='100%'>\n");
		Utils.writeLine(writer, htmlFile, "    <tr><td class='VHSep'></td></tr>\n");
		Utils.writeLine(writer, htmlFile, "    <tr>\n");
		Utils.writeLine(writer, htmlFile, "      <td class='HSep'></td>\n");
		Utils.writeLine(writer, htmlFile, "      <td class='ln'></td>\n");
		Utils.writeLine(writer, htmlFile, "      <td class='HSep'></td>\n");
		Utils.writeLine(writer, htmlFile, "    </tr>\n");
		Utils.writeLine(writer, htmlFile, "    <tr>\n");
		Utils.writeLine(writer, htmlFile, "      <td class='HSep'></td>\n");
		Utils.writeLine(writer, htmlFile, "      <td class='Footer'>Digital Proteomics LLC</td>\n");
		Utils.writeLine(writer, htmlFile, "    </tr>\n");
		Utils.writeLine(writer, htmlFile, "  </table>\n");
		Utils.writeLine(writer, htmlFile, "</div></body>\n");
				
		Utils.writeLine(writer, htmlFile, "<div id='navButtons' style='position: fixed; top: 0px; left: 0px;'></div>\n");
		Utils.writeLine(writer, htmlFile, "<ul id='myMenu' class='contextMenu'>\n");
		Utils.writeLine(writer, htmlFile, "  <li class='newtab'><a href='#newtab'>Open in new tab</a></li>\n");
		Utils.writeLine(writer, htmlFile, "  <li class='newwin'><a href='#newwin'>Open in new window</a></li>\n");
		Utils.writeLine(writer, htmlFile, "</ul>\n");
		Utils.writeLine(writer, htmlFile, "<div><a id='login-link' href='#login-box' class='login-window'></a></div>\n");
		Utils.writeLine(writer, htmlFile, "<div id='login-box' class='login-popup'>\n");
		Utils.writeLine(writer, htmlFile, "  <form id='pwdForm' method='post' class='signin' action='#'>\n");
		Utils.writeLine(writer, htmlFile, "    <fieldset class='textbox' id='pwdField'>\n");
		Utils.writeLine(writer, htmlFile, "      <label class='password'><span>Password</span><input id='password' name='password' value='' type='password' placeholder='Password'></label>\n");
		Utils.writeLine(writer, htmlFile, "      <button class='button' type='button' onclick='checkIds2();return false;'>Unlock</button>\n");
		Utils.writeLine(writer, htmlFile, "    </fieldset>\n");
		Utils.writeLine(writer, htmlFile, "  </form>\n");
		Utils.writeLine(writer, htmlFile, "</div>\n");
		Utils.writeLine(writer, htmlFile, "</html>\n");
		
		Utils.closeFileWriter(writer, htmlFile);
	}

	/**
	 * 
	 */
	private void populateTables() {
		
		boolean localDebug = true;
		System.out.println("[2] Populating data tables");
		HashSet<String> seenScans = new HashSet<String>(); //To make sure we only count each spectrum once (even if it has multiple similar candidates
		
		ArrayList<String[]> splitSequences = new ArrayList<String[]>();
		ArrayList<String> sequences = new ArrayList<String>();
		ArrayList <String[]> counts = new ArrayList<String[]>();
		
		
		//Populate these with the unmodified numbers
		long numRefProteins = this.refDB.getProteinCount();
		
		int[] specCounts = new int[(int) numRefProteins];
		specCounts = Utils.initializeIntArray(specCounts, 0);
		
		for(int i = 0; i < numRefProteins; ++i) {
			String seq = this.refDB.getProteinSequence(i);
			sequences.add(seq);
			splitSequences.add(Utils.splitSequenceWithModsToArray(seq));
			
			String[] currCount = new String[seq.length()];
			currCount = Utils.initializeStringArray(currCount, "0");
			counts.add(currCount);
		}
		
		//TableIndex;SpectrumIndex;Scan;ClusterIndex;ProteinName;PklbinFileIndex;
		//PklbinFilename;ReferenceSequence;HomologSequence;DeNovoSequence;UserSequence;
		//Mass;Charge;B_Per;Y_Per;BY_Int;OriginalFilename;ContigIndex;ProteinIndex;Tool;FragmentationModel
		
		int tableIdx = 0;
		
		
		String outputFileName = this.reportDataDir + File.separator + runValensMod.psmFileName;
		FileWriter writer = Utils.openFileWriter(outputFileName);
		Utils.writeLine(writer, outputFileName, "#TableIndex;SpectrumIndex;Scan;ClusterIndex;ProteinName;PklbinFileIndex;PklbinFilename;ReferenceSequence;HomologSequence;DeNovoSequence;UserSequence;Mass;Charge;B_Per;Y_Per;BY_Int;OriginalFilename;ContigIndex;ProteinIndex;Tool;FragmentationModel\n");
		
		//for(int i = 0; i < this.dbSearchOutputFiles.length; ++i) {
			if(localDebug)
				System.out.println(" - Parsing database search file " + this.sortedDBSearchFile);
			PeptideSpectrumMatch[] anns = MSGFPlusAnnotation.LoadResultsFileFilter(this.sortedDBSearchFile, MSGFPlusAnnotation.qValue, FDR_CUTOFF, true);
			if(localDebug)
				System.out.println(" - Loaded " + anns.length + " peptide spectrum matches");
			
			
			for(int j = 0; j < anns.length; ++j) {
				
				String unModPep = anns[j].getUnModifiedPeptide();
				String[] pepSeq = Utils.splitSequenceWithModsToArray(anns[j].getModifiedPeptide());
				
				String key = AntibodyUtils.CreateSpectrumKey(anns[j].getSpectrumFileName(), anns[j].getSpectrumIndex());
				if(seenScans.contains(key))
					continue;
			
				seenScans.add(key);
				
				int proteinIdx = -1;
				
				//Now add it to the coverage
				for(int k = 0; k < sequences.size(); ++k) {
					int index = sequences.get(k).indexOf(unModPep);
					if(index >= 0) {
						proteinIdx = k;
						specCounts[k] += 1;
						String[] currSplitSequence = splitSequences.get(k);
						String[] currCounts = counts.get(k);
						
						for(int pos = 0; pos < pepSeq.length; ++pos) {
							
							int aaIndex = Utils.getIndexOfStringInStringList(pepSeq[pos],currSplitSequence[index+pos],AA_DELIM);
							if(aaIndex < 0) //We have not seen this mod before
							{
								currSplitSequence[index+pos] += AA_DELIM + pepSeq[pos];
								currCounts[index+pos] += AA_DELIM + "1";
							}
							else {
								String[] bits = currCounts[index+pos].split(AA_DELIM);
								int currCount = Integer.parseInt(bits[aaIndex]) + 1;
								bits[aaIndex] = currCount + "";
								currCounts[index+pos] = Utils.JoinStringArray(bits, AA_DELIM);
							}
						}
						splitSequences.set(k, currSplitSequence);
						counts.set(k, currCounts);
						k = sequences.size(); //We only add it to the first one we find
					}
				}//end adding coverage
				
				int fileIdx = -1;
				for(int fNum = 0; fNum < this.spectrumFileNames.length; ++fNum) {
					if(anns[j].getSpectrumFileName().equalsIgnoreCase(Utils.GetBaseName(this.spectrumFileNames[fNum])))
						fileIdx = fNum;
				}
				if(fileIdx < 0) {
					ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to find '" + anns[j].getSpectrumFileName() + "' in spectrum file list");
				}
				
				//Add to spectrum table
				int specIdx = (anns[j].getSpectrumIndex()[0]+1);
				int scan = anns[j].getScanNumber()[0];
				if(scan < 0)
					scan = specIdx;
				String tableLine = tableIdx + ";" + specIdx + ";" + scan + ";";
				tableLine += scan + ";" + anns[j].getStringAttribute("Protein") + ";";
				tableLine += fileIdx + ";" + File.separator + "spectra" + File.separator + Utils.GetBaseName(this.pklbinFileNames[fileIdx]) + ";" + AntibodyUtils.ConvertMSGFPlusPeptide2SPS(anns[j].getModifiedPeptide()) + ";" + AntibodyUtils.ConvertMSGFPlusPeptide2SPS(anns[j].getModifiedPeptide()) + ";";
				tableLine += AntibodyUtils.ConvertMSGFPlusPeptide2SPS(anns[j].getModifiedPeptide()) + ";;" + anns[j].getDoubleAttribute("Precursor") + ";" + anns[j].getIntAttribute("Charge") + ";";
				tableLine += anns[j].getDoubleAttribute("PrecursorError(Da)") + ";" + anns[j].getDoubleAttribute("QValue") + ";0;" + Utils.GetBaseName(this.spectrumFileNames[fileIdx]) + ";-1;" + proteinIdx + ";ValensMod;" + this.frag[fileIdx].toUpperCase() + "\n";
				
				Utils.writeLine(writer, outputFileName, tableLine);
				tableIdx += 1;
			}
		//}
		
		//Write coverage table
		//Write protein table
		//ProteinIndex;ProteinName;ProteinDescription;ContigCount;SpectraCount;AAsCount;CoveragePercentage;ProteinSequence
		String coverageOutputFileName = this.reportDataDir + File.separator + runValensMod.coverageFileName;
		FileWriter coverageWriter = Utils.openFileWriter(coverageOutputFileName);
		Utils.writeLine(coverageWriter, coverageOutputFileName, "#ProteinIndex;ProteinName;ProteinLen;ProteinSequence;SpecCountData\n");
		
		String proteinFileName = this.reportDataDir + File.separator + runValensMod.proteinFileName;
		FileWriter proteinWriter = Utils.openFileWriter(proteinFileName);
		Utils.writeLine(proteinWriter, proteinFileName, "#ProteinIndex;ProteinName;ProteinDescription;ContigCount;SpectraCount;AAsCount;CoveragePercentage;ProteinSequence\n");
		
		for(int i = 0; i < numRefProteins; ++i) {
			String pName = this.refDB.getProteinName(i);
			
			String proteinTableLine = i + ";" + pName + ";;0;" + specCounts[i] + ";" + sequences.get(i).length() + ";" + this.computeCoverage(counts.get(i)) + ";";
			proteinTableLine += sequences.get(i) + "\n";
			Utils.writeLine(proteinWriter, proteinFileName, proteinTableLine);
			
			
			String coverageTableLine = i + ";" + pName + ";" + sequences.get(i).length() + ";" + Utils.JoinStringArray(splitSequences.get(i), TABLE_SEP_L1) + ";" + Utils.JoinStringArray(counts.get(i), TABLE_SEP_L2) + "\n";
			Utils.writeLine(coverageWriter, coverageOutputFileName, coverageTableLine);
			
			//Write the protein name
			//Utils.writeLine(coverageWriter, coverageOutputFileName, pName + "\n");
			
			//Write the split sequence
			//Utils.writeLine(coverageWriter, coverageOutputFileName, Utils.JoinStringArray(splitSequences.get(i), ";") + "\n");
			
			//Write the split coverage
			//Utils.writeLine(coverageWriter, coverageOutputFileName, Utils.JoinStringArray(counts.get(i), ";") + "\n");
			
		}
		Utils.closeFileWriter(proteinWriter, proteinFileName);
		Utils.closeFileWriter(coverageWriter, coverageOutputFileName);
		Utils.closeFileWriter(writer, outputFileName);
	}

	/**
	 * Calculates the percentage of amino acids with a non-zero count
	 * @param strings
	 * @return
	 */
	private double computeCoverage(String[] counts) {
		double fullLength = counts.length;
		double covered = 0;
		for(int i = 0; i < counts.length; ++i) {
			String[] bits = counts[i].split(AA_DELIM);
			
			for(int j = 0; j < bits.length; ++j) {
				if(Integer.parseInt(bits[j]) > 0) {
						covered += 1;
						j = bits.length;
				}
			}
		}
		
		return covered/fullLength;
	}

	/**
	 * Runs MSGFPlus on each spectrum file
	 */
	private void runMSGFPlus() {
		
		System.out.println("[1] Running MSGF+");
		
		this.sortedDBSearchFile = this.msgfplusOutputDir + File.separator + "combinedDBSearchResults.sorted.tsv";
		if(this.fastMode && Utils.IsFile(this.sortedDBSearchFile)) {
			System.out.println(" - Skipping MSGF+ since we are in fastmode and the files exist");
		}
		//iterate over each spectrum file
		this.dbSearchOutputFiles = new String[this.spectrumFileNames.length];
		
		for (int i = 0; i < this.spectrumFileNames.length; ++i) {
			
		
			this.dbSearchOutputFiles[i] = this.msgfplusOutputDir + File.separator +
									Utils.GetBaseNameNoExtension(this.spectrumFileNames[i]) + ".mzid";
			
			
			//Create the MSGFPlus arguments
			ArrayList<String> argList = new ArrayList<String>();

			argList.add("-s");
			argList.add(this.spectrumFileNames[i]);

			argList.add("-d");
			argList.add(this.searchDB);

			argList.add("-t");
			argList.add(this.pmTol + "Da");

			//argList.add("-thread");
			//argList.add("1");

			argList.add("-o");
			argList.add(dbSearchOutputFiles[i]);
			
			argList.add("-m");
			if(this.frag[i].equalsIgnoreCase("CID"))
				argList.add("1");
			else if(this.frag[i].equalsIgnoreCase("ETD"))
				argList.add("2");
			else if(this.frag[i].equalsIgnoreCase("HCD"))
				argList.add("3");
			else
				argList.add("0");
			

			argList.add("-inst");
			if (this.inst.equalsIgnoreCase("LTQ")) {
				argList.add("0");
			} else if (this.inst.equalsIgnoreCase("ORBI")) {
				argList.add("1");
			} else if (this.inst.equalsIgnoreCase("TOF")) {
				argList.add("2");
			}else if(this.inst.equalsIgnoreCase("QEXT"))
				argList.add("3");
			else
				argList.add("0");

			argList.add("-e");
			if(this.proteases[i].equalsIgnoreCase("trypsin"))
				argList.add("1");
			else if(this.proteases[i].equalsIgnoreCase("chymo"))
				argList.add("2");
			else if(this.proteases[i].equalsIgnoreCase("lysc"))
				argList.add("3");
			else if(this.proteases[i].equalsIgnoreCase("lysn"))
				argList.add("4");
			else if(this.proteases[i].equalsIgnoreCase("gluc"))
				argList.add("5");
			else if(this.proteases[i].equalsIgnoreCase("argc"))
				argList.add("6");
			else if(this.proteases[i].equalsIgnoreCase("aspn"))
				argList.add("7");
			else if(this.proteases[i].equalsIgnoreCase("alp"))
				argList.add("8");
			else
				argList.add("9");

			argList.add("-mod");
			argList.add(this.modFileName);
			
			argList.add("-tda");
			argList.add("1");

			argList.add("-ntt");
			argList.add("1");

			String[] argv = Utils.ConvertArraylistToStringArray(argList);
			System.out.println("  MSGFPlus args: "
								+ Utils.JoinStringArray(argv, " "));
			if(this.fastMode && Utils.IsFile(this.dbSearchOutputFiles[i])) {
				System.out.println(" - Skipping MSGF+ for " + this.spectrumFileNames[i] + " since we are in fastmode and file exists");
			}
			else {
				try {
					MSGFPlus.main(argv);
				} catch (Exception e) {
					Utils.writeStatusToFile(this.statusFileName, "Error");
						e.printStackTrace();
						ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
									"Unable to run MSGFPlus!");
				}
			}
			argList.clear();
			argList.add("-i");
			argList.add(this.dbSearchOutputFiles[i]);
			
			argv = Utils.ConvertArraylistToStringArray(argList);
			System.out.println("  MzIDToTsv args: "
					+ Utils.JoinStringArray(argv, " "));
			try {
				MzIDToTsv.main(argv);
			}catch (Exception e) {
				Utils.writeStatusToFile(this.statusFileName, "Error");
				e.printStackTrace();
				ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
							"Unable to run MzIDToTsv!");
			}
			this.dbSearchOutputFiles[i] = Utils.GetFileNameNoExtension(this.dbSearchOutputFiles[i]) + ".tsv";
		}
		
		//Concatenate all the files
		String combinedDBSearchOutputFile = this.msgfplusOutputDir + File.separator + "combinedDBSearchResults.tsv";
		String temp = this.msgfplusOutputDir + File.separator + "combinedDBSearchResults.reformat.tsv";
		String temp2 = this.msgfplusOutputDir + File.separator + "combinedDBSearchResults.reformat.sorted.tsv";
		
		
		String headerLine = Utils.concatenateFilesNoHeaders(dbSearchOutputFiles, combinedDBSearchOutputFile, '#');
		ArrayList<String> sedCmd = new ArrayList<String>();
		sedCmd.add("sed");
		sedCmd.add("s/index=//");
		sedCmd.add(combinedDBSearchOutputFile);
		//String cmd = "sed s/index=// " + combinedDBSearchOutputFile; 
		try {
			Utils.RunCommandWithStdOutToFile(sedCmd, temp, null, this.projectDir);
			//Utils.RunCommand(cmd, null, projectDir);
		} catch (Exception e) {
			Utils.writeStatusToFile(this.statusFileName, "Error");
			e.printStackTrace();
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to run sed command to process MSGF+ results");
		}
		
		//cmd = "sort -k 2n,2 -o " + sortedDBSearchFile;
		String cmd = "sort -k 2n,2 -o " + temp2 + " " + temp;
		try {
			Utils.RunCommand(cmd, null, projectDir);
		} catch (Exception e) {
			Utils.writeStatusToFile(this.statusFileName, "Error");
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to run sort command to process MSGF+ results");
		}
		
		Utils.addHeaderLineToFile(headerLine, temp2, this.sortedDBSearchFile);
	}

	/**
	 * Assumes the project directory itself already exists and is writeable
	 * 1.  Creates the subdirectories
	 *  - msgfplus
	 *  - ValensModData
	 *  - ValensModResults
	 * 2.  Creates pklbin files for input files
	 */
	private void setupProjectDirectory() {
		
		System.out.println("[0] Setting up the project directory");
		this.msgfplusOutputDir = this.projectDir + File.separator + "msgfplus";
		this.reportDataDir = this.projectDir + File.separator + "ValensModData";
		this.reportHTMLDir = this.projectDir + File.separator + "ValensModResults";
		
		Utils.MakeDir(this.msgfplusOutputDir);
		Utils.MakeDir(this.reportDataDir);
		Utils.MakeDir(this.reportHTMLDir);
		
		this.pklbinFileNames = new String[this.spectrumFileNames.length];
		for (int i = 0; i < this.spectrumFileNames.length; ++i) {
			this.pklbinFileNames[i] = Utils.GetFileNameNoExtension(this.spectrumFileNames[i]) + ".pklbin";
			if(this.fastMode && Utils.IsFile(this.pklbinFileNames[i])) {
				System.out.println(" - Running in fast mode, skipping conversion to pklbin because files exist!");
				continue;
			}
			System.out.println(" - converting " + this.spectrumFileNames[i] + " -> " + this.pklbinFileNames[i]);
			if(!PKLBINUtils.convertFileToPKLBIN(i, this.spectrumFileNames[i], this.pklbinFileNames[i], this.frag[i], this.inst)) {
				Utils.writeStatusToFile(this.statusFileName, "Error");
				ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR, "Unable to convert to pklbin");
			}
			
		}
		
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String[] options = { "-i","-f"};
		boolean[] values = { true,false};
		Hashtable<String,String> commandLineArgs = Utils.ParseCommandLine(args, options,
				values);

		if(!commandLineArgs.containsKey("-i"))
		{
			System.err.println("ERROR: Must specify an input parameter file");
			System.err.println(usageInfo);
			System.exit(1);
		}
		
		runValensMod runner = new runValensMod(commandLineArgs.get("-i"));
		runner.fastMode = commandLineArgs.containsKey("-f");
		runner.run();
		
		runner.cleanupDirectories();
		Utils.writeStatusToFile(runner.statusFileName, "Finished");

	}

	private void cleanupDirectories() {
		System.out.println("[4] Cleaning up project directory");
		//TODO fill this in
		
	}

}
