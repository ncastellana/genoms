package antibody;

import jabgraph.PartitionReads;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;

import errorUtils.ErrorThrower;

import basicUtils.Utils;


/**
 * This is the entry class for running Valens-DB, a tool for creating antibody sequence
 * databases for use with Valens.
 * Created 2015-08-17
 * @author natalie
 *
 */
public class runValensDB {

	
	public static String usageInfo = "antibody.runValensDB version 20150827.2\n" +
								"Options:\n" +
								" -i [FILE] File containing a list of files to be combined.  Each line is of the form\n" +
								"                  type filename (att1) (att2)\n" +
								"           where 'type' may be 'ref', 'base', 'seq' or 'template'.  Only one base database may be specified.  Sequence files\n" +
								"           should be preprocessed reads.  Template files contain addition columns containing the chain and the gene\n" +
								" -o [FILE] File name of the database to be created from merging the files in the file list\n";
	
	
	public static int K = 11;
	public static int TOTAL_CLASSES = 7;
	private String[] seqFileNames = null;
	private String[][] templateFileNames = null;
	private String[][] referenceFileNames = null;
	private String baseFileName = null;
	
	private String inputFileList = null;
	
	
	private String outputFileName;


	private String statusFileName;
	
	public runValensDB(String inputFileList, String outputFileName) {
		ArrayList<String> tempSeqFileNames = new ArrayList<String>();
		ArrayList<String[]> tempTemplateFileNames = new ArrayList<String[]>();
		ArrayList<String[]> tempReferenceFileNames = new ArrayList<String[]>();
		
		BufferedReader input = Utils.openBufferedReader(inputFileList);
		String line = Utils.readNextLine(input, inputFileList);
		this.inputFileList = inputFileList;
		
		
		while(line != null) {
			line = line.trim();
			if(line.length() == 0 || line.charAt(0) == '#')
			{
				line = Utils.readNextLine(input, inputFileList);
				continue;
			}
			String[] bits = line.split("\t");
			
			//Check if its the base file name
			if(bits[0].equalsIgnoreCase("base") && this.baseFileName == null) {
				if(!Utils.IsFile(bits[1])) {
					Utils.writeStatusToFile(statusFileName,"Error");
					ErrorThrower.ThrowError(1, bits[1]);
				}
				this.baseFileName = bits[1];
			}
			//Check if its a sequence file
			else if(bits[0].equalsIgnoreCase("seq")) {
				if(!Utils.IsFile(bits[1])) {
					Utils.writeStatusToFile(statusFileName,"Error");
					ErrorThrower.ThrowError(1, bits[1]);
				}
				tempSeqFileNames.add(bits[1]);
			}
			//Check if its a template file
			else if(bits[0].equalsIgnoreCase("template")) {
				if(bits.length != 4)
					ErrorThrower.ThrowError(12, "Template sequence lines must contain 4 columns: " + inputFileList);
				if(!Utils.IsFile(bits[1])) {
					Utils.writeStatusToFile(statusFileName,"Error");
					ErrorThrower.ThrowError(1, bits[1]);
				}
				String[] vals = new String[3];
				vals[0] = bits[1];
				vals[1] = bits[2].toLowerCase();
				if(vals[1].charAt(0) == 'h')
					vals[1] = "heavy";
				else
					vals[1] = "light";
				vals[2] = bits[3].toLowerCase().charAt(0) + "";
				
				tempTemplateFileNames.add(vals);
			}
			//Check if its a template file
			else if(bits[0].equalsIgnoreCase("ref")) {
				if(bits.length != 4) {
					Utils.writeStatusToFile(statusFileName,"Error");
					ErrorThrower.ThrowError(12, "Reference sequence lines must contain 4 columns: " + inputFileList);
				}if(!Utils.IsFile(bits[1])) {
					Utils.writeStatusToFile(statusFileName,"Error");
					ErrorThrower.ThrowError(1, bits[1]);
				}
				String[] vals = new String[3];
				vals[0] = bits[1];
				vals[1] = bits[2].toLowerCase();
				if(vals[1].charAt(0) == 'h')
					vals[1] = "heavy";
				else
					vals[1] = "light";
				vals[2] = bits[3].toLowerCase().charAt(0) + "";
				
				tempReferenceFileNames.add(vals);
			}
			else {
				
				ErrorThrower.ThrowWarning(12, "Ignoring file list line starting with '" + bits[0] + "'");
				
			}
			line = Utils.readNextLine(input, inputFileList);
		} //end while
		
		this.seqFileNames = Utils.ConvertArraylistToStringArray(tempSeqFileNames);
		
		this.templateFileNames = new String[tempTemplateFileNames.size()][3];
		for(int i = 0; i < this.templateFileNames.length; ++i) 
			this.templateFileNames[i] = tempTemplateFileNames.get(i);
		
		this.referenceFileNames = new String[tempReferenceFileNames.size()][3];
		for(int i = 0; i < this.referenceFileNames.length; ++i) 
			this.referenceFileNames[i] = tempReferenceFileNames.get(i);
		
		
		this.outputFileName = outputFileName;
		
		this.statusFileName = this.outputFileName + ".status.txt";
		//System.out.println("status file: " + this.statusFileName);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String[] options = { "-i", "-o"};
		boolean[] values = { true, true};
		Hashtable<String,String> commandLineArgs = Utils.ParseCommandLine(args, options,
				values);

		if(!commandLineArgs.containsKey("-i") && !commandLineArgs.containsKey("-o"))
		{
			System.err.println("ERROR: Must specify an input file list and an output file name");
			System.err.println(usageInfo);
			
			System.exit(0);
		}
		
		runValensDB runner = new runValensDB(commandLineArgs.get("-i"),commandLineArgs.get("-o"));
		
		runner.run();
		
	}

	private void run() {
		
		
		Utils.writeStatusToFile(statusFileName,"Running");
		
		String tempOutputDir = Utils.GetFileNameNoExtension(this.outputFileName);
		Utils.DeleteFile(tempOutputDir);
		Utils.MakeDir(tempOutputDir);
		
		//String currBaseFileName = this.baseFileName;
		/*if(this.baseFileName != null){
			currBaseFileName = tempOutputDir + File.separator + "tempBaseDB.fasta";
			
			Utils.copyFile(this.baseFileName, currBaseFileName);
			Utils.copyFile(FragmentDB.getConstraintFileName(this.baseFileName), FragmentDB.getConstraintFileName(currBaseFileName));
			
		}*/
		
		//String currOutputFileName = tempOutputDir + File.separator + "tempDB.fasta";

		
		String[] partitionedReadsAA = null;
		if(this.seqFileNames != null && this.seqFileNames.length > 0) {
			System.out.println("Processing sequencing reads...");
			partitionedReadsAA = this.processSequencingReads(statusFileName,tempOutputDir);
			
			for(int i = 0; i < partitionedReadsAA.length; ++i)
				System.out.println(" - ProcessedReads: " + partitionedReadsAA[i]);
		}
		else
			System.out.println("No sequencing reads found...");
		
		//String[] finalSeqs = new String[7];
		String[] geneFileNames = new String[TOTAL_CLASSES];
		int[] types = new int[TOTAL_CLASSES];
		for(int type = 0; type < TOTAL_CLASSES; ++type) {
			ArrayList<String> fileNames = new ArrayList<String>();
			
			//See if we have any reads of this type
			if(partitionedReadsAA != null && partitionedReadsAA[type] != null)
				fileNames.add(partitionedReadsAA[type]);
			if(this.templateFileNames != null) {
				for(int i = 0; i < this.templateFileNames.length; ++i)
				{
					int fileType = getFileType(this.templateFileNames[i]);
					if(fileType == type)
						fileNames.add(this.templateFileNames[i][0]);
				}
			}
			
			//System.out.println("All files (" + fileNames.size() + ") : " + Utils.JoinStringArray(Utils.ConvertArraylistToStringArray(fileNames), ","));
			//String geneFileName = null;
			
			if(fileNames.size() > 1) {
				geneFileNames[type] = tempOutputDir + File.separator + type + ".fasta";
				Utils.concatenateFiles(Utils.ConvertArraylistToStringArray(fileNames), geneFileNames[type]);
			}
			else if(fileNames.size() == 1)
				geneFileNames[type] = fileNames.get(0);
			else
				geneFileNames[type] = null;
			types[type] = type;
		}
				
		//System.out.println("Adding genes for " + FragmentDB.getSourceName(type) + " " + FragmentDB.getClassID(type));
		FragmentDB.augment(geneFileNames, types, this.baseFileName, this.outputFileName);
		/*if(currBaseFileName == null)
			currBaseFileName = tempOutputDir + File.separator + "tempBaseDB.fasta";
				
				System.out.println("Finished adding the genes");
				//Utils.WaitForEnter();
				Utils.copyFile(currOutputFileName, currBaseFileName);
				Utils.copyFile(FragmentDB.getConstraintFileName(currOutputFileName), FragmentDB.getConstraintFileName(currBaseFileName));
				
			}
			
		}
		Utils.copyFile(currOutputFileName, this.outputFileName);
		Utils.copyFile(FragmentDB.getConstraintFileName(currOutputFileName), FragmentDB.getConstraintFileName(this.outputFileName));
		*/
		//FragmentDB.dumpConstraints(FragmentDB.getConstraintFileName(this.outputFileName), FragmentDB.getConstraintFileName(this.outputFileName) + ".txt");
			
		//Utils.WaitForEnter();
		
		Utils.DeleteFile(tempOutputDir);
		Utils.DeleteFile(this.inputFileList);
		System.out.println("Updating program status to 'Finished'");
		Utils.writeStatusToFile(statusFileName,"Finished");
	}
	
	/*private String getSeqCode(int type) {
		// TODO Auto-generated method stub
		return null;
	}*/

	/**
	 * returns the integer corresponding to the chain-gene segment combination
	 * @param strings
	 * @return
	 */
	private int getFileType(String[] strings) {
		if(strings[1].equalsIgnoreCase("heavy")) {
			if(strings[2].equalsIgnoreCase("v")) return 0;
			if(strings[2].equalsIgnoreCase("d")) return 1;
			if(strings[2].equalsIgnoreCase("j")) return 2;
			if(strings[2].equalsIgnoreCase("c")) return 3;
		}
		else {
			if(strings[2].equalsIgnoreCase("v")) return 4;
			if(strings[2].equalsIgnoreCase("j")) return 5;
			if(strings[2].equalsIgnoreCase("c")) return 6;
		}
		return -1;
	}

	private String[] processSequencingReads(String statusFileName, String tempOutputDir) {
		
		
		//We assume that the reference files are in the correct order in the file
		ArrayList<String> args = new ArrayList<String>();
		
		//ArrayList<String> partitionedReads = new ArrayList<String>();
		String[] partitionedReads = new String[TOTAL_CLASSES];
		partitionedReads = Utils.initializeStringArray(partitionedReads, null);
		for(int i = 0; i < this.referenceFileNames.length; ++i) {
			args.add("-f"); //add the reference sequences
			args.add(this.referenceFileNames[i][0]);
			int type = this.getFileType(this.referenceFileNames[i]);
			partitionedReads[type] = tempOutputDir + File.separator + Utils.GetBaseName(this.referenceFileNames[i][0]);
			//partitionedReads.add(tempOutputDir + File.separator + Utils.GetBaseName(this.referenceFileNames[i][0]));
		}
		
		//Add the kmer length
		args.add("-k");
		args.add(runValensDB.K + "");
		
		//output directory
		args.add("--out");
		
		args.add(tempOutputDir);
		
		//input reads
		String tempReadFileName = tempOutputDir + File.separator + "allReads.txt";
		if(this.seqFileNames.length > 1)
			Utils.concatenateFiles(this.seqFileNames,tempReadFileName);
		else
			tempReadFileName = this.seqFileNames[0];
		args.add("--reads");
		args.add(tempReadFileName);
		
		
		//First we must create templates from the RNA file
		try {
			PartitionReads.main(Utils.ConvertArraylistToStringArray(args));
		} catch (Exception e) {
			Utils.writeStatusToFile(statusFileName,"Error");
			return null;
		}
		
		String[] partitionedReadsAA = new String[TOTAL_CLASSES];
		for(int i = 0; i < partitionedReads.length; ++i) {
			if(partitionedReads[i] == null) {
				partitionedReadsAA[i] = null;
				continue;
			}
			
			System.out.println(" - translating partitioned read file " + Utils.GetBaseName(partitionedReads[i]));
			String newFileName = Utils.GetFileNameNoExtension(partitionedReads[i]) + ".aa.fasta";
			Utils.translateSeqFile(partitionedReads[i],newFileName, 8);
			if(!Utils.IsFile(newFileName))
				System.out.println("FILE DOESN'T EXISTS!");
			partitionedReadsAA[i] = newFileName;
		}
		
		return partitionedReadsAA;
	}

}
