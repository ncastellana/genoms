package antibody;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;

import basicUtils.Utils;
import dbUtils.ProteinDB;
import errorUtils.ErrorThrower;

public class FragmentDB {

	
	/**
	 * Returns the class id from the type number
	 * 0 -> 1, 1-> 2, 2 -> 3, 3 -> 4 (heavy chain)
	 * 4 -> 1, 5 -> 2, 6 -> 3 (light chain)
	 * @param type
	 * @return
	 */
	public static int getClassID(int type) {
		if(type == 0)
			return 1; //HC V
		if(type == 1)
			return 2; //HC D
		if(type == 2)
			return 3; //HC J
		if(type == 3)
			return 4; //HC C
		if(type == 4)
			return 1; //LC V
		if(type == 5)
			return 2; //LC J
		if(type == 6)
			return 3; //LC C
		
		else { 
			ErrorThrower.ThrowWarningCustom(ErrorThrower.CUSTOM_WARNING, "Invalid template type '" + type + "'!");
			return -1;
		}
	}
	
	/**
	 * Returns the name of the chain for the template type (0-3) -> HC, (4-6) -> LC
	 * @param type
	 * @return Returns HC or LC or null
	 */
	public static String getSourceName(int type) {
		if(type >= 0 && type < 4)
			return "HC";
		if(type >= 4 && type < 7)
			return "LC";
		else { 
			ErrorThrower.ThrowWarningCustom(ErrorThrower.CUSTOM_WARNING, "Invalid template type '" + type + "'!");
			return null;
		}
	}
	
	/**
	 * Augments/creates a valens database with the given name
	 * @param geneFileName
	 * @param type
	 * @param baseFileName
	 * @param newDBFileName
	 */
	public static void augment(String geneFileName, int type, String baseFileName, String newDBFileName) {
		
		Template[] baseTemps = null;
		System.out.println("New templates are from " + geneFileName);
		
		int templateCount = 0;
		//Remove the file just in case
		Utils.DeleteFile(newDBFileName);
		int currProteinID = 0;
		
		//Load the constraints
		if(baseFileName != null) {
			System.out.println(" - loading old templates from " + FragmentDB.getConstraintFileName(baseFileName));
			baseTemps = FragmentDB.LoadAllTemplateFromConstraints(FragmentDB.getConstraintFileName(baseFileName));
			for(int i = 0; i < baseTemps.length; ++ i) {
				if(baseTemps[i].ProteinID >= currProteinID)
					currProteinID = baseTemps[i].ProteinID + 1;
			}
			templateCount += baseTemps.length;
			//Write the proteins to the new DB FileName
			Utils.copyFile(baseFileName, newDBFileName);
		}
		else
			System.out.println(" - no base database to load");
		
		FileWriter seqFile = Utils.openFileWriter(newDBFileName,true);
		
		ArrayList<Template> newTemps = new ArrayList<Template>();
		int classID = FragmentDB.getClassID(type);
		String source = FragmentDB.getSourceName(type);
		
		ProteinDB tempDB = new ProteinDB(geneFileName);
		System.out.println(" - Adding " + tempDB.getProteinCount() + " templates");
		for(int i = 0; i < tempDB.getProteinCount(); ++i) {
			Template t = new Template(currProteinID);
			t.classID = classID;
			t.SourceDB = source;
			newTemps.add(t);
			
			//Add the constraints around this new template
			for(int j = 0; baseTemps != null && j < baseTemps.length; ++j) {
				if(t.SourceDB.equalsIgnoreCase(baseTemps[j].SourceDB)) {
					//if(t.classID == baseTemps[j].classID) 
					//	t.AddForbiddenEdge(baseTemps[j]);
					if(t.classID == baseTemps[j].classID - 1)
						t.addForwardEdge(baseTemps[j]);
					else if(t.classID == baseTemps[j].classID + 1)
						baseTemps[j].addForwardEdge(t);
				}
			}
			
			
			//Write the template sequence to the file
			String header = ">" + t.classID + AntibodyUtils.CLASS_DELIM 
					+ t.SourceDB + AntibodyUtils.CLASS_DELIM + t.ProteinID + "\n";
			Utils.writeLine(seqFile, newDBFileName, header);
			String seq = tempDB.getProteinSequence(i);
			Utils.writeLine(seqFile, newDBFileName, seq + "\n");
			
			
			currProteinID += 1;
		}
		templateCount += newTemps.size();
		Template[] allTemplates = new Template[templateCount];
		int tempIdx = 0;
		if(baseTemps != null)
		{
			for(int i = 0; i < baseTemps.length; ++i) {
				allTemplates[tempIdx] = baseTemps[i];
				tempIdx += 1;
			}
		}
		for(int i = 0; i < newTemps.size(); ++ i) {
			allTemplates[tempIdx] = newTemps.get(i);
			tempIdx += 1;
		}
			
		FragmentDB.WriteConstraintsToFileBinary(allTemplates, FragmentDB.getConstraintFileName(newDBFileName));
		
		
		Utils.closeFileWriter(seqFile, newDBFileName);
	}
	
	/**
	 * Augments/creates a valens database with the given name
	 * @param geneFileName
	 * @param type
	 * @param baseFileName
	 * @param newDBFileName
	 */
	public static void augment(String[] geneFileName, int[] type, String baseFileName, String newDBFileName) {
		
		Template[] baseTemps = null;
		System.out.println("New templates are from " + Utils.JoinStringArray(geneFileName, ","));
		
		int templateCount = 0;
		//Remove the file just in case
		Utils.DeleteFile(newDBFileName);
		int currProteinID = 0;
		
		//Load the constraints
		if(baseFileName != null) {
			System.out.println(" - loading old templates from " + FragmentDB.getConstraintFileName(baseFileName));
			baseTemps = FragmentDB.LoadAllTemplateFromConstraints(FragmentDB.getConstraintFileName(baseFileName));
			for(int i = 0; i < baseTemps.length; ++ i) {
				if(baseTemps[i].ProteinID >= currProteinID)
					currProteinID = baseTemps[i].ProteinID + 1;
			}
			templateCount += baseTemps.length;
			//Write the proteins to the new DB FileName
			Utils.copyFile(baseFileName, newDBFileName);
		}
		else
			System.out.println(" - no base database to load");
		
		FileWriter seqFile = Utils.openFileWriter(newDBFileName,true);
		
		ArrayList<Template> newTemps = new ArrayList<Template>();
		for(int fIdx = 0; fIdx < geneFileName.length; ++fIdx) {
			int classID = FragmentDB.getClassID(type[fIdx]);
			String source = FragmentDB.getSourceName(type[fIdx]);
		
			if(geneFileName[fIdx] == null)
			{
				System.out.println(" - No new templates of type " + source + "." + classID + " to add");
				continue;
			}
			ProteinDB tempDB = new ProteinDB(geneFileName[fIdx]);
			System.out.println(" - Adding " + tempDB.getProteinCount() + " templates");
			for(int i = 0; i < tempDB.getProteinCount(); ++i) {
				Template t = new Template(currProteinID);
				t.classID = classID;
				t.SourceDB = source;
				newTemps.add(t);
			
				//Add the constraints around this new template
				for(int j = 0; baseTemps != null && j < baseTemps.length; ++j) {
					if(t.SourceDB.equalsIgnoreCase(baseTemps[j].SourceDB)) {
						//if(t.classID == baseTemps[j].classID) 
						//	t.AddForbiddenEdge(baseTemps[j]);
						if(t.classID == baseTemps[j].classID - 1)
							t.addForwardEdge(baseTemps[j]);
						else if(t.classID == baseTemps[j].classID + 1)
							baseTemps[j].addForwardEdge(t);
					}
				}
				
			
				//Write the template sequence to the file
				String header = ">" + t.classID + AntibodyUtils.CLASS_DELIM 
						+ t.SourceDB + AntibodyUtils.CLASS_DELIM + t.ProteinID + "\n";
				Utils.writeLine(seqFile, newDBFileName, header);
				String seq = tempDB.getProteinSequence(i);
				Utils.writeLine(seqFile, newDBFileName, seq + "\n");
			
			
				currProteinID += 1;
			}
			tempDB.destruct();
		}
		templateCount += newTemps.size();
		Template[] allTemplates = new Template[templateCount];
		int tempIdx = 0;
		if(baseTemps != null)
		{
			for(int i = 0; i < baseTemps.length; ++i) {
				allTemplates[tempIdx] = baseTemps[i];
				tempIdx += 1;
			}
		}
		for(int i = 0; i < newTemps.size(); ++ i) {
			allTemplates[tempIdx] = newTemps.get(i);
			tempIdx += 1;
		}
			
		FragmentDB.WriteConstraintsToFileBinary(allTemplates, FragmentDB.getConstraintFileName(newDBFileName));
		
		
		Utils.closeFileWriter(seqFile, newDBFileName);
	}
	
	public static boolean WriteConstraintsToFileBinary(Template[] Templates,
			String ConstraintFile) {

		System.out.println("Writing template constraints to " + ConstraintFile
				+ " for " + Templates.length + " templates ...");
		DataOutputStream f = null;
		try {
			f = new DataOutputStream(new FileOutputStream(ConstraintFile));
			// f.writeInt((contigCount));
			// f.writeInt(Integer.reverseBytes(contigCount));
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, ConstraintFile);
		}
		int ConstraintCount = 0;
		HashSet classes = new HashSet();

		// Write the constraints. There are three types
		// Information line: i pID cID (protein ID, classID)
		// Order line: o pID1 pID2 (protein ID 1 and protine ID2)
		// Forbidden line: f pID1 pID2 (protein ID1 and protein ID2)
		
		int idConstraints = 0;
		int fConstraints = 0;
		int oConstraints = 0;
		for (int i = 0; i < Templates.length; ++i) {
			// String FileLines = "";
			// FileLines += "i\t" + Templates[i].ProteinID + "\t" +
			// Templates[i].classID + "\n";
			if(Templates[i] == null)
				System.out.println("Template " + i + " is null");
			ConstraintCount += 1;
			idConstraints += 1;
			classes.add(Templates[i].SourceDB);
			
			for (int n = 0; Templates[i].NextTemplates != null
					&& n < Templates[i].NextTemplates.size(); ++n) {
				Template Next = (Template) (Templates[i].NextTemplates.get(n));
				// FileLines += "o\t" + Templates[i].ProteinID + "\t"
				// + Next.ProteinID + "\n";
				ConstraintCount += 1;
				oConstraints += 1;
			}
			for (int n = 0; Templates[i].ForbiddenTemplates != null
					&& n < Templates[i].ForbiddenTemplates.size(); ++n) {
				Template Next = (Template) (Templates[i].ForbiddenTemplates
						.get(n));
				if (Next.ProteinID > Templates[i].ProteinID) {
					// FileLines += "f\t" + Templates[i].ProteinID + "\t"
					// + Next.ProteinID + "\n";
					ConstraintCount += 1;
					fConstraints += 1;
				}
			}
		}
		
		//System.out.println("identity constraints: " + idConstraints);
		//System.out.println("forbidden constraints: " + fConstraints);
		//System.out.println("order constraints: " + oConstraints);
		
		String[] classNames = new String[classes.size()];
		Iterator<String> it = classes.iterator();
		int idx = 0;
		while(it.hasNext()) {
			classNames[idx] = it.next();
			idx += 1;
		}
			
		try {
			//write the constraints count
			f.writeInt(Integer.reverseBytes(ConstraintCount));
			//write the number of classes
			f.writeInt(Integer.reverseBytes(classes.size()));
			
			for(int i = 0; i < classNames.length; ++i)
				f.writeInt(Integer.reverseBytes(classNames[i].length()));
			
				
			
		} catch (IOException E) {
			ErrorThrower.ThrowError(7, ConstraintFile);
		}
		
		for(int i = 0; i < classNames.length; ++i)
			Utils.writeStringDataOutputStream(classNames[i], f, ConstraintFile);
		
		for (int i = 0; i < Templates.length; ++i) {
			if (i % 1000 == 0)
				System.out.println("Writing constraints for template " + i
						+ "/" + Templates.length);
			// String FileLines = "";
			// FileLines += "i\t" + Templates[i].ProteinID + "\t" +
			// Templates[i].classID + "\n";
			int classIdx = Utils.FindStringInArray(classNames, Templates[i].SourceDB);
			if(classIdx < 0 || classIdx >= classNames.length) {
				System.out.println("UNUSUAL CLASS NAME: " + Templates[i].SourceDB);
				//Utils.WaitForEnter();
			}
			try {
				f.writeChar(Character.reverseBytes('i'));
				f.writeInt(Integer.reverseBytes(Templates[i].ProteinID));
				f.writeInt(Integer.reverseBytes(Templates[i].classID));
				f.writeInt(Integer.reverseBytes(classIdx));
			} catch (IOException E) {
				ErrorThrower.ThrowError(7, ConstraintFile);
			}
			for (int n = 0; Templates[i].NextTemplates != null
					&& n < Templates[i].NextTemplates.size(); ++n) {
				Template Next = (Template) (Templates[i].NextTemplates.get(n));
				try {
					f.writeChar(Character.reverseBytes('o'));
					f.writeInt(Integer.reverseBytes(Templates[i].ProteinID));
					f.writeInt(Integer.reverseBytes(Next.ProteinID));
				} catch (IOException E) {
					ErrorThrower.ThrowError(7, ConstraintFile);
				}
			}
			for (int n = 0; Templates[i].ForbiddenTemplates != null
					&& n < Templates[i].ForbiddenTemplates.size(); ++n) {
				Template Next = (Template) (Templates[i].ForbiddenTemplates
						.get(n));
				if (Next.ProteinID > Templates[i].ProteinID) {
					try {
						f.writeChar(Character.reverseBytes('f'));
						f.writeInt(Integer.reverseBytes(Templates[i].ProteinID));
						f.writeInt(Integer.reverseBytes(Next.ProteinID));
					} catch (IOException E) {
						ErrorThrower.ThrowError(7, ConstraintFile);
					}
				}
			}
		}

		System.out.println("Wrote " + ConstraintCount + " constraints to "
				+ ConstraintFile);
		try {
			f.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(8, ConstraintFile);
		}

		return true;
	}
	
	/**
	 * Given a Valens database name, returns the constraints file name (assumed to be in the same directory as the database
	 * @param baseFileName
	 * @return
	 */
	public static String getConstraintFileName(String baseFileName) {
		
		String constraintFileName = Utils.GetFileNameNoExtension(baseFileName) + ".constraints.bin";
		
		return constraintFileName;
	}
	
	/**
	 * Loads the set of templates from a constraints file.  The sequences of the templates is not loaded, but the class
	 * numbers and the constraints are loaded
	 * @param FileName
	 * @return
	 */
	public static Template[] LoadAllTemplateFromConstraints(String FileName) {
		boolean LocalDebug = false;
		RandomAccessFile rIn = null;

		try {
			rIn = new RandomAccessFile(FileName, "r");

		} catch (IOException E) {
			ErrorThrower.ThrowError(5, FileName);
		}

		Hashtable<Integer, Template> ProteinIDHash = new Hashtable<Integer, Template>();

		String Line = "";
		HashSet<Integer> classes = new HashSet<Integer>();
		// Read the number of constraints (n)
		// Read the number of sources (m)
		// Read the lengths of the source names 
		int numConstraints = 0;
		int numSources = 0;
		boolean localDebug = false;
		
		int[] sourceNameLengths = null;
		String[] sourceNames = null;
		
		try {
			numConstraints = Integer.reverseBytes(rIn.readInt());
			numSources = Integer.reverseBytes(rIn.readInt());
			
		} catch (Exception E) {
			ErrorThrower.ThrowError(6, FileName);
		}
		if(localDebug)
			System.out.println("COnstraints: " + numConstraints + ", numSources: " + numSources);
		if(numConstraints < 0) 
			ErrorThrower.ThrowError(6, "Invalid number of constraints '" + numConstraints + "':" + FileName);
		if(numSources < 0) 
			ErrorThrower.ThrowError(6, "Invalid number of sources '" + numSources + "':" + FileName);
		
		sourceNameLengths = Utils.readIntArrayRandomAccessFile(rIn, FileName, numSources);
		sourceNames = new String[numSources];
		
		for(int i = 0; i < numSources; ++i) {
			sourceNameLengths[i] = Integer.reverseBytes(sourceNameLengths[i]);
			sourceNames[i] = Utils.readStringRandomAccessFile(rIn, FileName, sourceNameLengths[i]);
			if(localDebug)
				System.out.println("source[" + i + "]:" + sourceNameLengths[i] + ", " + sourceNames[i]);
		}
		
		
		for (int i = 0; i < numConstraints; ++i) {
			char type = 'x';
			int num1 = -1;
			int num2 = -1;
			int num3 = -1;
			try {
				type = Character.reverseBytes(rIn.readChar());
				num1 = Integer.reverseBytes(rIn.readInt());
				num2 = Integer.reverseBytes(rIn.readInt());
			} catch (Exception E) {
				ErrorThrower.ThrowError(6, FileName);
			}

			if (type == 'I' || type == 'i') {
				
				num3 = Integer.reverseBytes(Utils.readIntRandomAccessFile(rIn, FileName));
				int proteinID = num1;
				int classID = num2;
				if(num3 >= sourceNames.length)
					System.out.println("Weird line: " + num1 + " " + num2 + " " + num3);
				String sourceName = sourceNames[num3];
				
				classes.add(new Integer(classID));

				Template T1;
				if (ProteinIDHash.containsKey(new Integer(proteinID))) {
					T1 = (Template) ProteinIDHash.get(new Integer(proteinID));
				} else {
					T1 = new Template(proteinID);
				}

				T1.classID = classID;
				T1.SourceDB = sourceName;
				T1.ProteinID = proteinID;
				
				ProteinIDHash.put(new Integer(proteinID), T1);
			} else if (type == 'O' || type == 'o') {
				Template T1;
				Template T2;
				if (ProteinIDHash.containsKey(new Integer(num1))) {
					T1 = (Template) ProteinIDHash.get(new Integer(num1));
				} else {
					T1 = new Template(num1);
					ProteinIDHash.put(new Integer(num1), T1);
				}
				if (ProteinIDHash.containsKey(new Integer(num2))) {
					T2 = (Template) ProteinIDHash.get(new Integer(num2));
				} else {
					T2 = new Template(num2);
					ProteinIDHash.put(new Integer(num2), T2);
				}
				if (T1.NextTemplates == null) {
					T1.NextTemplates = new ArrayList<Template>();
				}
				boolean Found = false;
				for (int j = 0; j < T1.NextTemplates.size(); ++j) {
					Template CurrTemplate = (Template) (T1.NextTemplates.get(j));
					if (CurrTemplate.ProteinID == num2) {
						if (LocalDebug) {
							System.out
									.println("This relationship has laready been established");
							Utils.WaitForEnter();
						}
						Found = true;
						break;
					}
				}
				if (!Found) {

					T1.NextTemplates.add(T2);
					if (T2.PrevTemplates == null) {
						T2.PrevTemplates = new ArrayList<Template>();
					}
					T2.PrevTemplates.add(T1);
					if (LocalDebug) {
						System.out.println("This relationship is new!!");
						T1.DebugPrint();
						T2.DebugPrint();
						Utils.WaitForEnter();
					}
				}
			} else if (type == 'F' || type == 'f') {
				Template T1;
				Template T2;
				if (ProteinIDHash.containsKey(new Integer(num1))) {
					T1 = (Template) ProteinIDHash.get(new Integer(num1));
				} else {
					T1 = new Template(num1);
					ProteinIDHash.put(new Integer(num1), T1);
				}
				if (ProteinIDHash.containsKey(new Integer(num2))) {
					T2 = (Template) ProteinIDHash.get(new Integer(num2));
				} else {
					T2 = new Template(num2);
					ProteinIDHash.put(new Integer(num2), T2);
				}
				if (T1.ForbiddenTemplates == null) {
					T1.ForbiddenTemplates = new ArrayList<Template>();
				}
				boolean Found = false;
				for (int j = 0; j < T1.ForbiddenTemplates.size(); ++j) {
					Template CurrTemplate = (Template) (T1.ForbiddenTemplates
							.get(j));
					if (CurrTemplate.ProteinID == num2) {
						if (LocalDebug) {
							System.out
									.println("This relationship has laready been established");
							Utils.WaitForEnter();
						}
						Found = true;
						break;
					}
				}
				if (!Found) {
					T1.ForbiddenTemplates.add(T2);
					if (T2.ForbiddenTemplates == null) {
						T2.ForbiddenTemplates = new ArrayList<Template>();
					}
					T2.ForbiddenTemplates.add(T1);
					if (LocalDebug) {
						System.out.println("This relationship is new!!");
						T1.DebugPrint();
						T2.DebugPrint();
						Utils.WaitForEnter();
					}
				}
			} else {
				System.err.println("ERROR: Unacceptable character in " + Line);
			}
		}

		Utils.closeRandomAccessFile(rIn, FileName);
		
		//System.out.println("Total classes: " + classes.size());
		Template[] Ret = new Template[ProteinIDHash.size()];
		Enumeration<Integer> Keys = ProteinIDHash.keys();
		int Count = 0;
		while (Keys.hasMoreElements()) {
			Ret[Count] = (Template) (ProteinIDHash.get(Keys.nextElement()));
			Count += 1;
		}
		System.out.println("Loaded " + Count + " Templates with constraints");
		return Ret;
	}
	
	
	public static void dumpConstraints(String constraintFile, String humanReadableFile) {
		RandomAccessFile rIn = null;
		FileWriter f = Utils.openFileWriter(humanReadableFile);
		boolean LocalDebug = true;

		try {
			rIn = new RandomAccessFile(constraintFile, "r");

		} catch (IOException E) {
			ErrorThrower.ThrowError(5, constraintFile);
		}

		String Line = "";
		// Read the number of constraints (n)
		// Read the number of sources (m)
		// Read the lengths of the source names 
		int numConstraints = 0;
		int numSources = 0;
		boolean localDebug = false;
		
		int[] sourceNameLengths = null;
		String[] sourceNames = null;
		
		try {
			numConstraints = Integer.reverseBytes(rIn.readInt());
			numSources = Integer.reverseBytes(rIn.readInt());
			
		} catch (Exception E) {
			ErrorThrower.ThrowError(6, constraintFile);
		}
		
		Utils.writeLine(f, humanReadableFile, "constraints: " + numConstraints + "\n");
		Utils.writeLine(f, humanReadableFile, "num sources: " + numSources + "\n");;
		if(localDebug)
			System.out.println("COnstraints: " + numConstraints + ", numSources: " + numSources);
		if(numConstraints < 0) 
			ErrorThrower.ThrowError(6, "Invalid number of constraints '" + numConstraints + "':" + constraintFile);
		if(numSources < 0) 
			ErrorThrower.ThrowError(6, "Invalid number of sources '" + numSources + "':" + constraintFile);
		
		sourceNameLengths = Utils.readIntArrayRandomAccessFile(rIn, constraintFile, numSources);
		sourceNames = new String[numSources];
		
		for(int i = 0; i < numSources; ++i) {
			sourceNameLengths[i] = Integer.reverseBytes(sourceNameLengths[i]);
			sourceNames[i] = Utils.readStringRandomAccessFile(rIn, constraintFile, sourceNameLengths[i]);
			if(localDebug)
				System.out.println("source[" + i + "]:" + sourceNameLengths[i] + ", " + sourceNames[i]);
		}
		
		Utils.writeLine(f, humanReadableFile, Utils.JoinStringArray(sourceNames, ",") + "\n");
		
		for (int i = 0; i < numConstraints; ++i) {
			char type = 'x';
			int num1 = -1;
			int num2 = -1;
			int num3 = -1;
			try {
				type = Character.reverseBytes(rIn.readChar());
				num1 = Integer.reverseBytes(rIn.readInt());
				num2 = Integer.reverseBytes(rIn.readInt());
			} catch (Exception E) {
				ErrorThrower.ThrowError(6, constraintFile);
			}

			if (type == 'I' || type == 'i') {
				
				num3 = Integer.reverseBytes(Utils.readIntRandomAccessFile(rIn, constraintFile));
				int proteinID = num1;
				int classID = num2;
				String sourceName = sourceNames[num3];
				
				Utils.writeLine(f, humanReadableFile, "i " + num1 + " " + num2 + " " + num3 + "\n");
			} else if (type == 'O' || type == 'o') {
				Utils.writeLine(f, humanReadableFile, "o " + num1 + " " + num2 + "\n");
			} else if (type == 'F' || type == 'f') {
				Utils.writeLine(f, humanReadableFile, "f " + num1 + " " + num2 + "\n");
			} else {
				System.err.println("ERROR: Unacceptable character in " + Line);
			}
		}

		Utils.closeRandomAccessFile(rIn, constraintFile);
		Utils.closeFileWriter(f, humanReadableFile);
		
	}

	/**
	 * Creates a new Valens database from the db roots.  Simply updates the names of the sequences t oinclude the 
	 * class number and the rootname, and concatenates the sequences to a single file.
	 * @param projectDir
	 * @param DatabaseRootName
	 * @param Debug
	 * @return
	 */
	public static String createFASTADBandConstraintsFromRootName(
			String projectDir, String[] DatabaseRootName, boolean Debug) {
		BufferedReader Reader = null;
		String CombinedDBFasta = projectDir + File.separator;
		for (int i = 0; i < DatabaseRootName.length; ++i) {

			if (DatabaseRootName[i].charAt(0) != File.separatorChar)
				DatabaseRootName[i] = AntibodyUtils.makePathAbsolute(
						DatabaseRootName[i], projectDir);
			CombinedDBFasta += Utils.GetBaseName(DatabaseRootName[i]) + "_";
		}
		CombinedDBFasta += "combined.fa";
		BufferedWriter Writer = null;
		try {
			Writer = new BufferedWriter(new FileWriter(CombinedDBFasta));
		} catch (IOException E) {
			ErrorThrower.ThrowError(5, CombinedDBFasta);
		}
		if (Debug)
			System.out.println("Attempting to write: " + CombinedDBFasta);

		
		//Iterate over the different roots (e.g. chains)
		for (int j = 0; j < DatabaseRootName.length; ++j) {
			int DBIndex = 1; // Database index

			while (true) {
				String CurrFile = DatabaseRootName[j]
						+ AntibodyUtils.CLASS_DELIM + DBIndex;
				System.out.println("Looking for file " + CurrFile);
				try {
					Reader = new BufferedReader(new FileReader(CurrFile));
				} catch (IOException E) {
					// System.err.println(E.getMessage());
					break;
				}
				if (Debug)
					System.out.println("Loading DB " + CurrFile);
				String Line = null;
				try {
					Line = Reader.readLine();
				} catch (IOException E) {
					ErrorThrower.ThrowError(6, CurrFile);
				}
				while (Line != null) {
					// Write the line
					if (Line.length() > 0 && Line.charAt(0) == '>')
						Line = ">"
								+ DBIndex
								+ AntibodyUtils.CLASS_DELIM
								+ Utils.GetBaseNameNoExtension(DatabaseRootName[j])
								+ AntibodyUtils.CLASS_DELIM
								+ Line.substring(1, Line.length());
					else // This is a sequence line
					{
						// Make sure all letters are amino acids
						Line = Utils.AminoAcidify(Line);
					}
					if (Line.length() > 0) {
						try {
							Writer.write(Line);
							Writer.newLine();
						} catch (IOException E) {
							ErrorThrower.ThrowError(7, CombinedDBFasta);
						}
					}
					// Read another line
					try {
						Line = Reader.readLine();
					} catch (IOException E) {
						ErrorThrower.ThrowError(6, CurrFile);
					}

				}
				if (Line == null) {
					try {
						Reader.close();
					} catch (IOException E) {
						ErrorThrower.ThrowError(8, CurrFile);
					}
				}
				DBIndex += 1;
			}
		}
		try {
			Writer.close();
		} catch (IOException E) {
			ErrorThrower.ThrowError(8, CombinedDBFasta);
		}
		return CombinedDBFasta;

	}

}
