package antibody;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;

import trieUtils.TrieDB;
import basicUtils.Utils;
import errorUtils.ErrorThrower;

public class Template {

	public String TemplateName = null;
	public String SourceDB = null;

	public Seed anchors = null;
	public Seed midPoint = null;

	public String Sequence = null;
	public int ProteinID = -1;

	public int classID = -1;

	// public Template[] NextTemplates;
	// public Template[] PrevTemplates;
	public ArrayList<Template> NextTemplates = null;
	public ArrayList<Template> PrevTemplates = null;

	// public Template[] ForbiddenTemplates;
	public ArrayList<Template> ForbiddenTemplates = null;

	// Used when searching multiple paths, to mark that a template is used in a
	// path, or is forbidden with a template in the path
	boolean SelectedInPath = false;

	// For creating a linked list
	public Template NextInList = null;
	public Template PrevInList = null;

	public Template(int ProteinID) {
		this.ProteinID = ProteinID;
	}

	/**
	 * Loads a template constraint file, creating edges, etc. Constraint takes
	 * the form of edges of either a ordering or forbidding type Format: O
	 * ProteinID1 ProteinID2 F ProteinID1 ProteinID2
	 * 
	 * @param FileName
	 * @return
	 */
	/*public static Template[] LoadAllTemplateConstraints(String FileName) {
		boolean LocalDebug = false;
		BufferedReader Buf = Utils.openBufferedReader(FileName);

		Hashtable<Integer, Template> ProteinIDHash = new Hashtable<Integer, Template>();
		// ArrayList CurrTemplates = new ArrayList();

		String Line = "";
		try {
			Line = Buf.readLine();
		} catch (IOException E) {
			E.printStackTrace();
			return null;
		}
		while (Line != null) {
			Line = Line.trim();
			if (Line.length() == 0 || Line.charAt(0) == '#') {
				try {
					Line = Buf.readLine();
				} catch (IOException E) {
					E.printStackTrace();
					return null;
				}
				continue;
			}

			String[] Bits = Line.split("\t");
			if (Bits.length < 3) {
				System.err.println("Error parsing constraint file " + FileName
						+ ":" + Line);
				try {
					Line = Buf.readLine();
				} catch (IOException E) {
					E.printStackTrace();
					return null;
				}
				continue;
			}

			if (LocalDebug)
				System.out.println("Line: " + Line);

			if (Bits[0].compareTo("I") == 0 || Bits[0].compareTo("i") == 0) {
				int proteinID = Integer.parseInt(Bits[1]);
				int classID = Integer.parseInt(Bits[2]);

				Template T1;
				if (ProteinIDHash.containsKey(new Integer(proteinID))) {
					T1 = (Template) ProteinIDHash.get(new Integer(proteinID));
				} else {
					T1 = new Template(proteinID);
				}

				T1.classID = classID;
				ProteinIDHash.put(new Integer(proteinID), T1);

				try {
					Line = Buf.readLine();
				} catch (IOException E) {
					E.printStackTrace();
					return null;
				}
				continue;
			}

			int ProteinID1 = Integer.parseInt(Bits[1]);
			int ProteinID2 = Integer.parseInt(Bits[2]);

			Template T1;
			Template T2;
			if (ProteinIDHash.containsKey(new Integer(ProteinID1))) {
				T1 = (Template) ProteinIDHash.get(new Integer(ProteinID1));
			} else {
				T1 = new Template(ProteinID1);
				ProteinIDHash.put(new Integer(ProteinID1), T1);
			}
			if (ProteinIDHash.containsKey(new Integer(ProteinID2))) {
				T2 = (Template) ProteinIDHash.get(new Integer(ProteinID2));
			} else {
				T2 = new Template(ProteinID2);
				ProteinIDHash.put(new Integer(ProteinID2), T2);
			}

			// If the edge is of type, ordering
			if (Bits[0].compareTo("O") == 0 || Bits[0].compareTo("o") == 0) {
				if (T1.NextTemplates == null) {
					T1.NextTemplates = new ArrayList<Template>();
				}
				boolean Found = false;
				for (int i = 0; i < T1.NextTemplates.size(); ++i) {
					Template CurrTemplate = (Template) (T1.NextTemplates.get(i));
					if (CurrTemplate.ProteinID == ProteinID2) {
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

			} else if (Bits[0].compareTo("F") == 0
					|| Bits[0].compareTo("f") == 0) {
				if (T1.ForbiddenTemplates == null) {
					T1.ForbiddenTemplates = new ArrayList<Template>();
				}
				boolean Found = false;
				for (int i = 0; i < T1.ForbiddenTemplates.size(); ++i) {
					Template CurrTemplate = (Template) (T1.ForbiddenTemplates
							.get(i));
					if (CurrTemplate.ProteinID == ProteinID2) {
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
			try {
				Line = Buf.readLine();
			} catch (IOException E) {
				E.printStackTrace();
				return null;
			}
		}

		try {
			Buf.close();
		} catch (IOException E) {
			E.printStackTrace();
			return null;
		}

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
*/
	

	public void AddForbiddenEdge(Template B) {

		if (this.ForbiddenTemplates == null) {
			this.ForbiddenTemplates = new ArrayList<Template>();
		}
		if (B.ForbiddenTemplates == null)
			B.ForbiddenTemplates = new ArrayList<Template>();

		boolean Found = false;
		if (B.classID == this.classID && B.SourceDB.equalsIgnoreCase(this.SourceDB))
			Found = true;

		for (int i = 0; i < this.ForbiddenTemplates.size(); ++i) {
			Template CurrT = (Template) (this.ForbiddenTemplates.get(i));
			if (CurrT.ProteinID == B.ProteinID) {
				Found = true;
				break;
			}
		}
		if (!Found)
			this.ForbiddenTemplates.add(B);
		Found = false;
		for (int i = 0; i < B.ForbiddenTemplates.size(); ++i) {
			Template CurrT = (Template) (B.ForbiddenTemplates.get(i));
			if (CurrT.ProteinID == this.ProteinID) {
				Found = true;
				break;
			}
		}
		if (!Found)
			B.ForbiddenTemplates.add(this);

	}

	public void addAnchor(Peptide currAnn) {

		if (this.Sequence == null) {
			return;
		}
		boolean LocalDebug = false;

		String currUnModdedAnn = currAnn.getUnModdedPeptideSeq();
		int Start = this.Sequence.indexOf(currUnModdedAnn);
		if (Start < 0) {
			System.err.println("Peptide " + currUnModdedAnn
					+ " could not be found in template " + this.ProteinID + ":"
					+ this.Sequence);
			Utils.WaitForEnter();
			return;
		}
		int End = Start + currUnModdedAnn.length();
		if (LocalDebug)
			System.out.println("Adding new peptide: " + currUnModdedAnn + "["
					+ Start + "-" + End + "]");

		// If this is the first anchor, create a new Seed
		if (this.anchors == null) {
			this.anchors = new ProteinSeed(currAnn, Start, End, this.ProteinID,
					this);
			this.midPoint = this.anchors;
			return;
		}

		// If this is not the first anchor, search the list of current anchors
		// to see if the new peptide
		// overlaps an existing anchor
		Seed anchorPtr = this.anchors;
		while (anchorPtr != null) {

			// Check if this peptide overlaps the current anchor
			if (anchorPtr.hasOverlap(this.ProteinID, Start, End)) {

				// Add the peptide to the anchor
				anchorPtr.addNewAnnotation(this.ProteinID, Start, End, currAnn);

				// The boundaries of this anchor may have changed, so see if we
				// can merge with
				// an adjacent anchor
				anchorPtr = anchorPtr.mergeWithAdjacent();
				if (LocalDebug) {
					System.out.println("Seed [" + Start + "-" + End + "] on "
							+ this.ProteinID);
					System.out.println("Seed[" + anchorPtr.SeedStart + "-"
							+ anchorPtr.SeedEnd + "] on " + anchorPtr.ClassNum);
					System.out.println("Overlap!!");
					anchorPtr.DebugPrint();
					System.out.println("Template Seq: " + this.Sequence);
					Utils.WaitForEnter();

				}
				break;

				// If we've found an anchor that is greater than our peptide
				// range, then insert a new anchor
				// in front of this anchor.
			} else if (anchorPtr.compareTo(this.ProteinID, Start, End) > 0) {
				if (LocalDebug) {
					System.out.println("Inserting before ["
							+ anchorPtr.SeedStart + "-" + anchorPtr.SeedEnd
							+ "]");
					System.out.println("Template Seq: " + this.Sequence);
				}

				Seed NewSeed = new ProteinSeed(currAnn, Start, End,
						this.ProteinID, this);
				NewSeed.NextSeed = anchorPtr;
				if (anchorPtr.PrevSeed == null) // this is the first guy
				{
					this.anchors = NewSeed;
				} else
					anchorPtr.PrevSeed.NextSeed = NewSeed;
				NewSeed.PrevSeed = anchorPtr.PrevSeed;
				anchorPtr.PrevSeed = NewSeed;
				break;
			} else if (anchorPtr.NextSeed == null) // We've reached the end of
													// the
			// list
			{
				if (LocalDebug) {
					System.out.println("Inserting at end of list!");
					System.out.println("Template Seq: " + this.Sequence);
				}
				Seed NewSeed = new ProteinSeed(currAnn, Start, End,
						this.ProteinID, this);
				NewSeed.PrevSeed = anchorPtr;
				anchorPtr.NextSeed = NewSeed;

				break;
			} else if (LocalDebug) {
				System.out.println("Seed [" + Start + "-" + End + "] on "
						+ this.ProteinID);
				System.out.println("Seed[" + anchorPtr.SeedStart + "-"
						+ anchorPtr.SeedEnd + "] on " + anchorPtr.ClassNum);
				System.out.println("No overlap!!");
				System.out.println("Template Seq: " + this.Sequence);
				Utils.WaitForEnter();
			}
			anchorPtr = anchorPtr.NextSeed;

		}

	}

	/**
	 * Creates a linked list of ProteinSeeds from the Peptide objects
	 * 
	 * @param peptideObjs
	 *            An arraylist of Peptide objects mapping to this template
	 * @param trieDB
	 *            The trie formated database containing all template sequences
	 * @param isShuffledTrie
	 *            If true, then the TrieDB contains target and decoy sequences.
	 */
	public void createAnchors(ArrayList<Peptide> peptideObjs, TrieDB trieDB,
			boolean isShuffledTrie) {

		this.SourceDB = trieDB.GetDBFileName();

		// ASSUMPTION: ShuffledDB is of the form
		// Shuffled1*True1*Shuffled2*True2...
		if (isShuffledTrie) {
			this.Sequence = trieDB.getProteinSequence(this.ProteinID * 2 + 1);
			this.TemplateName = trieDB.getProteinName(this.ProteinID * 2 + 1);
		} else {
			this.Sequence = trieDB.getProteinSequence(this.ProteinID);
			this.TemplateName = trieDB.getProteinName(this.ProteinID);
		}

		// Make sure that we can get basic information about this template from
		// the trieDB
		if (this.TemplateName == null || this.Sequence == null) {
			System.out.println("TEMPLATE NAME OR SEQUENCE IS NULL!");
			System.out.println("ProteinID: " + this.ProteinID);
			System.out.println("Total proteins: " + trieDB.getNumProteins());
			for (int i = 0; i < peptideObjs.size(); ++i) {
				Peptide CurrAnn = (Peptide) (peptideObjs.get(i));
				System.out.println(CurrAnn.toString());
			}
			Utils.WaitForEnter();
		}
		boolean LocalDebug = false;

		// To speed things up, let's sort the Peptides by increasing position on
		// the template
		for (int i = 0; i < peptideObjs.size(); ++i) {

			String currUnModdedAnn = ((Peptide) peptideObjs.get(i))
					.getUnModdedPeptideSeq();

			int minStart = this.Sequence.indexOf(currUnModdedAnn);
			int minIndex = i;

			for (int j = i + 1; j < peptideObjs.size(); ++j) {
				currUnModdedAnn = ((Peptide) peptideObjs.get(j))
						.getUnModdedPeptideSeq();
				int newStart = this.Sequence.indexOf(currUnModdedAnn);
				if (newStart < minStart) {
					minStart = newStart;
					minIndex = j;
				}
			}
			if (minIndex != i) {
				Peptide temp = (Peptide) peptideObjs.remove(i);
				peptideObjs.add(i, peptideObjs.remove(minIndex - 1));
				peptideObjs.add(minIndex, temp);
			}

		}

		// Iterate over the Peptides and create anchors
		for (int i = 0; i < peptideObjs.size(); ++i) {
			Peptide currAnn = (Peptide) (peptideObjs.get(i));
			if (currAnn == null) {
				System.err.println("ERROR: WE ARE TRYING TO ADD A NULL ANN!");
			}
			this.addAnchor(currAnn);

			if (LocalDebug) {
				String CurrUnModdedAnn = currAnn.getUnModdedPeptideSeq();
				System.out.println("Adding " + CurrUnModdedAnn);
				this.DebugPrint();
				Utils.WaitForEnter();
			}
		}
	}

	/**
	 * Determines whether the template with the given ProteinID can be reached
	 * via a valid path from this template. Used for ordering anchors.
	 * 
	 * @param ProteinID
	 *            ID of protein to attempt to reach
	 * @return true, if ProteinID is reachable, and false otherwise
	 */
	public boolean CanReach(int ProteinID) {
		if (this.ProteinID == ProteinID)
			return true;

		for (int i = 0; this.NextTemplates != null
				&& i < this.NextTemplates.size(); ++i) {
			Template Next = (Template) (this.NextTemplates.get(i));
			if (Next.CanReach(ProteinID))
				return true;
		}
		return false;

	}

	/**
	 * @return Returns the number of amino acids covered by anchors
	 */
	public int GetOrigCoverage() {
		int Count = 0;
		Seed Ptr = this.anchors;
		while (Ptr != null) {
			Count += Ptr.AnchoredSeedSequence.length();
			Ptr = Ptr.NextSeed;
		}
		return Count;
	}

	
	/**
	 * Returns the sorted lists of templates based on order constraints. 
	 * @param AllTemplates
	 * @return
	 */
	private static Template[] TopologicalSort(Template[] AllTemplates) {

		ArrayList ret = new ArrayList();
		HashSet retIDs = new HashSet();
		HashSet parents = new HashSet();

		// Add parent nodes to parents
		for (int i = 0; i < AllTemplates.length; ++i) {
			if (AllTemplates[i].PrevTemplates == null
					|| AllTemplates[i].PrevTemplates.size() == 0) {
				parents.add(AllTemplates[i]);

			}
		}

		// System.out.println("Parents contains " + parents.size() +
		// " starting nodes");
		while (parents.size() > 0) {

			Iterator pIter = parents.iterator();
			Template n = (Template) (pIter.next());
			parents.remove(n);

			retIDs.add(new Integer(n.ProteinID));
			ret.add(n);
			// System.out.println("Popped " + n.ProteinID);

			for (int j = 0; n.NextTemplates != null
					&& j < n.NextTemplates.size(); ++j) {
				Template m = (Template) (n.NextTemplates.get(j));
				// System.out.println("Can we add child " + m.ProteinID);
				// Count how many remaining parents this template has
				int remainingParentCount = 0;
				for (int k = 0; k < m.PrevTemplates.size(); ++k) {
					Integer parentID = new Integer(
							((Template) (m.PrevTemplates.get(k))).ProteinID);
					if (retIDs.contains(parentID))
						continue;
					remainingParentCount += 1;
					// System.out.println("No");
				}
				if (remainingParentCount == 0) {
					parents.add(m);
					// System.out.println("Yes!");
				}
			}
			// ForEnter();
		}
		if (ret.size() != AllTemplates.length) {
			System.err.println("ERROR: We only sorted " + ret.size() + "/"
					+ AllTemplates.length);
			return null;
		}
		for (int i = 0; i < AllTemplates.length; ++i)
			AllTemplates[i] = (Template) (ret.get(i));

		return AllTemplates;
	}

	private static Template[] TopologicalSortOld(Template[] AllTemplates) {
		ArrayList NewTemplates = new ArrayList();
		Hashtable ProteinIDHash = new Hashtable();

		boolean LocalDebug = false;

		for (int i = 0; i < AllTemplates.length; ++i) {
			if (AllTemplates[i].PrevTemplates == null
					|| AllTemplates[i].PrevTemplates.size() == 0) {
				NewTemplates.add(AllTemplates[i]);
				ProteinIDHash.put(new Integer(AllTemplates[i].ProteinID),
						new Integer(1));
				if (LocalDebug)
					System.out.println("Adding " + AllTemplates[i].ProteinID
							+ " to no parents list");
			}
		}
		int Ptr = 0;
		while (Ptr < NewTemplates.size()) {
			Template CurrTemplate = (Template) (NewTemplates.get(Ptr));
			if (LocalDebug)
				System.out.println("[" + Ptr + "] - considering "
						+ CurrTemplate.ProteinID);
			ArrayList NextTemplates = CurrTemplate.NextTemplates;
			for (int i = 0; NextTemplates != null && i < NextTemplates.size(); ++i) {
				Template NextTemplate = (Template) (NextTemplates.get(i));
				if (LocalDebug)
					System.out.println("  looking at next template "
							+ NextTemplate.ProteinID);
				if (ProteinIDHash.containsKey(new Integer(
						NextTemplate.ProteinID))) {
					if (LocalDebug)
						System.out.println("  Skipping child of "
								+ CurrTemplate.ProteinID + " = "
								+ NextTemplate.ProteinID);
					continue;
				}
				// Check that there are no parents of NextTemplate not in
				// NewTemplates
				boolean Found = false;
				if (LocalDebug) {
					if (NextTemplate.PrevTemplates == null)
						System.out.println("Template " + NextTemplate.ProteinID
								+ " has no other parents!!");
					else
						System.out.println("Template " + NextTemplate.ProteinID
								+ " has " + NextTemplate.PrevTemplates.size()
								+ " other parents!!");
				}
				for (int j = 0; NextTemplate.PrevTemplates != null
						&& j < NextTemplate.PrevTemplates.size(); ++j) {
					Template PrevTemplate = (Template) (NextTemplate.PrevTemplates
							.get(j));
					if (LocalDebug)
						System.out.println("  Looking at parent of "
								+ NextTemplate.ProteinID + " = "
								+ PrevTemplate.ProteinID);
					if (!ProteinIDHash.containsKey(new Integer(
							PrevTemplate.ProteinID))) {
						Found = true;
						if (LocalDebug)
							System.out.println("  Cannot add child of "
									+ CurrTemplate.ProteinID + " = "
									+ NextTemplate.ProteinID + " since parent "
									+ PrevTemplate.ProteinID
									+ " has not been visited");
						break;
					}
				}
				if (!Found) {
					if (LocalDebug)
						System.out.println("  Adding child of "
								+ CurrTemplate.ProteinID + " = "
								+ NextTemplate.ProteinID);
					ProteinIDHash.put(new Integer(NextTemplate.ProteinID),
							new Integer(1));
					NewTemplates.add(NextTemplate);
					if (LocalDebug)
						System.out.println("Adding " + NextTemplate.ProteinID
								+ " to no parents list");
				}
			}
			Ptr += 1;
			if (LocalDebug)
				Utils.WaitForEnter();
		}

		// Sanity check: we didn't lose any templates
		if (NewTemplates.size() != AllTemplates.length) {
			System.err.println("ERROR: We only sorted " + NewTemplates.size()
					+ "/" + AllTemplates.length);
			return null;
		}
		for (int i = 0; i < AllTemplates.length; ++i)
			AllTemplates[i] = (Template) (NewTemplates.get(i));

		return AllTemplates;
	}

	/**
	 * Each template has a score, which is the raw coverage (int AAs). Find a
	 * path in the partial order graph which maximizes the score The path must
	 * not contains two templates which are connected by a forbidden edge
	 * 
	 * @param AllTemplates
	 *            The set of all templates, with order and forbiddeness
	 *            populated for each.
	 * @return Returns the templates in the path, in order and properly updated
	 *         (NextTemplates.size() == 1)
	 */
	public static Template[] ChooseBestTemplateSet(Template[] AllTemplates) {
		boolean LocalDebug = false;
		Template[] SortedAllTemplates = Template.TopologicalSort(AllTemplates);

		if (SortedAllTemplates == null)
			return null;

		// Try greedy
		// Find the highest scoring
		/*
		 * HashSet classIDs = new HashSet(); ArrayList bestIndexesTemp = new
		 * ArrayList(); for(int i = 0; i < SortedAllTemplates.length; ++i) {
		 * if(classIDs.contains(new Integer(SortedAllTemplates[i].classID)))
		 * continue; classIDs.add(new Integer(SortedAllTemplates[i].classID));
		 * int best =
		 * Template.findBestUnSelectedTemplate(SortedAllTemplates,SortedAllTemplates
		 * [i].classID); //System.out.println("Best score for class " +
		 * SortedAllTemplates[i].classID + ":" + best); if(best >= 0)
		 * bestIndexesTemp.add(new Integer(best)); } if(bestIndexesTemp.size() >
		 * 0) { int[] bestIndexes =
		 * Utils.ConvertArraylistToIntArray(bestIndexesTemp);
		 * 
		 * //Find out if we have a path //If this is sorted then the parentless
		 * template should be first //ArrayList paths = new ArrayList();
		 * ArrayList FinalSet = new ArrayList(); int bestScore = -1; int
		 * startIndex=-1, score=0, endIndex=-1;
		 * 
		 * for(int i = 0; i < bestIndexes.length; ++i) { Template t =
		 * SortedAllTemplates[bestIndexes[i]]; if(t.PrevTemplates == null ||
		 * t.PrevTemplates.size() == 0) { startIndex = i; score =
		 * t.GetOrigCoverage(); } if(t.NextTemplates == null ||
		 * t.NextTemplates.size() == 0) //End of a path { endIndex = i;
		 * if(startIndex >= 0 && score > bestScore) {
		 * //System.out.println("Found a path from " + startIndex + " thru " +
		 * endIndex + " with score " + score); bestScore = score;
		 * FinalSet.clear(); for(int j = startIndex; j <= endIndex; ++j)
		 * FinalSet.add(SortedAllTemplates[bestIndexes[j]]); } } else {
		 * if(!t.CanReach(SortedAllTemplates[bestIndexes[i+1]].ProteinID))
		 * startIndex = -1; else score +=
		 * SortedAllTemplates[bestIndexes[i+1]].GetOrigCoverage(); } }
		 * 
		 * if (FinalSet.size() == 0) return null; Template[] Ret = new
		 * Template[FinalSet.size()];
		 * 
		 * //Template.DeepCopy(FinalSet); for (int i = 0; i < Ret.length; ++i) {
		 * Ret[i] = (Template) (FinalSet.get(i)); Ret[i].SelectedInPath = true;
		 * if (i > 0) { Ret[i].PrevInList = Ret[i-1];
		 * 
		 * Ret[i - 1].NextInList = (Ret[i]); } if (i == 0) Ret[i].PrevInList =
		 * null; if (i == Ret.length - 1) Ret[i].NextInList = null; } return
		 * Ret; }
		 */

		/*
		 * if(LocalDebug) { String S = ""; for(int i = 0; i <
		 * SortedAllTemplates.length; ++i) S += " " +
		 * SortedAllTemplates[i].ProteinID; System.out.println("SortedList: " +
		 * S); Utils.WaitForEnter(); }
		 */
		int[] PathScores = new int[SortedAllTemplates.length];
		int[] BestPrev = new int[SortedAllTemplates.length];
		Hashtable ProteinIDHash = new Hashtable(); // Maps ProteinIDs to
		// position in the Sorted
		// template list

		/*
		 * Initialize scores of paths which are of length 1 (templates with no
		 * preceeders)
		 */
		int BestPathScore = 0;
		int BestPathEnd = -1;
		for (int i = 0; i < PathScores.length; ++i) {
			if (LocalDebug)
				System.out.println("Looking at path to "
						+ SortedAllTemplates[i].ProteinID + ":"
						+ SortedAllTemplates[i].TemplateName);

			BestPrev[i] = -1;
			ProteinIDHash.put(new Integer(SortedAllTemplates[i].ProteinID),
					new Integer(i));
			if (SortedAllTemplates[i].SelectedInPath) {
				if (LocalDebug)
					System.out.println("Already forbidden, continuing");
				continue;
			}
			if (SortedAllTemplates[i].PrevTemplates == null
					|| SortedAllTemplates[i].PrevTemplates.size() == 0) {
				// if(LocalDebug)
				// System.out.println("No prev templates, continuing");
				PathScores[i] = SortedAllTemplates[i].GetOrigCoverage();
				BestPrev[i] = -1;
			} else {
				int TemplateOrigScore = SortedAllTemplates[i].GetOrigCoverage();
				PathScores[i] = TemplateOrigScore;
				BestPrev[i] = -1;
				for (int j = 0; j < SortedAllTemplates[i].PrevTemplates.size(); ++j) {
					Template Prev = (Template) (SortedAllTemplates[i].PrevTemplates
							.get(j));
					if (Prev.SelectedInPath)
						continue;
					int Index = ((Integer) (ProteinIDHash.get(new Integer(
							Prev.ProteinID)))).intValue();
					// if(LocalDebug)
					// System.out.println("  considering path from " +
					// Prev.ProteinID);
					// Consider this better scoring path
					if (PathScores[Index] + TemplateOrigScore > PathScores[i]) {
						int BP = Index;
						boolean Forbidden = false;
						while (BP >= 0) {
							if (SortedAllTemplates[i]
									.HasForbiddenEdge(SortedAllTemplates[BP])) {
								Forbidden = true;
								break;
							}
							BP = BestPrev[BP];
						}
						if (!Forbidden) {
							PathScores[i] = PathScores[Index]
									+ TemplateOrigScore;
							BestPrev[i] = Index;
							if (LocalDebug) {
								System.out.println("New Best Score to "
										+ SortedAllTemplates[i].ProteinID
										+ " from " + Prev.ProteinID + ":"
										+ PathScores[i]);
							}
						} else if (LocalDebug)
							System.out.println("Path to " + Prev.ProteinID
									+ " contained forbidden node "
									+ SortedAllTemplates[BP].ProteinID);
					}
				}
			}
			if (PathScores[i] > BestPathScore) {
				BestPathScore = PathScores[i];
				BestPathEnd = i;
				if (LocalDebug) {
					System.out.println("New best path ending at " + i + ": "
							+ BestPathScore);
					System.out.println("Best Prev: " + BestPrev[i]);
					Utils.WaitForEnter();
				}
			}
			if (LocalDebug) {
				// System.out.println("PathScores: " +
				// AntibodyUtils.IntArrayToString(PathScores));
				// System.out.println("BestPrevs: " +
				// AntibodyUtils.IntArrayToString(BestPrev));
				// Utils.WaitForEnter();
				System.out.println("New path: " + PathScores[i]);
				System.out.println("Prev: " + BestPrev[i]);
				Utils.WaitForEnter();

			}

		}

		ArrayList FinalSet = new ArrayList();
		while (BestPathEnd >= 0) {
			FinalSet.add(0, SortedAllTemplates[BestPathEnd]);
			BestPathEnd = BestPrev[BestPathEnd];
		}

		if (FinalSet.size() == 0)
			return null;
		Template[] Ret = new Template[FinalSet.size()];

		// Template.DeepCopy(FinalSet);
		for (int i = 0; i < Ret.length; ++i) {
			Ret[i] = (Template) (FinalSet.get(i));
			Ret[i].SelectedInPath = true;
			if (i > 0) {
				Ret[i].PrevInList = Ret[i - 1];

				Ret[i - 1].NextInList = (Ret[i]);
			}
			if (i == 0)
				Ret[i].PrevInList = null;
			if (i == Ret.length - 1)
				Ret[i].NextInList = null;
		}

		return Ret;
	}

	private static int findBestUnSelectedTemplate(Template[] allTemplates,
			int classID2) {
		int best = -1;
		int bestScore = -1;
		for (int i = 0; i < allTemplates.length; ++i) {
			if (allTemplates[i].SelectedInPath)
				continue;
			if (allTemplates[i].classID != classID2)
				continue;
			int TemplateOrigScore = allTemplates[i].GetOrigCoverage();
			if (TemplateOrigScore > bestScore) {
				bestScore = TemplateOrigScore;
				best = i;
			}
		}

		return best;
	}

	/*
	 * public String TemplateName = null; public String SourceDB = null;
	 * 
	 * public Seed Anchors = null;
	 * 
	 * public String Sequence = null; public int ProteinID = -1;
	 * 
	 * public int classID = -1;
	 * 
	 * // public Template[] NextTemplates; // public Template[] PrevTemplates;
	 * public ArrayList NextTemplates = null; public ArrayList PrevTemplates =
	 * null;
	 * 
	 * // public Template[] ForbiddenTemplates; public ArrayList
	 * ForbiddenTemplates = null;
	 * 
	 * // Used when searching multiple paths, to mark that a template is used in
	 * a // path, or is forbidden with a template in the path boolean
	 * SelectedInPath = false;
	 * 
	 * // For creating a linked list public Template NextInList = null; public
	 * Template PrevInList = null;
	 */
	private static Template[] DeepCopy(ArrayList finalSet) {
		Template[] ret = new Template[finalSet.size()];

		// First just add new templates
		for (int i = 0; i < ret.length; ++i) {
			Template currTemplate = ((Template) (finalSet.get(i)));
			ret[i] = new Template(currTemplate.ProteinID);
			ret[i].classID = currTemplate.classID;
			ret[i].TemplateName = currTemplate.TemplateName;
			ret[i].Sequence = currTemplate.Sequence;
			ret[i].anchors = currTemplate.anchors;

			if (i > 0) {
				ret[i].PrevTemplates = new ArrayList();
				ret[i].PrevTemplates.add(ret[i - 1]);

				ret[i - 1].NextTemplates = new ArrayList();
				ret[i - 1].NextTemplates.add(ret[i]);
			}
			if (i == 0)
				ret[i].PrevTemplates = null;
			if (i == ret.length - 1)
				ret[i].NextTemplates = null;
		}
		return ret;

	}

	/**
	 * Each template has a score, which is the raw coverage (int AAs). Find all
	 * paths in the partial order graph which maximizes the score The paths must
	 * not contains two templates which are connected by a forbidden edge
	 * 
	 * @param AllTemplates
	 *            The set of all templates, with order and forbiddeness
	 *            populated for each.
	 * @return Returns the templates in the path, in order and properly updated
	 *         (NextTemplates.size() == 1)
	 */
	public static Template[] chooseBestTemplateSetMultiPath(
			Template[] AllTemplates, int MinPeptidesPerPath) {
		boolean LocalDebug = false;

		ArrayList<Template> SelectedTemplates = new ArrayList<Template>();
		if (LocalDebug) {
			System.out.println("Before selection:");
			for (int i = 0; i < AllTemplates.length; ++i)
				AllTemplates[i].DebugPrint();
		}
		while (true) {
			// Find the best path
			if (LocalDebug)
				System.out.println("Finding next best set!");
			Template[] Temp = Template.ChooseBestTemplateSet(AllTemplates);

			if (Temp == null) {
				if (LocalDebug)
					System.out.println("No path could be found!!");
				return Template
						.ConvertArrayListToTemplateArray(SelectedTemplates);
			}
			if (LocalDebug) {
				System.out.println("Set selected: ");
				for (int i = 0; i < Temp.length; ++i)
					Temp[i].DebugPrint();
				Utils.WaitForEnter();
			}
			if (LocalDebug) {
				System.out.println("AllTemplates: ");
				for (int i = 0; i < AllTemplates.length; ++i)
					AllTemplates[i].DebugPrint();
				Utils.WaitForEnter();
			}

			// Check that path has at least XXX Peptides
			int PeptideCount = 0;
			Hashtable<String, Integer> Peptides = new Hashtable<String, Integer>();
			for (int i = 0; i < Temp.length; ++i) {
				for (Seed S = Temp[i].anchors; S != null; S = S.NextSeed) {
					for (int j = 0; j < S.SupportingSpectra.size(); ++j) {

						SpectrumNode CurrAnn = (SpectrumNode) (S.SupportingSpectra
								.get(j));
						String CurrUnModdedAnn = CurrAnn.getUnModdedDBAnnotation();
						String NextPep = CurrUnModdedAnn;
						if (!Peptides.containsKey(NextPep)) {
							PeptideCount += 1;
							Peptides.put(NextPep, new Integer(1));
						}
					}
					if (PeptideCount >= MinPeptidesPerPath) {
						i = Temp.length;
						break;
					}
				}
			}
			if (LocalDebug)
				System.out.println("Path contains more than " + PeptideCount
						+ " peptides");
			if (PeptideCount < MinPeptidesPerPath) {
				if (LocalDebug)
					System.out.println("This is way too few!!");
				return Template
						.ConvertArrayListToTemplateArray(SelectedTemplates);
			}
			// Remove forbidden edge connected templates from the graph
			Template Head = Temp[0];

			for (int i = 0; i < Temp.length; ++i) {

				if (LocalDebug) {
					System.out.println("Removing forbidden edges of ");
					Temp[i].DebugPrint();
				}
				for (int j = 0; Temp[i].ForbiddenTemplates != null
						&& j < Temp[i].ForbiddenTemplates.size(); ++j) {
					((Template) (Temp[i].ForbiddenTemplates.get(j))).SelectedInPath = true;
					if (LocalDebug)
						System.out.println("["
								+ j
								+ "]:"
								+ ((Template) (Temp[i].ForbiddenTemplates
										.get(j))).TemplateName);
				}
				for (int j = 0; j < AllTemplates.length; ++j) {
					if (AllTemplates[j].classID == Temp[i].classID && AllTemplates[j].SourceDB.equalsIgnoreCase(Temp[i].SourceDB))
						AllTemplates[j].SelectedInPath = true;
				}
				if (LocalDebug) {
					System.out.println("For next round:");
					Temp[i].DebugPrint();
					Utils.WaitForEnter();
				}
				if (i < Temp.length - 1) {
					Temp[i].NextInList = Temp[i + 1];
					Temp[i + 1].PrevInList = Temp[i];
				}
			}
			SelectedTemplates.add(Head);
			if (LocalDebug) {
				System.out.println("Done with this path");
				Utils.WaitForEnter();
			}
			// Repeat
		}
	}

	public static Template[] ConvertArrayListToTemplateArray(
			ArrayList<Template> SelectedTemplates) {
		Template[] Ret = new Template[SelectedTemplates.size()];
		for (int i = 0; i < Ret.length; ++i) {
			Ret[i] = (Template) (SelectedTemplates.get(i));
		}
		return Ret;
	}

	/**
	 * Determines if anchors are too similar. Here similarity means that two
	 * anchors: 1. The alignment length is less than the length of both anchors
	 * 2. The alignment length is less than MIN_OVERLAPPING_AA_INTERNAL
	 * 
	 * @param A
	 * @param B
	 * @return
	 */
	public static boolean HasAnchorSimilarity(Template A, Template B,
			double Tolerance) {
		boolean LocalDebug = false;
		if (A.anchors == null || B.anchors == null)
			return false;

		Seed ASeed = A.anchors;
		while (ASeed != null) {
			Seed BSeed = B.anchors;
			while (BSeed != null) {
				if (ASeed.SeedSequence.indexOf(BSeed.SeedSequence) >= 0)
					return true;
				if (BSeed.SeedSequence.indexOf(ASeed.SeedSequence) >= 0)
					return true;
				if (LocalDebug) {
					System.out.println("Comparing anchors of templates "
							+ A.TemplateName + " and " + B.TemplateName);
					System.out.println("A: " + ASeed.SeedSequence);
					System.out.println("B: " + BSeed.SeedSequence);
				}

				int[] ScoreDet = AntibodyUtils.GetLongestAlignmentStrictSeqs(
						ASeed.SeedSequence, BSeed.SeedSequence, Tolerance);
				int Score = ScoreDet[0];
				if (LocalDebug) {
					System.out.println("Alignment Length: " + Score);
					Utils.WaitForEnter();
				}
				if (Score >= ASeed.SeedSequence.length()
						|| Score >= BSeed.SeedSequence.length()) {
					return true;
				}
				BSeed = BSeed.NextSeed;
			}

			ASeed = ASeed.NextSeed;
		}

		return false;
	}

	public boolean HasForbiddenEdge(Template Other) {
		int OtherID = Other.ProteinID;
		int otherClass = Other.classID;

		if (this.classID == otherClass && this.SourceDB.equalsIgnoreCase(Other.SourceDB))
			return true;
		for (int i = 0; this.ForbiddenTemplates != null
				&& i < this.ForbiddenTemplates.size(); ++i) {
			Template Test = (Template) (this.ForbiddenTemplates.get(i));
			if (Test.ProteinID == OtherID)
				return true;
		}
		return false;
	}

	/**
	 * If we are creating a set of templates/constraints from a trie file we
	 * assume that the proteins contain the correct indication of their order in
	 * the name. For example 1.DBFile.ProteinName is mutually exclusive with all
	 * other proteins with 1.DBFile.XXX names and directly precedes all proteins
	 * with 2.DBFile.XXX in the name. We do this in 2 passes, 1 to create the
	 * order constraints, and the second to create the forbidden edges To save
	 * space we only add forbidden edges between constraints at are reachable
	 * 
	 * @param TrieDB
	 * @return Returns a hashtable of hashtables OrigDBName->ClassNum->List of
	 *         Templates
	 */
	private static Hashtable CreateTemplatesFromDBOrderConstraints(TrieDB TrieDB) {

		Hashtable Classes = new Hashtable(); // Classes[OrigDBFile]->Hashtable[Index]->ArrayList)
		boolean LocalDebug = false;

		String DBFile = TrieDB.GetDBFileName();
		int NumProteins = TrieDB.getNumProteins();
		if (LocalDebug)
			System.out.println("Proteins in " + DBFile + ":" + NumProteins);
		for (int ProteinID = 0; ProteinID < NumProteins; ++ProteinID) {
			if (ProteinID % 1000 == 0)
				System.out.println("Built templates for " + ProteinID + "/"
						+ NumProteins);
			String Name = TrieDB.getProteinName(ProteinID);

			String[] Bits = Name.split(AntibodyUtils.CLASS_DELIM_REGEX);
			// int Index = Integer.parseInt(Bits[0]);
			String OrigDBFile = Bits[1];

			int AbsoluteOrder = Integer.parseInt(Bits[0]);

			Template NewTemplate = new Template(ProteinID);
			NewTemplate.classID = AbsoluteOrder;
			NewTemplate.TemplateName = Name;
			NewTemplate.SourceDB = DBFile;

			Hashtable OrigFileHash = null;
			ArrayList ClassList = null;
			if (!Classes.containsKey(OrigDBFile)) {
				if (LocalDebug) {
					System.out
							.println("This is the first time we've seen OrigDBFile: "
									+ OrigDBFile);
					Utils.WaitForEnter();
				}
				OrigFileHash = new Hashtable();
				ClassList = new ArrayList();
			} else {

				OrigFileHash = (Hashtable) (Classes.get(OrigDBFile));
				if (!OrigFileHash.containsKey(new Integer(AbsoluteOrder))) {
					if (LocalDebug) {
						System.out
								.println("This is the first time we've seen class "
										+ AbsoluteOrder
										+ "for file "
										+ OrigDBFile);
						Utils.WaitForEnter();
					}
					ClassList = new ArrayList();
				} else
					ClassList = (ArrayList) (OrigFileHash.get(new Integer(
							AbsoluteOrder)));
			}
			ClassList.add(NewTemplate);
			OrigFileHash.put(new Integer(AbsoluteOrder), ClassList);
			Classes.put(OrigDBFile, OrigFileHash);
		}

		Enumeration OrigDBFiles = Classes.keys();
		while (OrigDBFiles.hasMoreElements()) {
			String SubDB = (String) (OrigDBFiles.nextElement());
			Hashtable OrigFileHash = (Hashtable) (Classes.get(SubDB));
			Object[] ClassNums = (OrigFileHash.keySet().toArray());
			if (LocalDebug)
				System.out.println("Determining order constraints for SubDB: "
						+ SubDB);
			// Sort the class nums
			for (int i = 0; i < ClassNums.length; ++i) {
				int MinValue = ((Integer) (ClassNums[i])).intValue();
				int MinIndex = i;
				for (int j = i + 1; j < ClassNums.length; ++j) {
					if (((Integer) (ClassNums[j])).intValue() < MinValue) {
						MinValue = ((Integer) ClassNums[i]).intValue();
						MinIndex = j;
					}
				}
				if (i != MinIndex) {
					Integer Temp = ((Integer) ClassNums[MinIndex]);
					ClassNums[MinIndex] = ClassNums[i];
					ClassNums[i] = Temp;
				}
			}
			// System.out.println("Loaded " + ClassNums.length + " classes...");
			int Count = 0;
			for (int i = 0; i < ClassNums.length; ++i) {
				ArrayList CurrTemplates = (ArrayList) (OrigFileHash
						.get(ClassNums[i]));
				for (int j = 0; j < CurrTemplates.size(); ++j) {
					Template CurrTemplate = (Template) (CurrTemplates.get(j));
					if (LocalDebug)
						System.out.println("Considering [" + Count + "]:"
								+ CurrTemplate.TemplateName);

					if (i < ClassNums.length - 1) {
						ArrayList NextTemplates = (ArrayList) (OrigFileHash
								.get(ClassNums[i + 1]));
						for (int k = 0; k < NextTemplates.size(); ++k) {

							Template OtherTemplate = (Template) (NextTemplates
									.get(k));
							// if(LocalDebug)
							// System.out.println("  appears before " +
							// OtherTemplate.TemplateName);

							if (CurrTemplate.NextTemplates == null)
								CurrTemplate.NextTemplates = new ArrayList();
							if (OtherTemplate.PrevTemplates == null)
								OtherTemplate.PrevTemplates = new ArrayList();
							CurrTemplate.NextTemplates.add(OtherTemplate);
							OtherTemplate.PrevTemplates.add(CurrTemplate);
						}

					}
					Count += 1;
					if (LocalDebug)
						Utils.WaitForEnter();
				}
			}
		}

		return Classes;
	}

	private static Template[] CreateTemplatesFromDBForbiddenEdges(
			TrieDB TrieDB, Hashtable Classes) {

		boolean LocalDebug = false;

		// TrieDB TrieDB = new TrieDB(DBFile);
		String DBFile = TrieDB.GetDBFileName();
		int NumProteins = TrieDB.getNumProteins();
		if (LocalDebug)
			System.out.println("Proteins in " + DBFile + ":" + NumProteins);

		int Count = 0;
		if (LocalDebug)
			System.out.println("Creating forbidden edges within "
					+ Classes.size() + " classes");
		Enumeration Keys = Classes.keys();
		int LocalCount = 0;
		while (Keys.hasMoreElements()) {
			Hashtable OrigFileHash = (Hashtable) (Classes.get(Keys
					.nextElement()));
			Enumeration ClassKeys = OrigFileHash.keys();
			while (ClassKeys.hasMoreElements()) {
				ArrayList CurrTemplates = (ArrayList) (OrigFileHash
						.get(ClassKeys.nextElement()));
				Count += CurrTemplates.size();
				if (LocalDebug)
					System.out.println("Considering next class:"
							+ CurrTemplates.size());
				for (int i = 0; i < CurrTemplates.size(); ++i) {
					Template CurrTemplate = (Template) (CurrTemplates.get(i));
					// if(LocalDebug)
					// {
					// System.out.println("[" + i + "]:");
					// CurrTemplate.DebugPrint();
					// }
					if (LocalCount % 1000 == 0)
						System.out.println("Built forbidden edges for "
								+ LocalCount + "/" + NumProteins);
					LocalCount += 1;
					// Nothing is reachable from this template, so we won't be
					// adding any forbidden edges
					// if(CurrTemplate.NextTemplates == null)
					// continue;
					for (int j = i + 1; j < CurrTemplates.size(); ++j) {
						Template NextTemplate = (Template) (CurrTemplates
								.get(j));

						// if(CurrTemplate.CanReach(NextTemplate.ProteinID) ||
						// NextTemplate.CanReach(CurrTemplate.ProteinID))
						// {

						// if (CurrTemplate.ForbiddenTemplates == null)
						// CurrTemplate.ForbiddenTemplates = new ArrayList();
						// CurrTemplate.ForbiddenTemplates.add(NextTemplate);
						// if (NextTemplate.ForbiddenTemplates == null)
						// NextTemplate.ForbiddenTemplates = new ArrayList();
						// NextTemplate.ForbiddenTemplates.add(CurrTemplate);
					}
				}
			}
		}
		Enumeration OrigDBFiles = Classes.keys();
		Template[] Ret = new Template[Count];
		Count = 0;
		while (OrigDBFiles.hasMoreElements()) {
			Hashtable OrigFileHash = (Hashtable) (Classes.get(OrigDBFiles
					.nextElement()));
			Object[] ClassNums = (OrigFileHash.keySet().toArray());

			// Sort the class nums
			for (int i = 0; i < ClassNums.length; ++i) {
				int MinValue = ((Integer) (ClassNums[i])).intValue();
				int MinIndex = i;
				for (int j = i + 1; j < ClassNums.length; ++j) {
					if (((Integer) (ClassNums[j])).intValue() < MinValue) {
						MinValue = ((Integer) ClassNums[i]).intValue();
						MinIndex = j;
					}
				}
				if (i != MinIndex) {
					Integer Temp = ((Integer) ClassNums[MinIndex]);
					ClassNums[MinIndex] = ClassNums[i];
					ClassNums[i] = Temp;
				}
			}
			for (int i = 0; i < ClassNums.length; ++i) {
				ArrayList CurrTemplates = (ArrayList) (OrigFileHash
						.get(ClassNums[i]));
				for (int j = 0; j < CurrTemplates.size(); ++j) {
					Template CurrTemplate = (Template) (CurrTemplates.get(j));
					Ret[Count] = CurrTemplate;
					Count += 1;
				}
			}
		}
		// if(LocalDebug)
		// System.out.println("C");
		return Ret;
	}

	public static Template[] CreateTemplatesFromDB(String DBFile) {
		TrieDB TrieFile = new TrieDB(DBFile);
		Hashtable Classes = Template
				.CreateTemplatesFromDBOrderConstraints(TrieFile);
		// Utils.WaitForEnter();
		Template[] Ret = Template.CreateTemplatesFromDBForbiddenEdges(TrieFile,
				Classes);

		// Utils.WaitForEnter();
		return Ret;
	}

	/*
	public static boolean WriteConstraintsToFile(Template[] Templates,
			String ConstraintFile) {

		System.out.println("Writing template constraints to " + ConstraintFile
				+ "...");
		FileWriter Writer = null;
		int ConstraintCount = 0;
		try {
			Writer = new FileWriter(ConstraintFile);
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}

		// Write the constraints. There are three types
		// Information line: i pID cID (protein ID, classID)
		// Order line: o pID1 pID2 (protein ID 1 and protine ID2)
		// Forbidden line: f pID1 pID2 (protein ID1 and protein ID2)
		for (int i = 0; i < Templates.length; ++i) {
			String FileLines = "";
			FileLines += "i\t" + Templates[i].ProteinID + "\t"
					+ Templates[i].classID + "\n";
			for (int n = 0; Templates[i].NextTemplates != null
					&& n < Templates[i].NextTemplates.size(); ++n) {
				Template Next = (Template) (Templates[i].NextTemplates.get(n));
				FileLines += "o\t" + Templates[i].ProteinID + "\t"
						+ Next.ProteinID + "\n";
				ConstraintCount += 1;
			}
			for (int n = 0; Templates[i].ForbiddenTemplates != null
					&& n < Templates[i].ForbiddenTemplates.size(); ++n) {
				Template Next = (Template) (Templates[i].ForbiddenTemplates
						.get(n));
				if (Next.ProteinID > Templates[i].ProteinID) {
					FileLines += "f\t" + Templates[i].ProteinID + "\t"
							+ Next.ProteinID + "\n";
					ConstraintCount += 1;
				}
			}

			try {
				Writer.write(FileLines);
			} catch (IOException E) {
				E.printStackTrace();
				continue;
			}

		}
		System.out.println("Wrote " + ConstraintCount + " constraints to "
				+ ConstraintFile);
		try {
			Writer.close();
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}

		return true;
	}
*/
	

	public static boolean WriteAllSeedSpectraToPKLBIN_OldFormat(
			String outputFileName, Template[] SelectedTemplates) {
		boolean LocalDebug = false;

		if (LocalDebug)
			System.out.println("Template.WriteAllSeedSpectraToPKLBIN: "
					+ outputFileName);

		DataOutputStream f = null;
		int contigCount = 0;
		ArrayList peakCounts = new ArrayList();
		for (int i = 0; i < SelectedTemplates.length; ++i) {
			if (SelectedTemplates[i] == null)
				continue;
			Template currTemplate = SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					contigCount += 1;
					peakCounts
							.add(new Short((short) (currSeed.SeedPRMs.length)));
					currSeed = currSeed.NextSeed;
				}
				currTemplate = currTemplate.NextInList;
			}
		}

		try {
			f = new DataOutputStream(new FileOutputStream(outputFileName));
			// f.writeInt((contigCount));
			f.writeInt(Integer.reverseBytes(contigCount));
			if (LocalDebug)
				System.out.println("Total contigs: " + contigCount);
			for (int i = 0; i < contigCount; ++i) {
				f.writeShort(Short.reverseBytes(((Short) (peakCounts.get(i)))
						.shortValue()));
				if (LocalDebug)
					System.out.println("peakCount[" + i + "]: "
							+ ((Short) (peakCounts.get(i))).shortValue());
			}
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}

		for (int i = 0; i < SelectedTemplates.length; ++i) {
			if (SelectedTemplates[i] == null)
				continue;
			Template currTemplate = SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					try {
						// f.writeFloat((float) 0.0);
						// f.writeDouble(Utils.GetSeqMass(currSeed.SeedSequence));
						// f.writeFloat((float)(Utils.GetSeqMass(currSeed.SeedSequence)));
						f.writeInt(Integer.reverseBytes(Float.floatToIntBits((float) (Utils
								.GetSeqMass(currSeed.SeedSequence) + 19))));
						f.writeInt(Integer.reverseBytes(Float
								.floatToIntBits((float) (1))));
						if (LocalDebug)
							System.out
									.println("PM: "
											+ ((float) (Utils
													.GetSeqMass(currSeed.SeedSequence)))
											+ ", Z=" + ((float) 0.0));

						for (int j = 0; j < currSeed.SeedPRMs.length; ++j) {
							f.writeInt(Integer.reverseBytes(Float
									.floatToIntBits(((float) (currSeed.SeedPRMs[j]) / AntibodyUtils.MASS_SCALE))));
							f.writeInt(Integer.reverseBytes(Float
									.floatToIntBits(((float) (currSeed.SeedPRMScores[j])))));
							// f.writeFloat(((float)(currSeed.SeedPRMs[j])/AntibodyUtils.MASS_SCALE));
							// f.writeFloat(((float)(currSeed.SeedPRMScores[j])));
							if (LocalDebug)
								System.out
										.println((((float) (currSeed.SeedPRMs[j]) / AntibodyUtils.MASS_SCALE))
												+ " "
												+ (((float) (currSeed.SeedPRMScores[j]))));
						}
					} catch (IOException E) {
						E.printStackTrace();
						return false;
					}
					currSeed = currSeed.NextSeed;
				}
				currTemplate = currTemplate.NextInList;
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

	public static boolean writeAllSeedSpectraToPKLBIN(String outputFileName,
			Template[] SelectedTemplates) {
		boolean LocalDebug = false;

		if (LocalDebug)
			System.out.println("Template.WriteAllSeedSpectraToPKLBIN: "
					+ outputFileName);

		DataOutputStream f = null;
		int contigCount = 0;
		ArrayList peakCounts = new ArrayList();
		for (int i = 0; i < SelectedTemplates.length; ++i) {
			if (SelectedTemplates[i] == null)
				continue;
			Template currTemplate = SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					contigCount += 1;
					peakCounts
							.add(new Short((short) (currSeed.SeedPRMs.length)));
					currSeed = currSeed.NextSeed;
				}
				currTemplate = currTemplate.NextInList;
			}
		}

		try {
			f = new DataOutputStream(new FileOutputStream(outputFileName));

			ByteBuffer b = ByteBuffer.allocate(10);
			b.putInt(Integer.reverseBytes(0));
			b.put((byte) 1);
			b.put((byte) 0);
			// b.putChar(Character.reverseBytes((char)1));
			// b.putChar(Character.reverseBytes((char)0));
			b.putInt(Integer.reverseBytes(contigCount));

			// f.writeInt(Integer.reverseBytes(0)); //For the new format we
			// leave the first 4 bytes blank to indicate the new format
			// f.writeChar(Character.reverseBytes((char)1)); //version
			// f.writeChar(Character.reverseBytes((char)0)); //subversion
			// f.writeInt(Integer.reverseBytes(contigCount));
			f.write(b.array());
			if (LocalDebug)
				System.out.println("Total contigs: " + contigCount);

			// Write the scan Numbers
			for (int i = 0; i < contigCount; ++i) {
				f.writeInt(Integer.reverseBytes(i));
			}
			// Write the msLevels
			for (int i = 0; i < contigCount; ++i) {
				f.writeShort(Short.reverseBytes((short) 0));
			}
			// Write the peak counts
			for (int i = 0; i < contigCount; ++i) {
				f.writeShort(Short.reverseBytes(((Short) (peakCounts.get(i)))
						.shortValue()));
				if (LocalDebug)
					System.out.println("peakCount[" + i + "]: "
							+ ((Short) (peakCounts.get(i))).shortValue());
			}
		} catch (IOException E) {
			E.printStackTrace();
			return false;
		}

		for (int i = 0; i < SelectedTemplates.length; ++i) {
			if (SelectedTemplates[i] == null)
				continue;
			Template currTemplate = SelectedTemplates[i];
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					try {
						// f.writeFloat((float) 0.0);
						// f.writeDouble(Utils.GetSeqMass(currSeed.SeedSequence));
						// f.writeFloat((float)(Utils.GetSeqMass(currSeed.SeedSequence)));
						f.writeInt(Integer.reverseBytes(Float.floatToIntBits((float) (Utils
								.GetSeqMass(currSeed.SeedSequence) + 19))));
						f.writeFloat((float) 0.0);
						if (LocalDebug)
							System.out
									.println("PM: "
											+ ((float) (Utils
													.GetSeqMass(currSeed.SeedSequence)))
											+ ", Z=" + ((float) 0.0));

						for (int j = 0; j < currSeed.SeedPRMs.length; ++j) {
							f.writeInt(Integer.reverseBytes(Float
									.floatToIntBits(((float) (currSeed.SeedPRMs[j]) / AntibodyUtils.MASS_SCALE))));
							f.writeInt(Integer.reverseBytes(Float
									.floatToIntBits(((float) (currSeed.SeedPRMScores[j])))));
							// f.writeFloat(((float)(currSeed.SeedPRMs[j])/AntibodyUtils.MASS_SCALE));
							// f.writeFloat(((float)(currSeed.SeedPRMScores[j])));
							if (LocalDebug)
								System.out
										.println((((float) (currSeed.SeedPRMs[j]) / AntibodyUtils.MASS_SCALE))
												+ " "
												+ (((float) (currSeed.SeedPRMScores[j]))));
						}
					} catch (IOException E) {
						E.printStackTrace();
						return false;
					}
					currSeed = currSeed.NextSeed;
				}
				currTemplate = currTemplate.NextInList;
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

	/*
	 * public static boolean WriteAllSeedSpectraToBlankPKLBIN( String
	 * outputFileName, Template[] SelectedTemplates) { boolean LocalDebug =
	 * false;
	 * 
	 * if (LocalDebug)
	 * System.out.println("Template.WriteAllSeedSpectraToBlankPKLBIN: " +
	 * outputFileName);
	 * 
	 * DataOutputStream f = null; int contigCount = 0; // ArrayList peakCounts =
	 * new ArrayList(); for (int i = 0; i < SelectedTemplates.length; ++i) { if
	 * (SelectedTemplates[i] == null) continue; Template currTemplate =
	 * SelectedTemplates[i]; while (currTemplate != null) { Seed currSeed =
	 * currTemplate.anchors; while (currSeed != null) { contigCount += 1; //
	 * peakCounts.add(new Short(0)); currSeed = currSeed.NextSeed; }
	 * currTemplate = currTemplate.NextInList; } }
	 * 
	 * try { f = new DataOutputStream(new FileOutputStream(outputFileName)); //
	 * f.writeInt((contigCount)); f.writeInt(Integer.reverseBytes(contigCount));
	 * if (LocalDebug) System.out.println("Total contigs: " + contigCount); for
	 * (int i = 0; i < contigCount; ++i) {
	 * f.writeShort(Short.reverseBytes(((short) 0))); } } catch (IOException E)
	 * { E.printStackTrace(); return false; }
	 * 
	 * for (int i = 0; i < SelectedTemplates.length; ++i) { if
	 * (SelectedTemplates[i] == null) continue; Template currTemplate =
	 * SelectedTemplates[i]; while (currTemplate != null) { Seed currSeed =
	 * currTemplate.anchors; while (currSeed != null) { try { //
	 * f.writeFloat((float) 0.0); //
	 * f.writeDouble(Utils.GetSeqMass(currSeed.SeedSequence)); //
	 * f.writeFloat((float)(Utils.GetSeqMass(currSeed.SeedSequence)));
	 * f.writeInt(Integer.reverseBytes(Float .floatToIntBits((float) (0))));
	 * f.writeFloat((float) 0.0);
	 * 
	 * } catch (IOException E) { E.printStackTrace(); return false; } currSeed =
	 * currSeed.NextSeed; } currTemplate = currTemplate.NextInList; } }
	 * 
	 * try { f.close(); } catch (IOException E) { E.printStackTrace(); return
	 * false; }
	 * 
	 * return true;
	 * 
	 * }
	 */
	static Template[] CreateConstraintsFromDB(String TemplateConstraintFile,
			String CombinedTrieName) {

		System.out.println("Template Constraint File: "
				+ TemplateConstraintFile);
		System.out.println("CombinedDBName = " + CombinedTrieName);
		if (TemplateConstraintFile == null) {

			return null;
		}
		return Template.CreateTemplatesFromDB(CombinedTrieName);

	}

	public void DebugPrint() {
		System.out.println("------------------Template " + this.ProteinID
				+ "----------------------------");
		if (this.TemplateName != null) {
			System.out.println("Name: " + this.TemplateName);
			System.out.println("SourceDB: " + this.SourceDB);
			System.out.println("Sequence: " + this.Sequence);

		}
		System.out.println("PID: " + this.ProteinID);
		System.out.println("CID: " + this.classID);
		System.out.println("Selected: " + this.SelectedInPath);
		Seed Ptr = this.anchors;
		while (Ptr != null) {
			Ptr.DebugPrint();
			Ptr = Ptr.NextSeed;
		}
		System.out.println("Links: ");

		for (int i = 0; this.NextTemplates != null
				&& i < this.NextTemplates.size(); ++i) {
			Template T = (Template) (this.NextTemplates.get(i));
			System.out.println("  -> " + T.ProteinID);
		}
		for (int i = 0; this.PrevTemplates != null
				&& i < this.PrevTemplates.size(); ++i) {
			Template T = (Template) (this.PrevTemplates.get(i));
			System.out.println("  <- " + T.ProteinID);
		}
		for (int i = 0; this.ForbiddenTemplates != null
				&& i < this.ForbiddenTemplates.size(); ++i) {
			Template T = (Template) (this.ForbiddenTemplates.get(i));
			System.out.println("  X " + T.ProteinID);
		}

	}

	/**
	 * Creates a set of anchors on the template with ID Protein ID
	 * 
	 * @param AllTemplates
	 *            An array of all templates
	 * @param ProteinID
	 *            The ID for the template of interest
	 * @param peptides
	 *            The ArrayList of Peptide objects mapping to this template
	 * @param isInspect
	 *            If the DB search was performed by InsPecT
	 * @param TrieDatabase
	 *            The trie formated DB containing the template sequences
	 */
	public static void createAnchorsOnTemplate(Template[] AllTemplates,
			int ProteinID, ArrayList<Peptide> peptides, boolean isInspect,
			TrieDB TrieDatabase) {

		// Find the template
		Template currTemplate = null;
		for (int i = 0; i < AllTemplates.length; ++i) {
			if (ProteinID == AllTemplates[i].ProteinID) {
				currTemplate = AllTemplates[i];
				break;
			}
		}
		if (currTemplate == null) {
			ErrorThrower.ThrowErrorCustum(ErrorThrower.CUSTOM_ERROR,
					"ERROR: Unable to find template with Protein ID "
							+ ProteinID);

		}

		// If we are running Inspect then we expect the TrieDatabase to contain
		// shuffled proteins, otherwise the TrieDatabase
		// contains only target sequences
		if (isInspect)
			currTemplate.createAnchors(peptides, TrieDatabase, true);
		else
			currTemplate.createAnchors(peptides, TrieDatabase, false);

	}

	public boolean matchesForbiddenEdges(Template template) {
		if (this.ForbiddenTemplates == null) {
			if (template.ForbiddenTemplates == null)
				return true;
			return false;
		}
		if (this.ForbiddenTemplates.size() != template.ForbiddenTemplates
				.size())
			return false;

		int foundCount = 0;
		for (int i = 0; i < this.ForbiddenTemplates.size(); ++i) {
			Template currID = (Template) (this.ForbiddenTemplates.get(i));
			for (int j = 0; j < template.ForbiddenTemplates.size(); ++j) {
				Template otherID = (Template) (template.ForbiddenTemplates
						.get(j));
				if (currID.ProteinID == otherID.ProteinID) {
					foundCount += 1;
					break;
				}

			}
		}
		if (foundCount != this.ForbiddenTemplates.size())
			return false;
		return true;
	}

	public boolean matchesForwardEdges(Template template) {
		if (this.NextTemplates == null) {
			if (template.NextTemplates == null)
				return true;
			return false;
		}
		if (this.NextTemplates.size() != template.NextTemplates.size())
			return false;

		int foundCount = 0;
		for (int i = 0; i < this.NextTemplates.size(); ++i) {
			Template currID = (Template) (this.NextTemplates.get(i));
			for (int j = 0; j < template.NextTemplates.size(); ++j) {
				Template otherID = (Template) (template.NextTemplates.get(j));
				if (currID.ProteinID == otherID.ProteinID) {
					foundCount += 1;
					break;
				}

			}
		}
		if (foundCount != this.NextTemplates.size())
			return false;
		return true;
	}

	public static boolean writeAllSeedSpectraToMGF(String contigSpectrumFile,
			Template[] selectedTemplates) {

		FileWriter f = Utils.openFileWriter(contigSpectrumFile);

		int contigIndex = 1;
		// ArrayList peakCounts = new ArrayList();
		for (int i = 0; i < selectedTemplates.length; ++i) {
			System.out.println("Looking at seq: " + i);
			if (selectedTemplates[i] == null)
				continue;

			Template currTemplate = selectedTemplates[i];
			// currTemplate.DebugPrint();
			while (currTemplate != null) {
				Seed currSeed = currTemplate.anchors;
				while (currSeed != null) {
					// currSeed.DebugPrint();
					Utils.writeLine(f, contigSpectrumFile, "BEGIN IONS\n");
					Utils.writeLine(f,  contigSpectrumFile, "PEPMASS=" + AntibodyUtils.getPepMHFromSeq(currSeed.SeedSequence) + "\n");
					//Utils.writeLine(f,, string);
					//Utils.writeLine(f, contigSpectrumFile, "PEPMASS=0.00000\n");
					Utils.writeLine(f, contigSpectrumFile, "CHARGE=1+\n");
					Utils.writeLine(
							f,
							contigSpectrumFile,
							"FILENAME="
									+ Utils.GetBaseNameNoExtension(contigSpectrumFile)
									+ ".pklbin\n");
					Utils.writeLine(f, contigSpectrumFile,
							"TITLE=Scan Number: " + contigIndex + "\n");
					Utils.writeLine(f, contigSpectrumFile, "SCANS="
							+ contigIndex + "\n");
					Utils.writeLine(f, contigSpectrumFile, "ACTIVATION=CID\n");
					for (int j = 0; j < currSeed.SeedPRMs.length; ++j) {
						Utils.writeLine(f, contigSpectrumFile,
								((double) (currSeed.SeedPRMs[j]))
										/ AntibodyUtils.MASS_SCALE + " "
										+ currSeed.SeedPRMScores[j] + "\n");
					}

					Utils.writeLine(f, contigSpectrumFile, "END IONS\n");

					currSeed = currSeed.NextSeed;
					contigIndex += 1;
				}
				currTemplate = currTemplate.NextInList;
			}
		}
		Utils.closeFileWriter(f, contigSpectrumFile);
		return true;
	}

	/**
	 * Creates an ordered edge between this template and the next template (and a recirocal edge from the next template
	 * back to this template)
	 * @param template
	 */
	public void addForwardEdge(Template template) {
		if(this.NextTemplates == null)
			this.NextTemplates = new ArrayList<Template>();
		if(template.PrevTemplates == null)
			template.PrevTemplates = new ArrayList<Template>();
		
		boolean found = false;
		//make sure the edge doesn't already exists
		for(int i = 0; i < this.NextTemplates.size(); ++i) {
			if(template.ProteinID == this.NextTemplates.get(i).ProteinID)
			{
				found = true;
				break;
			}
		}
		if(!found)
			this.NextTemplates.add(template);
		found = false;
		for(int i = 0; i < template.PrevTemplates.size(); ++i) {
			if(this.ProteinID == template.PrevTemplates.get(i).ProteinID)
			{
				found = true;
				break;
			}
		}
		if(!found)
			template.PrevTemplates.add(this);
		
	}
}
