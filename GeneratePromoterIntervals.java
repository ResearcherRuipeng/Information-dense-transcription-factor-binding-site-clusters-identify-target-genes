import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

public class GeneratePromoterIntervals
{
	public static String chromosome = null;
	public static int TSS = 0;
	public static ArrayList<Integer> smallCoor = new ArrayList<Integer>();
	public static ArrayList<Integer> largeCoor = new ArrayList<Integer>();
	public static ArrayList<Integer> smallCoorTemp;
	public static ArrayList<Integer> largeCoorTemp;
	public static int promoterSmallCoor = 0;
	public static int promoterLargeCoor = 0;
	public static String strand = null;
	public static String geneName = null;
	public static String geneID = null;
	public static String geneType = null;
	public static int numberOfGenes = 0;
	public static String[] fields;
	public static PrintWriter output = null;	
	
	public static void main(String[] args)
	{
		String outputFolderPath = args[2];
		
		TreeMap<String, ArrayList<String>> allGenesToTranscripts = new TreeMap<String, ArrayList<String>>();
		
		try
		{
			Scanner input = new Scanner(new File(args[0]));
			input.nextLine();		
			while(input.hasNextLine())
			{
				String aLine = input.nextLine();
				String[] fields = aLine.split(",");
				if(!allGenesToTranscripts.containsKey(fields[0]))
				{
					ArrayList<String> transcripts = new ArrayList<String>();
					transcripts.add(aLine);
					allGenesToTranscripts.put(fields[0], transcripts);
				}
				else
				{
					ArrayList<String> transcripts = allGenesToTranscripts.get(fields[0]);
					transcripts.add(aLine);
				}
			}
			input.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open input file.");
			System.exit(1);
		}
		
		TreeMap<String, String> aSubsetOfGenes = readASubsetOfGenes(args[1]);		
		Set<Map.Entry<String, String>> setView = aSubsetOfGenes.entrySet();
		Iterator<Map.Entry<String, String>> it = setView.iterator();
		while(it.hasNext())
		{		
			smallCoor.clear();
			largeCoor.clear();
			
			Map.Entry<String, String> anEntry = it.next();
			String emsembleNoOfAGene = anEntry.getKey();
			ArrayList<String> allTranscriptsOfAGene = allGenesToTranscripts.get(emsembleNoOfAGene);
			if(allTranscriptsOfAGene == null)
			{
				System.out.println("The Ensembl No. " + emsembleNoOfAGene + " doesn't exist in the file containing all genes.");
			}
			else
			{
				for(int i = 0; i < allTranscriptsOfAGene.size(); i++)
				{
					fields = allTranscriptsOfAGene.get(i).split(",");					
					
					getPromoterOfFirstTranscriptOfOneGene();	
					
				}
				
				processOneGene(outputFolderPath);				
			}
		}		
		
		System.out.println("Done");
	}
	
	private static TreeMap<String, String> readASubsetOfGenes(String geneListFilePath)
	{
		TreeMap<String, String> aSubsetOfGenes = new TreeMap<String, String>();
		
		try
		{
			Scanner inputer = new Scanner(new File(geneListFilePath));
			while(inputer.hasNextLine())
			{
				String aLine = inputer.nextLine();
				String[] fields = aLine.split(",");
				aSubsetOfGenes.put(fields[0], fields[1]);
			}
			inputer.close();
		}
		catch(Exception e)
		{
			System.out.println("Unable to open input file:" + geneListFilePath);
			e.printStackTrace();
			System.exit(1);
		}
		
		return aSubsetOfGenes;
	}
	
	private static void getPromoterOfFirstTranscriptOfOneGene()
	{
		if(fields[6].equals("1"))
		{					
			TSS = Integer.parseInt(fields[3]);
			promoterLargeCoor = TSS - 1;
			promoterSmallCoor = promoterLargeCoor - 9999;
			if(promoterSmallCoor <= 0)
			{
				promoterSmallCoor = 1;
			}
			smallCoor.add(promoterSmallCoor);
			largeCoor.add(promoterLargeCoor);

			chromosome = fields[2];
			strand = fields[6];
			geneName = fields[4];
			geneID = fields[0];
			geneType = fields[12];
		}
		else
		{					
			TSS = Integer.parseInt(fields[3]);
			promoterSmallCoor = TSS + 1;
			promoterLargeCoor = promoterSmallCoor + 9999;
			smallCoor.add(promoterSmallCoor);
			largeCoor.add(promoterLargeCoor);

			chromosome = fields[2];
			strand = fields[6];
			geneName = fields[4];
			geneID = fields[0];
			geneType = fields[12];
		}
	}
	
	private static void getPromoterOfAnotherTranscriptOfOneGene()
	{
		if(fields[6].equals("1"))
		{	
			TSS = Integer.parseInt(fields[3]);
			promoterLargeCoor = TSS - 1;
			promoterSmallCoor = promoterLargeCoor - 9999;
			if(promoterSmallCoor <= 0)
			{
				promoterSmallCoor = 1;
			}
			smallCoor.add(promoterSmallCoor);
			largeCoor.add(promoterLargeCoor);					
		}
		else
		{
			TSS = Integer.parseInt(fields[3]);
			promoterSmallCoor = TSS + 1;
			promoterLargeCoor = promoterSmallCoor + 9999;
			smallCoor.add(promoterSmallCoor);
			largeCoor.add(promoterLargeCoor);
		}	
	}
	
	private static void getPromoterOfFirstTranscriptOfNextGene()
	{
		if(fields[6].equals("1"))
		{		
			smallCoor.clear();
			largeCoor.clear();

			TSS = Integer.parseInt(fields[3]);
			promoterLargeCoor = TSS - 1;
			promoterSmallCoor = promoterLargeCoor - 9999;
			if(promoterSmallCoor <= 0)
			{
				promoterSmallCoor = 1;
			}
			smallCoor.add(promoterSmallCoor);
			largeCoor.add(promoterLargeCoor);

			chromosome = fields[2];
			strand = fields[6];
			geneName = fields[4];
			geneID = fields[0];
			geneType = fields[12];						
		}
		else
		{			
			smallCoor.clear();
			largeCoor.clear();

			TSS = Integer.parseInt(fields[3]);
			promoterSmallCoor = TSS + 1;
			promoterLargeCoor = promoterSmallCoor + 9999;
			smallCoor.add(promoterSmallCoor);
			largeCoor.add(promoterLargeCoor);

			chromosome = fields[2];
			strand = fields[6];
			geneName = fields[4];
			geneID = fields[0];
			geneType = fields[12];
		}
	}
	
	private static void processOneGene(String outputFolderPath)
	{
		if(strand.equals("1"))
		{
			boolean done = false;						
			int j;
			int k;

			while(!done)
			{						
				done = true;
				smallCoorTemp = new ArrayList<Integer>(smallCoor);
				largeCoorTemp = new ArrayList<Integer>(largeCoor);

				for(j = 0; j < smallCoorTemp.size(); j++)
				{
					for(k = j + 1; k < smallCoorTemp.size(); k++)
					{
						if(smallCoorTemp.get(k) >= smallCoorTemp.get(j) && smallCoorTemp.get(k) <= largeCoorTemp.get(j))
						{
							if(largeCoorTemp.get(k) <= largeCoorTemp.get(j))
							{
								smallCoor.remove(k);
								largeCoor.remove(k);
								done = false;
								break;
							}
							else
							{
								largeCoor.set(j, largeCoor.get(k));
								smallCoor.remove(k);
								largeCoor.remove(k);
								done = false;
								break;
							}
						}

						if(largeCoorTemp.get(k) >= smallCoorTemp.get(j) && largeCoorTemp.get(k) <= largeCoorTemp.get(j))
						{
							if(smallCoorTemp.get(k) >= smallCoorTemp.get(j))
							{
								smallCoor.remove(k);
								largeCoor.remove(k);
								done = false;
								break;
							}
							else
							{
								smallCoor.set(j, smallCoor.get(k));
								smallCoor.remove(k);
								largeCoor.remove(k);
								done = false;
								break;
							}
						}
					}

					if(!done)
					{
						break;
					}								
				}									
			}						
		}
		else
		{
			boolean done = false;						
			int j;
			int k;

			while(!done)
			{						
				done = true;
				smallCoorTemp = new ArrayList<Integer>(smallCoor);
				largeCoorTemp = new ArrayList<Integer>(largeCoor);

				for(j = 0; j < smallCoorTemp.size(); j++)
				{
					for(k = j + 1; k < smallCoorTemp.size(); k++)
					{
						if(largeCoorTemp.get(k) >= smallCoorTemp.get(j) && largeCoorTemp.get(k) <= largeCoorTemp.get(j))
						{
							if(smallCoorTemp.get(k) >= smallCoorTemp.get(j))
							{
								smallCoor.remove(k);
								largeCoor.remove(k);
								done = false;
								break;
							}
							else
							{
								smallCoor.set(j, smallCoor.get(k));
								smallCoor.remove(k);
								largeCoor.remove(k);
								done = false;
								break;
							}
						}

						if(smallCoorTemp.get(k) >= smallCoorTemp.get(j) && smallCoorTemp.get(k) <= largeCoorTemp.get(j))
						{
							if(largeCoorTemp.get(k) <= largeCoorTemp.get(j))
							{
								smallCoor.remove(k);
								largeCoor.remove(k);
								done = false;
								break;
							}
							else
							{
								largeCoor.set(j, largeCoor.get(k));
								smallCoor.remove(k);
								largeCoor.remove(k);
								done = false;
								break;
							}
						}
					}

					if(!done)
					{
						break;
					}								
				}									
			}
		}
		
		for(int m = 0; m < smallCoor.size(); m++)
		{
			String currentGene = "chr" + chromosome + "\t" + (smallCoor.get(m) - 1) + "\t" + largeCoor.get(m) + "\t" + geneName + "\t" + geneID + "\t" + geneType + "\t" + strand + "\n";
			try
			{
				if(smallCoor.size() == 1)
				{
					output = new PrintWriter(outputFolderPath + "\\" + geneName + ".txt");
				}
				else
				{
					output = new PrintWriter(outputFolderPath + "\\" + geneName + "_" + m + ".txt");
				}							
				output.print(currentGene);
				output.close();							
			}
			catch(FileNotFoundException e)
			{
				System.out.println("Unable to open output file.");
				System.exit(1);
			}
		}					
		System.out.println("Number of genes processed: " + (++numberOfGenes));
	}
}