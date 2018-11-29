import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

public class GetTruePositivesAndTrueNegativesFromsiRNADataRevision_OnlyPValues 
{
	public static void main(String[] args)
	{
		String inputFilePath = args[0];
		String TFName = args[1];
		String truePositiveOutputFilePath = args[2];
		String trueNegativeOutputFilePath = args[3];
		String allGeneFilePath = args[4];
		
		ArrayList<String> geneEnsembles = new ArrayList<String>();
		ArrayList<String> geneNames = new ArrayList<String>();
		ArrayList<String> knockdownTFNames = new ArrayList<String>();
		ArrayList<String> chipseqTFNames = new ArrayList<String>();
		TreeMap<Double, String> trueNegativesPToGeneNames = new TreeMap<Double, String>();
		double[][] pValues = new double[8872][59];
		int[][] peakCounts = new int[8872][201];
		
		TreeMap<String, ArrayList<String>> allGenesToTranscripts = readAllGeneFile(allGeneFilePath);
		readInputFile(inputFilePath, pValues, peakCounts, geneEnsembles, geneNames, knockdownTFNames, chipseqTFNames);
		ArrayList<String> truePositivesForOutput = obtainTruePositivesForTFName(TFName, geneEnsembles, geneNames, knockdownTFNames, chipseqTFNames, pValues, peakCounts, trueNegativesPToGeneNames, allGenesToTranscripts);
		ArrayList<String> trueNegativesForOutput = obtainTrueNegativesForTFName(trueNegativesPToGeneNames, truePositivesForOutput.size());
		output(truePositivesForOutput, truePositiveOutputFilePath);
		output(trueNegativesForOutput, trueNegativeOutputFilePath);
		System.out.println("done");
	}
	
	private static TreeMap<String, ArrayList<String>> readAllGeneFile(String allGeneFilePath)
	{
		TreeMap<String, ArrayList<String>> allGenesToTranscripts = new TreeMap<String, ArrayList<String>>();
		
		try
		{
			Scanner input = new Scanner(new File(allGeneFilePath));
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
		
		return allGenesToTranscripts;
	}
	
	private static void output(ArrayList<String> truePositivesOrNegatives, String truePositiveOrNegativeOutputFilePath)
	{
		try
		{
			PrintWriter outputer = new PrintWriter(truePositiveOrNegativeOutputFilePath);
			for(int i = 0; i < truePositivesOrNegatives.size(); i++)
			{
				outputer.println(truePositivesOrNegatives.get(i));
			}
			outputer.close();
		}
		catch(Exception e)
		{
			System.out.println("Unable to open the file: " + truePositiveOrNegativeOutputFilePath);
			e.printStackTrace();
			System.exit(2);
		}
	}
	
	private static ArrayList<String> obtainTrueNegativesForTFName(TreeMap<Double, String> trueNegativesPToGeneNames, int numberOfTruePositives)
	{
		ArrayList<String> trueNegativesForOutput = new ArrayList<String>();
		Map<Double, String> descendingMap = trueNegativesPToGeneNames.descendingMap();
		Set<Map.Entry<Double, String>> entrySet = descendingMap.entrySet();
		Iterator<Map.Entry<Double, String>> it = entrySet.iterator();
		
		int i = 0;
		while(i < numberOfTruePositives)
		{
			Map.Entry<Double, String> anEntry = it.next();
			trueNegativesForOutput.add(anEntry.getValue() + "," + anEntry.getKey());
			i++;
		}
		
		return trueNegativesForOutput;
	}
	
	private static ArrayList<String> obtainTruePositivesForTFName(String TFName, ArrayList<String> geneEnsembles, ArrayList<String> geneNames, ArrayList<String> knockdownTFNames, ArrayList<String> chipseqTFNames, double[][] pValues, int[][] peakCounts, TreeMap<Double, String> trueNegativesPToGeneNames, TreeMap<String, ArrayList<String>> allGenesToTranscripts)
	{
		ArrayList<String> truePositivesForOutput = new ArrayList<String>();
		TreeMap<Double, String> truePositivesPToGeneNames = new TreeMap<Double, String>();
		int TFIndexInKnockdownTFNames = knockdownTFNames.indexOf(TFName);
		int TFIndexInChipseqTFNames = chipseqTFNames.indexOf(TFName);
		
		for(int i = 0; i < 8872; i++)
		{
			String ensembleNo = geneEnsembles.get(i);
			String aName = geneNames.get(i);
			if(pValues[i][TFIndexInKnockdownTFNames] >= 0.01)
			{
				if(DetermineIfThisGeneIsPC(allGenesToTranscripts, ensembleNo))
				{
					trueNegativesPToGeneNames.put(pValues[i][TFIndexInKnockdownTFNames], ensembleNo + "," + aName);
				}
			}
			else				
			{
				if(pValues[i][TFIndexInKnockdownTFNames] < 0.01 && pValues[i][TFIndexInKnockdownTFNames] >= 0 && peakCounts[i][TFIndexInChipseqTFNames] > 0)
				{
					if(DetermineIfThisGeneIsPC(allGenesToTranscripts, ensembleNo))
					{
						truePositivesPToGeneNames.put(pValues[i][TFIndexInKnockdownTFNames], ensembleNo + "," + aName + "," + peakCounts[i][TFIndexInChipseqTFNames]);
					}
				}
			}
		}
		
		Set<Map.Entry<Double, String>> entrySet = truePositivesPToGeneNames.entrySet();
		Iterator<Map.Entry<Double, String>> it = entrySet.iterator();
		while(it.hasNext())
		{
			Map.Entry<Double, String> anEntry = it.next();
			truePositivesForOutput.add(anEntry.getValue() + "," + anEntry.getKey());
		}
		
		return truePositivesForOutput;
	}
	
	private static boolean DetermineIfThisGeneIsPC(TreeMap<String, ArrayList<String>> allGenesToTranscripts, String ensembleNo)
	{
		ArrayList<String> allTranscriptsOfThisGene = allGenesToTranscripts.get(ensembleNo);
		if(allTranscriptsOfThisGene == null)
		{
			System.out.println("The Ensembl No. " + ensembleNo + " is absent in the file containing all genes.");
			return false;
		}
		else
		{
			String[] fields = allTranscriptsOfThisGene.get(0).split(",");
			if(fields[fields.length - 2].contains("protein_coding"))
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	
	private static void readInputFile(String inputFilePath, double[][] pValues, int[][] peakCounts, ArrayList<String> geneEnsembles, ArrayList<String> geneNames, ArrayList<String> knockdownTFNames, ArrayList<String> chipseqTFNames)
	{
		try
		{
			Scanner inputer = new Scanner(new File(inputFilePath));
			String firstLine = inputer.nextLine();
			String[] fields = firstLine.split("\t");
			for(int i = 3; i <= 61; i++)
			{
				knockdownTFNames.add(fields[i].split("_")[0]);
			}
			for(int i = 62; i <= 262; i++)
			{
				chipseqTFNames.add(fields[i].split("_")[0]);
			}
			
			System.out.println(knockdownTFNames.get(0));
			System.out.println(knockdownTFNames.get(knockdownTFNames.size() - 1));
			System.out.println(chipseqTFNames.get(0));
			System.out.println(chipseqTFNames.get(chipseqTFNames.size() - 1));
			
			int i = 0;
			while(inputer.hasNextLine())
			{
				String aLine = inputer.nextLine();
				fields = aLine.split("\t");
				geneEnsembles.add(fields[1]);
				geneNames.add(fields[2]);
				for(int j = 3; j <= 61; j++)
				{
					if(!fields[j].equals("NA"))
					{
						pValues[i][j-3] = Double.parseDouble(fields[j]);
					}
					else
					{
						pValues[i][j-3] = -1;
					}
				}
				for(int j = 62; j <= 262; j++)
				{
					peakCounts[i][j-62] = Integer.parseInt(fields[j]);
				}
				i++;
			}
			
			inputer.close();
		}
		catch(Exception e)
		{
			System.out.println("An exception in the method: readInputFile");
			e.printStackTrace();
			System.exit(1);
		}
	}
}
