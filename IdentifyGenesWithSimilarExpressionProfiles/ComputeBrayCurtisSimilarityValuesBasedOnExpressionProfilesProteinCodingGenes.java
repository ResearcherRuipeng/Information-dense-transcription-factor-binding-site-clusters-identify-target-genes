import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class ComputeBrayCurtisSimilarityValuesBasedOnExpressionProfilesProteinCodingGenes_Revision
{
	public static void main(String[] args)
	{
		String inputFilePath = args[0];		
		String outputFilePath = args[1];
		String specificGeneEnsemblID = args[2];	
		String allGeneFilePath = args[3];
				
		TreeMap<String, double[]> expressionValuesOfAllGenes = readExpressionValuesOfAllGenes(inputFilePath);
		TreeSet<String> allPCGeneEnsemblIDs = getAllPCGeneEnsemblIDs(allGeneFilePath);
		TreeMap<String, Double> allSimilarityValues = new TreeMap<String, Double>();
		TreeMap<String, double[]> expressionValuesOfAllPCGenes = getExpressionValuesOfAllPCGenes(expressionValuesOfAllGenes, allPCGeneEnsemblIDs);
		
		double[] expressionValuesOfSpecificGene = null;
		Set<Entry<String, double[]>> setView = expressionValuesOfAllPCGenes.entrySet();
		Iterator<Entry<String, double[]>> it = setView.iterator();
		while(it.hasNext())
		{
			Entry<String, double[]> anEntry = it.next();
			String key = anEntry.getKey();
			if(key.contains(specificGeneEnsemblID))
			{
				expressionValuesOfSpecificGene = anEntry.getValue();
				break;
			}
		}			
		
		Iterator<Entry<String, double[]>> it2 = setView.iterator();
		while(it2.hasNext())
		{
			Entry<String, double[]> anEntry = it2.next();
			String key = anEntry.getKey();
			double[] value = anEntry.getValue();
			
			if(!key.contains(specificGeneEnsemblID))
			{
				try
				{
					double sim = computeSimilarity(expressionValuesOfSpecificGene, value);
					allSimilarityValues.put(key, sim);					
				}
				catch(Exception e)
				{
					System.out.println(e.getMessage() + key);
				}
			}
		}
		
		output(outputFilePath, allSimilarityValues);
		System.out.println("Done");
	}
	
	private static TreeMap<String, double[]> getExpressionValuesOfAllPCGenes(TreeMap<String, double[]> expressionValuesOfAllGenes, TreeSet<String> allPCGeneEnsemblIDs)
	{
		TreeMap<String, double[]> expressionValuesOfAllPCGenes = new TreeMap<String, double[]>();
		Set<Entry<String, double[]>> setView = expressionValuesOfAllGenes.entrySet();
		Iterator<Entry<String, double[]>> it = setView.iterator();
		while(it.hasNext())
		{
			Entry<String, double[]> anEntry = it.next();
			if(allPCGeneEnsemblIDs.contains(anEntry.getKey().split(",")[0]))
			{
				expressionValuesOfAllPCGenes.put(anEntry.getKey(), anEntry.getValue());
			}
		}
		
		return expressionValuesOfAllPCGenes;
	}
	
	private static TreeSet<String> getAllPCGeneEnsemblIDs(String allGeneFilePath)
	{
		TreeSet<String> allPCGeneEnsemblIDs = new TreeSet<String>();
		
		try
		{
			Scanner input = new Scanner(new File(allGeneFilePath));
			input.nextLine();		
			while(input.hasNextLine())
			{
				String aLine = input.nextLine();
				String[] fields = aLine.split(",");
				if(fields[12].contains("protein_coding"))
				{
					allPCGeneEnsemblIDs.add(fields[0]);
				}
			}
			input.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open input file.");
			System.exit(1);
		}
		
		return allPCGeneEnsemblIDs;
	}	
	
	private static void output(String outputFilePath, TreeMap<String, Double> allSimilarityValues)
	{
		try
		{
			PrintWriter output = new PrintWriter(outputFilePath);
			
			Set<Entry<String, Double>> setView = allSimilarityValues.entrySet();
			Iterator<Entry<String, Double>> it = setView.iterator();
			while(it.hasNext())
			{
				Entry<String, Double> anEntry = it.next();				
				output.println(anEntry.getKey() + "," + anEntry.getValue());
			}
			output.close();
		}
		catch(Exception e)
		{
			System.out.println("Unable to open the file: " + outputFilePath);
			System.exit(1);
		}
	}
	
	private static double computeSimilarity(double[] expressionValuesOfSpecificGene, double[] expressionValuesOfaGene) throws Exception
	{
		if(expressionValuesOfaGene == null)
		{
			throw new Exception("Cannot find the expression values of the gene ");
		}
		double numerator = 0.0;
		double denominator = 0.0;
		for(int i = 0; i < expressionValuesOfSpecificGene.length; i++)
		{
			numerator += Math.abs(expressionValuesOfSpecificGene[i] - expressionValuesOfaGene[i]);
			denominator += expressionValuesOfSpecificGene[i] + expressionValuesOfaGene[i];
		}
		
		return 1 - ((denominator == 0.0)? 0.0 : (numerator / denominator));
	}
	
	private static TreeMap<String, double[]> readExpressionValuesOfAllGenes(String inputFilePath)
	{
		TreeMap<String, double[]> expressionValuesOfAllGenes = new TreeMap<String, double[]>();
		
		try
		{
			Scanner input = new Scanner(new File(inputFilePath));
			input.nextLine();
			input.nextLine();
			input.nextLine();
			while(input.hasNextLine())
			{
				double[] expressionValuesOfOneGene = new double[53];
				String aLine = input.nextLine();				
				String[] fields = aLine.split("\t");
				
				for(int i = 2; i < fields.length; i++)
				{
					expressionValuesOfOneGene[i - 2] = Double.parseDouble(fields[i]);
				}
				
				expressionValuesOfAllGenes.put(fields[0].split("\\.")[0] + "," + fields[1], expressionValuesOfOneGene);
			}
			input.close();
		}
		catch(IOException e)
		{
			System.out.println("Unable to open the file: " + inputFilePath);
			e.printStackTrace();
		}
		
		return expressionValuesOfAllGenes;
	}
}
