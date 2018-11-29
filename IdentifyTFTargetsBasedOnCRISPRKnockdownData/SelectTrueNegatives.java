import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.Map.Entry;

public class SelectTrueNegatives_Revision 
{
	public static void main(String[] args)
	{
		String geneExprValueFilePath = args[0];
		String truePosNameWithLargestDiffExpr = args[1];
		String trueNegFilePath = args[2];
		String outputFilePath = args[3];
		String allGeneFilePath = args[4];
		int numberOfTPs = Integer.parseInt(args[5]);
		
		TreeMap<String, double[]> expressionValuesOfAllGenes = readExpressionValuesOfAllGenes(geneExprValueFilePath);
		double[] exprValuesOfTruePosWithLargestDiffExpr = expressionValuesOfAllGenes.get(truePosNameWithLargestDiffExpr);
		ArrayList<String> trueNeg = readTrueNeg(trueNegFilePath);
		TreeMap<String, ArrayList<String>> allGenesToTranscripts = readAllGeneFile(allGeneFilePath);
		TreeMap<String, String> ensemblNos = getEnsemblNoForTrueNeg(allGenesToTranscripts, trueNeg);
		TreeMap<Double, String> simToGene = computeSim(exprValuesOfTruePosWithLargestDiffExpr, ensemblNos, expressionValuesOfAllGenes);
				
		output(outputFilePath, simToGene, numberOfTPs);
	}
	
	private static TreeMap<String, String> getEnsemblNoForTrueNeg(TreeMap<String, ArrayList<String>> allGenesToTranscripts, ArrayList<String> TrueNegGeneNames)
	{
		TreeMap<String, String> ensemblNos = new TreeMap<String, String>();
		
		for(int i = 0; i < TrueNegGeneNames.size(); i++)
		{
			String geneName = TrueNegGeneNames.get(i);
			ArrayList<String> allLinesOfGene = allGenesToTranscripts.get(geneName);
			if(allLinesOfGene == null)
			{
				System.out.println("The gene name " + geneName + " does not appear in the file containing all genes.");
				continue;
			}
			else
			{
				boolean ifMultipleEnsemblNos = false;
				String[] fields = allLinesOfGene.get(0).split(",");
				String ensemblNo = fields[0];
				for(int j = 1; j < allLinesOfGene.size(); j++)
				{
					String aTranscriptLine = allLinesOfGene.get(j);
					if(!ensemblNo.equals(aTranscriptLine.split(",")[0]))
					{
						System.out.println("The gene name " + geneName + " has multiple Ensembl numbers in the file containing all genes.");
						ifMultipleEnsemblNos = true;
						break;
					}
				}
				
				if(ifMultipleEnsemblNos)
				{
					continue;
				}
				else
				{
					if(fields[fields.length - 2].contains("protein_coding"))
					{
						ensemblNos.put(ensemblNo, geneName);
					}
				}
			}
		}
		
		return ensemblNos;
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
				if(!allGenesToTranscripts.containsKey(fields[4]))
				{
					ArrayList<String> transcripts = new ArrayList<String>();
					transcripts.add(aLine);
					allGenesToTranscripts.put(fields[4], transcripts);
				}
				else
				{
					ArrayList<String> transcripts = allGenesToTranscripts.get(fields[4]);
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
	
	private static void output(String outputFilePath, TreeMap<Double, String> simToGene, int numberOfTPs)
	{
		Set<Entry<Double, String>> setView = simToGene.entrySet();
		Iterator<Entry<Double, String>> it = setView.iterator();
		try
		{
			int counter = 0;
			PrintWriter pw = new PrintWriter(outputFilePath);
			Entry<Double, String> anEntry = it.next();
			if(anEntry.getKey() != 0.0)
			{
				pw.println(anEntry.getValue() + "," + anEntry.getKey());
				counter++;
			}
			while(counter < numberOfTPs)
			{
				anEntry = it.next();				
				pw.println(anEntry.getValue() + "," + anEntry.getKey());
				counter++;
			}
			pw.close();
		}
		catch(IOException e)
		{
			System.out.println("Cannot open the file: " + outputFilePath);
			System.exit(1);
		}
	}
	
	private static TreeMap<Double, String> computeSim(double[] exprValuesOfTruePosWithLargestDiffExpr, TreeMap<String, String> ensemblNos, TreeMap<String, double[]> expressionValuesOfAllGenes)
	{
		TreeMap<Double, String> simToGene = new TreeMap<Double, String>();
		Set<Entry<String, String>> setView = ensemblNos.entrySet();
		Iterator<Entry<String, String>> it = setView.iterator();
		while(it.hasNext())
		{
			Entry<String, String> anEntry = it.next();
			String geneName = anEntry.getValue();
			double[] exprValues = expressionValuesOfAllGenes.get(geneName);
			if(exprValues == null)
			{
				System.out.println("The expr values of the gene " + geneName + " doesn't exist in the file.");
				continue;
			}
			
			double sim = computeSimilarity(exprValuesOfTruePosWithLargestDiffExpr, exprValues);
			String old = simToGene.put(sim, anEntry.getKey() + "," + geneName);
			if(old != null)
			{
				System.out.println("The gene " + geneName + " has the same sim value: " + sim);
			}
		}
		
		return simToGene;
	}
	
	private static double computeSimilarity(double[] expressionValuesOfSpecificGene, double[] expressionValuesOfaGene)
	{		
		double numerator = 0.0;
		double denominator = 0.0;
		for(int i = 0; i < expressionValuesOfSpecificGene.length; i++)
		{
			numerator += Math.abs(expressionValuesOfSpecificGene[i] - expressionValuesOfaGene[i]);
			denominator += expressionValuesOfSpecificGene[i] + expressionValuesOfaGene[i];
		}
		
		return 1 - ((denominator == 0.0)? 0.0 : (numerator / denominator));
	}
	
	private static ArrayList<String> readTrueNeg(String trueNegFilePath)
	{
		ArrayList<String> trueNeg = new ArrayList<String>();
		
		try
		{
			Scanner sc = new Scanner(new File(trueNegFilePath));
			while(sc.hasNextLine())
			{
				trueNeg.add(sc.nextLine());
			}
			sc.close();
		}
		catch(IOException e)
		{
			System.out.println("Cannot open the file: " + trueNegFilePath);
			System.exit(1);
		}
		
		return trueNeg;
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
				
				expressionValuesOfAllGenes.put(fields[1], expressionValuesOfOneGene);
			}
			input.close();
		}
		catch(Exception e)
		{
			System.out.println("Unable to open the file: " + inputFilePath);			
			e.printStackTrace();
			System.exit(1);
		}
		
		return expressionValuesOfAllGenes;
	}
}
