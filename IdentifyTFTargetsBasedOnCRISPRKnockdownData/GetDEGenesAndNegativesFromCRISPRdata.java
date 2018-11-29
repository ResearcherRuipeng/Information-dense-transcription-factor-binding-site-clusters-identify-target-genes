import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

public class GetPCDEGenesAndTrueNegativesFromCRISPRdata_Revision
{
	public static void main(String[] args)
	{
		String inputFilePath = args[0];
		String TFName = args[1];
		String allGeneFilePath = args[2];		
		String outputFilePathTruePos = args[3];
		String outputFilePathTrueNeg = args[4];
		Double thresholdOnFoldChange = Double.parseDouble(args[5]);
		
		String[] sgRNAs = new String[26];
		String[] genes = new String[22026];
		TreeMap<String, Double> geneToCoeff = new TreeMap<String, Double>();
		
		double[][] values = readInputFile(sgRNAs, genes, inputFilePath);
		ArrayList<Integer> TFIndices = getTFIndices(sgRNAs, TFName);
		ArrayList<String> TruePosGeneNames = getTruePosWithChangesOfTheSameDirectionUnderAllGuideRNAs(values, TFIndices, genes, geneToCoeff, thresholdOnFoldChange);
		ArrayList<String> TrueNegGeneNames = getTrueNeg(values, TFIndices, genes);
		
		TreeMap<String, ArrayList<String>> allGenesToTranscripts = readAllGeneFile(allGeneFilePath);
		TreeMap<String, String> ensemblNos = getEnsemblNoForTruePos(allGenesToTranscripts, TruePosGeneNames);
		
		outputTP(ensemblNos, outputFilePathTruePos, geneToCoeff);		
		outputTN(TrueNegGeneNames, outputFilePathTrueNeg, allGenesToTranscripts);
		
		System.out.println("Done.");
	}	
	
	private static void outputTN(ArrayList<String> TrueNegGeneNames, String outputFilePathTrueNeg, TreeMap<String, ArrayList<String>> allGenesToTranscripts)
	{
		try
		{
			PrintWriter pw = new PrintWriter(outputFilePathTrueNeg);
			for(int i = 0; i < TrueNegGeneNames.size(); i++)
			{
				pw.println(TrueNegGeneNames.get(i));
			}
			pw.close();
		}
		catch(IOException e)
		{
			System.out.println("Cannot open the file: " + outputFilePathTrueNeg);
			System.exit(1);
		}
	}
	
	private static TreeMap<String, String> getEnsemblNoForTruePos(TreeMap<String, ArrayList<String>> allGenesToTranscripts, ArrayList<String> TruePosGeneNames)
	{
		TreeMap<String, String> ensemblNos = new TreeMap<String, String>();
		
		for(int i = 0; i < TruePosGeneNames.size(); i++)
		{
			String geneName = TruePosGeneNames.get(i);
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
	
	private static void outputTP(TreeMap<String, String> ensemblNos, String outputFilePathTruePos, TreeMap<String, Double> geneToCoeff)
	{
		try
		{
			PrintWriter outputer = new PrintWriter(outputFilePathTruePos);
			Set<Entry<String, String>> setView = ensemblNos.entrySet();
			Iterator<Entry<String, String>> it = setView.iterator();			
			while(it.hasNext())
			{
				Entry<String, String> anEntry = it.next();
				String geneName = anEntry.getValue();
				outputer.println(anEntry.getKey() + "," + geneName + "," + geneToCoeff.get(geneName));				
			}
			outputer.close();
		}
		catch(IOException e)
		{
			System.out.println("Cannot open the output file of true pos: " + outputFilePathTruePos);
			System.exit(1);
		}
	}
	
	private static ArrayList<String> getTrueNeg(double[][] values, ArrayList<Integer> TFIndices, String[] genes)
	{
		ArrayList<String> TrueNegGeneNames = new ArrayList<String>();
		
		for(int i = 0; i < values.length; i++)
		{
			boolean ifThisIsATrueNeg = true;
			for(int j = 0; j < TFIndices.size(); j++)
			{
				if(values[i][TFIndices.get(j)] != 0.0)
				{
					ifThisIsATrueNeg = false;
					break;
				}				
			}
			
			if(!ifThisIsATrueNeg)
			{
				continue;
			}
			else
			{
				TrueNegGeneNames.add(genes[i]);
			}
		}
		
		return TrueNegGeneNames;
	}
	
	private static ArrayList<String> getTruePosWithChangesOfTheSameDirectionUnderAllGuideRNAs(double[][] values, ArrayList<Integer> TFIndices, String[] genes, TreeMap<String, Double> geneToCoeff, Double thresholdOnFoldChange/*, String coeffToGeneOutputFilePathTruePos*/)
	{
		ArrayList<String> TruePosGeneNames = new ArrayList<String>();
		
		for(int i = 0; i < values.length; i++)
		{
			boolean ifThisIsATruePos = true;
			
			if(values[i][TFIndices.get(0)] == 0.0)
			{
				continue;
			}
			
			if(values[i][TFIndices.get(0)] > 0.0)
			{
				for(int j = 1; j < TFIndices.size(); j++)
				{
					if(values[i][TFIndices.get(j)] <= 0.0)
					{
						ifThisIsATruePos = false;
						break;
					}
				}
				
				if(!ifThisIsATruePos)
				{
					continue;
				}
				else
				{
					//set a threshold on the coefficients
					double sum = 0.0;
					for(int j = 0; j < TFIndices.size(); j++)
					{
						sum += values[i][TFIndices.get(j)];						
					}
					
					if(Math.pow(10, sum / TFIndices.size()) >= thresholdOnFoldChange/*0.041393*/)
					{						
						TruePosGeneNames.add(genes[i]);						
						geneToCoeff.put(genes[i], sum);
					}									
				}
			}
			
			if(values[i][TFIndices.get(0)] < 0.0)
			{
				for(int j = 1; j < TFIndices.size(); j++)
				{
					if(values[i][TFIndices.get(j)] >= 0.0)
					{
						ifThisIsATruePos = false;
						break;
					}
				}
				
				if(!ifThisIsATruePos)
				{
					continue;
				}
				else
				{
					//set a threshold on the coefficients
					double sum = 0.0;
					for(int j = 0; j < TFIndices.size(); j++)
					{
						sum += values[i][TFIndices.get(j)];						
					}
					
					if(Math.pow(10, sum / TFIndices.size()) <= (1 / thresholdOnFoldChange))
					{
						TruePosGeneNames.add(genes[i]);						
						geneToCoeff.put(genes[i], Math.abs(sum));
					}
				}
			}
		}		
		
		return TruePosGeneNames;
	}
	
	private static ArrayList<String> getTruePosWithChangesUnderAllGuideRNAs(double[][] values, ArrayList<Integer> TFIndices, String[] genes)
	{
		ArrayList<String> TruePosGeneNames = new ArrayList<String>();
		
		for(int i = 0; i < values.length; i++)
		{
			boolean ifThisIsATruePos = true;
			for(int j = 0; j < TFIndices.size(); j++)
			{
				if(values[i][TFIndices.get(j)] == 0.0)
				{
					ifThisIsATruePos = false;
					break;
				}
			}
			
			if(!ifThisIsATruePos)
			{
				continue;
			}
			else
			{
				TruePosGeneNames.add(genes[i]);
			}
		}
		
		return TruePosGeneNames;
	}
	
	private static ArrayList<Integer> getTFIndices(String[] sgRNAs, String TFName)
	{
		ArrayList<Integer> TFIndices = new ArrayList<Integer>();
		
		for(int i = 0; i < sgRNAs.length; i++)
		{
			if(sgRNAs[i].contains(TFName))
			{
				TFIndices.add(i);
			}
		}
		
		return TFIndices;
	}
	
	private static double[][] readInputFile(String[] sgRNAs, String[] genes, String inputFilePath)
	{
		double[][] values = new double[22026][26];
		
		try
		{
			Scanner inputer = new Scanner(new File(inputFilePath));
			
			String aLine = inputer.nextLine();
			String[] fields = aLine.split(",");
			for(int i = 1; i < fields.length; i++)
			{
				sgRNAs[i - 1] = fields[i];
			}
			
			int lineNumber = 0;
			while(inputer.hasNextLine())
			{
				aLine = inputer.nextLine();
				fields = aLine.split(",");
				genes[lineNumber] = fields[0];
				for(int i = 1; i < fields.length; i++)
				{
					values[lineNumber][i - 1] = Double.parseDouble(fields[i]);
				}
				lineNumber++;
			}
			inputer.close();
		}
		catch(IOException e)
		{
			System.out.println("Cannot open the file: " + inputFilePath);
			System.exit(1);
		}
		
		return values;
	}
}
