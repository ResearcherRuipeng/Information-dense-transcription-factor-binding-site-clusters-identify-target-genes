import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

public class Main
{
	public static void main(String[] args)
	{
		String featureFolderPath = args[0];
		String truePositiveSimFilePath = args[1];
		String trueNegativeSimFilePath = args[2];
		String outputFilePath = args[3];
		Boolean ifOneiPWM = Boolean.parseBoolean(args[4]);
		
		if(!ifOneiPWM)
		{
			ArrayList<Gene> featureSetsOfAllGenes = getFeatureSetsOfAllGenes(featureFolderPath);
			TreeMap<String, Gene> featureSetsOfLongestTranscripts = getFeatureSetsOfLongestTranscripts(featureSetsOfAllGenes);
			System.out.println("Number of genes with non-empty feature sets: " + featureSetsOfLongestTranscripts.size());
			int maxNumberOfClusters = obtainMaxNumberOfClusters(featureSetsOfLongestTranscripts);
			padZeroAtTheEnd(featureSetsOfLongestTranscripts, maxNumberOfClusters);
			ArrayList<String> geneNamesOfTruePositives = readSimFile(truePositiveSimFilePath);
			ArrayList<String> geneNamesOfTrueNegatives = readSimFile(trueNegativeSimFilePath);
			outputFeatureVectorsForTPTN(featureSetsOfLongestTranscripts, outputFilePath, maxNumberOfClusters, geneNamesOfTruePositives, geneNamesOfTrueNegatives);			
		}
		else
		{
			ArrayList<Gene_oneiPWM> featureSetsOfAllGenes = getFeatureSetsOfAllGenes_oneiPWM(featureFolderPath);
			TreeMap<String, Gene_oneiPWM> featureSetsOfLongestTranscripts = getFeatureSetsOfLongestTranscripts_oneiPWM(featureSetsOfAllGenes);
			System.out.println("Number of genes with non-empty feature sets: " + featureSetsOfLongestTranscripts.size());
			int maxNumberOfClusters = obtainMaxNumberOfClusters_oneiPWM(featureSetsOfLongestTranscripts);
			padZeroAtTheEnd_oneiPWM(featureSetsOfLongestTranscripts, maxNumberOfClusters);
			ArrayList<String> geneNamesOfTruePositives = readSimFile(truePositiveSimFilePath);
			ArrayList<String> geneNamesOfTrueNegatives = readSimFile(trueNegativeSimFilePath);
			outputFeatureVectorsForTPTN_oneiPWM(featureSetsOfLongestTranscripts, outputFilePath, maxNumberOfClusters, geneNamesOfTruePositives, geneNamesOfTrueNegatives);
		}
		System.out.println("done");
	}
	
	private static void padZeroForFeatures36(ArrayList<Double> featureVector, int numberOf0sToBePadded)
	{
		for(int i = 0; i < numberOf0sToBePadded; i++)
		{
			featureVector.add(0.0);
		}
	}
	
	private static void padZeroForFeatures1245(ArrayList<Integer> featureVector, int numberOf0sToBePadded)
	{
		for(int i = 0; i < numberOf0sToBePadded; i++)
		{
			featureVector.add(0);
		}
	}	
	
	private static void outputFeatureVectorsForTPTN_oneiPWM(TreeMap<String, Gene_oneiPWM> featureSetsOfLongestTranscripts, String outputFilePath, int maxNumberOfClusters, ArrayList<String> geneNamesOfTruePositives, ArrayList<String> geneNamesOfTrueNegatives)
	{
		try
		{
			PrintWriter outputer = new PrintWriter(outputFilePath);
			int numberOfColumns = 1 + maxNumberOfClusters * 6 + 1;
			for(int i = 0; i < numberOfColumns - 1; i++)
			{
				outputer.print((i + 1) + ",");
			}
			outputer.print("Label\n");
			
			for(int i = 0; i < geneNamesOfTruePositives.size(); i++)
			{
				String aGeneName = geneNamesOfTruePositives.get(i);
				Gene_oneiPWM anInstance = featureSetsOfLongestTranscripts.get(aGeneName);
				if(anInstance == null)
				{
					outputer.print(aGeneName + ",");
					for(int j = 1; j < numberOfColumns - 1; j++)
					{
						outputer.print(0 + ",");
					}
					outputer.print("1\n");
				}
				else
				{
					outputer.print(anInstance.getName() + ",");
					FeatureSetOfATSS_oneiPWM itsFeatureSet = anInstance.getFeatureSets().get(0);

					ArrayList<Integer> distances = itsFeatureSet.getFeature1().getDistances();
					for(int j = 0; j < distances.size(); j++)
					{
						outputer.print(distances.get(j) + ",");
					}

					ArrayList<Integer> lengths = itsFeatureSet.getFeature2().getLengths();
					for(int j = 0; j < lengths.size(); j++)
					{
						outputer.print(lengths.get(j) + ",");
					}

					ArrayList<Double> ICs = itsFeatureSet.getFeature3().getICs();
					for(int j = 0; j < ICs.size(); j++)
					{
						outputer.print(ICs.get(j) + ",");
					}

					ArrayList<Integer> numbersOfSites = itsFeatureSet.getFeature4().getNumbersOfSites();
					for(int j = 0; j < numbersOfSites.size(); j++)
					{
						outputer.print(numbersOfSites.get(j) + ",");
					}

					ArrayList<Integer> numbersOfStrongSites = itsFeatureSet.getFeature5().getNumbersOfStrongSites();
					for(int j = 0; j < numbersOfStrongSites.size(); j++)
					{
						outputer.print(numbersOfStrongSites.get(j) + ",");
					}
					
					ArrayList<Double> sumsOfICOfStrongSites = itsFeatureSet.getFeature6().getSumsOfICOfSites();
					for(int j = 0; j < sumsOfICOfStrongSites.size(); j++)
					{					
						outputer.print(sumsOfICOfStrongSites.get(j) + ",");					
					}

					String label = "1";				
					outputer.print(label + "\n");
				}
			}
						
			for(int i = 0; i < geneNamesOfTrueNegatives.size(); i++)
			{
				String aGeneName = geneNamesOfTrueNegatives.get(i);
				Gene_oneiPWM anInstance = featureSetsOfLongestTranscripts.get(aGeneName);
				if(anInstance == null)
				{
					outputer.print(aGeneName + ",");
					for(int j = 1; j < numberOfColumns - 1; j++)
					{
						outputer.print(0 + ",");
					}
					outputer.print("0\n");
				}
				else
				{
					outputer.print(anInstance.getName() + ",");
					FeatureSetOfATSS_oneiPWM itsFeatureSet = anInstance.getFeatureSets().get(0);

					ArrayList<Integer> distances = itsFeatureSet.getFeature1().getDistances();
					for(int j = 0; j < distances.size(); j++)
					{
						outputer.print(distances.get(j) + ",");
					}

					ArrayList<Integer> lengths = itsFeatureSet.getFeature2().getLengths();
					for(int j = 0; j < lengths.size(); j++)
					{
						outputer.print(lengths.get(j) + ",");
					}

					ArrayList<Double> ICs = itsFeatureSet.getFeature3().getICs();
					for(int j = 0; j < ICs.size(); j++)
					{
						outputer.print(ICs.get(j) + ",");
					}

					ArrayList<Integer> numbersOfSites = itsFeatureSet.getFeature4().getNumbersOfSites();
					for(int j = 0; j < numbersOfSites.size(); j++)
					{
						outputer.print(numbersOfSites.get(j) + ",");
					}

					ArrayList<Integer> numbersOfStrongSites = itsFeatureSet.getFeature5().getNumbersOfStrongSites();
					for(int j = 0; j < numbersOfStrongSites.size(); j++)
					{
						outputer.print(numbersOfStrongSites.get(j) + ",");
					}
					
					ArrayList<Double> sumsOfICOfStrongSites = itsFeatureSet.getFeature6().getSumsOfICOfSites();
					for(int j = 0; j < sumsOfICOfStrongSites.size(); j++)
					{					
						outputer.print(sumsOfICOfStrongSites.get(j) + ",");					
					}

					String label = "0";				
					outputer.print(label + "\n");
				}
			}			
			outputer.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open the file: " + outputFilePath);
			System.exit(1);
		}
	}
	
	private static void padZeroAtTheEnd_oneiPWM(TreeMap<String, Gene_oneiPWM> featureSetsOfAllGenes, int maxNumberOfClusters)
	{
		Set<Map.Entry<String, Gene_oneiPWM>> setView = featureSetsOfAllGenes.entrySet();
		Iterator<Map.Entry<String, Gene_oneiPWM>> it = setView.iterator();
		while(it.hasNext())
		{
			Map.Entry<String, Gene_oneiPWM> anEntry = it.next();
			FeatureSetOfATSS_oneiPWM featureSet = anEntry.getValue().getFeatureSets().get(0);
			ArrayList<Integer> featureVectorOfFeature1 = featureSet.getFeature1().getDistances();
			int numberOfClustersInThisFeatureSet = featureVectorOfFeature1.size();
			if(numberOfClustersInThisFeatureSet != maxNumberOfClusters)
			{
				int numberOf0sToBePadded = maxNumberOfClusters - numberOfClustersInThisFeatureSet;
				padZeroForFeatures1245(featureVectorOfFeature1, numberOf0sToBePadded);
				padZeroForFeatures1245(featureSet.getFeature2().getLengths(), numberOf0sToBePadded);
				padZeroForFeatures36(featureSet.getFeature3().getICs(), numberOf0sToBePadded);
				padZeroForFeatures1245(featureSet.getFeature4().getNumbersOfSites(), numberOf0sToBePadded);
				padZeroForFeatures1245(featureSet.getFeature5().getNumbersOfStrongSites(), numberOf0sToBePadded);
				padZeroForFeatures36(featureSet.getFeature6().getSumsOfICOfSites(), numberOf0sToBePadded);
			}
		}
	}
	
	private static int obtainMaxNumberOfClusters_oneiPWM(TreeMap<String, Gene_oneiPWM> featureSetsOfAllGenes)
	{
		int maxNumberOfClusters = Integer.MIN_VALUE;
		Set<Map.Entry<String, Gene_oneiPWM>> setView = featureSetsOfAllGenes.entrySet();
		Iterator<Map.Entry<String, Gene_oneiPWM>> it = setView.iterator();
		while(it.hasNext())
		{
			Map.Entry<String, Gene_oneiPWM> anEntry = it.next();
			FeatureSetOfATSS_oneiPWM featureSet = anEntry.getValue().getFeatureSets().get(0);
			int numberOfClustersInThisFeatureSet = featureSet.getFeature1().getDistances().size();
			if(numberOfClustersInThisFeatureSet > maxNumberOfClusters)
			{
				maxNumberOfClusters = numberOfClustersInThisFeatureSet;
			}
		}
		
		return maxNumberOfClusters;
	}
	
	private static TreeMap<String, Gene_oneiPWM> getFeatureSetsOfLongestTranscripts_oneiPWM(ArrayList<Gene_oneiPWM> featureSetsOfAllGenes)
	{
		TreeMap<String, Gene_oneiPWM> featureSetsOfLongestTranscripts = new TreeMap<String, Gene_oneiPWM>();
		
		int currentMaxLength = Integer.MIN_VALUE;
		FeatureSetOfATSS_oneiPWM featureSetOfLongestTranscriptOfAGene = null;
		Iterator<Gene_oneiPWM> it = featureSetsOfAllGenes.iterator();
		while(it.hasNext())
		{
			currentMaxLength = Integer.MIN_VALUE;
			featureSetOfLongestTranscriptOfAGene = null;
			
			Gene_oneiPWM aGene = it.next();
			ArrayList<FeatureSetOfATSS_oneiPWM> featureSetsOfAGene = aGene.getFeatureSets();
			for(int i = 0; i < featureSetsOfAGene.size(); i++)
			{
				ArrayList<Integer> lengths = featureSetsOfAGene.get(i).getLengths();
				for(int j = 0; j < lengths.size(); j++)
				{
					int lengthOfOneTranscript = lengths.get(j);
					if(lengthOfOneTranscript > currentMaxLength)
					{
						featureSetOfLongestTranscriptOfAGene = featureSetsOfAGene.get(i);
						currentMaxLength = lengthOfOneTranscript;
					}
				}
			}
			
			Gene_oneiPWM aGeneOnlyContainingFeatureSetOfLongestTranscript = new Gene_oneiPWM(aGene.getName());
			aGeneOnlyContainingFeatureSetOfLongestTranscript.addATranscript(featureSetOfLongestTranscriptOfAGene);
			featureSetsOfLongestTranscripts.put(aGene.getName(), aGeneOnlyContainingFeatureSetOfLongestTranscript);
		}
		
		return featureSetsOfLongestTranscripts;
	}
	
	private static ArrayList<Gene_oneiPWM> getFeatureSetsOfAllGenes_oneiPWM(String featureFolderPath)
	{
		ArrayList<Gene_oneiPWM> allInstances = new ArrayList<Gene_oneiPWM>();		
			
		File[] featureFiles = new File(featureFolderPath).listFiles();		
		for(int i = 0; i < featureFiles.length; i++)
		{
			/*String featureSetFilename = featureFiles[i].getName();
			String geneName = featureSetFilename.substring(0, featureSetFilename.length() - 4);	*/
			allInstances.add(restoreFeaturesOfOneGene_oneiPWM(featureFiles[i]));
		}			
		
		return allInstances;
	}
	
	private static TreeMap<String, Gene> getFeatureSetsOfLongestTranscripts(ArrayList<Gene> featureSetsOfAllGenes)
	{
		TreeMap<String, Gene> featureSetsOfLongestTranscripts = new TreeMap<String, Gene>();
		
		int currentMaxLength = Integer.MIN_VALUE;
		FeatureSetOfATSS featureSetOfLongestTranscriptOfAGene = null;
		Iterator<Gene> it = featureSetsOfAllGenes.iterator();
		while(it.hasNext())
		{
			currentMaxLength = Integer.MIN_VALUE;
			featureSetOfLongestTranscriptOfAGene = null;
			
			Gene aGene = it.next();
			ArrayList<FeatureSetOfATSS> featureSetsOfAGene = aGene.getFeatureSets();
			for(int i = 0; i < featureSetsOfAGene.size(); i++)
			{
				ArrayList<Integer> lengths = featureSetsOfAGene.get(i).getLengths();
				for(int j = 0; j < lengths.size(); j++)
				{
					int lengthOfOneTranscript = lengths.get(j);
					if(lengthOfOneTranscript > currentMaxLength)
					{
						featureSetOfLongestTranscriptOfAGene = featureSetsOfAGene.get(i);
						currentMaxLength = lengthOfOneTranscript;
					}
				}
			}
			
			Gene aGeneOnlyContainingFeatureSetOfLongestTranscript = new Gene(aGene.getName());
			aGeneOnlyContainingFeatureSetOfLongestTranscript.addATranscript(featureSetOfLongestTranscriptOfAGene);
			featureSetsOfLongestTranscripts.put(aGene.getName(), aGeneOnlyContainingFeatureSetOfLongestTranscript);
		}
		
		return featureSetsOfLongestTranscripts;
	}
	
	private static ArrayList<String> readSimFile(String simFilePath)
	{
		ArrayList<String> geneNamesOfTruePositivesOrNegatives = new ArrayList<String>();
		try
		{
			Scanner inputer = new Scanner(new File(simFilePath));
			while(inputer.hasNextLine())
			{
				String aLine = inputer.nextLine();
				String[] fields = aLine.split(",");
				geneNamesOfTruePositivesOrNegatives.add(fields[1]);
			}
			inputer.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open the file: " + simFilePath);
			System.exit(1);
		}
		
		return geneNamesOfTruePositivesOrNegatives;
	}	
	
	private static void outputFeatureVectorsForTPTN(TreeMap<String, Gene> featureSetsOfLongestTranscripts, String outputFilePath, int maxNumberOfClusters, ArrayList<String> geneNamesOfTruePositives, ArrayList<String> geneNamesOfTrueNegatives)
	{
		try
		{
			PrintWriter outputer = new PrintWriter(outputFilePath);
			int numberOfColumns = 1 + maxNumberOfClusters * 3 + 82 * maxNumberOfClusters * 4 + 1;
			for(int i = 0; i < numberOfColumns - 1; i++)
			{
				outputer.print((i + 1) + ",");
			}
			outputer.print("Label\n");
			
			for(int i = 0; i < geneNamesOfTruePositives.size(); i++)
			{
				String aGeneName = geneNamesOfTruePositives.get(i);
				Gene anInstance = featureSetsOfLongestTranscripts.get(aGeneName);
				if(anInstance == null)
				{
					outputer.print(aGeneName + ",");
					for(int j = 1; j < numberOfColumns - 1; j++)
					{
						outputer.print(0 + ",");
					}
					outputer.print("1\n");
				}
				else
				{
					outputer.print(anInstance.getName() + ",");
					FeatureSetOfATSS itsFeatureSet = anInstance.getFeatureSets().get(0);

					ArrayList<Integer> distances = itsFeatureSet.getFeature1().getDistances();
					for(int j = 0; j < distances.size(); j++)
					{
						outputer.print(distances.get(j) + ",");
					}

					ArrayList<Integer> lengths = itsFeatureSet.getFeature2().getLengths();
					for(int j = 0; j < lengths.size(); j++)
					{
						outputer.print(lengths.get(j) + ",");
					}

					ArrayList<Double> ICs = itsFeatureSet.getFeature3().getICs();
					for(int j = 0; j < ICs.size(); j++)
					{
						outputer.print(ICs.get(j) + ",");
					}

					ArrayList<int[]> numbersOfSites = itsFeatureSet.getFeature4().getNumbersOfSites();
					for(int j = 0; j < numbersOfSites.size(); j++)
					{
						int[] aCluster = numbersOfSites.get(j);
						for(int k = 0; k < aCluster.length; k++)
						{
							outputer.print(aCluster[k] + ",");
						}
					}

					ArrayList<int[]> numbersOfStrongSites = itsFeatureSet.getFeature5().getNumbersOfStrongSites();
					for(int j = 0; j < numbersOfStrongSites.size(); j++)
					{
						int[] aCluster = numbersOfStrongSites.get(j);
						for(int k = 0; k < aCluster.length; k++)
						{
							outputer.print(aCluster[k] + ",");
						}
					}

					ArrayList<double[]> sumsOfICOfSites = itsFeatureSet.getFeature6().getSumsOfICOfSites();
					for(int j = 0; j < sumsOfICOfSites.size(); j++)
					{
						double[] aCluster = sumsOfICOfSites.get(j);
						for(int k = 0; k < aCluster.length; k++)
						{
							outputer.print(aCluster[k] + ",");
						}
					}

					ArrayList<double[]> sumsOfICOfStrongSites = itsFeatureSet.getFeature7().getSumsOfICOfSites();
					for(int j = 0; j < sumsOfICOfStrongSites.size(); j++)
					{					
						double[] aCluster = sumsOfICOfStrongSites.get(j);
						for(int k = 0; k < aCluster.length; k++)
						{
							outputer.print(aCluster[k] + ",");
						}					
					}

					String label = "1";				
					outputer.print(label + "\n");
				}
			}
			
			for(int i = 0; i < geneNamesOfTrueNegatives.size(); i++)
			{
				String aGeneName = geneNamesOfTrueNegatives.get(i);
				Gene anInstance = featureSetsOfLongestTranscripts.get(aGeneName);
				if(anInstance == null)
				{
					outputer.print(aGeneName + ",");
					for(int j = 1; j < numberOfColumns - 1; j++)
					{
						outputer.print(0 + ",");
					}
					outputer.print("0\n");
				}
				else
				{
					outputer.print(anInstance.getName() + ",");
					FeatureSetOfATSS itsFeatureSet = anInstance.getFeatureSets().get(0);

					ArrayList<Integer> distances = itsFeatureSet.getFeature1().getDistances();
					for(int j = 0; j < distances.size(); j++)
					{
						outputer.print(distances.get(j) + ",");
					}

					ArrayList<Integer> lengths = itsFeatureSet.getFeature2().getLengths();
					for(int j = 0; j < lengths.size(); j++)
					{
						outputer.print(lengths.get(j) + ",");
					}

					ArrayList<Double> ICs = itsFeatureSet.getFeature3().getICs();
					for(int j = 0; j < ICs.size(); j++)
					{
						outputer.print(ICs.get(j) + ",");
					}

					ArrayList<int[]> numbersOfSites = itsFeatureSet.getFeature4().getNumbersOfSites();
					for(int j = 0; j < numbersOfSites.size(); j++)
					{
						int[] aCluster = numbersOfSites.get(j);
						for(int k = 0; k < aCluster.length; k++)
						{
							outputer.print(aCluster[k] + ",");
						}
					}

					ArrayList<int[]> numbersOfStrongSites = itsFeatureSet.getFeature5().getNumbersOfStrongSites();
					for(int j = 0; j < numbersOfStrongSites.size(); j++)
					{
						int[] aCluster = numbersOfStrongSites.get(j);
						for(int k = 0; k < aCluster.length; k++)
						{
							outputer.print(aCluster[k] + ",");
						}
					}

					ArrayList<double[]> sumsOfICOfSites = itsFeatureSet.getFeature6().getSumsOfICOfSites();
					for(int j = 0; j < sumsOfICOfSites.size(); j++)
					{
						double[] aCluster = sumsOfICOfSites.get(j);
						for(int k = 0; k < aCluster.length; k++)
						{
							outputer.print(aCluster[k] + ",");
						}
					}

					ArrayList<double[]> sumsOfICOfStrongSites = itsFeatureSet.getFeature7().getSumsOfICOfSites();
					for(int j = 0; j < sumsOfICOfStrongSites.size(); j++)
					{					
						double[] aCluster = sumsOfICOfStrongSites.get(j);
						for(int k = 0; k < aCluster.length; k++)
						{
							outputer.print(aCluster[k] + ",");
						}					
					}

					String label = "0";				
					outputer.print(label + "\n");
				}
			}
			
			outputer.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open file: " + outputFilePath);
			System.exit(1);
		}
	}
	
	private static void padZeroAtTheEnd(TreeMap<String, Gene> featureSetsOfAllGenes, int maxNumberOfClusters)
	{
		Set<Map.Entry<String, Gene>> setView = featureSetsOfAllGenes.entrySet();
		Iterator<Map.Entry<String, Gene>> it = setView.iterator();
		while(it.hasNext())
		{
			Map.Entry<String, Gene> anEntry = it.next();
			FeatureSetOfATSS featureSet = anEntry.getValue().getFeatureSets().get(0);
			ArrayList<Integer> featureVectorOfFeature1 = featureSet.getFeature1().getDistances();
			int numberOfClustersInThisFeatureSet = featureVectorOfFeature1.size();
			if(numberOfClustersInThisFeatureSet != maxNumberOfClusters)
			{
				int numberOf0sToBePadded = maxNumberOfClusters - numberOfClustersInThisFeatureSet;
				padZeroForFeatures12(featureVectorOfFeature1, numberOf0sToBePadded);
				padZeroForFeatures12(featureSet.getFeature2().getLengths(), numberOf0sToBePadded);
				padZeroForFeatures3(featureSet.getFeature3().getICs(), numberOf0sToBePadded);
				padZeroForFeatures45(featureSet.getFeature4().getNumbersOfSites(), numberOf0sToBePadded);
				padZeroForFeatures45(featureSet.getFeature5().getNumbersOfStrongSites(), numberOf0sToBePadded);
				padZeroForFeatures67(featureSet.getFeature6().getSumsOfICOfSites(), numberOf0sToBePadded);
				padZeroForFeatures67(featureSet.getFeature7().getSumsOfICOfSites(), numberOf0sToBePadded);
			}
		}
	}
	
	private static void padZeroForFeatures67(ArrayList<double[]> featureVector, int numberOf0sToBePadded)
	{
		double[] nullArrayToBePadded = new double[82];
		for(int i = 0; i < 82; i++)
		{
			nullArrayToBePadded[i] = 0.0;
		}
		
		for(int i = 0; i < numberOf0sToBePadded; i++)
		{
			featureVector.add(nullArrayToBePadded);
		}
	}
	
	private static void padZeroForFeatures45(ArrayList<int[]> featureVector, int numberOf0sToBePadded)
	{
		int[] nullArrayToBePadded = new int[82];
		for(int i = 0; i < 82; i++)
		{
			nullArrayToBePadded[i] = 0;
		}
		
		for(int i = 0; i < numberOf0sToBePadded; i++)
		{
			featureVector.add(nullArrayToBePadded);
		}
	}
	
	private static void padZeroForFeatures3(ArrayList<Double> featureVector, int numberOf0sToBePadded)
	{
		for(int i = 0; i < numberOf0sToBePadded; i++)
		{
			featureVector.add(0.0);
		}
	}
	
	private static void padZeroForFeatures12(ArrayList<Integer> featureVector, int numberOf0sToBePadded)
	{
		for(int i = 0; i < numberOf0sToBePadded; i++)
		{
			featureVector.add(0);
		}
	}
	
	private static int obtainMaxNumberOfClusters(TreeMap<String, Gene> featureSetsOfAllGenes)
	{
		int maxNumberOfClusters = Integer.MIN_VALUE;
		Set<Map.Entry<String, Gene>> setView = featureSetsOfAllGenes.entrySet();
		Iterator<Map.Entry<String, Gene>> it = setView.iterator();
		while(it.hasNext())
		{
			Map.Entry<String, Gene> anEntry = it.next();
			FeatureSetOfATSS featureSet = anEntry.getValue().getFeatureSets().get(0);
			int numberOfClustersInThisFeatureSet = featureSet.getFeature1().getDistances().size();
			if(numberOfClustersInThisFeatureSet > maxNumberOfClusters)
			{
				maxNumberOfClusters = numberOfClustersInThisFeatureSet;
			}
		}
		
		return maxNumberOfClusters;
	}
	
	private static ArrayList<Gene> getFeatureSetsOfAllGenes(String featureFolderPath)
	{
		ArrayList<Gene> allInstances = new ArrayList<Gene>();		
			
		File[] featureFiles = new File(featureFolderPath).listFiles();		
		for(int i = 0; i < featureFiles.length; i++)
		{
			allInstances.add(restoreFeaturesOfOneGene(featureFiles[i]));
		}			
		
		return allInstances;
	}
	
	private static Gene_oneiPWM restoreFeaturesOfOneGene_oneiPWM(File featureFileOfOneGene)
	{
		String geneName = null;
		ArrayList<FeatureSetOfATSS_oneiPWM> featureSets = new ArrayList<FeatureSetOfATSS_oneiPWM>();
		
		try
		{
			Scanner input = new Scanner(featureFileOfOneGene);
			geneName = input.nextLine().trim();
			
			ArrayList<String> linesInOneTranscript = new ArrayList<String>();
			boolean transcriptStart = false;
			while(input.hasNextLine())
			{
				String aLine = input.nextLine();
				if(!transcriptStart)
				{
					if(aLine.contains("***********************************"))
					{
						transcriptStart = true;
						continue;
					}
				}
				if(transcriptStart)
				{
					if(aLine.contains("***********************************"))
					{
						transcriptStart = false;
						FeatureSetOfATSS_oneiPWM aFeatureSet = restoreFeatureSetOfOneTSS_oneiPWM(linesInOneTranscript);
						featureSets.add(aFeatureSet);
												
						linesInOneTranscript.clear();
					}
					else
					{
						linesInOneTranscript.add(aLine);
					}
				}
			}			
			input.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open the file: " + featureFileOfOneGene.getAbsolutePath());
			e.printStackTrace();
			System.exit(1);
		}
		
		return new Gene_oneiPWM(geneName, featureSets);
	}
	
	private static FeatureSetOfATSS_oneiPWM restoreFeatureSetOfOneTSS_oneiPWM(ArrayList<String> linesInOneTranscript)
	{
		String transcripts = linesInOneTranscript.get(0).split("[ :]")[1].trim();
		String transcriptLengths = linesInOneTranscript.get(1).split("[ :]")[2].trim();		
		
		ArrayList<Integer> distances = new ArrayList<Integer>();
		ArrayList<Integer> lengths = new ArrayList<Integer>();
		ArrayList<Double> ICs = new ArrayList<Double>();
		ArrayList<Integer> numbersOfSites = new ArrayList<Integer>();
		ArrayList<Integer> numbersOfStrongSites = new ArrayList<Integer>();
		ArrayList<Double> sumsOfICOfStrongSites = new ArrayList<Double>();
		
		for(int i = 1; i < linesInOneTranscript.size(); i++)
		{
			String aLine = linesInOneTranscript.get(i);
			if(aLine.contains("Feature 1"))
			{
				String lineFeature1 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature1.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					distances.add(Integer.parseInt(fields[j].trim()));
				}				
			}
			else if(aLine.contains("Feature 2"))
			{
				String lineFeature2 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature2.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					lengths.add(Integer.parseInt(fields[j].trim()));
				}
			}
			else if(aLine.contains("Feature 3"))
			{
				String lineFeature3 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature3.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					ICs.add(Double.parseDouble(fields[j].trim()));
				}
			}
			else if(aLine.contains("Feature 4"))
			{
				String lineFeature4 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature4.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					numbersOfSites.add(Integer.parseInt(fields[j].trim()));
				}
			}
			else if(aLine.contains("Feature 5"))
			{
				String lineFeature5 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature5.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					numbersOfStrongSites.add(Integer.parseInt(fields[j].trim()));
				}
			}
			else if(aLine.contains("Feature 6"))
			{
				String lineFeature6 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature6.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					sumsOfICOfStrongSites.add(Double.parseDouble(fields[j].trim()));
				}
			}
		}
		
		return new FeatureSetOfATSS_oneiPWM(transcripts, transcriptLengths, new FeatureDistanceBetweenClusterAndTSS(distances), new FeatureClusterLength(lengths), new FeatureClusterIC(ICs), new FeatureNumberSitesTFCluster_oneiPWM(numbersOfSites), new FeatureNumberStrongSitesTFCluster_oneiPWM(numbersOfStrongSites), new FeatureSumICStrongSitesTFCluster_oneiPWM(sumsOfICOfStrongSites));
	}
	
	private static Gene restoreFeaturesOfOneGene(File featureFileOfOneGene)
	{
		String geneName = null;
		ArrayList<FeatureSetOfATSS> featureSets = new ArrayList<FeatureSetOfATSS>();
		
		try
		{
			Scanner input = new Scanner(featureFileOfOneGene);
			geneName = input.nextLine().trim();
			
			ArrayList<String> linesInOneTranscript = new ArrayList<String>();
			boolean transcriptStart = false;
			while(input.hasNextLine())
			{
				String aLine = input.nextLine();
				if(!transcriptStart)
				{
					if(aLine.contains("***********************************"))
					{
						transcriptStart = true;
						continue;
					}
				}
				if(transcriptStart)
				{
					if(aLine.contains("***********************************"))
					{
						transcriptStart = false;
						FeatureSetOfATSS aFeatureSet = restoreFeatureSetOfOneTSS(linesInOneTranscript);
						featureSets.add(aFeatureSet);
												
						linesInOneTranscript.clear();
					}
					else
					{
						linesInOneTranscript.add(aLine);
					}
				}
			}			
			input.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open the file: " + featureFileOfOneGene.getAbsolutePath());
			e.printStackTrace();
			System.exit(1);
		}
		
		return new Gene(geneName, featureSets);
	}
	
	private static FeatureSetOfATSS restoreFeatureSetOfOneTSS(ArrayList<String> linesInOneTranscript)
	{
		String transcripts = linesInOneTranscript.get(0).split("[ :]")[1].trim();
		String transcriptLengths = linesInOneTranscript.get(1).split("[ :]")[2].trim();		
		
		ArrayList<Integer> distances = new ArrayList<Integer>();
		ArrayList<Integer> lengths = new ArrayList<Integer>();
		ArrayList<Double> ICs = new ArrayList<Double>();
		ArrayList<int[]> numbersOfSites = new ArrayList<int[]>();
		ArrayList<int[]> numbersOfStrongSites = new ArrayList<int[]>();
		ArrayList<double[]> sumsOfICOfSites = new ArrayList<double[]>();
		ArrayList<double[]> sumsOfICOfStrongSites = new ArrayList<double[]>();
		
		for(int i = 1; i < linesInOneTranscript.size(); i++)
		{
			String aLine = linesInOneTranscript.get(i);
			if(aLine.contains("Feature 1"))
			{
				String lineFeature1 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature1.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					distances.add(Integer.parseInt(fields[j].trim()));
				}				
			}
			else if(aLine.contains("Feature 2"))
			{
				String lineFeature2 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature2.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					lengths.add(Integer.parseInt(fields[j].trim()));
				}
			}
			else if(aLine.contains("Feature 3"))
			{
				String lineFeature3 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature3.split("[\\[,\\]]");
				for(int j = 1; j < fields.length; j++)
				{
					ICs.add(Double.parseDouble(fields[j].trim()));
				}
			}
			else if(aLine.contains("Feature 4"))
			{
				String lineFeature4 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature4.split(";");
				for(int j = 0; j < fields.length; j++)
				{
					int[] aCluster = new int[82];
					String[] fieldsInOneCluster = fields[j].split("[\\[,\\]]");
					if(j == 0)
					{
						for(int k = 2; k < fieldsInOneCluster.length; k++)
						{
							aCluster[k - 2] = Integer.parseInt(fieldsInOneCluster[k].trim());
						}
					}
					else
					{
						for(int k = 1; k < fieldsInOneCluster.length; k++)
						{
							aCluster[k - 1] = Integer.parseInt(fieldsInOneCluster[k].trim());
						}
					}					
					
					numbersOfSites.add(aCluster);
				}
			}
			else if(aLine.contains("Feature 5"))
			{
				String lineFeature5 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature5.split(";");
				for(int j = 0; j < fields.length; j++)
				{
					int[] aCluster = new int[82];
					String[] fieldsInOneCluster = fields[j].split("[\\[,\\]]");
					if(j == 0)
					{
						for(int k = 2; k < fieldsInOneCluster.length; k++)
						{
							aCluster[k - 2] = Integer.parseInt(fieldsInOneCluster[k].trim());
						}
					}
					else
					{
						for(int k = 1; k < fieldsInOneCluster.length; k++)
						{
							aCluster[k - 1] = Integer.parseInt(fieldsInOneCluster[k].trim());
						}
					}
					
					numbersOfStrongSites.add(aCluster);
				}
			}
			else if(aLine.contains("Feature 6"))
			{
				String lineFeature6 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature6.split(";");
				for(int j = 0; j < fields.length; j++)
				{
					double[] aCluster = new double[82];
					String[] fieldsInOneCluster = fields[j].split("[\\[,\\]]");
					if(j == 0)
					{
						for(int k = 2; k < fieldsInOneCluster.length; k++)
						{
							aCluster[k - 2] = Double.parseDouble(fieldsInOneCluster[k].trim());
						}
					}
					else
					{
						for(int k = 1; k < fieldsInOneCluster.length; k++)
						{
							aCluster[k - 1] = Double.parseDouble(fieldsInOneCluster[k].trim());
						}
					}
					
					sumsOfICOfSites.add(aCluster);
				}
			}
			else if(aLine.contains("Feature 7"))
			{
				String lineFeature7 = linesInOneTranscript.get(i + 1);
				String[] fields = lineFeature7.split(";");
				for(int j = 0; j < fields.length; j++)
				{
					double[] aCluster = new double[82];
					String[] fieldsInOneCluster = fields[j].split("[\\[,\\]]");
					if(j == 0)
					{
						for(int k = 2; k < fieldsInOneCluster.length; k++)
						{
							aCluster[k - 2] = Double.parseDouble(fieldsInOneCluster[k].trim());
						}
					}
					else
					{
						for(int k = 1; k < fieldsInOneCluster.length; k++)
						{
							aCluster[k - 1] = Double.parseDouble(fieldsInOneCluster[k].trim());
						}
					}
					
					sumsOfICOfStrongSites.add(aCluster);
				}
			}			
		}
		
		return new FeatureSetOfATSS(transcripts, transcriptLengths, new FeatureDistanceBetweenClusterAndTSS(distances), new FeatureClusterLength(lengths), new FeatureClusterIC(ICs), new FeatureNumberSitesTFCluster(numbersOfSites), new FeatureNumberStrongSitesTFCluster(numbersOfStrongSites), new FeatureSumICSitesTFCluster(sumsOfICOfSites), new FeatureSumICStrongSitesTFCluster(sumsOfICOfStrongSites)/*, new FeatureOrderOfSitesCluster(orderOfSites), new FeatureOrderOfStrongSitesCluster(orderOfStrongSites)*/);
	}		
}
