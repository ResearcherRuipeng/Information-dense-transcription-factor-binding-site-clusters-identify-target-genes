import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;

public class Main
{
	public static void main(String[] args)
	{
		String clusterFolder = args[0];
		String allGenesFile = args[1];
		String promoterRegionsPath = args[2];
		String modelPath = args[3];		
		double threshold = Double.parseDouble(args[4]);
		double strongThreshold = Double.parseDouble(args[5]);
		String featureSetFolder = args[6];
		boolean ifAbsoluteCoor = Boolean.parseBoolean(args[7]);;
		boolean ifOneiPWM = Boolean.parseBoolean(args[8]);
		
		TreeSet<String> geneNames = getNamesOfGenesForWhichClustersWereObtained(clusterFolder);
		TreeMap<String, ArrayList<String>> geneAndTranscripts = getTranscriptsInfo(allGenesFile);
		TreeMap<String, ArrayList<File>> geneAndPromoterRegions = getPromoterRegions(promoterRegionsPath);
		TreeMap<String, Integer> TFsAndIndices = getTFsAndIndices(modelPath);
		HashMap<String, Double> modelToRseq = readRseqValues(modelPath);
		System.out.println("Size of the TreeSet geneNames: " + geneNames.size() + "; Size of the TreeMap geneAndPromoterRegions: " + geneAndPromoterRegions.size());
		
		String geneName = null;
		if(!ifOneiPWM)
		{
			Gene gene = null;
			
			Iterator<String> it = geneNames.iterator();		
			while(it.hasNext())
			{
				geneName = it.next();
				gene = getFeatureSetsForOneGene(geneName, clusterFolder, geneAndTranscripts, geneAndPromoterRegions, threshold, strongThreshold, TFsAndIndices, modelToRseq, ifAbsoluteCoor);
				
				try
				{
					PrintWriter output = new PrintWriter(featureSetFolder + "/" + geneName + ".txt");
					output.println(gene.toString());
					output.close();
				}
				catch(FileNotFoundException e)
				{
					System.out.println("The output file " + geneName + ".txt cannot be opened");
					System.exit(1);
				}
			}
		}
		else
		{
			Gene_oneiPWM gene = null;
			
			Iterator<String> it = geneNames.iterator();		
			while(it.hasNext())
			{
				geneName = it.next();
				gene = getFeatureSetsForOneGene_OneiPWM(geneName, clusterFolder, geneAndTranscripts, geneAndPromoterRegions, threshold, strongThreshold, TFsAndIndices, modelToRseq, ifAbsoluteCoor);
				
				try
				{
					PrintWriter output = new PrintWriter(featureSetFolder + "/" + geneName + ".txt");
					output.println(gene.toString());
					output.close();
				}
				catch(FileNotFoundException e)
				{
					System.out.println("The output file " + geneName + ".txt cannot be opened");
					System.exit(1);
				}
			}
		}
		
		System.out.println("Done");
	}
	
	private static FeatureSumICStrongSitesTFCluster_oneiPWM extractFeature6OfATranscript_OneiPWM(String strand, ArrayList<Cluster> clustersOfOneTranscript, double strongThreshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq)
	{
		ArrayList<Double> sumsOfICOfStrongSites = new ArrayList<Double>();		
		Cluster aCluster = null;
		HashMap<String, Double> TFsAndICs = null;
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				double sumsOfICOfStrongSitesInOneCluster = 0.0;
				aCluster = clustersOfOneTranscript.get(i);
				TFsAndICs = aCluster.getCenterSite().getTFsAndICs(modelToRseq, strongThreshold);
				if(TFsAndICs.size() == 1)
				{					
					sumsOfICOfStrongSitesInOneCluster += TFsAndICs.entrySet().iterator().next().getValue();
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFsAndICs = aSite.getTFsAndICs(modelToRseq, strongThreshold);
					if(TFsAndICs.size() == 1)
					{
						sumsOfICOfStrongSitesInOneCluster += TFsAndICs.entrySet().iterator().next().getValue();
					}
				}
				
				sumsOfICOfStrongSites.add(sumsOfICOfStrongSitesInOneCluster);
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				double sumsOfICOfStrongSitesInOneCluster = 0.0;
				aCluster = clustersOfOneTranscript.get(i);
				TFsAndICs = aCluster.getCenterSite().getTFsAndICs(modelToRseq, strongThreshold);
				if(TFsAndICs.size() == 1)
				{					
					sumsOfICOfStrongSitesInOneCluster += TFsAndICs.entrySet().iterator().next().getValue();
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFsAndICs = aSite.getTFsAndICs(modelToRseq, strongThreshold);
					if(TFsAndICs.size() == 1)
					{
						sumsOfICOfStrongSitesInOneCluster += TFsAndICs.entrySet().iterator().next().getValue();
					}
				}
				
				sumsOfICOfStrongSites.add(sumsOfICOfStrongSitesInOneCluster);
			}
		}
		
		return new FeatureSumICStrongSitesTFCluster_oneiPWM(sumsOfICOfStrongSites);
	}
	
	private static FeatureNumberStrongSitesTFCluster_oneiPWM extractFeature5OfATranscript_OneiPWM(String strand, ArrayList<Cluster> clustersOfOneTranscript, double strongThreshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq)
	{
		ArrayList<Integer> numbersOfStrongSites = new ArrayList<Integer>();
		int numberOfStrongSitesInOneCluster = 0;
		Cluster aCluster = null;
		TreeSet<String> TFs = null;
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				numberOfStrongSitesInOneCluster = 0;
				aCluster = clustersOfOneTranscript.get(i);
				TFs = aCluster.getCenterSite().getTFs(modelToRseq, strongThreshold);
				if(TFs.size() == 1)
				{
					numberOfStrongSitesInOneCluster++;
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFs = aSite.getTFs(modelToRseq, strongThreshold);
					if(TFs.size() == 1)
					{
						numberOfStrongSitesInOneCluster++;
					}
				}
				
				numbersOfStrongSites.add(numberOfStrongSitesInOneCluster);
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				numberOfStrongSitesInOneCluster = 0;
				aCluster = clustersOfOneTranscript.get(i);
				TFs = aCluster.getCenterSite().getTFs(modelToRseq, strongThreshold);
				if(TFs.size() == 1)
				{
					numberOfStrongSitesInOneCluster++;
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFs = aSite.getTFs(modelToRseq, strongThreshold);
					if(TFs.size() == 1)
					{
						numberOfStrongSitesInOneCluster++;
					}
				}
				
				numbersOfStrongSites.add(numberOfStrongSitesInOneCluster);
			}
		}		
		
		return new FeatureNumberStrongSitesTFCluster_oneiPWM(numbersOfStrongSites);
	}
	
	private static FeatureNumberSitesTFCluster_oneiPWM extractFeature4OfATranscript_OneiPWM(String strand, ArrayList<Cluster> clustersOfOneTranscript, double threshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq)
	{
		ArrayList<Integer> numbersOfSites = new ArrayList<Integer>();
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				numbersOfSites.add(clustersOfOneTranscript.get(i).getNumberOfSites());
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				numbersOfSites.add(clustersOfOneTranscript.get(i).getNumberOfSites());
			}
		}		
		
		return new FeatureNumberSitesTFCluster_oneiPWM(numbersOfSites);
	}
	
	private static FeatureSetOfATSS_oneiPWM ifTheTSSWasEncountered(int transcriptTSSNumbering, Gene_oneiPWM aGene)
	{
		ArrayList<FeatureSetOfATSS_oneiPWM> featureSets = aGene.getFeatureSets();
		for(int i = 0; i < featureSets.size(); i++)
		{
			FeatureSetOfATSS_oneiPWM aFeatureSet = featureSets.get(i);
			if(aFeatureSet.getTSS() == transcriptTSSNumbering)
			{
				return aFeatureSet;
			}
		}
		
		return null;
	}
	
	private static Gene_oneiPWM getFeatureSetsForOneGene_OneiPWM(String geneName, String clusterFolder, TreeMap<String, ArrayList<String>> geneAndTranscripts, TreeMap<String, ArrayList<File>> geneAndPromoterRegions, double threshold, double strongThreshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq, boolean ifAbsoluteCoor) /*throws Exception*/
	{
		Gene_oneiPWM aGene = new Gene_oneiPWM(geneName);
		Iterator<String> it = null;
		String oneTranscript = null;
		String[] fields = null;
		ArrayList<Cluster> allClustersOfTheGene = null;
		ArrayList<Cluster> clustersOfOneTranscript = null;
		int smallCoordinatePromoterRegion = -1;
		int transcriptTSSNumbering = -1;
		ArrayList<File> promoterRegions = geneAndPromoterRegions.get(geneName);
		ArrayList<String> transcripts = geneAndTranscripts.get(geneName);
		
		if(promoterRegions.size() == 1)
		{
			allClustersOfTheGene = restoreClustersFromOneFile(clusterFolder + "/" + promoterRegions.get(0).getName());
			smallCoordinatePromoterRegion = getSmallCoordinatePromoterRegion(promoterRegions.get(0));
			it = transcripts.iterator();			
			while(it.hasNext())
			{
				oneTranscript = it.next();
				fields = oneTranscript.split(",");
				transcriptTSSNumbering = Integer.parseInt(fields[3]);
				FeatureSetOfATSS_oneiPWM aFeatureSet = ifTheTSSWasEncountered(transcriptTSSNumbering, aGene);
				if(aFeatureSet == null)
				{
					clustersOfOneTranscript = getClustersOfOneTranscript(fields[6], transcriptTSSNumbering, smallCoordinatePromoterRegion, allClustersOfTheGene, ifAbsoluteCoor);
					if(clustersOfOneTranscript.size() == 0)
					{
						continue;
					}
					FeatureDistanceBetweenClusterAndTSS feature1 = extractFeature1OfATranscript(fields[6], clustersOfOneTranscript, smallCoordinatePromoterRegion, Integer.parseInt(fields[3]), ifAbsoluteCoor);
					FeatureClusterLength feature2 = extractFeature2OfATranscript(fields[6], clustersOfOneTranscript);
					FeatureClusterIC feature3 = extractFeature3OfATranscript(fields[6], clustersOfOneTranscript);
					FeatureNumberSitesTFCluster_oneiPWM feature4 = extractFeature4OfATranscript_OneiPWM(fields[6], clustersOfOneTranscript, threshold, TFsAndIndices, modelToRseq);
					FeatureNumberStrongSitesTFCluster_oneiPWM feature5 = extractFeature5OfATranscript_OneiPWM(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
					FeatureSumICStrongSitesTFCluster_oneiPWM feature6 = extractFeature6OfATranscript_OneiPWM(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
					/*FeatureSumICStrongSitesTFCluster feature7 = extractFeature7OfATranscript(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
					FeatureOrderOfSitesCluster feature8 = extractFeature8OfATranscript(fields[6], clustersOfOneTranscript, threshold, modelToRseq);
					FeatureOrderOfStrongSitesCluster feature9 = extractFeature9OfATranscript(fields[6], clustersOfOneTranscript, strongThreshold, modelToRseq);*/
					aGene.addATranscript(new FeatureSetOfATSS_oneiPWM(transcriptTSSNumbering, fields[5], Integer.parseInt(fields[9]), feature1, feature2, feature3, feature4, feature5, feature6));
				}
				else
				{
					aFeatureSet.addATranscript(fields[5]);
					aFeatureSet.addALength(Integer.parseInt(fields[9]));
				}
			}
		}		
		else
		{
			File promoterRegion = null;
			it = transcripts.iterator();
			while(it.hasNext())
			{
				oneTranscript = it.next();
				fields = oneTranscript.split(",");
				transcriptTSSNumbering = Integer.parseInt(fields[3]);
				FeatureSetOfATSS_oneiPWM aFeatureSet = ifTheTSSWasEncountered(transcriptTSSNumbering, aGene);
				if(aFeatureSet == null)
				{
					
					promoterRegion = determinePromoterRegion(transcriptTSSNumbering, promoterRegions);
					if(promoterRegion == null)
					{
						continue;
					}
					
					/*try
					{*/
						allClustersOfTheGene = restoreClustersFromOneFile(clusterFolder + "/" + promoterRegion.getName());
						if(allClustersOfTheGene == null)
						{
							continue;
						}
						smallCoordinatePromoterRegion = getSmallCoordinatePromoterRegion(promoterRegion);
						clustersOfOneTranscript = getClustersOfOneTranscript(fields[6], transcriptTSSNumbering, smallCoordinatePromoterRegion, allClustersOfTheGene, ifAbsoluteCoor);
						if(clustersOfOneTranscript.size() == 0)
						{
							continue;
						}
						FeatureDistanceBetweenClusterAndTSS feature1 = extractFeature1OfATranscript(fields[6], clustersOfOneTranscript, smallCoordinatePromoterRegion, Integer.parseInt(fields[3]), ifAbsoluteCoor);
						FeatureClusterLength feature2 = extractFeature2OfATranscript(fields[6], clustersOfOneTranscript);
						FeatureClusterIC feature3 = extractFeature3OfATranscript(fields[6], clustersOfOneTranscript);
						FeatureNumberSitesTFCluster_oneiPWM feature4 = extractFeature4OfATranscript_OneiPWM(fields[6], clustersOfOneTranscript, threshold, TFsAndIndices, modelToRseq);
						FeatureNumberStrongSitesTFCluster_oneiPWM feature5 = extractFeature5OfATranscript_OneiPWM(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
						FeatureSumICStrongSitesTFCluster_oneiPWM feature6 = extractFeature6OfATranscript_OneiPWM(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
						/*FeatureSumICStrongSitesTFCluster feature7 = extractFeature7OfATranscript(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
						FeatureOrderOfSitesCluster feature8 = extractFeature8OfATranscript(fields[6], clustersOfOneTranscript, threshold, modelToRseq);
						FeatureOrderOfStrongSitesCluster feature9 = extractFeature9OfATranscript(fields[6], clustersOfOneTranscript, strongThreshold, modelToRseq);*/
						aGene.addATranscript(new FeatureSetOfATSS_oneiPWM(transcriptTSSNumbering, fields[5], Integer.parseInt(fields[9]), feature1, feature2, feature3, feature4, feature5, feature6));
					/*}
					catch(Exception e)
					{
						System.out.println(e.getMessage());
						continue;
					}*/
				}
				else
				{
					aFeatureSet.addATranscript(fields[5]);
					aFeatureSet.addALength(Integer.parseInt(fields[9]));
				}
			}
		}
		
		return aGene;
	}
	
	private static File determinePromoterRegion(int TSS, ArrayList<File> promoterRegions)
	{
		Scanner input = null;
		String line = null;
		int smallCoor = -1;
		int largeCoor = -1;
		String[] fields = null;
		File currentRegion = null;
		
		try
		{
			for(int i = 0; i < promoterRegions.size(); i++)
			{
				currentRegion = promoterRegions.get(i);
				input = new Scanner(currentRegion);
				line = input.nextLine();
				input.close();
				fields = line.split("\t");
				smallCoor = Integer.parseInt(fields[1]);
				largeCoor = Integer.parseInt(fields[2]);

				if((TSS >= smallCoor) && (TSS <= largeCoor + 1))
				{
					return currentRegion;
				}			
			}
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open a promoter region file.");
			System.exit(1);
		}
		
		return null;
	}
	
	private static HashMap<String, Double> readRseqValues(String modelDirName)
	{
		HashMap<String, Double> modelToRseq = new HashMap<String, Double>();
		File modelDir = new File(modelDirName);
		File[] modelFiles = modelDir.listFiles();
		
		Scanner input;
		for(int i = 0; i < modelFiles.length; i++)
		{
			String filename = modelFiles[i].getName().substring(0, modelFiles[i].getName().length() - 4);
			try
			{
				input = new Scanner(modelFiles[i]);
				String[] fields = input.nextLine().split("\\s+");
				modelToRseq.put(filename, Double.valueOf(fields[fields.length - 1]));				
			}
			catch(FileNotFoundException e)
			{
				System.out.println("Unable to open the model file: " + filename);
				System.exit(1);
			}
		}
		
		return modelToRseq;
	}
	
	private static FeatureSetOfATSS ifTheTSSWasEncountered(int transcriptTSSNumbering, Gene aGene)
	{
		ArrayList<FeatureSetOfATSS> featureSets = aGene.getFeatureSets();
		for(int i = 0; i < featureSets.size(); i++)
		{
			FeatureSetOfATSS aFeatureSet = featureSets.get(i);
			if(aFeatureSet.getTSS() == transcriptTSSNumbering)
			{
				return aFeatureSet;
			}
		}
		
		return null;
	}
	
	private static Gene getFeatureSetsForOneGene(String geneName, String clusterFolder, TreeMap<String, ArrayList<String>> geneAndTranscripts, TreeMap<String, ArrayList<File>> geneAndPromoterRegions, double threshold, double strongThreshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq, boolean ifAbsoluteCoor) /*throws Exception*/
	{
		Gene aGene = new Gene(geneName);
		Iterator<String> it = null;
		String oneTranscript = null;
		String[] fields = null;
		ArrayList<Cluster> allClustersOfTheGene = null;
		ArrayList<Cluster> clustersOfOneTranscript = null;
		int smallCoordinatePromoterRegion = -1;
		int transcriptTSSNumbering = -1;
		ArrayList<File> promoterRegions = geneAndPromoterRegions.get(geneName);
		ArrayList<String> transcripts = geneAndTranscripts.get(geneName);
		if(transcripts == null)
		{
			System.out.println(geneName);
			System.exit(10);
		}
		
		if(promoterRegions.size() == 1)
		{
			allClustersOfTheGene = restoreClustersFromOneFile(clusterFolder + "/" + promoterRegions.get(0).getName());
			smallCoordinatePromoterRegion = getSmallCoordinatePromoterRegion(promoterRegions.get(0));
			it = transcripts.iterator();			
			while(it.hasNext())
			{
				oneTranscript = it.next();
				fields = oneTranscript.split(",");
				transcriptTSSNumbering = Integer.parseInt(fields[3]);
				FeatureSetOfATSS aFeatureSet = ifTheTSSWasEncountered(transcriptTSSNumbering, aGene);
				if(aFeatureSet == null)
				{
					clustersOfOneTranscript = getClustersOfOneTranscript(fields[6], transcriptTSSNumbering, smallCoordinatePromoterRegion, allClustersOfTheGene, ifAbsoluteCoor);
					if(clustersOfOneTranscript.size() == 0)
					{
						continue;
					}
					FeatureDistanceBetweenClusterAndTSS feature1 = extractFeature1OfATranscript(fields[6], clustersOfOneTranscript, smallCoordinatePromoterRegion, Integer.parseInt(fields[3]), ifAbsoluteCoor);
					FeatureClusterLength feature2 = extractFeature2OfATranscript(fields[6], clustersOfOneTranscript);
					FeatureClusterIC feature3 = extractFeature3OfATranscript(fields[6], clustersOfOneTranscript);
					FeatureNumberSitesTFCluster feature4 = extractFeature4OfATranscript(fields[6], clustersOfOneTranscript, threshold, TFsAndIndices, modelToRseq);
					FeatureNumberStrongSitesTFCluster feature5 = extractFeature5OfATranscript(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
					FeatureSumICSitesTFCluster feature6 = extractFeature6OfATranscript(fields[6], clustersOfOneTranscript, threshold, TFsAndIndices, modelToRseq);
					FeatureSumICStrongSitesTFCluster feature7 = extractFeature7OfATranscript(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
					aGene.addATranscript(new FeatureSetOfATSS(transcriptTSSNumbering, fields[5], Integer.parseInt(fields[9]), feature1, feature2, feature3, feature4, feature5, feature6, feature7));
				}
				else
				{
					aFeatureSet.addATranscript(fields[5]);
					aFeatureSet.addALength(Integer.parseInt(fields[9]));
				}
			}
		}		
		else
		{
			File promoterRegion = null;
			it = transcripts.iterator();
			while(it.hasNext())
			{
				oneTranscript = it.next();
				fields = oneTranscript.split(",");
				transcriptTSSNumbering = Integer.parseInt(fields[3]);
				FeatureSetOfATSS aFeatureSet = ifTheTSSWasEncountered(transcriptTSSNumbering, aGene);
				if(aFeatureSet == null)
				{
					
					promoterRegion = determinePromoterRegion(transcriptTSSNumbering, promoterRegions);
					if(promoterRegion == null)
					{
						continue;
					}
					
					
					allClustersOfTheGene = restoreClustersFromOneFile(clusterFolder + "/" + promoterRegion.getName());
					if(allClustersOfTheGene == null)
					{
						continue;
					}
					smallCoordinatePromoterRegion = getSmallCoordinatePromoterRegion(promoterRegion);
					clustersOfOneTranscript = getClustersOfOneTranscript(fields[6], transcriptTSSNumbering, smallCoordinatePromoterRegion, allClustersOfTheGene, ifAbsoluteCoor);
					if(clustersOfOneTranscript.size() == 0)
					{
						continue;
					}
					FeatureDistanceBetweenClusterAndTSS feature1 = extractFeature1OfATranscript(fields[6], clustersOfOneTranscript, smallCoordinatePromoterRegion, Integer.parseInt(fields[3]), ifAbsoluteCoor);
					FeatureClusterLength feature2 = extractFeature2OfATranscript(fields[6], clustersOfOneTranscript);
					FeatureClusterIC feature3 = extractFeature3OfATranscript(fields[6], clustersOfOneTranscript);
					FeatureNumberSitesTFCluster feature4 = extractFeature4OfATranscript(fields[6], clustersOfOneTranscript, threshold, TFsAndIndices, modelToRseq);
					FeatureNumberStrongSitesTFCluster feature5 = extractFeature5OfATranscript(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
					FeatureSumICSitesTFCluster feature6 = extractFeature6OfATranscript(fields[6], clustersOfOneTranscript, threshold, TFsAndIndices, modelToRseq);
					FeatureSumICStrongSitesTFCluster feature7 = extractFeature7OfATranscript(fields[6], clustersOfOneTranscript, strongThreshold, TFsAndIndices, modelToRseq);
					aGene.addATranscript(new FeatureSetOfATSS(transcriptTSSNumbering, fields[5], Integer.parseInt(fields[9]), feature1, feature2, feature3, feature4, feature5, feature6, feature7));
				}
				else
				{
					aFeatureSet.addATranscript(fields[5]);
					aFeatureSet.addALength(Integer.parseInt(fields[9]));
				}
			}
		}
		
		return aGene;
	}
	
	private static TreeMap<String, Integer> getTFsAndIndices(String modelFolderPath)
	{
		File modelFolder = new File(modelFolderPath);
		File[] models = modelFolder.listFiles();
		TreeSet<String> TFs = new TreeSet<String>();
		for(int i = 0; i < models.length; i++)
		{
			TFs.add(models[i].getName().split("[.-]")[0]);
		}
		
		TreeMap<String, Integer> TFsAndIndices = new TreeMap<String, Integer>();
		int index = 0;
		Iterator<String> it = TFs.iterator();
		while(it.hasNext())
		{
			TFsAndIndices.put(it.next(), index++);
		}
		
		return TFsAndIndices;
	}
	
	private static FeatureOrderOfStrongSitesCluster extractFeature9OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript, double strongThreshold, HashMap<String, Double> modelToRseq)
	{
		ArrayList<ArrayList<TreeSet<String>>> orderOfStrongSites = new ArrayList<ArrayList<TreeSet<String>>>();
		ArrayList<TreeSet<String>> orderOfStrongSitesInOneCluster = null;
		Cluster aCluster = null;
		TreeSet<BindingSite> allSitesInOneCluster = null;
		Iterator<BindingSite> it = null;
		BindingSite oneSite = null;
		TreeSet<String> TFs = null;
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				orderOfStrongSitesInOneCluster = new ArrayList<TreeSet<String>>();
				aCluster = clustersOfOneTranscript.get(i);
				allSitesInOneCluster = new TreeSet<BindingSite>(new BindingSiteComparator());
				
				allSitesInOneCluster.add(aCluster.getCenterSite());
				allSitesInOneCluster.addAll(aCluster.getOtherSites());
				it = allSitesInOneCluster.descendingIterator();
				while(it.hasNext())
				{
					oneSite = it.next();
					TFs = oneSite.getTFs(modelToRseq, strongThreshold);
					if(!TFs.isEmpty())
					{
						orderOfStrongSitesInOneCluster.add(TFs);
					}					
				}
				
				orderOfStrongSites.add(orderOfStrongSitesInOneCluster);
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				orderOfStrongSitesInOneCluster = new ArrayList<TreeSet<String>>();
				aCluster = clustersOfOneTranscript.get(i);
				allSitesInOneCluster = new TreeSet<BindingSite>(new BindingSiteComparator());
				
				allSitesInOneCluster.add(aCluster.getCenterSite());
				allSitesInOneCluster.addAll(aCluster.getOtherSites());
				it = allSitesInOneCluster.iterator();
				while(it.hasNext())
				{
					oneSite = it.next();
					TFs = oneSite.getTFs(modelToRseq, strongThreshold);
					if(!TFs.isEmpty())
					{
						orderOfStrongSitesInOneCluster.add(TFs);
					}
				}
				
				orderOfStrongSites.add(orderOfStrongSitesInOneCluster);
			}
		}
		
		return new FeatureOrderOfStrongSitesCluster(orderOfStrongSites);
	}
	
	private static FeatureOrderOfSitesCluster extractFeature8OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript, double threshold, HashMap<String, Double> modelToRseq)
	{
		ArrayList<ArrayList<TreeSet<String>>> orderOfSites = new ArrayList<ArrayList<TreeSet<String>>>();
		ArrayList<TreeSet<String>> orderOfSitesInOneCluster = null;
		Cluster aCluster = null;
		TreeSet<BindingSite> allSitesInOneCluster = null;
		Iterator<BindingSite> it = null;
		BindingSite oneSite = null;
		TreeSet<String> TFs = null;
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				orderOfSitesInOneCluster = new ArrayList<TreeSet<String>>();
				aCluster = clustersOfOneTranscript.get(i);
				allSitesInOneCluster = new TreeSet<BindingSite>(new BindingSiteComparator());
				
				allSitesInOneCluster.add(aCluster.getCenterSite());
				allSitesInOneCluster.addAll(aCluster.getOtherSites());
				it = allSitesInOneCluster.descendingIterator();
				while(it.hasNext())
				{
					oneSite = it.next();
					TFs = oneSite.getTFs(modelToRseq, threshold);
					if(!TFs.isEmpty())
					{
						orderOfSitesInOneCluster.add(TFs);
					}					
				}
				
				orderOfSites.add(orderOfSitesInOneCluster);
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				orderOfSitesInOneCluster = new ArrayList<TreeSet<String>>();
				aCluster = clustersOfOneTranscript.get(i);
				allSitesInOneCluster = new TreeSet<BindingSite>(new BindingSiteComparator());
				
				allSitesInOneCluster.add(aCluster.getCenterSite());
				allSitesInOneCluster.addAll(aCluster.getOtherSites());
				it = allSitesInOneCluster.iterator();
				while(it.hasNext())
				{
					oneSite = it.next();
					TFs = oneSite.getTFs(modelToRseq, threshold);
					if(!TFs.isEmpty())
					{
						orderOfSitesInOneCluster.add(TFs);
					}
				}
				
				orderOfSites.add(orderOfSitesInOneCluster);
			}
		}
		
		return new FeatureOrderOfSitesCluster(orderOfSites);
	}
	
	private static FeatureSumICStrongSitesTFCluster extractFeature7OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript, double strongThreshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq)
	{
		ArrayList<double[]> sumsOfICOfStrongSites = new ArrayList<double[]>();
		double[] sumsOfICOfStrongSitesInOneCluster = null;
		Cluster aCluster = null;
		HashMap<String, Double> TFsAndICs = null;
		Set<Entry<String, Double>> entries = null;
		Iterator<Entry<String, Double>> it = null;
		Entry<String, Double> entry = null;
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				sumsOfICOfStrongSitesInOneCluster = new double[TFsAndIndices.size()];
				for(int j = 0; j < sumsOfICOfStrongSitesInOneCluster.length; j++)
				{
					sumsOfICOfStrongSitesInOneCluster[j] = 0.0;
				}
				
				aCluster = clustersOfOneTranscript.get(i);
				TFsAndICs = aCluster.getCenterSite().getTFsAndICs(modelToRseq, strongThreshold);
				entries = TFsAndICs.entrySet();
				it = entries.iterator();
				while(it.hasNext())
				{
					entry = it.next();
					sumsOfICOfStrongSitesInOneCluster[TFsAndIndices.get(entry.getKey())] += entry.getValue();
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFsAndICs = aSite.getTFsAndICs(modelToRseq, strongThreshold);
					entries = TFsAndICs.entrySet();
					it = entries.iterator();
					while(it.hasNext())
					{
						entry = it.next();
						sumsOfICOfStrongSitesInOneCluster[TFsAndIndices.get(entry.getKey())] += entry.getValue();
					}
				}
				
				sumsOfICOfStrongSites.add(sumsOfICOfStrongSitesInOneCluster);
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				sumsOfICOfStrongSitesInOneCluster = new double[TFsAndIndices.size()];
				for(int j = 0; j < sumsOfICOfStrongSitesInOneCluster.length; j++)
				{
					sumsOfICOfStrongSitesInOneCluster[j] = 0.0;
				}
				
				aCluster = clustersOfOneTranscript.get(i);
				TFsAndICs = aCluster.getCenterSite().getTFsAndICs(modelToRseq, strongThreshold);
				entries = TFsAndICs.entrySet();
				it = entries.iterator();
				while(it.hasNext())
				{
					entry = it.next();
					sumsOfICOfStrongSitesInOneCluster[TFsAndIndices.get(entry.getKey())] += entry.getValue();
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFsAndICs = aSite.getTFsAndICs(modelToRseq, strongThreshold);
					entries = TFsAndICs.entrySet();
					it = entries.iterator();
					while(it.hasNext())
					{
						entry = it.next();
						sumsOfICOfStrongSitesInOneCluster[TFsAndIndices.get(entry.getKey())] += entry.getValue();
					}
				}
				
				sumsOfICOfStrongSites.add(sumsOfICOfStrongSitesInOneCluster);
			}
		}
		
		return new FeatureSumICStrongSitesTFCluster(sumsOfICOfStrongSites);
	}
	
	private static FeatureSumICSitesTFCluster extractFeature6OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript, double threshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq)
	{
		ArrayList<double[]> sumsOfICOfSites = new ArrayList<double[]>();
		double[] sumsOfICOfSitesInOneCluster = null;
		Cluster aCluster = null;
		HashMap<String, Double> TFsAndICs = null;
		Set<Entry<String, Double>> entries = null;
		Iterator<Entry<String, Double>> it = null;
		Entry<String, Double> entry = null;
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				sumsOfICOfSitesInOneCluster = new double[TFsAndIndices.size()];
				for(int j = 0; j < sumsOfICOfSitesInOneCluster.length; j++)
				{
					sumsOfICOfSitesInOneCluster[j] = 0.0;
				}
				
				aCluster = clustersOfOneTranscript.get(i);
				TFsAndICs = aCluster.getCenterSite().getTFsAndICs(modelToRseq, threshold);
				entries = TFsAndICs.entrySet();
				it = entries.iterator();
				while(it.hasNext())
				{
					entry = it.next();
					sumsOfICOfSitesInOneCluster[TFsAndIndices.get(entry.getKey())] += entry.getValue();
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFsAndICs = aSite.getTFsAndICs(modelToRseq, threshold);
					entries = TFsAndICs.entrySet();
					it = entries.iterator();
					while(it.hasNext())
					{
						entry = it.next();
						sumsOfICOfSitesInOneCluster[TFsAndIndices.get(entry.getKey())] += entry.getValue();
					}
				}
				
				sumsOfICOfSites.add(sumsOfICOfSitesInOneCluster);
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				sumsOfICOfSitesInOneCluster = new double[TFsAndIndices.size()];
				for(int j = 0; j < sumsOfICOfSitesInOneCluster.length; j++)
				{
					sumsOfICOfSitesInOneCluster[j] = 0.0;
				}
				
				aCluster = clustersOfOneTranscript.get(i);
				TFsAndICs = aCluster.getCenterSite().getTFsAndICs(modelToRseq, threshold);
				entries = TFsAndICs.entrySet();
				it = entries.iterator();
				while(it.hasNext())
				{
					entry = it.next();
					sumsOfICOfSitesInOneCluster[TFsAndIndices.get(entry.getKey())] += entry.getValue();
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFsAndICs = aSite.getTFsAndICs(modelToRseq, threshold);
					entries = TFsAndICs.entrySet();
					it = entries.iterator();
					while(it.hasNext())
					{
						entry = it.next();
						sumsOfICOfSitesInOneCluster[TFsAndIndices.get(entry.getKey())] += entry.getValue();
					}
				}
				
				sumsOfICOfSites.add(sumsOfICOfSitesInOneCluster);
			}
		}
		
		return new FeatureSumICSitesTFCluster(sumsOfICOfSites);
	}
	
	private static FeatureNumberStrongSitesTFCluster extractFeature5OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript, double strongThreshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq)
	{
		ArrayList<int[]> numbersOfStrongSites = new ArrayList<int[]>();
		int[] numberOfStrongSitesInOneCluster = null;
		Cluster aCluster = null;
		TreeSet<String> TFs = null;
		Iterator<String> it = null;
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				numberOfStrongSitesInOneCluster = new int[TFsAndIndices.size()];
				for(int j = 0; j < numberOfStrongSitesInOneCluster.length; j++)
				{
					numberOfStrongSitesInOneCluster[j] = 0;
				}
				
				aCluster = clustersOfOneTranscript.get(i);
				TFs = aCluster.getCenterSite().getTFs(modelToRseq, strongThreshold);
				it = TFs.iterator();
				while(it.hasNext())
				{
					numberOfStrongSitesInOneCluster[TFsAndIndices.get(it.next())]++;
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFs = aSite.getTFs(modelToRseq, strongThreshold);
					it = TFs.iterator();
					while(it.hasNext())
					{
						numberOfStrongSitesInOneCluster[TFsAndIndices.get(it.next())]++;
					}
				}
				
				numbersOfStrongSites.add(numberOfStrongSitesInOneCluster);
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				numberOfStrongSitesInOneCluster = new int[TFsAndIndices.size()];
				for(int j = 0; j < numberOfStrongSitesInOneCluster.length; j++)
				{
					numberOfStrongSitesInOneCluster[j] = 0;
				}
				
				aCluster = clustersOfOneTranscript.get(i);
				TFs = aCluster.getCenterSite().getTFs(modelToRseq, strongThreshold);
				it = TFs.iterator();
				while(it.hasNext())
				{
					numberOfStrongSitesInOneCluster[TFsAndIndices.get(it.next())]++;
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFs = aSite.getTFs(modelToRseq, strongThreshold);
					it = TFs.iterator();
					while(it.hasNext())
					{
						numberOfStrongSitesInOneCluster[TFsAndIndices.get(it.next())]++;
					}
				}
				
				numbersOfStrongSites.add(numberOfStrongSitesInOneCluster);
			}
		}		
		
		return new FeatureNumberStrongSitesTFCluster(numbersOfStrongSites);
	}
	
	private static FeatureNumberSitesTFCluster extractFeature4OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript, double threshold, TreeMap<String, Integer> TFsAndIndices, HashMap<String, Double> modelToRseq)
	{
		ArrayList<int[]> numbersOfSites = new ArrayList<int[]>();
		int[] numberOfSitesInOneCluster = null;
		Cluster aCluster = null;
		TreeSet<String> TFs = null;
		Iterator<String> it = null;
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				numberOfSitesInOneCluster = new int[TFsAndIndices.size()];
				for(int j = 0; j < numberOfSitesInOneCluster.length; j++)
				{
					numberOfSitesInOneCluster[j] = 0;
				}
				
				aCluster = clustersOfOneTranscript.get(i);
				TFs = aCluster.getCenterSite().getTFs(modelToRseq, threshold);
				it = TFs.iterator();
				while(it.hasNext())
				{
					numberOfSitesInOneCluster[TFsAndIndices.get(it.next())]++;
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFs = aSite.getTFs(modelToRseq, threshold);
					it = TFs.iterator();
					while(it.hasNext())
					{
						numberOfSitesInOneCluster[TFsAndIndices.get(it.next())]++;
					}
				}
				
				numbersOfSites.add(numberOfSitesInOneCluster);
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				numberOfSitesInOneCluster = new int[TFsAndIndices.size()];
				for(int j = 0; j < numberOfSitesInOneCluster.length; j++)
				{
					numberOfSitesInOneCluster[j] = 0;
				}
				
				aCluster = clustersOfOneTranscript.get(i);
				TFs = aCluster.getCenterSite().getTFs(modelToRseq, threshold);
				it = TFs.iterator();
				while(it.hasNext())
				{
					numberOfSitesInOneCluster[TFsAndIndices.get(it.next())]++;
				}
				
				TreeSet<BindingSite> otherSites = aCluster.getOtherSites();
				for(BindingSite aSite : otherSites)
				{
					TFs = aSite.getTFs(modelToRseq, threshold);
					it = TFs.iterator();
					while(it.hasNext())
					{
						numberOfSitesInOneCluster[TFsAndIndices.get(it.next())]++;
					}
				}
				
				numbersOfSites.add(numberOfSitesInOneCluster);
			}
		}		
		
		return new FeatureNumberSitesTFCluster(numbersOfSites);
	}
	
	private static FeatureClusterIC extractFeature3OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript)
	{
		ArrayList<Double> ICs = new ArrayList<Double>();
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				ICs.add(clustersOfOneTranscript.get(i).getIC());
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				ICs.add(clustersOfOneTranscript.get(i).getIC());
			}
		}		
		
		return new FeatureClusterIC(ICs);
	}
	
	private static FeatureClusterLength extractFeature2OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript)
	{
		ArrayList<Integer> lengths = new ArrayList<Integer>();
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				lengths.add(clustersOfOneTranscript.get(i).getLength());
			}
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				lengths.add(clustersOfOneTranscript.get(i).getLength());
			}
		}		
				
		return new FeatureClusterLength(lengths);
	}
	
	private static FeatureDistanceBetweenClusterAndTSS extractFeature1OfATranscript(String strand, ArrayList<Cluster> clustersOfOneTranscript, int smallCoordinatePromoterRegion, int transcriptTSSNumbering, boolean ifAbsoluteCoor)
	{
		int absoluteIndexTSS = transcriptTSSNumbering - 1;
		ArrayList<Integer> distances = new ArrayList<Integer>();
		
		if(strand.equals("1"))
		{
			for(int i = clustersOfOneTranscript.size() - 1; i >=0 ; i--)
			{
				int absoluteDistance = absoluteIndexTSS - clustersOfOneTranscript.get(i).getLargeCoordinate();
				if(!ifAbsoluteCoor)
				{
					absoluteDistance -= smallCoordinatePromoterRegion;
				}
				
				if(absoluteDistance < 0)
				{
					distances.add(0);
				}
				else
				{
					distances.add(absoluteDistance);
				}
			}
			
		}
		else
		{
			for(int i = 0; i < clustersOfOneTranscript.size(); i++)
			{
				int absoluteDistance = clustersOfOneTranscript.get(i).getSmallCoordinate() - absoluteIndexTSS;
				if(!ifAbsoluteCoor)
				{
					absoluteDistance += smallCoordinatePromoterRegion;
				}
				if(absoluteDistance < 0)
				{
					distances.add(0);
				}
				else
				{
					distances.add(absoluteDistance);
				}
			}
		}
		
		return new FeatureDistanceBetweenClusterAndTSS(distances);
	}
	
	private static int getSmallCoordinatePromoterRegion(File promoterRegion)
	{
		int smallCoordinatePromoterRegion = -1;
		
		try
		{
			Scanner input = new Scanner(promoterRegion);
			String[] fieldsInLineInBedFile = input.nextLine().split("\t");
			smallCoordinatePromoterRegion = Integer.parseInt(fieldsInLineInBedFile[1]);			
			input.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open the file: " + promoterRegion.getAbsolutePath());
			System.exit(1);
		}
		
		return smallCoordinatePromoterRegion;
	}
	
	private static ArrayList<Cluster> getClustersOfOneTranscript(String strand, int transcriptTSSNumbering, int smallCoordinatePromoterRegion, ArrayList<Cluster> allClustersOfTheRegion, boolean ifAbsoluteCoor)
	{
		ArrayList<Cluster> clustersOfOneTranscript = new ArrayList<Cluster>();		
		
		if(strand.equals("1"))
		{
			int smallAbsoluteIndex = transcriptTSSNumbering - 10001;
			int largeAbsoluteIndex = smallAbsoluteIndex + 9999;
			for(int i = 0; i < allClustersOfTheRegion.size(); i++)
			{
				Cluster aCluster = allClustersOfTheRegion.get(i);
				if(aCluster.intersect(smallAbsoluteIndex, largeAbsoluteIndex, ifAbsoluteCoor, smallCoordinatePromoterRegion))
				{
					clustersOfOneTranscript.add(aCluster);
				}
			}
		}
		else
		{
			int smallAbsoluteIndex = transcriptTSSNumbering;
			int largeAbsoluteIndex = smallAbsoluteIndex + 9999;
			for(int i = 0; i < allClustersOfTheRegion.size(); i++)
			{
				Cluster aCluster = allClustersOfTheRegion.get(i);
				if(aCluster.intersect(smallAbsoluteIndex, largeAbsoluteIndex, ifAbsoluteCoor, smallCoordinatePromoterRegion))
				{
					clustersOfOneTranscript.add(aCluster);
				}
			}
		}
		
		return clustersOfOneTranscript;
	}
	
	private static ArrayList<Cluster> restoreClustersFromOneFile(String clusterFile) 
	{
		ArrayList<Cluster> clusters = new ArrayList<Cluster>();
		try
		{
			ArrayList<String> linesInOneCluster = new ArrayList<String>();
			boolean clusterStart = false;
			
			Scanner input = new Scanner(new File(clusterFile));
			while(input.hasNextLine())
			{
				String aLine = input.nextLine();
				if(!clusterStart)
				{
					if(aLine.contains("********************************************************"))
					{
						clusterStart = true;
						continue;
					}
				}
				if(clusterStart)
				{
					if(aLine.contains("********************************************************"))
					{
						clusterStart = false;			
						clusters.add(restoreOneCluster(linesInOneCluster));
						linesInOneCluster.clear();
					}
					else
					{
						linesInOneCluster.add(aLine);
					}
				}				
			}
			input.close();
		}
		catch(FileNotFoundException e)
		{			
			return null;
		}
		
		return clusters;
	}
	
	private static Cluster restoreOneCluster(ArrayList<String> linesInOneCluster)
	{
		TreeSet<BindingSite> otherSites = new TreeSet<BindingSite>(new BindingSiteComparator());
		BindingSite centerSite = null;
		
		String[] fieldsFirstLine = linesInOneCluster.get(0).split("[:.\\-]");
		for(int i = 3; i < linesInOneCluster.size(); i++)
		{
			HashMap<String, Double> otherTFandRi = new HashMap<String, Double>();
			String[] fieldsOneSite = linesInOneCluster.get(i).split("\t");
			if(fieldsOneSite.length == 6)
			{
				String[] fieldsOtherTFs = fieldsOneSite[5].split("[{,}]");
				for(int j = 1; j < fieldsOtherTFs.length; j++)
				{
					String[] fieldsOneOtherTF = fieldsOtherTFs[j].trim().split("=");
					otherTFandRi.put(fieldsOneOtherTF[0], Double.valueOf(fieldsOneOtherTF[1]));
				}
				otherTFandRi.remove(fieldsOneSite[3]);
			}			
			
			if(fieldsOneSite[1].equals("+"))
			{
				if(i == 3)
				{
					centerSite = new BindingSite(Integer.parseInt(fieldsOneSite[0]), true, fieldsOneSite[2], fieldsOneSite[3], Double.parseDouble(fieldsOneSite[4]), otherTFandRi);
				}
				else
				{
					otherSites.add(new BindingSite(Integer.parseInt(fieldsOneSite[0]), true, fieldsOneSite[2], fieldsOneSite[3], Double.parseDouble(fieldsOneSite[4]), otherTFandRi));
				}
			}
			else
			{
				if(i == 3)
				{
					centerSite = new BindingSite(Integer.parseInt(fieldsOneSite[0]), false, fieldsOneSite[2], fieldsOneSite[3], Double.parseDouble(fieldsOneSite[4]), otherTFandRi);
				}
				else
				{
					otherSites.add(new BindingSite(Integer.parseInt(fieldsOneSite[0]), false, fieldsOneSite[2], fieldsOneSite[3], Double.parseDouble(fieldsOneSite[4]), otherTFandRi));
				}
			}
		}
		
		return new Cluster(centerSite, otherSites, Double.parseDouble(linesInOneCluster.get(2).split(":")[1].trim()), Integer.parseInt(fieldsFirstLine[1].trim()), Integer.parseInt(fieldsFirstLine[2].trim()), Integer.parseInt(fieldsFirstLine[4].trim()), Integer.parseInt(linesInOneCluster.get(1).split(":")[1].trim()));
	}
	
	private static TreeMap<String, ArrayList<File>> getPromoterRegions(String dir)
	{
		TreeMap<String, ArrayList<File>> geneAndPromoterRegions = new TreeMap<String, ArrayList<File>>();
		File promoterRegionFolder = new File(dir);
		File[] promoterRegionFiles = promoterRegionFolder.listFiles();
		String geneName = null;
		
		for(int i = 0; i < promoterRegionFiles.length; i++)
		{
			String aFilename = promoterRegionFiles[i].getName();
			if(aFilename.contains("__"))
			{
				String[] fields = aFilename.split("__");
				if(fields[1].contains("_"))
				{
					String[] subFields = fields[1].split("_");
					geneName = fields[0] + "__" + subFields[0];
				}
				else
				{
					geneName = aFilename.substring(0, aFilename.length() - 4);
				}
			}
			else
			{
				if(aFilename.contains("_"))
				{
					geneName = aFilename.split("_")[0];
				}
				else
				{
					geneName = aFilename.substring(0, aFilename.length() - 4);
				}
			}			
			
			if(!geneAndPromoterRegions.containsKey(geneName))
			{
				ArrayList<File> promoterRegions = new ArrayList<File>();
				promoterRegions.add(promoterRegionFiles[i]);
				geneAndPromoterRegions.put(geneName, promoterRegions);
			}
			else
			{
				ArrayList<File> promoterRegions = geneAndPromoterRegions.get(geneName);
				promoterRegions.add(promoterRegionFiles[i]);
			}
		}
		
		return geneAndPromoterRegions;
	}
	
	private static TreeMap<String, ArrayList<String>> getTranscriptsInfo(String filePath)
	{
		TreeMap<String, ArrayList<String>> geneAndTranscripts = new TreeMap<String, ArrayList<String>>();
		
		try
		{
			Scanner input = new Scanner(new File(filePath));
			input.nextLine();		
			while(input.hasNextLine())
			{
				String aLine = input.nextLine();
				String[] fields = aLine.split(",");
				if(!geneAndTranscripts.containsKey(fields[4]))
				{
					ArrayList<String> transcripts = new ArrayList<String>();
					transcripts.add(aLine);
					geneAndTranscripts.put(fields[4], transcripts);
				}
				else
				{
					ArrayList<String> transcripts = geneAndTranscripts.get(fields[4]);
					transcripts.add(aLine);
				}
			}
			input.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open the file: " + filePath);
			System.exit(1);
		}
		
		return geneAndTranscripts;
	}
	
	private static TreeSet<String> getNamesOfGenesForWhichClustersWereObtained(String folderPath)
	{
		File clusterFolder = new File(folderPath);
		File[] clusterFiles = clusterFolder.listFiles();
		
		TreeSet<String> geneNames = new TreeSet<String>();
		for(int i = 0; i < clusterFiles.length; i++)
		{
			String aFilename = clusterFiles[i].getName();
			if(aFilename.contains("__"))
			{
				String[] fields = aFilename.split("__");
				if(fields[1].contains("_"))
				{
					String[] subFields = fields[1].split("_");
					geneNames.add(fields[0] + "__" + subFields[0]);
				}
				else
				{
					geneNames.add(aFilename.substring(0, aFilename.length() - 4));
				}
			}
			else
			{
				if(aFilename.contains("_"))
				{
					geneNames.add(aFilename.split("_")[0]);
				}
				else
				{
					geneNames.add(aFilename.substring(0, aFilename.length() - 4));
				}
			}
		}
		
		return geneNames;
	}
}
