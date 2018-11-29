import java.util.ArrayList;

public class FeatureSetOfATSS 
{
	//private int TSS;
	//private String transcripts;
	private ArrayList<String> transcripts;
	private ArrayList<Integer> lengths;
	private FeatureDistanceBetweenClusterAndTSS feature1;
	private FeatureClusterLength feature2;
	private FeatureClusterIC feature3;
	private FeatureNumberSitesTFCluster feature4;
	private FeatureNumberStrongSitesTFCluster feature5;
	private FeatureSumICSitesTFCluster feature6;
	private FeatureSumICStrongSitesTFCluster feature7;
	/*private FeatureOrderOfSitesCluster feature8;
	private FeatureOrderOfStrongSitesCluster feature9;*/
	
	public FeatureSetOfATSS(/*int aTSS, */String someTranscripts, String someLengths)
	{
		//TSS = aTSS;
		transcripts = new ArrayList<String>();
		String[] transcriptNames = someTranscripts.split(",");
		for(int i = 0; i < transcriptNames.length; i++)
		{
			transcripts.add(transcriptNames[i].trim());
		}
		
		lengths = new ArrayList<Integer>();
		String[] transcriptLengths = someTranscripts.split(",");
		for(int i = 0; i < transcriptLengths.length; i++)
		{
			lengths.add(Integer.parseInt(transcriptLengths[i].trim()));
		}
		
		feature1 = new FeatureDistanceBetweenClusterAndTSS();
		feature2 = new FeatureClusterLength();
		feature3 = new FeatureClusterIC();
		feature4 = new FeatureNumberSitesTFCluster();
		feature5 = new FeatureNumberStrongSitesTFCluster();
		feature6 = new FeatureSumICSitesTFCluster();
		feature7 = new FeatureSumICStrongSitesTFCluster();
		/*feature8 = new FeatureOrderOfSitesCluster();
		feature9 = new FeatureOrderOfStrongSitesCluster();*/
	}
	
	public FeatureSetOfATSS(/*int aTSS, */String someTranscripts, String someLengths, FeatureDistanceBetweenClusterAndTSS aFeature1, FeatureClusterLength aFeature2, FeatureClusterIC aFeature3, FeatureNumberSitesTFCluster aFeature4, FeatureNumberStrongSitesTFCluster aFeature5, FeatureSumICSitesTFCluster aFeature6, FeatureSumICStrongSitesTFCluster aFeature7/*, FeatureOrderOfSitesCluster aFeature8, FeatureOrderOfStrongSitesCluster aFeature9*/)
	{
		//TSS = aTSS;
		transcripts = new ArrayList<String>();
		String[] transcriptNames = someTranscripts.split(",");
		for(int i = 0; i < transcriptNames.length; i++)
		{
			transcripts.add(transcriptNames[i].trim());
		}
		
		lengths = new ArrayList<Integer>();
		String[] transcriptLengths = someLengths.split(",");
		for(int i = 0; i < transcriptLengths.length; i++)
		{
			lengths.add(Integer.parseInt(transcriptLengths[i].trim()));
		}
		/*transcripts.add(aTranscript);*/
		feature1 = aFeature1;
		feature2 = aFeature2;
		feature3 = aFeature3;
		feature4 = aFeature4;
		feature5 = aFeature5;
		feature6 = aFeature6;
		feature7 = aFeature7;
		/*feature8 = aFeature8;
		feature9 = aFeature9;*/
	}
	
	public ArrayList<String> getTranscripts()
	{
		return transcripts;
	}
	
	public ArrayList<Integer> getLengths()
	{
		return lengths;
	}
	
	public FeatureDistanceBetweenClusterAndTSS getFeature1()
	{
		return feature1;
	}
	
	public FeatureClusterLength getFeature2()
	{
		return feature2;
	}
	
	public FeatureClusterIC getFeature3()
	{
		return feature3;
	}
	
	public FeatureNumberSitesTFCluster getFeature4()
	{
		return feature4;
	}
	
	public FeatureNumberStrongSitesTFCluster getFeature5()
	{
		return feature5;
	}
	
	public FeatureSumICSitesTFCluster getFeature6()
	{
		return feature6;
	}
	
	public FeatureSumICStrongSitesTFCluster getFeature7()
	{
		return feature7;
	}
	
	/*public void addATranscript(String aTranscript)
	{
		transcripts.add(aTranscript);
	}*/
	
	/*public int getTSS()
	{
		return TSS;
	}*/
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		sb.append("***************************************************************************\n");
		/*String allTranscriptNames = "";
		for(int i = 0; i < transcripts.size() - 1; i++)
		{
			allTranscriptNames += transcripts.get(i) + ",";
		}
		allTranscriptNames += transcripts.get(transcripts.size() - 1);*/
		sb.append("Transcripts " + transcripts + ":\n");
		sb.append("Transcript length " + lengths + ":\n");
		sb.append("Feature 1:\n");
		sb.append(feature1.toString() + "\n");
		sb.append("Feature 2:\n");
		sb.append(feature2.toString() + "\n");
		sb.append("Feature 3:\n");
		sb.append(feature3.toString() + "\n");
		sb.append("Feature 4:\n");
		sb.append(feature4.toString() + "\n");
		sb.append("Feature 5:\n");
		sb.append(feature5.toString() + "\n");
		sb.append("Feature 6:\n");
		sb.append(feature6.toString() + "\n");
		sb.append("Feature 7:\n");
		sb.append(feature7.toString() + "\n");
		/*sb.append("Feature 8:\n");
		sb.append(feature8.toString() + "\n");
		sb.append("Feature 9:\n");
		sb.append(feature9.toString() + "\n");*/
		sb.append("***************************************************************************\n");
		
		return sb.toString();
	}
}
