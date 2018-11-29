import java.util.ArrayList;

public class FeatureSetOfATSS_oneiPWM
{
	private ArrayList<String> transcripts;
	private ArrayList<Integer> lengths;
	private FeatureDistanceBetweenClusterAndTSS feature1;
	private FeatureClusterLength feature2;
	private FeatureClusterIC feature3;
	private FeatureNumberSitesTFCluster_oneiPWM feature4;
	private FeatureNumberStrongSitesTFCluster_oneiPWM feature5;
	private FeatureSumICStrongSitesTFCluster_oneiPWM feature6;
	
	public FeatureSetOfATSS_oneiPWM(String someTranscripts, String someLengths)
	{
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
		feature4 = new FeatureNumberSitesTFCluster_oneiPWM();
		feature5 = new FeatureNumberStrongSitesTFCluster_oneiPWM();
		feature6 = new FeatureSumICStrongSitesTFCluster_oneiPWM();
	}
	
	public FeatureSetOfATSS_oneiPWM(String someTranscripts, String someLengths, FeatureDistanceBetweenClusterAndTSS aFeature1, FeatureClusterLength aFeature2, FeatureClusterIC aFeature3, FeatureNumberSitesTFCluster_oneiPWM aFeature4, FeatureNumberStrongSitesTFCluster_oneiPWM aFeature5, FeatureSumICStrongSitesTFCluster_oneiPWM aFeature6)
	{
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
		
		feature1 = aFeature1;
		feature2 = aFeature2;
		feature3 = aFeature3;
		feature4 = aFeature4;
		feature5 = aFeature5;
		feature6 = aFeature6;
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
	
	public FeatureNumberSitesTFCluster_oneiPWM getFeature4()
	{
		return feature4;
	}
	
	public FeatureNumberStrongSitesTFCluster_oneiPWM getFeature5()
	{
		return feature5;
	}
		
	public FeatureSumICStrongSitesTFCluster_oneiPWM getFeature6()
	{
		return feature6;
	}
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		sb.append("***************************************************************************\n");		
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
		sb.append("***************************************************************************\n");
		
		return sb.toString();
	}
}

