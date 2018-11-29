import java.util.ArrayList;

public class FeatureSetOfATSS_oneiPWM
{
	private int TSS;
	private ArrayList<String> transcripts;
	private ArrayList<Integer> transcriptLengths;
	private FeatureDistanceBetweenClusterAndTSS feature1;
	private FeatureClusterLength feature2;
	private FeatureClusterIC feature3;
	private FeatureNumberSitesTFCluster_oneiPWM feature4;
	private FeatureNumberStrongSitesTFCluster_oneiPWM feature5;
	private FeatureSumICStrongSitesTFCluster_oneiPWM feature6;
	
	public FeatureSetOfATSS_oneiPWM(int aTSS, String aTranscript, int aLength)
	{
		TSS = aTSS;
		transcripts = new ArrayList<String>();
		transcriptLengths = new ArrayList<Integer>();
		transcripts.add(aTranscript);
		transcriptLengths.add(aLength);
		feature1 = new FeatureDistanceBetweenClusterAndTSS();
		feature2 = new FeatureClusterLength();
		feature3 = new FeatureClusterIC();
		feature4 = new FeatureNumberSitesTFCluster_oneiPWM();
		feature5 = new FeatureNumberStrongSitesTFCluster_oneiPWM();
		feature6 = new FeatureSumICStrongSitesTFCluster_oneiPWM();
	}
	
	public FeatureSetOfATSS_oneiPWM(int aTSS, String aTranscript, int aLength, FeatureDistanceBetweenClusterAndTSS aFeature1, FeatureClusterLength aFeature2, FeatureClusterIC aFeature3, FeatureNumberSitesTFCluster_oneiPWM aFeature4, FeatureNumberStrongSitesTFCluster_oneiPWM aFeature5, FeatureSumICStrongSitesTFCluster_oneiPWM aFeature6)
	{
		TSS = aTSS;
		transcripts = new ArrayList<String>();
		transcriptLengths = new ArrayList<Integer>();
		transcripts.add(aTranscript);
		transcriptLengths.add(aLength);
		feature1 = aFeature1;
		feature2 = aFeature2;
		feature3 = aFeature3;
		feature4 = aFeature4;
		feature5 = aFeature5;
		feature6 = aFeature6;
	}
	
	public void addATranscript(String aTranscript)
	{
		transcripts.add(aTranscript);
	}
	
	public void addALength(int aLength)
	{
		transcriptLengths.add(aLength);
	}
	
	public int getTSS()
	{
		return TSS;
	}
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		sb.append("***************************************************************************\n");
		
		String allTranscriptNames = "";
		for(int i = 0; i < transcripts.size() - 1; i++)
		{
			allTranscriptNames += transcripts.get(i) + ",";
		}
		allTranscriptNames += transcripts.get(transcripts.size() - 1);
		
		String allTranscriptLengths = "";
		for(int i = 0; i < transcriptLengths.size() - 1; i++)
		{
			allTranscriptLengths += transcriptLengths.get(i) + ",";
		}
		allTranscriptLengths += transcriptLengths.get(transcriptLengths.size() - 1);
		
		sb.append("Transcripts " + allTranscriptNames + ":\n");
		sb.append("Transcript lengths " + allTranscriptLengths + ":\n");
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
