import java.util.ArrayList;

public class Gene_oneiPWM
{
	private String name;
	private ArrayList<FeatureSetOfATSS_oneiPWM> featureSets;
	
	public Gene_oneiPWM(String aName)
	{
		name = aName;
		featureSets = new ArrayList<FeatureSetOfATSS_oneiPWM>();
	}
	
	public Gene_oneiPWM(String aName, ArrayList<FeatureSetOfATSS_oneiPWM> someFeatureSets)
	{
		name = aName;
		featureSets = someFeatureSets;
	}
	
	public String getName()
	{
		return name;
	}
	
	public void addATranscript(FeatureSetOfATSS_oneiPWM aTranscript)
	{
		featureSets.add(aTranscript);
	}
	
	public ArrayList<FeatureSetOfATSS_oneiPWM> getFeatureSets()
	{
		return featureSets;
	}
	
	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder();
		sb.append(name + "\n");
		int numberOfTranscripts = featureSets.size();
		for(int i = 0; i < numberOfTranscripts; i++)
		{
			sb.append(featureSets.get(i).toString());
			if(i != numberOfTranscripts - 1)
			{
				sb.append("\n");
			}
		}
		
		return sb.toString();
	}
}
