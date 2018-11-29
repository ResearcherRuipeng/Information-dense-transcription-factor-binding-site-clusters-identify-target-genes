import java.util.ArrayList;

public class Gene 
{
	private String name;
	private ArrayList<FeatureSetOfATSS> featureSets;
	
	public Gene(String aName)
	{
		name = aName;
		featureSets = new ArrayList<FeatureSetOfATSS>();
	}
	
	public Gene(String aName, ArrayList<FeatureSetOfATSS> someFeatureSets)
	{
		name = aName;
		featureSets = someFeatureSets;
	}
	
	public String getName()
	{
		return name;
	}
	
	public void addATranscript(FeatureSetOfATSS aTranscript)
	{
		featureSets.add(aTranscript);
	}
	
	public ArrayList<FeatureSetOfATSS> getFeatureSets()
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
