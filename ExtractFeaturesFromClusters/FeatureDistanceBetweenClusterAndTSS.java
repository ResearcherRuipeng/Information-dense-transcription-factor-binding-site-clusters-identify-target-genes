import java.util.ArrayList;

public class FeatureDistanceBetweenClusterAndTSS
{
	private ArrayList<Integer> distances;
	
	public FeatureDistanceBetweenClusterAndTSS()
	{
		distances = new ArrayList<Integer>();
	}
	
	public FeatureDistanceBetweenClusterAndTSS(ArrayList<Integer> aList)
	{
		distances = aList;
	}
	
	public void addADistance(int aDistance)
	{
		distances.add(aDistance);
	}
	
	public ArrayList<Integer> getDistances()
	{
		return distances;
	}
	
	@Override	
	public String toString()
	{
		return distances.toString();
	}
}
