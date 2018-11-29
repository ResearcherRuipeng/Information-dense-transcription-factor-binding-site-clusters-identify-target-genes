import java.util.ArrayList;

public class FeatureNumberSitesTFCluster_oneiPWM
{
	private ArrayList<Integer> numbersOfSites;
	
	public FeatureNumberSitesTFCluster_oneiPWM()
	{
		numbersOfSites = new ArrayList<Integer>();
	}
	
	public FeatureNumberSitesTFCluster_oneiPWM(ArrayList<Integer> aList)
	{
		numbersOfSites = aList;
	}
	
	public void addACluster(int aCluster)
	{
		numbersOfSites.add(aCluster);
	}
	
	public ArrayList<Integer> getNumbersOfSites()
	{
		return numbersOfSites;
	}
	
	@Override	
	public String toString()
	{
		return numbersOfSites.toString();		
	}
}
