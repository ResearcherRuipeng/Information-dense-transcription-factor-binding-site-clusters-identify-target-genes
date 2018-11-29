import java.util.ArrayList;

public class FeatureNumberStrongSitesTFCluster_oneiPWM 
{
	private ArrayList<Integer> numbersOfStrongSites;
	
	public FeatureNumberStrongSitesTFCluster_oneiPWM()
	{
		numbersOfStrongSites = new ArrayList<Integer>();
	}
	
	public FeatureNumberStrongSitesTFCluster_oneiPWM(ArrayList<Integer> aList)
	{
		numbersOfStrongSites = aList;
	}
	
	public void addACluster(int aCluster)
	{
		numbersOfStrongSites.add(aCluster);
	}
	
	public ArrayList<Integer> getNumbersOfStrongSites()
	{
		return numbersOfStrongSites;
	}
	
	@Override	
	public String toString()
	{		
		return numbersOfStrongSites.toString();
	}
}
