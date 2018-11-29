import java.util.ArrayList;
import java.util.Arrays;

public class FeatureNumberStrongSitesTFCluster 
{
	private ArrayList<int[]> numbersOfStrongSites;
	
	public FeatureNumberStrongSitesTFCluster()
	{
		numbersOfStrongSites = new ArrayList<int[]>();
	}
	
	public FeatureNumberStrongSitesTFCluster(ArrayList<int[]> aList)
	{
		numbersOfStrongSites = aList;
	}
	
	public void addACluster(int[] aCluster)
	{
		numbersOfStrongSites.add(aCluster);
	}
	
	public ArrayList<int[]> getNumbersOfStrongSites()
	{
		return numbersOfStrongSites;
	}
	
	@Override	
	public String toString()
	{
		StringBuilder sb = new StringBuilder("[");
		int[] numbersOfStrongSitesOneCluster = null;
		int numberOfClusters = numbersOfStrongSites.size();
		for(int i = 0; i < numberOfClusters; i++)
		{
			numbersOfStrongSitesOneCluster = numbersOfStrongSites.get(i);
			sb.append(Arrays.toString(numbersOfStrongSitesOneCluster));
			if(i != numberOfClusters - 1)
			{
				sb.append("; ");
			}			
		}
		sb.append("]");
		return sb.toString();
	}
}
