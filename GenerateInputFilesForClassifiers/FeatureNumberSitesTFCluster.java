import java.util.ArrayList;
import java.util.Arrays;

public class FeatureNumberSitesTFCluster 
{
	private ArrayList<int[]> numbersOfSites;
	
	public FeatureNumberSitesTFCluster()
	{
		numbersOfSites = new ArrayList<int[]>();
	}
	
	public FeatureNumberSitesTFCluster(ArrayList<int[]> aList)
	{
		numbersOfSites = aList;
	}
	
	public void addACluster(int[] aCluster)
	{
		numbersOfSites.add(aCluster);
	}
	
	public ArrayList<int[]> getNumbersOfSites()
	{
		return numbersOfSites;
	}
	
	@Override	
	public String toString()
	{
		StringBuilder sb = new StringBuilder("[");
		int[] numbersOfSitesOneCluster = null;
		int numberOfClusters = numbersOfSites.size();
		for(int i = 0; i < numberOfClusters; i++)
		{
			numbersOfSitesOneCluster = numbersOfSites.get(i);
			sb.append(Arrays.toString(numbersOfSitesOneCluster));
			if(i != numberOfClusters - 1)
			{
				sb.append("; ");
			}			
		}
		sb.append("]");
		return sb.toString();
	}
}
