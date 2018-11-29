import java.util.ArrayList;
import java.util.TreeSet;

public class FeatureOrderOfSitesCluster 
{
	private ArrayList<ArrayList<TreeSet<String>>> orderOfSites;
	
	public FeatureOrderOfSitesCluster()
	{
		orderOfSites = new ArrayList<ArrayList<TreeSet<String>>>();
	}
	
	public FeatureOrderOfSitesCluster(ArrayList<ArrayList<TreeSet<String>>> aList)
	{
		orderOfSites = aList;
	}
	
	public void addACluster(ArrayList<TreeSet<String>> aCluster)
	{
		orderOfSites.add(aCluster);
	}
	
	public ArrayList<ArrayList<TreeSet<String>>> getOrderOfSites()
	{
		return orderOfSites;
	}
	
	@Override	
	public String toString()
	{
		StringBuilder sb = new StringBuilder("[");
		ArrayList<TreeSet<String>> orderOfSitesOneCluster = null;
		TreeSet<String> aBindingSite = null;
		int numberOfClusters = orderOfSites.size();		
		int numberOfSitesOneCluster = -1;
		
		for(int i = 0; i < numberOfClusters; i++)
		{
			orderOfSitesOneCluster = orderOfSites.get(i);
			sb.append("[");
			numberOfSitesOneCluster = orderOfSitesOneCluster.size();
			for(int j = 0; j < numberOfSitesOneCluster; j++)
			{
				aBindingSite = orderOfSitesOneCluster.get(j);
				sb.append(aBindingSite.toString());
				if(j != numberOfSitesOneCluster - 1)
				{
					sb.append("; ");
				}				
			}
			sb.append("]");
			
			if(i != numberOfClusters - 1)
			{
				sb.append(": ");
			}			
		}
		sb.append("]");
		return sb.toString();
	}
}
