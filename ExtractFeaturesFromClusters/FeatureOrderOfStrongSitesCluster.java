import java.util.ArrayList;
import java.util.TreeSet;

public class FeatureOrderOfStrongSitesCluster 
{
	private ArrayList<ArrayList<TreeSet<String>>> orderOfStrongSites;
	
	public FeatureOrderOfStrongSitesCluster()
	{
		orderOfStrongSites = new ArrayList<ArrayList<TreeSet<String>>>();
	}
	
	public FeatureOrderOfStrongSitesCluster(ArrayList<ArrayList<TreeSet<String>>> aList)
	{
		orderOfStrongSites = aList;
	}
	
	public void addACluster(ArrayList<TreeSet<String>> aCluster)
	{
		orderOfStrongSites.add(aCluster);
	}
	
	public ArrayList<ArrayList<TreeSet<String>>> getOrderOfStrongSites()
	{
		return orderOfStrongSites;
	}
	
	@Override	
	public String toString()
	{
		StringBuilder sb = new StringBuilder("[");
		ArrayList<TreeSet<String>> orderOfStrongSitesOneCluster = null;
		TreeSet<String> aBindingSite = null;
		int numberOfClusters = orderOfStrongSites.size();		
		int numberOfStrongSitesOneCluster = -1;
		
		for(int i = 0; i < numberOfClusters; i++)
		{
			orderOfStrongSitesOneCluster = orderOfStrongSites.get(i);
			sb.append("[");
			numberOfStrongSitesOneCluster = orderOfStrongSitesOneCluster.size();
			for(int j = 0; j < numberOfStrongSitesOneCluster; j++)
			{
				aBindingSite = orderOfStrongSitesOneCluster.get(j);
				sb.append(aBindingSite.toString());
				if(j != numberOfStrongSitesOneCluster - 1)
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
