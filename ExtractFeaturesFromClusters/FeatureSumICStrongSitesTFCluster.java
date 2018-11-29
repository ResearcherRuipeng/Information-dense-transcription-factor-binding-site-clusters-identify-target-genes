import java.util.ArrayList;
import java.util.Arrays;

public class FeatureSumICStrongSitesTFCluster 
{
	private ArrayList<double[]> sumsOfICOfStrongSites;
	
	public FeatureSumICStrongSitesTFCluster()
	{
		sumsOfICOfStrongSites = new ArrayList<double[]>();
	}
	
	public FeatureSumICStrongSitesTFCluster(ArrayList<double[]> aList)
	{
		sumsOfICOfStrongSites = aList;
	}
	
	public void addACluster(double[] aCluster)
	{
		sumsOfICOfStrongSites.add(aCluster);
	}
	
	public ArrayList<double[]> getSumsOfICOfSites()
	{
		return sumsOfICOfStrongSites;
	}
	
	@Override	
	public String toString()
	{
		StringBuilder sb = new StringBuilder("[");
		double[] sumsOfICOfStrongSitesOneCluster = null;
		int numberOfClusters = sumsOfICOfStrongSites.size();
		for(int i = 0; i < numberOfClusters; i++)
		{
			sumsOfICOfStrongSitesOneCluster = sumsOfICOfStrongSites.get(i);
			sb.append(Arrays.toString(sumsOfICOfStrongSitesOneCluster));
			if(i != numberOfClusters - 1)
			{
				sb.append("; ");
			}			
		}
		sb.append("]");
		return sb.toString();
	}
}
