import java.util.ArrayList;
import java.util.Arrays;

public class FeatureSumICSitesTFCluster 
{
	private ArrayList<double[]> sumsOfICOfSites;
	
	public FeatureSumICSitesTFCluster()
	{
		sumsOfICOfSites = new ArrayList<double[]>();
	}
	
	public FeatureSumICSitesTFCluster(ArrayList<double[]> aList)
	{
		sumsOfICOfSites = aList;
	}
	
	public void addACluster(double[] aCluster)
	{
		sumsOfICOfSites.add(aCluster);
	}
	
	public ArrayList<double[]> getSumsOfICOfSites()
	{
		return sumsOfICOfSites;
	}
	
	@Override	
	public String toString()
	{
		StringBuilder sb = new StringBuilder("[");
		double[] sumsOfICOfSitesOneCluster = null;
		int numberOfClusters = sumsOfICOfSites.size();
		for(int i = 0; i < numberOfClusters; i++)
		{
			sumsOfICOfSitesOneCluster = sumsOfICOfSites.get(i);
			sb.append(Arrays.toString(sumsOfICOfSitesOneCluster));
			if(i != numberOfClusters - 1)
			{
				sb.append("; ");
			}			
		}
		sb.append("]");
		return sb.toString();
	}
}
