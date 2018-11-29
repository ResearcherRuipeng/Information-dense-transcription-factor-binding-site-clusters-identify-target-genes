import java.util.ArrayList;

public class FeatureSumICStrongSitesTFCluster_oneiPWM
{
	private ArrayList<Double> sumsOfICOfSites;
	
	public FeatureSumICStrongSitesTFCluster_oneiPWM()
	{
		sumsOfICOfSites = new ArrayList<Double>();
	}
	
	public FeatureSumICStrongSitesTFCluster_oneiPWM(ArrayList<Double> aList)
	{
		sumsOfICOfSites = aList;
	}
	
	public void addACluster(double aCluster)
	{
		sumsOfICOfSites.add(aCluster);
	}
	
	public ArrayList<Double> getSumsOfICOfSites()
	{
		return sumsOfICOfSites;
	}
	
	@Override	
	public String toString()
	{		
		return sumsOfICOfSites.toString();
	}
}
