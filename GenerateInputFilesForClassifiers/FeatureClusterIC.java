import java.util.ArrayList;

public class FeatureClusterIC 
{
	private ArrayList<Double> ICs;
	
	public FeatureClusterIC()
	{
		ICs = new ArrayList<Double>();
	}
	
	public FeatureClusterIC(ArrayList<Double> aList)
	{
		ICs = aList;
	}
	
	public void addAnIC(double anIC)
	{
		ICs.add(anIC);
	}
	
	public ArrayList<Double> getICs()
	{
		return ICs;
	}
	
	@Override	
	public String toString()
	{
		return ICs.toString();
	}
}
