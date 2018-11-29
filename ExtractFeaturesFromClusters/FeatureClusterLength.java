import java.util.ArrayList;

public class FeatureClusterLength 
{
	private ArrayList<Integer> lengths;
	
	public FeatureClusterLength()
	{
		lengths = new ArrayList<Integer>();
	}
	
	public FeatureClusterLength(ArrayList<Integer> aList)
	{
		lengths = aList;
	}
	
	public void addALength(int aLength)
	{
		lengths.add(aLength);
	}
	
	public ArrayList<Integer> getLengths()
	{
		return lengths;
	}
	
	@Override	
	public String toString()
	{
		return lengths.toString();
	}
}
