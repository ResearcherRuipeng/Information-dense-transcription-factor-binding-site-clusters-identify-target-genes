import java.util.Comparator;

public class ClusterComparator implements Comparator<Cluster>
{
	public int compare(Cluster cluster1, Cluster cluster2)
	{
		if(cluster1.getCenterSite().getPosition() < cluster2.getCenterSite().getPosition())
		{
			return -1;
		}
		else if(cluster1.getCenterSite().getPosition() > cluster2.getCenterSite().getPosition())
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}
