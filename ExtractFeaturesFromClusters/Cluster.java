import java.util.Iterator;
import java.util.TreeSet;

public class Cluster implements Comparable<Cluster> {
	private BindingSite centerSite;
	private TreeSet<BindingSite> otherSites;
	private double totalIC;
	private int smallCoordinate;
	private int largeCoordinate;
	private int length;
	private int numberOfSites;
	
	public Cluster(BindingSite aCenterSite, TreeSet<BindingSite> someSites, double NIC) {
		centerSite = aCenterSite;
		otherSites = someSites;
		totalIC = NIC;
	}
	
	public Cluster(BindingSite aCenterSite, TreeSet<BindingSite> someSites, double NIC, int aMinCoordinate, int aMaxCoordinate, int aLength, int aNumberOfSites) {
		centerSite = aCenterSite;
		otherSites = someSites;
		totalIC = NIC;
		smallCoordinate = aMinCoordinate;
		largeCoordinate = aMaxCoordinate;
		length = aLength;
		numberOfSites = aNumberOfSites;
	}
	
	public boolean intersect(int smallAbsoluteIndexOfAPromoter, int largeAbsoluteIndexOfAPromoter, boolean ifAbsoluteCoor, int smallCoordinatePromoterRegion)
	{
		/*if((smallCoordinate >= smallRelativeIndexOfAPromoter && smallCoordinate < largeRelativeIndexOfAPromoter) || (largeCoordinate <= largeRelativeIndexOfAPromoter && largeCoordinate > smallRelativeIndexOfAPromoter))*/
		int middlePointCoordinate = (smallCoordinate + largeCoordinate) / 2;
		if(!ifAbsoluteCoor)
		{
			middlePointCoordinate += smallCoordinatePromoterRegion;
		}
		if(middlePointCoordinate >= smallAbsoluteIndexOfAPromoter && middlePointCoordinate <= largeAbsoluteIndexOfAPromoter)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	
	public int getSmallCoordinate()
	{
		return smallCoordinate;
	}
	
	public int getLargeCoordinate()
	{
		return largeCoordinate;
	}
	
	public int getLength()
	{
		return length;
	}
	
	public int getNumberOfSites()
	{
		return numberOfSites;
	}
	
	public BindingSite getCenterSite()
	{
		return centerSite;
	}
	
	public TreeSet<BindingSite> getOtherSites()
	{
		return otherSites;
	}
	
	public double getIC()
	{
		return totalIC;
	}
	
	public double getICWhereEachSiteIsCountedOnce()
	{
		double totalICWhereEachSiteIsCountedOnce = centerSite.getMajorRi();
		for(BindingSite aSite : otherSites) {
			totalICWhereEachSiteIsCountedOnce += aSite.getMajorRi();
		}
		
		return totalICWhereEachSiteIsCountedOnce;
	}
	
	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();
		str.append("***************************************************************************\n");
		int[] range = getRange();
		str.append("The coordinates of the cluster: " + range[0] + "-" + range[1] + ". Length: " + range[2] + "\n");
		str.append("The total number of binding sites in this cluster: " + (otherSites.size() + 1) + "\n");
		str.append("IC: " + totalIC + "\n");
		str.append(centerSite.toString());
		Iterator<BindingSite> it = otherSites.iterator();
		while(it.hasNext()) {
			str.append(it.next().toString());
		}
		str.append("***************************************************************************\n");
		
		return str.toString();
	}
	
	public int[] getRange() {
		int minCoordinate = centerSite.getPosition();
		int maxCoordinate = centerSite.getPosition();		
		
		Iterator<BindingSite> it = otherSites.iterator();
		while(it.hasNext()) {
			BindingSite next = it.next();
			if(next.getPosition() < minCoordinate) {
				minCoordinate = next.getPosition();
			}
			else if(next.getPosition() > maxCoordinate) {
				maxCoordinate = next.getPosition();
			}
		}
		int length = maxCoordinate - minCoordinate + 1;
		
		int[] range = {minCoordinate, maxCoordinate, length};
		return range;
	}
	
	@Override
	public int compareTo(Cluster other) {
		return (totalIC < other.totalIC) ? -1 : ((totalIC > other.totalIC) ? 1 : 0);
	}
	
	@Override public boolean equals(Object otherObject) {
		if(this == otherObject)
			return true;
		if(otherObject == null)
			return false;
		if(!(otherObject instanceof Cluster))
			return false;
		
		Cluster other = (Cluster) otherObject;
		return (centerSite.equals(other.centerSite)) && (totalIC == other.totalIC) && (otherSites.size() == other.otherSites.size());
	}
	
	@Override
	public int hashCode() {
		return centerSite.hashCode() + Double.hashCode(totalIC) + Integer.hashCode(otherSites.size());
	}
	
	public void recomputeTotalIC() {
		totalIC = 0.0;
		for(BindingSite aSite : otherSites) {
			totalIC += centerSite.getMajorRi() + aSite.getMajorRi();
		}
	}
	
	public BindingSite newCenterSite() {
		double maxIC = 0.0;
		BindingSite siteWithMaxIC = null;
		Iterator<BindingSite> it = otherSites.iterator();
		while(it.hasNext()) {
			BindingSite next = it.next();
			if(next.getMajorRi() > maxIC) {
				siteWithMaxIC = next;
				maxIC = next.getMajorRi();
			}
		}
		
		return siteWithMaxIC;
	}
	
	public static int merge(Cluster cluster1, Cluster cluster2)
	{
		if(cluster1.centerSite.getMajorRi() >= cluster2.centerSite.getMajorRi())
		{
			cluster1.otherSites.add(cluster2.centerSite);
			cluster1.otherSites.addAll(cluster2.otherSites);
			cluster1.otherSites.remove(cluster1.centerSite);
			cluster1.recomputeTotalIC();
			
			return 2;
		}
		else
		{
			cluster2.otherSites.add(cluster1.centerSite);
			cluster2.otherSites.addAll(cluster1.otherSites);
			cluster2.otherSites.remove(cluster2.centerSite);
			cluster2.recomputeTotalIC();
			
			return 1;
		}
	}
	
	/*public Cluster merge(Cluster cluster2) {
		BindingSite newCenterSite = null;
		TreeSet<BindingSite> newOtherSites = new TreeSet<BindingSite>(new BindingSiteComparator());
		if(centerSite.IC > cluster2.centerSite.IC) {
			newCenterSite = new BindingSite(centerSite);
			newOtherSites.addAll(otherSites);
			Iterator<BindingSite> it = otherSites.iterator();
			while(it.hasNext()) {
				newOtherSites.add(it.next());
			}
			Iterator<BindingSite> it2 = cluster2.otherSites.iterator();
			while(it2.hasNext()) {
				newOtherSites.add(it2.next());
			}
			newOtherSites.addAll(cluster2.otherSites);
			newOtherSites.add(cluster2.centerSite);
			newOtherSites.remove(centerSite);
		}
		else if(cluster2.centerSite.IC > centerSite.IC) {
			newCenterSite = new BindingSite(cluster2.centerSite);
			newOtherSites.addAll(otherSites);
			newOtherSites.addAll(cluster2.otherSites);
			Iterator<BindingSite> it = otherSites.iterator();
			while(it.hasNext()) {
				newOtherSites.add(it.next());
			}
			Iterator<BindingSite> it2 = cluster2.otherSites.iterator();
			while(it2.hasNext()) {
				newOtherSites.add(it2.next());
			}
			newOtherSites.add(centerSite);
			newOtherSites.remove(cluster2.centerSite);
		}
		else {
			if(otherSites.size() >= cluster2.otherSites.size()) {
				newCenterSite = new BindingSite(centerSite);
				newOtherSites.addAll(otherSites);
				newOtherSites.addAll(cluster2.otherSites);
				Iterator<BindingSite> it = otherSites.iterator();
				while(it.hasNext()) {
					newOtherSites.add(it.next());
				}
				Iterator<BindingSite> it2 = cluster2.otherSites.iterator();
				while(it2.hasNext()) {
					newOtherSites.add(it2.next());
				}
				newOtherSites.add(cluster2.centerSite);
				newOtherSites.remove(centerSite);
			}
			else {
				newCenterSite = new BindingSite(cluster2.centerSite);
				newOtherSites.addAll(otherSites);
				newOtherSites.addAll(cluster2.otherSites);
				Iterator<BindingSite> it = otherSites.iterator();
				while(it.hasNext()) {
					newOtherSites.add(it.next());
				}
				Iterator<BindingSite> it2 = cluster2.otherSites.iterator();
				while(it2.hasNext()) {
					newOtherSites.add(it2.next());
				}
				newOtherSites.add(centerSite);
				newOtherSites.remove(cluster2.centerSite);
			}
		}
		
		double newTotalIC = 0.0;
		for(BindingSite aSite : newOtherSites) {
			newTotalIC += newCenterSite.IC + aSite.IC;
		}
		
		return new Cluster(newCenterSite, newOtherSites, newTotalIC);
	}*/
}
