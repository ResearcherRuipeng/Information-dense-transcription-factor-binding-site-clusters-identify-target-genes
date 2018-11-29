
public class BindingSite {
	public int position;
	public double IC;	
	public boolean ifOnPlusStrand;
	public String sequence;
	public String TF;
	public String TFandRi;
	
	public BindingSite(int aPos, double anIC, boolean aBool, String aSequence, String aTF, String someTFsandRis) {
		position = aPos;
		IC = anIC;
		ifOnPlusStrand = aBool;
		sequence = aSequence;
		TF = aTF;
		TFandRi = someTFsandRis;
	}
	
	public BindingSite(BindingSite otherSite) {
		position = otherSite.position;
		IC = otherSite.IC;
		ifOnPlusStrand = otherSite.ifOnPlusStrand;
		sequence = otherSite.sequence;
		TF = otherSite.TF;
		TFandRi = otherSite.TFandRi;
	}
	
	public boolean IsInOneDNaseRegion(int[] aDNaseRegion)
	{
		if((position >= aDNaseRegion[0]) && (position/* + sequence.length() */<= aDNaseRegion[1]))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	
	public void changeToAbsoluteCoordinate(int smallCoorOfPromoterRegion)
	{
		position += smallCoorOfPromoterRegion;
	}
	
	@Override
	public String toString() {
		if(ifOnPlusStrand) {
			if(TFandRi.equals(""))
			{
				return position + "\t+\t" + sequence + "\t" + TF + "\t" + IC + "\n";
			}
			else
			{
				return position + "\t+\t" + sequence + "\t" + TF + "\t" + IC + "\t" + TFandRi.toString() + "\n";
			}									
		}
		else {						
			if(TFandRi.equals(""))
			{
				return position + "\t-\t" + sequence + "\t" + TF + "\t" + IC + "\n";
			}
			else
			{
				return position + "\t-\t" + sequence + "\t" + TF + "\t" + IC + "\t" + TFandRi.toString() + "\n";
			}						
		}
	}
	
	@Override 
	public boolean equals(Object otherObject) {
		if(this == otherObject)
			return true;
		if(otherObject == null)
			return false;
		if(!(otherObject instanceof BindingSite))
			return false;
		
		BindingSite other = (BindingSite) otherObject;
		return (position == other.position) && /*(IC == other.IC) && (ifOnPlusStrand == other.ifOnPlusStrand) && (TF.equals(other.TF)) && */(sequence.length() == other.sequence.length());
	}
	
	@Override
	public int hashCode() {
		return Integer.hashCode(position) + /*Double.hashCode(IC) + Boolean.hashCode(ifOnPlusStrand) + */new Integer(sequence.length()).hashCode()/* + TF.hashCode()*/;
	}
}
