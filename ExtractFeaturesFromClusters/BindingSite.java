import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

public class BindingSite 
{
	private int position;
	private boolean ifOnPositiveStrand;
	private String sequence;
	private String majorTF;
	private double majorRi;
	private HashMap<String, Double> TFandRi;
	
	public BindingSite(int aPosition, boolean aBoolean, String aSequence, String TF, Double Ri)
	{
		position = aPosition;
		ifOnPositiveStrand = aBoolean;
		sequence = aSequence.toUpperCase();
		majorTF = TF;
		majorRi = Ri;
		TFandRi = new HashMap<String, Double>();
		TFandRi.put(TF, Ri);
	}
	
	public BindingSite(int aPosition, boolean aBoolean, String aSequence, String TF, Double Ri, HashMap<String, Double> someTFandRi)
	{
		position = aPosition;
		ifOnPositiveStrand = aBoolean;
		sequence = aSequence.toUpperCase();
		majorTF = TF;
		majorRi = Ri;
		TFandRi = someTFandRi;
	}
	
	public HashMap<String, Double> getTFsAndICs(HashMap<String, Double> modelToRseq, double threshold)
	{
		HashMap<String, Double> TFsAndICs = new HashMap<String, Double>();
		if(threshold == 0.0)
		{
			TFsAndICs.put(majorTF.split("-")[0], majorRi);
		}
		else if(majorRi >= modelToRseq.get(majorTF) * threshold)
		{
			TFsAndICs.put(majorTF.split("-")[0], majorRi);
		}
		
		Set<Entry<String, Double>> otherTFs = TFandRi.entrySet();
		Iterator<Entry<String, Double>> it = otherTFs.iterator();
		while(it.hasNext())
		{
			Entry<String, Double> anEntry = it.next();
			String key = anEntry.getKey();
			Double value = anEntry.getValue();
			
			if(threshold == 0.0)
			{
				String TFName = key.split("-")[0];
				Double oldRi = TFsAndICs.put(TFName, value);
				if(oldRi != null)
				{
					if(oldRi > value)
					{
						TFsAndICs.put(TFName, oldRi);
					}
				}
			}
			else if(value >= modelToRseq.get(key) * threshold)
			{
				String TFName = key.split("-")[0];
				Double oldRi = TFsAndICs.put(TFName, value);
				if(oldRi != null)
				{
					if(oldRi > value)
					{
						TFsAndICs.put(TFName, oldRi);
					}
				}
			}
		}
		
		return TFsAndICs;
	}
	
	public TreeSet<String> getTFs(HashMap<String, Double> modelToRseq, double threshold)
	{
		TreeSet<String> TFs = new TreeSet<String>();
		if(threshold == 0.0)
		{
			TFs.add(majorTF.split("-")[0]);
		}
		else if(majorRi >= modelToRseq.get(majorTF) * threshold)
		{
			TFs.add(majorTF.split("-")[0]);
		}
		/*TFs.add(majorTF.split("-")[0]);*/
		Set<Entry<String, Double>> otherTFs = TFandRi.entrySet();
		Iterator<Entry<String, Double>> it = otherTFs.iterator();
		while(it.hasNext())
		{
			Entry<String, Double> anEntry = it.next();
			String key = anEntry.getKey();
			
			if(threshold == 0.0)
			{
				TFs.add(key.split("-")[0]);
			}
			else if(anEntry.getValue() >= modelToRseq.get(key) * threshold)
			{
				TFs.add(key.split("-")[0]);
			}
		}
		
		return TFs;
	}
	
	@Override
	public String toString()
	{		
		if(ifOnPositiveStrand)
		{
			return position + "\t+\t" + sequence + "\t" + majorTF + "\t" + majorRi + "\t" + TFandRi.toString() + "\n";
		}
		else
		{
			return position + "\t-\t" + sequence + "\t" + majorTF + "\t" + majorRi + "\t" + TFandRi.toString() + "\n";
		}
	}
	
	public void setPosition(int aPosition)
	{
		position = aPosition;
	}
	
	public void setStrand(boolean bool)
	{
		ifOnPositiveStrand = bool;
	}
	
	public void setSequence(String aSequence)
	{
		sequence = aSequence;
	}
	
	public String getMajorTF()
	{
		return majorTF;
	}
	
	public double getMajorRi()
	{
		return majorRi;
	}
	
	public void setMajorTF(String newMajorTF)
	{
		majorTF = newMajorTF;
	}
	
	public void setMajorRi(double newMajorRi)
	{
		majorRi = newMajorRi;
	}
	
	public HashMap<String, Double> getTFandRi()
	{
		return TFandRi;
	}
	
	public int getPosition()
	{
		return position;
	}
	
	public String getSequence()
	{
		return sequence;
	}
	
	/*public String getReverseSequence()
	{
		StringBuilder reverseSequence = new StringBuilder();
		for(int i = sequence.length() - 1; i >= 0; i--)
		{
			switch(sequence.charAt(i))
			{
			case 'A':
				reverseSequence.append('T');
			case 'C':
				reverseSequence.append('G');
			case 'G':
				reverseSequence.append('C');
			case 'T':
				reverseSequence.append('A');
			}
		}
		
		return reverseSequence.toString();
	}*/
	
	@Override
	public boolean equals(Object otherObject)
	{
		if(this == otherObject)
			return true;
		if(otherObject == null)
			return false;
		if(!(otherObject instanceof BindingSite))
			return false;
		
		BindingSite other = (BindingSite) otherObject;
		return (position == other.position) && (sequence.length() == other.getSequence().length())/* && (sequence.equals(other.getSequence())) && (sequence.equals(other.getReverseSequence()))*/;
	}
}
