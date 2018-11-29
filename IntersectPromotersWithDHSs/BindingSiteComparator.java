import java.util.Comparator;

public class BindingSiteComparator implements Comparator<BindingSite> {
	public int compare(BindingSite site1, BindingSite site2) {
		if(site1.position < site2.position) {
			return -1;
		}
		else if(site1.position > site2.position) {
			return 1;
		}
		else {
			if(site1.sequence.length() < site2.sequence.length())
			{
				return -1;
			}
			else if(site1.sequence.length() > site2.sequence.length())
			{
				return 1;
			}
			else
			{
				return 0;
			}
			/*if(site1.ifOnPlusStrand == true && site2.ifOnPlusStrand == false) {
				return -1;
			}
			else if(site1.ifOnPlusStrand == false && site2.ifOnPlusStrand == true) {
				return 1;
			}
			else {
				int sequenceComparison = site1.sequence.compareTo(site2.sequence);
				if(sequenceComparison != 0) {
					return sequenceComparison;
				}
				else {
					if(site1.IC < site2.IC) {
						return -1;
					}
					else if(site1.IC > site2.IC) {
						return 1;
					}
					else {
						if(site1.TF < site2.TF) {
							return -1;
						}
						else if(site1.TF > site2.TF) {
							return 1;
						}
						else {
							return 0;
						}
					}
				}
			}*/
		}
	}
}
