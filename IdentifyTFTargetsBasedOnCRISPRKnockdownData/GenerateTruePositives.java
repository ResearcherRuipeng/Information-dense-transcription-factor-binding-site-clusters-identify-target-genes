import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;

public class GenerateTruePositives 
{
	public static void main(String[] args)
	{
		String truePosWithChipseq = args[0];
		String outputFilePath = args[1];
		
		TreeMap<Double, String> coeffToNameAndEnsemTruePos = readTruePosWithChipseqFile(truePosWithChipseq);		
		output(coeffToNameAndEnsemTruePos, outputFilePath);
	}
	
	private static void output(TreeMap<Double, String> coeffToNameAndEnsemTruePos, String outputFilePath)
	{
		try
		{
			PrintWriter pw = new PrintWriter(outputFilePath);
			Map<Double, String> descendingMapView = coeffToNameAndEnsemTruePos.descendingMap();
			Set<Entry<Double, String>> setView = descendingMapView.entrySet();
			Iterator<Entry<Double, String>> it = setView.iterator();
			while(it.hasNext())
			{
				Entry<Double, String> anEntry = it.next();
				pw.println(anEntry.getKey() + "," + anEntry.getValue());
			}
			pw.close();
		}
		catch(IOException e)
		{
			System.out.println("Cannot open the file: " + outputFilePath);
			System.exit(1);
		}
	}
	
	
	private static TreeMap<Double, String> readTruePosWithChipseqFile(String path)
	{
		TreeMap<Double, String> coeffToNameAndEnsemTruePos = new TreeMap<Double, String>();
		
		try
		{
			Scanner sc = new Scanner(new File(path));
			while(sc.hasNextLine())
			{
				String aLine = sc.nextLine();
				String[] fields = aLine.split("\t");
				coeffToNameAndEnsemTruePos.put(Double.parseDouble(fields[5]), fields[3] + "," + fields[4]);
			}
			sc.close();
		}
		catch(IOException e)
		{
			System.out.println("Cannot open the file: " + path);
			System.exit(1);
		}
		
		return coeffToNameAndEnsemTruePos;
	}
}
