import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;

public class Main
{
	public static void main(String[] args)
	{
		String bindingSiteFolderPath = args[0];
		String promoterRegionFolderPath = args[1];
		String outputFolderPath = args[2];		
		
		ArrayList<TreeMap<String, ArrayList<int[]>>> DNaseRegionFiles = new ArrayList<TreeMap<String, ArrayList<int[]>>>();
		for(int i = 3; i < args.length; i++)
		{
			TreeMap<String, ArrayList<int[]>> DNaseRegionsInOneFile = new TreeMap<String, ArrayList<int[]>>(); 
			DNaseRegionsInOneFile = readDNaseTrackFile(args[i]);
			//test
			/*Set<String> keySet = DNaseRegionsInOneFile.keySet();
			Iterator<String> it = keySet.iterator();
			while(it.hasNext())
			{
				System.out.println(it.next());
			}*/
			DNaseRegionFiles.add(DNaseRegionsInOneFile);
		}		
		
		File[] allBindingSiteFiles = new File(bindingSiteFolderPath).listFiles();
		for(int i = 0; i < allBindingSiteFiles.length; i++)
		{
			String bindingSiteFilename = allBindingSiteFiles[i].getName();
			ArrayList<BindingSite> allBindingSitesInOneFile = readOneBindingSiteFile(allBindingSiteFiles[i]);
			TreeMap<String, String> promoterRegionToChr = changeToAbsoluteCoordinates(bindingSiteFilename, allBindingSitesInOneFile, promoterRegionFolderPath);
						
			TreeSet<BindingSite> allBindingSitesInDNaseRegionsInOneFile = new TreeSet<BindingSite>(new BindingSiteComparator());
			for(int j = 0; j < DNaseRegionFiles.size(); j++)
			{				
				intersectBindingSitesWithDNaseRegions(bindingSiteFilename, allBindingSitesInOneFile, DNaseRegionFiles.get(j), promoterRegionToChr, allBindingSitesInDNaseRegionsInOneFile);
			}
			
			output(outputFolderPath, bindingSiteFilename, allBindingSitesInDNaseRegionsInOneFile);			
		}		
		
		System.out.println("Done");
	}
	
	private static void output(String outputFolderPath, String bindingSiteFilename, TreeSet<BindingSite> allBindingSitesInDNaseRegionsInOneFile)
	{		
		if(allBindingSitesInDNaseRegionsInOneFile.size() == 0)
		{
			System.out.println("The binding site file " + bindingSiteFilename + " does not overlap with DNase sites");
			return;
		}

		try
		{
			PrintWriter output = new PrintWriter(outputFolderPath + "/" + bindingSiteFilename);
			Iterator<BindingSite> it = allBindingSitesInDNaseRegionsInOneFile.iterator();
			while(it.hasNext())
			{
				output.print(it.next().toString());
			}
			output.close();
		}
		catch(Exception e)
		{
			System.out.println("The output file cannot be opened: " + bindingSiteFilename);
			System.exit(1);
		}		
	}
	
	private static void intersectBindingSitesWithDNaseRegions(String bindingSiteFilename, ArrayList<BindingSite> allBindingSitesInOneFile, TreeMap<String, ArrayList<int[]>> DNaseRegions, TreeMap<String, String> promoterRegionToChr, TreeSet<BindingSite> allBindingSitesInDNaseRegionsInOneFile)
	{				
		int pointerToFirstDNaseRegionForOneSite = 0;
		int pointerToCurrentDNaseRegion = 0;
		ArrayList<int[]> DNaseRegionsInOneChr = DNaseRegions.get(promoterRegionToChr.get(bindingSiteFilename));
		if(DNaseRegionsInOneChr == null)
		{
			return;
		}
		boolean ifSiteIsInOneDNaseRegion = false;
		
		Iterator<BindingSite> it = allBindingSitesInOneFile.iterator();
		while(it.hasNext())
		{
			BindingSite aSite = it.next();
			for(int i = pointerToFirstDNaseRegionForOneSite; i < DNaseRegionsInOneChr.size(); i++)
			{
				pointerToCurrentDNaseRegion = i;
				ifSiteIsInOneDNaseRegion = aSite.IsInOneDNaseRegion(DNaseRegionsInOneChr.get(pointerToCurrentDNaseRegion));
				if(ifSiteIsInOneDNaseRegion)
				{
					allBindingSitesInDNaseRegionsInOneFile.add(aSite);
					break;
				}
			}

			if(ifSiteIsInOneDNaseRegion)
			{
				pointerToFirstDNaseRegionForOneSite = pointerToCurrentDNaseRegion;
			}
		}		
		
		return;
	}
	
	private static TreeMap<String, String> changeToAbsoluteCoordinates(String bindingSiteFilename, ArrayList<BindingSite> allBindingSitesInOneFile, String promoterRegionFolderPath)
	{
		TreeMap<String, String> promoterRegionToChr = new TreeMap<String, String>();		
		
		try
		{
			Scanner input = new Scanner(new File(promoterRegionFolderPath + "/" + bindingSiteFilename));
			String promoterRegion = input.nextLine();
			String[] fields = promoterRegion.split("\t");
			promoterRegionToChr.put(bindingSiteFilename, fields[0]);
			int smallCoor = Integer.parseInt(fields[1]);

			Iterator<BindingSite> it = allBindingSitesInOneFile.iterator();
			while(it.hasNext())
			{
				BindingSite aSite = it.next();
				aSite.changeToAbsoluteCoordinate(smallCoor);
			}				
			input.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("The promoter region file cannot be opened: " + bindingSiteFilename);
			System.exit(1);
		}		
		
		return promoterRegionToChr;
	}	
	
	private static ArrayList<BindingSite> readOneBindingSiteFile(File oneBindingSiteFile)
	{
		ArrayList<String> allSiteLines = new ArrayList<String>();
		ArrayList<BindingSite> allBindingSitesInOneFile = null;
		
		try 
		{
			Scanner input = new Scanner(oneBindingSiteFile);
			while(input.hasNextLine()) 
			{
				allSiteLines.add(input.nextLine());
			}
			input.close();
			
			allBindingSitesInOneFile = decomposeAllSites(allSiteLines);
		}
		catch(FileNotFoundException e) 
		{
			System.out.println("The binding site file cannot be found: " + oneBindingSiteFile.getAbsolutePath());
			System.exit(1);
		}
		
		return allBindingSitesInOneFile;
	}
	
	private static ArrayList<BindingSite> decomposeAllSites(ArrayList<String> siteLines) 
	{
		ArrayList<BindingSite> allBindingSitesInOneFile = new ArrayList<BindingSite>();
		Iterator<String> it = siteLines.iterator();
		BindingSite aSite;
		String aSiteLine;
		String[] fields;		
		String someTFsandRis = "";
		
		while(it.hasNext()) 
		{
			someTFsandRis = "";
			aSiteLine = it.next();
			fields = aSiteLine.split("\t");			
			
			String TF = fields[3];
			if(fields.length == 6)
			{
				someTFsandRis += fields[5];
			}								

			if(fields[1].equals("+")) 
			{					
				aSite = new BindingSite(Integer.parseInt(fields[0]), Double.parseDouble(fields[4]), true, fields[2], TF, someTFsandRis);
			} 
			else 
			{
				aSite = new BindingSite(Integer.parseInt(fields[0]), Double.parseDouble(fields[4]), false, fields[2], TF, someTFsandRis);
			}
			allBindingSitesInOneFile.add(aSite);			
		}
		
		return allBindingSitesInOneFile;
	}
	
	private static TreeMap<String, ArrayList<int[]>> readDNaseTrackFile(String DNaseTrackFilePath)
	{
		TreeMap<String, ArrayList<int[]>> DNaseRegions = new TreeMap<String, ArrayList<int[]>>();
		String chr = null;
		ArrayList<int[]> DNaseRegionsInOneChr = null;
		int[] aDNaseRegion = null;
		
		try
		{
			Scanner input = new Scanner(new File(DNaseTrackFilePath));
			while(input.hasNext())
			{
				String aLine = input.nextLine();
				aDNaseRegion = new int[2];
				String[] fields = aLine.split("\t");
				
				if(!fields[0/*1*/].equals(chr))
				{
					if(chr != null)
					{
						DNaseRegions.put(chr, DNaseRegionsInOneChr);						
					}
					
					chr = fields[0/*1*/];
					DNaseRegionsInOneChr = new ArrayList<int[]>();
				}				
				
				aDNaseRegion[0] = Integer.parseInt(fields[1/*2*/]);
				aDNaseRegion[1] = Integer.parseInt(fields[2/*3*/]);
				DNaseRegionsInOneChr.add(aDNaseRegion);
			}
			DNaseRegions.put(chr, DNaseRegionsInOneChr);
			input.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println("Unable to open the file: " + DNaseTrackFilePath);
			System.exit(1);
		}
		
		return DNaseRegions;
	}
}

