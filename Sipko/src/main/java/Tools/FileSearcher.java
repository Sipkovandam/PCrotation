package Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import PCA.MatrixStruct;
import STAR.StringFilter;

/**
 * @author Sipko
 *
 */
public class FileSearcher extends Script<FileSearcher> {
	// This script searches a folder (and all subdirectories) all for all files
	// containing a certain "String" and
	// outputs all filenames (including pathname) that contain this string into
	// a .txt file
	// same can be achieved with shell command : find . | grep
	// "whateverYouAreSearchingFor" | grep -v "whateverYouWantToExclude"

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String folderComment = "/root/directory/,/root/directory2/; MANDATORY UNLESS//Comma separated list of input folders to search files in";
	String folders = null;// "E:/Groningen/Test/STAR/STAR/,E:/Groningen/Test/PCA/Test2/";
	private String writeNameComment = "/root/directory/fastqs.txt; OPTIONAL //filename of the output file";//semitransient
	String writeName = null;// "E:/Groningen/Test/STAR/STAR/textsearch.txt";//semitransient
	private String searchStringsComment = "[\".txt\",\".fq\"]; MANDATORY UNLESS (fastQfilesFN) IS DEFINED //Comma separated list of strings that define which files are fastq files. Any files containing this string will be input files for the pipeline except those excluded by (forbiddenSearchStrings)";
	String[] searchStrings = new String[] { ".fq", ".fastq" };
	private String requiredStringFNComment = "/root/directory/requiredStrings.txt; OPTIONAL //filename of a file containing enter separated strings. One of these strings needs to be in the filename of searched files in order to be included in the results";
	String requiredStringFN = null;// "E:/Groningen/Test/STAR/STAR/RequiredStrings.txt";
	private String sampleNamesFNComment = "/root/directory/sampleNames.txt; OPTIONAL //filename of a file containing enter separated strings. This needs to be the name of the file (excluding root) after removing (removeBits)";
	String sampleNamesFN = null;// "E:/Groningen/Test/STAR/STAR/RequiredStrings.txt";
	private String sampleCheckRemoveAfterStringsComment = "[\"_R1_\",\".fq\"]; OPTIONAL //bits to remove from filename before checking if the sample should be included";
	private String[] sampleCheckRemoveAfterStrings = new String[]{"_R1.fq","_R2.fq","_1.fq","_2.fq"};
	private String forbiddenStringsComment = "[\"md5\",\"DISCARDED\"]; OPTIONAL //Comma separated list of strings. Any fastq file, defined by fastQSearchStringsComment, containing this string is not included in the analysis";
	String[] forbiddenStrings = new String[]{".md5","DISCARDED"};

	private class RequiredStringChecker implements StringFilter {
		private String[] requiredString = null;
		
		RequiredStringChecker(String[] requiredString)
		{
			this.requiredString=requiredString;
		}
		
		/**
		 * @param fileName
		 * @return True, if filter is pass, false if not
		 */
		@Override
		public boolean isPass(String fileName) {
			if (requiredString == null)
				return true;
			else// check if the other required string (any of them) is also present
				// in the filename
				for (int r = 0; r < requiredString.length; r++)
					if (requiredString[r] != null && fileName.contains(requiredString[r]))
						return true;
			return false;
		}
	}
	
	private class ForbiddenStringChecker implements StringFilter {
		private String[] forbiddenString = null;
		
		ForbiddenStringChecker(String[] forbiddenString)
		{
			this.forbiddenString=forbiddenString;
		}
		/**
		 * @param fileName
		 * @return True, if filter is pass, false if not
		 */
		@Override
		public boolean isPass(String fileName) {
			if (forbiddenString != null)
				for (int s = 0; s < forbiddenString.length; s++)
				{
					if (forbiddenString[s] != null && fileName.contains(forbiddenString[s]))
						return false;
				}
			return true;
		}
	}
	
	private class SampleNameChecker implements StringFilter {
		private Set<String> includeList = null;
		private String[] replaceElements = null;
		
		SampleNameChecker(Set<String> includeList,String[] replaceElements)
		{
			this.includeList=includeList;
			this.replaceElements=replaceElements;
		}
		/**
		 * @param fileName
		 * @return True, if filter is pass, false if not
		 */
		@Override
		public boolean isPass(String fileName) {
			
			if(includeList==null)
				return true;
		
			for (String replaceElement : replaceElements) {
				if(replaceElement!= null)
					fileName = fileName.replaceAll(replaceElement+".*", "");
			}
			if (includeList.contains(fileName))
			{
				return true;
			}
			return false;
		}
	}

	public FileSearcher() {
	}

	public FileSearcher(String folderName) {
		this.jsonFN = FileUtils.makeFolderNameEndWithSlash(folderName) + this.getNewJsonFN();
	}

	public static void main(String[] args) throws Exception {
		new FileSearcher().run(args);
	}

	public void run(String[] args) throws Exception {
		// java -jar -Xmx1g SearchFilesInDirectories.jar
		// folderName=/groups/umcg-pub/tmp04/public-rna-seq/
		// searchString=.fastq.gz
		// writefn=/local/groups/umcg-bios/scr01/Sipko/Juha/samples.txt

		if (args.length == 0)
			checkArgs(args);
		for (int a = 0; a < args.length; a++) {
			String arg = args[a].split("=")[0];
			String value = args[a].split("=")[1];
			switch (arg.toLowerCase()) {
			case "foldername":
				folders = parseString(value);
				break;
			case "folders":
				folders = parseString(value);
				break;
			case "searchstrings":
				searchStrings = value.split(",");
				break;
			case "writefn":
				writeName = value;
				break;
			case "forbiddenstrings":
				forbiddenStrings = value.split(",");
				break;
			case "requiredstringfn":
				requiredStringFN = parseString(value);
				break;
			default:
				System.out.println("Incorrect argument supplied: " + args[a] + "\n" + "exiting");
				checkArgs(args);
				System.exit(1);
			}
		}
		
		this.run();
		// String folderName =
		// "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated
		// gluten specific Tcell clones";
		// String writeName =
		// "/groups/umcg-wijmenga/prm02/data_projects/Celiac_Disease/gene_expression/RNASEQ/RNAseq_NG_data_08032016/HiSeq/stimulated
		// gluten specific Tcell clones/allSamples.txt";

	}

	public void run() {
		try {
			//init
			if (writeName == null)
				writeName = folders.split(",")[0] + this.getClassName() + "_result.txt";
			if (this.jsonFN == null)
				this.jsonFN = new File(writeName).getParent() + getJsonFN();
			writeConfig();
			String[] folderNames = folders.split(",");
			BufferedWriter writer = FileUtils.createWriter(writeName);
			BufferedWriter writerFailed = FileUtils.createWriter(FileUtils.removeExtention(writeName)+"_otherFiles.txt");
			
			Set<String> sampleNames = parseSampleNames(sampleNamesFN);
			SampleNameChecker sampleNameChecker = new SampleNameChecker(sampleNames,sampleCheckRemoveAfterStrings);
			
			ForbiddenStringChecker forbiddenStringChecker = new ForbiddenStringChecker(forbiddenStrings);
			List<StringFilter> checks = new ArrayList<>();
			
			RequiredStringChecker searchStringChecker = new RequiredStringChecker(searchStrings);
			
			final String[] requiredStringsFinal = readRequiredStrings(requiredStringFN);
			RequiredStringChecker requiredStringChecker = new RequiredStringChecker(requiredStringsFinal);
			
			checks.add(forbiddenStringChecker);
			checks.add(searchStringChecker);
			checks.add(requiredStringChecker);
			checks.add(sampleNameChecker);
			
			Stream.of(folderNames).forEach(folderName -> searchDirectory(new File(folderName), writer, writerFailed, searchStrings,
					requiredStringsFinal, forbiddenStrings,checks));
			
			writerFailed.close();
			writer.close();
			p("File writen to:" + writeName);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private Set<String> parseSampleNames(String sampleNamesFn) throws FileNotFoundException, IOException {
		if(sampleNamesFn==null)
			return null;
		Set<String> sampleNames = new HashSet<>();
		BufferedReader sampleNamesReader = FileUtils.createReader(sampleNamesFn);
		sampleNamesReader.lines().forEach(line -> sampleNames.add(line));
		return sampleNames;
	}

	private String[] readRequiredStrings(String requiredStringFN) throws FileNotFoundException, IOException {
		if (requiredStringFN == null)
			return null;
		BufferedReader reader = FileUtils.createReader(requiredStringFN);
		String[] requiredStrings = reader.lines().map(line -> {
			return line;
		}).collect(Collectors.toList()).toArray(new String[] {});

		return requiredStrings;
	}



	public void searchDirectory(File directory, BufferedWriter writer, BufferedWriter writerFailed, String[] searchString,
			String[] requiredString, String[] forbiddenString, List<StringFilter> checks){

		try {
			File[] files = directory.listFiles();
			for (File file : files) {
				if (file.isDirectory())
					searchDirectory(file, writer, writerFailed, searchString, requiredString, forbiddenString, checks);
				else {
					String fileName = file.getName();
					boolean include = checkInclude(searchString, fileName, requiredString, forbiddenString, checks);

					if (include)
						writer.write(file.getAbsolutePath() + "\n");
					else
						writerFailed.write(file.getAbsolutePath() + "\n");
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private boolean checkInclude(String[] searchString, String fileName, String[] requiredString,
			String[] forbiddenString, List<StringFilter> checks) throws Exception {

		boolean include = true;
		//-1 = failed a mandatory test
		//1 = passed but continue other checks to see if any other mandatory check are failing
		//2 = check passed return true (break)
		for(StringFilter check : checks)
		{
			include=(boolean) check.isPass(fileName);
			if(include ==false)
				return false;
		}
		return true;
	}

	

	public ArrayList<String> searchDirectory(File directory, String searchString, String forbiddenString)
			throws IOException {
		File[] files = directory.listFiles();
		ArrayList<String> fastqFiles = new ArrayList<String>();
		for (File file : files) {
			if (file.isDirectory())
				searchDirectory(file, searchString, forbiddenString);
			else {
				if (file.getName().contains(searchString) && !file.getName().contains(forbiddenString)) {
					fastqFiles.add(file.getAbsolutePath());
				}
			}
		}
		return fastqFiles;
	}

	private String parseString(String value) {
		if (value.equals("null"))
			return null;
		return value;
	}

	public String getWriteName() {
		return this.writeName;
	}

	public void setWriteName(String writeName) {
		this.writeName = writeName;
	}

	public String getFolders() {
		return folders;
	}

	public void setFolders(String folders) {
		this.folders = folders;
	}

	public String[] getSearchStrings() {
		return searchStrings;
	}

	public void setSearchStrings(String[] searchStrings) {
		this.searchStrings = searchStrings;
	}

	public String getRequiredStringFN() {
		return requiredStringFN;
	}

	public void setRequiredStringFN(String requiredStringFN) {
		this.requiredStringFN = requiredStringFN;
	}

	public String[] getForbiddenStrings() {
		return forbiddenStrings;
	}

	public void setForbiddenStrings(String[] forbiddenStrings) {
		this.forbiddenStrings = forbiddenStrings;
	}

	public void checkArgs(String[] args) {
		if (System.getProperty("user.dir").contains("C:\\Users\\Sipko\\git\\PCrotation\\Sipko"))
			return;
		System.out.println("This script takes the following arguments:\n"
				+ "1. folderName=<foldername>, name of the folder to search through\n"
				+ "2. searchStrings=<searchString,searchString2>, a comma separated list of strings for which to include the files\n"
				+ "3. forbiddenStrings=<forbiddenString,forbiddenString2>, a comma separated list of strings that are not allowed to be in the filenames (default=.md5)\n"
				+ "4. writeFN=<writeFN>, name of the file to write (defaul=<searchString>.txt)\n");
		System.exit(1);
	}
}
