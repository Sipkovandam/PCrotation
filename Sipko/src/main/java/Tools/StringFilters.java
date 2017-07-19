package Tools;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class StringFilters
{
	private List<StringFilter> filters = new ArrayList<>();
	
	public List<StringFilter> addFilter(StringFilter filter)
	{
		filters.add(filter);
		return this.filters;
	}
	
	public class RequiredStringChecker implements StringFilter {
		private String[] requiredString = null;
		
		public RequiredStringChecker(String[] requiredString)
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
	
	public class ForbiddenStringChecker implements StringFilter {
		private String[] forbiddenString = null;
		
		public ForbiddenStringChecker(String[] forbiddenString)
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
	
	public class SampleNameChecker implements StringFilter {
		private Set<String> includeList = null;
		private String[] replaceElements = null;
		
		public SampleNameChecker(Set<String> includeList,String[] replaceElements)
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
		
			if(replaceElements != null)
			{
				for (String replaceElement : replaceElements) {
					if(replaceElement!= null)
						fileName = fileName.replaceAll(replaceElement, "");
				}
			}
			
			if (includeList.contains(fileName))
			{
				return true;
			}
			return false;
		}
	}

	public List<StringFilter> getFilters()
	{
		return filters;
	}

	public void setFilters(List<StringFilter> filters)
	{
		this.filters = filters;
	}

	public boolean check(String string)
	{
		boolean include = true;
		for(StringFilter check : this.filters)
		{
			include=(boolean) check.isPass(string);
			if(include ==false)
				return false;
		}
		return true;
	}
}
