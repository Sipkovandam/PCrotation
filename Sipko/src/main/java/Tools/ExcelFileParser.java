package Tools;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
 
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.CellType;
import org.apache.poi.ss.usermodel.CreationHelper;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

public class ExcelFileParser
{
	//I never tested reading an excel file, only writing...
	
	Iterator<Row> iterator = null;
	XSSFWorkbook workbook = null;
	FileInputStream inputStream = null;
	HashMap<String,Integer> rowCount = new HashMap<String,Integer>();
	HashMap<String, Sheet> sheets = new HashMap<String, Sheet>();
	String writeFn= null;
	CellStyle dateStyle=null;
	
	public ExcelFileParser getReader(String fileName, int sheetNumber) throws IOException {
		String excelFilePath = fileName;
        inputStream = new FileInputStream(new File(excelFilePath));
		this.workbook = new XSSFWorkbook(this.inputStream);
            
        Sheet sheetToRead = workbook.getSheetAt(sheetNumber);
        iterator = sheetToRead.iterator();
        
        return this;
	}
	
	public ArrayList<String> getSheetNames(String fn) throws IOException
	{
		inputStream = new FileInputStream(new File(fn));
		
		if(workbook==null)
			workbook = new XSSFWorkbook(inputStream);    
		
		ArrayList<String> sheetNames = new ArrayList<String>();
		for (int i=0; i<workbook.getNumberOfSheets(); i++) {
		    sheetNames.add( workbook.getSheetName(i) );
		}
		return sheetNames;
	}
	
	public void closeReader() throws IOException
	{
		if(workbook!=null)
			workbook.close();
		if(inputStream!=null)
			inputStream.close();
	}
	public String readLine() throws IOException
	{
		return readLine(-1);
	}
	
	public String readLine(int totalCols) throws IOException {
        if (iterator.hasNext()) 
        {
            Row nextRow = iterator.next();
            StringBuilder line = new StringBuilder(); 
            if(nextRow.getLastCellNum()>0)
            	line = addCell(line, nextRow.getCell(0));
            
            if(totalCols==-1)
            	totalCols=nextRow.getLastCellNum();
            for(int cn=1; cn<totalCols; cn++) 
            {
            	
            	Cell cell = nextRow.getCell(cn);
            	line.append("\t");
        		line = addCell(line, cell);
            }
            return line.toString(); 
        }
        return null;       
    }

	private StringBuilder addCell(StringBuilder line, Cell cell)
	{
		if(cell==null)
			return line;
        
         switch (cell.getCellTypeEnum()) 
         {
             case STRING:
            	 line.append(cell.getStringCellValue());
                 break;
             case BOOLEAN:
            	 line.append(cell.getBooleanCellValue());
                 break;
             case NUMERIC:
            	 line.append(cell.getNumericCellValue());
                 break;
			 default:
				break;
         }
         return line;
    }
	
	public void createWriter(String writeFn) throws IOException
	{
		this.writeFn=writeFn;
		workbook=new XSSFWorkbook();  
	}
	
	public Sheet createSheet(String sheetName)
	{
		//if sheet already exists;
		if(sheets.containsKey(sheetName))
			return sheets.get(sheetName);
		//create sheet
		Sheet sheet = workbook.createSheet(sheetName);

		sheets.put(sheetName, sheet);
		rowCount.put(sheetName, 0);
		return sheet;
	}
	
	public void closeWriter() throws IOException
	{

		FileOutputStream outputStream =new FileOutputStream(this.writeFn); 
		if(this.workbook!=null)
			this.workbook.write(outputStream);
		this.workbook.close();
		if(outputStream!=null)
			outputStream.close();
	}
	
	public void write(String writeLine, String sheetName)
	{
		write(writeLine.split("\t"), sheetName);
	}
	
	public void write(String[] fields, String sheetName)
	{
		
		Sheet sheet = sheets.get(sheetName);
		if (sheet == null)
			 sheet = createSheet(sheetName);
		
		int currentRow = this.rowCount.get(sheetName);
		Row row = sheet.createRow(currentRow);
		incrementRowCount(currentRow+1,sheetName);
		
		int columnCount = 0;
	     
	    for (String field : fields) {//I think I can just write everything as string, but i ll find out...
	        Cell cell = row.createCell(columnCount);
	        columnCount++;
	        if (field instanceof String) {
	        	String stringValue=(String) field;
	        	if(stringValue.trim().length()==0)
	        		continue;
	            cell.setCellValue(stringValue);
//	            if(field.matches("..-..-...."))
//	            {
//	            	if(this.dateStyle==null)//"dd-MM-yyyy"
//	            		createStyle("yyyy/MM/dd");
//	            	cell.setCellStyle(dateStyle);
//	            }
	        }
	    }
	    //sheet.getRow(rownum)
	}

	private void createStyle(String style)
	{
		this.dateStyle = this.workbook.createCellStyle();
		CreationHelper createHelper = workbook.getCreationHelper();
		this.dateStyle.setDataFormat(
		    createHelper.createDataFormat().getFormat(style));
	}

	private void incrementRowCount(	int currentRow,
									String sheetName)
	{
		this.rowCount.put(sheetName, currentRow);
	}
	
	
    
    
	
}
