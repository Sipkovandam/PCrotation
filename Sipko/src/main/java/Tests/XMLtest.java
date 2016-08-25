package Tests;

import org.w3c.dom.*;
import org.xml.sax.SAXException;

import javax.xml.parsers.*;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import java.io.*;

public class XMLtest 
{
	static String writeFN = "E:\\Groningen\\Test\\XML\\test.xml";
	
	public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException, TransformerException 
	{	
//		writeXML();
				
		File inputFile = new File(writeFN);
		
		inputFile = new File("E:/Groningen/Data/PublicSamples/100SamplesTest/Rsample/TESTexpression/config.xml");
		
        DocumentBuilderFactory dbFactory = 
           DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        Document doc = dBuilder.parse(inputFile);
        
        Element config = doc.getDocumentElement();
        
        System.out.println(config.getNodeName());
        
        System.out.println("res=" + config.getElementsByTagName("writefolder").item(0).getTextContent());
        
//        System.out.println(config.getNodeName());
//        //doc.getDocumentElement().normalize();
//        System.out.print("Root element: ");
//        System.out.println(doc.getDocumentElement().getNodeName());
//        
//        System.out.println("----------------------------");
//        Element test = doc.createElement("test");
//        Text content = doc.createTextNode("Content of text");
//        test.appendChild(content);
//        doc.getElementsByTagName("config").item(0).appendChild(test);
//      
//        TransformerFactory tff = TransformerFactory.newInstance();
//        Transformer tf = tff.newTransformer();
//        DOMSource source = new DOMSource(doc);
//        
//        String FN2= writeFN.replace(".xml", "_2.xml");
//        System.out.println(FN2);
//        tf.setOutputProperty(OutputKeys.INDENT, "yes");
//		tf.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
//        StreamResult result = new StreamResult(new File(FN2).getAbsolutePath());
//        tf.transform(source,  result);

 
//        for (int temp = 0; temp < nList.getLength(); temp++) {
//           Node nNode = nList.item(temp);
//           System.out.println("\nCurrent Element :");
//           System.out.print(nNode.getNodeName());
//           if (nNode.getNodeType() == Node.ELEMENT_NODE) {
//              Element eElement = (Element) nNode;
//              System.out.print("company : ");
//              System.out.println(eElement.getAttribute("company"));
//              NodeList carNameList = 
//                 eElement.getElementsByTagName("carname");
//              for (int count = 0; 
//                 count < carNameList.getLength(); count++) {	 
//                 Node node1 = carNameList.item(count);
//                 if (node1.getNodeType() ==
//                    node1.ELEMENT_NODE) {
//                    Element car = (Element) node1;
//                    System.out.print("car name : ");
//                    System.out.println(car.getTextContent());
//                    System.out.print("car type : ");
//                    System.out.println(car.getAttribute("type"));
//                 }
//              }
//           }
//        }
		System.out.println("Done!");
	}

	private static void writeXML() throws TransformerException, ParserConfigurationException {
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder builder = factory.newDocumentBuilder();
		
		Document doc = builder.newDocument();
		Element config = doc.createElement("config");
		doc.appendChild(config);
		
		Element ele = doc.createElement("wdir");	
		config.appendChild(ele);
		ele.appendChild(
				doc.createTextNode("/Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/"));
		
		Element ele2 = doc.createElement("wdir2");	
		config.appendChild(ele2);
		ele2.appendChild(
				doc.createTextNode("/Volumes/Promise_RAID/GeneNetwork/Sipko/06-2016/Juha_QuantNorm_Covariance/Prediction/2"));
		
		
		//writing it to a file
		TransformerFactory tff = TransformerFactory.newInstance();
		Transformer tf = tff.newTransformer();
		DOMSource source = new DOMSource(doc);
		StreamResult result = new StreamResult(new File(writeFN).getAbsolutePath());
		tf.setOutputProperty(OutputKeys.INDENT, "yes");
		tf.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
		tf.transform(source, result);
	}

}
