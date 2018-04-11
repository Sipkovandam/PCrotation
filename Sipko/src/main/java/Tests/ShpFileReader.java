package Tests;

import java.io.File;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.FeatureSource;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
//import org.geotools.data.shapefile.shp;
import org.opengis.feature.Feature;
import org.opengis.feature.GeometryAttribute;

import Tools.Script;

public class ShpFileReader extends Script<ShpFileReader>
{
	String fileName = "E:/Groningen/Data/Juha/31-01-2018/ICA/C++/testResults/test.sph";

	@Override
	public void run()
	{
		File file = new File(fileName);//"mayshapefile.shp");

		try {
		  Map<String, String> connect = new HashMap();
		  log("File:" + file);
		  log("Exists:" + file.exists());
		  connect.put("url", file.toURI().toString());
		  log("Exists1:\t " + file.toURI().toString());
		  
		  Iterator availableStores =  DataStoreFinder.getAvailableDataStores();

          while (availableStores.hasNext()) {
              log(availableStores.next().toString());
          }

		  
		  DataStore dataStore = DataStoreFinder.getDataStore(connect);
		  
		  //ShapefileDataStore dataStore = new ShapefileDataStore(file.toURI().toURL());
		  
		  String[] typeNames = dataStore.getTypeNames();
		  String typeName = typeNames[0];

		  System.out.println("Reading content " + typeName);

		  FeatureSource featureSource = dataStore.getFeatureSource(typeName);
		  FeatureCollection collection = featureSource.getFeatures();
		  FeatureIterator iterator = collection.features();


		  try {
		    while (iterator.hasNext()) {
		      Feature feature = iterator.next();
		      GeometryAttribute sourceGeometry = feature.getDefaultGeometryProperty();
		    }
		  } finally {
		    iterator.close();
		  }

		} catch (Throwable e) {e.printStackTrace();}
		
	}
}
