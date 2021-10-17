### Subtask 1: Find areas where the elevation is below the flood height level of 3.7m
## Functional Unit 1
# 1) Import QGIS frameworks
from qgis.core import *
from qgis.gui import *
import processing
import os

# 2) Set the filepath location which stores the sample data required for running the application
# The end user may have to change this depending on the location where they store the sample data
filepath = "/Users/Karen Yuenying Leung/Documents/RMIT/GEOM2157/Major Project/Data/"

# 3) Declare a variable for each input file (i.e., sample data)
# Digital elevation model (DEM) for Southeast Queensland (QLD):
seqDemFile = "DP_SEQ_DEM_25M_100K/DP_SEQ_DEM_25M_100K/seq/w001001.adf"
# QLD suburb boundaries:
qldSuburbFile = "QLD_Locality_Boundaries/QSC_Extracted_Data_20210927_151039049000-33408/Locality_Boundaries.shp"
# QLD roads & tracks:
qldRoadFile = "QLD_Roads_Tracks/QSC_Extracted_Data_20210927_153310676000-29952/Baseline_roads_and_tracks.shp"
# Brisbane (BNE) bus stops - this is a .csv file and it needs to be converted into a point .shp file later
bneBusTable = "Bus_Stop_locations.csv"

# 4) Load input layers into QGIS
seqDemLayer = iface.addRasterLayer((filepath + seqDemFile), seqDemFile[:-4])
qldSuburbLayer = iface.addVectorLayer((filepath + qldSuburbFile), qldSuburbFile[:4], "ogr")
qldRoadLayer = iface.addVectorLayer((filepath + qldRoadFile), qldRoadFile[:4], "ogr")

# 5) Convert bus stop .csv file into a .shp file & add it to interface
# Codes below adapted from sample codes on howtoinqgis.wordpress.com & lecture notes
# Set the names for fields storing longitudes and latitudes
long_field = "LONGITUDE"
lat_field = "LATITUDE"
crs = 28356 # BNE is in GDA94/MGA Zone 56
# Declare a variable for the output .shp file & set a filepath for it
outputBneBusFileName = "Converted_Bus_Stop_locations.shp"
outputBneBusFile = filepath + outputBneBusFileName
# Define CRS for the shapefile
spatRef = QgsCoordinateReferenceSystem(crs, QgsCoordinateReferenceSystem.EpsgCrsId)
# Set the input csv file (i.e. table) as data provider of the new vector layer
inpTable = QgsVectorLayer((filepath + bneBusTable), bneBusTable[:4], "ogr")
prov = inpTable.dataProvider()
fields = inpTable.fields() # return fields of the table
# Generate a temporal point layer and add new attributes (ie fields in csv file) to it
tmpLayer = QgsVectorLayer("Point", "temp", "memory")
tmpLayerData = tmpLayer.dataProvider()
tmpLayerData.addAttributes(fields)
tmpLayer.updateFields()
# Assign XY coordinates to the points
pt = QgsPointXY()
outFeature = QgsFeature()
for feat in inpTable.getFeatures():
    attrs = feat.attributes()
    pt.setX(float(feat[long_field])) # X for longitude
    pt.setY(float(feat[lat_field]))  # Y for latitude
    outFeature.setAttributes(attrs)  # Set attributes for each feature
    outFeature.setGeometry(QgsGeometry.fromPointXY(pt)) # Set geometry of features - point features
    tmpLayerData.addFeature(outFeature) # Add the output point features defined above to the temporal layer
# Write info above to a vector layer and create this layer
QgsVectorFileWriter.writeAsVectorFormat(tmpLayer, outputBneBusFile, "utf-8", spatRef, "ESRI Shapefile" )
del tmpLayer

# 6) Load the newly created bus stop shapefile to QGIS
bneBusLayer = iface.addVectorLayer(outputBneBusFile, "Converted_Bus_Stop_locations", "ogr")

# 7) Reproject the DEM layer
# No need to reproject the other layers since the original files are already in EPSG:28356 
# Firstly create a dictionary to store parameters required for the projection tool
reprojRasFile = "DP_SEQ_DEM_25M_100K/DP_SEQ_DEM_25M_100K/SEQ_DEM.tif" # declare a variable for the output reprojected raster
reprojParaDict = {"INPUT": seqDemLayer, "TARGET_CRS": "EPSG:28356", "RESAMPLING": 0, "DATA_TYPE": 0, "MULTITHREADING": False, "OUTPUT": (filepath + reprojRasFile)}
processing.run("gdal:warpreproject", reprojParaDict)
# 8) Add the reprojected raster layer into QGIS
reprojSeqDemLayer = iface.addRasterLayer((filepath + reprojRasFile), reprojRasFile[:-4])


## Functional Unit 2
# 1) Create a new layer that only contains suburbs in BNE LGA
# Select suburbs in BNE LGA from the existing QLD suburb layer
qldSuburbLayer.selectByExpression("\"LGA\" = 'Brisbane City'") 
# Write info of selected features into the new vector layer
bneSuburbFileName = "QLD_Locality_Boundaries/QSC_Extracted_Data_20210927_151039049000-33408/BNE_Suburbs.shp" # declare a variable for the output layer (path)
bneSuburbFile = filepath + bneSuburbFileName
writer = QgsVectorFileWriter.writeAsVectorFormat(qldSuburbLayer, bneSuburbFile, "utf-8", spatRef, "ESRI Shapefile", onlySelected = True) 
# load the newly created layer into QGIS
bneSuburbLayer = iface.addVectorLayer(bneSuburbFile, bneSuburbFile[:4], "ogr") # load the newly created layer into QGIS
del (writer)

# 2) Create a new layer that only shows DEM of BNE
# Clip the SEQ DEM with bneSuburbLayer 
# Create a dictionary to store the parameters required for the clip tool
bneDemFile = "DP_SEQ_DEM_25M_100K/DP_SEQ_DEM_25M_100K/BNE_DEM.tif" # Declare a ariable for the output clipped DEM layer
clipDemParaDict = {"INPUT": reprojSeqDemLayer, "MASK": bneSuburbLayer, "TARGET_CRS": "EPSG:28356", "ALPHA_BAND": False, "CROP_TO_CUTLINE": True, "KEEP_RESOLUTION": True, "SET_RESOLUTION": False, "MULTITHREADING": False, "DATA_TYPE": 0, "OUTPUT": (filepath + bneDemFile)}
processing.run("gdal:cliprasterbymasklayer", clipDemParaDict)
# Add the clipped DEM layer into QGIS
bneDemLayer = iface.addRasterLayer((filepath + bneDemFile), "BNE_DEM")

# 3) Create a new layer that only contains roads passing through BNE suburbs
# Clip the QLD road layer with bneSuburbLayer
# Create a dictionary to store the parameters required for the clip tool
bneRoadFile = "QLD_Roads_Tracks/QSC_Extracted_Data_20210927_153310676000-29952/BNE_roads.shp" # Declare a ariable for the output clipped roads layer
clipRoadParaDict = {"INPUT": qldRoadLayer, "OVERLAY": bneSuburbLayer, "OUTPUT": (filepath + bneRoadFile)}
processing.run("native:clip", clipRoadParaDict)
# Add the clipped roads layer into QGIS
bneRoadLayer = iface.addVectorLayer((filepath + bneRoadFile), "BNE_roads", "ogr")


## Functional Unit 3
# 1) Find areas where the elevation is below 3.7m (i.e., flood-prone areas)
# Use raster calculator to identify pixels < 3.7m
# Create a dictionary to store the parameters required for the raster calculator
floodDemFile = "DP_SEQ_DEM_25M_100K/DP_SEQ_DEM_25M_100K/flood_DEM.tif" # Declare a ariable for the output flood layer
rasCalParaDict = {"EXPRESSION": "\"BNE_DEM@1\"< 3.7", "LAYERS": bneDemLayer, "CRS": "EPSG:28356", "OUTPUT": (filepath + floodDemFile)}
processing.run("qgis:rastercalculator", rasCalParaDict)
# Add the flood layer into QGIS
floodDemLayer = iface.addRasterLayer((filepath + floodDemFile), "flood_DEM")

# 2) Extract flood-prone areas from the layer above
# The layer above includes both flood-prone (band 1) and non flood areas (band 0)
# Therefore need to extract the flood-prone ones and create a layer containing these areas only
# same method as above
fldOnlyDemFile = "DP_SEQ_DEM_25M_100K/DP_SEQ_DEM_25M_100K/flood_only_DEM.tif" # Declare a ariable for the output flood layer
rasCalParaDict2 = {"EXPRESSION": "\"flood_DEM@1\" / ( \"flood_DEM@1\" > 0 )", "LAYERS": floodDemLayer, "CRS": "EPSG:28356", "OUTPUT": (filepath + fldOnlyDemFile)}
processing.run("qgis:rastercalculator", rasCalParaDict2)
# Add the flood-prone area layer into QGIS
fldOnlyDemLayer = iface.addRasterLayer((filepath + fldOnlyDemFile), "flood_only_DEM")

# 3) Convert the raster layer for the flood-prone areas to vector layer
# This step makes Subtasks 2 and 3 easier
# Create a dictionary to store the parameters required for the Polygonize tool
fldOnlyVecFile = "DP_SEQ_DEM_25M_100K/DP_SEQ_DEM_25M_100K/flood_only_vec.shp" # Declare a ariable for the output vector layer
rasToVecDict = {"INPUT": fldOnlyDemLayer, "BAND": 1, "FIELD": "DN", "EIGHT_CONNECTEDNESS": False, "OUTPUT": (filepath + fldOnlyVecFile)}
processing.run("gdal:polygonize", rasToVecDict)
# Add the vector layer into QGIS
fldOnlyVecLayer = iface.addVectorLayer((filepath + fldOnlyVecFile), "flood_only_vec", "ogr")

# 4) Fix the invalid geometries in the layer above
# processing tools would stop running when there are invalid geometries
fixGeoFldOnlyVecFile = "DP_SEQ_DEM_25M_100K/DP_SEQ_DEM_25M_100K/fixGeo_flood_only_vec.shp"
fixGeoDict = {"INPUT": fldOnlyVecLayer, "OUTPUT": (filepath + fixGeoFldOnlyVecFile)}
processing.run("native:fixgeometries", fixGeoDict)
fixGeoFldOnlyVecLayer = iface.addVectorLayer((filepath + fixGeoFldOnlyVecFile), "fixGeo_flood_only_vec", "ogr")

# 5) Dissolve all features in the layer above into one single feature (i.e., polygon)
# The attribute table of the vector layer converted from the raster layer has a lot of rows/features, corresponding to the pixels in the original raster layer
# having lots of features in the layer makes it time-consuming to check point(bus stop)-polygon(flood-prone area) intersection
# Thus, much more efficient to make it a single polygon (as these flood features all share the same value - elevation <3.7m)
# Create a dict to store parameters required for the Dissolve tool
dsvFldOnlyVecFile = "DP_SEQ_DEM_25M_100K/DP_SEQ_DEM_25M_100K/dsv_flood_only.shp"
dissolveDict = {"INPUT": fixGeoFldOnlyVecLayer, "FIELD": "DN", "OUTPUT": (filepath + dsvFldOnlyVecFile)}
processing.run("native:dissolve", dissolveDict)
# Add the dissolved layer into QGIS
dsvFldOnlyVecLayer = iface.addVectorLayer((filepath + dsvFldOnlyVecFile), "dsv_flood_only", "ogr")



## Functional Unit 4
# 1) count # of flooded pixelx in each suburb polygon
# Create a dict to store parameters required for the zonal histogram tool
pixInSubFile = "QLD_Locality_Boundaries/QSC_Extracted_Data_20210927_151039049000-33408/fldpix_in_suburb.shp"
zonalParaDict = {"INPUT_RASTER": fldOnlyDemLayer, "RASTER_BAND": 1, "INPUT_VECTOR": bneSuburbLayer, "OUTPUT": (filepath + pixInSubFile)}
processing.run("native:zonalhistogram", zonalParaDict)
# Add the resulting layer into QGIS
# HISTO_1 in the attribute table represents # of flooded pixels, while HISTO_NODA includes nodata and non-flooded pixels within a suburb
pixInSubLayer = iface.addVectorLayer((filepath + pixInSubFile), "fldpix_in_suburb", "ogr")


## Functional Unit 5
# 1) Add a new field containing % of flooded area in each suburb to the fld_in_suburb layer
pixInSubLayer.startEditing() # firstly make the attribute table editable
pIS_data_provider = pixInSubLayer.dataProvider() # data provider is the layer itself
pIS_data_provider.addAttributes([QgsField("PercFlood", QVariant.Double)]) # Specify the data type of the new field - %
pixInSubLayer.updateFields() 
pixInSubLayer.commitChanges() # Save all changes and stop editing > new field will appear in the table


## Functional Unit 6
# Calculate the % flooded areas within a suburb and assign the value to the newly created field
# 1)Get all the features from the layer
suburbs = pixInSubLayer.getFeatures()

# 2) then create a for loop to loop through all these suburbs and apply the formula below
# % flooded areas = (# flooded pixels)/(# flooded pixels + # other pixels)
for aSbrb in suburbs:
    fldAreaPercentage = (aSbrb["HISTO_1"])/(aSbrb["HISTO_1"] + aSbrb["HISTO_NODA"])
    pixInSubLayer.startEditing() # Start editing the table
    aSbrb["PercFlood"] = fldAreaPercentage # Write the value into the field
    pixInSubLayer.updateFeature(aSbrb) # Update features in the attribute table
pixInSubLayer.commitChanges() # save all changes and stop editing


### Subtask 2: Identify bus stops that are located within the flood-prone areas.
## Functional Unit 7
# Add a new field to the attribute table of the bus stop layer
# this tells whether a bus stop is at risk of being inundated - namely located in the flood-prone areas
# same procedures as those in Functional Unit 5
bneBusLayer.startEditing()
bneBus_data_provider = bneBusLayer.dataProvider()
bneBus_data_provider.addAttributes([QgsField("FloodRisk", QVariant.String)]) # value: Yes/No(= at risk or not)
bneBusLayer.updateFields() 
bneBusLayer.commitChanges()


## Functional Unit 8
# 1) Use a for loop to evaluate if a bus stop is in a flood-prone area
# This is done by checking if the a bus stop point intersects with the flood-prone polygon
# If contains > assign Yes to the FloodRisk field. If not > No to the field
busStops = bneBusLayer.getFeatures()
for aBusStop in busStops:
    floodAreas = dsvFldOnlyVecLayer.getFeatures() # after retrieving features in the point layer, also get those in the polygon layer
    for aFloodArea in floodAreas:
        if aBusStop.geometry().intersects(aFloodArea.geometry()):
            risk = "Yes"
        else:
            risk = "No"
    
#2) Write the classification into the new field
# same procedure as those in Functional Unit 6
    bneBusLayer.startEditing()
    aBusStop["FloodRisk"] = risk
    bneBusLayer.updateFeature(aBusStop)
bneBusLayer.commitChanges()
    

### Subtask 3: Identify roads that pass through the flood-prone areas
## Functional Unit 9
# 1) add a new field to the attribute table of the bne road layer
# This field gives info on whether the road is at risk of being flooded
# same procedure as creating a new field for the bus stop layer
bneRoadLayer.startEditing()
bneRoad_data_provider = bneRoadLayer.dataProvider()
bneRoad_data_provider.addAttributes([QgsField("FloodRisk", QVariant.String)]) # value: Yes/No(= at risk or not)
bneRoadLayer.updateFields() 
bneRoadLayer.commitChanges()

## Functional Unit 10
# 1) Merge road segments with the same name into one road
# the original file separates a road into multiple segments. If not combined, a road may appear multiple times in the .txt file generated in the last functional unit
# same procedure as dissolving flood-prone areas
dsvBneRoadFile = "QLD_Roads_Tracks/QSC_Extracted_Data_20210927_153310676000-29952/dsv_BNE_roads.shp"
dissolveDict2 = {"INPUT": bneRoadLayer, "FIELD": "STREET", "OUTPUT": (filepath + dsvBneRoadFile)}
processing.run("native:dissolve", dissolveDict2)
dsvBneRoadLayer = iface.addVectorLayer((filepath + dsvBneRoadFile), "dsv_BNE_roads", "ogr")

# 2) Create a for loop to check if a road passes through the flood-prone area
# check: if a road (line) intersects with the area (polygone) - if yes, assign "Yes", if no, assign "No"
roads = dsvBneRoadLayer.getFeatures()
for aRoad in roads:
    floodAreas2 = dsvFldOnlyVecLayer.getFeatures()
    for aFloodArea2 in floodAreas2:
        if aRoad.geometry().intersects(aFloodArea2.geometry()):
            risk = "Yes"
        else:
            risk = "No"
    
#3) Write the classification into the new field
    dsvBneRoadLayer.startEditing()
    aRoad["FloodRisk"] = risk
    dsvBneRoadLayer.updateFeature(aRoad)
dsvBneRoadLayer.commitChanges()


## Functional Unit 11
# a) For flood-prone suburbs:
# 1) Create layers that only contains the flood-prone suburbs
pixInSubLayer.selectByExpression("\"PercFlood\" > 0") 
# Write info of selected features into the new vector layer
fldSuburbOnlyFileName = "QLD_Locality_Boundaries/QSC_Extracted_Data_20210927_151039049000-33408/flooded_suburbs_only.shp" # declare a variable for the output layer (path)
fldSuburbOnlyFile = filepath + fldSuburbOnlyFileName
writer_fldSub = QgsVectorFileWriter.writeAsVectorFormat(pixInSubLayer, fldSuburbOnlyFile, "utf-8", spatRef, "ESRI Shapefile", onlySelected = True) 
# load the newly created layer into QGIS
fldSuburbOnlyLayer = iface.addVectorLayer(fldSuburbOnlyFile, "flooded_suburbs_only", "ogr") # load the newly created layer into QGIS
del (writer_fldSub)

# 2) Import suburb names and percentage of flooded areas of each suburb into a .txt file
# Set up the filepath for the output .txt file
txtFldSuburbFile = "Output_txt_files/flooded_suburbs.txt"
output_txtFldSuburb = open((filepath + txtFldSuburbFile), "wb")
# Create a list to store suburb names and percentage of flooded areas
featSuburbs = []
# Create a for loop to retrieve all attributes of interest
for aFeatSuburb in fldSuburbOnlyLayer.getFeatures():
    msgout1 = '%s, %s\n' % (aFeatSuburb["LOCALITY"], round((aFeatSuburb["PercFlood"])*100, 2)) # convert float to percentage and round it to 2 decimal places 
    unicode_message1 = msgout1.encode('utf-8')
    featSuburbs.append(unicode_message1)
featSuburbs.sort() # sort suburb in an alphabetically order

# Write info into the .txt file
for suburbItem in featSuburbs:
    output_txtFldSuburb.write(suburbItem)
output_txtFldSuburb.close()


# b) For bus stops:
# same procedures as a)
# 1) Create layers that only contains bus stops located in flood prone areas
bneBusLayer.selectByExpression("\"FloodRisk\" = 'Yes'") 
fldBusFileName = "flooded_bus.shp"
fldBusFile = filepath + fldBusFileName
writer_fldBus = QgsVectorFileWriter.writeAsVectorFormat(bneBusLayer, fldBusFile, "utf-8", spatRef, "ESRI Shapefile", onlySelected = True) 
fldBusLayer = iface.addVectorLayer(fldBusFile, "flooded_bus", "ogr")
del (writer_fldBus)

# 2) Import information about those stops into a .txt file - stop name, the street and suburb it's in
txtFldBusFile = "Output_txt_files/flooded_bus.txt"
output_txtFldBus = open((filepath + txtFldBusFile), "wb")
featBuses = []
for aFeatBus in fldBusLayer.getFeatures():
    msgout2 = '%s, %s, %s\n' % (aFeatBus["DESCRIPTIO"], aFeatBus["STREETNAME"], aFeatBus["SUBURB"])
    unicode_message2 = msgout2.encode('utf-8')
    featBuses.append(unicode_message2)
featBuses.sort()
for busItem in featBuses:
    output_txtFldBus.write(busItem)
output_txtFldBus.close()

# c) For roads:
# same procedures as a) and b)
# 1) Create layers that only contains roads passing through flood prone areas
dsvBneRoadLayer.selectByExpression("\"FloodRisk\" = 'Yes'") 
fldRoadFileName = "flooded_roads.shp"
fldRoadFile = filepath + fldRoadFileName
writer_fldRoad = QgsVectorFileWriter.writeAsVectorFormat(dsvBneRoadLayer, fldRoadFile, "utf-8", spatRef, "ESRI Shapefile", onlySelected = True) 
fldRoadLayer = iface.addVectorLayer(fldRoadFile, "flooded_roads", "ogr")
del (writer_fldRoad)

# 2) Import names of the flooded roads into a .txt file
txtFldRoadFile = "Output_txt_files/flooded_roads.txt"
output_txtFldRoad = open((filepath + txtFldRoadFile), "wb")
featRoads = []
for aFeatRoad in fldRoadLayer.getFeatures():
    msgout3 = '%s\n' % (aFeatRoad["STREET"])
    unicode_message3 = msgout3.encode('utf-8')
    featRoads.append(unicode_message3)
featRoads.sort()
for roadItem in featRoads:
    output_txtFldRoad.write(roadItem)
output_txtFldRoad.close()



## Functional Unit 12
# This is a newly added functional unit in order to export the results into a map as mentioned in the proposal
# I missed writing this step into the implementation plan - now add it back to the script

# 1) Symbolise layers of a) all suburbs in BNE, b) flood-prone areas
# a) all suburbs in BNE - layer "fldpix_in_suburb" created in Functional Unit 4
allSuburb_style_renderer = pixInSubLayer.renderer()
allSuburb_symbol = allSuburb_style_renderer.symbol()
allSuburb_symbol.setColor(QColor.fromRgb(245, 255, 230,80)) # set the fill colour
allSuburb_symbol.symbolLayer(0).setStrokeColor(QColor.fromRgb(149,158,131,60)) # set the border colour. The last value is opacity.
allSuburb_symbol.symbolLayer(0).setStrokeWidth(0.45)
pixInSubLayer.triggerRepaint()

# b) flood-prone areas - layer "dsv_flood_only" created in Functional Unit 3
fldArea_style_renderer = dsvFldOnlyVecLayer.renderer()
fldArea_symbol = fldArea_style_renderer.symbol()
fldArea_symbol.setColor(QColor.fromRgb(153,204,255))
fldArea_symbol.symbolLayer(0).setStrokeColor(QColor.fromRgb(153,204,255))
dsvFldOnlyVecLayer.triggerRepaint()

# 2) Export BNE suburbs & flood-prone areas into a .png file
# This image shows which part of Brisbane LGA is at risk of being flooded.
# Codes below are adapted from opensourceoptions.com/blog/pyqgis-render-print-save-a-layer-as-an-image/
# Declare variable for the image size and background colour
image1 = QImage(QSize(800, 800), QImage.Format_ARGB32_Premultiplied)
colour1 = QColor(255, 255, 255, 255)
image1.fill (colour1.rgba())
# Set up a painter to paint the layer on a image
painter1 = QPainter()
painter1.begin(image1)
painter1.setRenderHint(QPainter.Antialiasing)
# Set up map setting & set up background colour
mapSetting1 = QgsMapSettings()
mapSetting1.setBackgroundColor(colour1)
# Retrieve the layers to be exported
# Note: make sure that the flood-prone area layer is on the top of the suburb layer in the layer panel. If not, manually move them on the panel.
layer1_1 = QgsProject.instance().mapLayersByName("fldpix_in_suburb")
layer1_2 = QgsProject.instance().mapLayersByName("dsv_flood_only")
mapSetting1.setLayers([layer1_1[0], layer1_2[0]])
# Set extent of the output png file 
rectangle1 = QgsRectangle(mapSetting1.fullExtent()) # export full extent of the layer.
rectangle1.scale(1.1) # increase rectangle size by factor of 1.1
mapSetting1.setExtent(rectangle1)
mapSetting1.setOutputSize(image1.size())
# Paint layers on the image file
render1 = QgsMapRendererCustomPainterJob(mapSetting1, painter1)
render1.start()
render1.waitForFinished()
painter1.end()
# Save the output file into the working folder
image1File = "/Output_img_files/all_suburbs_and_flood_area.png"
image1.save ((filepath + image1File))

### End of the script! ###
