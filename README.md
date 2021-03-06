# GEOM2157 - Major Project
Demonstration site for Semester 2 2021 Geospatial Programming Course


This site stores files created for the major project. 

This project creates a python script to evaluate flood risks in Brisbane, Queensland. Any area that is below 3.7m will be considered as flood-prone.


- How to use the script?
    1. Download the script named 'project_code.py' from this repository 
    2. Open QGIS and open the Python Console (from Plugins) and Editor
    3. Open the script in the Editor and run! (Scan the entire codes before running. Remember to change the filepath of the input and output layers to where you store the sample data/your own data.)


- Files stored in this repository are...?
    1. test.docx: just a empty word document for testing Git
    2. project_code.py: the python script for automating the task
    3. README.md: the file you are now reading - a brief description of things in this repository
    4. Bus_Stop_locations.csv: a excel file storing information about bus stops in Brisbane, e.g., stop name, suburb, longitude & latitude etc.
    5. QLD_Locality_Boundaries: suburb boundaries of Queensland. Use the shapefile 'Locality_Boundaries.shp' to run the script
    6. DP_SEQ_DEM_25M_100K: digital elevation model (DEM) of Southeast Queensland (SEQ). Use the adf file 'w001001x.adf' to run the script
    7. QLD_Roads_Tracks: road networks in Queensland. Use the shapefile 'Baseline_roads_and_tracks.shp' to run the script
   
    ** Unfortunately both DP_SEQ_DEM_25M_100K and QLD_Roads_Tracks folders exceed 25MB that Github/Git does not allow me to upload them (or it simply uploads part of the folder into Github such as the DP_SEQ_DEM_25M_100K one which is actually imcomplete! w001001x.adf is not there because it is large.)
   
    ** To download the DEM and roads data, please go to 
          Queensland Spatial Catalogue (Qspatial) https://qldspatial.information.qld.gov.au/catalogue/custom/index.page > 
          
        (for DEM data) Type 'digital elevation model southeast queeeensland' in the search box and find 'Digital Elevation Model??? 25 metre ??? Southeast Queensland'.
    ![image](https://user-images.githubusercontent.com/91411718/137626620-77153347-9d47-48f0-bb29-25e6db7501f4.png)
          
        It should be in the first page of the search results.
        
        Click Download datasets > enter your email and it will be sent to your email minutes later.
          
           'w001001x.adf' is in the subfolder 'seq'.
          

           (for road data) Search 'Baseline Roads and Tracks' and find 'Baseline Roads and Tracks - Queensland'
          
     ![image](https://user-images.githubusercontent.com/91411718/137626784-e52d75e7-ceeb-4156-b3da-9644551fdc23.png)
          
           It should be in the first page as well.
           
           Download datasets > select file type and coordinate reference system and enter your email address. It will be sent to your email minutes later.
          
     ![image](https://user-images.githubusercontent.com/91411718/137626956-109df3a3-02b0-45f9-90e2-21aa92526366.png)
          
          'Baseline_roads_and_tracks.shp' is in the folder.
          
   
    8. Output_img_files: 'all_suburbs_and_flood_areas.png' is a sample output image (or call it a map?) generated by the script. It shows the areas that are prone to flooding (in blue) in Brisbane.
   
    9. Output_txt_files: all the three files are sample output text files generated by the script.

             a).'flooded_bus.txt': list bus stops that are located in the flood-prone areas (contents: bus stop name, road, suburb)                    
             b).'flooded_roads.txt': roads that pass through the flood-prone areas (contents: road name)                     
             c). 'flooded_suburbs.txt': suburbs that have parts (or the entire suburb) at risk of flooding, i.e., elevation < 3.7m (contents: suburb name, % of risky areas)



That's it! Enjoy exploring Brisbane :)

   
