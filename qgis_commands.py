import os
import sys
import qgis
from qgis.core import QgsApplication
from qgis.core import QgsCoordinateReferenceSystem
from qgis.core import QgsProcessing
from qgis.core import QgsProcessingAlgorithm
from qgis.core import QgsProcessingMultiStepFeedback
from qgis.core import QgsProcessing
from qgis.core import QgsProcessingParameterVectorDestination
import geopandas as gpd
#start QGIS instance without GUI
#Make sure the prefix is correct. Even though QGIS is in '/usr/share/qgis',
#the prefix needs to be '/usr' (assumes Debian OS)

QgsApplication.setPrefixPath('/usr', True)
myqgs = QgsApplication([], False)
myqgs.initQgis()

######### INITIALISE THE PROCESSING FRAMEWORK ################

# Append the path where processing plugin can be found (assumes Debian OS)
sys.path.append('/usr/share/qgis/python/plugins')

#import modules needed
import processing
from processing.core.Processing import Processing

#start the processing module
processing.core.Processing.Processing.initialize()

def voronoi_polygons(buffer,input_layer,output_layer,):
    """Creates a layer of polygons from a layer of points. The layer extent is divides in such a way 
    that everything inside a polygon is closer to its point than to any other of the points.
    See documentation https://docs.qgis.org/3.16/en/docs/user_manual/processing_algs/qgis/vectorgeometry.html?highlight=voronoi#voronoi-polygons

    Args:
        buffer (float): Percentage of the point layer extent that the new layer must cover outside the point layer extent
        input_layer (str): Path to input point layer
        output_layer (str): Path to the new layer

    Returns:
        
    """
    return  processing.run("qgis:voronoipolygons",{ 'BUFFER' : buffer, 'INPUT' : input_layer, 'OUTPUT' : output_layer })

def rasterize(extent,field,height,width,input_vector,output_raster,init=None,invert=False,nodata=-9999.0,options="",units=0,data_type=5,extra="",burn=0.0):
    """From vector to raster See documentation https://docs.qgis.org/3.16/en/docs/user_manual/processing_algs/gdal/vectorconversion.html#rasterize-vector-to-raster . 
    
        Derived from gdal rasterize utility. See docs https://gdal.org/programs/gdal_rasterize.html                 

    Args:           
        burn (float): Fixed value to burn into a band for all features. Default is 0.0
        extent (str): lonmin, lonmax, latmin, latmax [EPSG:code]. If the extent is not specified, the minimum extent that covers the selected reference layer(s) will be used.
        field (str): Name of the band to burn values into. This is, the name of the column of the shapefile whose values we want to have in the raster layer.
        height (int): Height of the raster layer in pixels if units=0.  Horizontal resolution in georeferenced units if units=1
        width (int): Width of the raster layer in pixels if units=0. Vertical resolution in georeferenced units if units=1
        input_vector (str): path to vector input layer          
        output_raster (str): Path to output raster layer (with file extension)_
        init (_type_, optional): Pre-initialize the output image bands with these values. However, it is not marked as the nodata value in the output file. 
            If only one value is given, the same value is used in all the bands.. Defaults to None.
        invert (bool, optional): Invert rasterization?  Defaults to False.
        nodata (float, optional): Assign a specified nodata value to output bands. Defaults to -9999.0.
        options (str, optional): _description_. Defaults to "".
        units (int, optional): 0 if pixels, 1 if georeferenced units. Defaults to 0.
        data_type (int, optional): Defines the format of the output raster file. 0 (byte), 1(Int16), 2 (Uint16), 3 (Uint32), 4 (Int32), 5 (Float32), 6 (Float64), 7 (Cint16), 8 (CInt32), 9 (CFloat32), 10 (CFloat64). Defaults to 5.
        extra (str, optional): _description_. Defaults to "".

    Returns:
        _type_: _description_
    """
    return processing.run("gdal:rasterize",{ 'BURN' : burn, 'DATA_TYPE' : data_type, 'EXTENT' : extent, 'EXTRA' : extra, 'FIELD' : field, 'HEIGHT' : height, 'INIT' : None, 'INPUT' : input_vector, 'INVERT' : False, 'NODATA' : nodata, 'OPTIONS' : options, 'OUTPUT' : output_raster, 'UNITS' : units, 'WIDTH' : width })
def dissolve_geometry(input_layer,output_layer):
    gdf=gpd.read_file(input_layer)
    gdf["column_for_dissolving"]="abc"
    gdf.to_file(input_layer)
    return processing.run("native:dissolve",{ 'FIELD' : ['abc'], 'INPUT' : input_layer, 'OUTPUT' : output_layer})

def merge_vector_layers(epsg_code,input_layers,output_layer):
    return processing.run("native:mergevectorlayers",{ 'CRS' : QgsCoordinateReferenceSystem(epsg_code), 'LAYERS' : input_layers , 'OUTPUT' : output_layer })

def intersection(input_file,input_fields,overlay_file,overlay_fields,output_file):
    return processing.run("qgis:intersection",{ 'INPUT' : input_file, 'INPUT_FIELDS' : input_fields, 'OUTPUT' : output_file, 'OVERLAY' : overlay_file, 'OVERLAY_FIELDS' : overlay_fields, 'OVERLAY_FIELDS_PREFIX' : 'overlay' })

def difference(input_file,overlay_file,output_file):
    return processing.run("qgis:difference",{"INPUT":input_file, "OUTPUT":output_file,"OVERLAY":overlay_file})
def reclass(input_file,rules_file,output_file):
    """Function to change categoric values of a raster file. It is intended to reduce the number of categories or
    apply modifications.       

    Args:
        input_file (str): Path to the input raster file
        rules_file (str): Path to the txt file with the reclassification rules. Best is to look at an example
        at /home/usuario/OneDrive/playas/InVEST/Habitat Quality/reclass_rules.txt'
        output_file (str): Path to the output raster file.

    Returns:
        
    """
    return processing.run("grass7:r.reclass",{"input":input_file,"rules":rules_file,"output":output_file})

def linealization(input_file,rules_file,output_file):
    processing.run("qgis:creategrid",{})
    outputs = {}
    outputs["reclassified"]=processing.run("grass7:r.reclass",{"input": input_file,"rules":rules_file, "output":"./temporary_output.tif"})
    outputs["thined"]=processing.run("grass7:r.thin", {"input":outputs["reclassified"]["output"],"output":"./temporary_output.tif"})
    outputs["linealized"]=processing.run("grass7:r.to.vect", {"input": outputs["thined"]["output"],"type":0,"column":"value","output":output_file})
    os.remove("./temporary_output.tif")
    os.remove("./temporary_output.tfw")
    os.remove("./temporary_output.tif.aux.xml")
    return 


def create_grid(input_file,type,hspacing,vspacing,crs,output_file,hoverlay=0,voverlay=0):
    gdf=gpd.read_file(input_file)
    
    minx=str(float(gdf.bounds["minx"]))
    miny=str(float(gdf.bounds["miny"]))
    maxx=str(float(gdf.bounds["maxx"]))
    maxy=str(float(gdf.bounds["maxy"]))

    extent=str(minx+","+maxx+","+miny+","+maxy+" ["+crs+"]")

    processing.run("qgis:creategrid", {"TYPE":type,"EXTENT":extent,"HSPACING":hspacing,"VSPACING":vspacing,"HOVERLAY":hoverlay,"VOVERLAY":voverlay,"CRS":crs,"OUTPUT":"delete.shp"})
    processing.run("qgis:clip",{"INPUT":"delete.shp","OVERLAY":input_file,"OUTPUT":output_file})
    os.remove("borrar.shp")
    return

def length_of_lines(polygon_file,lines_file,length_field,count_field,output):
    return processing.run("qgis:sumelinelengths",{"LINES":lines_file,"POLYGONS":polygon_file,"LEN_FIELD":length_field,"COUNT_FIELD":count_field,"OUTPUT":output})


def polygonize(input_raster,band,field,output_vector):
    return processing.run("qgis:polygonize",{"INPUT":input_raster, "BAND":band, "FIELD":field,"OUTPUT":output_vector})
if __name__ == "__main__":
    x=1
    #create voroni polygons form point data
    #voronoi_polygons(10,'/home/usuario/OneDrive/geo_data/POINTS_TA_AVG_1.5m.shp','/home/usuario/OneDrive/InVEST/Recreation/borrar_ordenador_nuevo.shp')
    #voronoi_polygons(10,'/home/usuario/OneDrive/geo_data/POINTS_HSOL_SUM_1.5m.shp','/home/usuario/OneDrive/InVEST/Recreation/VORONOI_HSOL_SUM_1.5m.shp')
   
    #rasterization
    #rasterize(extent="-10.0,-6.0,41.5,44 [EPSG:4258]",field="valor",height=500,width= 500, input_vector="/home/usuario/OneDrive/InVEST/Recreation/VORONOI_HSOL_SUM_1.5m.shp",output_raster="/home/usuario/OneDrive/InVEST/Recreation/borrar.tif")
    #rasterize(extent="-10.0,-6.0,41.5,44 [EPSG:4258]",field="valor",height=500,width= 500, input_vector="/home/usuario/OneDrive/InVEST/Recreation/VORONOI_TA_AVG_1.5m.shp",output_raster="/home/usuario/OneDrive/InVEST/Recreation/VORONOI_TA_AVG_1.5mheight500width500.tif")

    #dissolving
    #dissolve_geometry("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/galician_coastline_buffer10km_nosea.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/galician_coastline_buffer10km_nosea.shp")

    #merging
    # 

    #intersection no funciona, no es capaz de encontrar el algoritmo
    #intersection("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/provincias_costeras.shp",["COUNTRY"],'/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/galician_coastline_buffer10km.shp', ['fid','cat','coastline','label'],'/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/borrar.shp')

    #reclassification
    #reclass('/home/usuario/OneDrive/playas/InVEST/lulcRCP 8.52081-2100.tif','/home/usuario/OneDrive/playas/InVEST/Habitat Quality/reclass_rules.txt','/home/usuario/OneDrive/playas/InVEST/Habitat Quality/prueba_reclass.tif')

    #linealizado
    #linealization(input_file='/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/linea_de_costa.tif',rules_file='/home/usuario/OneDrive/InVEST/Recreation/reclass_rules_coastline.txt',output_file='/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/prueba_linealization.shp')
   
    #create a grid
    #create_grid(input_file="/home/usuario/OneDrive/recreation/qgis/dissolvedaoi.shp",type=4,hspacing=10000,vspacing=10000,crs="EPSG:3035",output_file="/home/usuario/OneDrive/recration/qgis/prueba_gridding.shp")

    #compute_length_of lines
    #length_of_lines("home/usuario/OneDrive/recreation/qgis/aoi_5km_grid.shp","/home/usuario/OneDrive/geo_data/caminos_naturales/all_trial_network.shp","trial_len","trial_n","/home/usuario/OneDrive/recreation/qgis/trial_length.shp")

    # #polygonize
    # polygonize(input_raster="/home/usuario/OneDrive/recreation/InVEST/CZ_2012_galicia_reclass.tif",band=1,field="LULC_codes",output_vector="/home/usuario/OneDrive/recreation/InVEST/borrar.shp")

    #try for difference function 
    difference(input_file='/home/usuario/OneDrive/recreation/InVEST/lulc/arable_land.shp',overlay_file='/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_4.5_2026-2045.shp',output_file="home/usuario/OneDrive/recreation/InVEST/scenario_data/prueba_difference.shp")