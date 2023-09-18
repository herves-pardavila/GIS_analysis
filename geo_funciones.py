import geopandas as gpd
import pandas as pd
import rasterio
from rasterio import features
from rasterio.enums import MergeAlg
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
import time
from shapely import geometry
from tqdm import tqdm
def crop_to_bounds(original,mask):
    """Function to crop the spatial extension of a layer to the extension of another

    It does not work as we would like

    Args:
        original (geopandas.GeoDataFrame): _Layer (read from shapefile) whose spatial extension we need to crop
        mask ( geopandas.GeoDataFrame): Other layer (read from shapefile ) whose spatial extension will be used to crop original

    Returns:
        geopandas.GeoDataFrame: The cropped layer
    """
    xmin, ymin, xmax, ymax = mask.total_bounds
    new=original.cx[xmin:xmax,ymin:ymax]
    return new
def transform_lulc_vector(code_table,columns_in_table,vector,column_in_vector,path,save=True):
    """This function groups several land use and land cover types in just one. In order to perform this
    grouping you need to create/edit and excel file with (at least) these columns: the lulc codes, the name of the habitat that each code
    represents, the new lucl codes and their new names.

    Args:
        code_table (str): Path to the excel file where each lucl type has its new equivalent.
        columns_in_table (list): list of names of relevant columns of the code_table: old lulc codes, old names
        new lulc codes, new names. 
        vector (geopandas.GeoDataFrame): vector data of the land use and land cover 
        column_in_vector (str): Name of the column in vector with the lulc codes.
        path (str): parth to the output file
        save(bool,optional) Save the result. Default True

    Returns:
        geopandas.GeoDataFrame: Vector with new columns yielding the new habitat classification. Old columns are not
        deleted so you can also acces the old classification.
    """
        
    df=pd.read_excel(code_table)
    df=df[columns_in_table]
    newvec=pd.merge(vector,df,how="inner",left_on=column_in_vector,right_on=columns_in_table[0])
    if save==True:
        newvec.to_file(path)
    return  newvec
def rasterize_vector_with_template(input_file,template,variable,output_CRS,formato,output_file):
    """This functions rasterizes a shapefile layer by using another existing
    raster as s template. You needt to verify (easiest way to do it is in QGIS or other
    GIS platform) that vector file and raster template have the same extension when you set
    the same crs for both

    Args:
        input_file (str): Path to vector file
        template (str): Path to template raster
        variable (str): Name of the column in the vector whose values we want to burn
        output_file (str): Path to the output raster file.

    Returns:
        numpy.ndarray: raster values in matrix format
    """
    raster=rasterio.open(template)
    shapefile=gpd.read_file(input_file)
    geom_value = ((geom,value) for geom, value in zip(shapefile.geometry, shapefile[variable]))
    #we create the array    
    rasterized = features.rasterize(geom_value,
                                out_shape = raster.shape,
                                transform = raster.transform,
                                all_touched = True,
                                fill = 0,   # background value for no data
                                merge_alg = MergeAlg.replace,
                                dtype = np.float64)
    
    #we save the array as a raster file.
    with rasterio.open(output_file, "w",
    driver = "GTiff",
    transform = raster.transform,
    crs=output_CRS,
    dtype = formato,
    count = 1,
    width = raster.width,
    height = raster.height) as dst:
        dst.write(rasterized, indexes = 1)
    
    return rasterized
def set_raster_projection(raster_file,prj):
    """Change the crs of a raster layer

    Args:
        raster_file (path): path to the layer_
        prj (str): new projection in EPSG code. Example EPSG:3035
    """
    raster=gdal.Open(raster_file, gdal.GA_Update)
    print("Antigua Proyección",raster.GetProjection())
    #nueva proyeccion
    raster.SetProjection(prj)
    print("Nueva Proyección",raster.GetProjection())
    del raster #cerrar y guardar
    return
def see_raster_projection(raster_file):
    """Prints the CRS of a raster layer

    Args:
        raster_file (str): Path to raster layer
    """
    raster=gdal.Open(raster_file, gdal.GA_Update)
    print(raster.GetProjection())
    del raster #cerrar y guardar
    return
def see_vector_projection(vector_file):
    """Prints CRS of a vector layer

    Args:
        vector_file (str): Path to the vector layer
    """
    gdf=gpd.read_file(vector_file)
    print(gdf.crs)
    del gdf #cerrar y guardar
    return
def visualize_raster(input_file,layer=1):
    """Plots raster values

    Args:
        input_file (str): Path to input raster
        output_file (str) Path to output file for saving an image
        layer (int, optional): Band. Defaults to 1.
        
    """
    f=rasterio.open(input_file)
    x=f.read(layer)
    fig=plt.figure()
    plt.imshow(x)
    plt.colorbar()
    plt.show(True)

    return 
def save_array_as_raster(array,input_template,output_file):
    """Save a numpy array as a tif raster using another 
    raster template. The new raster will have same crs, transformation
    and dimensions as template!. Be careful when chosing template!

    Args:
        array (numpy.ndarray): Raster values we want to save
        input_template (str): Path to template
        output_file (str): Output path
    """
    
    template=rasterio.open(input_template)
        
    with rasterio.open(output_file, "w",
    driver = "GTiff",
    transform = template.transform,
    crs=template.crs,
    dtype = np.int32,
    count = 1,
    width = template.width,
    height = template.height) as dst:
        dst.write(array, indexes = 1)
    return 
def compute_coast_line(dem_file,sea_value):
    """This function computes the coastline from a digital elevation
    model in raster format. The DEM needs to have some no-data value that 
    represents sea, or, some unique value for only the sea

    Args:
        dem_file (str): Path to DEM file
        sea_value (float): The value in the DEM raster used to represent pixels where there
        is sea

    Returns:
        numpy.ndarray: Matrix with 1's representing the coastline and 0 for everything else.
    """
    f=rasterio.open(dem_file) #load dem file
    dem=f.read(1) #create array 
    dem[dem!=sea_value]=1 #land pixels are all set to one
    dem[dem==sea_value]=0 #sea pixels are all set to zero
    for fila in range(1,np.shape(dem)[0]-1): 
        for columna in range(1,np.shape(dem)[1]-1):
            #we check, for the land values if the surrounding pizels are sea, the pixels that verify
            #such conditions are the coast
            arriba,abajo,izq,der=dem[fila-1,columna],dem[fila+1,columna],dem[fila,columna-1],dem[fila,columna+1]
            if dem[fila,columna]==1 and any(np.array([arriba,abajo,izq,der])==0):
                #if any of the surroinfing pixels is sea, then that pixel is coastline
                #print("Si")
                dem[fila,columna]=2
    dem[dem!=2]=0
    dem[dem==2]=1
    return dem
def light_path(noriginal, moriginal,n,m, nmax, mmax, dem_matrix,path_matrix,visibility_matrix,steps,step=1):
    steps=steps+step
    distance=np.sqrt((m-moriginal)**2+(n-noriginal)**2) 

    path_matrix[n][m] = distance # le pongo un 1 al nodo  por el que estoy caminando
    #miro si ya estoy en el mar
    if steps > 900 :
        pass
    elif distance > np.sqrt((7/3)*(6378*1000/25)*(dem_matrix[noriginal,moriginal]+1.6/25)):
        pass
    
    elif dem_matrix[n][m]==dem_matrix.min(): #he llegado al mar

        print("Ha llegado al mar")
        print(noriginal,moriginal)
        #fig=plt.figure()
        #plt.imshow(path_matrix)
        #plt.show(True) 
        #visibility_matrix[noriginal,moriginal]=1
           
    # miro si puedo caminar hacia abajo
    
    elif 0 <  n < nmax - 1 and 0 < m < mmax - 1 and dem_matrix[n+1][m] < dem_matrix[n][m]+1.6/25: #si no he llegado a los bordes y el nodo inferior tiene menor altura:
        light_path(noriginal,moriginal,n+1,m, nmax,mmax,dem_matrix,path_matrix,visibility_matrix,steps) #vuelvo a llamar la función pero en el nuevo punto al que voy
    
    elif n + 1 == nmax or m + 1 == mmax: #si llegué a los bordes, salgo de la función
        pass
    
    # si no puedo caminar hacia abajo y no terminé , pruebo a la izquierda
    
    elif m>0 and dem_matrix[n][m - 1] < dem_matrix[n][m]+1.6/25: #si el nodo  de mi izquierda está a menor altura vuelvo a empezar pero en el nuevo nodo 
        light_path(noriginal,moriginal,n, m-1,nmax,mmax, dem_matrix, path_matrix, visibility_matrix,steps)
        
        #pruebo a la derecha
    
    elif m < mmax - 1 and dem_matrix[n][m + 1] < dem_matrix[n,m]+1.6/25: 
        light_path(noriginal,moriginal,n, m+1,nmax,mmax ,dem_matrix,path_matrix,visibility_matrix,steps)
    
    elif n>0 and dem_matrix[n-1][m] < dem_matrix[n][m]+1.6/25: #a menos que esté en la primera fila me muevo hacia el nodo que tengo encima 
        light_path(noriginal,moriginal,n-1, m,nmax,mmax, dem_matrix, path_matrix,visibility_matrix,steps)

    return visibility_matrix, path_matrix,steps
def visibility_of_sea(dem_file,coast_file,dem_layer=1,coast_layer=1):
    """This function computes the areas from where the coast is visible. The function needs
    previous treatment on a Digital Elevation Model (DEM) file to get the regions of interest. The fucntion is time-consuming
    for a 25 m resolution raster for Galicia it lasts several days when the area being evaluated as coast was 
    a 2 - 10 km buffer from coastline. 

    Args:
        dem_file (path): Path to a DEM file of the region of interest. Elevation are assumed to be in meters and 
            to be 
        coastline_file (str): Path to raster file with the area that coud be coast. It must have 1 for pixels that 
            could be coast and 0 otherwise. It must have same extension, resolution and CRS as the dem_file.
        dem_layer (int, optional): Layer of the dem_file where elevation data is. Default to 1
        coast_layer (int, optional): Layer of the coast_file where possible coast pixels are. Default to 1

    Returns:
        numpy.ndarray : Array with same dimension as the raster files. 1 values are pixels where coastline is visible
        and 0 values pixels where is not.
    """

    #dem_file with elevarion data
    f_dem=rasterio.open(dem_file)
    dem=f_dem.read(dem_layer) #we assume 
    dem[:,8500:]=dem.min()
    #file with the region to perform the visibility analysis
    f_coast=rasterio.open(coast_file)
    coast=f_coast.read(coast_layer)
    coastal_values=np.argwhere(coast==1)
    
    #limito los valores de la costa para hacer una prueba
    coastal_values=coastal_values[:700000+2290000]
    
    #this is array is what function returns, lo hemos cambiado temporalmente para ahorarnos tiempo de computacion
    #visibility=np.zeros(np.shape(dem))
    f_visibility=rasterio.open("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/visibility.tif")
    visibility=f_visibility.read(1)
    print("Interactions=",len(coastal_values)) #total number of pixels to be analyzed
    
    Rt=6378000 #earth radius in meters
    
    #raster dimensions
    mmax=f_dem.width
    nmax=f_dem.height
    
    #we will print the maximum elevation in the studied area
    max_dem=0

    t=time.time()
    
    for i in range(len(coastal_values)):
         
        if i%1000==0: #time for 1000 iterations 
            print(i)
            print(time.time()-t)
            t=time.time()
        
        row=coastal_values[i][0]
        column=coastal_values[i][1]
        
        if dem[row,column] < -1.6: 
        # if elevation is lower than a persons height (m) horizont is negative. No sense to perform this calculations
            continue
        #here we save the maximun recorded elevation
        max_dem=max(max_dem,dem[row,column])
       
        
        #distance of the horizont (in pixels) as a function of observer's elevation.
        hmax1357=np.sqrt((7/3)*(Rt/25)*(dem[row,column]+1.6)/25)
        hmax1357=int(hmax1357)

        #visibility towards north
        
        #all way of doing it
        #hmax=max(0,row-hmax1357) #vemos si el horizonte de visión está mas lejos que el final de la capa
        #vision=dem[hmax:row,column]
        
        vision=dem[0:row,column]
        sea=int(max(np.argwhere(vision==dem.min()))) #closest sea pixel in north direction
        vision=vision[sea:]

        #height reduction due to earth's curvature and refraction in meters
        deltah=[((25*i)**2)/(2*Rt) for i in range(1,len(vision)+1)]
        deltah=np.array(deltah)

        #if obstacles are not as tall as observer  (meters) and the sea is before the horizont (in pixels) I can see 
        #the sea from tha pixel
        if (vision -deltah < dem[row,column]+1.6).all()==True and row-sea < hmax1357:
            visibility[row,column]=1
        else:
            #visibility towards east, analogous as towards north
            
            #oldway
            #hmax=min(mmax,column+hmax1357)
            #vision=dem[row,column:hmax]

            vision=dem[row,column:mmax]
            sea=int(min(np.argwhere(vision==dem.min())))
            vision=vision[:sea]
            
            deltah=[((25*i)**2)/(2*Rt) for i in range(1,len(vision)+1)]
            deltah=np.array(deltah)

            if (vision - deltah < dem[row,column]+1.6).all()==True and sea-column < hmax1357:
                visibility[row,column]=1
            else:
                #visibility towards east, analogous as towards west

                #oldway
                #hmax=max(0,column-hmax1357)
                #vision=dem[row,hmax:column]

                vision=dem[row,0:column]
                sea=int(max(np.argwhere(vision==dem.min())))
                vision=vision[sea:]
         
                deltah=[((25*i)**2)/(2*Rt) for i in range(1,len(vision)+1)]
                deltah=np.array(deltah)

                if (vision-deltah < dem[row,column]+1.6).all()==True and column-sea< hmax1357:
                    visibility[row,column]=1
                else:
                    #visibility towards south

                    #oldway
                    #hmax=min(nmax,row+hmax1357)
                    #vision=dem[row:hmax,column]
                    
                    vision=dem[row:nmax,column]
                    sea=int(min(np.argwhere(vision==dem.min())))
                    vision=vision[:sea]
                    
                    deltah=[((25*i)**2)/(2*Rt) for i in range(1,len(vision)+1)]
                    deltah=np.array(deltah)
                    if (vision - deltah < dem[row,column]+1.6).all()==True and sea-row <- hmax1357:
                        visibility[row,column]=1
      
        
           
    print("Altura máxima=",max_dem)

    return visibility      
def create_square_grid(input_file,size,output_file):

    """Given an area of interest it builds an square within its sizes

    Args:
        input_file (str): Path to the input area of interest. Must be shapefile
        size (float): Size of the squared cells in meters 
        output_file (str): Path to output file
    """
    gdf=gpd.read_file(input_file)
    # Get the extent of the shapefile
    total_bounds = gdf.total_bounds
    # Get minX, minY, maxX, maxY
    minX, minY, maxX, maxY = total_bounds
    # Create a fishnet
    x, y = (minX, minY)
    geom_array = []
    # Polygon Size
    square_size = size
    print("Entra en el bucle")
    while y <= maxY:
        while x <= maxX:
            geom = geometry.Polygon([(x,y), (x, y+square_size), (x+square_size, y+square_size), (x+square_size, y), (x, y)])
            geom_array.append(geom)
            x += square_size
        
        x = minX
        y += square_size
    fishnet = gpd.GeoDataFrame(geom_array, columns=['geometry']).set_crs(gdf.crs)
    fishnet=gpd.clip(fishnet,gdf)
    fishnet.to_crs("EPSG:3035",inplace=True)
    fishnet.to_file(output_file)


    return
def length_of_lines(lines_file,polygon_file,spatial_indicator,name,output_file,save=False):
    """Computes the length of a layer of lines in a layer of polygons. 

    Args:
        lines_file (str): Path to lines layer. Must be shapefile.
        polygon_file (str): Path to the polygons layer. Must be shapefile.
        spatial_indicator (str): Name of the column that identifies the different polygons in the polygon layer
        name (str): Name of the new column where information of length of lines will be stored.
        output_file (str): Path to the output file.

    Returns:
        grid (geopandas.GeoDataFrame): gridded area of interest
    """
    lines=gpd.read_file(lines_file)
    #we only need index and geometry
    lines["index"]=lines.index
    lines=lines[["geometry","index"]]
    #polygons
    if type(polygon_file)==str:
        grid=gpd.read_file(polygon_file)
    else:
        grid=polygon_file

    if grid.crs != lines.crs:
        lines.to_crs(grid.crs,inplace=True)

    overlay=gpd.overlay(lines,grid,how="intersection")
    overlay[name]=overlay.geometry.length
    for i in tqdm(range(len(overlay))):
        grid.loc[grid[spatial_indicator]==overlay.iloc[i][spatial_indicator],name]=overlay.iloc[i][name]
    grid[pd.isna(grid.length)==True][name]=0.0
    if save==True:
        grid.to_file(output_file,index=False)
    return grid
def distance_to_points(point_layer,other_layer,field_name,output_file,save=False):
    """Computes, for each polygon in a polygons layer, their distance to the closest point 
    of a layer of points    

    Args:
        point_layer (str): Path to a the point layers. Must be shapefile.    
        other_layer (str): Path to the polygons layer. Must be shapefile.
        field_name (str): Name of the new column that must be created in order to store the information.    
        output_file (str): _Path to the output file

    Returns:
        str: Layer with the new column
    """

    if type(other_layer)==str:
        layer=gpd.read_file(other_layer)
    else:
        layer=other_layer
    point_layer=gpd.read_file(point_layer)
    if point_layer.crs != layer.crs:
        point_layer.to_crs(layer.crs,inplace=True)
    distances=[point_layer.distance(layer.iloc[i]["geometry"]).min() for i in tqdm(range(len(layer)))]
    layer[field_name]=distances
    if save==True:
        layer.to_file(output_file,index=False)
    return layer
def polygonize(input_raster,output_vector=None,save=False):
    """From raster (tif) to vector (shapefile). It will automatically burn the first band. I did not 
    test it with rasters of more than one band.

    Args:
        input_raster (str): Path to the input raster
        output_vector (str, optional): Path to the output vector if you want to save.  Defaults to None.
        save (bool, optional): Save the file. Defaults to False.

    Returns:
        geopandas.GeoDataFrame : Polygon layer
    """
    mask = None
    with rasterio.Env():
        with rasterio.open(input_raster) as src:
            image = src.read(1) # first band
            results = ({'properties': {'raster_val': v}, 'geometry': s} for i, (s, v) in enumerate(features.shapes(image, mask=mask, transform=src.transform)))
            geoms = list(results)
            gdf  = gpd.GeoDataFrame.from_features(geoms,crs=src.crs)
            
    
    if save==True:
        gdf.to_file(output_vector)
    return gdf
def spatial_overlay(layer1,layer2,path_voronoi1,path_voronoi2,crs,spatial_indicator):
    """Given two polygon layers and their respective voronoi polygons this functions assigns to each polygon in the
    layer 1, its correspondent polygon in layer 2 with all its information. The criteria to choose which poltygon of layer 2 is assigned to each
    polygon of layer 1 is the spatial overlay, i.e the polygon of layer 2 whose intersection with the polygon of layer1 is the
    largest. All information of layer 2 will be included in the polygons of layer 1. 

    Args:
        layer1 (pandas.DataFrame): Information of Layer 1
        layer2 (pandas.DataFrame): Information of Layer 2
        path_voronoi1 (str): Path to the voronoi polygons of layer 1
        path_voronoi2 (_type_): Path to the voronoi polygons of layer 2
        crs (str): EPSG code. For example EPSG:3035
        spatial_indicator (str): The name of the column in layer 1 with the names of the spatial units

    Returns:
        pandas.DataFrame: layer 1 with the information of layer 2
    """

    #loading voronoi files with only the spatial information
    voronoi1=gpd.read_file(path_voronoi1)
    voronoi2=gpd.read_file(path_voronoi2)
    voronoi1.to_crs(crs,inplace=True)
    voronoi2.to_crs(crs,inplace=True)

    #spatial overlay
    overlay=gpd.overlay(voronoi1,voronoi2,how="intersection")
    overlay["area"]=overlay["geometry"].area/10**6
    #we sort this area
    overlay.sort_values(by=["area"],inplace=True)
    #and we keep the largest overlap for each voronoi polygon of the meteogalicia data
    overlay.drop_duplicates(subset=["idStation"],keep="last",inplace=True)
    overlay=overlay[[spatial_indicator,"lat_2","lon_2"]]
    #we start recovering the temporal information.
    layer1=layer1.merge(overlay,how="inner",on=spatial_indicator)

    #merge information from both sources
    layer1.time=pd.to_datetime(layer1.time).dt.strftime("%Y/%m")
    layer2.time=pd.to_datetime(layer2.time).dt.strftime("%Y/%m")
    layer1=layer1.merge(layer2,how="inner",left_on=["time","lat_2","lon_2"],right_on=["time","lat","lon"])

    return layer1
def extend_polygons_to_line(polygon_file,lines_file,buffer_distance,variable,output_path,save=False,plot=False):
    """This function is intended for the following purpose: Imagin you have a shapefile of a coastline (geometry are lines)
    and a shapefile of the sea level rise in some beaches (geoemtries are polygons). You want to to extrapolate the sea level 
    rise to all the coastline. Every piece of the coast is assigned a sea level rise value according to the closest beach. The 
    coastline shapefile has to be previously converted to polygons using a buffer. This buffer has to be as small as possible so 
    the geometries are as similar to a line as possible (5 m for example). Output is a shapefile of polygons that represents
    the coast with sea level rise values along all of it. The function is valid for any calculation that operates the same way,
    extending polygon values to a whole line. Polygons have to overlap with the line.
    Args:
        polygon_file (str): Path to th
        lines_file (_type_): _description_
        buffer_distance (_type_): _description_
        variable (_type_): _description_
        output_path (_type_): _description_
        save (bool, optional): _description_. Defaults to False.
        plot (bool, optional): _description_. Defaults to False.

    Returns:
        _type_: _description_
    """

    lines=gpd.read_file(lines_file)
    polygons=gpd.read_file(polygon_file)
    buffer=lines.buffer(5)
    lines["geometry"]=buffer.geometry
    lines=lines.sjoin(polygons,how="left")
    gdf1=lines[lines[variable].notna()]
    gdf2=lines[lines[variable].isna()]
    gdf1=gdf1.reset_index()
    gdf2=gdf2.reset_index()
    matrix=gdf2.geometry.apply(lambda g: gdf1.distance(g))
    for i in tqdm(range(len(gdf2))):
        columna=matrix.iloc[i].argmin()
        gdf2.loc[i,variable]=gdf1.iloc[columna][variable]
    joined=pd.concat([gdf1,gdf2])
    if save==True:
        joined.to_file(output_path)
    if plot==True:
        joined.plot(column=variable)
        plt.show()
    return joined

def point_count(polygon_layer,point_layer,spatial_indicator,column):

    sjoin=gpd.sjoin(polygon_layer,point_layer)
    pivot=pd.pivot_table(sjoin,index=spatial_indicator,columns=column,aggfunc={column:len})
    pivot.columns=pivot.columns.droplevel()
    df=pd.merge(polygon_layer,pivot,how="left",on="FID")
    gdf=gpd.GeoDataFrame(df,crs=polygon_layer.crs,geometry=df.geometry)
    return gdf

if __name__ == "__main__":

    #===== EXAMPLE FOR crop_to_bounds FUNCTION===========================
    # gdf2012=gpd.read_file("Coastal_Zones_CORINE_Land_Use/CZ_2012_DU004_3035_V010_fgdb/CZ_2012_DU004_3035_V010/Data/CZ_2012_galicia.shp")
    # gdf2018=gpd.read_file("Coastal_Zones_CORINE_Land_Use/CZ_2018_DU004_3035_V010_fgdb/CZ_2018_DU004_3035_V010/Data/CZ_2018_galicia.shp")
    # print("Originalmente el archivo del 2012 tiene la extensión",gdf2012.geometry.total_bounds)
    # print("Originalmente el archivo del 2018 tiene la extensión",gdf2018.geometry.total_bounds)
    # print("Cambiamos la extension")
    # newgdf2012=crop_to_bounds(gdf2012,gdf2018)
    # print("Ahora, el archivo del 2012 tiene la extensión",newgdf2012.geometry.total_bounds)



    #======= EJEMPLO PARA TRANS    #============================ POLYGONIZE ========================
    pol=polygonize("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_4.5_2026-2045.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_4.5_2026-2045.shp",save=True)
    print(pol)
    #vectorlulc=transform_lulc_vector(code_table="/home/david/OneDrive/geo_data/Coastal_Zones_CORINE_Land_Use/tabla_codigos.xlsx",columns_in_table=["LULC","NAME","Capital","Type of Capital"],vector=vectorlulc,column_in_vector="CODE_5_12",path="/home/david/OneDrive/geo_data//Coastal_Zones_CORINE_Land_Use/CZ_2012_DU004_3035_V010_fgdb/CZ_2012_DU004_3035_V010/Data/CZ_2012_galicia_natural_built_capital.shp")
    #print(vectorlulc)

    #==== EJEMPLO PARA TRANSFORMAR SRC=======

    #see_vector_projection("/home/david/OneDrive/geo_data/Coastal_Zones_CORINE_Land_Use/CZ_2012_DU004_3035_V010_fgdb/CZ_2012_DU004_3035_V010/Data/CZ_2012_galicia_natural_capital.shp")
    #see_vector_projection("/home/david/OneDrive/geo_data/Coastal_Zones_CORINE_Land_Use/CZ_2012_DU004_3035_V010_fgdb/CZ_2012_DU004_3035_V010/Data/CZ_2012_galicia_built_capital.shp")
    

    #=============== EJEMPLO PARA RASTERIZAR ================================
    #raster=rasterize_vector_with_template("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/galician_coastline_buffer2-10km.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","fid","EPSG:3035","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/coastline_buffer2-10km.tif")
    #visualize_raster("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/coastline_buffer2-10km.tif","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/coastline_buffer2-10km.png",save=True)
    #visualize_raster("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.png",save=True)


    #====================== VISIBILITY ANALYSIS =====================================================
    #visibility=visibility_of_sea("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/coastline_buffer2-10km.tif")
    # fig=plt.figure()
    # plt.imshow(visibility)
    # plt.show(True)
    # save_array_as_raster(visibility,"/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/visibility.tif")


    #======================0 CREATE SQUARE GRID======================================
    #create_square_grid(input_file="/home/usuario/OneDrive/recreation/qgis/dissolvedaoi.shp",size=5000,output_file="/home/usuario/OneDrive/recreation/qgis/prueba_grid5km.shp")

    # #================================= DISTANCE TO POINTS =============
    # grid=distance_to_points("/home/usuario/OneDrive/geo_data/transport/comunication_nodes.shp","/home/usuario/OneDrive/recreation/qgis/recreation_sdata.shp","Com-nodes",output_file="/home/usuario/OneDrive/recreation/qgis/recreation_sdata.shp")
    # print(grid)

    #============================ POLYGONIZE ========================
    pol=polygonize("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_4.5_2026-2045.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_4.5_2026-2045.shp",save=True)
    print(pol)