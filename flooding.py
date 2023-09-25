import sys
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from scipy import ndimage 
from tqdm import tqdm
sys.path.append("/home/usuario/OneDrive/geo_data")
from geo_funciones import rasterize_vector_with_template,visualize_raster,save_array_as_raster, extend_polygons_to_line

def plot_array(array):
    fig=plt.figure()
    plt.imshow(array)
    plt.colorbar()
    plt.show(True)
    return

def find_clusters(array):
    clustered = np.empty_like(array)
    unique_vals = np.unique(array)
    cluster_count = 0
    for val in unique_vals:
        labelling, label_count = ndimage.label(array == val)
        for k in range(1, label_count + 1):
            clustered[labelling == k] = cluster_count
            cluster_count += 1
    return clustered, cluster_count

def flooding(dem_raster,rise_raster,output_raster):
    # DIGITAL ELEVATION MODEL
    dem=rasterio.open(dem_raster)
    dem_array=dem.read(1)
    dem_array[dem_array==dem_array[0,0]]=1e9

    # SEA LEVEL RISE
    
    beaches=rasterio.open(rise_raster)
    beaches_array=beaches.read(1)

    flooded=np.zeros(np.shape(dem_array))

    slr_values=np.unique(beaches_array)
    slr_values=slr_values[slr_values != 0]

    for rise in tqdm(slr_values):
        print(rise)
        possibly_flooded=dem_array <= rise
        possibly_flooded[possibly_flooded==True]=1.0
        possibly_flooded[possibly_flooded==False]=0.0
        
        possibly_flooded=possibly_flooded+beaches_array
        possibly_flooded[possibly_flooded != 0.0]=1.0
    
        clusters,cluster_count=find_clusters(possibly_flooded)


        for cluster in range(1,cluster_count+1):
            if np.any(beaches_array[clusters==cluster]==rise):
                flooded[clusters==cluster]=1
    #plot_array(flooded)
    save_array_as_raster(flooded,dem_raster,output_raster)
    return flooded


if __name__ == "__main__":

    
    #x=extend_polygons_to_line("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2026-2045.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/galician_coastline.shp",1,"Aumento (m","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2026-2045_all_coastline.shp",save=True,plot=True)
    #x=extend_polygons_to_line("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2081-2100.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/galician_coastline.shp",1,"Aumento (m","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2081-2100_all_coastline.shp",save=True,plot=True)
    #x=extend_polygons_to_line("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2026-2045.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/galician_coastline.shp",1,"Aumento (m","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2026-2045_all_coastline.shp",save=True,plot=True)
    #x=extend_polygons_to_line("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2081-2100.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/galician_coastline.shp",1,"Aumento (m","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2081-2100_all_coastline.shp",save=True,plot=True)

    #rasterize_vector_with_template("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2026-2045_all_coastline.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","Aumento (m","EPSG:3035",np.float64,"/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2026-2045.tif")
    #rasterize_vector_with_template("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2081-2100_all_coastline.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","Aumento (m","EPSG:3035",np.float64,"/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2081-2100.tif")
    #rasterize_vector_with_template("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2026-2045_all_coastline.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","Aumento (m","EPSG:3035",np.float64,"/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2026-2045.tif")
    #rasterize_vector_with_template("/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2081-2100_all_coastline.shp","/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","Aumento (m","EPSG:3035",np.float64,"/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2081-2100.tif")

    flood1=flooding("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2026-2045.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_4.5_2026-2045.tif")
    # flood2=flooding("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_4.5_2081-2100.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_4.5_2081-2100.tif")
    # flood3=flooding("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2026-2045.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_8.5_2026-2045.tif")
    # flood4=flooding("/home/usuario/OneDrive/geo_data/DEM/eu_dem_v11_E20N20/Galicia_DEM.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/vector_data/GEAMA_BRUNN_aumento_medio_RCP_8.5_2081-2100.tif","/home/usuario/OneDrive/geo_data/2021_03_26_GEAMA_PIMA_Universidades/2021_03_26_GEAMA_PIMA_Universidades/GEAMA_BRUNN_floods_RCP_8.5_2081-2100.tif")

    plot_array(flood1)
    # plot_array(flood2)
    # plot_array(flood3)
    # plot_array(flood4)