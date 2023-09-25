import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import netCDF4
import datetime
import xarray
from meteogalicia import Meteogalicia
from qgis_commands import voronoi_polygons
from geo_funciones import spatial_overlay
import matplotlib
import random
class ESGF():
    """Class for working with downloaded files of the ESGF server. From my experience is better to work with downloaded files
    because the different servers it supports can stop working for some days. It has been tested only with nc files.
    """
    
    def __init__(self,path):
        """Initialize the class

        Args:
            path (str): Path to an nc file. Other file formats have not been yet tested.
        """
        
        
        self.Data=netCDF4.Dataset(path)
        self.variables=self.Data.variables
        self.tunits=self.variables["time"].units
        
    
    def get_time_days_since_epoch(self,epoch):
        """Returns a date index for the file when time units are days since a given 
        epoch. This epoch can be known by checking the self.tunits

        Args:
            epoch (datetime.datetime): epoch in the format datetime.datetime(yyyy,m,d,m,s) 

        Returns:
            dates (datetime64[ns]): date index in pandas format
        """
        fun1= lambda x  : epoch+datetime.timedelta(x)
        dates=list(map(fun1,np.array(self.variables["time"])))
        dates=pd.to_datetime(dates)
        return dates

    def get_time_bounds(self,dates,start_year,end_year):
        """Returns the indices of the time dimensions that correspond to a given start and end years.

        Args:
            dates (self.get_time_days_since_epoch): Time vector extracted from the file.    
            start_year (int): 
            end_year (int): _description_

        Returns:
            arg_start (int): Index in the dimensions that corresponds to the first elemt of the starting year.
            arg_finish (int): Index in the dimensions that corresponds to the last elemt of the last year.
        """
        arg_start=np.argwhere(dates.year==start_year)[0]
        arg_finish=np.argwhere(dates.year==end_year)[-1]  
        return arg_start,arg_finish
    
    def bbox(self,grid_type,latname,lonname,latmin,lonmin,latmax,lonmax,variable_name,epoch):
        """Function that selects the data of a given bounding box of coordinates. Before directly applying this 
        function you should check your type of grid, short name of the variable in the nc file, which short names 
        are set for latitude and longitude and if the fucntion is getting the coordinates well.

        Args:
            grid_type (str): Type of the grid. For example primarily when lat and lon are monotously increasing.
            latname (str): Short name for latitude in the nc file.
            lonname (_type_): Short name for longitude in the nc file.
            latmin (float): Bounding box
            lonmin (float):  Bounding box
            latmax (float):  Bounding box
            lonmax (_type_):  Bounding box
            variable_name (str):  Short name of the variable.
            epoch (datetime.datetime): Current class version only supports nc data where time is understood as days since an epoch.
                                        This epoch must be provided in the format datetime.datetime(year,month,day,hour,minute)

        Returns:
            pandas.DataFrame: Data bounded to the box limits
        """

        lats=self.Data.variables[latname]
        lons=self.Data.variables[lonname]
        #get indices for the box bounds
        if grid_type == "primarily":           
            dtheta=2*lats[-1]/len(lats)
            dphi=lons[-1]/len(lons)
            jmin=int(0.5*len(lats)+int(latmin/dtheta))
            jmax=int(0.5*len(lats)+int(latmax/dtheta))
            imin=int((lonmin-lons[0])/dphi)
            imax=int((lonmax-lons[0])/dphi)
        #clip data to this bounding box
        data=self.Data.variables[variable_name][:,jmin:jmax+1,imin:imax+1]
        #get a dataframe
        lats=lats[jmin:jmax+1]
        lons=lons[imin:imax+1]
        data=xarray.DataArray(data)
        #get dates to add them afterwards
        t=self.get_time_days_since_epoch(epoch)
        arrays=[t,lats,lons]
        multi_index=pd.MultiIndex.from_product(arrays,names=("time","lat","lon"))
        df=data.to_dataframe(name=variable_name)
        newdf=df.reindex(multi_index,copy=False)
        newdf=newdf.reset_index(level=["time","lat","lon"])
        newdf[variable_name]=list(df[variable_name])
        #Always check if this last line is  neccesary
        newdf["lon"]=newdf.lon-180.0
        return newdf
    


if __name__ =="__main__":

    #comparison of historical of the model with meteogalicia data. Example for the pressure
        
    #meteogalicia observational data
    met=Meteogalicia()
    datos=met.get_time_series(variable_code="PR_AVG_1.5m",starting_year=2000,ending_year=2014,freq="monthly")    
    datos=datos[datos.altitude < 100] #we are not interested in high altitude stations

    #model hist-data gfdl
    con_GFDL=ESGF("/home/usuario/OneDrive/geo_data/CMIP6/CMIP/ps_Emon_GFDL-CM4_historical_r1i1p1f1_gr1_195001-201412.nc")
    df_gfdl=con_GFDL.bbox("primarily","lat","lon",41.5,170,44,174,"ps",datetime.datetime(1850,1,1,0,0))

    #model mpi data
    files=["/home/usuario/OneDrive/geo_data/CMIP6/CMIP/ps_Amon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200001-200412.nc","/home/usuario/OneDrive/geo_data/CMIP6/CMIP/ps_Amon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_200501-200912.nc","/home/usuario/OneDrive/geo_data/CMIP6/CMIP/ps_Amon_MPI-ESM1-2-HR_historical_r1i1p1f1_gn_201001-201412.nc"]
    df_mpi=pd.DataFrame()
    for file in files:
        con_MPI=ESGF(file)
        df=con_MPI.bbox("primarily","lat","lon",41.5,170,44,174,"ps",datetime.datetime(1850,1,1,0,0))
        df_mpi=pd.concat([df_mpi,df])
    
    #model HiRAM-SIT-HR
    years=np.arange(2000,2014+1,1)
    df_hiram=pd.DataFrame()
    for year in years:
        con_hiram=ESGF("/home/usuario/OneDrive/geo_data/CMIP6/HighResMIP/ps_Amon_HiRAM-SIT-HR_hist-1950_r1i1p1f1_gn_%s01-%s12.nc" %(year,year))
        df=con_hiram.bbox("primarily","lat","lon",41.5,170,44,174,"ps",datetime.datetime(1948,1,1,0,0))
        df_hiram=pd.concat([df_hiram,df])
    
    
    #create shapefiles
    meteogalicia_data=datos
    df_hiram.time=df_hiram.time.astype(str)
    df_mpi.time=df_mpi.time.astype(str)
    df_gfdl.time=df_gfdl.time.astype(str)
    
    gdf_hiram=gpd.GeoDataFrame(df_hiram,crs="EPSG:4326",geometry=gpd.points_from_xy(df_hiram.lon,df_hiram.lat))
    gdf_mpi=gpd.GeoDataFrame(df_mpi,crs="EPSG:4326",geometry=gpd.points_from_xy(df_mpi.lon,df_mpi.lat))
    gdf_gfdl=gpd.GeoDataFrame(df_gfdl,crs="EPSG:4326",geometry=gpd.points_from_xy(df_gfdl.lon,df_gfdl.lat))

    gdf_meteogalicia=gpd.GeoDataFrame(meteogalicia_data,crs="EPSG:4326",geometry=gpd.points_from_xy(meteogalicia_data.lon,meteogalicia_data.lat))
    
    #create voronoi polygons for computing later the spatial overlay
    gdf_hiram.to_file("/home/usuario/OneDrive/geo_data/CMIP6/presion_hiram.shp",index=False)
    gdf_mpi.to_file("/home/usuario/OneDrive/geo_data/CMIP6/presion_mpi.shp",index=False)
    gdf_gfdl.to_file("/home/usuario/OneDrive/geo_data/CMIP6/presion_gfdl.shp",index=False)
    gdf_meteogalicia.to_file("/home/usuario/OneDrive/geo_data/CMIP6/presion_meteogalicia.shp",index=False)
    voronoi_polygons(10,"/home/usuario/OneDrive/geo_data/CMIP6/presion_hiram.shp","/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_hiram.shp")
    voronoi_polygons(10,"/home/usuario/OneDrive/geo_data/CMIP6/presion_mpi.shp","/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_mpi.shp")
    voronoi_polygons(10,"/home/usuario/OneDrive/geo_data/CMIP6/presion_gfdl.shp","/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_gfdl.shp")
    voronoi_polygons(10,"/home/usuario/OneDrive/geo_data/CMIP6/presion_meteogalicia.shp","/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_meteogalicia.shp")


    meteogalicia_data.rename(columns={"Dates":"time"},inplace=True)
    
    hiram=spatial_overlay(meteogalicia_data,df_hiram,"/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_meteogalicia.shp","/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_hiram.shp","EPSG:3035","idStation")
    hiram.time=pd.to_datetime(hiram.time)
    hiram["Month"]=hiram.time.dt.month

    mpi=spatial_overlay(meteogalicia_data,df_mpi,"/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_meteogalicia.shp","/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_mpi.shp","EPSG:3035","idStation")
    mpi.time=pd.to_datetime(mpi.time)
    mpi["Month"]=mpi.time.dt.month

    gfdl=spatial_overlay(meteogalicia_data,df_gfdl,"/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_meteogalicia.shp","/home/usuario/OneDrive/geo_data/CMIP6/voronoi_presion_gfdl.shp","EPSG:3035","idStation")
    gfdl.time=pd.to_datetime(gfdl.time)
    gfdl["Month"]=gfdl.time.dt.month

    # grouping by

    hiram=hiram[hiram.Value != -9999.0]
    hiram=hiram.groupby(["Month","idStation"],as_index=False).mean()
    #hiram=hiram[hiram.Value < 1040]
    hiram.ps=hiram.ps/100
    print(hiram.corr())
    hiram_medias=hiram.groupby("Month",as_index=False).mean()
    hiram_std=hiram.groupby("Month",as_index=False).std()
    hiram_meses=pd.merge(hiram_medias,hiram_std[["Month","ps","Value"]],how="inner",on="Month",suffixes=["_mean","_std"])


    mpi=mpi[mpi.Value != -9999.0]
    mpi=mpi.groupby(["Month","idStation"],as_index=False).mean()
    #mpi=mpi[mpi.Value < 1040]
    mpi.ps=mpi.ps/100
    print(mpi.corr())
    mpi_medias=mpi.groupby("Month",as_index=False).mean()
    mpi_std=mpi.groupby("Month",as_index=False).std()
    mpi_meses=pd.merge(mpi_medias,mpi_std[["Month","ps","Value"]],how="inner",on="Month",suffixes=["_mean","_std"])
  

    gfdl=gfdl[gfdl.Value != -9999.0]
    gfdl=gfdl.groupby(["Month","idStation"],as_index=False).mean()
    gfdl=gfdl[gfdl.Value < 1040]
    gfdl.ps=gfdl.ps/100
    print(gfdl.corr())
    gfdl_medias=gfdl.groupby("Month",as_index=False).mean()
    gfdl_std=gfdl.groupby("Month",as_index=False).std()
    gfdl_meses=pd.merge(gfdl_medias,gfdl_std[["Month","ps","Value"]],how="inner",on="Month",suffixes=["_mean","_std"])
    

    #plotting

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_xlabel("Mes")
    ax.set_ylabel("Presión (hPa)")
    ax.errorbar(hiram_meses.Month,hiram_meses.ps_mean,yerr=hiram_meses.ps_std,fmt="-o",label="HiRAM-SIT-HR")
    ax.errorbar(mpi_meses.Month,mpi_meses.ps_mean,yerr=mpi_meses.ps_std,fmt="-o",label="MPI-ESM1.2-HR")
    ax.errorbar(gfdl_meses.Month,gfdl_meses.ps_mean,yerr=gfdl_meses.ps_std,fmt="-o",label="GFDL-CM4")
    ax.errorbar(gfdl_meses.Month,gfdl_meses.Value_mean,yerr=gfdl_meses.Value_std,fmt="-o",label="Meteogalicia")
    fig.legend(loc="upper center",ncol=4,mode="expand")
    plt.show()
    fig.savefig("/home/usuario/OneDrive/recreation/InVEST/presion.pdf")



   



















            





        