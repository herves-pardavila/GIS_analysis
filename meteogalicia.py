import requests
import json
import pandas as pd
import numpy as np
import geopandas as gpd
import utm
import netCDF4
import time
import matplotlib.pyplot as plt
import contextily as ctx

class Meteogalicia:
    """Class for accesing Metogalicia's data. The class covers the THREDDS server for only one of the models that can be found there. However the function
    for thredds data can be easily modified to acces  """
    def __init__(self):
        self.basic_url_mensuais="https://servizos.meteogalicia.gal/mgrss/observacion/datosMensuaisEstacionsMeteo.action"
    
    def stations_info(self):
        """Returns information about all the meteorological stations

        Returns:
            pandas.DataFrame: Altitude, name, county and coordinates of the stations
        """
        r=requests.get("https://servizos.meteogalicia.gal/mgrss/observacion/listaEstacionsMeteo.action")
        myjson=r.json()
        return pd.DataFrame.from_dict(myjson["listaEstacionsMeteo"])
    
    def get_stations_data(self,variable_code,start_date,finish_date,frequency="daily"):
        """Gets the data from MeteoGalicia

            To check the available variables, go to https://www.meteogalicia.gal/observacion/rede/parametrosIndex.action and choose
            the desired requency: 10-minute, hourly, daily, monthly. The function is only built for the two most common cases: daily and monthly
        Args:
            variable_code (str): Check them in the link above
            start_date (str): Starting date as dd/mm/yyyy
            finish_date (str): Finish date as dd/mm/yyyy
            frequency (str, optional): _description_. Defaults to "daily".

        Returns:
            _type_: _description_
        """
        if frequency== "daily":
            freq="Diarios"
        elif frequency == "monthly":
            freq="Mensuais"
            
        r=requests.get("https://servizos.meteogalicia.gal/mgrss/observacion/datos"+freq+"EstacionsMeteo.action?idParam=%s&dataIni=%s&dataFin=%s" %(variable_code,start_date,finish_date))
        myjson=r.json()
        #print(myjson)
        dates=pd.json_normalize(myjson,record_path=["listDatos"+freq]) #dates
        #print(dates)
        for i in range(len(dates)):
            stations=[dates.listaEstacions.iloc[i][j]["estacion"] for j in range(len(dates.listaEstacions.iloc[i]))]
            concellos=[dates.listaEstacions.iloc[i][j]["concello"] for j in range(len(dates.listaEstacions.iloc[i]))]
            idestaciones=[dates.listaEstacions.iloc[i][j]["idEstacion"] for j in range(len(dates.listaEstacions.iloc[i]))]
            provincia=[dates.listaEstacions.iloc[i][j]["provincia"] for j in range(len(dates.listaEstacions.iloc[i]))]
            utmx=[dates.listaEstacions.iloc[i][j]["utmx"] for j in range(len(dates.listaEstacions.iloc[i]))]
            utmy=[dates.listaEstacions.iloc[i][j]["utmy"] for j in range(len(dates.listaEstacions.iloc[i]))]
            codigoparam=[dates.listaEstacions.iloc[i][j]["listaMedidas"][0]["codigoParametro"] for j in range(len(dates.listaEstacions.iloc[i]))]
            lncodval=[dates.listaEstacions.iloc[i][j]["listaMedidas"][0]["lnCodigoValidacion"] for j in range(len(dates.listaEstacions.iloc[i]))]
            nomeparametro=[dates.listaEstacions.iloc[i][j]["listaMedidas"][0]["nomeParametro"] for j in range(len(dates.listaEstacions.iloc[i]))]
            units=[dates.listaEstacions.iloc[i][j]["listaMedidas"][0]["unidade"] for j in range(len(dates.listaEstacions.iloc[i]))]
            value=[dates.listaEstacions.iloc[i][j]["listaMedidas"][0]["valor"] for j in range(len(dates.listaEstacions.iloc[i]))]
            date=[dates.data.iloc[i]]*len(stations)

            date=[dates.data.iloc[i]]*len(stations)
            if i== 0 :
                df=pd.DataFrame({"Dates":date,"Stations":stations,"Counties":concellos,"idStation":idestaciones,"Province":provincia,"utmx":utmx,"utmy":utmy,"CodeParameter":codigoparam,"CodeValidation":lncodval,"NameParameter":nomeparametro,"Units":units,"Value":value})
            else:
                df=pd.concat([df,pd.DataFrame({"Dates":date,"Stations":stations,"Counties":concellos,"idStation":idestaciones,"Province":provincia,"utmx":utmx,"utmy":utmy,"CodeParameter":codigoparam,"CodeValidation":lncodval,"NameParameter":nomeparametro,"Units":units,"Value":value})])
        #adding more information about stations: altitude and coordinates in degrees
        info_stations=self.stations_info()
        df=df.merge(info_stations[["altitude","idEstacion","lat","lon"]],how="inner",left_on=["idStation"],right_on=["idEstacion"])
        
        return df
    
    def get_time_series(self,variable_code,starting_year,ending_year,freq):

        df=pd.DataFrame()
        for year in range(starting_year,ending_year+1):
            print(year)
            for month in range(1,13):
                print(month)
                newdf=self.get_stations_data(variable_code=variable_code,start_date="01/"+"%02d" %month+"/"+str(year),finish_date="31/"+"%02d" %month+"/"+str(year),frequency=freq)
                df=pd.concat([df,newdf])
                time.sleep(2)
            time.sleep(1)
        return df
    def find_indices(self,Lats,Lons,latitude,longitude): 
        """
        Función usada en archivos.nc donde la latitud y longitud vienen dadas como
        una malla[i,j]. Devuelve los i,j correspondientes a la coordenada deseada.

        Parameters
        ----------
        Lats : clase (depende de con qué librería lo abras) 
            Clase con las latitudes del archivo.nc
        Lons : clase (depende de con qué librería lo abras)
            Clase con las longitudes del archivo.nc
        latitude : float
            Latitud cuyos indices queremos saber
        longitude : float
            Longitud cuyos índices queremos saber

        Returns
        -------
        x: float
            Indice i tal que i,j corresponde a las coordenadas de interés
        y: float
            Indice j tal que i,j corresponde a las coordenadas de interés

        """
        abslat = np.abs(np.asarray(Lats)-latitude)
        abslon = np.abs(np.asarray(Lons)-longitude)
        c = np.maximum(abslon, abslat)
        #print(np.where(c == np.min(c)))
        a = np.where(c == np.min(c))
        x,y=a[0][0],a[1][0]

        return x,y

    def thredds_data(self,variable_name,model="WRF_HIST/d02",start_year=2010,end_year=2021):
        
        latmin,lonmin,latmax,lonmax=41.0,-10.0,44.0,-7.0
        f=netCDF4.Dataset("http://mandeo.meteogalicia.gal/thredds/dodsC/modelos/"+str(model)+"/"+str(2020)+"/"+str(10)+"/"+"wrf_arw_det_history_d02_"+str(2020)+str(10)+str(10)+"_1200.nc4")
        lat=f.variables["lat"]
        lon=f.variables["lon"]

        imin,jmin=self.find_indices(lat,lon,latmin,lonmin)
        imax,jmax=self.find_indices(lat,lon,latmax,lonmax)
        
        lat=lat[imin:imax,jmin:jmax]
        lon=lon[imin:imax,jmin:jmax]

        years=np.arange(start_year,end_year,1)
        months=np.arange(1,13,1)
        days=[31,28,31,30,31,30,31,31,30,31,30,31]
        thredds_df=pd.DataFrame()
        for year in years:
            print(year)
            for month in months:
                print(month)
                total_days=days[month-1]
                for day in range(1,total_days+1):
                    try:
                        f=netCDF4.Dataset("http://mandeo.meteogalicia.gal/thredds/dodsC/modelos/"+model+"/"+str(year)+"/"+"%02d"%month+"/"+"wrf_arw_det_history_d02_"+str(year)+"%02d"%month+"%02d"%day+"_1200.nc4")
                        var=f.variables[variable_name][:,imin:imax,jmin:jmax]
                        var1=np.mean(var,axis=0)
                    except OSError:
                        continue
                    try:
                        f=netCDF4.Dataset("http://mandeo.meteogalicia.gal/thredds/dodsC/modelos/"+model+"/"+str(year)+"/"+"%02d"%month+"/"+"wrf_arw_det_history_d02_"+str(year)+"%02d"%month+"%02d"%day+"_0000.nc4")
                        var=f.variables[variable_name][:,imin:imax,jmin:jmax]
                        var2=np.mean(var,axis=0)
                    except OSError:
                        continue
                    var=0.5*(var1+var2) #daily average of the variable
                    temporal_df=pd.DataFrame({"Date":" ", "Latitude":np.ravel(lat),"Longitude":np.ravel(lon),"value":np.ravel(var)})
                    temporal_df.Date=pd.Timestamp(year,month,day,12)
                    temporal_df["Name"]=f.variables[variable_name].long_name
                    temporal_df["Units"]=f.variables[variable_name].units
                    temporal_df["Grid"]=f.variables[variable_name].grid_mapping
                    thredds_df=pd.concat([thredds_df,temporal_df])



        return thredds_df
    
if __name__ == "__main__":

    #initiate the class
    con=Meteogalicia()

    #see info about stations
    print(con.stations_info())
    
    #map stations with their altitude
    df=pd.DataFrame(con.stations_info())
    gdf=gpd.GeoDataFrame(data=df,crs="EPSG:4326",geometry=gpd.points_from_xy(df.lon,df.lat))
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.title.set_text("Altitude (m) and location of the stations")
    gdf.plot("altitude",ax=ax,legend=True)
    ctx.add_basemap(ax=ax, crs=gdf.crs, source=ctx.providers.OpenStreetMap.DE.url)
    plt.show()

    # #get maximum daily temperature at 15 m height for 2018 in the geojson services
    # data=con.get_stations_data("TA_MAX_15m","01/01/2018","31/12/2018","daily")
    # print(data)

    #get monthly average wind speed  at 2m in the geojson services during years 2020-2022 and plot mean at each station
    data=con.get_time_series("HSOL_SUM_1.5m",2020,2022,freq="monthly")
    print(data)
    data=data.groupby(by="idStation",as_index=False).mean(numeric_only=True)
    new_gdf=gpd.GeoDataFrame(data=data,crs="EPSG:4326",geometry=gpd.points_from_xy(data.lon,data.lat))
    print(new_gdf)
    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.title.set_text("Mean wind speed (m/s)")
    ctx.add_basemap(ax=ax2, crs=new_gdf.crs, source=ctx.providers.OpenStreetMap.DE.url)
    new_gdf.plot("Value",ax=ax2,legend=True)
    plt.show()
