import http.client
import json
import requests
import pandas as pd
import geopandas as gpd
import time
import matplotlib.pyplot as plt
import contextily as ctx
class Aemet():
    """
    Class for accesing State Meteorological Agency of Spain (AEMET) data. 
    This class was written for historic daily climate data. For other type of data you must
    change some of the functions below. A good starting point is to check the link  :
    https://opendata.aemet.es/dist/index.html?#!/valores-climatologicos/Climatolog%C3%ADas_diarias
    
        
    """
    def __init__(self, api_key):
        """
        Starts the conection. Obtain your API in the link: https://opendata.aemet.es/centrodedescargas/altaUsuario?

        Parameters
        ----------
        api_key : str
            API key
        Returns
        -------
        None.

        """
        self.conn = http.client.HTTPSConnection("opendata.aemet.es")
        self.headers = {'cache-control': "no-cache", 'api_key':api_key}

    def get_stations(self):
        """
        Makes a request and obtains information of metereological stations.

        Returns
        -------
        stations : pandas.DataFrame
            Coordinates, province, altitude and other metadata of metereological stations.

        """
        self.conn.request("GET", "/opendata/api/valores/climatologicos/inventarioestaciones/todasestaciones", headers=self.headers)
        data = self.conn.getresponse().read()
        response = requests.get(json.loads(data.decode("utf-8"))['datos'])
        response=pd.json_normalize(response.json())
        response=self.transform_coordinates(response)
        return response

    def get_data(self, initial_date, final_date):
        """
        List where each element is a json with climate data of an station in a given day. This
        type of request is for a maximum of 30 days.

        Parameters
        ----------
        initial_date : str
            Start of the time series in format {year}-{month}-{day}T{hour}%3A{minute}%3A{second}UTC
        final_date : str
            End of the time series in format {year}-{month}-{day}T{hour}%3A{minute}%3A{second}UTC

        Returns
        -------
        list
            List where each element is a json with station name, date, hour of measurement and climate data.

        """
        self.conn.request("GET", "https://opendata.aemet.es/opendata/api/valores/climatologicos/diarios/datos/fechaini/{}/fechafin/{}/todasestaciones".format(initial_date, final_date), headers=self.headers)
        data = self.conn.getresponse().read()
        response = requests.get(json.loads(data.decode("utf-8"))['datos'])
        return response.json()
    
    def get_time_series(self,first_year,last_year):
        """
        The client allows for requests of 31 days maximum. As a result, this functions calls and iterates the 
        previous function to obtain the desired time series. Daily data only.
  
        Parameters
        ---------
            
        first_year: int
            First year of the time series
        last_year : int
            End of the time series

        Returns
        -------
        res : list
            List where each element is a json with the name of the station, date, hour of observation
            and values of climate variables.

        """
        res=[]
        for year in range(first_year,last_year+1):
            #print(year)
            for month in range(1,13):
                #print(month)


                if month==1:
                    start_date = '{}-{}-01T00%3A00%3A00UTC'.format(year, month)
                else:
                    start_date = '{}-{}-02T00%3A00%3A00UTC'.format(year, month)
                    
                if month==12:
                    finish_date= '{}-{}-31T00%3A00%3A00UTC'.format(year, month)

                else:
                    finish_date = '{}-{}-01T00%3A00%3A00UTC'.format(year, month+1)
                res.append(self.get_data(start_date, finish_date))
                time.sleep(2)
        res=pd.DataFrame(sum(res,[]))
        #add coordinates and altitude of the stations
        stations=self.get_stations()
        res=pd.merge(res,stations[["latitud","longitud","nombre"]],on="nombre",how="inner")
        return res
    

    def transform_latitude(self,lat):
        """AEMET cordinate format is weird. 42.03567 North is given as "4203567N". This function
        removes the letter so coordinates are changed to a much more friendly format

        Args:
            lat (str): latitude in AEMET's format

        Returns:
            str: Latitude in degrees north
        """
        return lat[0:2]+"."+lat[2:-1]

    def transform_longitude(self,lon):
        """AEMET cordinate format is weird. -0.89765 longitude appears as 0089765W. his function
        removes the letter so coordinates are changed to a much more friendly format

        Args:
            lon (str): longitude in AEMET's format

        Returns:
           newlon (str): longitude in degrees east
        """

        if lon[-1]=="E":
            newlon=lon[0:2]+"."+lon[2:-1]
        else:
            newlon="-"+lon[0:2]+"."+lon[2:-1]
        return newlon
    def transform_coordinates(self,dataframe):
        """Applies the transformation of cordinates to pandas dataframe

        Args:
            dataframe (pd.DataFrame): AEMET data

        Returns:
            dataframe (pd.DataFrame): AEMET data with latitude and longitude in the new format
        """

        dataframe["latitud"]=list(map(self.transform_latitude,dataframe["latitud"])) 
        dataframe["longitud"]=list(map(self.transform_longitude,dataframe["longitud"])) 
        dataframe["latitud"]=dataframe["latitud"].astype(float) 
        dataframe["longitud"]=dataframe["longitud"].astype(float)
        return dataframe




if __name__ == '__main__':  

    #start the connection. Get your api key here: https://opendata.aemet.es/centrodedescargas/altaUsuario?
    api_key = "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJkYXZpZC5oZXJ2ZXMucGFyZGF2aWxhQHVzYy5lcyIsImp0aSI6ImVhYmRmNWMwLWY1MGQtNGNkMC1iYTUwLTdkNzcwNzg4NjZhMCIsImlzcyI6IkFFTUVUIiwiaWF0IjoxNjYzNTc5NDE2LCJ1c2VySWQiOiJlYWJkZjVjMC1mNTBkLTRjZDAtYmE1MC03ZDc3MDc4ODY2YTAiLCJyb2xlIjoiIn0.yabR_VNglnoQx8L4hbZP--8cRzpLXZXF5z5Yocws6cI"
    con = Aemet(api_key)
    
    #get and plot stations in a map
    stations = con.get_stations() 
    print(stations)
    gdf=gpd.GeoDataFrame(data=stations,crs="EPSG:4326",geometry=gpd.points_from_xy(stations.longitud,stations.latitud))
    fig=plt.figure()
    ax=fig.add_subplot(111)
    gdf.plot(ax=ax)
    ctx.add_basemap(ax=ax, crs=gdf.crs, source= ctx.providers.OpenStreetMap.DE.url)
    plt.show()


    res=con.get_time_series(2015, 2016) #datos de todas las estaciones desde 2015 hasta 2019
    print(res)
    print(res.info())


    

