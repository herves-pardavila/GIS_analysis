import http.client
import json
import requests
import pandas as pd
import geopandas as gpd
import time


class Aemet():
    """
    Clase para acceder a los datos de la Aemet. En concreto, a los valores 
    climátológicos históricos diarios de todas las estaciones. Enlace :
    https://opendata.aemet.es/dist/index.html?#!/valores-climatologicos/Climatolog%C3%ADas_diarias
    
    Para acceder a otro tipo de datos, habría que modificar en las funciones de
    abajo las url's e ir haciendo pruebas corrigiendo los errores
        
    """
    def __init__(self, api_key):
        """
        Establece los parámetros para realizar la conexión, como es la url de la
        Aemet y la api. Para solicitar API: https://opendata.aemet.es/centrodedescargas/altaUsuario?

        Parameters
        ----------
        api_key : str
            API que debemos solicitar previamente

        Returns
        -------
        None.

        """
        self.conn = http.client.HTTPSConnection("opendata.aemet.es")
        self.headers = {'cache-control': "no-cache", 'api_key':api_key}

    def get_estaciones(self):
        """
        Establece conexión con el cliente haciendo una request y obtiene una 
        lista de dicts con información de cada estación

        Returns
        -------
        response.json() : list
            Lista de jsons con información sobre coordenadas, provincia y 
            otros parámetros de las estaciones meteorológicas

        """
        self.conn.request("GET", "/opendata/api/valores/climatologicos/inventarioestaciones/todasestaciones", headers=self.headers)
        data = self.conn.getresponse().read()
        response = requests.get(json.loads(data.decode("utf-8"))['datos'])
        return response.json()

    def get_temperaturas_todas_estaciones(self, fechaini, fechafin):
        """
        Establece conexión y obtiene una lista donde cada elemento es un json 
        con el valor de las variables climatológicas de cada estación en un 
        instante de tiempo. La request solo permite solicitar 30 dias de 
        información en cada vez.

        Parameters
        ----------
        fechaini : str
            Inicio de la serie temporal en formato {año}-{mes}-{dia}T{hora}%3A{minuto}%3A{segundo}UTC
        fechafin : str
            Final de la serie temporal en formato {año}-{mes}-{dia}T{hora}%3A{minuto}%3A{segundo}UTC

        Returns
        -------
        list
            Lista donde cada elemto es un json con el nombre de la estación
            le fecha y hora de la obsrevación y los valores de las variables
            climatológicas.

        """
        self.conn.request("GET", "https://opendata.aemet.es/opendata/api/valores/climatologicos/diarios/datos/fechaini/{}/fechafin/{}/todasestaciones".format(fechaini, fechafin), headers=self.headers)
        data = self.conn.getresponse().read()
        print(data)
        response = requests.get(json.loads(data.decode("utf-8"))['datos'])
        return response.json()
    
    def get_time_series(get_temperaturas_todas_estaciones,añoini,añofin):
        """
        El cliente sólo permite obtener 31 días de información en cada conexión.
        Por lo tanto, esta fucnión itera sobre la función anterior para obtener
        de forma agrupada toda la serie temporal deseada. Tal y como está escrita
        la clase actualmente obtiene datos diarios (una observación al dia)

        Parameters
        ----------
        get_temperaturas_todas_estaciones : function
            Definida en la clase
            
        añoini : int
            Inicio de la serie temporal
        añofin : int
            Fin de la serie temporal

        Returns
        -------
        res : list
            Lista donde cada elemto es un json con el nombre de la estación
            le fecha y hora de la obsrevación y los valores de las variables
            climatológicas.

        """
        res=[]
        for año in range(añoini,añofin+1):
            print(año)
            for mes in range(1,13):
                print(mes)


                if mes==1:
                    fechaini = '{}-{}-01T00%3A00%3A00UTC'.format(año, mes)
                else:
                    fechaini = '{}-{}-02T00%3A00%3A00UTC'.format(año, mes)
                    
                if mes==12:
                    fechafin = '{}-{}-31T00%3A00%3A00UTC'.format(año, mes)

                else:
                    fechafin = '{}-{}-01T00%3A00%3A00UTC'.format(año, mes+1)
                res.append(con.get_temperaturas_todas_estaciones(fechaini, fechafin))
                time.sleep(2)
        return res
    
    def add_coordinates(dataframe,stations):
        """
        

        Añade a la serie temporal de observaciones las coordenadas de cada
        estación.
        ----------
        dataframe : pandas.DataFrame
            Observaciones de las estaciones 
        stations : list
            output de la función get_estaciones, parte de esta clase.

        Returns
        -------
        dataframe : pandas.DataFrame
            Dataframe con dos nuevas columnas, latitud y longitud de la 
            estación que realiza la observacion.

        """
        
        for station in dataframe["nombre"].unique():
            info = next( item for item in stations if item["nombre"]==station)
            dataframe.loc[dataframe["nombre"]==station,"latitud"]=info["latitud"]
            dataframe.loc[dataframe["nombre"]==station,"longitud"]=info["longitud"]
        return dataframe
                    
    def transform_latitude(lat):
        """_La Aemet devuelve las latitudes en un formato raro. De tal forma que 42.03567 Norte
        aparece como 4203567N. En esta función se modifica la latitud para convertirla a un formato
        más amigable. Actua sobre una latitud sólo, para cambiar todas habría que combinarla con un
        map.

        Args:
            lat (str): Latitud en el formato dado por la AEMET

        Returns:
            str: Latitud en formato decimal
        """
        return lat[0:2]+"."+lat[2:-1]

    def transform_longitude(lon):
        """La Aemet devuelve las longitudes en un formato raro. De tal forma que -0.89765 
        aparece como 0089765W. En esta función se modifica la longitud para convertirla a un formato
        más amigable. Actua sobre una longitud sólo, para cambiar todas habría que combinarla con un
        map.

        Args:
            lon (str): longitud en el formato dado por la AEMET

        Returns:
            str: longitud en el nuevo formato
        """

        if lon[-1]=="E":
            newlon=lon[0:2]+"."+lon[2:-1]
        else:
            newlon="-"+lon[0:2]+"."+lon[2:-1]
        return newlon

def check_provinces(object1,name1,object2,name2):
    """Printea las provincias de un objeto1 y un objeto2 para comprobar que tienen los mismos nombres.

    Args:
        object1 (_type_): ob
        name1 (str): nombre de la columna con las provincias en el objeto1.
        object2 (_type_): _description_
        name2 (str): nombre de la columna con las provincias en el objeto2.
    """
    for i in range(len(object1[name1].unique())):
        print(sorted(object1[name1].unique())[i],sorted(object2[name2].unique())[i])
    return



if __name__ == '__main__':  

    api_key = "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJkYXZpZC5oZXJ2ZXMucGFyZGF2aWxhQHVzYy5lcyIsImp0aSI6ImVhYmRmNWMwLWY1MGQtNGNkMC1iYTUwLTdkNzcwNzg4NjZhMCIsImlzcyI6IkFFTUVUIiwiaWF0IjoxNjYzNTc5NDE2LCJ1c2VySWQiOiJlYWJkZjVjMC1mNTBkLTRjZDAtYmE1MC03ZDc3MDc4ODY2YTAiLCJyb2xlIjoiIn0.yabR_VNglnoQx8L4hbZP--8cRzpLXZXF5z5Yocws6cI"
    con = Aemet(api_key)
    estaciones = con.get_estaciones() #todas las estaciones
    res=con.get_time_series(2015, 2019) #datos de todas las estaciones desde 2015 hasta 2019
    datos = pd.DataFrame(sum(res, [])) #lo guardamos como dataframe
    datos=datos[["fecha","nombre","provincia","altitud","tmed","tmin","tmax"]] #columnas que me interesan
    datos=Aemet.add_coordinates(datos,estaciones) #añadimos a cada estación sus coordendadas
    datos["latitud"]=list(map(Aemet.transform_latitude,datos["latitud"])) #cambiamos formato de latitudes
    datos["longitud"]=list(map(Aemet.transform_longitude,datos["longitud"])) #cambiamos formato de longitudes
    datos["latitud"]=datos["latitud"].astype(float) 
    datos["longitud"]=datos["longitud"].astype(float)
    #cambiamos comas por puntos como separador decimal en los datos 
    datos['tmed'] = datos['tmed'].str.replace(',', '.').astype(float) 
    datos['tmin'] = datos['tmin'].str.replace(',', '.').astype(float)
    datos['tmax'] = datos['tmax'].str.replace(',', '.').astype(float)
    print(datos)
    #creo  un geodataframe
    gdf=gpd.GeoDataFrame(datos,crs="EPSG:4258",geometry=gpd.points_from_xy(datos.longitud,datos.latitud))
    #como ya tengo la geometría, quito las columnas de latitud y longitud
    gdf=gdf.drop(["latitud","longitud"],axis=1)
    print(gdf)
    gdf.to_file("aemet_todas.shp")

    # #======PARTE PARA ELEGIR LA ESTACIÓN MÁS CERCANA A LA CAPITAL DE PROVINCIA====
    gdf_aemet=gdf
    del(gdf)
    gdf_capitales=gpd.read_file("capitales_de_provincia.shp")

    #Primero tengo que igualar los nombres de las provincias
    gdf_aemet.provincia=gdf_aemet.provincia.str.lower() #pongo minúsculas
    gdf_capitales.Provincia=gdf_capitales.Provincia.str.lower() #pongo minúsculas
    
    #tildes
    gdf_aemet.replace("araba/alava","araba",inplace=True)
    gdf_aemet.replace("avila","ávila",inplace=True)
    gdf_aemet.replace("castellon","castelló",inplace=True)
    gdf_aemet.replace("caceres","cáceres",inplace=True)
    gdf_aemet.replace("cadiz","cádiz",inplace=True)
    gdf_aemet.replace("cordoba","córdoba",inplace=True)
    gdf_aemet.replace("jaen","jaén",inplace=True)
    gdf_aemet.replace("leon","león",inplace=True)
    gdf_aemet.replace("malaga","málaga",inplace=True)
    gdf_aemet.replace("almeria","almería",inplace=True)
    gdf_aemet.replace("sta. cruz de tenerife","santa cruz de tenerife",inplace=True)
    #cambio de idioma
    gdf_capitales.replace("vizcaya","bizkaia",inplace=True)
    gdf_capitales.replace("álava","araba",inplace=True)

    #ya he conseguido que salgan en el mismo orden, hago un replace para que sean exactamente igual en ambos objetos
    for i in range(len(gdf_aemet.provincia.unique())):
        print(sorted(gdf_aemet.provincia.unique())[i],sorted(gdf_capitales.Provincia.unique())[i])
        gdf_capitales.replace(sorted(gdf_capitales.Provincia.unique())[i],sorted(gdf_aemet.provincia.unique())[i],inplace=True)
    #compruebo que efectivamente tienen los mismos nombres
    check_provinces(gdf_aemet,"provincia",gdf_capitales,"Provincia")
    estaciones_a_borrar=["SABADELL AEROPUERTO", "GUADALAJARA, EL SERRANILLO","LA ALDEA DE SAN NICOLÁS",
    "PUEBLA DE DON RODRIGO","ANTEQUERA","PORQUERES","ANDÚJAR","REUS AEROPUERTO","VILLAFRANCA DEL CID/VILLAFRANCA",
    "CASTELLFORT","VINARÒS","JACA","ARAGÜÉS DEL PUERTO","BIELSA","TORLA","IZAÑA","TENERIFE NORTE AEROPUERTO",
    "FORONDA-TXOKIZA"]
    mis_estaciones=[] #lista con los nombres de las estaciones más cercanas a cada capital de provincia
    gdf_aemet=gdf_aemet[~gdf_aemet.nombre.isin(estaciones_a_borrar)]
    for i in range(len(gdf_aemet.provincia.unique())):
        provincia=sorted(gdf_aemet.provincia.unique())[i] #provincia i
        print(provincia)
        estaciones=gdf_aemet[gdf_aemet.provincia==provincia].nombre.unique() #estaciones en esa provincia i
        print(estaciones)
        #coordenadas de todas las estaciones de una misma provincia
        s1=gpd.GeoSeries(gdf_aemet[gdf_aemet.nombre.isin(estaciones)].geometry.unique())
        #coordenadas de la capital de provincia
        point=gdf_capitales[gdf_capitales.Provincia==provincia].geometry.iloc[0] 
        #estaciones con la menor distancia a la capital
        mis_estaciones+=[estaciones[s1.distance(point).argmin()]]
        del(point)
    print(mis_estaciones)
    print(len(mis_estaciones))
    #selecciono en el geodataframe sólo las estaciones más próximas 
    gdf_aemet=gdf_aemet[gdf_aemet.nombre.isin(mis_estaciones)]
    print(gdf_aemet)
    gdf_aemet.to_file("aemet_capitales_provincia.shp",index=False)
    df=gdf_aemet
    df=df.drop("geometry",axis=1)
    df.to_csv("aemet_capitales_provincia.csv",index=False)

