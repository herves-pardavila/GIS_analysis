"""Importing request from OpenStreetMaps. 
Read: https://towardsdatascience.com/loading-data-from-openstreetmap-with-python-and-the-overpass-api-513882a27fd0"""


import requests
import json
import pandas as pd
import geopandas as gpd
import time
overpass_url = "http://overpass-api.de/api/interpreter"

def overpass_node_query(south,west,north,east,key,value):
    
    overpass_query= """
    [out:json];
    node[%s=%s]
    (%f,%f,%f,%f);
    (._;>;);
    out center;
    
    """%(key,value,south,west,north,east)

    response=requests.get(overpass_url,params={"data":overpass_query})
    res=response.json()
    version=res["version"]
    generator=res["generator"]
    time_of_query=res["osm3s"]["timestamp_osm_base"]
    lats=[elem["lat"] for elem in res["elements"]]
    lons=[elem["lon"] for elem in res["elements"]]
    ids=[elem["id"] for elem in res["elements"]]
    df=pd.DataFrame({"id":ids,"lat":lats,"lon":lons})
    df["Version"]=version
    df["Generator"]=generator
    df["time"]=time_of_query
    df["type"]="node"
    df["key"]=key
    df["Value"]=value
    gdf=gpd.GeoDataFrame(df,crs="EPSG:4326",geometry=gpd.points_from_xy(df.lon,df.lat))
    print(gdf)
    return gdf

def overpass_way_query(south,west,north,east,key,value):

    overpass_query= """
    [out:json];
    way[%s=%s]
    (%f,%f,%f,%f);
    (._;>;);
    out center;
    
    """%(key,value,south,west,north,east)
    response=requests.get(overpass_url,params={"data":overpass_query})
    res=response.json()
    version=res["version"]
    generator=res["generator"]
    time_of_query=res["osm3s"]["timestamp_osm_base"]
    ids,lats,lons=[],[],[]
    for elem in res["elements"]:
        if elem["type"]=="way":
            ids+=[elem["id"]]
            lats+=[elem["center"]["lat"]]
            lons+=[elem["center"]["lon"]]
    df=pd.DataFrame({"id":ids,"lat":lats,"lon":lons})
    df["Version"]=version
    df["Generator"]=generator
    df["time"]=time_of_query
    df["type"]="node"
    df["key"]=key
    df["Value"]=value
    gdf=gpd.GeoDataFrame(df,crs="EPSG:4326",geometry=gpd.points_from_xy(df.lon,df.lat))
    print(gdf)
    return gdf




    return res

if __name__ =="__main__":


    x=1
    # overpass_query = """
    # [out:json];
    # area["ISO3166-1"="DE"][admin_level=2];
    # (node["amenity"="biergarten"](area);
    #  way["amenity"="biergarten"](area);
    #  rel["amenity"="biergarten"](area);
    # );
    # out center;
    # """
    # response = requests.get(overpass_url, 
    #                         params={'data': overpass_query})
    # data = response.json()

    # restaurant=overpass_node_query(41.5,-10.0,44,-7,"amenity","restaurant")
    # restaurant.to_crs("EPSG:3035",inplace=True)
    # restaurant.to_file("/home/usuario/OneDrive/geo_data/OpenSteetMaps/restaurants_galicia.shp",index=False)
    
    # cities=overpass_node_query(42,-10.0,44,-7,"place","city")
    # cities.to_crs("EPSG:3035",inplace=True)
    # print(cities)
    # cities.to_file("/home/usuario/OneDrive/geo_data/OpenStreetMaps/cities_galicia.shp",index=False)

    # towns=overpass_node_query(41.5,-10.0,44,-7,"place","town")
    # towns.to_crs("EPSG:3035",inplace=True)
    # towns.to_file("/home/usuario/OneDrive/geo_data/OpenSteetMaps/towns_galicia.shp",index=False)
    
    
    # lighthouse=overpass_node_query(41.5,-10.0,44,-7,"man_made","lighthouse")
    # lighthouse.to_crs("EPSG:3035",inplace=True)
    # lighthouse.to_file("/home/usuario/OneDrive/geo_data/OpenSteetMaps/lighthouse_galicia.shp",index=False)

    # viewpoint=overpass_node_query(41.5,-10.0,44,-7,"tourism","viewpoint")
    # viewpoint.to_crs("EPSG:3035",inplace=True)
    # viewpoint.to_file("/home/usuario/OneDrive/geo_data/OpenSteetMaps/viewpoints_galicia.shp",index=False)
    

    # tourism=pd.concat([hotels,motels,hostels,camp_pitch,camp_site,caravan_site])
    # tourism=gpd.GeoDataFrame(tourism,crs="EPSG:4326",geometry=tourism.geometry)
    # tourism.to_crs("EPSG:3035",inplace=True)
    # tourism.to_file("/home/usuario/OneDrive/recreation/qgis/OSMtourism.shp")

    tnodes=overpass_way_query(41.5,-10,44,-7,"public_transport","station")
    airports=overpass_way_query(41.5,-10,44,-7,"aeroway","terminal")
    tnodes=pd.concat([tnodes,airports])
    tnodes=gpd.GeoDataFrame(tnodes,crs="EPSG:4326",geometry=tnodes.geometry)
    tnodes.to_crs("EPSG:3035",inplace=True)
    tnodes.to_file("/home/usuario/OneDrive/geo_data/OpenStreetMaps/tnodes.shp",index=False)
    print(tnodes)