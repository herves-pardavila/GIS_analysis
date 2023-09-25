import pandas as pd
import numpy as np

import requests
import inspect
import operator
import itertools
from functools import reduce

"""
Created on Fri Oct 28 12:38:13 2022

Script para acceder a los datos de la red eleéctrica de España. Está inspirado en una clase creada en:
https://github.com/AyrtonB/RED-Electricity-API-Wrapper.git pero nos valió porque debe haber un error en la clase 
que hace que sólo funcione para dos combinaciones de categoría-widget.
@author: david
"""

def request(categoria,wiget,inicio,fin,lugar="peninsular",ts="month"):

    API_stream_url = f'https://apidatos.ree.es/en/datos/{categoria}/{widget}'
    params = {
        'start_date' : inicio,
        'end_date' : fin,
        'time_trunc' : ts,
        "geo_limit": lugar
    }
    r = requests.get(API_stream_url, params=params)
    r_json = r.json()
    precio_total=r_json["included"][-1]
    df = pd.json_normalize(precio_total["attributes"]["values"])

    return df


if __name__ == "__main__":

    RED_documentation_url = 'https://www.ree.es/en/apidatos'
    r = requests.get(RED_documentation_url)
    documentation_tables = pd.read_html(r.content)

    df_category_widgets = (documentation_tables
                        [4]
                        .drop(columns='lang')
                        .dropna(how='all')
                        .ffill()
                        .reset_index(drop=True)
                        )

    category_widget_combos = df_category_widgets.apply(lambda s: tuple(s.tolist()), axis=1).tolist()

    #print(category_widget_combos)

    category = 'mercados'
    widget="componentes-precio"
    
    años=np.arange(2015,2021,2) #la API sólo permite obtener un máximo de 24 meses en cada consulta si el paso temporal es el mes o de 366 días si el paso temporal son los dias


    for año in años:
        if año==años[0]:
            dataframe=request(category,widget,str(año)+"-01-01T00:00",str(año+1)+"-12-31T22:00")
        else:
           dataframe=dataframe.append(request(category,widget,str(año)+"-01-01T00:00",str(año+1)+"-12-31T22:00"))
    print(dataframe)
    dataframe.to_csv("REprecios.csv",index=False)


 




