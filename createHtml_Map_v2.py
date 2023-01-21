# -*- coding: utf-8 -*-

"""
Spyder Editor

This is a temporary script file.
"""
import sys
import os
import fiona
import pandas as pd
import math
import geopandas as gpd
import matplotlib.pyplot as plt
import folium
from folium.features import DivIcon
from shapely import geometry
from geopy.distance import distance
from shapely.geometry import Polygon, Point, LineString, LinearRing
import numpy as np
import csv
import descartes
# import geojson
import json
from geojson import Feature, Point, FeatureCollection, dumps
import tkinter as tk
from tkinter import filedialog as fd
import interpolation
from scipy.interpolate import griddata, interp2d
import seaborn as sns
import webbrowser as web

iLS_Dict = {} #global
list_pointname =[] #Liste
list_latitude = []
list_longitude = []
list_profondeur = []
grp_name = []
geo_df_list = {}
length = {}
ecart = {}
archi = {}
XY = {}
LS_inter = {}
LineS = {}
poly_table = {}

def searchFolder(file):
    return os.path.join(os.getcwd(), 'files' + file)

def read_df(df):
    data = {}
    for row in df.itertuples():
        values = row[1].split(";")
        data[values[0]]={"coordinates":{"latitude": float(values[1]),
                                        "longitude": float(values[2])},
                         "complements":{
                                        "profondeur": float(values[3])}}
    return data

def haversine_distance(pt1, pt2):
    lat1, lon1 = pt1
    lat2, lon2 = pt2
    
    # Convertir les degrés en radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    # Calculer les différences
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Appliquer la formule
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))

    # Rayon de la Terre (moyen) en kilomètres
    R = 6371

    return c * R

def closest_poly_table(data):
    poly_table = {}
    for key1, value1 in data.items():
        if isinstance(key1, str) and "PZ" in key1:
            for key2, value2 in value1.items():
                if key2 == "coordinates":
                    pt1 = (value2['longitude'], value2['latitude'])
                    poly_table[key1] = {'first_point': pt1, 'points':[]}
                    distances = {}
                    for key3, value3 in data.items():
                        if isinstance(key3, str) and "PZ" in key3 and key3 != key1:
                            for key4, value4 in value3.items():
                                if key4 == "coordinates":
                                    pt2 = (value4['longitude'], value4['latitude'])
                                    distance = haversine_distance(pt1, pt2)
                                    distances[key3] = (pt2, distance)
                    distances = dict(sorted(distances.items(), key=lambda x: x[1][1])[:2])
                    for key, value in distances.items():
                        poly_table[key1]['points'].append({'name': key,'coordinates':value[0], 'distance':value[1]})
    return poly_table

# def closest_poly_table(data):
#     poly_table = {}
#     for k1, v1 in data.items():
#         if(type(k1) == str and "PZ" in k1):
#             for k2, v2 in v1.items():
#                 if k2 == "coordinates":
#                     pt1 = (v2['longitude'], v2['latitude'])
#                     poly_table[k1]={'firstPoint': pt1, 'points':[]}
#                     distances = {}
#                     for k1_bis, v1_bis in data.items():
#                         if(type(k1_bis) == str and "PZ" in k1_bis and k1_bis != k1):
#                             for k2_bis, v2_bis in v1_bis.items():
#                                 if k2_bis == "coordinates":
#                                     pt2 = (v2_bis['longitude'], v2_bis['latitude'])
#                                     distance = haversine_distance(pt1, pt2)
#                                     distances[k1_bis] = (pt2, distance)
#                     distances = dict(sorted(distances.items(), key=lambda x: x[1][1])[:2])
#                     for key, value in distances.items():
#                         poly_table[k1]['points'].append({'name': key,'coordinates':value[0], 'distance':value[1]})
    
#     polygon_dict = {}
#     for k2, v2 in poly_table.items():
#         poly_name = k2
#         for k3, v3 in v2.items():
#             if k3 == 'points':
#                 for elt in v3:
#                     poly_name = poly_name + elt['name']
#         for k4, v4 in poly_table.items():
#             polygon_dict[poly_name] = Polygon([v4['firstPoint']] +
#                                               [p['coordinates'] for p in v4['points']])
        
#         if check_polygon_intersection_adjacent(polygon_dict):
#             coordinates ={}


        # poly_name = k + v['name']
        # print(str('[----> ') + str(k)+str(v))
    # print(poly_name)

    # polygon_dict = {}
    # for k, v in poly_table.items():
    #     polygon_dict[k] = Polygon([v['firstPoint']] + 
    #                               [p['coordinates'] for p in v['points']])
        
    # if check_polygon_intersection_adjacent(polygon_dict):
    #     coordinates = {}
    #     for k, v in poly_table.items():
    #         coordinates[k] = [v['firstPoint']]
    #         for point in v['points']:
    #             coordinates[k].append(point['coordinates'])
    #     return coordinates
    # else:
    #     return {}

def check_polygon_intersection_adjacent(poly_table):
    for k1, v1 in poly_table.items():
        first_polygon = Polygon([poly_table[k1]['firstPoint']] + [p['coordinates'] for p in poly_table[k1]['points']])
        for k2, v2 in poly_table.items():
            if k1 == k2:
                continue
            second_polygon = Polygon([v2['firstPoint']] + [p['coordinates'] for p in v2['points']])
            if first_polygon.intersects(second_polygon) or not first_polygon.touches(second_polygon):
                return False
    return True

def htmlCreator():
    folderFiles2 = searchFolder('\\test.html')
    folderFiles3 = searchFolder('\\out.geojson')

    currentFolder = os.getcwd()
    if os.path.isdir(currentFolder) == True:
        filename = frame(currentFolder) 
    
    data = {}
    df = pd.DataFrame()
    
    if os.path.isfile(filename) == True:
        df = pd.read_csv(filename)
    data = read_df(df)
    p_table = closest_poly_table(data)
    print(p_table)
    # return p_table
    
    # for k1, v1 in data.items():
    #     if(type(k1) == str and "PZ" in k1):
    #         for k2, v2 in v1.items():
    #             my_list = {}
    #             if k2 == "coordinates":
    #                 pt1 = (v2['longitude'], v2['latitude'])
    #                 poly_table[k1]={'firstPoint': pt1}
    #                 for k1_bis, v1_bis in data.items():
    #                     if(type(k1_bis) == str and "PZ" in k1_bis and k1_bis != k1):
    #                         for k2_bis, v2_bis in v1_bis.items():
    #                             ect = 0
    #                             if k2_bis == "coordinates":
    #                                 pt2 = (v2_bis['longitude'], v2_bis['latitude'])
    #                                 ect = haversine_distance(pt1, pt2)
    #                                 if ect < ect_temp or not my_list:
    #                                     if len(my_list) >= 2:
    #                                         my_list.pop(max(my_list, key=my_list.get))
    #                                     ect_temp = ect
    #                                     my_list[k1_bis] = pt2
    #                                     poly_table['Points'] = my_list


    
    
    # for k1, v1 in data.items():
    #     if(type(k1) == str and "PZ" in k1):
    #         for k2, v2 in v1.items():
    #             my_list = {}
    #             if k2 == "coordinates":
    #                 pt1 = (v2['longitude'], v2['latitude'])
    #                 poly_table[k1]={'firstPoint': pt1}
    #                 for k1_bis, v1_bis in data.items():
    #                     if(type(k1_bis) == str and 
    #                        "PZ" in k1 and 
    #                        k1_bis != k1):

    #                         for k2_bis, v2_bis in v1_bis.items():
    #                             ect = 0
    #                             print(v1_bis)
    #                             if k2_bis == "coordinates":
    #                                 pt2 = (v2_bis['longitude'], v2_bis['latitude'])
    #                                 ect = haversine_distance(pt1, pt2)
    #                                 if (my_list):
    #                                     print(len(my_list.keys()))
    #                                     if (ect < ect_temp):
    #                                         print(my_list)
    #                                 elif not my_list:
    #                                     ect_temp = ect
    #                                     my_list[k1_bis] = {}
                                        

                                    
                                    

                            


                
    

        

        
    # if os.path.isfile(filename)==True:
    #     df = pd.read_csv(filename) #df=>dataframe
    
    #     for row in df.itertuples():

    #         values = row[1].split(";")
    #         data[values[0]]={"Coordinates":{"Latitude": float(values[1]),
    #                                         "Longitude": float(values[2]),
    #                                         "Profondeur": float(values[3])}}

    #     for ks1, vs1 in data.items():
    #         list_pointname.append(ks1)
    #         grp_name.append('set1')
    #         geo_df_list[ks1]={}
    #         for ks2, vs2 in vs1.items():
    #             for ks3, vs3 in vs2.items():
    #                 if ks3 == "Latitude":
    #                     lat = vs3
    #                     list_latitude.append(vs3)
    #                 elif ks3 == "Longitude":
    #                     long = vs3
    #                     list_longitude.append(vs3)
    #                 elif  ks3 == "Profondeur":
    #                     prof = vs3
    #                     list_profondeur.append(vs3)
                        
    #             coord = {'coordinates': [long, lat]}
    #             pro = {'profondeur': prof}
    #             coord.update(pro)
    #             geo_df_list[ks1] = coord

    #     vals1=[] 
    #     for k1, v1 in geo_df_list.items():
    #         for k2, v2 in v1.items():
    #             if k2 == 'coordinates':
    #                 vals1.append(v2)
        
    #     number_of_polygon = int(np.ceil(len(vals1)/3)) #calcul du nombre de polygone
    #     LS = {}
    #     iLS = {}
    #     vals_temp = {}

    #     for index_polygon in range(number_of_polygon):
    #         vals_temp = create_vals2(vals1)
        
    #     j = 0
    #     k = 0
    #     for cle, val in vals_temp.items():
    #         if cle == 'coordinates':
    #             for elt in val:
    #                 LineS[k] = {}
    #                 for i in range(len(elt) - 1):
    #                     if j < 3:
    #                       LineS[k][j] = LineString([elt[i], elt[i + 1]])
    #                       j += 1
    #                     else:
    #                       k += 1
    #                       LineS[k] = {}
    #                       j = 0
    #                       m = 0
    #                       LineS[k][j] = LineString([elt[-1], elt[0]])
    #                       for l in range(len(elt) - 1):
    #                           if m < 3:
    #                               LineS[k][m] = LineString([elt[l], elt[l + 1]])
    #                               m += 1
        
    #     for k1, v1 in LineS.items():
    #         for k2, v2 in v1.items():
    #             for pt1, pt2 in zip(v2.coords,
    #                                 v2.coords[1:]):
    #                 LS[j] = v2
    #                 j += 1
                    
    #                 df1 = pd.DataFrame(
    #                     {'Site':list_pointname,
    #                      'Latitude':list_latitude,
    #                      'Longitude':list_longitude,
    #                      'Profondeur':list_profondeur,
    #                     })
                    
    #                 x = df1['Latitude'].tolist()
    #                 y = df1['Longitude'].tolist()
            
    #                 poly = Polygon(list(zip(y, x)))
    #                 geo_j = folium.GeoJson(data = poly,
    #                                         style_function = lambda x:{'color': 'blue'})
                
    #                 text_poly = str(poly.centroid)
                    
    #                 for r in text_poly.split(" "):
    #                     if '(' in r:
    #                         y_centroid = r.split('(')[1]
    #                     if ')' in r:
    #                         x_centroid = r.split(')')[0]
                        
    #                 y_centroid_around = np.around(float(y_centroid), decimals=6)
    #                 x_centroid_around = np.around(float(x_centroid), decimals=6)
            
    #                 gdf = gpd.GeoDataFrame(df1,
    #                                         geometry = gpd.points_from_xy((df1.Longitude), 
    #                                                                      (df1.Latitude)),
    #                                         crs='epsg:4326')
            
    #                 gdf_size = gdf['Site'].size
    #                 gdf_index = gdf_size + 1
    #                 gdf.loc[gdf_index - 1] = ['PZCentroid',
    #                                           x_centroid_around, 
    #                                           y_centroid_around, 
    #                                           '', 
    #                                           poly.centroid]
                    
    #             iLS = construct_iLSDict(k2, v2, gdf)
                
    #             coor = {'coordinates': [x_centroid_around, 
    #                                     y_centroid_around]}
                
    #             pr = {'profondeur':None}
    #             coor.update(pr)
    #             geo_df_list['PZCentroid'] = coor
                
    #             coordonnees = (geo_df_list[ks1]['coordinates'][1], 
    #                             geo_df_list[ks1]['coordinates'][0])
            
    #             m = folium.Map(location=coordonnees,
    #                             zoom_start=16,
    #                             zoom_control=(True))

    #             for index, rows in gdf.iterrows():
    #                 html=f"""<h1>{rows['Site']}</h1></n>
    #                 Coordonées (Lat, Long): {rows['Latitude'], rows['Longitude']}</br>
    #                 Profondeur mesurée (m): {rows['Profondeur']}</br>
    #                 Profondeur crépines (m):</br>
    #                 Diamètre int. (cm):</br>
    #                 Type:</br>
    #                 Masse d'eau:</br>
    #                 BDLISAv3:</br>"""
    #                 marker_coordinates = (rows['Latitude'], rows['Longitude'])
    #                 iframe = folium.IFrame(html=html, width=400, height=250)
    #                 popup = folium.Popup(iframe, max_width=500)
    #                 m.add_child(folium.Marker(location=marker_coordinates,
    #                                           popup = popup,
    #                                           tooltip = rows['Profondeur']
    #                                           )
    #                             )
    #         deltaZ = {}

    #     for k3, v3 in iLS.items():
    #         if (type(k3)== int):
    #             cle = k3
    #             print(str('k3---->') + str(k3) + str(v3))
    #             if (type(v3) != LineString):
    #                 print(v3)
    #                 delta = abs(v3['z1'] - v3['z2'])
    #                 deltaZ[str(cle) + '_' + 'ecart'] = delta
    #             # if(type(k3 == str)):
    #             #     print(v3)
    #             #     delta = abs(v3['z1'] - v3['z2'])
    #             #     deltaZ[str(cle) + '_' + 'ecart'] = delta
    
    #     iLS.update(deltaZ)

    #     for k, v in iLS.items():
    #         if type(k) == int:
    #             my_Values = v
    #             longueur = my_Values.length
    #             archi[k] = {'longueur': longueur,
    #                         'ecart':'',
    #                         'distance':'',
    #                         'valeur_de_depart':'',
    #                         'intervalle':''}

        
    #     for k, v in archi.items():
    #         for k2, v2 in iLS.items():
    #             rslt = 0
    #             print(v2)
    #             if (type(k2) == str and 'PZ' not in k2):
    #                 if (k == int(k2[0])):
    #                     archi[k]['ecart'] = v2
    #                     dist = archi[k]['longueur']/archi[k]['ecart']
    #                     archi[k]['distance'] = dist
                        
    #             if(type(k2) == str and 'PZ' in k2 and type(v2) == dict):
    #                 if (k == int(k2[0])):
    #                     archi[k]['valeur_de_depart'] = v2['z1']
    #                     rslt = v2['z1'] - v2['z2']
    #                     if (rslt > 0):
    #                         archi[k]['intervalle'] = -1
    #                     elif (rslt < 0):
    #                         archi[k]['intervalle'] = 1
                            
    #     print(str('archi---> ' +str(archi)))
        
    #     x_list = []
    #     y_list = []

    #     df2 = pd.DataFrame.from_dict(XY)
    #     for i in df2:
    #         print(i)
    #         k=0
    #         points = []
    #         temp = df2[i]['prof']
    #         pt1 = (df2[i]['long'], df2[i]['lat'])
    #         points.append(tuple([df2[i]['lat'], df2[i]['long']]))
    #         for j in df2:
    #             if j > i and temp == df2[j]['prof'] and k <= 3:
    #                 k += 1
    #                 pt2 = (df2[j]['long'], df2[j]['lat'])
    #                 points.append(tuple([df2[j]['lat'], df2[j]['long']]))
    #                 ect = haversine_distance(pt1, pt2)
    #                 print(str('ect---> ') + str(ect))

    #         for x in points:
    #             draw_polyligne(m, x) #goto draw polyligne module                  
           
    ##########################################
    # geo_j.add_to(m)
    # m.save(folderFiles2)    
    # web.open(folderFiles2)
    # gdf.to_file(folderFiles3)
    ##########################################
     
w = tk.Tk()
w.geometry("400x100")
w.grid_bbox(column=2, row=3)

def create_vals2(vals1):    
  # Vérifiez que les coordonnées forment un anneau fermé
  if (vals1[0] != vals1[-1]):
    vals1.append(vals1[0])

  # Créez un polygone à partir des coordonnées
  try:
    polygon = Polygon(vals1)
  except AssertionError:
    # Si l'erreur d'assertion se produit, cela signifie que les coordonnées sont dans le mauvais ordre
    # Inversez l'ordre des coordonnées et créez à nouveau le polygone
    vals1 = vals1[::-1]
    polygon = Polygon(vals1)

  # Créez un dictionnaire à partir du polygone
  data = {
    "type": "Polygon",
    "coordinates": [[c for c in ring.coords] for ring in list(polygon.interiors) + [polygon.exterior]]
  }
  return data    

def draw_polyligne(m, points):
    temp_tuple=tuple(points)
    m.add_child(folium.CircleMarker(location=temp_tuple,
                                    fill='true',
                                    radius=6,
                                    fill_color='green',
                                    color='clear',
                                    fill_opacity=1))
    
    folium.PolyLine(locations=points, 
                    color='blue',
                    weight=2.5, 
                    opacity=1).add_to(m) 

def geo_json(lat, lon, value, step):
    return {
      "type": "FeatureCollection",
      "features": [
        {
          "type": "Feature",
          "properties": {
            'color': 'black',
            'weight': 1,
            'fillColor': 'green',
            'fillOpacity': 0.5,
          },
          "geometry": {
            "type": "Polygon",
            "coordinates": [
                [
                    [lon - step, lat - step],
                    [lon - step, lat + step],
                    [lon + step, lat + step],
                    [lon + step, lat - step],
                    [lon - step, lat - step],
                ]          
                            ]      
                        }
          }         ]    
      }

def construct_iLSDict(k, LS, gdf):
    i = k
    idxMax = int(gdf.index.max(-1))
    iLS_Dict[k] = LS
    iLS = {}
    if (i <= idxMax):
        chaine = str(gdf['Site'][i])
        chaine_plus_1 = str(gdf['Site'][i + 1])
        if('PZCentroid' not in chaine):
            if ('PZCentroid' not in chaine_plus_1):
                iLS_name = str(i) + "_" + gdf['Site'][i] + "_" + gdf['Site'][i + 1]
                iLS[iLS_name] = {'z1': gdf['Profondeur'][i],
                                      'z2': gdf['Profondeur'][i + 1]}
            elif ('PZCentroid' == str(chaine_plus_1)):
                iLS_name = str(i) + "_" + gdf['Site'][i] + "_" + gdf['Site'][0]
                iLS[iLS_name] = {'z1': gdf['Profondeur'][i],
                                      'z2': gdf['Profondeur'][0]}
        iLS_Dict.update(iLS)
        i += 1
    
    return iLS_Dict

def frame(initDir):
  
    buttonKO = tk.Button(w, text="Quitter", width=10)
    buttonKO.grid(sticky="e", column=2, row=1, ipadx=5, ipady=5, padx=5, pady=5)
    if os.path.isdir(initDir) == True:
        filename = file_open(initDir)
        
    return filename
    w.mainloop()
    
def file_open(initDir):
    frame
    filename = fd.askopenfilename(initialdir=initDir, 
                                  filetypes = (
                                            ('csv files', '*.csv'),
                                            ('All files', '*.*')))
    return filename


htmlCreator()
w.destroy()