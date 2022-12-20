# -*- coding: utf-8 -*-

"""
Spyder Editor

This is a temporary script file.
"""
import sys
import os
import fiona
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import folium
from folium.features import DivIcon
from shapely import geometry
from shapely.geometry import Polygon, Point, LineString
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

def htmlCreator():
    folderFiles = os.path.join(os.getcwd(), 'files') # définition du dossier contenant les fichiers
    folderFiles2 = folderFiles + "\\test.html"
    folderFiles3 = folderFiles + "\\out.geojson"


    currentFolder = os.getcwd() # définition du dossier courant
    if os.path.isdir(currentFolder)==True:
        filename = frame(currentFolder) #filename contient en retour le chemin d'accès au fichier devant être lu
    
    data = {}
    df = pd.DataFrame()
        
    if os.path.isfile(filename)==True:
        df = pd.read_csv(filename) #df=>dataframe
    
        for row in df.itertuples():
            
            values = row[1].split(";")
            data[values[0]]={"Coordinates":{"Latitude":float(values[1]),
                                            "Longitude":float(values[2]),
                                            "Profondeur":float(values[3])}}
        # list_pointname =[] #Liste
        # list_latitude = []
        # list_longitude = []
        # list_profondeur = []
        # grp_name = []
        # geo_df_list = {}
        # length = {}
        # ecart = {}
        # archi = {}
        # XY = {}
        # LS_inter = {}

        for ks1, vs1 in data.items():
            list_pointname.append(ks1)
            grp_name.append('set1')
            geo_df_list[ks1]={}
            for ks2, vs2 in vs1.items():
                for ks3, vs3 in vs2.items():
                    if ks3 == "Latitude":
                        lat = vs3
                        list_latitude.append(vs3)
                    elif ks3 == "Longitude":
                        long = vs3
                        list_longitude.append(vs3)
                    elif  ks3 == "Profondeur":
                        prof = vs3
                        list_profondeur.append(vs3)
                        
                coord = {'coordinates': [long, lat]}
                pro = {'profondeur': prof}
                coord.update(pro)
                geo_df_list[ks1] = coord

        vals1=[] 
        for k1, v1 in geo_df_list.items():
            for k2, v2 in v1.items():
                if k2 == 'coordinates':
                    vals1.append(v2)
        
        number_of_polygon = int(np.ceil(len(vals1)/3)) #calcul du nombre de polygone
        print(number_of_polygon)
        vals2 = {}
        LS = {}
        iLS = {}
        vals_temp = {}
        ##################################################
        for i in range(0, number_of_polygon, 1):
            vals_temp = create_vals2(vals1, i) #goto create_vals2 module
            vals2[i] = vals_temp
            line = LineString(vals2[i])
            j = 0
            for pt1, pt2 in zip(line.coords, 
                                line.coords[1:]):
                LS[j] = LineString([pt1, pt2])
                j += 1

            df1 = pd.DataFrame(
                {'Site':list_pointname,
                 'Latitude':list_latitude,
                 'Longitude':list_longitude,
                 'Profondeur':list_profondeur,
                })
            
            x = df1['Latitude'].tolist()
            y = df1['Longitude'].tolist()
        
            poly = Polygon(list(zip(y, x)))
            geo_j = folium.GeoJson(data = poly,
                                   style_function = lambda x:{'color': 'blue'})
            
            text_poly = str(poly.centroid)
            for r in text_poly.split(" "):
                if '(' in r:
                    y_centroid = r.split('(')[1]
                if ')' in r:
                    x_centroid = r.split(')')[0]
                    
            y_centroid_around = np.around(float(y_centroid), decimals=6)
            x_centroid_around = np.around(float(x_centroid), decimals=6)
            
            gdf = gpd.GeoDataFrame(df1,
                                   geometry=gpd.points_from_xy((df1.Longitude), 
                                                                (df1.Latitude)),
                                   crs='epsg:4326')
            
            gdf_size = gdf['Site'].size
            gdf_index = gdf_size + 1
            gdf.loc[gdf_index - 1] = ['PZCentroid',
                                      x_centroid_around, 
                                      y_centroid_around, 
                                      '', 
                                      poly.centroid]
        
            for k, Line_String in LS.items():
                iLS = construct_iLSDict(k, Line_String, gdf)
    
            coor = {'coordinates': [x_centroid_around, 
                                    y_centroid_around]}
            pr = {'profondeur':None}
            coor.update(pr)
            geo_df_list['PZCentroid'] = coor
            
            coordonnees = (geo_df_list[ks1]['coordinates'][1], 
                           geo_df_list[ks1]['coordinates'][0])
    
            m = folium.Map(location=coordonnees,
                           zoom_start=16,
                           zoom_control=(True))

            for index, rows in gdf.iterrows():
                html=f"""<h1>{rows['Site']}</h1></n>
                Coordonées (Lat, Long): {rows['Latitude'], rows['Longitude']}</br>
                Profondeur mesurée (m): {rows['Profondeur']}</br>
                Profondeur crépines (m):</br>
                Diamètre int. (cm):</br>
                Type:</br>
                Masse d'eau:</br>
                BDLISAv3:</br>"""
                marker_coordinates = (rows['Latitude'], rows['Longitude'])
                iframe = folium.IFrame(html=html, width=400, height=250)
                popup = folium.Popup(iframe, max_width=500)
                m.add_child(folium.Marker(location=marker_coordinates,
                                          popup = popup,
                                          tooltip = rows['Profondeur']
                                          )
                            )
            deltaZ = {}
            ##################################################
            for key, values in iLS.items():
                if type(key)== int:
                    cle = key
                if(type(key) == str):
                    delta = abs(values['z1'] - values['z2'])
                    deltaZ[str(cle) + '_' + 'ecart'] = delta

        iLS.update(deltaZ)

        for k, v in iLS.items():
            if type(k) == int:
                my_Values = v
                
                longueur = my_Values.length
                archi[k] = {'longueur': longueur,
                            'ecart':'',
                            'distance':'',
                            'valeur_de_depart':'',
                            'intervalle':''}
        
        for k, v in archi.items():
            for k2, v2 in iLS.items():
                rslt=0

                if (type(k2) == str and 'PZ' not in k2):
                    if (k == int(k2[0])):
                        archi[k]['ecart'] = v2
                        dist = archi[k]['longueur']/archi[k]['ecart']
                        archi[k]['distance'] = dist
                        
                if(type(k2) == str and 'PZ' in k2 and type(v2) == dict):
                    if (k == int(k2[0])):
                        archi[k]['valeur_de_depart'] = v2['z1']
                        rslt = v2['z1'] - v2['z2']
                        if (rslt > 0):
                            archi[k]['intervalle'] = -1
                        elif (rslt < 0):
                            archi[k]['intervalle'] = 1
        
        i = 0
        for k, v in iLS.items():
            x = 0
            y = 0
            for k2, v2 in archi.items():
                prof_du_point_pour_0 = 0
                if k == k2:
                    line = v
                    tot_length = float(line.length)
                    prof_du_point = float(archi[k2]['valeur_de_depart'])
                    distance = 0
                    while distance < tot_length:
                        if (distance > 0):

                            name = 'inter: ' + str(i)
                            XY[name] = {'lat':[],
                                        'long':[],
                                        'prof':[]}
                            prof_du_point += float(archi[k2]['intervalle'])
                            new_point = line.interpolate(distance)
                            lat = np.around(float(new_point.xy[1][0]),
                                            decimals=6)
                            long = np.around(float(new_point.xy[0][0]),
                                              decimals=6)
                            
                            marker_coordinates = (lat, long) #inversion entre le html 
                            marker_web = (long, lat) #et le web
                            point_coordinates = Point(marker_web, precision=6)
                            #===================================
                            gdf_size = gdf['Site'].size
                            gdf_index = gdf_size + 1
                            gdf.loc[gdf_index - 1] = [name, lat, long,
                                                      prof_du_point,
                                                      point_coordinates]
                            XY[name]['lat'] = lat
                            XY[name]['long'] = long
                            XY[name]['prof'] = prof_du_point
                            #===================================   
                            
                            html_inter = f"""<h1>{name}</h1></n>
                            Coordonnées (Lat, Long): {lat, long}</br>
                            Profondeur théorique (m): {prof_du_point}</br>"""
                            
                            iframe_inter = folium.IFrame(html = html_inter,
                                                         width = 400,
                                                         height = 250)
                            
                            popup = folium.Popup(iframe_inter,
                                                 max_width = 500)
                            
                            geo_j.add_child(folium.Marker(location=marker_coordinates,
                                                          popup=popup,
                                                          tooltip=prof_du_point,
                                                          icon=folium.Icon(icon='glyphicon-screenshot',
                                                                           icon_color='#9b59b6',
                                                                           color='green'
                                                                           )
                                                          )
                                            )
                            i += 1
                        elif (distance == 0):
                            name = 'inter: ' + str(i)
                            XY[name] = {'lat':[],
                                        'long':[],
                                        'prof':[]}
                            prof_du_point_pour_0 = prof_du_point
                            new_point = line.interpolate(distance)
                            lat = np.around(float(new_point.xy[1][0]),
                                            decimals=6)
                            long = np.around(float(new_point.xy[0][0]),
                                            decimals=6)
                            
                            XY[name]['lat'] = lat
                            XY[name]['long'] = long
                            XY[name]['prof'] = prof_du_point_pour_0
                            # i += 1

                        distance += float(v2['distance'])

        ##########################################
        x_list = []
        y_list = []

        df2 = pd.DataFrame.from_dict(XY)
        for i in df2:
            k=0
            points = []
            temp = df2[i]['prof']
            pt1 = (df2[i]['long'], df2[i]['lat'])
            points.append(tuple([df2[i]['lat'], df2[i]['long']]))
            for j in df2:
                if j > i and temp == df2[j]['prof'] and k <= 3:
                    k += 1
                    pt2 = (df2[j]['long'], df2[j]['lat'])
                    points.append(tuple([df2[j]['lat'], df2[j]['long']]))

            for x in points:
                draw_polyligne(m, x) #goto draw polyligne module                  
                    
    ##########################################
    geo_j.add_to(m)
    m.save(folderFiles2)    
    web.open(folderFiles2)
    gdf.to_file(folderFiles3)
    ##########################################
     
w = tk.Tk()
w.geometry("400x100")
w.grid_bbox(column=2, row=3)

def create_vals2(vals1, i):
    vals2 = []

    for x in vals1:
        if (i <= 3):
            vals2.append(x)
            #Si la valeur de len(vals2) est supérieure à 3 alors il faut
            #créer une autre liste "vals2" qui va permettre de créer 
            #un autre polygone qui devra être accolé au premier
        if (i > 3):
            vals2.append(x)
    return vals2
        
        # if (i > 3):
        #     vals2 = []
    
def draw_polyligne(m, points):
    temp_tuple=tuple(points)
    m.add_child(folium.CircleMarker(location=temp_tuple,
                                    fill='true',
                                    radius=6,
                                    fill_color='green',
                                    color='clear',
                                    fill_opacity=1))
    
    # folium.PolyLine(locations=points, 
    #                 color='blue',
    #                 weight=2.5, 
    #                 opacity=1).add_to(m) 

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