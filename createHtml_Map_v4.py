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
import haversine
import geopandas as gpd
import matplotlib.pyplot as plt
import folium
import haversine
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

    # Convert degrees to radians
    lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])

    # Calculate differences
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Apply haversine formula
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))

    # Earth's radius (mean) in kilometers
    R = 6371

    return c * R

def closest_poly_table(data):
    poly_table = {}
    for key1, value1 in data.items():
        pt1 = (value1['coordinates']['longitude'], value1['coordinates']['latitude'])
        poly_table[key1] = {'coordinates':(value1['coordinates']['longitude'], 
                                           value1['coordinates']['latitude']), 
                            'profondeur':value1['complements']['profondeur'],
                            'points':[]}
        for key2, value2 in data.items():
            if key2 != key1:
                pt2 = (value2['coordinates']['longitude'], 
                       value2['coordinates']['latitude'])
                poly_table[key1]['points'].append((key2, pt2))
                
    return poly_table

def keep_adjacent_polygons(polygons):
    # Utilisez la fonction de distance de point à point pour trouver les points les plus proches pour chaque point donné.
    # Utilisez un dictionnaire pour stocker les informations de chaque polygone, y compris les coordonnées du premier point, les points les plus proches.
    # Parcourez les polygones et vérifiez s'il y a des recouvrements en utilisant la distance entre les points.
    # Si un recouvrement est détecté, ne pas ajouter le polygone au dictionnaire de sortie.
    my_dict = {}
    output_list = []
    for polygon_key, polygon_data in polygons.items():
        my_list = []
        for k, v in polygon_data.items():
            if k == 'coordinates':
                pt1 = v
                my_dict[polygon_key] = {'coordinates': v, 
                                        'profondeur': '',
                                        'points': ''}
            if k == 'profondeur':
                my_dict[polygon_key]['profondeur'] = v
                
            if k == 'points':
                distances = []
                for elt in v:
                    pt2_temp = (elt[1][0], elt[1][1])
                    name = elt[0]
                    pt2 = pt2_temp
                    distance = haversine_distance(pt1, pt2)
                    distances.append((distance, name, pt2))
                distances = sorted(distances, key=lambda x: x[0])[:2]
                my_list.append([(x[1], x[2]) for x in distances])
                my_dict[polygon_key]['points'] = my_list
    return my_dict

def remove_overlapping_polygons(polygons):
    output_dict = {}
    processed_points = []
    # print("polygons---> " + str(polygons))
    for polygon_key, polygon_data in polygons.items():
        pt1_1 = polygon_data['coordinates']
        pt2_1 = polygon_data['points'][0][0][1]
        pt3_1 = polygon_data['points'][0][1][1]
        polygon1 = Polygon([pt1_1, pt2_1, pt3_1])
        is_overlapping = False
        if not output_dict:
            output_dict[polygon_key] = polygon_data

        for polygon_key_2, polygon_data_2 in polygons.items():
            # print(polygon_data_2)
            if polygon_key != polygon_key_2:
                pt1_2 = polygon_data_2['coordinates']
                pt2_2 = polygon_data_2['points'][0][0][1]
                pt3_2 = polygon_data_2['points'][0][1][1]
                polygon2 = Polygon([pt1_2, pt2_2, pt3_2])
                if polygon2.covers(polygon1) or polygon1.covers(polygon2):
                    is_overlapping = True
        if not is_overlapping:
            output_dict[polygon_key] = polygon_data
    # print('output-dict---> ' + str(output_dict))
    return output_dict

def extract_centroid(centroid):
    point_string = str(centroid)
    point_parts = point_string.split(' ')
    latitude = float(point_parts[2][:-1])
    longitude = float(point_parts[1][1:])
    coordinates = (latitude, longitude)
    return coordinates

def create_Map(p_non_overlapping):
    output_polygon = []
    output_polygon_dict = {}
    first_key = next(iter(p_non_overlapping)) #Récupère la première clé du dict
    
    for k, v in p_non_overlapping.items():
        if k == first_key:
            coordonnees = (v['coordinates'][1],
                           v['coordinates'][0])
            
            m = folium.Map(location=coordonnees,
                            zoom_start=16,
                            zoom_control=(True),
                            tiles = 'OpenStreetMap')
    return m

def create_PointsPolygons(p_non_overlapping, m):
    output_polygon_dict =  {}
    output_polygon = []
    points = []
    Dict = {}
    for k, v in p_non_overlapping.items():
        i = 1
        pt1 = (v['coordinates'][0],
                v['coordinates'][1])
        
        points.append(pt1)
        name = k
        ls_Dict = {}
        profondeur_pt1 = float(p_non_overlapping[k]['profondeur'])
        for k1, v1 in v.items():
            if k1 == 'profondeur':
                Dict[k] = v1
            if k1 == 'points':
                for elt in v1:
                    pt2 = ((elt[0][1][0]), (elt[0][1][1]))
                    points.append(pt2)
                    ls1 = LineString([pt1, pt2])

                    profondeur_pt2 = float(p_non_overlapping[(elt[0][0])]['profondeur'])
                    ecart_pt1pt2 = profondeur_pt1 - profondeur_pt2
                    pt3 = ((elt[1][1][0]), (elt[1][1][1]))
                    points.append(pt3)
                    profondeur_pt3 = float(p_non_overlapping[(elt[1][0])]['profondeur'])
                    ecart_pt2pt3 = profondeur_pt2 - profondeur_pt3
                    ecart_pt3pt1 = profondeur_pt3 - profondeur_pt1
                    ls2 = LineString([pt2, pt3])
                    ls3 = LineString([pt3, pt1])
                    
                    name = name + '+' + str(v1[0][0][0]) + '+' + str(v1[0][1][0])
                    ls_Dict['LS_' + name] = [ls1, ls2, ls3]
                    ls_Dict['LS_' + name + '_ecart'] = [ecart_pt1pt2, 
                                                        ecart_pt2pt3, 
                                                        ecart_pt3pt1]
                    
                    
                    my_Polygon = Polygon([pt1, pt2, pt3])
                    is_overlapping = False
                    if output_polygon:
                        for p in output_polygon:
                            if (my_Polygon.overlaps(p)) or (p.overlaps(my_Polygon)):
                                is_overlapping = True
                                
                        if not is_overlapping:
                            output_polygon.append(my_Polygon)
                            output_polygon_dict[name] = my_Polygon
                            output_polygon_dict['LS_' + name] = ls_Dict['LS_' + name]
                            output_polygon_dict['LS_' + name + '_ecart'] = ls_Dict['LS_' + name + '_ecart']
                    else:
                        output_polygon.append(my_Polygon)
                        output_polygon_dict[name] = my_Polygon
                        output_polygon_dict['LS_' + name] = ls_Dict['LS_' + name]
                        output_polygon_dict['LS_' + name + '_ecart'] = ls_Dict['LS_' + name + '_ecart']
                                    
    return output_polygon_dict, points, Dict

def create_html(k, v):
    html1 = """Polygon: """ + k + """
    </br>Area: """ + str(round(v.area, 6))
    return html1

def create_iframe(html):
    return folium.IFrame(html = html,
                            width = 200,
                            height = 125)

def create_popup(iframe):
    return folium.Popup(iframe, max_width = 250)
    
def create_marker(centroid, popup, geo_p):
    return folium.Marker(location = centroid, 
                         popup = popup).add_to(geo_p)
    
def htmlCreator():
    folderFiles2 = searchFolder('\\test.html')
    folderFiles3 = searchFolder('\\out.geojson')
    existing_polylines = []

    currentFolder = os.getcwd()
    if os.path.isdir(currentFolder) == True:
        filename = frame(currentFolder) 
    
    data = {}
    df = pd.DataFrame()
    
    if os.path.isfile(filename) == True:
        df = pd.read_csv(filename)
        
    data = read_df(df) #GoTo 'read_df' with a dataframe argument and return a data dict with this information 'latitude', 'longitude' and 'profondeur'
    p_table = closest_poly_table(data) #GoTo 'closest_poly_table' with data and return dict
    p_adjacent = keep_adjacent_polygons(p_table)
    p_non_overlapping = remove_overlapping_polygons(p_adjacent)

    output_polygon = []
    first_key = next(iter(p_non_overlapping)) #Récupère la première clé du dict
    
    m = create_Map(p_non_overlapping)
    output_polygon_dict, points, Dict = create_PointsPolygons(p_non_overlapping, m)

    for k_poly, v_poly in output_polygon_dict.items():
        processed_polylines = []
        if '_' not in k_poly:
            geo_p = folium.GeoJson(data = v_poly,
                                   style_function = (lambda x: {'fillColor': 'blue', 
                                                                'weight':0.5}))
            geo_p.add_to(m)
            
            centroid = extract_centroid(v_poly.centroid)
            html1 = create_html(k_poly, v_poly)
            iframe1 = create_iframe(html1)
            popup1 = create_popup(iframe1)
            marker1 = create_marker(centroid, popup1, geo_p)
        
        elif '_' in k_poly:
            if ('ecart' not in k_poly and 'LS' in k_poly):
                i = 0
                for elt in v_poly:
                    i += 1
                    line = elt
                    line.coords = [(y, x) for x, y in line.coords]
                    fp = folium.PolyLine(locations = line.coords[:],
                                         color = '#FF1100',
                                         dash_array = '10, 10',
                                         weight = 2)
                    if fp.locations not in processed_polylines:
                        processed_polylines.append(fp.locations)
                        fp.add_to(m)
                    else:
                        v_poly.pop[i]
    print(output_polygon_dict)
                        
    #         cle = k
    #         for k2, v2 in output_polygon_dict.items():
    #             if '_' in k2 and cle in k2 and 'LS' in k2:
    #                 if ('ecart' not in k2):
    #                     i = 0
    #                     for elt in v2:
    #                         i += 1
    #                         line = elt
    #                         line.coords = [(y, x) for x, y in line.coords]
    #                         fp = folium.PolyLine(locations = line.coords[:], 
    #                                         color = '#FF1100', 
    #                                         dash_array = '10, 10',
    #                                         weight = 2)
    #                         if fp.locations not in existing_polylines:
    #                             existing_polylines.append(fp.locations)
    #                             fp.add_to(m)
    #                         else:
    #                             v2.pop(i)
                                
                        
                                


            
                            # for k3, v3 in output_polygon_dict.items():
                            #     tot_length = 0
                            #     i = 0
                            #     print(k3)
                            #     if '_' in k3 and cle in k3 and 'LS' in k3:
                            #         if ('ecart' in k3):
                            #             c = str(k3.split('_')[1][0:3])
                            #             #Valable uniquement pour des point nommés de 3 caractères
                            #             depthValue = Dict[c]
                            #             ecart_between_points = float(line.length)/float(abs(v3[i]))
                                        
                            #             while float(tot_length) < float(line.length):
                            #                 if tot_length > 0:
                            #                     if (v3[i]) > 0:
                            #                         depthValue += -1
                            #                     elif (v3[i]) < 0:
                            #                         depthValue += +1
                                                    
                            #                     html2 = """Profondeur: """ + str(depthValue)
                            #                     iframe2 = folium.IFrame(html = html2,
                            #                                             width = 200,
                            #                                             height = 125)
                            #                     popup2 = folium.Popup(iframe2, max_width = 250)
                            #                     new_point = line.interpolate(tot_length, False)
                            #                     coordins = (new_point.x, new_point.y)
                            #                     i_con = folium.Icon(color = 'green')
                            #                     folium.Marker(location = coordins,
                            #                                   icon = i_con, 
                            #                                   popup = popup2).add_to(geo_p)
                                                
                            #                 tot_length += ecart_between_points
                            #             i += 1
                                    
            
                    

                            # for k3, v3 in output_polygon_dict.items():
                            #     tot_length = 0
                            #     if '_' in k3 and k2 in k3 and 'LS' in k3:
                            #         if('ecart' in k3):

                            #             length_of_line = line.length
                            #             ecart_between_points = length_of_line/float(abs(v3[i]))
                            #             if tot_length > 0  or tot_length < length_of_line:
                            #                 while tot_length < length_of_line:
                            #                     new_point = line.interpolate(tot_length, False)
                            #                     coordins = (new_point.x, new_point.y)
                            #                     i_con = folium.Icon(color='green')
                            #                     folium.Marker(location = coordins, 
                            #                                   icon = i_con).add_to(geo_p)
                            #                     tot_length += ecart_between_points
                                
                            #             i += 1
   





                                

            # for elt in points:

            #     html = """Polygône: """ + k + """
            #     </br>
            #     Coordonnées: """ + str(v) + """
            #     </br>
            #     Aire : """ + str(v.area) + """
            #     </br>
            #     Coord. centroïd: """ + str(centroid) + """
            #     </br>
            #     """ + str(v.boundary.wkt)
            #     iframe = folium.IFrame(html = html, 
            #                             width = 400, 
            #                             height = 250)
            #     popup = folium.Popup(iframe, max_width = 500)
            #     coord = (elt[1], elt[0])
            #     folium.Marker(location = coord,  
            #                   popup = popup).add_to(geo_p)

    

    


    

    # for k1, v1 in p_non_overlapping.items():
    #     pt1 = (v1['coordinates'][0],
    #             v1['coordinates'][1])
    #     name = k1
    #     # is_overlapping = False
    #     for k2, v2 in v1.items():
    #         if k2 == 'points':
    #             for elt in v2:
    #                 pt2 = ((elt[0][1][0]), (elt[0][1][1]))
    #                 pt3 = ((elt[1][1][0]), (elt[1][1][1]))
    #                 name = name + str(v2[0][0][0]) + str(v2[0][1][0])
    #                 my_Polygon = Polygon([pt1, pt2, pt3])
    #                 is_overlapping = False
    #                 if output_polygon:
    #                     for p in output_polygon:
    #                         # print(type(p), type(my_Polygon))
    #                         if (my_Polygon.overlaps(p)) or (p.overlaps(my_Polygon)):
    #                             is_overlapping = True
                                
    #                     if not is_overlapping:
    #                         output_polygon.append(my_Polygon)
    #                         output_polygon_dict[name] = my_Polygon
    #                 else:
    #                     output_polygon.append(my_Polygon)
    #                     output_polygon_dict[name] = my_Polygon
                
    #                 # geo_p = folium.GeoJson(data=my_Polygon, style_function=(lambda x: {'fillColor': 'orange'}))
    #                 # centroid = extract_centroid(my_Polygon.centroid)
    #                 # folium.Marker(location=centroid, tooltip=name).add_to(geo_p)
    #                 # geo_p.add_to(m)
                    
    # # print(output_polygon_dict)
    #         for k, v in output_polygon_dict.items():
                
    #             geo_p = folium.GeoJson(data = v, 
    #                                    style_function = (lambda x: {'fillColor': 'orange'}))
    #             centroid = extract_centroid(v.centroid)
    #             html = """Polygône: """ + k + """
    #             </br>
    #             Coordonnées: """ + str(v) + """
    #             </br>
    #             Aire : """ + str(v.area) + """
    #             </br>
    #             Coord. centroïd: """ + str(centroid) + """
    #             </br>
    #             """ + str(v.boundary.wkt)
                
    #             print(output_polygon_dict['PZ1PZ2PZ5'])
                
    #             Line_String = v.boundary
    #             print(Line_String.length)
    #             pt = Line_String.interpolate(0.01, True)
    #             print(pt.bounds)
                
                

            
                
                
                
                
    #             iframe = folium.IFrame(html = html, 
    #                                    width = 400, 
    #                                    height = 250)
    #             popup = folium.Popup(iframe, max_width = 500)
    #             folium.Marker(location = centroid, 
    #                           tooltip = k, 
    #                           popup = popup).add_to(geo_p)
    #             geo_p.add_to(m)
            

    
    
                
                # for k2, v2 in v1.items():
                #     if k2 == 'points':
                #         for elt in v2:
                #             is_overlapping = False
                #             pt2 = ((elt[0][1][0]), (elt[0][1][1]))
                #             pt3 = ((elt[1][1][0]), (elt[1][1][1]))
                #             name = name + str(v2[0][0][0]) + str(v2[0][1][0])
                #             my_Polygon = Polygon([pt1, pt2, pt3])
                #             if output_polygon:
                #                 output_polygon.append(my_Polygon)
                #             else:
                #                 for elt in output_polygon:
                #                     if (my_Polygon.covers(elt)) or (elt.covers(my_Polygon)):
                #                         is_overlapping = True
                            
                #                 if not is_overlapping:
                #                     output_polygon.append(my_Polygon)
                                        
                #             for elt in output_polygon:
                #                 geo_p = folium.GeoJson(data = elt,
                #                                        style_function=(lambda x: {'fillColor': 'orange'})
                #                                        )
                            
                #                 centroid = extract_centroid(elt.centroid)
                            
                #                 folium.Marker(location = centroid,
                #                               tooltip = name).add_to(geo_p)

                            
                #             geo_p.add_to(m)
                    

    
        
        
        

        
            

        


    
    

    
    
    #p_non_overlapping---> {
    #'PZ1': {
    #           'coordinates': (4.989822, 45.774408), 
    #           'profondeur': 10.0, 
    #           'points': [[('PZ2', (4.989062, 45.780709)), 
    #                       ('PZ3', (4.980835, 45.778293))]]}, 
    #'PZ4': {
    #           'coordinates': (4.958782, 45.770784), 
    #           'profondeur': 13.0, 
    #           'points': [[('PZ3', (4.980835, 45.778293)), 
    #                       ('PZ1', (4.989822, 45.774408))]]}}
    
    

    ######################################
    # for k1, v1 in p_adjacent.items():
    #     if(type(k1) == str and 'PZ' in k1):
    #         for k2, v2 in v1.items():
    #             my_list = {}
    #             if (k2 == 'first_point'):
                    
            


    
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
    #     # df = pd.read_csv(filename) #df=>dataframe
    
    #     for row in df.itertuples():

    #         values = row[1].split(";")
    #         data[values[0]]={"Coordinates":{"Latitude": float(values[1]),
    #                                         "Longitude": float(values[2]),
    #                                         "Profondeur": float(values[3])}}
            
    #         print('data' + str(data))


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
    #                       'Latitude':list_latitude,
    #                       'Longitude':list_longitude,
    #                       'Profondeur':list_profondeur,
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
    #                                                                       (df1.Latitude)),
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
                
                # coordonnees = (geo_df_list[ks1]['coordinates'][1], 
                #                 geo_df_list[ks1]['coordinates'][0])
            
                # m = folium.Map(location=coordonnees,
                #                 zoom_start=16,
                #                 zoom_control=(True))

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
           
    # #########################################
    # geo_j.add_to(m)
    m.save(folderFiles2)    
    web.open(folderFiles2)
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