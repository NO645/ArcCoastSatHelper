import urllib.request
import numpy as np
import pickle
import matplotlib.pyplot as plt
from osgeo import gdal
import os, glob, shutil, subprocess, copy
import shapefile
import scipy.stats as st
import datetime
from scipy.spatial.transform import Rotation as R
from shapely.geometry import LineString
from scipy.spatial.distance import directed_hausdorff
import matplotlib.image as mpimg

import time
import numpy as np
import pickle
import warnings
from coastsat import SDS_download, SDS_preprocess, SDS_shoreline, SDS_tools, SDS_transects, SDS_slope
import shapefile
import datetime
import pandas as pd
import pytz
import cc_functionsv2
from pyproj import CRS
import os
import copy
import glob

import matplotlib as mpl
import matplotlib.cm as cm










def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))





def compute_intersection(output, transects):
    """
    Computes the intersection between the 2D shorelines and the shore-normal.
    transects. It returns time-series of cross-shore distance along each transect.
    
    KV WRL 2018       

    Arguments:
    -----------
    output: dict
        contains the extracted shorelines and corresponding metadata
    transects: dict
        contains the X and Y coordinates of each transect
    settings: dict with the following keys
        'along_dist': int
            alongshore distance considered caluclate the intersection
              
    Returns:    
    -----------
    cross_dist: dict
        time-series of cross-shore distance along each of the transects. 
        Not tidally corrected.        
    """    
    
    # loop through shorelines and compute the median intersection    
    intersections = np.zeros((len(output['shorelines']),len(transects)))
    for i in range(len(output['shorelines'])):

        sl = output['shorelines'][i]
        
        for j,key in enumerate(list(transects.keys())): 
            
            # compute rotation matrix
            X0 = transects[key][0,0]
            Y0 = transects[key][0,1]
            temp = np.array(transects[key][-1,:]) - np.array(transects[key][0,:])
            phi = np.arctan2(temp[1], temp[0])
            Mrot = np.array([[np.cos(phi), np.sin(phi)],[-np.sin(phi), np.cos(phi)]])
    
            # calculate point to line distance between shoreline points and the transect
            p1 = np.array([X0,Y0])
            p2 = transects[key][-1,:]
            d_line = np.abs(np.cross(p2-p1,sl-p1)/np.linalg.norm(p2-p1))
            # calculate the distance between shoreline points and the origin of the transect
            d_origin = np.array([np.linalg.norm(sl[k,:] - p1) for k in range(len(sl))])
            # find the shoreline points that are close to the transects and to the origin
            # the distance to the origin is hard-coded here to 1 km 
            idx_dist = np.logical_and(d_line <= 25, d_origin <= 1000)
            # find the shoreline points that are in the direction of the transect (within 90 degrees)
            temp_sl = sl - np.array(transects[key][0,:])
            phi_sl = np.array([np.arctan2(temp_sl[k,1], temp_sl[k,0]) for k in range(len(temp_sl))])
            diff_angle = (phi - phi_sl)
            idx_angle = np.abs(diff_angle) < np.pi/2
            # combine the transects that are close in distance and close in orientation
            idx_close = np.where(np.logical_and(idx_dist,idx_angle))[0]     
            
            # in case there are no shoreline points close to the transect 
            if len(idx_close) == 0:
                intersections[i,j] = np.nan
            else:
                # change of base to shore-normal coordinate system
                xy_close = np.array([sl[idx_close,0],sl[idx_close,1]]) - np.tile(np.array([[X0],
                                   [Y0]]), (1,len(sl[idx_close])))
                xy_rot = np.matmul(Mrot, xy_close)
                # compute the median of the intersections along the transect
                intersections[i,j] = np.nanmedian(xy_rot[0,:])
    
    # fill the a dictionnary
    cross_dist = dict([])
    for j,key in enumerate(list(transects.keys())): 
        cross_dist[key] = intersections[:,j]   
    
#    # save a .csv file for Excel users
#    out_dict = dict([])
#    out_dict['dates'] = output['dates']
#    for key in transects.keys():
#        out_dict['Transect '+ key] = cross_dist[key]
#    df = pd.DataFrame(out_dict)
#    fn = os.path.join(settings['inputs']['filepath'],settings['inputs']['sitename'],
#                      'transect_time_series.csv')
#    df.to_csv(fn, sep=',')
#    print('Time-series of the shoreline change along the transects saved as:\n%s'%fn)
    
    return cross_dist












def Timestack(base_folder,beforeDate,eventDate,trackthroughDate,eventName,sitename,normalize,normalizing_transect,fat,InterestedPolygon):

    TC = True
    
    
    
    filepath = os.path.join(base_folder, 'data', sitename)
    
        
        
    with open(filepath+"\\metadata.pkl", 'rb') as f:
        metadata = pickle.load(f) 
    with open(filepath+"\\output.pkl", 'rb') as f:
        output = pickle.load(f) 
    with open(filepath+"\\cross_distance.pkl", 'rb') as f:
        cross_distance = pickle.load(f) 
#    with open(filepath+"\\cross_distance_tidally_corrected.pkl", 'rb') as f:
#        cross_distance_tidally_corrected = pickle.load(f) 
    
    
    
    
    if not fat:
        with open(filepath+"\\transects.pkl", 'rb') as f:
            transects = pickle.load(f) 
            
            
            
            
            
    if fat:
        sf = shapefile.Reader(filepath+"\\transects.shp")
        shaperecs = sf.shapeRecords()
        transects = {}
        k = -1
        for shaperec in shaperecs:
            k=k+1
            transects['Transect '+str(k)]=np.array([list(shaperec.shape.points[0]),list( shaperec.shape.points[1])])
        sf.close()
        
    
    
    
    
    

    
    
    
    
    
    #cross_distance_tidally_corrected    
    tz = output['dates'][0].tzinfo
    output2={}
    shores=[]
    dates = []
    
    
    try:
        sf = shapefile.Reader(filepath+"\\"+sitename+"_outputTidalC.shp")
        shapeRecs = sf.shapeRecords()
        for shapeRec in shapeRecs:
            shores.append(np.array(shapeRec.shape.points))
            dates.append(datetime.datetime(year=shapeRec.record[0].year,month=shapeRec.record[0].month, day=shapeRec.record[0].day, tzinfo=tz))
    except:
        print('TIDALLY CORRECTED LOAD FAILED - FALLING BACK TO TRIMMED SHORELINES')
        sf = shapefile.Reader(filepath+"\\"+sitename+"_outputTrim.shp")
        shapeRecs = sf.shapeRecords()
        for shapeRec in shapeRecs:
            shores.append(np.array(shapeRec.shape.points))
            dates.append(datetime.datetime(year=shapeRec.record[0].year,month=shapeRec.record[0].month, day=shapeRec.record[0].day, tzinfo=tz))
        
        
        
    output2['shorelines']=shores#################################################################################
    output2['dates']=dates#######################################################################################
    cross_distance_tidally_corrected = compute_intersection(output2, transects) ###################################################
    
    
    
    
    
    
    
    utc = output['dates'][0].tzinfo
    eventDatetime = datetime.datetime(year= int(eventDate[0:4]),month=int(eventDate[5:7]) ,day= int(eventDate[-2:]),tzinfo=utc)
    trackthroughDatetime = datetime.datetime(year= int(trackthroughDate[0:4]),month=int(trackthroughDate[5:7]) ,day= int(trackthroughDate[-2:]),tzinfo=utc)
    
    WWextentLine = None
    
    
    
    
    
    
    # BEFORE BEFORE BEFORE BEFORE
    
    eventName2 = eventName+'BEFORE'
    beforeDatetime = datetime.datetime(year= int(beforeDate[0:4]),month=int(beforeDate[5:7]) ,day= int(beforeDate[-2:]),tzinfo=utc)
    before_list = np.where((np.array(output['dates'])<eventDatetime)&(np.array(output['dates'])>beforeDatetime))[0]
    outputNew = {}
    for key in output:
        listNew = []
        n=-1
        for entry in output[key]:
            n=n+1
            if n in before_list:
                listNew.append(entry)
        outputNew[key]=listNew
    
    cross_distance_tidally_correctedNew = {}
    for key in cross_distance_tidally_corrected:
        listNew = []
        n=-1
        for entry in list(cross_distance_tidally_corrected[key]):
            n=n+1
            if n in before_list:
                listNew.append(entry)
        cross_distance_tidally_correctedNew[key]=np.array(listNew)
        
    cross_distanceNew = {}
    for key in cross_distance:
        listNew = []
        n=-1
        for entry in list(cross_distance[key]):
            n=n+1
            if n in before_list:
                listNew.append(entry)
        cross_distanceNew[key]=np.array(listNew)
        
        
    
    if TC:
        avg_positions = []
        for key in cross_distance_tidally_correctedNew:
            distances = cross_distance_tidally_correctedNew[key]
#            avg_positions.append(np.nanmean(distances[-3:])) # last 3 shorelineslines to create historic mean for comparison
            avg_positions.append(np.nanmean(distances)) # all before shorelines to create historic mean for comparison
    else:
        avg_positions = []
        for key in cross_distanceNew:
            distances = cross_distanceNew[key]
#            avg_positions.append(np.nanmean(distances[-3:]))
            avg_positions.append(np.nanmean(distances))
    
    
    # AFTER AFTER AFTER AFTER
    
    
    eventName2 = eventName+'AFTER'
    
#    after_list = np.where(np.array(output['dates'])>eventDatetime)[0]
    after_before_list = np.where((np.array(output['dates'])<trackthroughDatetime)&(np.array(output['dates'])>beforeDatetime))[0]
    outputNew = {}
    for key in output:
        listNew = []
        n=-1
        for entry in output[key]:
            n=n+1
            if n in after_before_list:
                listNew.append(entry)
        outputNew[key]=listNew
    
    cross_distance_tidally_correctedNew = {}
    for key in cross_distance_tidally_corrected:
        listNew = []
        n=-1
        for entry in list(cross_distance_tidally_corrected[key]):
            n=n+1
            if n in after_before_list:
                listNew.append(entry)
        cross_distance_tidally_correctedNew[key]=np.array(listNew)
        
    cross_distanceNew = {}
    for key in cross_distance:
        listNew = []
        n=-1
        for entry in list(cross_distance[key]):
            n=n+1
            if n in after_before_list:
                listNew.append(entry)
        cross_distanceNew[key]=np.array(listNew)
        
    
    
    slopes = []
    slopedates = []
    end = np.where(np.array(outputNew['dates'])==nearest(outputNew['dates'], trackthroughDatetime))[0][0]
    for n in range(0,end):
        if TC:
            new_positions = []
            for key in cross_distance_tidally_correctedNew:
                distances = cross_distance_tidally_correctedNew[key]
                new_positions.append(distances[n])
        else:
            new_positions = []
            for key in cross_distanceNew:
                distances = cross_distanceNew[key]
                new_positions.append(distances[n])
        
        
        
        slopes.append(np.array(new_positions)-np.array(avg_positions))
        slopedates.append(outputNew['dates'][n])
    
    
    slopedates2=[]
    for slopedate in slopedates:
        slopedates2.append(str(slopedate)[0:10])
        
    if normalize:
        slopes1=[]
        for slope in slopes:
            normal=slope[normalizing_transect]
            slop1=[]
            for slop in slope:
                slop1.append(slop-normal)
            slopes1.append(slop1)
            slopes=copy.copy(slopes1)
            
        
        
        
    
    vmin = -1*np.nanmax(np.abs(slopes))
    vmax = np.nanmax(np.abs(slopes))
    nnn=-1
    
    
    
    



    # plot timestack 
    fig,ax=plt.subplots(figsize=(25,25))
    slopes2 = np.array(slopes)
    plt.imshow(slopes2, cmap=mpl.cm.seismic_r, vmin=vmin, vmax=vmax)
    plt.xticks([])
    plt.yticks([])
    #    ax.set_yticks(np.arange(0,len(slopedates2)))
    #    ax.set_yticklabels(slopedates2)
    #    plt.tick_params(labelsize=40)
    date = str(outputNew['dates'][0])[0:10]
    plt.text(-12,0,date,size=20)
    date = str(outputNew['dates'][-1])[0:10]
    plt.text(-12,len(slopes2),date,size=20)
    td = (outputNew['dates'][-1] - outputNew['dates'][0]).days
    tdd = (eventDatetime - outputNew['dates'][0]).days
    tddd = tdd/td
    tdddd = tddd*len(slopes2)

    plt.text(-12,tdddd,str(eventDatetime)[0:10],size=20)
    # plt.show()
    plt.savefig(filepath+"\\"+eventName+"_timestack.png")
    plt.close('all')
        
        

    # plot colorbar
    cmapp = mpl.cm.seismic_r
    normm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.80, 0.9, 0.1])
    cb = mpl.colorbar.ColorbarBase(ax, orientation='horizontal', cmap=cmapp,norm=normm,label='Shoreline Change (m)')
    plt.tick_params(labelsize=40)
    # plt.show()
    plt.savefig(filepath+"\\"+eventName+"_timestackCbar.png")
    plt.close('all')



    # plot basemap
    # last_s2 = np.where(np.array(output['satname'])=='S2')[0][-1]
    fig = plt.figure(figsize=(15,15))
    # plt.imshow(output['imRGBs'][last_s2], extent = output['Extents'][last_s2])
    j=-1
    for key in transects.keys():
        j=j+1
        if j==0:
            plt.scatter(transects[key][0][0], transects[key][0][1],c='red', s= 100)
        else:
            plt.scatter(transects[key][0][0], transects[key][0][1],c='black', s= 100)
            
            
            
    if InterestedPolygon != 'na':
        sf = shapefile.Reader(InterestedPolygon)
        shaperecs = sf.shapeRecords()
        for shaperec in shaperecs:
            for point in shaperec.shape.points:
                plt.scatter(point[0],point[1],c='magenta', s= 100)
        sf.close()
            
    plt.savefig(filepath+"\\"+eventName+"_timestackMap.png")
    plt.close('all')



    print('Done!')
        












 
USERNAME = os.getlogin()
ccSettings = pickle.load(open("C:\\Users\\"+USERNAME+"\\Desktop\\coastsatv2\\temp\\ccpickle.pkl", 'rb'))


base_folder = ccSettings['base_folder']
sitename = ccSettings['sitename']
beforeDate = ccSettings['beforeDate']
eventDate = ccSettings['eventDate']
trackthroughDate = ccSettings['trackthroughDate']
normalize = ccSettings['normalize']
normalizing_transect = ccSettings['normalizing_transect']
fat = ccSettings['fat']
InterestedPolygon = ccSettings['InterestedPolygon']
eventName = ccSettings['eventName']
ccSettings = None

Timestack(base_folder, beforeDate,eventDate,trackthroughDate,eventName,sitename,normalize,normalizing_transect,fat,InterestedPolygon)














