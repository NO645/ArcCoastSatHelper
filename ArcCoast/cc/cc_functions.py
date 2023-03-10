import time
import numpy as np
import pickle
from netCDF4 import Dataset
import xarray as xr
import warnings
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from coastsat import SDS_download, SDS_preprocess, SDS_shoreline, SDS_tools, SDS_transects, SDS_slope
from osgeo import gdal
from osgeo import osr
import os, sys, glob, shutil, subprocess, copy
from osgeo import ogr
import shapefile
import scipy.stats as st
import datetime
from scipy import spatial
from scipy import interpolate
from scipy import signal
from scipy import optimize
import statsmodels as sm
import statsmodels.api as sm
from statsmodels.graphics.api import qqplot
import netCDF4 as nc
import matplotlib as mpl
import matplotlib.cm as cm
from scipy import stats
import pylab as pl
from scipy import ndimage
import dask.array as da
from scipy.spatial.transform import Rotation as R
from osgeo.gdalconst import GA_Update
import urllib.request
import pytz
import urllib.request
from shapely.geometry import LineString
import geopandas
import matplotlib.image as mpimg
from scipy.spatial.distance import directed_hausdorff
from sklearn import neighbors, datasets
from matplotlib import gridspec



import cc_save






def UserTransectt(UserTransect):
    transects = {}
    sf = shapefile.Reader(UserTransect)
    shaperecs = sf.shapeRecords()
    n=0
    for shaperec in shaperecs:
        n=n+1
        transects['Transect ' + str(n)] = np.array(shaperec.shape.points)
    sf.close()
    


def SortImagesBig(filepath, imageSortPkl):
    site = filepath+"\\jpg_files\detection"
    with open(imageSortPkl, 'rb') as f:
        clf = pickle.load(f) 
    
    # delete folders
    a = glob.glob(site+"\\*.png")
    b = glob.glob(site+"\\*")
    aa = set(a)
    bb=set(b)
    cc = bb-aa
    for c in cc:
        shutil.rmtree(c)
        
        
    # make bad folder
    try:
        os.mkdir(site+"\\bad")
    except:
        pass
    

    
    sitename = site.split('\\')[-3]
    shapefilename = filepath+'\\'+sitename+'_output.shp'
    shapefilee = shapefile.Reader(shapefilename)
    shapeRecs = shapefilee.shapeRecords()
    shapedates = []
    for recordd in shapeRecs:
        y=str(recordd.record[0].year)
        m=str(recordd.record[0].month)
        d=str(recordd.record[0].day)
        if len(m)<2:
            m = '0'+m
        if len(d)<2:
            d = '0'+d
        shapedates.append(y+m+d+recordd.record[1])
        
    shapedates = np.array(shapedates)
    
    
    
    # go through images
    l5length=[]
    l5center=[]
    l5lenBreaks=[]
    l5lines=[]
    checkin=0
    for aaa in a:
        checkin=checkin+1
        
        shapedate = aaa.split('\\')[-1].split('-')[0]+aaa.split('\\')[-1].split('-')[1]+aaa.split('\\')[-1].split('-')[2]+aaa.split('\\')[-1].split('-')[-1].split('_')[-1].split('.')[0]
        
        
        
        try:
            h = np.where(shapedates ==shapedate)[0][0]
        except:
            shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])
            continue
            
        
        recordd = shapeRecs[h]
        points = recordd.shape.points
        
        
        
        
        
        
        pointDistances = []
        for n in range(0,len(points)-1):
            dd = distancePP(points[n],points[n+1])
            pointDistances.append(dd)
            
        
        breakpoints = [0] + list(np.where(np.array(pointDistances)>100)[0])+[len(points)]
        lengths = []
        for n in np.flip(np.arange(1,len(breakpoints))):
            lengths.append(breakpoints[n]-breakpoints[n-1])
        lengths = np.flip(np.array(lengths))
        g = np.where(lengths==np.nanmax(lengths))[0][0]
        newShore = points[breakpoints[g]:breakpoints[g+1]]
        
        points=copy.copy(newShore)
        breakpoints=len(breakpoints)
        
        
        
        length = LineString(points).length
        centroid = (LineString(points).centroid.x,LineString(points).centroid.y)
        
        
        
        if len(l5length)<2: #want at least x images
            key = GiveUserChoice2(aaa, site)
            if key ==1:
                l5length.append(length)
                l5center.append(centroid)
                l5lenBreaks.append(breakpoints)
                l5lines.append(points)
            if key ==0:
                shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])
            
            
            
            
            
            
            
        elif checkin>=100: # check every once in a while
            key = GiveUserChoice2(aaa, site)
            if key ==1:
                l5length.append(length)
                l5center.append(centroid)
                l5lenBreaks.append(breakpoints)
                l5lines.append(points)
                checkin = 0
            if key ==0:
                shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])




            
        else: #if its NOT the first few or the checkin....
            testResults = TEST(clf, length, l5length, centroid,l5center,breakpoints,l5lenBreaks,points,l5lines)
            if testResults==2: # if it passes test
                pass
            
            elif testResults==0: # if it fails test
                shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])
                
            else: # if test is unsure
                key = GiveUserChoice2(aaa, site)
                if key ==1:
                    l5length.append(length)
                    l5center.append(centroid)
                    l5lenBreaks.append(breakpoints)
                    l5lines.append(points)
                if key ==0:
                    shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])




        if len(l5length)>5: #only 5 most recent images
            l5length=l5length[-5:]
            l5center=l5center[-5:]
            l5lenBreaks=l5lenBreaks[-5:]
            l5lines=l5lines[-5:]

    shapefilee.close
    
    

















def SortImagesSmall(filepath, imageSortPkl):
    site = filepath+"\\jpg_files\detection"
    with open(imageSortPkl, 'rb') as f:
        clf = pickle.load(f) 
    
    # reset from prior sorts
    b = glob.glob(site+"\\bad\\*.jpg")
    for bb in b:
        shutil.move(bb)
        
        
    # make bad folder
    try:
        os.mkdir(site+"\\bad")
    except:
        pass
    
    a = glob.glob(site+"\\*.jpg")

    
    sitename = site.split('\\')[-3]
    shapefilename = filepath+'\\'+sitename+'_output.shp'
    shapefilee = shapefile.Reader(shapefilename)
    shapeRecs = shapefilee.shapeRecords()
    shapedates = []
    for recordd in shapeRecs:
        y=str(recordd.record[0].year)
        m=str(recordd.record[0].month)
        d=str(recordd.record[0].day)
        if len(m)<2:
            m = '0'+m
        if len(d)<2:
            d = '0'+d
        shapedates.append(y+m+d+recordd.record[1])
        
    shapedates = np.array(shapedates)
    
    
    
    # go through images
    l5length=[]
    l5center=[]
    l5lenBreaks=[]
    l5lines=[]
    checkin=0
    for aaa in a:
        checkin=checkin+1
        
        shapedate = aaa.split('\\')[-1].split('-')[0]+aaa.split('\\')[-1].split('-')[1]+aaa.split('\\')[-1].split('-')[2]+aaa.split('\\')[-1].split('-')[-1].split('_')[-1].split('.')[0]
        
        
        
        try:
            h = np.where(shapedates ==shapedate)[0][0]
        except:
            shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])
            continue
            
        
        recordd = shapeRecs[h]
        points = recordd.shape.points
        
        
        
        
        
        
        pointDistances = []
        for n in range(0,len(points)-1):
            dd = distancePP(points[n],points[n+1])
            pointDistances.append(dd)
            
        
        breakpoints = [0] + list(np.where(np.array(pointDistances)>100)[0])+[len(points)]
        lengths = []
        for n in np.flip(np.arange(1,len(breakpoints))):
            lengths.append(breakpoints[n]-breakpoints[n-1])
        lengths = np.flip(np.array(lengths))
        g = np.where(lengths==np.nanmax(lengths))[0][0]
        newShore = points[breakpoints[g]:breakpoints[g+1]]
        
        points=copy.copy(newShore)
        breakpoints=len(breakpoints)
        
        
        
        length = LineString(points).length
        centroid = (LineString(points).centroid.x,LineString(points).centroid.y)
        
        
        
        if len(l5length)<5: #want at least x images
            key = GiveUserChoice2(aaa, site)
            if key ==1:
                l5length.append(length)
                l5center.append(centroid)
                l5lenBreaks.append(breakpoints)
                l5lines.append(points)
            if key ==0:
                shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])
            
            
            
            
            
            
            
        elif checkin>=50: # check every once in a while
            key = GiveUserChoice2(aaa, site)
            if key ==1:
                l5length.append(length)
                l5center.append(centroid)
                l5lenBreaks.append(breakpoints)
                l5lines.append(points)
                checkin = 0
            if key ==0:
                shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])




            
        else: #if its NOT the first few or the checkin....
            testResults = TEST(clf, length, l5length, centroid,l5center,breakpoints,l5lenBreaks,points,l5lines)
            if testResults==2: # if it passes test
                pass
            
            elif testResults==0: # if it fails test
                shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])
                
            else: # if test is unsure
                key = GiveUserChoice2(aaa, site)
                if key ==1:
                    l5length.append(length)
                    l5center.append(centroid)
                    l5lenBreaks.append(breakpoints)
                    l5lines.append(points)
                if key ==0:
                    shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])




        if len(l5length)>5: #only 5 most recent images
            l5length=l5length[-5:]
            l5center=l5center[-5:]
            l5lenBreaks=l5lenBreaks[-5:]
            l5lines=l5lines[-5:]

    shapefilee.close
    
    



def lengthTest(length, l5length):
    l5length = np.array(l5length)
    lengthpercents = (length-l5length)/l5length
    return lengthpercents
    
    
    
    
def centroidTest(centroid,l5center):
    l5center = np.array(l5center)
    distances=(np.sqrt((centroid[0]-l5center[:,0])**2+(centroid[1]-l5center[:,1])**2))
    return distances
    
def hausdorffTest(points, l5lines):
    u = points
    distances = []
    for v in l5lines:
        distances.append(np.nanmax([directed_hausdorff(u, v)[0], directed_hausdorff(v, u)[0]]))
    return np.array(distances)
    
def breaksTest(breakpoints,lenBreaks):
    #test for good
    breakmean= np.nanmean(lenBreaks)
    breaksd = np.std(lenBreaks) # if true its good
    breaks = (breakmean*np.ones(len(lenBreaks))+2*breaksd>breakpoints) & (breakmean*np.ones(len(lenBreaks))-2*breaksd<breakpoints) | (np.nanmax(lenBreaks)*np.ones(len(lenBreaks))==breakpoints) | (np.nanmin(lenBreaks)*np.ones(len(lenBreaks))==breakpoints)
    breaks=np.array(breaks)
    
    return breaks.all()




def TEST(clf, length, l5length, centroid,l5center,breakpoints,l5lenBreaks,points,l5lines):
    lengthResult = lengthTest(length, l5length)
    centroidResult = centroidTest(centroid,l5center)
    breakResult = breaksTest(breakpoints,l5lenBreaks)
    hausdorffResult = hausdorffTest(points,l5lines)
    
#    print(lengthResult)
#    print(centroidResult)
#    print(breakResult)
#    print(hausdorffResult)
    
    percentgood = 0.019
    percentbad = 0.25
    
    cdistancegood = 75
    cdistancebad = 400
    
    hdistancegood = 50
    hdistancebad = 500
    
    lengthgeneralgood = (lengthResult<percentgood).any()
    centroidgeneralgood = (centroidResult<cdistancegood).any()
    hausdorffgeneralgood =(hausdorffResult<hdistancegood).any()
    
    lengthgeneralbad = (lengthResult>percentbad).all()
    centroidgeneralbad = (centroidResult>cdistancebad).all()
    hausdorffgeneralbad = (hausdorffResult>hdistancebad).all()
    
    
#    ZZ = []
#    for n in range(0, len(lengthResult)):
#        zz = clf.predict(np.array([lengthResult[n], hausdorffResult[n]]))
#        ZZ.append(zz)
        
    zz = clf.predict(np.dstack((lengthResult, hausdorffResult))[0])
    
    
    if np.nansum(zz) == len(zz):
        result = 2
    elif np.nansum(zz) == 0:
        result = 0
    else:
        result = 1
    
    return result
    
    

def GiveUserChoice2(aaa, site):    
    img = mpimg.imread(aaa)




    fig = plt.figure()
    fig.set_size_inches([18, 9])
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()



    # create image 1 (RGB)
    plt.imshow(img)
    plt.title(site, fontweight='bold', fontsize=16)



    # set a key event to accept/reject the detections (see https://stackoverflow.com/a/15033071)
    # this variable needs to be immuatable so we can access it after the keypress event
    key_event = {}
    def press(event):
        # store what key was pressed in the dictionary
        key_event['pressed'] = event.key
    # let the user press a key, right arrow to keep the image, left arrow to skip it
    # to break the loop the user can press 'escape'
    while True:
        btn_keep = plt.text(1.1, 0.9, 'keep ⇨', size=12, ha="right", va="top",   bbox=dict(boxstyle="square", ec='k',fc='w'))
        btn_skip = plt.text(-0.1, 0.9, '⇦ trash', size=12, ha="left", va="top",   bbox=dict(boxstyle="square", ec='k',fc='w'))
        btn_esc = plt.text(0.5, 0, '<esc> to quit', size=12, ha="center", va="top",  bbox=dict(boxstyle="square", ec='k',fc='w'))
        plt.draw()
        fig.canvas.mpl_connect('key_press_event', press)
        plt.waitforbuttonpress()
        # after button is pressed, remove the buttons
        btn_skip.remove()
        btn_keep.remove()
        btn_esc.remove()
        
 # keep/skip image according to the pressed key, 'escape' to break the loop
        if key_event.get('pressed') == 'left':
            skip_image = 0
            break
        elif key_event.get('pressed') == 'right':
            skip_image = 1
            break
        else:
            plt.waitforbuttonpress()


    plt.close('all')
    
    
    return skip_image



    
    
    
    
    
    
    
    






def TrimOutputShorelines2(output,filepath,sitename):
    shoresShapefile = filepath+"\\"+sitename+"_output.shp"
    shoresShapefileC = filepath+"\\"+sitename+"_outputTrim.shp"
    
    pointDistances = []
    
    
    shoreshape = shapefile.Reader(shoresShapefile)
    shoreRecs = shoreshape.shapeRecords()
    
    w = shapefile.Writer(shoresShapefileC)
    w.fields = shoreshape.fields[1:]    
    
    
    newshores = []
    for shoreRec in shoreRecs:
        shorepoints = shoreRec.shape.points
        pointDistances = []
        for n in range(0,len(shorepoints)-1):
            d = distancePP(shorepoints[n],shorepoints[n+1])
            pointDistances.append(d)
            
        
        # plt.scatter(np.array(shorepoints)[:,0], np.array(shorepoints)[:,1])
        # plt.show()
        # plt.close('all')
        
        # plt.scatter(np.arange(0,len(pointDistances)), pointDistances)
        # plt.show()
        # plt.close('all')
        
        breakpoints = [0] + list(np.where(np.array(pointDistances)>100)[0])+[len(shorepoints)]
        lengths = []
        for n in np.flip(np.arange(1,len(breakpoints))):
            lengths.append(breakpoints[n]-breakpoints[n-1])
        lengths = np.flip(np.array(lengths))
        a = np.where(lengths==np.nanmax(lengths))[0][0]
        newShore = shorepoints[breakpoints[a]:breakpoints[a+1]]
#        print(newShore)
        w.record(shoreRec.record[0],shoreRec.record[1],shoreRec.record[2],shoreRec.record[3])
        w.line([newShore])
        
        newshores.append(newShore)
        
        

    newshores2 = []
    for newshore in newshores:
        newshore2 = []
        for points in newshore:
            newshore2.append(np.array(list(points)))
        newshore2 = np.array(newshore2)
        newshores2.append(newshore2)



        
    output['TrimmedShorelines']=newshores2
        
    w.close()
    shutil.copy(filepath+"\\"+sitename+"_output.prj", filepath+"\\"+sitename+"_outputTrim.prj")
    
    return output





def AdjustOutputShapefileByTides4(slope_ests,output,filepath,sitename,Contour):
    shoresShapefile = filepath+"\\"+sitename+"_outputTrim.shp"
    shoresShapefileC = filepath+"\\"+sitename+"_outputTidalC.shp"
    
    
#    with open(filepath+"\\output.pkl", 'rb') as f:
#        output = pickle.load(f) 
#    with open(filepath+"\\slope_ests.pkl", 'rb') as f:
#        slope_ests = pickle.load(f) 


    shoreshape = shapefile.Reader(shoresShapefile)
    shoreRecs = shoreshape.shapeRecords()
    
    
    
    # find beach slope
    slopes = []
    for slope_est in slope_ests:
        slopes.append(slope_ests[slope_est])
    slope = np.nanmean(np.array(slopes))
    
    
    w = shapefile.Writer(shoresShapefileC)
    w.fields = shoreshape.fields[1:]    
    
    lrx, uly, lry, ulx, array_gray = NewTransectsBuilderUtils(sitename, filepath, output)
    
    newshores=[]
    for shoreRec in shoreRecs:
        shorepoints = shoreRec.shape.points
        
        NewTransectsBuilder2(sitename, filepath, output, shorepoints, lrx, uly, lry, ulx, array_gray)
        transectsShapefile = filepath+"\\transectsTideMove.shp"
        transectshape = shapefile.Reader(transectsShapefile)
        transectRecs = transectshape.shapeRecords()
        
        # find tide position
        date = shoreRec.record[0]
        days = []
        months = []
        years = []
        for fdgfd in output['dates']:
            days.append(fdgfd.day==date.day)
            months.append(fdgfd.month==date.month)
            years.append(fdgfd.year==date.year)
            
        for n in range(0,len(days)):
            if days[n] and months[n] and years[n]:
                break
        
        newShore = []
        shorenum=0
        for transectRec in transectRecs:
            shorenum=shorenum+1
            transectpoints = transectRec.shape.points
            a = shorepoints[shorenum]
                
            try:
            
                v = np.array([transectpoints[1][0],transectpoints[1][1]])-np.array([a[0],a[1]])
                u = v/(v[0]**2+v[1]**2)**0.5
                
                distancee1 = (output['tides'][n])/slope
                distancee2 = (-Contour)/slope
                distancee=distancee1+distancee2
                newPoint = np.array([a[0],a[1]])+distancee*u
                
#                plt.scatter(newPoint[0],newPoint[1],c='red')
#                plt.scatter(a[0],a[1],c='black')
#                plt.plot(np.array(transectpoints)[:,0],np.array(transectpoints)[:,1],c='magenta')
#                plt.plot(np.array(shorepoints)[:,0],np.array(shorepoints)[:,1],c='blue')
#                plt.axis([a[0]-50,a[0]+50,a[1]-50,a[1]+50])
#                plt.show()
#                plt.close('all')
                
                newShore.append(list(newPoint))
            except:
                pass
        
        w.record(shoreRec.record[0],shoreRec.record[1],shoreRec.record[2],shoreRec.record[3])
        w.line([newShore])
        
        newshores.append(newShore)
    output['TrimmedCorrectedShorelines']=newshores
        
        
        
        
#         import matplotlib.pyplot as plt
#         plt.figure(figsize=(25,25))
#         plt.scatter(np.array(newShore)[:,0],np.array(newShore)[:,1],c='red')
#         plt.plot(np.array(shorepoints)[:,0],np.array(shorepoints)[:,1],c='black')
#         plt.show()
#         plt.close('all')
        
        
    w.close()
    shutil.copy(filepath+"\\"+sitename+"_output.prj", filepath+"\\"+sitename+"_outputTidalC.prj")
    shoreshape.close()
    transectshape.close()
    
    return output



def distancePointPoint(point1, point2):
    distance = (np.sqrt((point1[0]-point2[0])**2+(point1[1]-point2[1])**2))
    return distance



def AdjustOutputShapefileByTides5(slope_ests,output,filepath,sitename,Contour,transects):
    
    # for the chosen point along the reference shoreline find a slope for best fit line
    shoreshape = shapefile.Reader(filepath+"\\"+sitename+"_reference_shoreline.shp")
    shoreRece = shoreshape.shapeRecords()[0]
    shorepoints=shoreRece.shape.points
    shoreshape.close()
    
    deriv_points = np.array(shorepoints)
    x = deriv_points[:,0]
    y = deriv_points[:,1]
    slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
    
    # find a slope perpindicular to shorelne slope for transect
    pslope = -1*(1/slope)
    angle = np.arctan(pslope)
    
    # find beach slope
    slopes = []
    for slope_est in slope_ests:
        slopes.append(slope_ests[slope_est])
    slope = np.nanmean(np.array(slopes))
    
    # Get point from ocean
    shoreRece = list(transects.keys())[int(len(transects)/2)]
    ocean_lat_long=transects[shoreRece][0]
    
    # open trimmed shorelines and start tidally corrected shorelines file
    shoreshape = shapefile.Reader(filepath+"\\"+sitename+"_outputTrim.shp")
    shoreRecs = shoreshape.shapeRecords()
    w = shapefile.Writer(filepath+"\\"+sitename+"_outputTidalC.shp")
    w.fields = shoreshape.fields[1:]    
    
    newshores=[]
    for shoreRec in shoreRecs:
        shorepoints = shoreRec.shape.points
        
        # find tide position
        date = shoreRec.record[0]
        days = []
        months = []
        years = []
        for fdgfd in output['dates']:
            days.append(fdgfd.day==date.day)
            months.append(fdgfd.month==date.month)
            years.append(fdgfd.year==date.year)
            
        for n in range(0,len(days)):
            if days[n] and months[n] and years[n]:
                break

        Vdistance = -output['tides'][n] + Contour
        Hshift = Vdistance/slope

        if Hshift>0:
            outTOsea=False
        if Hshift<0:
            outTOsea=True
            
        shorepoints1 = np.array(shorepoints)
        shorepoints2 = np.dstack((shorepoints1[:,0] + np.cos(angle) * Hshift, shorepoints1[:,1] + np.sin(angle) * Hshift))[0]
        shorepoints3 = np.dstack((shorepoints1[:,0] + np.cos(angle) * -Hshift, shorepoints1[:,1] + np.sin(angle) * -Hshift))[0]

        testpoint1 = shorepoints2[int(len(shorepoints2)/2)]
        testpoint2 = shorepoints[int(len(shorepoints)/2)]
        if (distancePointPoint(testpoint1,ocean_lat_long) < distancePointPoint(testpoint2,ocean_lat_long)) and outTOsea:
            NEWshorepoints = shorepoints2
        elif (distancePointPoint(testpoint1,ocean_lat_long) < distancePointPoint(testpoint2,ocean_lat_long)) and not outTOsea:
            NEWshorepoints = shorepoints3
        elif (distancePointPoint(testpoint1,ocean_lat_long) > distancePointPoint(testpoint2,ocean_lat_long)) and outTOsea:
            NEWshorepoints = shorepoints3
        elif (distancePointPoint(testpoint1,ocean_lat_long) > distancePointPoint(testpoint2,ocean_lat_long)) and not outTOsea:
            NEWshorepoints = shorepoints2
        else:
            NEWshorepoints = shorepoints1



        w.record(shoreRec.record[0],shoreRec.record[1],shoreRec.record[2],shoreRec.record[3])
        w.line([NEWshorepoints])
        
        newshores.append(NEWshorepoints)
    output['TrimmedCorrectedShorelines']=newshores
        
        
    w.close()
    shutil.copy(filepath+"\\"+sitename+"_output.prj", filepath+"\\"+sitename+"_outputTidalC.prj")
    shoreshape.close()
    
    return output








def NewTransectsBuilderUtils(sitename, filepath, output):
    
    array_gray = copy.copy(output['imClassifs'][0])
    
        
    top = array_gray[1,:]
    right = array_gray[:,-2]
    bottom = array_gray[-2,:]
    left = array_gray[:,1]
    
    top2 = len(np.where(top==0.4284314)[0])/len(top)
    left2 = len(np.where(left==0.4284314)[0])/len(left)
    bottom2 = len(np.where(bottom==0.4284314)[0])/len(bottom)
    right2 = len(np.where(right==0.4284314)[0])/len(right)
    
    test = [top2, left2, bottom2, right2]
    
    open_water_position = np.where(test == np.nanmax(test))[0][0]
    
    

    
    # open a rster to get extent for transects
    try:
        rasters = glob.glob(filepath+"\\S2\\10m\\*_"+sitename+"_10m.tif")
        src = gdal.Open(rasters[0])
    except:
        rasters = glob.glob(filepath+"\\L8\\pan\\*.tif")
        src = gdal.Open(rasters[0])
    
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    src = None
    
    tryagain = 1
    while lrx-ulx < 1000:
        try:
            src = gdal.Open(rasters[tryagain])
            ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
            lrx = ulx + (src.RasterXSize * xres)
            lry = uly + (src.RasterYSize * yres)
            src = None
            tryagain = tryagain+1
        except:
            break
    return lrx, uly, lry, ulx, array_gray





    
    
    
def TrimOutputShorelines(filepath,sitename):
    shoresShapefile = filepath+"\\"+sitename+"_output.shp"
    shoresShapefileC = filepath+"\\"+sitename+"_outputTrim.shp"
    
    pointDistances = []
    
    
    shoreshape = shapefile.Reader(shoresShapefile)
    shoreRecs = shoreshape.shapeRecords()
    
    w = shapefile.Writer(shoresShapefileC)
    w.fields = shoreshape.fields[1:]    
    
    for shoreRec in shoreRecs:
        shorepoints = shoreRec.shape.points
        pointDistances = []
        for n in range(0,len(shorepoints)-1):
            d = distancePP(shorepoints[n],shorepoints[n+1])
            pointDistances.append(d)
            
        
        # plt.scatter(np.array(shorepoints)[:,0], np.array(shorepoints)[:,1])
        # plt.show()
        # plt.close('all')
        
        # plt.scatter(np.arange(0,len(pointDistances)), pointDistances)
        # plt.show()
        # plt.close('all')
        
        breakpoints = [0] + list(np.where(np.array(pointDistances)>100)[0])+[len(shorepoints)]
        lengths = []
        for n in np.flip(np.arange(1,len(breakpoints))):
            lengths.append(breakpoints[n]-breakpoints[n-1])
        lengths = np.flip(np.array(lengths))
        a = np.where(lengths==np.nanmax(lengths))[0][0]
        newShore = shorepoints[breakpoints[a]:breakpoints[a+1]]
#        print(newShore)
        w.record(shoreRec.record[0],shoreRec.record[1],shoreRec.record[2],shoreRec.record[3])
        w.line([newShore])
        
    w.close()
    shutil.copy(filepath+"\\"+sitename+"_output.prj", filepath+"\\"+sitename+"_outputTrim.prj")
    
    



def NewTransectsBuilder2(sitename, filepath, output, shorepoints, lrx, uly, lry, ulx, array_gray):
    

    # add prpidicular points to new transect shapefile
    
    w = shapefile.Writer(filepath+"\\transectsTideMove.shp")
    w.field('name', 'C')
    transects = dict([])
    transectsFull = dict([])
    transectsRotations = []
    transectsSlopes = []
    COUNTERtransectsfull=0
    COUNTERtransectssmall=0
    COUNTERtransectssmall2=0
    for n in range(1,len(shorepoints)-2):######
        try:
            COUNTERtransectsfull=COUNTERtransectsfull+1
            COUNTERtransectssmall=COUNTERtransectssmall+1
            
            # for the chosen point along the reference shoreline find a slope for best fit line
            deriv_points = np.array(shorepoints[n-1:n+1])######
            x = deriv_points[:,0]
            y = deriv_points[:,1]
            X = np.array(shorepoints)[n,0]
            Y = np.array(shorepoints)[n,1]
            slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
            
            # find a slope perpindicular to shorelne slope for transect
            pslope = -1*(1/slope)
            
            # start transect with first value
            XX = [X]
            YY = [Y + pslope*(XX[-1]-X)]
            
            # grow transect in forward X direction
            waterPointsForward = 0
            PointsForward = 0
            i = 0
            # greyClassification2 = [0]######################################################################
            while XX[-1]<lrx+100 and YY[-1]<uly+100 and YY[-1]>lry+100:
                i=i+1
                XX.append(XX[-1]+1)
                YY.append(Y + pslope*(XX[-1]-X))
                greyClassification = getArrayValueScaled(XX[-1], YY[-1], array_gray, lry,uly,ulx,lrx)
                # greyClassification2.append(greyClassification)######################################################################
                if ((greyClassification==2) or (greyClassification==3)):
                    waterPointsForward = waterPointsForward+1
                PointsForward = PointsForward+1
            # grow transect in backward X direction
            waterPointsReverse = 0
            PointsReverse = 0
            i = 0
            while XX[0]>ulx+100 and YY[0]<uly+100 and YY[0]>lry+100:
                i=i+1
                XX.insert(0, XX[0]-1)
                YY.insert(0, Y + pslope*(XX[0]-X))
                greyClassification = getArrayValueScaled(XX[0], YY[0], array_gray, lry,uly,ulx,lrx)
                # greyClassification2.insert(0, greyClassification)######################################################################
                if ((greyClassification==2) or (greyClassification==3)):
                    waterPointsReverse = waterPointsReverse+1
                PointsReverse = PointsReverse+1
                

            if (waterPointsReverse/PointsReverse) > (waterPointsForward/PointsForward):
                # cccc = 'black'
                XX = np.flip(np.array(XX))
                YY = np.flip(np.array(YY))
                

            
            raw_transect_line = []
            nn=-1
            for all in YY:
                nn=nn+1
                raw_transect_line.append((XX[nn],YY[nn]))
            
            
            distances = distance(raw_transect_line, (X,Y)) 
            raw_transect_line = np.array(raw_transect_line)
            transect_start_point = np.where(distances==nearest(distances,100))[0][0]
            x1 = raw_transect_line[transect_start_point][0]
            y1 = raw_transect_line[transect_start_point][1]
            x2 = raw_transect_line[-1][0]
            y2 = raw_transect_line[-1][1]
            accept = True
        except:
            COUNTERtransectsfull=COUNTERtransectsfull-1
            COUNTERtransectssmall=COUNTERtransectssmall-1
            accept = False



        if accept:
            w.line([[ [x1, y1], [x2, y2] ]])
            
            
#            plt.figure(figsize=(25,25))
#            plt.plot([x1,x2],[y1,y2])
##            plt.plot(raw_transect_line[:,0], raw_transect_line[:,1])
#            plt.scatter(X,Y, c='red',s=25)
##            plt.scatter(XX,YY, c='black',s=25,alpha=.25)
#            plt.show()
#            plt.close('all')
#            break
            
            
            
            transectsFull['Transect '+str(COUNTERtransectsfull)] = np.array([ [x1, y1], [x2, y2] ])
            
            if slope < 0:
                angle = -1*(np.arctan(slope) - (np.pi/2))
            else:    
                angle = (np.pi/2) - np.arctan(slope)
            rotation = R.from_euler('z', [angle])
            transectsRotations.append(rotation)
            transectsSlopes.append(pslope)
            
            
            
            if COUNTERtransectssmall == 70:
                COUNTERtransectssmall=0
                COUNTERtransectssmall2=COUNTERtransectssmall2+1
                transects['Transect '+str(COUNTERtransectssmall2)] = np.array([ [x1, y1], [x2, y2] ])
            w.record(str(n))
        

    w.close()
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transectsTideMove.prj")
    
    





def AdjustOutputShapefileByTides3(filepath,sitename):
    shoresShapefile = filepath+"\\"+sitename+"_outputTrim.shp"
    shoresShapefileC = filepath+"\\"+sitename+"_outputTidalC.shp"
    
    
    with open(filepath+"\\output.pkl", 'rb') as f:
        output = pickle.load(f) 
    with open(filepath+"\\slope_ests.pkl", 'rb') as f:
        slope_ests = pickle.load(f) 


    shoreshape = shapefile.Reader(shoresShapefile)
    shoreRecs = shoreshape.shapeRecords()
    
    
    
    # find beach slope
    slopes = []
    for slope_est in slope_ests:
        slopes.append(slope_ests[slope_est])
    slope = np.nanmean(np.array(slopes))
    
    
    w = shapefile.Writer(shoresShapefileC)
    w.fields = shoreshape.fields[1:]    
    
    lrx, uly, lry, ulx, array_gray = NewTransectsBuilderUtils(sitename, filepath, output)
    
    for shoreRec in shoreRecs:
        shorepoints = shoreRec.shape.points
        
        NewTransectsBuilder2(sitename, filepath, output, shorepoints, lrx, uly, lry, ulx, array_gray)
        transectsShapefile = filepath+"\\transectsTideMove.shp"
        transectshape = shapefile.Reader(transectsShapefile)
        transectRecs = transectshape.shapeRecords()
        
        # find tide position
        date = shoreRec.record[0]
        days = []
        months = []
        years = []
        for fdgfd in output['dates']:
            days.append(fdgfd.day==date.day)
            months.append(fdgfd.month==date.month)
            years.append(fdgfd.year==date.year)
            
        for n in range(0,len(days)):
            if days[n] and months[n] and years[n]:
                break
        
        newShore = []
        shorenum=0
        for transectRec in transectRecs:
            shorenum=shorenum+1
            transectpoints = transectRec.shape.points
            a = shorepoints[shorenum]
                
            try:
            
                v = np.array([transectpoints[1][0],transectpoints[1][1]])-np.array([a[0],a[1]])
                u = v/(v[0]**2+v[1]**2)**0.5
                
                distancee = output['tides'][n]/slope
                newPoint = np.array([a[0],a[1]])+distancee*u
                
#                plt.scatter(newPoint[0],newPoint[1],c='red')
#                plt.scatter(a[0],a[1],c='black')
#                plt.plot(np.array(transectpoints)[:,0],np.array(transectpoints)[:,1],c='magenta')
#                plt.plot(np.array(shorepoints)[:,0],np.array(shorepoints)[:,1],c='blue')
#                plt.axis([a[0]-15,a[0]+15,a[1]-15,a[1]+15])
#                plt.show()
#                plt.close('all')
                
                newShore.append(list(newPoint))
            except:
                pass
        
        w.record(shoreRec.record[0],shoreRec.record[1],shoreRec.record[2],shoreRec.record[3])
        w.line([newShore])
        
        
#         import matplotlib.pyplot as plt
#         plt.figure(figsize=(25,25))
#         plt.scatter(np.array(newShore)[:,0],np.array(newShore)[:,1],c='red')
#         plt.plot(np.array(shorepoints)[:,0],np.array(shorepoints)[:,1],c='black')
#         plt.show()
#         plt.close('all')
        
        
    w.close()
    shutil.copy(filepath+"\\"+sitename+"_output.prj", filepath+"\\"+sitename+"_outputTidalC.prj")
    shoreshape.close()
    transectshape.close()







def SlopesEstimates(sitename,filepath,transects,dates,output):
    # 2. Load 2D shorelines and transects
    
    # load the sitename_output.pkl generated by CoastSat
    
    
        
    utc = output['dates'][0].tzinfo
        
        
    
    # # remove S2 shorelines (the slope estimation algorithm needs only Landsat)
    # if 'S2' in output['satname']:
    #     idx_S2 = np.array([_ == 'S2' for _ in output['satname']])
    #     for key in output.keys():
    #         output[key] = [output[key][_] for _ in np.where(~idx_S2)[0]]
    
    bad_list = np.where(np.array(output['satname'])!='L8')[0]
    outputNew = {}
    for key in output:
        listNew = []
        n=-1
        for entry in output[key]:
            n=n+1
            if n not in bad_list:
                listNew.append(entry)
    
        outputNew[key]=listNew
            
            
    
    # remove duplicates 
    output = SDS_slope.remove_duplicates(outputNew)
    # remove shorelines from images with poor georeferencing (RMSE > 10 m)
    output = SDS_slope.remove_inaccurate_georef(output, 10)
    
    # # plot shorelines and transects
    # fig,ax = plt.subplots(1,1,figsize=[12,  8])
    # fig.set_tight_layout(True)
    # ax.axis('equal')
    # ax.set(xlabel='Eastings', ylabel='Northings', title=sitename)
    # ax.grid(linestyle=':', color='0.5')
    # for i in range(len(output['shorelines'])):
    #     coords = output['shorelines'][i]
    #     date = output['dates'][i]
    #     ax.plot(coords[:,0], coords[:,1], '.', label=date.strftime('%d-%m-%Y'))
    # for key in transects.keys():
    #     ax.plot(transects[key][:,0],transects[key][:,1],'k--',lw=2)
    #     ax.text(transects[key][-1,0], transects[key][-1,1], key)
    
    # a more robust method to compute intersection is needed here to avoid outliers
    # as these can affect the slope detection algorithm
    settings_transects = { # parameters for shoreline intersections
                          'along_dist':         25,             # along-shore distance to use for intersection
                          'max_std':            15,             # max std for points around transect
                          'max_range':          30,             # max range for points around transect
                          'min_val':            -100,           # largest negative value along transect (landwards of transect origin)
                          # parameters for outlier removal
                          'nan/max':            'auto',         # mode for removing outliers ('auto', 'nan', 'max')
                          'prc_std':            0.1,            # percentage to use in 'auto' mode to switch from 'nan' to 'max'
                          'max_cross_change':   40,        # two values of max_cross_change distance to use
                          }
    # compute intersections [advanced version]
    cross_distance = SDS_slope.compute_intersection(output, transects, settings_transects) 
    # remove outliers [advanced version]
    cross_distance = SDS_slope.reject_outliers(cross_distance,output,settings_transects)       
    
    
    # # plot time-series
    # SDS_slope.plot_cross_distance(output['dates'],cross_distance)
        
    # slope estimation settings
    days_in_year = 365.2425
    seconds_in_day = 24*3600
    settings_slope = {'slope_min':        0.035,
                      'slope_max':        0.2, 
                      'delta_slope':      0.005,
                      'date_range':       [1999,2020],            # range of dates over which to perform the analysis
                      'n_days':           8,                      # sampling period [days]
                      'n0':               50,                     # for Nyquist criterium
                      'freqs_cutoff':     1./(seconds_in_day*30), # 1 month frequency
                      'delta_f':          100*1e-10,              # deltaf for buffer around max peak                                           # True to save some plots of the spectrums
                      }
    # settings_slope['date_range'] = [pytz.utc.localize(datetime(settings_slope['date_range'][0],5,1)),
    #                                 pytz.utc.localize(datetime(settings_slope['date_range'][1],1,1))]
    date1 = pytz.utc.localize(datetime.datetime(int(dates[0][0:4]),int(dates[0][5:7]),int(dates[0][8:10])))
    date2 = pytz.utc.localize(datetime.datetime(int(dates[1][0:4]),int(dates[1][5:7]),int(dates[1][8:10])))
    settings_slope['date_range'] = [date1,date2]
    
    
    beach_slopes = SDS_slope.range_slopes(settings_slope['slope_min'], settings_slope['slope_max'], settings_slope['delta_slope'])
    
    # clip the dates between 1999 and 2020 as we need at least 2 Landsat satellites 
    idx_dates = [np.logical_and(_>settings_slope['date_range'][0],_<settings_slope['date_range'][1]) for _ in output['dates']]
    dates_sat = [output['dates'][_] for _ in np.where(idx_dates)[0]]
    for key in cross_distance.keys():
        cross_distance[key] = cross_distance[key][idx_dates]
    
    # 3. Tide levels
        
    
    # Option 2. otherwise load tide levels associated with "dates_sat" from a file
    # with open(os.path.join('example_data', sitename + '_tide' + '.pkl'), 'rb') as f:
    #     tide_data = pickle.load(f) 
    # tide_sat = tide_data['tide']
    tide_sat = np.array(output['tides'])
    
    
    
    
    # # plot time-step distribution
    # t = np.array([_.timestamp() for _ in dates_sat]).astype('float64')
    # delta_t = np.diff(t)
    # fig, ax = plt.subplots(1,1,figsize=(12,3), tight_layout=True)
    # ax.grid(which='major', linestyle=':', color='0.5')
    # bins = np.arange(np.min(delta_t)/seconds_in_day, np.max(delta_t)/seconds_in_day+1,1)-0.5
    # ax.hist(delta_t/seconds_in_day, bins=bins, ec='k', width=1);
    # ax.set(xlabel='timestep [days]', ylabel='counts',
    #         xticks=settings_slope['n_days']*np.arange(0,20),
    #         xlim=[0,50], title='Timestep distribution');
    
    # find tidal peak frequency
    settings_slope['freqs_max'] = SDS_slope.find_tide_peak(dates_sat,tide_sat,settings_slope)
    
    # 4. Estimate beach slopes along the transects
    
    slope_est = dict([])
    for key in cross_distance.keys():
        # remove NaNs
        idx_nan = np.isnan(cross_distance[key])
        dates = [dates_sat[_] for _ in np.where(~idx_nan)[0]]
        tide = tide_sat[~idx_nan]
        composite = cross_distance[key][~idx_nan]
        # apply tidal correction
        tsall = SDS_slope.tide_correct(composite,tide,beach_slopes)
        # SDS_slope.plot_spectrum_all(dates,composite,tsall,settings_slope)
        slope_est[key] = SDS_slope.integrate_power_spectrum(dates,tsall,settings_slope)
        # print('Beach slope at transect %s: %.3f'%(key, slope_est[key]))
    
    return slope_est
    
    
    
    
    
    


def FindEPSG(polygon):
    point = polygon[0][1]
    
    a = np.array([[32601,-180.0000,0.0000,-174.0000,84.0000],
    [32602,-174.0000,0.0000,-168.0000,84.0000],
    [32603,-168.0000,0.0000,-162.0000,84.0000],
    [32604,-162.0000,0.0000,-156.0000,84.0000],
    [32605,-156.0000,0.0000,-150.0000,84.0000],
    [32606,-150.0000,0.0000,-144.0000,84.0000],
    [32607,-144.0000,0.0000,-138.0000,84.0000],
    [32608,-138.0000,0.0000,-132.0000,84.0000],
    [32609,-132.0000,0.0000,-126.0000,84.0000],
    [32610,-126.0000,0.0000,-120.0000,84.0000],
    [32611,-120.0000,0.0000,-114.0000,84.0000],
    [32612,-114.0000,0.0000,-108.0000,84.0000],
    [32613,-108.0000,0.0000,-102.0000,84.0000],
    [32614,-102.0000,0.0000,-96.0000,84.0000],
    [32615,-96.0000,0.0000,-90.0000,84.0000],
    [32616,-90.0000,0.0000,-84.0000,84.0000],
    [32617,-84.0000,0.0000,-78.0000,84.0000],
    [32618,-78.0000,0.0000,-72.0000,84.0000],
    [32619,-72.0000,0.0000,-66.0000,84.0000],
    [32620,-66.0000,0.0000,-60.0000,84.0000],
    [32621,-60.0000,0.0000,-54.0000,84.0000],
    [32622,-54.0000,0.0000,-48.0000,84.0000],
    [32623,-48.0000,0.0000,-42.0000,84.0000],
    [32624,-42.0000,0.0000,-36.0000,84.0000],
    [32625,-36.0000,0.0000,-30.0000,84.0000],
    [32626,-30.0000,0.0000,-24.0000,84.0000],
    [32627,-24.0000,0.0000,-18.0000,84.0000],
    [32628,-18.0000,0.0000,-12.0000,84.0000],
    [32629,-12.0000,0.0000,-6.0000,84.0000],
    [32630,-6.0000,0.0000,0.0000,84.0000],
    [32631,0.0000,0.0000,6.0000,84.0000],
    [32632,6.0000,0.0000,12.0000,84.0000],
    [32633,12.0000,0.0000,18.0000,84.0000],
    [32634,18.0000,0.0000,24.0000,84.0000],
    [32635,24.0000,0.0000,30.0000,84.0000],
    [32636,30.0000,0.0000,36.0000,84.0000],
    [32637,36.0000,0.0000,42.0000,84.0000],
    [32638,42.0000,0.0000,48.0000,84.0000],
    [32639,48.0000,0.0000,54.0000,84.0000],
    [32640,54.0000,0.0000,60.0000,84.0000],
    [32641,60.0000,0.0000,66.0000,84.0000],
    [32642,66.0000,0.0000,72.0000,84.0000],
    [32643,72.0000,0.0000,78.0000,84.0000],
    [32644,78.0000,0.0000,84.0000,84.0000],
    [32645,84.0000,0.0000,90.0000,84.0000],
    [32646,90.0000,0.0000,96.0000,84.0000],
    [32647,96.0000,0.0000,102.0000,84.0000],
    [32648,102.0000,0.0000,108.0000,84.0000],
    [32649,108.0000,0.0000,114.0000,84.0000],
    [32650,114.0000,0.0000,120.0000,84.0000],
    [32651,120.0000,0.0000,126.0000,84.0000],
    [32652,126.0000,0.0000,132.0000,84.0000],
    [32653,132.0000,0.0000,138.0000,84.0000],
    [32654,138.0000,0.0000,144.0000,84.0000],
    [32655,144.0000,0.0000,150.0000,84.0000],
    [32656,150.0000,0.0000,156.0000,84.0000],
    [32657,156.0000,0.0000,162.0000,84.0000],
    [32658,162.0000,0.0000,168.0000,84.0000],
    [32659,168.0000,0.0000,174.0000,84.0000],
    [32660,174.0000,0.0000,180.0000,84.0000],
    [32701,-180.0000,-80.0000,-174.0000,0.0000],
    [32702,-174.0000,-80.0000,-168.0000,0.0000],
    [32703,-168.0000,-80.0000,-162.0000,0.0000],
    [32704,-162.0000,-80.0000,-156.0000,0.0000],
    [32705,-156.0000,-80.0000,-150.0000,0.0000],
    [32706,-150.0000,-80.0000,-144.0000,0.0000],
    [32707,-144.0000,-80.0000,-138.0000,0.0000],
    [32708,-138.0000,-80.0000,-132.0000,0.0000],
    [32709,-132.0000,-80.0000,-126.0000,0.0000],
    [32710,-126.0000,-80.0000,-120.0000,0.0000],
    [32711,-120.0000,-80.0000,-114.0000,0.0000],
    [32712,-114.0000,-80.0000,-108.0000,0.0000],
    [32713,-108.0000,-80.0000,-102.0000,0.0000],
    [32714,-102.0000,-80.0000,-96.0000,0.0000],
    [32715,-96.0000,-80.0000,-90.0000,0.0000],
    [32716,-90.0000,-80.0000,-84.0000,0.0000],
    [32717,-84.0000,-80.0000,-78.0000,0.0000],
    [32718,-78.0000,-80.0000,-72.0000,0.0000],
    [32719,-72.0000,-80.0000,-66.0000,0.0000],
    [32720,-66.0000,-80.0000,-60.0000,0.0000],
    [32721,-60.0000,-80.0000,-54.0000,0.0000],
    [32722,-54.0000,-80.0000,-48.0000,0.0000],
    [32723,-48.0000,-80.0000,-42.0000,0.0000],
    [32724,-42.0000,-80.0000,-36.0000,0.0000],
    [32725,-36.0000,-80.0000,-30.0000,0.0000],
    [32726,-30.0000,-80.0000,-24.0000,0.0000],
    [32727,-24.0000,-80.0000,-18.0000,0.0000],
    [32728,-18.0000,-80.0000,-12.0000,0.0000],
    [32729,-12.0000,-80.0000,-6.0000,0.0000],
    [32730,-6.0000,-80.0000,0.0000,0.0000],
    [32731,0.0000,-80.0000,6.0000,0.0000],
    [32732,6.0000,-80.0000,12.0000,0.0000],
    [32733,12.0000,-80.0000,18.0000,0.0000],
    [32734,18.0000,-80.0000,24.0000,0.0000],
    [32735,24.0000,-80.0000,30.0000,0.0000],
    [32736,30.0000,-80.0000,36.0000,0.0000],
    [32737,36.0000,-80.0000,42.0000,0.0000],
    [32738,42.0000,-80.0000,48.0000,0.0000],
    [32739,48.0000,-80.0000,54.0000,0.0000],
    [32740,54.0000,-80.0000,60.0000,0.0000],
    [32741,60.0000,-80.0000,66.0000,0.0000],
    [32742,66.0000,-80.0000,72.0000,0.0000],
    [32743,72.0000,-80.0000,78.0000,0.0000],
    [32744,78.0000,-80.0000,84.0000,0.0000],
    [32745,84.0000,-80.0000,90.0000,0.0000],
    [32746,90.0000,-80.0000,96.0000,0.0000],
    [32747,96.0000,-80.0000,102.0000,0.0000],
    [32748,102.0000,-80.0000,108.0000,0.0000],
    [32749,108.0000,-80.0000,114.0000,0.0000],
    [32750,114.0000,-80.0000,120.0000,0.0000],
    [32751,120.0000,-80.0000,126.0000,0.0000],
    [32752,126.0000,-80.0000,132.0000,0.0000],
    [32753,132.0000,-80.0000,138.0000,0.0000],
    [32754,138.0000,-80.0000,144.0000,0.0000],
    [32755,144.0000,-80.0000,150.0000,0.0000],
    [32756,150.0000,-80.0000,156.0000,0.0000],
    [32757,156.0000,-80.0000,162.0000,0.0000],
    [32758,162.0000,-80.0000,168.0000,0.0000],
    [32759,168.0000,-80.0000,174.0000,0.0000],
    [32760,174.0000,-80.0000,180.0000,0.0000]])
    
    for all in a:
        if point[0]>all[1] and point[0]<all[3] and point[1]>all[2] and point[1]<all[4]:
            epsg = all[0]
    
    return int(epsg)


def FindEPSG2(filepath, output):
    
    if output['satname'][0]=="L5":
        a = glob.glob(filepath + "\\L5\\30m\\" + output['filename'][0])
    elif output['satname'][0]=="L8":
        a = glob.glob(filepath + "\\L8\\pan\\" + output['filename'][0])
    elif output['satname'][0]=="S2":
        a = glob.glob(filepath + "\\S2\\10m\\" + output['filename'][0])
    image_epsg = int(gdal.Info(a[0], format='json')['coordinateSystem']['wkt'].rsplit('"EPSG","', 1)[-1].split('"')[0])
    return image_epsg







def Get_Tidal_Data(base_folder, sitename,datetimeDates,lat,long):

    # get tides
    
    
    hycom_list = [
    'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_56.3?lat,lon,time,surf_el',
    'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.2?lat,lon,time,surf_el',
    'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.8?lat,lon,time,surf_el',
    'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.7?lat,lon,time,surf_el',
    'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.9?lat,lon,time,surf_el',
    'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_93.0?lat,lon,time,surf_el',
    'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0?lat,lon,time,surf_el'
    ]
    
    # hycom date ranges
    #2014-07-01T12:00:00 ... 2016-09-30T09:00:00    
    #2016-05-01T12:00:00 ... 2017-02-01T09:00:00    
    # 2017-02-01T12:00:00 ... 2017-06-01T09:00:00
    # 2017-06-01T12:00:00 ... 2017-10-01T09:00:00
    # 2017-10-01T12:00:00 ... 2018-03-20T09:00:00
    # 2018-09-19T12:00:00 ... 2018-12-09T09:00:00
    # 2018-01-01T12:00:00 ... 2020-02-19T09:00:00

    hycom_date_list = [
    [datetime.datetime(year=2014, month=7,  day=1),    datetime.datetime(year=2016, month=9, day=30)],  
    [datetime.datetime(year=2016, month=5,  day=1),    datetime.datetime(year=2017, month=2, day=1)],
    [datetime.datetime(year=2017, month=2,  day=1),    datetime.datetime(year=2017, month=6, day=1)],
    [datetime.datetime(year=2017, month=6,  day=1),     datetime.datetime(year=2017, month=10, day=1)],
    [datetime.datetime(year=2017, month=10, day=1),    datetime.datetime(year=2018, month=3, day=20)],
    [datetime.datetime(year=2018, month=9,  day=19),    datetime.datetime(year=2018, month=12, day=9)],
    [datetime.datetime(year=2018, month=1, day=1),    datetime.datetime(year=2020, month=2, day=19)]
    ]
    
    FIRST = True
    hycom_list_flag = []
    for n in [0,1,2,3,4,5,6]:
        if hycom_date_list[n][0] <= datetimeDates[0] <= hycom_date_list[n][1] and FIRST:
            hycom_list_flag.append(n)
            FIRST = False
        if hycom_date_list[n][0] <= datetimeDates[1] <= hycom_date_list[n][1] and not FIRST:
            hycom_list_flag.append(n)
            break
    
    list_range = np.arange(hycom_list_flag[0],hycom_list_flag[1]+1)
    
    
    
    
    elevations = np.array([])
    elevations_dates = np.array([], dtype='datetime64[ns]')
    
    for n in list_range:
        ds = xr.open_dataset(hycom_list[n])
        if long < 0:
            lon = -1*((-1*long)-180-180)
        else:
            lon = long
        
        lat2 = ds['lat'].sel(lat=lat,method='nearest')
        long2 = ds['lon'].sel(lon=lon,method='nearest')
        
        elevations = np.append(elevations, np.array(ds['surf_el'].sel(lat=lat2, lon=long2) ))
        elevations_dates= np.append(elevations_dates, np.array(ds['time'].sel()) )
    
    unique_dates = np.unique(elevations_dates,return_index=True)[1]
    elevations = elevations[unique_dates]
    elevations_dates= elevations_dates[unique_dates]
    
    
    cc_save.Save_Tides(base_folder, sitename, elevations, elevations_dates)
    return elevations,elevations_dates










def Get_Tidal_DataNOAA(base_folder, sitename,output,lat,long,utc):
    
    datetimeDates = output['dates']
    
    gages = np.genfromtxt(base_folder+'gages.txt',dtype=[('S50'),('S12'),('f8'),('f8')],delimiter=',')
    
    distances = []
    for line in gages:
        distances.append(distancePP([line[2],line[3]],[lat,long]))
    a = np.where(np.array(distances)== np.nanmin(distances))[0][0]
    gage = gages[a]
    print('Is gage at ' + str(gage[0]) + ' close enough?')
    station = str(gage[1])[2:-1]
    
    
    datums = ['msl','mhhw','mhw','mtl','msl','mlw','mllw','navd']
    n=-1
    datumm = False
    while not datumm:
        n=n+1
        try:
            datum = datums[n]
        except IndexError:
            break
        tides = []
        tidedates = []
        
        for date in datetimeDates:
            # station = '8454000'
            rangee=False
            if rangee==True:
                beginDateStr = '20130101%2000:00'
                endDateStr = '20130102%2000:00'
                dateStr = 'begin_date='+beginDateStr+'&end_date='+endDateStr
            else:
                
                onedayless = ((date-datetime.timedelta(days=1)))
                onedaymore = ((date+datetime.timedelta(days=1)))
                
                # less
                if onedayless.month<10:
                    month = '0'+str(onedayless.month)
                else:
                    month = str(onedayless.month)
                    
                    
                if onedayless.day<10:
                    day = '0'+str(onedayless.day)
                else:
                    day = str(onedayless.day)
                beginDateStr = str(onedayless.year)+month+day+'%2000:00'
                
                # more
                
                
                if onedaymore.month<10:
                    month = '0'+str(onedaymore.month)
                else:
                    month = str(onedaymore.month)
                    
                    
                if onedaymore.day<10:
                    day = '0'+str(onedaymore.day)
                else:
                    day = str(onedaymore.day)
                endDateStr = str(onedaymore.year)+month+day+'%2000:00'
                    
                    
                    
                dateStr = 'begin_date='+beginDateStr+'&end_date='+endDateStr
                
                
                
                
            noaaCode = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?'+dateStr+'&station='+station+'&product=water_level&datum='+datum+'&units=metric&time_zone=gmt&application=web_services&format=xml'
            dateStr = 'begin_date='+beginDateStr[:-8]+'&end_date='+endDateStr[:-8]
            noaaCode = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=predictions&application=usace&'+dateStr+'&station='+station+'&product=water_level&datum='+datum+'&units=metric&time_zone=gmt&application=web_services&format=xml'


            try:    
                with urllib.request.urlopen(noaaCode) as f:
                    html = f.read().decode('utf-8')
                datumm = True
            except:
                break
                
        
            html2 = html.split('pr')
            dates =[]
            vs = []
            for line in html2[1:-1]:
                parts = line.split('  ')
                for part in parts:
                    if 't=' in part:
                        time = ''
                        for i in part:
                            if i.isdigit():
                                time=time+i
                    if 'v=' in part:
                        v = ''
                        for i in part:
                            if i.isdigit() or i=='.' or i=='-':
                                v=v+i
                dates.append(time)
                vs.append(float(v))


            datetimes = []
            for datee in dates:
                datetimes.append(datetime.datetime(year=int(datee[0:4]),month=int(datee[4:6]),day=int(datee[6:8]),hour=int(datee[8:10]),minute=int(datee[10:12]),tzinfo=utc))
        
        
        
            
            # html2 = html.split('<wl')
            # html3 = html2[1:]
            # html3[-1] = html3[-1].split('</observations>')[0]
            # datees = []
            # datetimes = []
            # vs = []
            # ss = []
            # for line in html3:
            #     a = line.split('\" ')
            #     datee = a[0][4:]
            #     v = a[1][4:]
            #     s = a[2][3:]
            #     datees.append(datee)
                
            #     try:
            #         vs.append(float(v))
            #     except:
            #         vs.append(0)
                
            #     try:
            #         ss.append(float(s))
            #     except:
            #         ss.append(0)
            #     datetimes.append(datetime.datetime(year=int(datee.split('-')[0]),month=int(datee.split('-')[1]),day=int(datee.split('-')[2].split(' ')[0]),hour=int(datee.split('-')[2].split(' ')[1].split(':')[0]),minute=int(datee.split('-')[2].split(' ')[1].split(':')[1]),tzinfo=utc))
        
            a = np.where(np.array(datetimes)==nearest(np.array(datetimes),date))[0][0]
            tides.append(vs[a])
            tidedates.append(date)
    
    
    output['tides'] = tides
    output['tidedates'] = tidedates
    
    
    return output










def Get_Tidal_DataNOAA2(base_folder, sitename,output,lat,long,utc):
    
    datetimeDates = output['dates']
    
    gages = np.genfromtxt(base_folder+'gages.txt',dtype=[('S50'),('S12'),('f8'),('f8')],delimiter=',')
    
    distances = []
    for line in gages:
        distances.append(distancePP([line[2],line[3]],[lat,long]))
    a = np.where(np.array(distances)== np.nanmin(distances))[0][0]
    gage = gages[a]
    print('Is gage at ' + str(gage[0]) + ' close enough?')
    station = str(gage[1])[2:-1]
    
    
    datum = 'msl'
    tides = []
    tidedates = []
    
    for date in datetimeDates:
        # station = '8454000'
        rangee=False
        if rangee==True:
            beginDateStr = '20130101%2000:00'
            endDateStr = '20130102%2000:00'
            dateStr = 'begin_date='+beginDateStr+'&end_date='+endDateStr
        else:
            
            onedayless = ((date-datetime.timedelta(days=1)))
            onedaymore = ((date+datetime.timedelta(days=1)))
            
            # less
            if onedayless.month<10:
                month = '0'+str(onedayless.month)
            else:
                month = str(onedayless.month)
                
                
            if onedayless.day<10:
                day = '0'+str(onedayless.day)
            else:
                day = str(onedayless.day)
            beginDateStr = str(onedayless.year)+month+day+'%2000:00'
            
            # more
            
            
            if onedaymore.month<10:
                month = '0'+str(onedaymore.month)
            else:
                month = str(onedaymore.month)
                
                
            if onedaymore.day<10:
                day = '0'+str(onedaymore.day)
            else:
                day = str(onedaymore.day)
            endDateStr = str(onedaymore.year)+month+day+'%2000:00'
                
                
                
            dateStr = 'begin_date='+beginDateStr+'&end_date='+endDateStr
            
            
            
            
        noaaCode = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?'+dateStr+'&station='+station+'&product=water_level&datum='+datum+'&units=metric&time_zone=gmt&application=web_services&format=xml'
        dateStr = 'begin_date='+beginDateStr[:-8]+'&end_date='+endDateStr[:-8]
        noaaCode = 'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=predictions&application=usace&'+dateStr+'&station='+station+'&product=water_level&datum='+datum+'&units=metric&time_zone=gmt&application=web_services&format=xml'


        vs = np.genfromtxt(noaaCode, skip_header=4, skip_footer=10,delimiter='\"', usecols=3)            
        
        madeupdates = [datetime.datetime(year=onedayless.year,month=onedayless.month,day=onedayless.day,hour=0,minute=0,tzinfo=utc)]
        for all in vs[1:]:
            madeupdates.append(madeupdates[-1] + datetime.timedelta(minutes=6))
    

        a = np.where(np.array(madeupdates)==nearest(np.array(madeupdates),date))[0][0]
        tides.append(vs[a])
        tidedates.append(date)
    
    
    output['tides'] = tides
    output['tidedates'] = tidedates
    
    
    return output






    
    key = input('1-GOOD    0-BAD')
    if key =='1':
        pass
    elif key =='0':
        shutil.move(aaa, site+"\\bad\\"+aaa.split('\\')[-1])
    return key



def MapTideGagesAllExistingSites(base_folder, sitename,output,lat,long):
    datetimeDates = output['dates']
    
    gages = np.genfromtxt(base_folder+'gages.txt',dtype=[('S50'),('S12'),('f8'),('f8')],delimiter=',')
    
    distances = []
    gagelats = []
    gagelongs = []
    for line in gages:
        distances.append(distancePP([line[2],line[3]],[lat,long]))
        gagelats.append(line[2])
        gagelongs.append(line[3])
        
    a = np.where(np.array(distances)== np.nanmin(distances))[0][0]
    
    plt.figure(figsize=(10,10))
    
    
    sf = shapefile.Reader(base_folder+"\\classification\\NOAAUSAShoreline\\usaShore.shp")
    shapes = sf.shapes()
    for shape in shapes:
        plt.scatter(np.array(shape.points)[:,0], np.array(shape.points)[:,1], c='black',s=1)
    plt.scatter(gagelongs, gagelats, s=30,c='cyan')
    
    n=-1
    for all in gagelongs:
        n=n+1
        if (gagelongs[n]>long-0.5 and gagelongs[n]<long+0.5) and (gagelats[n]>lat-0.5 and gagelats[n]<lat+0.5):
            plt.text(gagelongs[n], gagelats[n], str(n))
    
    
    
    sf.close()
    
    plt.scatter(gages[a][3],gages[a][2], s=70,c='blue')
    plt.scatter(long,lat,s=70,c='red')
    plt.axis([long-0.5,long+0.5,lat-0.5,lat+0.5])
#        plt.show()
    plt.savefig(base_folder + 'data\\' + sitename + '\\TideStationSelection.png')
    plt.close('all')
    
    
    
    print('Is th highlighted gage close enough?')
    key = input('Type \"yes\" or \"no\"')
    if key =='yes':
        pass
    elif key =='no':
        a = input('Which gage should be used?')
        
    gage = gages[a]
    print('Using gage', a,',',gage[0])





def Get_Tidal_DataNOAA3(base_folder, sitename,output,lat,long,utc, a,tideType, USshape, gagesloc):
    datetimeDates = output['dates']
    
    gages = np.genfromtxt(gagesloc,dtype=[('S50'),('S12'),('f8'),('f8')],delimiter=',')
    
    distances = []
    gagelats = []
    gagelongs = []
    for line in gages:
        distances.append(distancePP([line[2],line[3]],[lat,long]))
        gagelats.append(line[2])
        gagelongs.append(line[3])
        
        
        
        
    if int(a) == -1:
        a = np.where(np.array(distances)== np.nanmin(distances))[0][0]
        
        plt.figure(figsize=(10,10))
        
        
        sf = shapefile.Reader(USshape)
        shapes = sf.shapes()
        for shape in shapes:
            plt.scatter(np.array(shape.points)[:,0], np.array(shape.points)[:,1], c='black',s=1)
        plt.scatter(gagelongs, gagelats, s=30,c='cyan')
        
        n=-1
        for all in gagelongs:
            n=n+1
            if (gagelongs[n]>long-0.5 and gagelongs[n]<long+0.5) and (gagelats[n]>lat-0.5 and gagelats[n]<lat+0.5):
                plt.text(gagelongs[n], gagelats[n], str(n))
        
        
        
        sf.close()
        
        plt.scatter(gages[a][3],gages[a][2], s=70,c='blue')
        plt.scatter(long,lat,s=70,c='red')
        plt.axis([long-0.5,long+0.5,lat-0.5,lat+0.5])
#        plt.show()
        plt.savefig(base_folder + 'data\\' + sitename + '\\TideStationSelection.png')
        plt.close('all')
        
        
        
        # print('Is the highlighted gage close enough?')
        # key = input('Type \"yes\" or \"no\"')
        # if key =='yes':
            # pass
        # elif key =='no':
            # a = input('Which gage should be used?')
        
        
        
        
        
    gage = gages[a]
    print('Using gage', a,',',gage[0])
    station = str(gage[1])[2:-1]
    
    
    datum = 'msl'
    tides = []
    tidedates = []
    alltides = []
    alltidedates = []
    
    for date in datetimeDates:
        # station = '8454000'
        rangee=False
        if rangee==True:
            beginDateStr = '20130101%2000:00'
            endDateStr = '20130102%2000:00'
            dateStr = 'begin_date='+beginDateStr+'&end_date='+endDateStr
        else:
            
            onedayless = ((date-datetime.timedelta(days=1)))
            onedaymore = ((date+datetime.timedelta(days=1)))
            
            # less
            if onedayless.month<10:
                month = '0'+str(onedayless.month)
            else:
                month = str(onedayless.month)
                
                
            if onedayless.day<10:
                day = '0'+str(onedayless.day)
            else:
                day = str(onedayless.day)
            beginDateStr = str(onedayless.year)+month+day+'%2000:00'
            
            # more
            
            
            if onedaymore.month<10:
                month = '0'+str(onedaymore.month)
            else:
                month = str(onedaymore.month)
                
                
            if onedaymore.day<10:
                day = '0'+str(onedaymore.day)
            else:
                day = str(onedaymore.day)
            endDateStr = str(onedaymore.year)+month+day+'%2000:00'
                
                
                
            dateStr = 'begin_date='+beginDateStr+'&end_date='+endDateStr
            
            
        dateStr = 'begin_date='+beginDateStr[:-8]+'&end_date='+endDateStr[:-8]
        noaaCode = 'https://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&'+dateStr+'&datum='+tideType+'&station='+station+'&time_zone=GMT&units=metric&interval=h&format=xml'
    
        
        data = urllib.request.urlopen(noaaCode)    
        data2 = []
        for line in data:
            data2.append(line)
        data3 = data2[4:-7]
        vs = []
        for line in data3:
            bla = str(line)[31:-7]
            vs.append(float(bla))
        
        
        madeupdates = [datetime.datetime(year=onedayless.year,month=onedayless.month,day=onedayless.day,hour=0,minute=0,tzinfo=utc)]
        for all in vs[1:]:
            madeupdates.append(madeupdates[-1] + datetime.timedelta(minutes=60))
    
        a = np.where(np.array(madeupdates)==nearest(np.array(madeupdates),date))[0][0]
        
        
        tides.append(vs[a])
        tidedates.append(date)
        
#        print(vs[0])
    
        alltides = alltides +vs
        alltidedates = alltidedates + madeupdates
    
    
    
    
    
    if lat > 40 and lat < 50 and long > -92 and long < -75:
        tides = np.zeros(len(tides))
    
    
    
    output['tides'] = tides
    output['tidedates'] = tidedates
    
#    plt.figure(figsize=(20,20))
#    plt.plot(alltidedates, alltides)
#    plt.show()
#    plt.close('all')
    
    return output













def Get_Wave_Data(base_folder, sitename, datetimeDates,long,lat,utc):

    print('check database for updates???')
    print('Last available is 5/2019')
    
    yearList = np.arange(datetimeDates[0].year, datetimeDates[1].year+1)
    
    noaaList = []
    for year in yearList:
        if year <=2016:
            for month in ['01','02','03','04','05','06','07','08','09','10','11','12']:
                noaaList.append('https://data.nodc.noaa.gov/thredds/dodsC/ncep/nww3/'+str(year)+'/'+month+'/glo_30m/multi_1.glo_30m.hs.'+str(year)+month+'.grb2?lat,lon,time,Significant_height_of_combined_wind_waves_and_swell_surface')
        if year ==2017:
            for month in ['01','02','03','04','05','06','07','08','09','10','11','12']:
                noaaList.append('https://data.nodc.noaa.gov/thredds/dodsC/ncep/nww3/'+str(year)+'/'+month+'/gribs/multi_1.glo_30m.hs.'+str(year)+month+'.grb2?lat,lon,time,Significant_height_of_combined_wind_waves_and_swell_surface')
        if year ==2018:
            for month in ['01','02','03','04','05','06','07','08','09','10','11','12']:
                noaaList.append('https://data.nodc.noaa.gov/thredds/dodsC/ncep/nww3/'+str(year)+'/'+month+'/gribs/multi_1.glo_30m.hs.'+str(year)+month+'.grb2?lat,lon,time,Significant_height_of_combined_wind_waves_and_swell_surface')
        if year ==2019:
            for month in ['01','02','03','04','05']:
                noaaList.append('https://data.nodc.noaa.gov/thredds/dodsC/ncep/nww3/'+str(year)+'/'+month+'/gribs/multi_1.glo_30m.hs.'+str(year)+month+'.grb2?lat,lon,time,Significant_height_of_combined_wind_waves_and_swell_surface')
    
    
    if long < 0:
        lon = -1*((-1*long)-180-180)
    else:
        lon = long
            
        
    hsTS = np.array([])
    tpTS = np.array([])
    dpTS = np.array([])
    wave_timeTS = np.array([], dtype='datetime64[ns]')
    for noaa in noaaList:
        # hs, time
        ds = xr.open_dataset(noaa)
        lat2 = ds['lat'].sel(lat=lat,method='nearest')
        long2 = ds['lon'].sel(lon=lon,method='nearest')
        hs = np.array(ds['Significant_height_of_combined_wind_waves_and_swell_surface'].sel(lon=long2, lat=lat2))
        hsTS = np.append(hsTS, hs)
        wave_time = np.array(ds['time'].sel())
        wave_timeTS = np.append(wave_timeTS, wave_time)
        # tp
        noaa2 = noaa.replace('hs','tp')
        noaa2 = noaa2.replace('Significant_height_of_combined_wind_waves_and_swell_surface','Primary_wave_mean_period_surface')
        ds = xr.open_dataset(noaa2)
        lat2 = ds['lat'].sel(lat=lat,method='nearest')
        long2 = ds['lon'].sel(lon=lon,method='nearest')
        tp = np.array(ds['Primary_wave_mean_period_surface'].sel(lon=long2, lat=lat2))
        tpTS = np.append(tpTS, tp)
        # dp
        noaa3 = noaa2.replace('.tp.','.dp.')
        noaa3 = noaa3.replace('Primary_wave_mean_period_surface','Primary_wave_direction_surface')
        ds = xr.open_dataset(noaa3)
        lat2 = ds['lat'].sel(lat=lat,method='nearest')
        long2 = ds['lon'].sel(lon=lon,method='nearest')
        dp = np.array(ds['Primary_wave_direction_surface'].sel(lon=long2, lat=lat2))
        dpTS = np.append(dpTS, dp)
        

        
    wave_timeTS2 = []
    for date in wave_timeTS:
        ts = (date - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
        # wave_timeTS2.append( datetime.datetime.utcfromtimestamp(ts).replace(tzinfo=output['dates'][0].tzinfo) )
        wave_timeTS2.append( datetime.datetime.utcfromtimestamp(ts).replace(tzinfo=utc) )
    wave_timeTS = np.array(wave_timeTS2)    
        
        
    
    
    cc_save.Save_Waves(base_folder, sitename, hsTS, tpTS, dpTS, wave_timeTS)
    return hsTS, tpTS, dpTS, wave_timeTS



def Process_Wave_Data(base_folder, sitename, hsTS, tpTS, dpTS, wave_timeTS):    
    
    # hsTS = np.loadtxt(base_folder+'data\\'+sitename+'\\waveheights.txt')
    # tpTS = np.loadtxt(base_folder+'data\\'+sitename+'\\waveperiods.txt')
    # dpTS = np.loadtxt(base_folder+'data\\'+sitename+'\\wavedirections.txt')
    # with open(base_folder+'data\\'+sitename+'\\wavedates.pkl', 'rb') as f:
    #     wave_timeTS = pickle.load(f) 

    
    yeears = []
    mnths = []
    for date in wave_timeTS:
        yeears.append(date.year)
        mnths.append(date.month)
    yeears = np.array(yeears)
    mnths = np.array(mnths)
    
    # Ensure clean data
    a = np.where((hsTS!=np.nan)&(tpTS!=np.nan))
    heightst_clean = hsTS[a]
    periodst_clean = tpTS[a]
    directionst_clean = dpTS[a]
    timest2_clean = wave_timeTS[a[0]]
    years_clean = yeears[a]
    months_clean = mnths[a]
    
    
    #  #plot of monthly mean wave heights, stdevs, periods
    wave_heights_monthly_means = []
    wave_heights_monthly_stdevs = []
    wave_periods_monthly_means = []
    
    for month in range(1,13):
        a = np.where(months_clean == month)
        wave_heights_monthly_means.append(np.nanmean(heightst_clean[a]))
        wave_heights_monthly_stdevs.append(np.std(heightst_clean[a]))
        wave_periods_monthly_means.append(np.nanmean(periodst_clean[a]))
    
    
    
    
    # plot of monthly 1% & 10% exceedance
    
    wave_heights_monthly_exceed2 = []
    wave_heights_monthly_exceed10 = []
    for month in range(1,13):
        a = np.where(months_clean == month)
        monthwaves = heightst_clean[a]
        a = st.rayleigh
        params = a.fit(monthwaves)
        x = np.arange(0,6,.01)
        probs = a.pdf(x, loc = params[0], scale = params[1])
        cumprobs = a.cdf(x, loc = params[0], scale = params[1])
        
        prob2 = nearest(cumprobs, .99)
        prob10 = nearest(cumprobs, .9)
        Hs2 = x[np.where(cumprobs == prob2)]
        Hs10 = x[np.where(cumprobs == prob10)]
        wave_heights_monthly_exceed2.append(Hs2)
        wave_heights_monthly_exceed10.append(Hs10)
        
        
    
    
    
    
    # Seasonal Wave Frequency Compass
    
    
    
    n=-1
    
    WaveFrequencyByDirectionBandsSeasonally = []
    meanWaveHeightsByDirectionBandsSeasonally = []
    
    
    for month in [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]:
        n=n+1
        a = np.where(((months_clean == month[0])|(months_clean == month[1])|(months_clean == month[2])))
        monthwaves = heightst_clean[a]
        monthdirections = directionst_clean[a]
        bandwidth=10
        directionalBands = np.arange(0,360+bandwidth,bandwidth)
        
        WaveFrequencyByDirectionBands = []
        meanWaveHeightsByDirectionBands = []
        for a in range(0, len(directionalBands)-1):
            bandLow = directionalBands[a]
            bandHigh = directionalBands[a+1]
            a = np.where((monthdirections>=bandLow)&(monthdirections<bandHigh))
            WaveFrequencyByDirectionBands.append(len(monthwaves[a]))
            meanWaveHeightsByDirectionBands.append(np.nanmean(monthwaves[a]))
        WaveFrequencyByDirectionBandsSeasonally.append(WaveFrequencyByDirectionBands)
        meanWaveHeightsByDirectionBandsSeasonally.append(meanWaveHeightsByDirectionBands)
            
        
        
    
    
    
    
    # Seasonal BIG Wave Frequency Compass
    BigThreshold = np.nanmean(wave_heights_monthly_exceed10)
    
    n=-1
    
    WaveFrequencyByDirectionBandsSeasonallyBIG = []
    meanWaveHeightsByDirectionBandsSeasonallyBIG = []
    
    
    for month in [[1,2,3],[4,5,6],[7,8,9],[10,11,12]]:
        n=n+1
        a = np.where(((months_clean == month[0])|(months_clean == month[1])|(months_clean == month[2]))&(heightst_clean>=BigThreshold))
        monthwaves = heightst_clean[a]
        monthdirections = directionst_clean[a]
        bandwidth=10
        directionalBands = np.arange(0,360+bandwidth,bandwidth)
        
        WaveFrequencyByDirectionBands = []
        meanWaveHeightsByDirectionBands = []
        for a in range(0, len(directionalBands)-1):
            bandLow = directionalBands[a]
            bandHigh = directionalBands[a+1]
            a = np.where((monthdirections>bandLow)&(monthdirections<=bandHigh))
            WaveFrequencyByDirectionBands.append(len(monthwaves[a]))
            meanWaveHeightsByDirectionBands.append(np.nanmean(monthwaves[a]))
            
            
            
            
        WaveFrequencyByDirectionBandsSeasonallyBIG.append(WaveFrequencyByDirectionBands)
        meanWaveHeightsByDirectionBandsSeasonallyBIG.append(meanWaveHeightsByDirectionBands)


    Wave_Data = {'meanWaveHeightsByDirectionBandsSeasonallyBIG':meanWaveHeightsByDirectionBandsSeasonallyBIG}
    Wave_Data['wave_heights_monthly_means']=wave_heights_monthly_means
    Wave_Data['wave_heights_monthly_stdevs']=wave_heights_monthly_stdevs
    Wave_Data['wave_heights_monthly_exceed2']=wave_heights_monthly_exceed2
    Wave_Data['wave_heights_monthly_exceed10']=wave_heights_monthly_exceed10
    Wave_Data['wave_periods_monthly_means']=wave_periods_monthly_means
    Wave_Data['WaveFrequencyByDirectionBandsSeasonally']=WaveFrequencyByDirectionBandsSeasonally
    Wave_Data['directionalBands']=directionalBands
    Wave_Data['meanWaveHeightsByDirectionBandsSeasonally']=meanWaveHeightsByDirectionBandsSeasonally
    Wave_Data['WaveFrequencyByDirectionBandsSeasonallyBIG']=WaveFrequencyByDirectionBandsSeasonallyBIG


    # return meanWaveHeightsByDirectionBandsSeasonallyBIG,wave_heights_monthly_means,wave_heights_monthly_stdevs,wave_heights_monthly_exceed2,wave_heights_monthly_exceed10,wave_periods_monthly_means,WaveFrequencyByDirectionBandsSeasonally,directionalBands,meanWaveHeightsByDirectionBandsSeasonally,WaveFrequencyByDirectionBandsSeasonallyBIG
    return Wave_Data









def RemoveBadImageryLater(base_folder, sitename, output, metadata):
    

    
    
    
    bad_list = glob.glob(base_folder+'data\\'+sitename+'\\jpg_files\\detection\\bad\\*')
    newbads = []
    for bad in bad_list:
        newbad = bad.split('\\')[-1][:16]
        newbads.append(newbad)
    locations = []
    n=-1
    for file in output['filename']:
        n=n+1
        compare = file[:16]
        if not compare in newbads:
            locations.append(n)
    
        
    outputNew = {}
    for key in output:
        listNew = []
        n=-1
        for entry in output[key]:
            n=n+1
            if n in locations:
                listNew.append(entry)
    
        outputNew[key]=listNew
        
        
        
    metadataNew = {}
    for key in metadata:
        satKey = {}
        filenames=[]
        acc_georef=[]
        epsg=[]
        dates=[]
        n=-1
        for date in metadata[key]['dates']:
            n=n+1
            if date in outputNew['dates']:
                filenames.append(metadata[key]['filenames'][n])
                acc_georef.append(metadata[key]['acc_georef'][n])
                epsg.append(metadata[key]['epsg'][n])
                dates.append(metadata[key]['dates'][n])
        outputNew[key]=listNew
        
        satKey['filenames']=filenames
        satKey['acc_georef']=acc_georef
        satKey['epsg']=epsg
        satKey['dates']=dates
    metadataNew[key]=satKey
    
    
    
    
    return outputNew, metadataNew
    
    














def geojsonTOshp_lines(geojson, shp):
    a = "ogr2ogr -nlt LINESTRING -skipfailures -overwrite "
    b = shp + " "
    c = geojson
    # d = " "
    # d = "OGRGeoJSON"
    clipString = a + b + c# + d + e
    subprocess.call(clipString)
    


def shapefileTOgeoJSON_lines(geojson, shp):
    # a = "ogr2ogr -f \"GeoJSON\" -t_srs "
    # b = shp + " "
    # c = geojson + " "
    # # d = " "
    # # d = "ESRI Shapefile"
    # clipString = a + c + b# + d + e
    # # print(clipString)
    # subprocess.call(clipString)    
    a = "ogr2ogr -nlt LINESTRING -skipfailures "
    b = geojson + " "
    c = shp
    # d = " "
    # d = "OGRGeoJSON"
    clipString = a + b + c# + d + e
    subprocess.call(clipString)



def grayscale(img):
    grayImage = np.zeros(img.shape)
    R = np.array(img[:, :, 0])
    G = np.array(img[:, :, 1])
    B = np.array(img[:, :, 2])
    a = np.array(img[:, :, 2])

    R = (R *0.75)
    G = (G *0.5)
    B = (B *0.25)
    a = (a *0)

    Avg = (R+G+B)
    return Avg


def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))



def distancePP(point1, point2):
    distances=(np.sqrt((point1[0]-point2[0])**2+(point1[1]-point2[1])**2))
    return distances




def distance(line, point):
    distances = []
    for all in line:
        distances.append(np.sqrt((all[0]-point[0])**2+(all[1]-point[1])**2))
    return distances


def getArrayValueScaled(x, y, array, lry,uly,ulx,lrx):
    xpercent = (x-ulx)/(lrx-ulx)
    ypercent = (y-lry)/(uly-lry)
    ypercent = 1-ypercent
    C = xpercent*array.shape[1]
    R = ypercent*array.shape[0]
    C = int(C)
    R = int(R)
    try:
        arrayValue = array[R,C]
    except:
        arrayValue = 0
    return arrayValue




def imageExtracter(im):
    a = np.shape(im)
    b = int(a[0]/2)

    middleRow = im[b,:]
    
    oneThirdColumn = int((1/3)*a[1])
    
    im2 = im[:,oneThirdColumn:oneThirdColumn*2,:]
    
    comparison = im2[:,0]
    im3=[]
    for n in range(0,im2.shape[1]):
        col = im2[:,n]
        if not np.array_equal(col, comparison):
            im3.append(col)
    im3 = np.array(im3)            
    # im3 = im3.reshape((im2.shape[0],-1,4))
    
    comparison = im3[:,0]
    im4=[]
    for n in range(0,im3.shape[1]):
        col = im3[:,n]
        if not np.array_equal(col, comparison):
            im4.append(col)
    im4 = np.array(im4)            
    # im4 = im4.reshape((im3.shape[0],-1,4))
    
    
    return im4    




def imageExtracterS2(im2):
    
    comparison = im2[:,0]
    im3=[]
    for n in range(0,im2.shape[1]):
        col = im2[:,n]
        if not np.array_equal(col, comparison):
            im3.append(col)
    im3 = np.array(im3)            
    # im3 = im3.reshape((im2.shape[0],-1,4))
    
    comparison = im3[:,0]
    im4=[]
    for n in range(0,im3.shape[1]):
        col = im3[:,n]
        if not np.array_equal(col, comparison):
            im4.append(col)
    im4 = np.array(im4)            
    # im4 = im4.reshape((im3.shape[0],-1,4))
    
    
    return im4    
    






def imageExtracter_thirdImage(im):
    a = np.shape(im)
    b = int(a[0]/2)

    middleRow = im[b,:]
    
    oneThirdColumn = int((1/3)*a[1])
    twoThirdColumn = int((2/3)*a[1])
    
    im2 = im[:,twoThirdColumn:oneThirdColumn*3,:]
    
    comparison = im2[:,0]
    im3=[]
    for n in range(0,im2.shape[1]):
        col = im2[:,n]
        if not np.array_equal(col, comparison):
            im3.append(col)
    im3 = np.array(im3)            
    # im3 = im3.reshape((im2.shape[0],-1,4))
    
    comparison = im3[:,0]
    im4=[]
    for n in range(0,im3.shape[1]):
        col = im3[:,n]
        if not np.array_equal(col, comparison):
            im4.append(col)
    im4 = np.array(im4)            
    # im4 = im4.reshape((im3.shape[0],-1,4))
    
    
    return im4    










def NewTransectsBuilder_StartParallel(sitename, base_folder):
    
    filepath = base_folder+'data\\'+sitename
    
    
    
    
    # make shorelines into shapefile
    shores_geo_filename = filepath+"\\"+sitename+"_output.geojson"
    shores_shp_filename = filepath+"\\"+sitename+"_output.shp"
    geojsonTOshp_lines(shores_geo_filename, shores_shp_filename)
    
    
    
    # get image
    images = glob.glob(filepath+"\\jpg_files\\detection\\*.png")
    im = plt.imread(images[0])
    
    c = imageExtracter(im)    
    
    array_gray = grayscale(c)
    
    
    # print(c.shape)
    # print(array_gray.shape)
    # print(array_gray)
    # plt.figure(figsize=(30,30))
    # plt.imshow(array_gray)
    # plt.show()
    # plt.close('all')
    
    
    
    # np.savetxt(r"C:\Users\RDCHLNRO\Desktop\grey.txt", array_gray)######################################################################
    
    
    top = array_gray[1,:]
    right = array_gray[:,-2]
    bottom = array_gray[-2,:]
    left = array_gray[:,1]
    
    top2 = len(np.where(top==0.4284314)[0])/len(top)
    left2 = len(np.where(left==0.4284314)[0])/len(left)
    bottom2 = len(np.where(bottom==0.4284314)[0])/len(bottom)
    right2 = len(np.where(right==0.4284314)[0])/len(right)
    
    test = [top2, left2, bottom2, right2]
    
    open_water_position = np.where(test == np.nanmax(test))[0][0]
    
    

    
    # open a rster to get extent for transects
    try:
        rasters = glob.glob(filepath+"\\S2\\10m\\*_"+sitename+"_10m.tif")
        src = gdal.Open(rasters[0])
    except:
        rasters = glob.glob(filepath+"\\L8\\pan\\*.tif")
        src = gdal.Open(rasters[0])
    
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    src = None
    
    tryagain = 1
    while lrx-ulx < 1000:
        try:
            src = gdal.Open(rasters[tryagain])
            ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
            lrx = ulx + (src.RasterXSize * xres)
            lry = uly + (src.RasterYSize * yres)
            src = None
            tryagain = tryagain+1
        except:
            break
    
    
    
    
    # make reference shoreline into shapefile
    shoreRef_geo_Filename = filepath+"\\"+sitename+"_reference_shoreline.geojson"
    shoreRef_shp_filename = filepath+"\\"+sitename+"_reference_shoreline.shp"
    geojsonTOshp_lines(shoreRef_geo_Filename, shoreRef_shp_filename)
    
    # get reference shapefile points and then get a subset of those
    sf = shapefile.Reader(shoreRef_shp_filename)
    shapes = sf.shapes()
    points = shapes[0].points
    # less_points = points[1::70]
    sf.close()
    
    # add prpidicular points to new transect shapefile
    
    w = shapefile.Writer(filepath+"\\transects.shp")
    w.field('name', 'C')
    transects = dict([])
    transectsFull = dict([])
    transectsRotations = []
    transectsSlopes = []
    COUNTERtransectsfull=0
    COUNTERtransectssmall=0
    COUNTERtransectssmall2=0
    for n in range(5,len(points)-5):
        try:
            COUNTERtransectsfull=COUNTERtransectsfull+1
            COUNTERtransectssmall=COUNTERtransectssmall+1
            
            # for the chosen point along the reference shoreline find a slope for best fit line
            deriv_points = np.array(points[n-2:n+2])
            x = deriv_points[:,0]
            y = deriv_points[:,1]
            X = np.array(points)[n,0]
            Y = np.array(points)[n,1]
            slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
            
            # find a slope perpindicular to shorelne slope for transect
            pslope = -1*(1/slope)
            
            # start transect with first value
            XX = [X]
            YY = [Y + pslope*(XX[-1]-X)]
            
            # grow transect in forward X direction
            waterPointsForward = 0
            PointsForward = 0
            i = 0
            # greyClassification2 = [0]######################################################################
            while XX[-1]<lrx-3 and YY[-1]<uly-3 and YY[-1]>lry-3:
                i=i+1
                XX.append(XX[-1]+1)
                YY.append(Y + pslope*(XX[-1]-X))
                greyClassification = getArrayValueScaled(XX[-1], YY[-1], array_gray, lry,uly,ulx,lrx)
                # greyClassification2.append(greyClassification)######################################################################
                if ((greyClassification > 0.4284313 and greyClassification < 0.4284315) or (greyClassification > 1.34 and greyClassification < 1.36)):
                    waterPointsForward = waterPointsForward+1
                PointsForward = PointsForward+1
            # grow transect in backward X direction
            waterPointsReverse = 0
            PointsReverse = 0
            i = 0
            while XX[0]>ulx-3 and YY[0]<uly-3 and YY[0]>lry-3:
                i=i+1
                XX.insert(0, XX[0]-1)
                YY.insert(0, Y + pslope*(XX[0]-X))
                greyClassification = getArrayValueScaled(XX[0], YY[0], array_gray, lry,uly,ulx,lrx)
                # greyClassification2.insert(0, greyClassification)######################################################################
                if ((greyClassification > 0.4284313 and greyClassification < 0.4284315) or (greyClassification > 1.34 and greyClassification < 1.36)):
                    waterPointsReverse = waterPointsReverse+1
                PointsReverse = PointsReverse+1
                

            if (waterPointsReverse/PointsReverse) > (waterPointsForward/PointsForward):
                # cccc = 'black'
                XX = np.flip(np.array(XX))
                YY = np.flip(np.array(YY))
                

            
            raw_transect_line = []
            n=-1
            for all in YY:
                n=n+1
                raw_transect_line.append((XX[n],YY[n]))
            
            
            distances = distance(raw_transect_line, (X,Y)) 
            raw_transect_line = np.array(raw_transect_line)
            transect_start_point = np.where(distances==nearest(distances,100))[0][0]
            transect_processed = raw_transect_line[transect_start_point:]
            x1 = raw_transect_line[transect_start_point][0]
            y1 = raw_transect_line[transect_start_point][1]
            x2 = raw_transect_line[-1][0]
            y2 = raw_transect_line[-1][1]
            accept = True
        except:
            COUNTERtransectsfull=COUNTERtransectsfull-1
            COUNTERtransectssmall=COUNTERtransectssmall-1
            accept = False



        if accept:
            w.line([[ [x1, y1], [x2, y2] ]])
            
            transectsFull['Transect '+str(COUNTERtransectsfull)] = np.array([ [x1, y1], [x2, y2] ])
            
            if slope < 0:
                angle = -1*(np.arctan(slope) - (np.pi/2))
            else:    
                angle = (np.pi/2) - np.arctan(slope)
            rotation = R.from_euler('z', [angle])
            transectsRotations.append(rotation)
            transectsSlopes.append(pslope)
            
            
            
            if COUNTERtransectssmall == 70:
                COUNTERtransectssmall=0
                COUNTERtransectssmall2=COUNTERtransectssmall2+1
                transects['Transect '+str(COUNTERtransectssmall2)] = np.array([ [x1, y1], [x2, y2] ])
            w.record(str(n))
        
        

    w.close()
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transects.prj")
    
    # convert transect shapefile into geoJSON
    shapefileTOgeoJSON_lines(filepath+"\\transects.geojson", filepath+"\\transects.shp")
    
    return filepath+"\\transects.geojson", transects, transectsFull, transectsRotations, transectsSlopes
        
    
    
    



def NewTransectsBuilder_StartParallel2(sitename, base_folder, output):
    
    filepath = base_folder+'data\\'+sitename
    
    
    
    
    # make shorelines into shapefile
    shores_geo_filename = filepath+"\\"+sitename+"_output.geojson"
    shores_shp_filename = filepath+"\\"+sitename+"_output.shp"
    geojsonTOshp_lines(shores_geo_filename, shores_shp_filename)
    
    
    
    # get image
    
    array_gray = copy.copy(output['imClassifs'][0])
    
    
    # print(c.shape)
    # print(array_gray.shape)
    # print(array_gray)
    # plt.figure(figsize=(30,30))
    # plt.imshow(array_gray)
    # plt.show()
    # plt.close('all')
    
    
    
    # np.savetxt(r"C:\Users\RDCHLNRO\Desktop\grey.txt", array_gray)######################################################################
    
    
    top = array_gray[1,:]
    right = array_gray[:,-2]
    bottom = array_gray[-2,:]
    left = array_gray[:,1]
    
    top2 = len(np.where(top==0.4284314)[0])/len(top)
    left2 = len(np.where(left==0.4284314)[0])/len(left)
    bottom2 = len(np.where(bottom==0.4284314)[0])/len(bottom)
    right2 = len(np.where(right==0.4284314)[0])/len(right)
    
    test = [top2, left2, bottom2, right2]
    
    open_water_position = np.where(test == np.nanmax(test))[0][0]
    
    

    
    # open a rster to get extent for transects
    try:
        rasters = glob.glob(filepath+"\\S2\\10m\\*_"+sitename+"_10m.tif")
        src = gdal.Open(rasters[0])
    except:
        rasters = glob.glob(filepath+"\\L8\\pan\\*.tif")
        src = gdal.Open(rasters[0])
    
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    src = None
    
    tryagain = 1
    while lrx-ulx < 1000:
        try:
            src = gdal.Open(rasters[tryagain])
            ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
            lrx = ulx + (src.RasterXSize * xres)
            lry = uly + (src.RasterYSize * yres)
            src = None
            tryagain = tryagain+1
        except:
            break
    
    
    
    
    # make reference shoreline into shapefile
    shoreRef_geo_Filename = filepath+"\\"+sitename+"_reference_shoreline.geojson"
    shoreRef_shp_filename = filepath+"\\"+sitename+"_reference_shoreline.shp"
    geojsonTOshp_lines(shoreRef_geo_Filename, shoreRef_shp_filename)
    
    # get reference shapefile points
    sf = shapefile.Reader(shoreRef_shp_filename)
    shapes = sf.shapes()
    points = shapes[0].points
    # less_points = points[1::70]
    sf.close()
    
    
    
    
    # add prpidicular points to new transect shapefile
    w = shapefile.Writer(filepath+"\\transects.shp")
    w.field('name', 'C')
    
    
    transects = dict([])
    transectsFull = dict([])
    transectsRotations = []
    transectsSlopes = []
    COUNTERtransectsfull=0
    COUNTERtransectssmall=0
    COUNTERtransectssmall2=0
    for n in range(5,len(points)-5):
        try:
            COUNTERtransectsfull=COUNTERtransectsfull+1
            COUNTERtransectssmall=COUNTERtransectssmall+1
            
            # for the chosen point along the reference shoreline find a slope for best fit line
            deriv_points = np.array(points[n-2:n+2])
            x = deriv_points[:,0]
            y = deriv_points[:,1]
            X = np.array(points)[n,0]
            Y = np.array(points)[n,1]
            slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
            
            # find a slope perpindicular to shorelne slope for transect
            pslope = -1*(1/slope)
            
            # start transect with first value
            XX = [X]
            YY = [Y + pslope*(XX[-1]-X)]
            
            # grow transect in forward X direction
            waterPointsForward = 0
            PointsForward = 0
            i = 0
            # greyClassification2 = [0]######################################################################
            while XX[-1]<lrx-3 and YY[-1]<uly-3 and YY[-1]>lry-3:
                i=i+1
                XX.append(XX[-1]+1)
                YY.append(Y + pslope*(XX[-1]-X))
                greyClassification = getArrayValueScaled(XX[-1], YY[-1], array_gray, lry,uly,ulx,lrx)
                # greyClassification2.append(greyClassification)######################################################################
                if ((greyClassification==2) or (greyClassification==3)):
                    waterPointsForward = waterPointsForward+1
                PointsForward = PointsForward+1
            # grow transect in backward X direction
            waterPointsReverse = 0
            PointsReverse = 0
            i = 0
            while XX[0]>ulx-3 and YY[0]<uly-3 and YY[0]>lry-3:
                i=i+1
                XX.insert(0, XX[0]-1)
                YY.insert(0, Y + pslope*(XX[0]-X))
                greyClassification = getArrayValueScaled(XX[0], YY[0], array_gray, lry,uly,ulx,lrx)
                # greyClassification2.insert(0, greyClassification)######################################################################
                if ((greyClassification==2) or (greyClassification==3)):
                    waterPointsReverse = waterPointsReverse+1
                PointsReverse = PointsReverse+1
                

            if (waterPointsReverse/PointsReverse) > (waterPointsForward/PointsForward):
                # cccc = 'black'
                XX = np.flip(np.array(XX))
                YY = np.flip(np.array(YY))
                

            
            raw_transect_line = []
            n=-1
            for all in YY:
                n=n+1
                raw_transect_line.append((XX[n],YY[n]))
            
            
            distances = distance(raw_transect_line, (X,Y)) 
            raw_transect_line = np.array(raw_transect_line)
            transect_start_point = np.where(distances==nearest(distances,100))[0][0]
            transect_processed = raw_transect_line[transect_start_point:]
            x1 = raw_transect_line[transect_start_point][0]
            y1 = raw_transect_line[transect_start_point][1]
            x2 = raw_transect_line[-1][0]
            y2 = raw_transect_line[-1][1]
            accept = True
        except:
            COUNTERtransectsfull=COUNTERtransectsfull-1
            COUNTERtransectssmall=COUNTERtransectssmall-1
            accept = False



        if accept:
            w.line([[ [x1, y1], [x2, y2] ]])
            transectsFull['Transect '+str(COUNTERtransectsfull)] = np.array([ [x1, y1], [x2, y2] ])
            
            
            if slope < 0:
                angle = -1*(np.arctan(slope) - (np.pi/2))
            else:    
                angle = (np.pi/2) - np.arctan(slope)
            rotation = R.from_euler('z', [angle])
            transectsRotations.append(rotation)
            transectsSlopes.append(pslope)
            
            if COUNTERtransectssmall == 70:
                COUNTERtransectssmall=0
                COUNTERtransectssmall2=COUNTERtransectssmall2+1
                transects['Transect '+str(COUNTERtransectssmall2)] = np.array([ [x1, y1], [x2, y2] ])
                
                
                
            w.record(str(n))
        
        

    w.close()
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transects.prj")
    
    # convert transect shapefile into geoJSON
    shapefileTOgeoJSON_lines(filepath+"\\transects.geojson", filepath+"\\transects.shp")
    
    return filepath+"\\transects.geojson", transects, transectsFull, transectsRotations, transectsSlopes
            
    
    




def NewTransectsBuilder_StartParallel3(sitename, base_folder, output):
    
    filepath = base_folder+'data\\'+sitename
    
    
    
    
    # make shorelines into shapefile
    shores_geo_filename = filepath+"\\"+sitename+"_output.geojson"
    shores_shp_filename = filepath+"\\"+sitename+"_output.shp"
    geojsonTOshp_lines(shores_geo_filename, shores_shp_filename)
    
    
    
    # get image
    
    array_gray = copy.copy(output['imClassifs'][0])
    
    
    # print(c.shape)
    # print(array_gray.shape)
    # print(array_gray)
    # plt.figure(figsize=(30,30))
    # plt.imshow(array_gray)
    # plt.show()
    # plt.close('all')
    
    
    
    # np.savetxt(r"C:\Users\RDCHLNRO\Desktop\grey.txt", array_gray)######################################################################
    
    
    top = array_gray[1,:]
    right = array_gray[:,-2]
    bottom = array_gray[-2,:]
    left = array_gray[:,1]
    
    top2 = len(np.where(top==0.4284314)[0])/len(top)
    left2 = len(np.where(left==0.4284314)[0])/len(left)
    bottom2 = len(np.where(bottom==0.4284314)[0])/len(bottom)
    right2 = len(np.where(right==0.4284314)[0])/len(right)
    
    test = [top2, left2, bottom2, right2]
    
    open_water_position = np.where(test == np.nanmax(test))[0][0]
    
    

    
    # open a rster to get extent for transects
    try:
        rasters = glob.glob(filepath+"\\S2\\10m\\*_"+sitename+"_10m.tif")
        src = gdal.Open(rasters[0])
    except:
        rasters = glob.glob(filepath+"\\L8\\pan\\*.tif")
        src = gdal.Open(rasters[0])
    
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    src = None
    
    tryagain = 1
    while lrx-ulx < 1000:
        try:
            src = gdal.Open(rasters[tryagain])
            ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
            lrx = ulx + (src.RasterXSize * xres)
            lry = uly + (src.RasterYSize * yres)
            src = None
            tryagain = tryagain+1
        except:
            break
    
    
    
    
    # make reference shoreline into shapefile
    shoreRef_geo_Filename = filepath+"\\"+sitename+"_reference_shoreline.geojson"
    shoreRef_shp_filename = filepath+"\\"+sitename+"_reference_shoreline.shp"
    geojsonTOshp_lines(shoreRef_geo_Filename, shoreRef_shp_filename)
    
    # get reference shapefile points
    sf = shapefile.Reader(shoreRef_shp_filename)
    shapes = sf.shapes()
    points = shapes[0].points
    # less_points = points[1::70]
    sf.close()
    
    
    
    
    
    # add prpidicular points to new transect shapefile
    w = shapefile.Writer(filepath+"\\transects.shp")
    w.field('name', 'C')
    w2 = shapefile.Writer(filepath+"\\transectsroots.shp")
    w2.field('name', 'C')
    
    
    transects = dict([])
    transectsFull = dict([])
    transectsRotations = []
    transectsSlopes = []
    COUNTERtransectsfull=0
    COUNTERtransectssmall=0
    COUNTERtransectssmall2=0
    XXs = []
    YYs = []
    Xs = []
    Ys = []
    
    for n in range(5,len(points)-5):
        COUNTERtransectsfull=COUNTERtransectsfull+1
        COUNTERtransectssmall=COUNTERtransectssmall+1
        
        # for the chosen point along the reference shoreline find a slope for best fit line
        deriv_points = np.array(points[n-2:n+2])
        x = deriv_points[:,0]
        y = deriv_points[:,1]
        X = np.array(points)[n,0]
        Y = np.array(points)[n,1]
        slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
        if slope ==0:
            slope = .00000001
        if np.isnan(slope):
            continue
        
        
#            more = 0
#            while np.isnan(slope):
#                more = more+1
#                deriv_points = np.array(points[n-2-more:n+2+more])
#                x = deriv_points[:,0]
#                y = deriv_points[:,1]
#                X = np.array(points)[n,0]
#                Y = np.array(points)[n,1]
#                slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
#            deriv_points = np.array(points[n-2-more:n+2+more])
#            x = deriv_points[:,0]
#            y = deriv_points[:,1]
#            X = np.array(points)[n,0]
#            Y = np.array(points)[n,1]
#            slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
        
        # find a slope perpindicular to shorelne slope for transect
        pslope = -1*(1/slope)
        
        # start transect with first value
        XX = [X]
#            YY = [Y + pslope*(XX[-1]-X)]
        YY = [Y]
        
        
        # grow transect in forward X direction
        waterPointsForward = []
        while distancePP([XX[0],YY[0]], [XX[-1],YY[-1]]) < 700:
#                XX.append(XX[-1]+1)
#                YY.append(Y + pslope*(XX[-1]-X))
            XX.append(XX[-1] + 1*((1)/(np.sqrt(1+pslope**2))))
            YY.append(YY[-1] + 1*((pslope)/(np.sqrt(1+pslope**2))))
            greyClassification = getArrayValueScaled(XX[-1], YY[-1], array_gray, lry,uly,ulx,lrx)
            if ((greyClassification==2) or (greyClassification==3)):
                waterPointsForward.append(1)
            else:
                waterPointsForward.append(0)
            
            
            
            
            
        # grow transect in backward X direction
        waterPointsReverse = []
        while distancePP([XX[0],YY[0]], [XX[-1],YY[-1]]) < 1400:
#                XX.insert(0, XX[0]-1)
#                YY.insert(0, Y + pslope*(XX[0]-X))
            XX.insert(0,XX[0] - 1*((1)/(np.sqrt(1+pslope**2))))
            YY.insert(0,YY[0] - 1*((pslope)/(np.sqrt(1+pslope**2))))
            greyClassification = getArrayValueScaled(XX[0], YY[0], array_gray, lry,uly,ulx,lrx)
            if ((greyClassification==2) or (greyClassification==3)):
                waterPointsReverse.append(1)
            else:
                waterPointsReverse.append(0)
            
            
        # which direction reaches water first
        if np.nansum(waterPointsForward[0:200])>np.nansum(waterPointsReverse[0:200]):
            REVERSE = False
        if np.nansum(waterPointsForward[0:200])<np.nansum(waterPointsReverse[0:200]):
            REVERSE = True
        else:
            REVERSE=False
        if REVERSE:
            XX = np.flip(np.array(XX))
            YY = np.flip(np.array(YY))
            
            
        XX2 = XX[500:]
        YY2 = YY[500:]

            
            
        XXs.append(XX2)
        YYs.append(YY2)
        Xs.append(X)
        Ys.append(Y)

          
        
        

        
        
    
    if len(XXs)==0:
        print('There is something wrong and no transects were created.')
    
    
    c = -1
    COUNTERtransectssmall = -1
    COUNTERtransectssmall2 = 0
    for n in XXs:
        COUNTERtransectssmall = COUNTERtransectssmall+1
        c=c+1
        X = Xs[c]
        XX = XXs[c]
        Y = Ys[c]
        YY = YYs[c]

        x1 = XX[0]
        y1 = YY[0]
        x2 = XX[-1]
        y2 = YY[-1]
        
        
        if np.all(n!=XXs[0]) and distancePP([x1,y1], [XXs[c-1][0],YYs[c-1][1]]) > 20:
            continue
        
        w.line([[ [x1, y1], [x2, y2] ]])
        
        transectsFull['Transect '+str(COUNTERtransectsfull)] = np.array([ [x1, y1], [x2, y2] ])
        
        if slope < 0:
            angle = -1*(np.arctan(slope) - (np.pi/2))
        else:    
            angle = (np.pi/2) - np.arctan(slope)
        rotation = R.from_euler('z', [angle])
        transectsRotations.append(rotation)
        transectsSlopes.append(pslope)
        
        
        
        if COUNTERtransectssmall == 70:
            COUNTERtransectssmall=0
            COUNTERtransectssmall2=COUNTERtransectssmall2+1
            transects['Transect '+str(COUNTERtransectssmall2)] = np.array([ [x1, y1], [x2, y2] ])
            
            
            
        w.record(str(X)+','+str(Y))
        w2.point(x1, y1)
        w2.record(str(X)+','+str(Y))
        
        
        
        
        

    w.close()
    w2.close()
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transects.prj")
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transectsroots.prj")
    
    
    # convert transect shapefile into geoJSON
    shapefileTOgeoJSON_lines(filepath+"\\transects.geojson", filepath+"\\transects.shp")
    
    
    
    
    
    
    return filepath+"\\transects.geojson", transects, transectsFull, transectsRotations, transectsSlopes
            
    










def NewTransectsBuilder_StartParallel4(sitename, base_folder, output):
    
    filepath = base_folder+'\\data\\'+sitename
    
        
    
    # get image
    
    array_gray = copy.copy(output['imClassifs'][0])
    
#    
#    plt.figure(figsize=(10,10))
#    plt.imshow(array_gray)
#    plt.colorbar()
#    plt.show()
#    plt.close('all')
#    
    

    
    

    
    # open a rster to get extent for transects
    try:
        rasters = glob.glob(filepath+"\\S2\\10m\\*_"+sitename+"_10m.tif")
        src = gdal.Open(rasters[0])
    except:
        rasters = glob.glob(filepath+"\\L8\\pan\\*.tif")
        src = gdal.Open(rasters[0])
    
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    src = None
    
    tryagain = 1
    while lrx-ulx < 1000:
        try:
            src = gdal.Open(rasters[tryagain])
            ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
            lrx = ulx + (src.RasterXSize * xres)
            lry = uly + (src.RasterYSize * yres)
            src = None
            tryagain = tryagain+1
        except:
            break
    
    
    
    
    shoreRef_shp_filename = filepath+"\\"+sitename+"_reference_shoreline.shp"
    # get reference shapefile points
    sf = shapefile.Reader(shoreRef_shp_filename)
    shapes = sf.shapes()
    points = shapes[0].points
    # less_points = points[1::70]
    sf.close()
    
    
    
    
    
    # add prpidicular points to new transect shapefile
    w = shapefile.Writer(filepath+"\\transects.shp")
    w.field('name', 'C')
    w.field('num', 'C')
    
    w2 = shapefile.Writer(filepath+"\\transectsroots.shp")
    w2.field('name', 'C')
    w2.field('num', 'C')
    
    w3 = shapefile.Writer(filepath+"\\transectsSmall.shp")
    w3.field('name', 'C')
    w3.field('num', 'C')
    
    slope, intercept, r_value, p_value, std_err = st.linregress(np.array(points)[:,0],np.array(points)[:,1])
    if 1 >= slope <= -1 :
        VERT = True
    else:
        VERT = False
    
    
    
    transects = dict([])
    transectsFull = dict([])
    transectsRotations = []
    transectsSlopes = []
    COUNTERtransectsfull=0
    COUNTERtransectssmall=0
    COUNTERtransectssmall2=0
    XXs = []
    YYs = []
    XXXs = []
    YYYs = []
    Xs = []
    Ys = []
    REVERSES = []
    pslopes = []
    skips = 0
    for n in range(5,len(points)-5):
        COUNTERtransectsfull=COUNTERtransectsfull+1
        COUNTERtransectssmall=COUNTERtransectssmall+1
        
        # for the chosen point along the reference shoreline find a slope for best fit line
        deriv_points = np.array(points[n-2:n+2])
        x = deriv_points[:,0]
        y = deriv_points[:,1]
        X = np.array(points)[n,0]
        Y = np.array(points)[n,1]
        slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
        if slope ==0:
            slope = .00000001
        if np.isnan(slope):
            skips=skips+1
            continue
        
        

        
        # find a slope perpindicular to shorelne slope for transect
        pslope = -1*(1/slope)
        
        # start transect with first value
        XX = [X]
#            YY = [Y + pslope*(XX[-1]-X)]
        YY = [Y]
        
        
        # grow transect in forward X direction
        waterPointsForward = []
        while distancePP([XX[0],YY[0]], [XX[-1],YY[-1]]) < 700:
#                XX.append(XX[-1]+1)
#                YY.append(Y + pslope*(XX[-1]-X))
            XX.append(XX[-1] + 1*((1)/(np.sqrt(1+pslope**2))))
            YY.append(YY[-1] + 1*((pslope)/(np.sqrt(1+pslope**2))))
            greyClassification = getArrayValueScaled(XX[-1], YY[-1], array_gray, lry,uly,ulx,lrx)
            if ((greyClassification==2) or (greyClassification==3)):
                waterPointsForward.append(1)
            else:
                waterPointsForward.append(0)
            
            
            
            
            
        # grow transect in backward X direction
        waterPointsReverse = []
        while distancePP([XX[0],YY[0]], [XX[-1],YY[-1]]) < 1400:
#                XX.insert(0, XX[0]-1)
#                YY.insert(0, Y + pslope*(XX[0]-X))
            XX.insert(0,XX[0] - 1*((1)/(np.sqrt(1+pslope**2))))
            YY.insert(0,YY[0] - 1*((pslope)/(np.sqrt(1+pslope**2))))
            greyClassification = getArrayValueScaled(XX[0], YY[0], array_gray, lry,uly,ulx,lrx)
            if ((greyClassification==2) or (greyClassification==3)):
                waterPointsReverse.append(1)
            else:
                waterPointsReverse.append(0)
            
            
        # which direction reaches water first
        if np.nansum(waterPointsForward[0:200])>np.nansum(waterPointsReverse[0:200]):
            REVERSE = False
        if np.nansum(waterPointsForward[0:200])<np.nansum(waterPointsReverse[0:200]):
            REVERSE = True
        else:
            REVERSE=False
            
        if not VERT:
            if pslope < 0:
                XX = np.flip(np.array(XX))
                YY = np.flip(np.array(YY))
                REVERSE = not REVERSE
            
        REVERSES.append(REVERSE)
        XXXs.append(XX)
        YYYs.append(YY)
        Xs.append(X)
        Ys.append(Y)
        pslopes.append(pslope)
            
            
    if REVERSES.count(False) > REVERSES.count(True):
        REVERSE_CONSENSUS = False
    else:
        REVERSE_CONSENSUS = True
        
        
        
        
#    plt.figure(figsize=(20,20))
#    nn = -1
#    for n in range(5,len(points)-5):
#        nn=nn+1
#        XX = XXXs[nn]
#        YY = YYYs[nn]
#        XXa=XX[0]
#        XXb=XX[-1]
#        YYa=YY[0]
#        YYb=YY[-1]
#        plt.scatter(XXa,YYa,s=1,c='red')
#        plt.scatter(XXb,YYb,s=1,c='black')
#    plt.plot(np.array(points)[:,0],np.array(points)[:,1],c='red')
#    plt.gca().set_aspect('equal', adjustable='box')
#    plt.show()
#    plt.close('all')
        
        
        
        
        
    nn = -1
    for n in range(5,len(points)-5-skips):
        nn=nn+1
        XX = XXXs[nn]
        YY = YYYs[nn]
            
        if REVERSE_CONSENSUS:
            XX = np.flip(np.array(XX))
            YY = np.flip(np.array(YY))
            
        XX2 = XX[500:]
        YY2 = YY[500:]
            
        XXs.append(XX2)
        YYs.append(YY2)




        
#    plt.figure(figsize=(20,20))
#    nn = -1
#    for n in range(5,len(points)-5):
#        nn=nn+1
#        XX = XXs[nn]
#        YY = YYs[nn]
#        plt.scatter(XX,YY,s=1)
#    plt.plot(np.array(points)[:,0],np.array(points)[:,1],c='red')
#    plt.gca().set_aspect('equal', adjustable='box')
#    plt.show()
#    plt.close('all')        
        
        
        
        
    
    if len(XXs)==0:
        print('There is something wrong and no transects were created.')
    
    
    c = -1
    COUNTERtransectssmall = -1
    COUNTERtransectssmall2 = 0
    for n in XXs:
        COUNTERtransectssmall = COUNTERtransectssmall+1
        c=c+1
        X = Xs[c]
        XX = XXs[c]
        Y = Ys[c]
        YY = YYs[c]

        x1 = XX[0]
        y1 = YY[0]
        x2 = XX[-1]
        y2 = YY[-1]
        
        
        if np.all(n!=XXs[0]) and distancePP([x1,y1], [XXs[c-1][0],YYs[c-1][1]]) > 20:
            continue
        
        w.line([[ [x1, y1], [x2, y2] ]])
        
        transectsFull['Transect '+str(COUNTERtransectsfull)] = np.array([ [x1, y1], [x2, y2] ])
        
        if slope < 0:
            angle = -1*(np.arctan(slope) - (np.pi/2))
        else:    
            angle = (np.pi/2) - np.arctan(slope)
        rotation = R.from_euler('z', [angle])
        transectsRotations.append(rotation)
        transectsSlopes.append(pslope)
        
        
        
        if COUNTERtransectssmall == 70:
            COUNTERtransectssmall=0
            COUNTERtransectssmall2=COUNTERtransectssmall2+1
            transects['Transect '+str(COUNTERtransectssmall2)] = np.array([ [x1, y1], [x2, y2] ])
            w3.line([[ [x1, y1], [x2, y2] ]])
            n2='Transect '+str(COUNTERtransectssmall2)
            c2=str(X)+','+str(Y)
            r = [c2, n2]
            w3.record(*r)
        
        
        n2='Transect '+str(COUNTERtransectsfull)
        c2=str(X)+','+str(Y)
        r = [c2, n2]
        w.record(*r)
        w2.point(x1, y1)
        n2='Transect '+str(COUNTERtransectsfull)
        c2=str(X)+','+str(Y)
        r = [c2, n2]
        w2.record(*r)
        
        
        
        
        

    w.close()
    w2.close()
    w3.close()
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transects.prj")
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transectsroots.prj")
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transectsSmall.prj")
    
    # convert transect shapefile into geoJSON
    shapefileTOgeoJSON_lines(filepath+"\\transects.geojson", filepath+"\\transects.shp")
    
    # convert transectsmall shapefile into geoJSON
    shapefileTOgeoJSON_lines(filepath+"\\transectsSmall.geojson", filepath+"\\transectsSmall.shp")
    
    
    
    
    return filepath+"\\transects.geojson", filepath+"\\transectsSmall.geojson", transects, transectsFull, transectsRotations, transectsSlopes
            
    












def NewTransectsBuilder_StartParallel5(sitename, base_folder, output):
    
        
    filepath = base_folder+'data\\'+sitename
    
    
    
    
    
    
    
    # get image
    
    array_gray = copy.copy(output['imClassifs'][0])
    
#    
#    plt.figure(figsize=(10,10))
#    plt.imshow(array_gray)
#    plt.colorbar()
#    plt.show()
#    plt.close('all')
#    
    

    
    

    
    # open a rster to get extent for transects
    try:
        rasters = glob.glob(filepath+"\\S2\\10m\\*_"+sitename+"_10m.tif")
        src = gdal.Open(rasters[0])
    except:
        rasters = glob.glob(filepath+"\\L8\\pan\\*.tif")
        src = gdal.Open(rasters[0])
    
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    src = None
    
    tryagain = 1
    while lrx-ulx < 1000:
        try:
            src = gdal.Open(rasters[tryagain])
            ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
            lrx = ulx + (src.RasterXSize * xres)
            lry = uly + (src.RasterYSize * yres)
            src = None
            tryagain = tryagain+1
        except:
            break
    
    
    
    
    shoreRef_shp_filename = filepath+"\\"+sitename+"_reference_shoreline.shp"
    # get reference shapefile points
    sf = shapefile.Reader(shoreRef_shp_filename)
    shapes = sf.shapes()
    points = shapes[0].points
    # less_points = points[1::70]
    sf.close()
    
    
    
    
    
    # add prpidicular points to new transect shapefile
    w = shapefile.Writer(filepath+"\\transects.shp")
    w.field('name', 'C')
    w2 = shapefile.Writer(filepath+"\\transectsroots.shp")
    w2.field('name', 'C')
    
    
    
    slope, intercept, r_value, p_value, std_err = st.linregress(np.array(points)[:,0],np.array(points)[:,1])
    if 1 >= slope <= -1 :
        VERT = True
    else:
        VERT = False
    
    
    
    transects = dict([])
    transectsFull = dict([])
    transectsRotations = []
    transectsSlopes = []
    COUNTERtransectsfull=0
    COUNTERtransectssmall=0
    COUNTERtransectssmall2=0
    XXs = []
    YYs = []
    XXXs = []
    YYYs = []
    Xs = []
    Ys = []
    REVERSES = []
    pslopes = []
    for n in range(5,len(points)-5):
        COUNTERtransectsfull=COUNTERtransectsfull+1
        COUNTERtransectssmall=COUNTERtransectssmall+1
        
        # for the chosen point along the reference shoreline find a slope for best fit line
        deriv_points = np.array(points[n-2:n+2])
        x = deriv_points[:,0]
        y = deriv_points[:,1]
        X = np.array(points)[n,0]
        Y = np.array(points)[n,1]
        slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
        if slope ==0:
            slope = .00000001
        if np.isnan(slope):
            continue
        
        

        
        # find a slope perpindicular to shorelne slope for transect
        pslope = -1*(1/slope)
        
        # start transect with first value
        XX = [X]
#            YY = [Y + pslope*(XX[-1]-X)]
        YY = [Y]
        
        
        # grow transect in forward X direction
        waterPointsForward = []
        while distancePP([XX[0],YY[0]], [XX[-1],YY[-1]]) < 700:
#                XX.append(XX[-1]+1)
#                YY.append(Y + pslope*(XX[-1]-X))
            XX.append(XX[-1] + 1*((1)/(np.sqrt(1+pslope**2))))
            YY.append(YY[-1] + 1*((pslope)/(np.sqrt(1+pslope**2))))
            greyClassification = getArrayValueScaled(XX[-1], YY[-1], array_gray, lry,uly,ulx,lrx)
            if ((greyClassification==2) or (greyClassification==3)):
                waterPointsForward.append(1)
            else:
                waterPointsForward.append(0)
            
            
            
            
            
        # grow transect in backward X direction
        waterPointsReverse = []
        while distancePP([XX[0],YY[0]], [XX[-1],YY[-1]]) < 1400:
#                XX.insert(0, XX[0]-1)
#                YY.insert(0, Y + pslope*(XX[0]-X))
            XX.insert(0,XX[0] - 1*((1)/(np.sqrt(1+pslope**2))))
            YY.insert(0,YY[0] - 1*((pslope)/(np.sqrt(1+pslope**2))))
            greyClassification = getArrayValueScaled(XX[0], YY[0], array_gray, lry,uly,ulx,lrx)
            if ((greyClassification==2) or (greyClassification==3)):
                waterPointsReverse.append(1)
            else:
                waterPointsReverse.append(0)
            
            
        # which direction reaches water first
        if np.nansum(waterPointsForward[0:200])>np.nansum(waterPointsReverse[0:200]):
            REVERSE = False
        if np.nansum(waterPointsForward[0:200])<np.nansum(waterPointsReverse[0:200]):
            REVERSE = True
        else:
            REVERSE=False
            
        if not VERT:
            if pslope < 0:
                XX = np.flip(np.array(XX))
                YY = np.flip(np.array(YY))
                REVERSE = not REVERSE
            
        REVERSES.append(REVERSE)
        XXXs.append(XX)
        YYYs.append(YY)
        Xs.append(X)
        Ys.append(Y)
        pslopes.append(pslope)
            
            
    if REVERSES.count(False) > REVERSES.count(True):
        REVERSE_CONSENSUS = False
    else:
        REVERSE_CONSENSUS = True
        
        
        
        
    plt.figure(figsize=(20,20))
    nn = -1
    for n in range(5,len(points)-5):
        nn=nn+1
        XX = XXXs[nn]
        YY = YYYs[nn]
        XXa=XX[0]
        XXb=XX[-1]
        YYa=YY[0]
        YYb=YY[-1]
        plt.scatter(XXa,YYa,s=1,c='red')
        plt.scatter(XXb,YYb,s=1,c='black')
    plt.plot(np.array(points)[:,0],np.array(points)[:,1],c='red')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
    plt.close('all')
        
        
        
        
        
    nn = -1
    for n in range(5,len(points)-5):
        nn=nn+1
        XX = XXXs[nn]
        YY = YYYs[nn]
            
        if REVERSE_CONSENSUS:
            XX = np.flip(np.array(XX))
            YY = np.flip(np.array(YY))
            
        XX2 = XX[500:]
        YY2 = YY[500:]
            
        XXs.append(XX2)
        YYs.append(YY2)




        
#    plt.figure(figsize=(20,20))
#    nn = -1
#    for n in range(5,len(points)-5):
#        nn=nn+1
#        XX = XXs[nn]
#        YY = YYs[nn]
#        plt.scatter(XX,YY,s=1)
#    plt.plot(np.array(points)[:,0],np.array(points)[:,1],c='red')
#    plt.gca().set_aspect('equal', adjustable='box')
#    plt.show()
#    plt.close('all')        
        
        
        
        
    
    if len(XXs)==0:
        print('There is something wrong and no transects were created.')
    
    
    c = -1
    COUNTERtransectssmall = -1
    COUNTERtransectssmall2 = 0
    for n in XXs:
        COUNTERtransectssmall = COUNTERtransectssmall+1
        c=c+1
        X = Xs[c]
        XX = XXs[c]
        Y = Ys[c]
        YY = YYs[c]

        x1 = XX[0]
        y1 = YY[0]
        x2 = XX[-1]
        y2 = YY[-1]
        
        
        if np.all(n!=XXs[0]) and distancePP([x1,y1], [XXs[c-1][0],YYs[c-1][1]]) > 20:
            continue
        
        w.line([[ [x1, y1], [x2, y2] ]])
        
        transectsFull['Transect '+str(COUNTERtransectsfull)] = np.array([ [x1, y1], [x2, y2] ])
        
        if slope < 0:
            angle = -1*(np.arctan(slope) - (np.pi/2))
        else:    
            angle = (np.pi/2) - np.arctan(slope)
        rotation = R.from_euler('z', [angle])
        transectsRotations.append(rotation)
        transectsSlopes.append(pslope)
        
        
        
        if COUNTERtransectssmall == 70:
            COUNTERtransectssmall=0
            COUNTERtransectssmall2=COUNTERtransectssmall2+1
            transects['Transect '+str(COUNTERtransectssmall2)] = np.array([ [x1, y1], [x2, y2] ])
            
            
            
        w.record(str(X)+','+str(Y))
        w2.point(x1, y1)
        w2.record(str(X)+','+str(Y))
        
        
        
        
        

    w.close()
    w2.close()
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transects.prj")
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\transectsroots.prj")
    
    
    # convert transect shapefile into geoJSON
    shapefileTOgeoJSON_lines(filepath+"\\transects.geojson", filepath+"\\transects.shp")
    
    
    
    
    
    
    return filepath+"\\transects.geojson", transects, transectsFull, transectsRotations, transectsSlopes
            
    



    return filepath+"\\transects.geojson", transects, transectsFull, transectsRotations, transectsSlopes
















    
    
def BeachSlopes(sitename, base_folder,dates_sat,cross_distance_tidally_corrected,cross_distance,TC,output):
    
    
    slopes = []
    p_values = []
    
    if TC:
        
        dates = dates_sat
        dates2 = [0]
        d = dates[0]
        for date in dates[1:]:
            dates2.append((date-d).total_seconds())
        
        dates2=np.array(dates2)
        
        
        for i in cross_distance_tidally_corrected:
            if i != 'dates':
                ts = cross_distance_tidally_corrected[i]
                
                tsGood = ts[~np.isnan(ts)]
                dates2Good = dates2[~np.isnan(ts)]
                
                slope, intercept, r_value, p_value, std_err = st.linregress(dates2Good, tsGood)
                slopes.append(slope)
                p_values.append(p_value)
            
    else:
        
        dates = output['dates']
        dates2 = [0]
        d = dates[0]
        for date in dates[1:]:
            dates2.append((date-d).total_seconds())
        
        dates2=np.array(dates2)
        
        
        for i in cross_distance:
            if i != 'dates':

                ts = cross_distance[i]
                
                tsGood = ts[~np.isnan(ts)]
                dates2Good = dates2[~np.isnan(ts)]
                
                slope, intercept, r_value, p_value, std_err = st.linregress(dates2Good, tsGood)
                slopes.append(slope)
                p_values.append(p_value)

    return slopes,p_values

def Average_Shoreline2(sitename, base_folder, refPoints,shoreOutputPoints):
    
    filepath = base_folder+'data\\'+sitename
    
    

    shoreOutputPoints = np.array(shoreOutputPoints)
    
    

    
    
    # add prpidicular points to new transect shapefile
    avg_shoreline_points = []
    transect_rotations = []
    
    for n in range(5,len(refPoints)-5):
        # for the chosen point along the reference shoreline find a slope for best fit line
        deriv_points = np.array(refPoints[n-2:n+2])
        x = deriv_points[:,0]
        y = deriv_points[:,1]
        X = np.array(refPoints)[n,0]
        Y = np.array(refPoints)[n,1]
        XY = np.array([X,Y,0])
        slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
        if slope < 0:
            angle = -1*(np.arctan(slope) - (np.pi/2))
        else:    
            angle = (np.pi/2) - np.arctan(slope)
        rotation = R.from_euler('z', [angle])
        transect_rotations.append(rotation)
        
        transect_intersections = []
        for shoreline2 in shoreOutputPoints:
            shoreline = np.append(shoreline2,np.zeros([len(shoreline2),1]),1)
            rotatedShoreline = rotation.apply(shoreline)
            rotatedXY = rotation.apply(XY)
            rotatedY = rotatedXY[0][1]
            rotatedShoreline[:,1] = rotatedShoreline[:,1] - rotatedY
            
            rotatedCoordinate = rotatedShoreline[np.abs(rotatedShoreline[:,1])==np.nanmin(np.abs(rotatedShoreline[:,1])),:]
                    
            rotatedCoordinate[0][1] = rotatedCoordinate[0][1] + rotatedY
            coordinate = rotation.inv().apply(rotatedCoordinate)
            transect_intersections.append(coordinate[0])
            
        transect_intersections = np.array(transect_intersections)
        avg_shoreline_points.append( [ np.nanmean(transect_intersections[:,0]), np.nanmean(transect_intersections[:,1]) ] )    
        
        
    cc_save.Save_Average_Shoreline(filepath,sitename,avg_shoreline_points)
    
    
    
    
    return avg_shoreline_points








        
def WhiteWaterByPixel2(sitename, base_folder, output):        
    
    whitewaterArrays = []
    classArrays = []
    classificationArraysExtents = []
    SRSs = []
    SGTs = []
    
    # with open(os.path.join(base_folder+ "data\\"+sitename, sitename + '_output' + '.pkl'), 'rb') as f:
    #     output = pickle.load(f) 
    
    
    shore_feature = -1
    for image in output['filename']:
        
        shore_feature = shore_feature+1
        
        # get classified image    
        array_whitewater = output['imClassifs'][shore_feature]
        
        # find image
        aa = glob.glob(base_folder+ "data\\"+sitename+"\\*\\*\\"+image)
        
        # open image
        src = gdal.Open(aa[0])
        SRS = src.GetProjection()
        
        # get extent
        ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
        lrx = ulx + (src.RasterXSize * xres)
        lry = uly + (src.RasterYSize * yres)
        extent = [ulx, lrx, lry, uly]
        
        # save extent
        classificationArraysExtents.append(extent)
               
        # save classification raster
        classArrays.append(copy.copy(array_whitewater))
        
        # turn classified image into whitewater only
        array_whitewater[array_whitewater!=2]=np.nan
        array_whitewater[array_whitewater==2]=1
        
        # save WW raster
        whitewaterArrays.append(array_whitewater)
        
        
        # get projection info
        SGT = list(src.GetGeoTransform())
        aa = np.shape(array_whitewater)
        bb = np.shape(np.array(src.GetRasterBand(1).ReadAsArray()))
        SGT[1] = SGT[1]* (bb[1]/aa[1])
        SGT[5] = SGT[5]*(bb[0]/aa[0])
        SGT = tuple(SGT)
        SRS = src.GetProjection()
        
        # save arrays as rasters
        cc_save.Save_Classification_Rasters(sitename, base_folder, image, whitewaterArrays[-1], classArrays[-1], SRS, SGT)
        
        SRSs.append(SRS)
        SGTs.append(SGT)
        
        
    output['SRS'] = SRSs
    output['SGT'] = SGTs
    output['Extents'] = classificationArraysExtents
    output['WhiteWaterMask'] = whitewaterArrays
    output['ClassificationMask'] = classArrays
        
        
    return output



def GetProjectionInfo(sitename, base_folder, output):        
    
    classificationArraysExtents = []
    SRSs = []
    SGTs = []
    
    # with open(os.path.join(base_folder+ "data\\"+sitename, sitename + '_output' + '.pkl'), 'rb') as f:
    #     output = pickle.load(f) 
    
    
    shore_feature = -1
    for image in output['filename']:
        
        shore_feature = shore_feature+1
        
        # get classified image    
        array_whitewater = output['imClassifs'][shore_feature]
        
        # find image
        aa = glob.glob(base_folder+ "data\\"+sitename+"\\*\\*\\"+image)
        
        # open image
        src = gdal.Open(aa[0])
        SRS = src.GetProjection()
        
        # get extent
        ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
        lrx = ulx + (src.RasterXSize * xres)
        lry = uly + (src.RasterYSize * yres)
        extent = [ulx, lrx, lry, uly]
        
        # save extent
        classificationArraysExtents.append(extent)
               
        
        
        # get projection info
        SGT = list(src.GetGeoTransform())
        aa = np.shape(array_whitewater)
        bb = np.shape(np.array(src.GetRasterBand(1).ReadAsArray()))
        SGT[1] = SGT[1]* (bb[1]/aa[1])
        SGT[5] = SGT[5]*(bb[0]/aa[0])
        SGT = tuple(SGT)
        SRS = src.GetProjection()
        
        # save arrays as rasters
        # cc_save.Save_Classification_Rasters(sitename, base_folder, image, whitewaterArrays[-1], classArrays[-1], SRS, SGT)
        
        SRSs.append(SRS)
        SGTs.append(SGT)
        
        
    output['SRS'] = SRSs
    output['SGT'] = SGTs
    output['Extents'] = classificationArraysExtents
        
        
    return output

        






def get_basemap(output):
    # classed_images = glob.glob(base_folder+ "data\\"+sitename+"\\jpg_files\\detection\\*.png")
    # dateWanted = classed_images[last_s2].split('\\')[-1][:-4]
    # imageWanted = base_folder+ "data\\"+sitename+"\\jpg_files\\preprocessed\\"+dateWanted+".png"
    # im = plt.imread(imageWanted)
    
    # # find image
    # c = imageExtracterS2(im)
    
    try:
        last_s2 = np.where(np.array(output['satname'])=='S2')[0][-1]
    except:
        last_s2 = np.where(np.array(output['satname'])=='L8')[0][-1]
    
    im = output['imRGBs'][last_s2]
    extent = output['Extents'][last_s2]
    
    return im,extent















def WarpArrayAsRaster(Array1,SRS1,SGT1,Extent1,Array2,SRS2,SGT2,Extent2):
    Raster1 = array2raster_InMemory(Array1,SRS1,SGT1)
    
    Raster2 = array2raster_InMemory(Array2,SRS2,SGT2)
    
    outputbound = (Extent2[0],Extent2[2],Extent2[1],Extent2[3])
    Array2Shape = Array2.shape
    
    Raster3 = gdal.Warp('', Raster1, format='MEM',targetAlignedPixels = False, width= Array2Shape[1], height= Array2Shape[0], outputBounds= outputbound)
    outputArray = np.array(Raster3.GetRasterBand(1).ReadAsArray())
    
    return outputArray



def WhiteWaterExtentSumAll2(sitename, base_folder, transectsFull, transectsSlopes, output):
    filepath = base_folder+'data\\'+sitename
    
    try:
        last_s2 = np.where(np.array(output['satname'])=='S2')[0][-1]
    except:
        last_s2 = np.where(np.array(output['satname'])=='L8')[0][-1]
    lasts2SRS = output['SRS'][last_s2]
    lasts2SGT = output['SGT'][last_s2]
    lasts2Extent = output['Extents'][last_s2]
    lasts2Shape = output['WhiteWaterMask'][last_s2].shape
    
    WWsum = np.zeros(lasts2Shape)
    for n in range(0,len(output['WhiteWaterMask'])):
        if output['WhiteWaterMask'][n].shape != lasts2Shape:
            WW_raster = WarpArrayAsRaster(output['WhiteWaterMask'][n],output['SRS'][n],output['SGT'][n],output['Extents'][n],output['WhiteWaterMask'][last_s2],output['SRS'][last_s2],output['SGT'][last_s2],output['Extents'][last_s2]) 
            
            # WW_raster[WW_raster!=0]=1
            # plt.imshow(WW_raster)
            # plt.title(np.nansum(WW_raster))
            # plt.show()
            # plt.close('all')
            # print(WW_raster)
            WWsum = WWsum + np.nan_to_num(WW_raster)
            
        elif output['WhiteWaterMask'][n].shape == lasts2Shape:
            WWsum = WWsum + np.nan_to_num(output['WhiteWaterMask'][n])
    
    
    
    # plt.figure(figsize=(30,30))
    # plt.imshow(WWsum)
    # plt.title(np.nansum(WWsum))
    # plt.show()
    # plt.close('all')

    
    
    ly= lasts2Extent[2]
    uy= lasts2Extent[3]
    lx= lasts2Extent[0]
    ux= lasts2Extent[1]
    
    
    
    # plt.figure(figsize=(30,30))
    # plt.imshow(WWsum, extent = lasts2Extent)
    # plt.scatter(lrx,lry,s=1000,c='black')
    # plt.scatter(ulx,uly,s=1000,c='red')
    # plt.title(np.nansum(WWsum))
    # plt.show()
    # plt.close('all')    
    
    
    
    n=-1
    transectt = []
    positions = []
    for key, transect in transectsFull.items():
        n=n+1
        
        pslope = transectsSlopes[n]
        X = transect[0][0]
        Y = transect[0][1]
        
        # generate transect values
        XX = np.linspace(X, transect[1][0], 1000)
        YY = Y + pslope*(XX-X)

        arrayValuesAlongTransect = []
        for nn in range(0, len(XX)):
            arrayValuesAlongTransect.append(getArrayValueScaled2(XX[nn], YY[nn], WWsum, ly,uy,lx,ux))
        arrayValuesAlongTransect = np.array(arrayValuesAlongTransect)
        try:
            arrayValuesAlongTransect[arrayValuesAlongTransect == 0] = np.nan
        except ValueError:
            arrayValuesAlongTransect = np.nan*np.ones(len(arrayValuesAlongTransect))
        try:
            a = np.where(arrayValuesAlongTransect>1)[0]
            position = a[-1]
            if position == 999:
                x = np.nan
                y = np.nan
            else:
                x = XX[position]
                y = YY[position]
            
        except:
            position = 999
            x = np.nan
            y = np.nan
        transectt.append([x,y])
        positions.append(position)
        
        
    # plt.figure(figsize=(30,30))
    # plt.imshow(WWsum, extent = lasts2Extent)
    # # plt.scatter(lrx,lry,s=1000,c='black')
    # # plt.scatter(ulx,uly,s=1000,c='red')
    # plt.scatter(np.array(transectt)[:,0],np.array(transectt)[:,1],s=100, c='green')
    # plt.title(np.nansum(WWsum))
    # plt.show()
    # plt.close('all')
        
    try:
        cc_save.MaxOuterWaveBreak(filepath,sitename,transectt)
    except:
        print('Couldnt Save MaxOuterWaveBreak as shapefile')
    return transectt
    
        
    
def getArrayValueScaled2(x, y, array, ly,uy,lx,ux):
    xpercent = (x-lx)/(ux-lx)
    ypercent = (y-ly)/(uy-ly)
    ypercent = 1-ypercent
    C = xpercent*array.shape[1]
    R = ypercent*array.shape[0]
    C = int(C)
    R = int(R)
    try:
        arrayValue = array[R,C]
    except:
        arrayValue = 0
    return arrayValue    
    
    
    
    


def array2raster_InMemory(array,SRS,SGT):
    driver = gdal.GetDriverByName("MEM")
    raster = driver.Create("myraster", array.shape[1], array.shape[0], 1, gdal.GDT_Float32)
    raster.GetRasterBand(1).SetNoDataValue(0)
    raster.SetGeoTransform(SGT)
    raster.SetProjection(SRS)
    raster.GetRasterBand(1).WriteArray(array)
    return raster




def WhiteWaterExtentIndividualandDEPTH3(sitename,base_folder,utc,gamma,hsTS,tpTS,dpTS,wave_timeTS,avg_shoreline_points,output):
    
    filepath = base_folder+'data\\'+sitename
    
    

    shoreline = np.array(avg_shoreline_points)
    
    
    

    
    whitewater_extents = []
    whitewater_depth_estimates = []
    whitewater_wave_heights1 = []
    whitewater_wave_heights2 = []
    whitewater_crossShoreDistances = []
    n=-1
    for whitewaterArray in output['WhiteWaterMask']:
        n=n+1
        
        
        
        
        wwRaster = array2raster_InMemory(whitewaterArray,output['SRS'][n],output['SGT'][n]) 


        # get spatial ref
        prj=wwRaster.GetProjection()
        srs=osr.SpatialReference(wkt=prj)
        wwRasterBand = wwRaster.GetRasterBand(1)


        # prep for clipping
        clipinput = base_folder+'data\\'+sitename+'\\transects.shp'
        clipoutput = base_folder+"data\\"+sitename+"\\jpg_files\\detection\\"+output['filename'][n][:-4]+'\\transects_clip.shp'

        # transects to clip
        driverSHP = ogr.GetDriverByName("ESRI Shapefile")
        # shape to clip
        inDataSource = driverSHP.Open(clipinput, 0)
        inLayer = inDataSource.GetLayer()
        inSpatialRef = inLayer.GetSpatialRef()
                
        #  create whitewater polygon to clip transects
        dst_layername = base_folder+ "data\\"+sitename+"\\jpg_files\\detection\\"+output['filename'][n][:-4]+'\\vectorWhiteWater'
        driverMEM = ogr.GetDriverByName("Memory")
        wwPolyDatasource = driverMEM.CreateDataSource( dst_layername + ".shp" )
        wwPolyLayer = wwPolyDatasource.CreateLayer(dst_layername, srs = inSpatialRef )
        gdal.Polygonize( wwRasterBand, wwRasterBand, wwPolyLayer, -1, [], callback=None )
                
        ## clip output
        outDataSource = driverMEM.CreateDataSource(clipoutput)
        outLayer = outDataSource.CreateLayer('FINAL', srs = inSpatialRef, geom_type=ogr.wkbMultiLineString)
        
        # ogr.Layer.Clip(inLayer, wwPolyLayer, outLayer)
        inLayer.Clip(wwPolyLayer, outLayer)
        inDataSource = None

        


        extents = []
        for i in range(0, outLayer.GetFeatureCount()):
            wwclippedTransectFeature = outLayer.GetFeature(i)
            wwclippedTransectGeometry = wwclippedTransectFeature.GetGeometryRef()
            extents.append([wwclippedTransectGeometry.GetPoint(1)[0],wwclippedTransectGeometry.GetPoint(1)[1]])
        crossShoreDistances = []
        for extent in extents:
            cross_shore_distance = np.nanmin(distance(shoreline, extent))
            crossShoreDistances.append(cross_shore_distance)
        whitewater_crossShoreDistances.append(crossShoreDistances)
            
            
        outDataSource = None
        whitewater_extents.append(extents)
        
        
        datestr = output['filename'][n][:13].split('-')
        date = datetime.datetime(year=int(datestr[0]),month=int(datestr[1]),day=int(datestr[2]),hour=int(datestr[3])).replace(tzinfo=utc)
        nearest_date = nearest(wave_timeTS, date)
        a = np.where(wave_timeTS==nearest_date)[0][0]
        waveheight = hsTS[a]
        waveperiod = tpTS[a]
        wavedirection = dpTS[a]
        waveheightI = hsTS[a]
            
            
        depth = 4
            
            
            
        # SHOALING WAVE CODE
        # SHOALING WAVE CODE
        # SHOALING WAVE CODE
        if False:#output['filename'][n].split('-')[5][3:4] == 'S':
            tryAngle = True
            image = ''
            for all in wwRasters[0].split('\\')[:-1]:
                image =os.path.join(image, all)
            image = image+'.png'
            try:
                shoaled_angle = waveangle(sitename, base_folder, image)
                wavedirection2 = angleWithShore(sitename, base_folder, np.tan(np.deg2rad(wavedirection)), base_folder+'data\\'+sitename+"\\"+sitename+"mean_shore_unproj.shp")
                waveheight2 = ShoalWave(waveheight, depth, waveperiod,wavedirection2,shoaled_angle)
                depths = (waveheight2/gamma)
                for n in range(1,10):
                    depth = depths[3]
                    waveheight2 = ShoalWave(waveheight, depth, waveperiod,wavedirection2,shoaled_angle)
                    depths = (waveheight2/gamma)
            except:
                waveheight2 = ShoalWave(waveheight, depth, waveperiod)
                depths = (waveheight2/gamma)
                for n in range(1,10):
                    depth = depths[3]
                    waveheight2 = ShoalWave(waveheight, depth, waveperiod)
                    depths = (waveheight2/gamma)
        # SHOALING WAVE CODE
        # SHOALING WAVE CODE
        # SHOALING WAVE CODE
                    
                   
        else:
            waveheight2 = ShoalWave(waveheight, depth, waveperiod)
            depths = (waveheight2/gamma)
    
            for nnn in range(1,10):
                depth = depths[3]
                waveheight2 = ShoalWave(waveheight, depth, waveperiod)
                depths = (waveheight2/gamma)
        
        
        
        whitewater_depth_estimates.append(depths)
        whitewater_wave_heights1.append(waveheight)
        whitewater_wave_heights2.append(waveheight2)
        
        
        
        
        #     folder = base_folder+"data\\"+sitename+"\\jpg_files\\detection\\"+output['filename'][n][:-4]
        #     cc_save.Save_WhiteWater_Extents(sitename, base_folder,folder+'\\whitewater_point_extent.shp',extents,depths,waveheight,waveheight2,folder)
            
            
        # except:
        #     pass


    return whitewater_extents,whitewater_depth_estimates,whitewater_wave_heights1,whitewater_wave_heights2,whitewater_crossShoreDistances




def collinear(a, b, c):
    "Return true iff a, b, and c all lie on the same line."
    wtv = ((b[0] - a[0] ) * (c[1] - a[1]) - (c[0]  - a[0] ) * (b[1] - a[1]))
    if wtv < 1 and wtv > -1:
        WTV = True
    else:
        WTV = False
    return WTV


def WhiteWaterExtentIndividualandDEPTH4(sitename,base_folder,utc,gamma,hsTS,tpTS,dpTS,wave_timeTS,avg_shoreline_points,output,transectsFull):
    
    filepath = base_folder+'data\\'+sitename
    
    

    shoreline = np.array(avg_shoreline_points)
    
    
    

    
    whitewater_extents = []
    whitewater_depth_estimates = []
    whitewater_wave_heights1 = []
    whitewater_wave_heights2 = []
    whitewater_crossShoreDistances = []
    n=-1
    for whitewaterArray in output['WhiteWaterMask']:
        n=n+1
        
        
        
        
        wwRaster = array2raster_InMemory(whitewaterArray,output['SRS'][n],output['SGT'][n]) 


        # get spatial ref
        prj=wwRaster.GetProjection()
        srs=osr.SpatialReference(wkt=prj)
        wwRasterBand = wwRaster.GetRasterBand(1)


        # prep for clipping
        clipinput = base_folder+'data\\'+sitename+'\\transects.shp'
        clipoutput = base_folder+"data\\"+sitename+"\\jpg_files\\detection\\"+output['filename'][n][:-4]+'\\transects_clip.shp'

        # transects to clip
        driverSHP = ogr.GetDriverByName("ESRI Shapefile")
        # shape to clip
        inDataSource = driverSHP.Open(clipinput, 0)
        inLayer = inDataSource.GetLayer()
        inSpatialRef = inLayer.GetSpatialRef()
                
        #  create whitewater polygon to clip transects
        dst_layername = base_folder+ "data\\"+sitename+"\\jpg_files\\detection\\"+output['filename'][n][:-4]+'\\vectorWhiteWater'
        driverMEM = ogr.GetDriverByName("Memory")
        wwPolyDatasource = driverMEM.CreateDataSource( dst_layername + ".shp" )
        wwPolyLayer = wwPolyDatasource.CreateLayer(dst_layername, srs = inSpatialRef )
        gdal.Polygonize( wwRasterBand, wwRasterBand, wwPolyLayer, -1, [], callback=None )
                
        ## clip output
        outDataSource = driverMEM.CreateDataSource(clipoutput)
        outLayer = outDataSource.CreateLayer('FINAL', srs = inSpatialRef, geom_type=ogr.wkbMultiLineString)
        
        # ogr.Layer.Clip(inLayer, wwPolyLayer, outLayer)
        inLayer.Clip(wwPolyLayer, outLayer)
        inDataSource = None

        

        # extents = []
        # for i in range(0, outLayer.GetFeatureCount()):
        #     wwclippedTransectFeature = outLayer.GetFeature(i)
        #     wwclippedTransectGeometry = wwclippedTransectFeature.GetGeometryRef()
        #     extent_point = [wwclippedTransectGeometry.GetPoint(1)[0],wwclippedTransectGeometry.GetPoint(1)[1]]
        #     for transect in transectsFull:
        #         point1 = transectsFull[transect][0]
        #         point2 = transectsFull[transect][1]
        #         which_transect  = collinear(point1, point2, extent_point)
        #         if which_transect:
        #             extents.append(extent_point)
        #             break
        #         else:
        #             extents.append([np.nan,np.nan])
                    
                    
        extents = np.nan*np.zeros((len(transectsFull),2))
        transectsFull2 = copy.copy(transectsFull)
        for i in range(0, outLayer.GetFeatureCount()):
            wwclippedTransectFeature = outLayer.GetFeature(i)
            wwclippedTransectGeometry = wwclippedTransectFeature.GetGeometryRef()
            extent_point = [wwclippedTransectGeometry.GetPoint(1)[0],wwclippedTransectGeometry.GetPoint(1)[1]]
            for transect in transectsFull2:
                point1 = transectsFull[transect][0]
                point2 = transectsFull[transect][1]
                which_transect  = collinear(point1, point2, extent_point)
                if which_transect:
                    extents[int(transect[9:])-1] = extent_point
                    del transectsFull2[transect]
                    break
        whitewater_extents.append(extents)
                    
                
        crossShoreDistances = []
        for extent in extents:
            if not np.isnan(extent[0]):
                cross_shore_distance = np.nanmin(distance(shoreline, extent))
                crossShoreDistances.append(cross_shore_distance)
            else:
                crossShoreDistances.append(np.nan)
        whitewater_crossShoreDistances.append(crossShoreDistances)
            
            
        outDataSource = None
        
        
        datestr = output['filename'][n][:13].split('-')
        date = datetime.datetime(year=int(datestr[0]),month=int(datestr[1]),day=int(datestr[2]),hour=int(datestr[3])).replace(tzinfo=utc)
        nearest_date = nearest(wave_timeTS, date)
        a = np.where(wave_timeTS==nearest_date)[0][0]
        waveheight = hsTS[a]
        waveperiod = tpTS[a]
        wavedirection = dpTS[a]
        waveheightI = hsTS[a]
            
            
        depth = 4
            
            
            
        # SHOALING WAVE CODE
        # SHOALING WAVE CODE
        # SHOALING WAVE CODE
        if False:#output['filename'][n].split('-')[5][3:4] == 'S':
            tryAngle = True
            image = ''
            for all in wwRasters[0].split('\\')[:-1]:
                image =os.path.join(image, all)
            image = image+'.png'
            try:
                shoaled_angle = waveangle(sitename, base_folder, image)
                wavedirection2 = angleWithShore(sitename, base_folder, np.tan(np.deg2rad(wavedirection)), base_folder+'data\\'+sitename+"\\"+sitename+"mean_shore_unproj.shp")
                waveheight2 = ShoalWave(waveheight, depth, waveperiod,wavedirection2,shoaled_angle)
                depths = (waveheight2/gamma)
                for n in range(1,10):
                    depth = depths[3]
                    waveheight2 = ShoalWave(waveheight, depth, waveperiod,wavedirection2,shoaled_angle)
                    depths = (waveheight2/gamma)
            except:
                waveheight2 = ShoalWave(waveheight, depth, waveperiod)
                depths = (waveheight2/gamma)
                for n in range(1,10):
                    depth = depths[3]
                    waveheight2 = ShoalWave(waveheight, depth, waveperiod)
                    depths = (waveheight2/gamma)
        # SHOALING WAVE CODE
        # SHOALING WAVE CODE
        # SHOALING WAVE CODE
                    
                   
        else:
            waveheight2 = ShoalWave(waveheight, depth, waveperiod)
            depths = (waveheight2/gamma)
    
            for nnn in range(1,10):
                depth = depths[3]
                waveheight2 = ShoalWave(waveheight, depth, waveperiod)
                depths = (waveheight2/gamma)
        
        
        
        whitewater_depth_estimates.append(depths)
        whitewater_wave_heights1.append(waveheight)
        whitewater_wave_heights2.append(waveheight2)
        
        
        
        
        #     folder = base_folder+"data\\"+sitename+"\\jpg_files\\detection\\"+output['filename'][n][:-4]
        #     cc_save.Save_WhiteWater_Extents(sitename, base_folder,folder+'\\whitewater_point_extent.shp',extents,depths,waveheight,waveheight2,folder)
            
            
        # except:
        #     pass




    WhiteWater_Data = {'whitewater_extents':whitewater_extents}
    WhiteWater_Data['whitewater_depth_estimates']=whitewater_depth_estimates
    WhiteWater_Data['whitewater_wave_heights1']=whitewater_wave_heights1
    WhiteWater_Data['whitewater_wave_heights2']=whitewater_wave_heights2
    WhiteWater_Data['whitewater_crossShoreDistances']=whitewater_crossShoreDistances

    # return whitewater_extents,whitewater_depth_estimates,whitewater_wave_heights1,whitewater_wave_heights2,whitewater_crossShoreDistances
    return WhiteWater_Data










def ShoalWave(waveheight, depth, waveperiod, angle1=90,angle2=90):
    g = 9.81
    T = waveperiod
    d = depth
    L = ((g*T**2)/(2*np.pi)) * np.sqrt(np.tanh((4*np.pi**2*d)/(g*T**2)))
    k = (2*np.pi)/L
    n = 0.5 * (1 + ((2*k*d)/(np.sinh(2*k*d))) )
    Ks = np.sqrt((2*n*np.tanh(k*d))**-1)
    Kr = np.sqrt((np.cos(np.deg2rad(angle1)))/(np.cos(np.deg2rad(angle2))))
    newWaveHeight = waveheight * Ks * Kr
    return newWaveHeight



# image ='C:\\Users\\RDCHLNRO\\Desktop\\wave2014\\data\\DUCK1\\jpg_files\\detection\\2016-12-17-15-56-46_S2\\raster_whitewater.tif'
# imagesplit = image.split('\\')
# imagesplit[0] = imagesplit[0] + '\\'
# dat = imagesplit[-2]
# imagesplit = imagesplit[:-2]
# imagesplit.append(dat+'.png')
# image = os.path.join(*imagesplit)



def extractDataSquare(wavesOnly):
    for corner in range(0,4):
        if corner == 0:
            corner0 = wavesOnly[0:2,0:2]
            n=3
            while not np.any(np.isnan(corner0)):
                corner0 = wavesOnly[0:n,0:n]
                n=n+1
        if corner == 1:
            corner1 = wavesOnly[0:2,wavesOnly.shape[1]-2:wavesOnly.shape[1]]
            n=3
            while not np.any(np.isnan(corner1)):
                corner1 = wavesOnly[0:n,wavesOnly.shape[1]-n:wavesOnly.shape[1]]
                n=n+1
        if corner == 2:
            corner2 = wavesOnly[wavesOnly.shape[0]-2:wavesOnly.shape[0],wavesOnly.shape[1]-2:wavesOnly.shape[1]]
            n=3
            while not np.any(np.isnan(corner1)):
                corner2 = wavesOnly[wavesOnly.shape[0]-n:wavesOnly.shape[0],wavesOnly.shape[1]-n:wavesOnly.shape[1]]
                n=n+1
        if corner == 3:
            corner3 = wavesOnly[wavesOnly.shape[0]-2:wavesOnly.shape[0],0:2]
            n=3
            while not np.any(np.isnan(corner0)):
                corner3 = wavesOnly[wavesOnly.shape[0]-n:wavesOnly.shape[0],0:n]
                n=n+1
    sizes = np.array([corner0.size,corner1.size,corner2.size,corner3.size])
    winner = np.where(sizes == np.nanmax(sizes))[0]
    if len(winner)==1:
        if winner ==0:
            return corner0
        if winner ==1:
            return corner1
        if winner ==2:
            return corner2
        if winner ==3:
            return corner3
    else:
        return False




def waveangle(sitename, base_folder, image):
    im = plt.imread(image)
    
    im2 = imageExtracter_thirdImage(im)
    plt.imshow(im2)
    plt.show()
    plt.close()

    array, extent = waterMask(sitename, base_folder)


    d1list = []
    for row in range(0,im2.shape[0]):
        for col in range(0,im2.shape[1]):
            d1list.append(np.sum(im2[row,col]))
    d2array = np.array(d1list).reshape(im2.shape[0],im2.shape[1])
    plt.imshow(d2array)
    plt.show()
    plt.close()
    
    
    plt.imshow(array)
    plt.show()
    plt.close()


    array[array!=array[-1,-1]]=np.nan
    array[array==array[-1,-1]]=1
    
    
    wavesOnly = array*d2array
    
    plt.imshow(wavesOnly)
    plt.show()
    plt.close()

    subsetWavesOnly = extractDataSquare(wavesOnly)
    plt.imshow(subsetWavesOnly)
    plt.show()
    plt.close()




    theta = np.linspace(0., 180., max(subsetWavesOnly.shape), endpoint=False)
    sinogram = radon(subsetWavesOnly, theta=theta, circle=True)



    
    
    slope = theta[np.where(sinogram == np.nanmax(sinogram))[1][0]]
    riseoverrun = np.tan(np.deg2rad(slope))
    xs = np.arange(0,500)
    ys = xs*riseoverrun*-1


    a = angleWithShore(sitename, base_folder, slope, base_folder+'data\\'+sitename+"\\"+sitename+"mean_shore_unproj.shp")


    
    return a




def angleWithShore(sitename, base_folder, slope, shore):
    # get reference shapefile points and then get a subset of those
    sf = shapefile.Reader(shore)
    shapes = sf.shapes()
    points = np.array(shapes[0].points)
    # less_points = points[1::70]
    sf.close()
    
    shoreSlope, intercept, r_value, p_value, std_err = st.linregress(points[:,0], points[:,1])
    
    
    a = np.abs(np.rad2deg(np.arctan((slope-shoreSlope)/(1+slope*shoreSlope))))
    b = 180-a
    
    if b<a:
        a = b

    return a







def waterMask(sitename, base_folder):        
    
    grayarrays = glob.glob(base_folder+ "data\\"+sitename+"\\jpg_files\\detection\\*\\raster_gray.tif")       
    for grayarray in [grayarrays[-1]]:

        src = gdal.Open(grayarray)
        SRS = src.GetProjection()
        ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
        lrx = ulx + (src.RasterXSize * xres)
        lry = uly + (src.RasterYSize * yres)
        extent = [ulx, lrx, lry, uly]
        
        array = np.array(src.GetRasterBand(1).ReadAsArray())
        array2 = np.where((array==2)|(array==3),np.nan,array)
        
        break

    return array2, extent

















def AdjustOutputShapefileByWaves1(slope_ests,output,filepath,sitename,buoy,dissipative,transects):
    
    years = []
    for date in output['dates']:
        years.append(date.year)
    Uyears = np.unique(np.array(years))




    utc = output['dates'][0].tzinfo
    waveheights = []
    waveperiods = []
    waveVdates = []
    for Uyear in Uyears:
        try:
            try:
                data = urllib.request.urlopen("https://www.ndbc.noaa.gov/view_text_file.php?filename="+buoy+"h"+str(Uyear)+".txt.gz&dir=data/historical/stdmet/")        
            except:
                time.sleep(30)
                data = urllib.request.urlopen("https://www.ndbc.noaa.gov/view_text_file.php?filename="+buoy+"h"+str(Uyear)+".txt.gz&dir=data/historical/stdmet/")        
            nn=-1
            heights = []
            periods = []
            ms = []
            ds = []
            hrs=[]
            mins=[]
            wavedatetimes = []
            for line in data:
                nn=nn+1
                if nn>1:
                    line = str(line).replace('  ', ' ')
                    line =      line.replace('   ', ' ')
                    line =      line.replace('    ', ' ')
                    line = line.split(' ')
                    heights.append( float(line[8]) )
                    periods.append( float(line[10]) ) # APD?????????????????
                    ms.append(int(line[1]))
                    ds.append(int(line[2]))
                    hrs.append(int(line[3]))
                    mins.append(int(line[4]))
                    wavedatetimes.append(datetime.datetime(year=Uyear,month=int(line[1]),day=int(line[2]),hour=int(line[3]),minute=int(line[4]),tzinfo=utc))
            heights = np.array(heights)
            heights[heights>30]=np.nan
            heights[heights<0]=np.nan
            periods = np.array(periods)
            periods[periods>90]=np.nan
            periods[periods<0]=np.nan
            
            
            
            for date in output['dates']:
                if date.year == Uyear:
                    near = nearest(wavedatetimes, date)
                    a = np.where(np.array(wavedatetimes)==near)[0][0]
                    waveheights.append(heights[a])
                    waveperiods.append(periods[a])
                    waveVdates.append(wavedatetimes[a])
        except:
            for date in output['dates']:
                if date.year == Uyear:
                    waveheights.append(np.nan)
                    waveperiods.append(np.nan)
                    waveVdates.append(date)
    
    
    
    
    
    # for the chosen point along the reference shoreline find a slope for best fit line
    shoreshape = shapefile.Reader(filepath+"\\"+sitename+"_reference_shoreline.shp")
    shoreRece = shoreshape.shapeRecords()[0]
    shorepoints=shoreRece.shape.points
    shoreshape.close()
    
    deriv_points = np.array(shorepoints)
    x = deriv_points[:,0]
    y = deriv_points[:,1]
    slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
    
    # find a slope perpindicular to shorelne slope for transect
    pslope = -1*(1/slope)
    angle = np.arctan(pslope)
    
    # find beach slope
    slopes = []
    for slope_est in slope_ests:
        slopes.append(slope_ests[slope_est])
    slope = np.nanmean(np.array(slopes))
    
    # Get point from ocean
    shoreRece = list(transects.keys())[int(len(transects)/2)]
    ocean_lat_long=transects[shoreRece][0]
    
    # open trimmed shorelines and start tidally corrected shorelines file
    shoreshape = shapefile.Reader(filepath+"\\"+sitename+"_outputTidalCC.shp")
    shoreRecs = shoreshape.shapeRecords()
    w = shapefile.Writer(filepath+"\\"+sitename+"_outputTidalCRunupC.shp")
    w.fields = shoreshape.fields[1:]    
    
    
    
    
    
    newshores = []
    Vdistances = []
    Vdates = []
    for shoreRec in shoreRecs:
        shorepoints = shoreRec.shape.points
        
        # find tide position
        date = shoreRec.record[0]
        sat = shoreRec.record[1]
        days = []
        months = []
        years = []
        sats = []
        for fdgfd in output['dates']:
            days.append(fdgfd.day==date.day)
            months.append(fdgfd.month==date.month)
            years.append(fdgfd.year==date.year)
        for fdgfd in output['satname']:
            sats.append(fdgfd==sat)
        for n in range(0,len(days)):
            if days[n] and months[n] and years[n] and sats[n]:
                break





    #    %% Stockdon et al., 2006 Model
        
    #    %R2=1.1[0.35Bf(HoLo)^1/2 + {[HoLo(0.563Bf^2 + 0.004]^1/2}/2]
    #    %R2_dissipative=0.043(sqrt(Ho.*Lo))
        
    #    %Bf = foreshore beach slope
    #    %Ho = wave height
    #    %T = wave period
    
        if np.isfinite(waveheights[n]) and np.isfinite(waveperiods[n]):
            T = waveperiods[n]
            Ho = waveheights[n]
            Bf = copy.copy(slope)
            Lo=9.81/(2*np.pi)*T**2
            IB=Bf/np.sqrt(Ho/Lo)
            setup=(0.35*Bf)*np.sqrt(Ho*Lo);
            swash=np.sqrt((Ho*Lo)*(0.563*(Bf**2) + 0.004))
            swashIN=0.75*Bf*np.sqrt(Ho*Lo)
            swashIG=0.06*np.sqrt(Ho*Lo)
            R2=1.1*(setup+(swash/2))
            R2_dissipative=0.043*np.sqrt(Ho*Lo)
            if dissipative:
                Vdistance = -1*R2_dissipative
            else:
                Vdistance = -1*R2
            Vdate = waveVdates[n]
        else:
            Vdistance = 0
            Vdate = date
        Hshift = Vdistance/slope
        outTOsea=True
            
        shorepoints1 = np.array(shorepoints)
        shorepoints2 = np.dstack((shorepoints1[:,0] + np.cos(angle) * Hshift, shorepoints1[:,1] + np.sin(angle) * Hshift))[0]
        shorepoints3 = np.dstack((shorepoints1[:,0] + np.cos(angle) * -Hshift, shorepoints1[:,1] + np.sin(angle) * -Hshift))[0]

        testpoint1 = shorepoints2[int(len(shorepoints2)/2)]
        testpoint2 = shorepoints[int(len(shorepoints)/2)]
        if (distancePointPoint(testpoint1,ocean_lat_long) < distancePointPoint(testpoint2,ocean_lat_long)) and outTOsea:
            NEWshorepoints = shorepoints2
        if (distancePointPoint(testpoint1,ocean_lat_long) < distancePointPoint(testpoint2,ocean_lat_long)) and not outTOsea:
            NEWshorepoints = shorepoints3
        if (distancePointPoint(testpoint1,ocean_lat_long) > distancePointPoint(testpoint2,ocean_lat_long)) and outTOsea:
            NEWshorepoints = shorepoints3
        if (distancePointPoint(testpoint1,ocean_lat_long) > distancePointPoint(testpoint2,ocean_lat_long)) and not outTOsea:
            NEWshorepoints = shorepoints2
        if Vdistance == 0 or Hshift == 0 or Hshift < -100:
            NEWshorepoints = shorepoints1
            Vdistance = 0
            
            
            
        
        w.record(shoreRec.record[0],shoreRec.record[1],shoreRec.record[2],shoreRec.record[3])
        w.line([NEWshorepoints])
        
        Vdistances.append(-1*Vdistance)
        Vdates.append(Vdate)
        newshores.append(NEWshorepoints)
    output['TrimmedCorrectedShorelines+Runup']=newshores
    output['WaveVertical']=Vdistances
    output['WaveDates']=Vdates
        
    w.close()
    shutil.copy(filepath+"\\"+sitename+"_output.prj", filepath+"\\"+sitename+"_outputTidalCRunupC.prj")
    shoreshape.close()
    
    return output






























def AdjustOutputShapefileByWaves2(slope_ests,output,filepath,sitename,buoy,dissipative):
    
    # attempt to use ML to pick when to wave adjust
    
    with open("C:\\Users\\RDCHLNRO\\Desktop\\waves2021\\runup_classifier.pkl", 'rb') as f:
        clf = pickle.load(f) 
    
    
    years = []
    for date in output['dates']:
        years.append(date.year)
    Uyears = np.unique(np.array(years))




    utc = output['dates'][0].tzinfo
    waveheights = []
    waveperiods = []
    for Uyear in Uyears:
        try:
            data = urllib.request.urlopen("https://www.ndbc.noaa.gov/view_text_file.php?filename="+buoy+"h"+str(Uyear)+".txt.gz&dir=data/historical/stdmet/")        
            nn=-1
            heights = []
            periods = []
            ms = []
            ds = []
            hrs=[]
            mins=[]
            wavedatetimes = []
            for line in data:
                nn=nn+1
                if nn>1:
                    line = str(line).replace('  ', ' ')
                    line =      line.replace('   ', ' ')
                    line =      line.replace('    ', ' ')
                    line = line.split(' ')
                    heights.append( float(line[8]) )
                    periods.append( float(line[10]) ) # APD?????????????????
                    ms.append(int(line[1]))
                    ds.append(int(line[2]))
                    hrs.append(int(line[3]))
                    mins.append(int(line[4]))
                    wavedatetimes.append(datetime.datetime(year=Uyear,month=int(line[1]),day=int(line[2]),hour=int(line[3]),minute=int(line[4]),tzinfo=utc))
            heights = np.array(heights)
            heights[heights>30]=np.nan
            heights[heights<0]=np.nan
            periods = np.array(periods)
            periods[periods>90]=np.nan
            periods[periods<0]=np.nan
            
            
            
            for date in output['dates']:
                if date.year == Uyear:
                    near = nearest(wavedatetimes, date)
                    a = np.where(np.array(wavedatetimes)==near)[0][0]
                    waveheights.append(heights[a])
                    waveperiods.append(periods[a])
        except:
            for date in output['dates']:
                if date.year == Uyear:
                    waveheights.append(np.nan)
                    waveperiods.append(np.nan)
    
    
    
    
    
    
    # for the chosen point along the reference shoreline find a slope for best fit line
    shoreshape = shapefile.Reader(filepath+"\\"+sitename+"_reference_shoreline.shp")
    shoreRece = shoreshape.shapeRecords()[0]
    shorepoints=shoreRece.shape.points
    shoreshape.close()
    
    deriv_points = np.array(shorepoints)
    x = deriv_points[:,0]
    y = deriv_points[:,1]
    slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
    
    # find a slope perpindicular to shorelne slope for transect
    pslope = -1*(1/slope)
    angle = np.arctan(pslope)
    
    # find beach slope
    slopes = []
    for slope_est in slope_ests:
        slopes.append(slope_ests[slope_est])
    slope = np.nanmean(np.array(slopes))
    
    # Get point from ocean
    transects = shapefile.Reader(filepath+"\\transects.shp")
    shoreRece = transects.shapeRecords()[int(len(transects.shapeRecords())/2)]
    ocean_lat_long=shoreRece.shape.points[1]
    transects.close()
    
    # open trimmed shorelines and start tidally corrected shorelines file
    shoreshape = shapefile.Reader(filepath+"\\"+sitename+"_outputTidalCC.shp")
    shoreRecs = shoreshape.shapeRecords()
    w = shapefile.Writer(filepath+"\\"+sitename+"_outputTidalCRunupC.shp")
    w.fields = shoreshape.fields[1:]    
    
    
    
    
    
    newshores=[]
    for shoreRec in shoreRecs:
        shorepoints = shoreRec.shape.points
        
        # find tide position
        date = shoreRec.record[0]
        sat = shoreRec.record[1]
        days = []
        months = []
        years = []
        sats = []
        for fdgfd in output['dates']:
            days.append(fdgfd.day==date.day)
            months.append(fdgfd.month==date.month)
            years.append(fdgfd.year==date.year)
        for fdgfd in output['satname']:
            sats.append(fdgfd==sat)
        for n in range(0,len(days)):
            if days[n] and months[n] and years[n] and sats[n]:
                break





    #    %% Stockdon et al., 2006 Model
        
    #    %R2=1.1[0.35Bf(HoLo)^1/2 + {[HoLo(0.563Bf^2 + 0.004]^1/2}/2]
    #    %R2_dissipative=0.043(sqrt(Ho.*Lo))
        
    #    %Bf = foreshore beach slope
    #    %Ho = wave height
    #    %T = wave period
    
        try:
            zz = clf.predict(np.dstack((waveheights[n], waveperiods[n]))[0])
            if zz[0] != 0:
                print('REJECTED!')
        except:
            zz = [0]
        if np.isfinite(waveheights[n]) and np.isfinite(waveperiods[n]) and zz[0] != 0:
            T = waveperiods[n]
            Ho = waveheights[n]
            Bf = copy.copy(slope)
            Lo=9.81/(2*np.pi)*T**2
            IB=Bf/np.sqrt(Ho/Lo)
            setup=(0.35*Bf)*np.sqrt(Ho*Lo);
            swash=np.sqrt((Ho*Lo)*(0.563*(Bf**2) + 0.004))
            swashIN=0.75*Bf*np.sqrt(Ho*Lo)
            swashIG=0.06*np.sqrt(Ho*Lo)
            R2=1.1*(setup+(swash/2))
            R2_dissipative=0.043*np.sqrt(Ho*Lo)
            if dissipative:
                Vdistance = -1*R2_dissipative
            else:
                Vdistance = -1*R2
        else:
            Vdistance = 0
        
        Hshift = Vdistance/slope
        outTOsea=True
            
        shorepoints1 = np.array(shorepoints)
        shorepoints2 = np.dstack((shorepoints1[:,0] + np.cos(angle) * Hshift, shorepoints1[:,1] + np.sin(angle) * Hshift))[0]
        shorepoints3 = np.dstack((shorepoints1[:,0] + np.cos(angle) * -Hshift, shorepoints1[:,1] + np.sin(angle) * -Hshift))[0]

        testpoint1 = shorepoints2[int(len(shorepoints2)/2)]
        testpoint2 = shorepoints[int(len(shorepoints)/2)]
        if (distancePointPoint(testpoint1,ocean_lat_long) < distancePointPoint(testpoint2,ocean_lat_long)) and outTOsea:
            NEWshorepoints = shorepoints2
        if (distancePointPoint(testpoint1,ocean_lat_long) < distancePointPoint(testpoint2,ocean_lat_long)) and not outTOsea:
            NEWshorepoints = shorepoints3
        if (distancePointPoint(testpoint1,ocean_lat_long) > distancePointPoint(testpoint2,ocean_lat_long)) and outTOsea:
            NEWshorepoints = shorepoints3
        if (distancePointPoint(testpoint1,ocean_lat_long) > distancePointPoint(testpoint2,ocean_lat_long)) and not outTOsea:
            NEWshorepoints = shorepoints2
        if Vdistance == 0 or Hshift == 0:
            NEWshorepoints = shorepoints1
            
            
        
        if Hshift < -100:
            break
        
        w.record(shoreRec.record[0],shoreRec.record[1],shoreRec.record[2],shoreRec.record[3])
        w.line([NEWshorepoints])
        
        newshores.append(NEWshorepoints)
    output['TrimmedCorrectedShorelines+Runup']=newshores
        
        
    w.close()
    shutil.copy(filepath+"\\"+sitename+"_output.prj", filepath+"\\"+sitename+"_outputTidalCRunupC.prj")
    shoreshape.close()
    
    return output





















def AdjustOutputShapefileByWaves3(slope_ests,output,filepath,sitename,buoy,dissipative):
    #1.5m wave cutoff for duck only
        
    years = []
    for date in output['dates']:
        years.append(date.year)
    Uyears = np.unique(np.array(years))




    utc = output['dates'][0].tzinfo
    waveheights = []
    waveperiods = []
    for Uyear in Uyears:
        try:
            data = urllib.request.urlopen("https://www.ndbc.noaa.gov/view_text_file.php?filename="+buoy+"h"+str(Uyear)+".txt.gz&dir=data/historical/stdmet/")        
            nn=-1
            heights = []
            periods = []
            ms = []
            ds = []
            hrs=[]
            mins=[]
            wavedatetimes = []
            for line in data:
                nn=nn+1
                if nn>1:
                    line = str(line).replace('  ', ' ')
                    line =      line.replace('   ', ' ')
                    line =      line.replace('    ', ' ')
                    line = line.split(' ')
                    heights.append( float(line[8]) )
                    periods.append( float(line[10]) ) # APD?????????????????
                    ms.append(int(line[1]))
                    ds.append(int(line[2]))
                    hrs.append(int(line[3]))
                    mins.append(int(line[4]))
                    wavedatetimes.append(datetime.datetime(year=Uyear,month=int(line[1]),day=int(line[2]),hour=int(line[3]),minute=int(line[4]),tzinfo=utc))
            heights = np.array(heights)
            heights[heights>30]=np.nan
            heights[heights<0]=np.nan
            periods = np.array(periods)
            periods[periods>90]=np.nan
            periods[periods<0]=np.nan
            
            
            
            for date in output['dates']:
                if date.year == Uyear:
                    near = nearest(wavedatetimes, date)
                    a = np.where(np.array(wavedatetimes)==near)[0][0]
                    waveheights.append(heights[a])
                    waveperiods.append(periods[a])
        except:
            for date in output['dates']:
                if date.year == Uyear:
                    waveheights.append(np.nan)
                    waveperiods.append(np.nan)
    
    
    
    
    
    
    # for the chosen point along the reference shoreline find a slope for best fit line
    shoreshape = shapefile.Reader(filepath+"\\"+sitename+"_reference_shoreline.shp")
    shoreRece = shoreshape.shapeRecords()[0]
    shorepoints=shoreRece.shape.points
    shoreshape.close()
    
    deriv_points = np.array(shorepoints)
    x = deriv_points[:,0]
    y = deriv_points[:,1]
    slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
    
    # find a slope perpindicular to shorelne slope for transect
    pslope = -1*(1/slope)
    angle = np.arctan(pslope)
    
    # find beach slope
    slopes = []
    for slope_est in slope_ests:
        slopes.append(slope_ests[slope_est])
    slope = np.nanmean(np.array(slopes))
    
    # Get point from ocean
    transects = shapefile.Reader(filepath+"\\transects.shp")
    shoreRece = transects.shapeRecords()[int(len(transects.shapeRecords())/2)]
    ocean_lat_long=shoreRece.shape.points[1]
    transects.close()
    
    # open trimmed shorelines and start tidally corrected shorelines file
    shoreshape = shapefile.Reader(filepath+"\\"+sitename+"_outputTidalCC.shp")
    shoreRecs = shoreshape.shapeRecords()
    w = shapefile.Writer(filepath+"\\"+sitename+"_outputTidalCRunupC.shp")
    w.fields = shoreshape.fields[1:]    
    
    
    
    
    
    newshores=[]
    for shoreRec in shoreRecs:
        shorepoints = shoreRec.shape.points
        
        # find tide position
        date = shoreRec.record[0]
        sat = shoreRec.record[1]
        days = []
        months = []
        years = []
        sats = []
        for fdgfd in output['dates']:
            days.append(fdgfd.day==date.day)
            months.append(fdgfd.month==date.month)
            years.append(fdgfd.year==date.year)
        for fdgfd in output['satname']:
            sats.append(fdgfd==sat)
        for n in range(0,len(days)):
            if days[n] and months[n] and years[n] and sats[n]:
                break





    #    %% Stockdon et al., 2006 Model
        
    #    %R2=1.1[0.35Bf(HoLo)^1/2 + {[HoLo(0.563Bf^2 + 0.004]^1/2}/2]
    #    %R2_dissipative=0.043(sqrt(Ho.*Lo))
        
    #    %Bf = foreshore beach slope
    #    %Ho = wave height
    #    %T = wave period
    
        if np.isfinite(waveheights[n]) and np.isfinite(waveperiods[n]):
            T = waveperiods[n]
            Ho = waveheights[n]
            Bf = copy.copy(slope)
            Lo=9.81/(2*np.pi)*T**2
            IB=Bf/np.sqrt(Ho/Lo)
            setup=(0.35*Bf)*np.sqrt(Ho*Lo);
            swash=np.sqrt((Ho*Lo)*(0.563*(Bf**2) + 0.004))
            swashIN=0.75*Bf*np.sqrt(Ho*Lo)
            swashIG=0.06*np.sqrt(Ho*Lo)
            R2=1.1*(setup+(swash/2))
            R2_dissipative=0.043*np.sqrt(Ho*Lo)
            if dissipative:
                Vdistance = -1*R2_dissipative
            else:
                Vdistance = -1*R2
        else:
            Ho=0
            Vdistance = 0
        
        Hshift = Vdistance/slope
        outTOsea=True
            
        shorepoints1 = np.array(shorepoints)
        shorepoints2 = np.dstack((shorepoints1[:,0] + np.cos(angle) * Hshift, shorepoints1[:,1] + np.sin(angle) * Hshift))[0]
        shorepoints3 = np.dstack((shorepoints1[:,0] + np.cos(angle) * -Hshift, shorepoints1[:,1] + np.sin(angle) * -Hshift))[0]

        testpoint1 = shorepoints2[int(len(shorepoints2)/2)]
        testpoint2 = shorepoints[int(len(shorepoints)/2)]
        if (distancePointPoint(testpoint1,ocean_lat_long) < distancePointPoint(testpoint2,ocean_lat_long)) and outTOsea:
            NEWshorepoints = shorepoints2
        if (distancePointPoint(testpoint1,ocean_lat_long) < distancePointPoint(testpoint2,ocean_lat_long)) and not outTOsea:
            NEWshorepoints = shorepoints3
        if (distancePointPoint(testpoint1,ocean_lat_long) > distancePointPoint(testpoint2,ocean_lat_long)) and outTOsea:
            NEWshorepoints = shorepoints3
        if (distancePointPoint(testpoint1,ocean_lat_long) > distancePointPoint(testpoint2,ocean_lat_long)) and not outTOsea:
            NEWshorepoints = shorepoints2
        if Vdistance == 0 or Hshift == 0:
            NEWshorepoints = shorepoints1
        if Ho<1.5:
            NEWshorepoints = shorepoints1
            
        
        if Hshift < -100:
            break
        
        
        
        
        w.record(shoreRec.record[0],shoreRec.record[1],shoreRec.record[2],shoreRec.record[3])
        w.line([NEWshorepoints])
        
        newshores.append(NEWshorepoints)
    output['TrimmedCorrectedShorelines+Runup']=newshores
        
        
    w.close()
    shutil.copy(filepath+"\\"+sitename+"_output.prj", filepath+"\\"+sitename+"_outputTidalCRunupC.prj")
    shoreshape.close()
    
    return output






def SyncOutputDictShape(output, TCC_File):
    
    checkUniques = []
    
    originalkeys = output.keys()
    outputNew = {}
    for originalkey in originalkeys:
        outputNew[originalkey]=[]

    sf = shapefile.Reader(TCC_File)
    shaperecs = sf.shapeRecords()
    for shaperec in shaperecs:
        date = shaperec.record[0]
        sat = shaperec.record[1]
        
        
        done = []
        a = -1
        for outdate in output['dates']:
            a=a+1
            outsat = output['satname'][a]
            code = str(date.year) + str(date.month) + str(date.day) + sat
            if outdate.year==date.year and outdate.month==date.month and outdate.day==date.day and outsat == sat and np.nanmean(output['shorelines'][a]) not in checkUniques and code not in done:
                done.append(code)
                checkUniques.append(np.nanmean(output['shorelines'][a]))
                for originalkey in originalkeys:
                    outputNew[originalkey].append(output[originalkey][a])
    sf.close()
    return outputNew

    
    
    
    