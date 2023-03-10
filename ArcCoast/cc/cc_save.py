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
import numpy as np
import os, sys, glob, shutil, subprocess, copy
import matplotlib.pyplot as plt
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
from osgeo import gdal
from osgeo import osr
import numpy as np
import os, sys, glob, shutil
import matplotlib.pyplot as plt
from osgeo import ogr
import shapefile
import subprocess
import copy
from osgeo import gdal
from osgeo import osr
import numpy as np
import os, sys, glob, shutil
import matplotlib.pyplot as plt
from osgeo import ogr
import shapefile
import subprocess
import copy
import scipy.stats as st
from scipy import ndimage
import dask.array as da
from scipy.spatial.transform import Rotation as R
import pandas as pd
from datetime import timedelta
import pytz
from osgeo.gdalconst import GA_Update


import cc_functions




def array2raster(newRasterfn,array,SRS,SGT):
    
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, array.shape[1], array.shape[0], 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform(SGT)
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRaster.SetProjection(SRS)
    outband.FlushCache()


def Save_Tides(base_folder, sitename, elevations, elevations_dates):
    return None
    

def Save_Waves(base_folder, sitename, hsTS, tpTS, dpTS, wave_timeTS):
    np.savetxt(base_folder+'data\\'+sitename+'\\waveheights.txt',hsTS)
    np.savetxt(base_folder+'data\\'+sitename+'\\waveperiods.txt',tpTS)
    np.savetxt(base_folder+'data\\'+sitename+'\\wavedirections.txt',dpTS)
    pickle.dump( wave_timeTS, open( base_folder+'data\\'+sitename+'\\wavedates.pkl', "wb" ) )


def Save_Classification_Rasters(sitename, base_folder, image, wwArray, classArray, SRS, SGT):

    path = base_folder+ "data\\"+sitename+"\\jpg_files\\detection"

    # make folder for shapefiles
    new_fold = image.split('\\')[-1][:-4]
    new_fold = path+'\\'+new_fold
    os.mkdir(new_fold)

    newRasterfn = new_fold+"\\raster_gray.tif"
    array2raster(newRasterfn,classArray,SRS,SGT)

    newRasterfn = new_fold+"\\raster_whitewater.tif"
    array2raster(newRasterfn,wwArray,SRS,SGT)


def Save_WhiteWater_Extents(sitename, base_folder,shapefile_name,extents,depths,waveheight,waveheight2,folder):

    # w = shapefile.Writer(shapefile_name, shapeType=1)
    # w.field('cross_shore_distances', 'F', decimal=10)
    # w.field('waveheight initial', 'F', decimal=10)
    # w.field('waveheight shoaled', 'F', decimal=10)
    # w.field('depths5', 'F', decimal=10)
    # w.field('depths6', 'F', decimal=10)
    # w.field('depths7', 'F', decimal=10)
    # w.field('depths78', 'F', decimal=10)
    # w.field('depths8', 'F', decimal=10)
    # w.field('depths9', 'F', decimal=10)
    # w.field('depths1', 'F', decimal=10)
    
        
    # for extent in extents:
    #     print(extent)
    #     cross_shore_distance = np.nanmin(distance(shoreline, extent))
    #     w.point(extent[0],extent[1]) 
    #     w.record(cross_shore_distance,waveheight,waveheight2,depths[0],depths[1],depths[2],depths[3],depths[4],depths[5],depths[6])
    # w.close()
    # shutil.copyfile(base_folder+'data\\'+sitename+"\\"+sitename+"_reference_shoreline.prj",  folder+'\\whitewater_point_extent.prj')
    return 0



def Save_Average_Shoreline(filepath,sitename,avg_shoreline_points):
    w = shapefile.Writer(filepath+"\\"+sitename+"_MeanShoreline.shp", shapeType=8)
    w.field('name', 'C')
    w.multipoint(avg_shoreline_points) 
    w.record('1')
    w.close()
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\"+sitename+"_MeanShoreline.prj")

def MaxOuterWaveBreak(filepath,sitename,transectt):
    w = shapefile.Writer(filepath+"\\"+sitename+"_MaxOuterWaveBreak.shp", shapeType=8)
    w.field('name', 'C')
    w.multipoint(transectt) 
    w.record('1')
    w.close()
    shutil.copyfile(filepath+"\\"+sitename+"_reference_shoreline.prj", filepath+"\\"+sitename+"_MaxOuterWaveBreak.prj")




def SaveOutput(base_folder, sitename, output):
    pickle.dump( output, open( base_folder+'data\\'+sitename+'\\output.pkl', "wb" ) )
    # pickle.dump( WhiteWater_Data, open( base_folder+'data\\'+sitename+'\\WhiteWater_Data.pkl', "wb" ) )
    # pickle.dump( Wave_Data, open( base_folder+'data\\'+sitename+'\\Wave_Data.pkl', "wb" ) )

