
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












def Wave_Plot(base_folder,sitename,Wave_Data):
        
    colorHighBIG=np.nanmax(np.array(Wave_Data['meanWaveHeightsByDirectionBandsSeasonallyBIG']))
    colorLowBIG=np.nanmin(np.array(Wave_Data['meanWaveHeightsByDirectionBandsSeasonallyBIG']))
    
    cmapBIG = mpl.cm.cool
    normBIG = mpl.colors.Normalize(vmin=colorLowBIG, vmax=colorHighBIG)



    colorHigh=np.nanmax(np.array(Wave_Data['meanWaveHeightsByDirectionBandsSeasonally']))
    colorLow=np.nanmin(np.array(Wave_Data['meanWaveHeightsByDirectionBandsSeasonally']))
    
    cmap = mpl.cm.cool
    norm = mpl.colors.Normalize(vmin=colorLow, vmax=colorHigh)

    season_names = ['Jan Feb Mar','Apr May Jun','Jul Aug Sep','Oct Nov Dec']
    
    import matplotlib
    matplotlib.rcParams.update({'font.size': 22})
    
    
    ################################# FORMATTING THE FIGURE
    ################################# FORMATTING THE FIGURE
    ################################# FORMATTING THE FIGURE
    ################################# FORMATTING THE FIGURE
    ################################# FORMATTING THE FIGURE
    ################################# FORMATTING THE FIGURE
    ################################# FORMATTING THE FIGURE
    ################################# FORMATTING THE FIGURE
    ################################# FORMATTING THE FIGURE
    
    
    gs=GridSpec(8,4) # 5 rows, 4 columns
    gs.update(hspace = .75)
    
    fig= plt.figure(figsize=(20,25))
    ax0 = plt.subplot(gs[0:2,:])
    ax1 = plt.subplot(gs[2:4,0], projection='polar' )
    ax2 = plt.subplot(gs[2:4,1], projection='polar' )
    ax3 = plt.subplot(gs[2:4,2], projection='polar' )
    ax4 = plt.subplot(gs[2:4,3], projection='polar')
    ax5 = plt.subplot(gs[4,1:3])
    ax6 = plt.subplot(gs[5:7,0], projection='polar' )
    ax7 = plt.subplot(gs[5:7,1], projection='polar' )
    ax8 = plt.subplot(gs[5:7,2], projection='polar' )
    ax9 = plt.subplot(gs[5:7,3], projection='polar')
    ax10 = plt.subplot(gs[7,1:3])
    
    
    
    
    # fig.suptitle('Latitude: ' + str(lat) + ', Longitude: '+str(long))
    
    # graph1
    ax0.set_xlabel('Month')
    ax0.set_ylabel('Mean Hs (m)')
    ax0.plot(range(1,13), Wave_Data['wave_heights_monthly_means'], label = "Mean Wave Height", color= 'red')
    a = np.array(Wave_Data['wave_heights_monthly_means'])-np.array(Wave_Data['wave_heights_monthly_stdevs'])
    b = np.array(Wave_Data['wave_heights_monthly_means'])+np.array(Wave_Data['wave_heights_monthly_stdevs'])
    ax0.fill_between(range(1,13),a ,b, color= 'red', alpha=.1, label='1 St.Dev.')
    
    ax0.plot(range(1,13), Wave_Data['wave_heights_monthly_exceed2'], label = "2% Exceedance", color= 'black',alpha=0.25)
    ax0.plot(range(1,13), Wave_Data['wave_heights_monthly_exceed10'], label = "10% Exceedance", color= 'black',alpha=0.5)
    
    ax0_2 = ax0.twinx()
    ax0_2.set_xlabel('Month')
    ax0_2.set_ylabel('Mean Period (s)')
    ax0_2.plot(range(1,13), Wave_Data['wave_periods_monthly_means'], label = "Mean Wave Period", color='blue')
    ax0_2.tick_params(axis='y', labelcolor = 'blue')
    
    
    #graph2
    rmax = np.nanmax(np.array(Wave_Data['WaveFrequencyByDirectionBandsSeasonally']))
    n=-1
    for month in range(0,4):
        n=n+1
        if month==0:
            for a in range(0, len(Wave_Data['directionalBands'])-1):
                bandLow = Wave_Data['directionalBands'][a]
                Rr=Wave_Data['meanWaveHeightsByDirectionBandsSeasonally'][month][a]
                ax1.bar(np.deg2rad(bandLow), Wave_Data['WaveFrequencyByDirectionBandsSeasonally'][month][a],width=.3, color=cmap(norm(Rr)))
            ax1.set_theta_direction(-1)
            ax1.set_theta_offset(np.pi/2.0)
            ax1.set_title(season_names[n])
            ax1.set_rmax(rmax)
            ax1.set_rmin(0)
            ax1.set_xticklabels([])
            # if (month-1)==0:
            #     cax = f.add_axes()
            #f.colorbar(im)
        if month==1:
            for a in range(0, len(Wave_Data['directionalBands'])-1):
                bandLow = Wave_Data['directionalBands'][a]
                Rr=Wave_Data['meanWaveHeightsByDirectionBandsSeasonally'][month][a]
                ax2.bar(np.deg2rad(bandLow), Wave_Data['WaveFrequencyByDirectionBandsSeasonally'][month][a],width=.3, color=cmap(norm(Rr)))
            ax2.set_theta_direction(-1)
            ax2.set_theta_offset(np.pi/2.0)
            ax2.set_title(season_names[n])
            ax2.set_rmax(rmax)
            ax2.set_rmin(0)
            ax2.set_xticklabels([])
            ax2.set_yticklabels([])
            # if (month-1)==0:
            #     cax = f.add_axes()
            #f.colorbar(im)
        if month==2:
            for a in range(0, len(Wave_Data['directionalBands'])-1):
                bandLow = Wave_Data['directionalBands'][a]
                Rr=Wave_Data['meanWaveHeightsByDirectionBandsSeasonally'][month][a]
                ax3.bar(np.deg2rad(bandLow), Wave_Data['WaveFrequencyByDirectionBandsSeasonally'][month][a],width=.3, color=cmap(norm(Rr)))
            ax3.set_theta_direction(-1)
            ax3.set_theta_offset(np.pi/2.0)
            ax3.set_title(season_names[n])
            ax3.set_rmax(rmax)
            ax3.set_rmin(0)
            ax3.set_xticklabels([])
            ax3.set_yticklabels([])
            # if (month-1)==0:
            #     cax = f.add_axes()
            #f.colorbar(im)
        if month==3:
            for a in range(0, len(Wave_Data['directionalBands'])-1):
                bandLow = Wave_Data['directionalBands'][a]
                Rr=Wave_Data['meanWaveHeightsByDirectionBandsSeasonally'][month][a]
                ax4.bar(np.deg2rad(bandLow), Wave_Data['WaveFrequencyByDirectionBandsSeasonally'][month][a],width=.3, color=cmap(norm(Rr)))
            ax4.set_theta_direction(-1)
            ax4.set_theta_offset(np.pi/2.0)
            ax4.set_title(season_names[n])
            ax4.set_rmax(rmax)
            ax4.set_rmin(0)
            ax4.set_xticklabels([])
            ax4.set_yticklabels([])
            # if (month-1)==0:
            #     cax = f.add_axes()
            #f.colorbar(im)
    
    
    
    # graph3
    cb1=mpl.colorbar.ColorbarBase(ax5,cmap=cmap,norm=norm,orientation='horizontal')
    cb1.set_label('Mean Significant Wave Height (m)')
    
    
    
    
    #graph4
    rmax = np.nanmax(np.array(Wave_Data['WaveFrequencyByDirectionBandsSeasonallyBIG']))
    n=-1
    for month in range(0,4):
        n=n+1
        if month==0:
            for a in range(0, len(Wave_Data['directionalBands'])-1):
                bandLow = Wave_Data['directionalBands'][a]
                Rr=Wave_Data['meanWaveHeightsByDirectionBandsSeasonallyBIG'][month][a]
                ax6.bar(np.deg2rad(bandLow), Wave_Data['WaveFrequencyByDirectionBandsSeasonallyBIG'][month][a],width=.3, color=cmapBIG(normBIG(Rr)))
            ax6.set_theta_direction(-1)
            ax6.set_theta_offset(np.pi/2.0)
            ax6.set_title(season_names[n])
            ax6.set_rmax(rmax)
            ax6.set_rmin(0)
            ax6.set_xticklabels([])
            # if (month-1)==0:
            #     cax = f.add_axes()
            #f.colorbar(im)
        if month==1:
            for a in range(0, len(Wave_Data['directionalBands'])-1):
                bandLow = Wave_Data['directionalBands'][a]
                Rr=Wave_Data['meanWaveHeightsByDirectionBandsSeasonallyBIG'][month][a]
                ax7.bar(np.deg2rad(bandLow), Wave_Data['WaveFrequencyByDirectionBandsSeasonallyBIG'][month][a],width=.3, color=cmapBIG(normBIG(Rr)))
            ax7.set_theta_direction(-1)
            ax7.set_theta_offset(np.pi/2.0)
            ax7.set_title(season_names[n])
            ax7.set_rmax(rmax)
            ax7.set_rmin(0)
            ax7.set_xticklabels([])
            ax7.set_yticklabels([])
            # if (month-1)==0:
            #     cax = f.add_axes()
            #f.colorbar(im)
        if month==2:
            for a in range(0, len(Wave_Data['directionalBands'])-1):
                bandLow = Wave_Data['directionalBands'][a]
                Rr=Wave_Data['meanWaveHeightsByDirectionBandsSeasonallyBIG'][month][a]
                ax8.bar(np.deg2rad(bandLow), Wave_Data['WaveFrequencyByDirectionBandsSeasonallyBIG'][month][a],width=.3, color=cmapBIG(normBIG(Rr)))
            ax8.set_theta_direction(-1)
            ax8.set_theta_offset(np.pi/2.0)
            ax8.set_title(season_names[n])
            ax8.set_rmax(rmax)
            ax8.set_rmin(0)
            ax8.set_xticklabels([])
            ax8.set_yticklabels([])
            # if (month-1)==0:
            #     cax = f.add_axes()
            #f.colorbar(im)
        if month==3:
            for a in range(0, len(Wave_Data['directionalBands'])-1):
                bandLow = Wave_Data['directionalBands'][a]
                Rr=Wave_Data['meanWaveHeightsByDirectionBandsSeasonallyBIG'][month][a]
                ax9.bar(np.deg2rad(bandLow), Wave_Data['WaveFrequencyByDirectionBandsSeasonallyBIG'][month][a],width=.3, color=cmapBIG(normBIG(Rr)))
            ax9.set_theta_direction(-1)
            ax9.set_theta_offset(np.pi/2.0)
            ax9.set_title(season_names[n])
            ax9.set_rmax(rmax)
            ax9.set_rmin(0)
            ax9.set_xticklabels([])
            ax9.set_yticklabels([])
            # if (month-1)==0:
            #     cax = f.add_axes()
            #f.colorbar(im)
    
    
    
    # graph5
    cb1=mpl.colorbar.ColorbarBase(ax10,cmap=cmapBIG,norm=normBIG,orientation='horizontal')
    cb1.set_label('Mean Significant Large Wave Height (m)')
    plt.savefig(base_folder+'data\\'+sitename+'\\01waveClimate.png')
    
    plt.close('all')
    
    
    
    
    
    



def Shore_OverLays_Plot(base_folder,sitename,TC,transects,slopes,output,eventName='',basemap=1):

    
    plt.figure(figsize=(25,25))
    
    
    if basemap==1:
        # plot basemap
        im,extent = cc_functions.get_basemap(output)
        plt.imshow(im, extent = extent)
    elif basemap!=1:
        plt.imshow(basemap[0], extent = basemap[1],alpha=0.5)
    
    
    # # plot water mask
    # watermask = np.where((output['ClassificationMask'][last_s2]==1)|(output['ClassificationMask'][last_s2]==0),np.nan,copy.copy(output['ClassificationMask'][last_s2]))
    # extentwatermask = output['Extents'][last_s2]
    # plt.imshow(watermask, extent = extentwatermask, cmap ='coolwarm_r' ,vmin=0,vmax=.0001)
    
    
    # # plot whitewaters
    # alphaVal = 3/len(output['WhiteWaterMask'])
    # n=-1
    # for whitewaterArray in output['WhiteWaterMask']:
    #     n=n+1
    #     whitewaterArray[whitewaterArray == 0] = np.nan
    #     plt.imshow(whitewaterArray, alpha=alphaVal, extent = output['Extents'][n], cmap='spring_r',vmin=0,vmax=.0001)
    

    # # plot mean shoreline
    # plt.plot(np.array(avg_shoreline_points)[:,0],np.array(avg_shoreline_points)[:,1],c='black',lw=5)
    
    
    # # # plot mean outer wave breaking
    # cc_functions.Average_Shoreline2(filepath_data, sitename, base_folder)
    # sf = shapefile.Reader(filepath+"\\"+sitename+"_MeanOuterWaveBreak.shp")
    # shapes = sf.shapes()
    # for shape in shapes:
    #     points = np.array(shape.points)
    #     plt.plot(points[:,0],points[:,1],c='red',lw=5)
    # sf.close()
    
    
    # # # plot max outer wave breaking
    # if len(WWextentLine)!=0:
    #     plt.scatter(np.array(WWextentLine)[:,0],np.array(WWextentLine)[:,1],c='red',s=10)
    
    
    
    # # # plot max outer wave breaking 10th degree polynomial
    # sf = shapefile.Reader(filepath+"\\"+sitename+"_MaxOuterWaveBreak.shp")
    # shapes = sf.shapes()
    # points = np.array(shapes[0].points)
    # sf.close()
    # points2 = points[~np.isnan(points).any(axis=1)]
    # poly = np.polyfit(points2[:,0], points2[:,1], 20)
    # polynomial = np.poly1d(poly)
    # plt.plot(points2[:,0], polynomial(points2[:,0]), c='red')
    
    
    # # # plot max outer wave breaking interpolation
    # from scipy import interpolate
    # sf = shapefile.Reader(filepath+"\\"+sitename+"_MaxOuterWaveBreak.shp")
    # shapes = sf.shapes()
    # points = np.array(shapes[0].points)
    # sf.close()
    # points2 = points[~np.isnan(points).any(axis=1)]
    # interp = interpolate.interp1d(points2[::100][:,0], points2[::100][:,1],kind='cubic')
    # points = np.linspace(np.nanmin(points2[::100][:,0]),np.nanmax(points2[::100][:,0]),1000)
    # interpolation = interp(points)
    # plt.plot(points, interpolation, c='red')
    
    
    # # plot max outer wave breaking connected
    # sf = shapefile.Reader(base_folder+'data\\'+sitename+"\\"+sitename+"_MaxOuterWaveBreak.shp")
    # shapes = sf.shapes()
    # points = np.array(shapes[0].points)
    # sf.close()
    # points2 = points[~np.isnan(points).any(axis=1)]
    # plt.plot(points2[:,0],points2[:,1] , c='red')
    
    
    # # plot all shorelines
    # sf = shapefile.Reader(base_folder+'data\\'+sitename+'\\'+sitename+'_output.shp')
    # shapes = sf.shapes()
    # for shape in shapes:
    #     points = np.array(shape.points)
    #     plt.plot(points[:,0],points[:,1])
    # sf.close()
    
    
    # plot transect points
    n = -1
    x=[]
    y=[]
    for point in slopes:
        n=n+1
        x.append(transects['Transect '+str(n+1)][0][0])
        y.append(transects['Transect '+str(n+1)][0][1])
    plt.scatter(x, y, s=200, c=slopes, cmap=mpl.cm.jet_r, vmin=-1*np.nanmax(np.abs(slopes)), vmax=np.nanmax(np.abs(slopes)))
    
    # plt.show()
    if basemap == 1:
        plt.axis([extent[0], extent[1], extent[2], extent[3]])
    elif basemap!=1:
        plt.axis([basemap[1][0], basemap[1][1], basemap[1][2], basemap[1][3]])
    plt.savefig(base_folder+'data\\'+sitename+'\\03shoreOverlays'+eventName+'.png')
    plt.close('all')
    
    
    
    
    




def Shore_OverLays_Colorbar(base_folder,sitename,slopes,eventName=''): 

    cmapp = mpl.cm.jet_r
    normm = mpl.colors.Normalize(vmin=-31536000*np.nanmax(np.abs(slopes)), vmax=31536000*np.nanmax(np.abs(slopes)))
    
    
    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.80, 0.9, 0.1])
    cb = mpl.colorbar.ColorbarBase(ax, orientation='horizontal', cmap=cmapp,norm=normm,label='Shoreline Change m/yr')
    plt.savefig(base_folder+'data\\'+sitename+'\\03color'+eventName+'.png')
    plt.close('all')






def Make_PDF(base_folder, sitename):
    
    from PIL import Image
    
    im1 = Image.open(base_folder+'data\\'+sitename+'\\01waveClimate.png')
    im1rgb = Image.new('RGB', im1.size, (255, 255, 255))  # white background
    im1rgb.paste(im1, mask=im1.split()[3])               # paste using alpha channel as mask
    
    im2 = Image.open(base_folder+'data\\'+sitename+'\\03color.png')
    im2rgb = Image.new('RGB', im2.size, (255, 255, 255))  # white background
    im2rgb.paste(im2, mask=im2.split()[3])               # paste using alpha channel as mask
    
    im3 = Image.open(base_folder+'data\\'+sitename+'\\03shoreOverlays.png')
    im3rgb = Image.new('RGB', im3.size, (255, 255, 255))  # white background
    im3rgb.paste(im3, mask=im3.split()[3])               # paste using alpha channel as mask
    
    # im4 = Image.open(base_folder+'data\\'+sitename+'\\04bathyMap.png')
    # im4rgb = Image.new('RGB', im4.size, (255, 255, 255))  # white background
    # im4rgb.paste(im4, mask=im4.split()[3])               # paste using alpha channel as mask
    
    im5 = Image.open(base_folder+'data\\'+sitename+'\\05bathySlopeEstimates.png')
    im5rgb = Image.new('RGB', im5.size, (255, 255, 255))  # white background
    im5rgb.paste(im5, mask=im5.split()[3])               # paste using alpha channel as mask
    
    
    
    im_list = [im2rgb,im3rgb,im4rgb,im5rgb]
    im_list = [im2rgb,im3rgb,im5rgb]
    
    pdf1_filename = base_folder+'data\\'+sitename+'\\'+sitename+'_Report.pdf'
    
    im1rgb.save(pdf1_filename, "PDF" ,resolution=100.0, save_all=True, append_images=im_list)
    
    
