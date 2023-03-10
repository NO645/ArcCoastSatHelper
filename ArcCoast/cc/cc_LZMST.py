#     %matplotlib inline  
#     %matplotlib qt  
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
np.set_printoptions(linewidth=10000)
warnings.filterwarnings("ignore")


# collects the active user's username for use in creating filepaths

USERNAME = os.getlogin()



# this file was saved by the portion of the code in the arc toolbox. it is reposnsible for delivering the user's choices to the main code contained in this file
ccSettings = pickle.load(open("C:\\Users\\"+USERNAME+"\\Desktop\\coastsatv2\\temp\\ccpickle.pkl", 'rb'))

base_folder="C:\\Users\\"+USERNAME+"\\Desktop\\coastsatv2\\output"

sitename=ccSettings['sitename']
polygon=ccSettings['polygon'] # 4 points, clockwise, start from top left
Contour=float(ccSettings['Contour']) # elevation contour to push the shoreline to during tidal correction, based off of NOAA's estimate for MSL at tidal station
UserSlope = float(ccSettings['UserSlope']) # the user given beach slope
if ccSettings['UserSlope']=='na':
    UserSlope=ccSettings['UserSlope']
elif float(ccSettings['UserSlope'])<0:
    UserSlope='na'
else :
    UserSlope=float(ccSettings['UserSlope'])
    
tgage=ccSettings['tgage'] # user slected tidal gage. Usually only given on a rerun because the user needs to look at the tidal map
if ccSettings['RERUN'] == 'false': # rerun option, if true skips the download and few other steps. Also allows for user selected bad shorelines
    RERUN = False
else:
    RERUN = True
dates = [ccSettings['date1'], ccSettings['date2']] # date range
Tspace = int(ccSettings['Tspace'])

# grab the polygon shapefile and make into format acceptable by coastsat
sf = shapefile.Reader(polygon)
polygon1 = sf.shapeRecords()[0].shape.points
sf.close()
polygon2 = []
for point in polygon1:
    polygon2.append(list(point))
polygon=[polygon2]



# build inputs dictionary
inputs = {}
inputs['polygon'] = polygon
inputs['dates'] = dates
inputs['sat_list'] = ['L5', 'L8', 'S2'] # satellite missions ['L5','L7','L8','L9','S2']
inputs['sitename'] = sitename
inputs['filepath'] = os.path.join(base_folder, 'data')
inputs['landsat_collection'] = 'C02'
# inputs['include_T2'] = True # In case you need to access Tier 2 images for qualitative analysis, you need to set True



# build settings dictionary
settings = {}
settings['cloud_thresh'] = 0.5        # threshold on maximum cloud cover
settings['dist_clouds'] = 300,         # ditance around clouds where shoreline can't be mapped
settings['check_detection'] = False    # if True, shows each shoreline detection to the user for validation
settings['adjust_detection'] = False  # if True, allows user to adjust the postion of each shoreline by changing the threhold
settings['save_figure'] = True        # if True, saves a figure showing the mapped shoreline for each image
settings['inputs'] = inputs
settings['min_beach_area'] = 4500     # minimum area (in metres^2) for an object to be labelled as a beach
settings['min_length_sl'] = 200       # minimum length (in metres) of shoreline perimeter to be valid
settings['cloud_mask_issue'] = True  # switch this parameter to True if sand pixels are masked (in black) on many images  
settings['sand_color'] = 'default'   # 'default', 'dark' (for grey/black sand beaches) or 'bright' (for white sand beaches)
settings['pan_off'] = False   # True to switch pansharpening off for Landsat 7/8/9 imagery
# settings['output_epsg'] = 4326 we add this later


# OLD SETTINGS that I don't use anymore
# settings['max_dist_ref'] = 200  
# settings['along_dist'] = 25
# settings['buffer_size'] = 150         # radius (in metres) of the buffer around sandy pixels considered in the shoreline detection




# build filepath that gets used a lot later to make full path filenames
filepath = os.path.join(inputs['filepath'], sitename)


# location of image sorting algorithm training results
imageSortPkl = "C:\\Users\\"+USERNAME+"\\Desktop\\coastsatv2\\cc\\sort_classifier.pkl"


# tidal correction datum type from NOAA station
tideType = 'MSL'


# not used anymore, was once used to id buoy for runup corrections
buoy = 'na'
dissipative = False



# not used, was intended as flag for users providing their own transects
UserTransect = 'na'




# erasing the pickled file frm ARC from memory
ccSettings=None





# getting a point from the polygon file to use for some distance calculations, ie: which tidal station
long = polygon[0][0][0]
lat = polygon[0][0][1]





if not RERUN:

    SDS_download.check_images_available(inputs) # before downloading the images, check how many images are available for your inputs
    metadata = SDS_download.retrieve_images(inputs) # retrives the satellite images from Google Earth Engine
    
    pickle.dump( metadata, open( filepath+"\\metadataO.pkl", "wb" ) ) # the 'O' files are the original ones, keeping these intact helps with the reruns, they are not modified by reruns
    pickle.dump( metadata, open( filepath+"\\metadata.pkl", "wb" ) ) # the files without 'O' get modified if there is a rerun
    pickle.dump( settings, open( filepath+"\\settingsO.pkl", "wb" ) )
    
    os.chdir(base_folder)
    
    # SDS_preprocess.save_jpg(metadata, settings) # Saves .jpg files of the preprocessed satellite images (cloud masking + pansharpening/down-sampling) under ./data/sitename/jpeg_files\preprocessed
    
    
    
    # everything should be in the same UTM zone this chunk figures out what that should be from the images from GEE or looks up the location of the coordinates and the UTM boundaries
    try:
        epsg = cc_functionsv2.FindEPSG2(filepath , output)
        print('Got EPSG from image...')
    except:
        epsg = cc_functionsv2.FindEPSG(polygon)
        print('Got EPSG from polygon...')
    settings['output_epsg'] = epsg
    
    # coastsat's prompt for a reference shore. Does not run if the file already exists
    settings['reference_shoreline'] = SDS_preprocess.get_reference_sl(metadata, settings) # get reference shore from user
    settings['max_dist_ref'] = 200 # max distance (in meters) allowed from the reference shoreline

    # get all the shorelines for all the images, not the same as usual coastsat process which makes the user pick along the way. I do that later!
    output = SDS_shoreline.extract_shorelines(metadata, settings)
    
    output = SDS_tools.remove_duplicates(output) # removes duplicates (images taken on the same date by the same satellite)
    output = SDS_tools.remove_inaccurate_georef(output, 10) # remove inaccurate georeferencing (set threshold to 10 m)
    
    # save an original and modifiable version of the output dictionary
    pickle.dump( output, open( filepath+"\\outputO.pkl", "wb" ) )
    pickle.dump( output, open( filepath+"\\output.pkl", "wb" ) )
    


    
    
if RERUN:
    # since this was run before we load up the original metadata and outputs
    with open(filepath+"\\metadataO.pkl", 'rb') as f:
        metadata = pickle.load(f) 
    with open(filepath+"\\outputO.pkl", 'rb') as f:
        output = pickle.load(f) 
            
    # we reset the shorelines based on what's in the 'bad' folder at the start of the run
    output,metadata = cc_functionsv2.RemoveBadImageryLater(base_folder, sitename, output, metadata)
    pickle.dump( output, open( filepath+"\\output.pkl", "wb" ) )
    pickle.dump( metadata, open( filepath+"\\metadata.pkl", "wb" ) )
    
    # everything should be in the same UTM zone this chunk figures out what that should be from the images from GEE or looks up the location of the coordinates and the UTM boundaries
    try:
        epsg = cc_functionsv2.FindEPSG2(filepath)
        print('Got EPSG from image...')
    except:
        epsg = cc_functionsv2.FindEPSG(polygon)
        print('Got EPSG from polygon...')
    settings['output_epsg'] = epsg
    
    # coastsat's prompt for a reference shore. Does not run if the file already exists
    # get reference shore from user
    settings['reference_shoreline'] = SDS_preprocess.get_reference_sl(metadata, settings)



# save geojson of shorelines history
geomtype = 'lines' # choose 'points' or 'lines' for the layer geometry
gdf = SDS_tools.output_to_gdf(output, geomtype)
if gdf is None:
    raise Exception("output does not contain any mapped shorelines")
gdf.crs = CRS(settings['output_epsg']) # set layer projection
gdf.to_file(os.path.join(inputs['filepath'], inputs['sitename'], '%s_output.geojson'%(sitename)), driver='GeoJSON', encoding='utf-8') # save GEOJSON layer to file


# save shoreline history geojson as shapefile as well
shores_geo_filename = filepath+"\\"+sitename+"_output.geojson"
shores_shp_filename = filepath+"\\"+sitename+"_output-.shp"
cc_functionsv2.geojsonTOshp_lines(shores_geo_filename, shores_shp_filename)




if not RERUN:
    # this is the inital runs ML-issisted shoreline sorting
    cc_functionsv2.SortImagesSmall(filepath, imageSortPkl)
    # after sorting bad images are removed from the output and metadats dictionaries
    output,metadata = cc_functionsv2.RemoveBadImageryLater(base_folder, sitename, output, metadata)
    # the dictionaries are saved
    pickle.dump( output, open( filepath+"\\output.pkl", "wb" ) )
    pickle.dump( metadata, open( filepath+"\\metadata.pkl", "wb" ) )



# 1.1 Sync Output Dict and Output Shapefile to make sure the same shorelines are in each
shores_shp_filename = filepath+"\\"+sitename+"_output-.shp"
output = cc_functionsv2.SyncOutputDictShape(output, shores_shp_filename)
pickle.dump( output, open( filepath+"\\output.pkl", "wb" ) )



# 2.1. Remove Hanging Shorelines | this function attempts to trim shorelines that jump inland from pools, shiny roofs, lagoons, etc
output = cc_functionsv2.TrimOutputShorelines2(output,filepath,sitename)
# a key is added to the output dictionary for these trimmed shorelines
output['shorelines'] = copy.copy(output['TrimmedShorelines'])




# this deletes some of the intermediate output shapefilesfiles
deletes = glob.glob(filepath+"\\"+sitename+"_output-.*")
deletes = deletes + [shores_geo_filename]
for file in deletes:
    os.remove(file)


# save new geojson of trimmed shorelines history
geomtype = 'lines' # choose 'points' or 'lines' for the layer geometry
gdf = SDS_tools.output_to_gdf(output, geomtype)
if gdf is None:
    raise Exception("output does not contain any mapped shorelines")
gdf.crs = CRS(settings['output_epsg']) # set layer projection
gdf.to_file(os.path.join(inputs['filepath'], inputs['sitename'], '%s_output.geojson'%(sitename)), driver='GeoJSON', encoding='utf-8') # save GEOJSON layer to file

# save new trimmed shoreline history geojson as shapefile
shores_geo_filename = filepath+"\\"+sitename+"_output.geojson"
shores_shp_filename = filepath+"\\"+sitename+"_output-.shp"
cc_functionsv2.geojsonTOshp_lines(shores_geo_filename, shores_shp_filename)

# make reference shoreline into shapefile
shoreRef_geo_Filename = filepath+"\\"+sitename+"_reference_shoreline.geojson"
shoreRef_shp_filename = filepath+"\\"+sitename+"_reference_shoreline.shp"
cc_functionsv2.geojsonTOshp_lines(shoreRef_geo_Filename, shoreRef_shp_filename)







# 2. Build Transects
transects_filename, transectsSmall_filename, transects, transectsFull, transectsRotations, transectsSlopes = cc_functionsv2.NewTransectsBuilder_StartParallel4(Tspace, sitename, base_folder,output)
# save transect data as pickles
pickle.dump( transects, open( filepath+"\\transects.pkl", "wb" ) )
pickle.dump( transectsFull, open( filepath+"\\transectsFull.pkl", "wb" ) )
pickle.dump( transectsRotations, open( filepath+"\\transectsRotations.pkl", "wb" ) )
pickle.dump( transectsSlopes, open( filepath+"\\transectsSlopes.pkl", "wb" ) )
if len(transects)==0:
    print('There is something wrong and no transects were created.')





# Use User-given Transects
if UserTransect != 'na':
    transects = cc_functionsv2.UserTransectt(sitename, base_folder, UserTransect)
    pickle.dump( transects, open( filepath+"\\transects.pkl", "wb" ) )
    if len(transects)==0:
        print('There is something wrong and no transects were created.')







# 2.2. Shore Position Output

# get shoreline positions at transects
settings_transects = {}
settings_transects['along_dist'] = 25 #: (in metres), alongshore distance to caluclate the intersection (median of points within this distance).
# newer options not used at the moment:
# settings_transects['min_points'] #: minimum number of shoreline points to calculate an intersection.
# settings_transects['max_std'] #: (in metres) maximum STD for the shoreline points within the alongshore range, if STD is above this value a NaN is returned for this intersection.
# settings_transects['max_range'] #: (in metres) maximum RANGE for the shoreline points within the alongshore range, if RANGE is above this value a NaN is returned for this intersection.
# settings_transects['min_chainage'] #: (in metres) furthest distance landward of the transect origin that an intersection is accepted, beyond this point a NaN is returned.
# settings_transects['multiple_inter'] #: ('auto','nan','max') defines how to deal with multiple shoreline intersections
# settings_transects['auto_prc'] #: (value between 0 and 1) by default 0.1, percentage of the time that a multiple intersection needs to be present to use the max in auto mode

# calculate non-tidally corrected intersections
cross_distance = SDS_transects.compute_intersection(output, transects, settings_transects) 

# save non-tidally corrected intersections as pickle
pickle.dump( cross_distance, open( filepath+"\\cross_distance.pkl", "wb" ) )

# process non-tidally corrected intersections into a CSV file
out_dict = dict([])
out_dict['dates'] = output['dates']
for key in transects.keys():
    out_dict[key] = cross_distance[key]
    
df = pd.DataFrame(out_dict)
fn = os.path.join(settings['inputs']['filepath'],settings['inputs']['sitename'],'transect_time_series.csv')
df.to_csv(fn, sep=',')
print('Time-series of the shoreline change along the transects saved as:\n%s'%fn)





# 2.3. Get Tidal Time Series from NOAA station online

utc = output['dates'][0].tzinfo
USshape = base_folder+"\\classification\\NOAAUSAShoreline\\usaShore.shp"
gagesloc = base_folder+'\\gages.txt'
try:
    output = cc_functionsv2.Get_Tidal_DataNOAA3(base_folder+'\\', sitename,output,lat,long,utc,tgage,tideType, USshape, gagesloc)
except:
# if it doesnt work the first time it tries again 30 seconds later, sometimes NOAA site is flakey
    time.sleep(30)
    try:
        output = cc_functionsv2.Get_Tidal_DataNOAA3(base_folder+'\\', sitename,output,lat,long,utc,tgage,tideType, USshape, gagesloc)
    except:
        print('Cannot reach NOAA Tide Data!')
# saves tidal elevations in output pickle
pickle.dump( output, open( filepath+"\\output.pkl", "wb" ) )







# 3. Beach Slopes from coastsat.slope if user slope was not submitted in Arc

if UserSlope == 'na':
    try:
        days_in_year = 365.2425
        seconds_in_day = 24*3600
        settings_slope = {'slope_min':        0.035,                  # minimum slope to trial
                          'slope_max':        0.2,                    # maximum slope to trial
                          'delta_slope':      0.005,                  # slope increment
                          'date_range':       [1999,2022],            # range of dates over which to perform the analysis
                          'n_days':           8,                      # sampling period [days]
                          'n0':               50,                     # parameter for Nyquist criterium in Lomb-Scargle transforms
                          'freqs_cutoff':     1./(seconds_in_day*30), # 1 month frequency
                          'delta_f':          100*1e-10,              # deltaf for identifying peak tidal frequency band
                          'prc_conf':         0.05,                   # percentage above minimum to define confidence bands in energy curve
                          }
        
        settings_slope['date_range'] = [pytz.utc.localize(datetime.datetime(settings_slope['date_range'][0],5,1)), pytz.utc.localize(datetime.datetime(settings_slope['date_range'][1],1,1))]
        beach_slopes = SDS_slope.range_slopes(settings_slope['slope_min'], settings_slope['slope_max'], settings_slope['delta_slope'])
        
        # clip the dates between 1999 and 2020 as we need at least 2 Landsat satellites 
        idx_dates = [np.logical_and(_>settings_slope['date_range'][0],_<settings_slope['date_range'][1]) for _ in output['dates']]
        dates_sat = [output['dates'][_] for _ in np.where(idx_dates)[0]]
        tide_sat = [output['tides'][_] for _ in np.where(idx_dates)[0]]
        tide_sat = np.array(tide_sat)
        Tcross_distance = {}
        for key in cross_distance.keys():
            Tcross_distance[key] = cross_distance[key][idx_dates]
        
        # find tidal peak frequency
        settings_slope['n_days'] = 8
        settings_slope['freqs_max'] = SDS_slope.find_tide_peak(dates_sat,tide_sat,settings_slope)
        
        # estimate beach-face slopes along the transects
        slope_est, cis = dict([]), dict([])
        for key in cross_distance.keys():
            # remove NaNs
            idx_nan = np.isnan(cross_distance[key])
            Tdates = [dates_sat[_] for _ in np.where(~idx_nan)[0]]
            tide = tide_sat[~idx_nan]
            composite = cross_distance[key][~idx_nan]
            # apply tidal correction
            tsall = SDS_slope.tide_correct(composite,tide,beach_slopes)
            title = 'Transect %s'%key
            SDS_slope.plot_spectrum_all(Tdates,composite,tsall,settings_slope, title)
            slope_est[key],cis[key] = SDS_slope.integrate_power_spectrum(Tdates,tsall,settings_slope)
            print('Beach slope at transect %s: %.3f'%(key, slope_est[key]))
        TC = True
    except:
        TC = False
        print('Slope Estimation Failed!')
else:
    TC = True
    slope_est = {}
    for key in cross_distance.keys():
        slope_est[key] = UserSlope
        
        
if TC:
    # save slopes as pickle
    pickle.dump( slope_est, open( filepath+"\\slope_ests.pkl", "wb" ) )
    slope_ests_list = []
    for key in slope_est.keys():
        slope_ests_list.append(slope_est[key])
        
    # also save slopes as txt file
    with open(filepath+"\\slope_ests.txt", 'w') as fp:
        yyy = -1
        for item in slope_ests_list:
            yyy = yyy + 1
            fp.write(str(list(slope_est.keys())[yyy]) + ' : ' + "%s\n" % item)



    
# 4. Tidal Correction of Shore Position Output Table
# IMPORTANT: the mean of all the slopes along the region are used, not the slope corresponding to the transect along which the shore is moved
if TC:
    beach_slope = np.nanmean(slope_ests_list)
    reference_elevation = 0 # elevation at which you would like the shoreline time-series to be
    cross_distance_tidally_corrected = {}
    for key in cross_distance.keys():
        correction = (np.array(output['tides'])-reference_elevation)/beach_slope
        cross_distance_tidally_corrected[key] = cross_distance[key] + correction
        
    
    # store the tidally-corrected intersection series time-series in a .csv file
    out_dict = dict([])
    out_dict['dates'] = copy.copy(output['dates'])
    for key in cross_distance_tidally_corrected.keys():
        out_dict['Transect '+ key] = cross_distance_tidally_corrected[key]
    df = pd.DataFrame(out_dict)
    fn = os.path.join(settings['inputs']['filepath'], settings['inputs']['sitename'], 'transect_time_series_T_corrected.csv')
    df.to_csv(fn, sep=',')
    print('Tidally-corrected time-series of the shoreline change along the transects saved as:\n%s'%fn)
 
 
 
 
 
output = cc_functionsv2.GetProjectionInfo(sitename, base_folder+'\\', output)


# 4. Tidal Correction of Shore Position Shapefile
# IMPORTANT: the shores are moved perpendicular to the slope of the entire AOI, not to the small area around the point getting moved

if TC:
    output = cc_functionsv2.AdjustOutputShapefileByTides5(slope_est,output,filepath,sitename,Contour,transects)
    # output = cc_functions.AdjustOutputShapefileByWaves1(slope_est,output,filepath,sitename,buoy,dissipative,transects)
    # if you ever have to do runup split apart the horizontal corrections and the actual shapefile change
    # instead get the tide correction, then get the runup correction, add them, and then adjust shapefile
else:
    pass
    
    
    
try:
    output = cc_functionsv2.SyncOutputDictShape(output,  filepath+"\\"+sitename+"_outputTidalC.shp")
    pickle.dump( output, open( filepath+"\\output.pkl", "wb" ) )
except:
    print('Sync Failed')
    pickle.dump( output, open( filepath+"\\output.pkl", "wb" ) )

    
    