#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Version: v1.2
Date: 2021-04-01
Authors: Mullissa A., Vollrath A., Braun, C., Slagter B., Balling J., Gou Y., Gorelick N.,  Reiche J.
Description: A wrapper function to derive the Sentinel-1 ARD
"""

import ee
import geemap
import os
import border_noise_correction as bnc
import speckle_filter as sf
import terrain_flattening as trf
import helper

# VPN port
geemap.set_proxy(port=5188)
try:
    ee.Initialize()
except:
    ee.Authenticate()
    ee.Initialize()


###########################################
# DO THE JOB
###########################################

def s1_preproc(params):
    """
    Applies preprocessing to a collection of S1 images to return an analysis ready sentinel-1 data.

    Parameters
    ----------
    params : Dictionary
        These parameters determine the data selection and image processing parameters.

    Raises
    ------
    ValueError


    Returns
    -------
    ee.ImageCollection
        A processed Sentinel-1 image collection

    """

    APPLY_BORDER_NOISE_CORRECTION = params['APPLY_BORDER_NOISE_CORRECTION']
    APPLY_TERRAIN_FLATTENING = params['APPLY_TERRAIN_FLATTENING']
    APPLY_SPECKLE_FILTERING = params['APPLY_SPECKLE_FILTERING']
    POLARIZATION = params['POLARIZATION']
    PLATFORM_NUMBER = params['PLATFORM_NUMBER']
    ORBIT = params['ORBIT']
    ORBIT_NUM = params['ORBIT_NUM']
    SPECKLE_FILTER_FRAMEWORK = params['SPECKLE_FILTER_FRAMEWORK']
    SPECKLE_FILTER = params['SPECKLE_FILTER']
    SPECKLE_FILTER_KERNEL_SIZE = params['SPECKLE_FILTER_KERNEL_SIZE']
    SPECKLE_FILTER_NR_OF_IMAGES = params['SPECKLE_FILTER_NR_OF_IMAGES']
    TERRAIN_FLATTENING_MODEL = params['TERRAIN_FLATTENING_MODEL']
    DEM = params['DEM']
    TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = params['TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER']
    FORMAT = params['FORMAT']
    START_DATE = params['START_DATE']
    STOP_DATE = params['STOP_DATE']
    ROI = params['ROI']
    CLIP_TO_ROI = params['CLIP_TO_ROI']
    SAVE_ASSET = params['SAVE_ASSET']
    ASSET_ID = params['ASSET_ID']
    SAVE_LOCAL = params['SAVE_LOCAL']
    EXPORT_CRS = params['EXPORT_CRS']
    EXPORT_SCALE = params['EXPORT_SCALE']
    VISUALIZATION = params['VISUALIZATION']
    RENDER_SCALE = params['RENDER_SCALE']
    LOCAL_DIR = params['LOCAL_DIR']

    ###########################################
    # 0. CHECK PARAMETERS
    ###########################################

    if APPLY_BORDER_NOISE_CORRECTION is None:
        APPLY_BORDER_NOISE_CORRECTION = True
    if APPLY_TERRAIN_FLATTENING is None:
        APPLY_TERRAIN_FLATTENING = True
    if APPLY_SPECKLE_FILTERING is None:
        APPLY_SPECKLE_FILTERING = True
    if POLARIZATION is None:
        POLARIZATION = 'VVVH'
    if ORBIT is None:
        ORBIT = 'BOTH'
    if SPECKLE_FILTER_FRAMEWORK is None:
        SPECKLE_FILTER_FRAMEWORK = 'MULTI BOXCAR'
    if SPECKLE_FILTER is None:
        SPECKLE_FILTER = 'GAMMA MAP'
    if SPECKLE_FILTER_KERNEL_SIZE is None:
        SPECKLE_FILTER_KERNEL_SIZE = 7
    if SPECKLE_FILTER_NR_OF_IMAGES is None:
        SPECKLE_FILTER_NR_OF_IMAGES = 10
    if TERRAIN_FLATTENING_MODEL is None:
        TERRAIN_FLATTENING_MODEL = 'VOLUME'
    if TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER is None:
        TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER = 0
    if FORMAT is None:
        FORMAT = 'DB'
    if ORBIT is None:
        ORBIT = 'DESCENDING'
    if EXPORT_CRS is None:
        EXPORT_CRS = 'EPSG:4326'
    if EXPORT_SCALE is None:
        EXPORT_SCALE = 10
    if RENDER_SCALE is None:
        RENDER_SCALE = 100

    pol_required = ['VV', 'VH', 'VVVH', 'HH', 'HV', 'HHHV']
    if POLARIZATION not in pol_required:
        raise ValueError("ERROR!!! Parameter POLARIZATION not correctly defined")

    orbit_required = ['ASCENDING', 'DESCENDING', 'BOTH']
    if ORBIT not in orbit_required:
        raise ValueError("ERROR!!! Parameter ORBIT not correctly defined")

    model_required = ['DIRECT', 'VOLUME']
    if TERRAIN_FLATTENING_MODEL not in model_required:
        raise ValueError("ERROR!!! Parameter TERRAIN_FLATTENING_MODEL not correctly defined")

    format_required = ['LINEAR', 'DB']
    if FORMAT not in format_required:
        raise ValueError("ERROR!!! FORMAT not correctly defined")

    frame_needed = ['MONO', 'MULTI']
    if SPECKLE_FILTER_FRAMEWORK not in frame_needed:
        raise ValueError("ERROR!!! SPECKLE_FILTER_FRAMEWORK not correctly defined")

    format_sfilter = ['BOXCAR', 'LEE', 'GAMMA MAP'
        , 'REFINED LEE', 'LEE SIGMA']
    if SPECKLE_FILTER not in format_sfilter:
        raise ValueError("ERROR!!! SPECKLE_FILTER not correctly defined")

    if TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER < 0:
        raise ValueError("ERROR!!! TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER not correctly defined")

    if SPECKLE_FILTER_KERNEL_SIZE <= 0:
        raise ValueError("ERROR!!! SPECKLE_FILTER_KERNEL_SIZE not correctly defined")

    ###########################################
    # 1. DATA SELECTION
    ###########################################

    # select S-1 image collection
    s1 = ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT') \
        .filter(ee.Filter.eq('instrumentMode', 'IW')) \
        .filter(ee.Filter.eq('resolution_meters', 10)) \
        .filterDate(START_DATE, STOP_DATE) \
        .filterBounds(ROI)

    if PLATFORM_NUMBER == 'A' or PLATFORM_NUMBER == 'B':
        s1 = s1.filter(ee.Filter.eq('platform_number', PLATFORM_NUMBER))

    if ORBIT_NUM is not None:
        s1 = s1.filter(ee.Filter.eq('relativeOrbitNumber_start', ORBIT_NUM))
        # .filter(ee.Filter.eq('relativeOrbitNumber_start',None))

    # select orbit
    if ORBIT != 'BOTH':
        s1 = s1.filter(ee.Filter.eq('orbitProperties_pass', ORBIT))

    # select polarization
    polar_bands = []
    if POLARIZATION == 'VV':
        polar_bands = ['VV', 'angle']
    elif POLARIZATION == 'VH':
        polar_bands = ['VH', 'angle']
    elif POLARIZATION == 'VVVH':
        polar_bands = ['VV', 'VH', 'angle']
    elif POLARIZATION == 'HH':
        polar_bands = ['HH', 'angle']
    elif POLARIZATION == 'HV':
        polar_bands = ['HV', 'angle']
    elif POLARIZATION == 'HHHV':
        polar_bands = ['HHHV', 'angle']

    s1 = s1.select(polar_bands)
    print('Number of images in collection: ', s1.size().getInfo())

    # Get ImageCollection footprint list
    sizeRaw = s1.size().getInfo()
    imlistRaw = s1.toList(sizeRaw)
    footprintList = []
    for idx in range(0, sizeRaw):
        img = imlistRaw.get(idx)
        img = ee.Image(img)
        footprint = ee.Geometry(img.get('system:footprint'))
        footprintList.append(footprint)

    ###########################################
    # 2. ADDITIONAL BORDER NOISE CORRECTION
    ###########################################

    if APPLY_BORDER_NOISE_CORRECTION:
        s1_1 = s1.map(bnc.f_mask_edges)
        print('Additional border noise correction is completed')
    else:
        s1_1 = s1
    ########################
    # 3. SPECKLE FILTERING
    #######################

    if APPLY_SPECKLE_FILTERING:
        if SPECKLE_FILTER_FRAMEWORK == 'MONO':
            s1_1 = ee.ImageCollection(sf.MonoTemporal_Filter(s1_1, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER))
            print('Mono-temporal speckle filtering is completed')
        else:
            s1_1 = ee.ImageCollection(
                sf.MultiTemporal_Filter(s1_1, SPECKLE_FILTER_KERNEL_SIZE, SPECKLE_FILTER, SPECKLE_FILTER_NR_OF_IMAGES))
            print('Multi-temporal speckle filtering is completed')

    ########################
    # 4. TERRAIN CORRECTION
    #######################

    if APPLY_TERRAIN_FLATTENING:
        s1_1 = (trf.slope_correction(s1_1
                                     , TERRAIN_FLATTENING_MODEL
                                     , DEM
                                     , TERRAIN_FLATTENING_ADDITIONAL_LAYOVER_SHADOW_BUFFER))
        print('Radiometric terrain normalization is completed')

    ########################
    # 5. OUTPUT
    #######################

    if FORMAT == 'DB':
        s1_1 = s1_1.map(helper.lin_to_db)

    # clip to roi
    if CLIP_TO_ROI:
        s1_1 = s1_1.map(lambda image: image.clip(ROI))

    # save to asset
    if SAVE_ASSET:

        size = s1_1.size().getInfo()
        imlist = s1_1.toList(size)
        for idx in range(0, size):
            img = imlist.get(idx)
            img = ee.Image(img)
            name = str(img.id().getInfo())
            description = name
            assetId = ASSET_ID + '/' + name

            task = ee.batch.Export.image.toAsset(image=img,
                                                 assetId=assetId,
                                                 description=description,
                                                 region=s1_1.geometry(),
                                                 scale=EXPORT_SCALE,
                                                 maxPixels=1e13)
            task.start()
            print('Exporting {} to {}'.format(name, assetId))

    # save to local
    if SAVE_LOCAL:

        size = s1_1.size().getInfo()
        imlist = s1_1.toList(size)
        for idx in range(0, size):
            img = imlist.get(idx)
            img = ee.Image(img)
            name = str(img.id().getInfo())

            # save raw images to local
            if not os.path.exists(LOCAL_DIR):
                os.makedirs(LOCAL_DIR)
            filenameRaw = os.path.join(LOCAL_DIR, name + '.tif')
            print('Downloading Raw Image: {} to {}'.format(name, filenameRaw))
            if CLIP_TO_ROI:
                geemap.download_ee_image(img,
                                         filenameRaw,
                                         region=ROI,
                                         crs=EXPORT_CRS,
                                         scale=EXPORT_SCALE)
            else:
                geemap.download_ee_image(img,
                                         filenameRaw,
                                         region=footprintList[idx],
                                         crs=EXPORT_CRS,
                                         scale=EXPORT_SCALE)

            # save visualization images to local
            if VISUALIZATION:
                if len(POLARIZATION) > 2:
                    # raise ValueError("ERROR!!! Only can convert single band image into an 32-int gray image")
                    img = img.select(polar_bands)
                    visImage = img.visualize(**{
                        'bands': polar_bands,
                        'min': [-25, -25, 30],
                        'max': [0, 0, 45],
                        'gamma': 1.0
                    })
                else:
                    img = img.select([POLARIZATION])
                    visImage = img.visualize(**{
                        'bands': [POLARIZATION],
                        'min': -25,
                        'max': 0.0,
                        'gamma': 1.0
                    })
                filename_render = os.path.join(LOCAL_DIR, name + '_render.tif')
                print('Downloading Visualization Image to {}'.format(filename_render))
                if CLIP_TO_ROI:
                    # geemap.download_ee_image(visImage, filename_render, region=ROI, crs=EXPORT_CRS, scale=RENDER_SCALE, dtype='int32')
                    geemap.ee_export_image(visImage, filename=filename_render, crs=EXPORT_CRS, scale=RENDER_SCALE, region=ROI)
                else:
                    # geemap.download_ee_image(visImage, filename_render, region=footprintList[idx], crs=EXPORT_CRS, scale=RENDER_SCALE, dtype='int32')
                    geemap.ee_export_image(visImage, filename=filename_render, crs=EXPORT_CRS, scale=RENDER_SCALE, region=footprintList[idx])

    return s1_1
