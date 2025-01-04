#%%
from pathlib import Path
import os
from astropy import units
from astropy.nddata import CCDData
from astropy.stats import mad_std
from ccdproc import ImageFileCollection, combine, subtract_bias, subtract_dark, flat_correct
import numpy as np

# PATH INFO
WORKING_DIR   = os.path.join(os.getcwd(), 'test', 'LOT20220724')

BIAS_DIR      = Path(WORKING_DIR, 'bias')
DARK_DIR      = Path(WORKING_DIR, 'dark')
FLAT_DIR      = Path(WORKING_DIR, 'flat')
DARK_FLAT_DIR = Path(WORKING_DIR, 'dark flat')
LIGHT_DIR     = Path(WORKING_DIR, 'cngeow')

MEM_LIM = 10e9  # bytes


def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=0.5):

    dark_exposures = np.array(list(dark_exposure_times))
    idx = np.argmin(np.abs(dark_exposures - image.header['exptime']))
    closest_dark_exposure = dark_exposures[idx]

    if (tolerance is not None and 
        np.abs(image.header['exptime'] - closest_dark_exposure) > tolerance):
        
        raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           'time {}.'.format(closest_dark_exposure, a_flat.header['exptime']))    
    
    return closest_dark_exposure


print('CURRENT WORKING DIRECTORY:', WORKING_DIR)

Path(WORKING_DIR, 'calibrated_data').mkdir(exist_ok=True)
calibrated_data_path = Path(WORKING_DIR, 'calibrated_data')

# BIAS COMBINATION
if os.path.isdir(BIAS_DIR) is True:
    bias_files = ImageFileCollection(BIAS_DIR)

    bias_path_list = bias_files.files_filtered(imagetyp = 'bias',
                                            include_path = True)

    combined_bias = combine(bias_path_list,
                            method = 'average',
                            sigma_clip = True,
                            sigma_clip_low_thresh = 5,
                            sigma_clip_high_thresh = 5,
                            sigma_clip_func = np.ma.median,
                            sigma_clip_dev_func = mad_std,
                            mem_limit = MEM_LIM,
                            unit = 'adu')
                                
    combined_bias.meta['combined'] = True
    combined_bias.write(calibrated_data_path / 'combined_bias.fits', overwrite = True)

else:
    print('\033[1m\033[91mBIAS NOT FOUND\033[0m')
    pass


# DARK CALIBRATION & COMBINATION
if os.path.isdir(FLAT_DIR) is True:
    dark_files = ImageFileCollection(DARK_DIR)
    dark_exptime_set = set(dark_files.summary['exptime'])

    for exptime in dark_exptime_set:
        
        calibrated_dark_list = []
        for dark_ccd, dark_name in dark_files.ccds(imagetyp = 'dark',
                                                   return_fname = True,
                                                   exptime = exptime,
                                                   ccd_kwargs = dict(unit='adu')):
            
            calibrated_dark_list.append(subtract_bias(dark_ccd, combined_bias))

        combined_dark = combine(calibrated_dark_list,
                                method = 'average',
                                sigma_clip = True,
                                sigma_clip_low_thresh = 5,
                                sigma_clip_high_thresh = 5,
                                sigma_clip_func = np.ma.median,
                                signma_clip_dev_func = mad_std,
                                mem_limit = MEM_LIM,
                                unit = 'adu')

        combined_dark.meta['combined'] = True               
        combined_dark.write(calibrated_data_path / f'combined_dark_{exptime}.fits', overwrite = True)

else:
    print('\033[1m\033[91mDARK NOT FOUND\033[0m')
    pass


# DARK FLAT CALIBRATION & COMBINATION
if os.path.isdir(DARK_FLAT_DIR) is True:
    dark_flat_files = ImageFileCollection(DARK_FLAT_DIR)
    dark_falt_exptime_set = set(dark_flat_files.summary['exptime'])

    for exptime in dark_falt_exptime_set:

        calibrated_dark_flat = []
        for dark_flat_ccd, dark_flat_name in dark_flat_files.ccds(imagetyp = 'dark',
                                                                  return_fname = True,
                                                                  exptime = exptime,
                                                                  ccd_kwargs = dict(unit='adu')):
            
            calibrated_dark_flat.append(subtract_bias(dark_flat_ccd, combined_bias))

    combined_dark_flat = combine(calibrated_dark_flat,
                                 method = 'average',
                                 sigma_clip = True,
                                 sigma_clip_low_thresh = 5,
                                 sigma_clip_high_thresh = 5,
                                 sigma_clip_func = np.ma.median,
                                 signma_clip_dev_func = mad_std,
                                 mem_limit = MEM_LIM,
                                 unit = 'adu')
                                    
    combined_dark_flat.meta['combined'] = True
    combined_dark_flat.write(calibrated_data_path / f'combined_dark_flat_{exptime}.fits', overwrite = True)

else:
    print('\033[1m\033[91mDARK FLAT NOT FOUND\033[0m')
    pass


# FLAT CALIBRATION & COMBINATION
if os.path.isdir(FLAT_DIR) is True:
    flat_files = ImageFileCollection(FLAT_DIR)
    a_flat = CCDData.read(flat_files.files_filtered(imagetyp = 'flat', include_path = True)[0], unit = 'adu')

    def inv_median(a):
        return 1 / np.median(a)

    calibrated_files = ImageFileCollection(calibrated_data_path)
    combined_dark_flat_dict = {ccd.header['exptime']: ccd for ccd in calibrated_files.ccds(imagetyp='dark', combined=True)}

    flat_filter_set = set(flat_files.summary['filter'])

    for filter in flat_filter_set:

        calibrated_flat_list = []
        for flat_ccd, flat_name in flat_files.ccds(imagetyp = 'flat',
                                                   return_fname = True,
                                                   filter = filter,
                                                   ccd_kwargs = dict(unit='adu')):
        
            calibrated_flat = subtract_bias(flat_ccd, combined_bias)

            closest_dark = find_nearest_dark_exposure(calibrated_flat, combined_dark_flat_dict.keys())

            calibrated_flat_list.append(subtract_dark(calibrated_flat,
                                                      combined_dark_flat_dict[closest_dark], 
                                                      exposure_time = 'exptime',
                                                      exposure_unit = units.second,
                                                      scale = True))

        combined_flat = combine(calibrated_flat_list,
                                method = 'average',
                                scale = inv_median,
                                sigma_clip = True,
                                sigma_clip_low_thresh = 5,
                                sigma_clip_high_thresh = 5,
                                sigma_clip_func = np.ma.median,
                                signma_clip_dev_func = mad_std,
                                mem_limit = MEM_LIM,
                                unit = 'adu')
                                    
        combined_flat.meta['combined'] = True
        combined_flat.write(calibrated_data_path / f'combined_flat_{filter}.fits', overwrite = True)

else:
    print('\033[1m\033[91mFLAT NOT FOUND\033[0m')
    pass


# LIGHT REDUCTION
Path(WORKING_DIR, 'calibrated_data', 'reduced_light_lights').mkdir(exist_ok=True)
calibrated_lights_path = Path(calibrated_data_path, 'reduced_light_lights')

light_files = ImageFileCollection(LIGHT_DIR)

combined_dark_dict = {ccd.header['exptime']: ccd for ccd in calibrated_files.ccds(imagetyp='dark', combined=True)}

for light_ccd, light_name in light_files.ccds(imagetyp='light',
                                              return_fname=True,
                                              ccd_kwargs=dict(unit='adu')):

    if os.path.isdir(BIAS_DIR) is True:
        reduced_light = subtract_bias(light_ccd, combined_bias)
    else:
        pass

    if os.path.isdir(DARK_DIR) is True:
        closest_dark = find_nearest_dark_exposure(reduced_light, combined_dark_dict.keys())

        reduced_light = subtract_dark(reduced_light,
                                      combined_dark_dict[closest_dark],
                                      exposure_time = 'exptime',
                                      exposure_unit = units.second,
                                      scale = False)
    else:
        pass
    
    if os.path.isdir(FLAT_DIR) is True:
        calibrated_flat_dict = {ccd.header['filter']: ccd for ccd in calibrated_files.ccds(imagetyp='flat', combined=True)}
        print(calibrated_flat_dict.keys())
        print(reduced_light.header['filter'])
        combined_flat = calibrated_flat_dict[reduced_light.header['filter']]
        reduced_light = flat_correct(reduced_light, combined_flat)
    else:
        pass
    
    reduced_light.write(calibrated_lights_path / f'{light_name}_calibrated', overwrite = True)

print(f'{WORKING_DIR} \033[1m\033[32mDONE REDUCTION\033[0m')

# LOG
light_exp_time   = set(light_files.summary['exptime'])
light_filter     = set(light_files.summary['filter'])
flat_exptime_set = set(flat_files.summary['exptime'])

print('\033[1mLOG:\033[0m')
print('light_exp:    ', light_exp_time)
print('light_filter: ', light_filter)
print('dark_exp:     ', dark_exptime_set)
print('flat_filter:  ', flat_filter_set)
print('flat_exp:     ', flat_exptime_set)
print('dark_flat_exp:', dark_falt_exptime_set)