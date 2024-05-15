#%%
from pathlib import Path
import os
from astropy import units
from astropy.nddata import CCDData
from astropy.stats import mad_std
import ccdproc as ccdp
import numpy as np

FOLDER          = 'LOT20220910'
CALIBRATED_DIR  = FOLDER + '/calibrated_data'
BIAS_DIR        = FOLDER + '/bias'
DARK_DIR        = FOLDER + '/dark'
FLAT_DIR        = FOLDER + '/flat'
DARK_FLAT_DIR   = FOLDER + '/dark flat'
LIGHT_DIR       = FOLDER + '/cngeow'
RDCED_FLATS_DIR = CALIBRATED_DIR + '/reduced_flats'
REDUCED_LIGHT   = CALIBRATED_DIR + '/reduced_lights'

MEM_LIM = 10e9


Path('.', CALIBRATED_DIR).mkdir(exist_ok = True)
calibrated_path = Path(CALIBRATED_DIR)

Path('.', REDUCED_LIGHT).mkdir(exist_ok = True)
lights_path = Path(REDUCED_LIGHT)


# Bias
bias_files = ccdp.ImageFileCollection(Path(BIAS_DIR))

bias_list = bias_files.files_filtered(imagetyp = 'bias',
                                      include_path = True)

combined_bias = ccdp.combine(bias_list,
                             method = 'average',
                             sigma_clip = True,
                             sigma_clip_low_thresh = 5,
                             sigma_clip_high_thresh = 5,
                             sigma_clip_func = np.ma.median,
                             sigma_clip_dev_func = mad_std,
                             mem_limit = MEM_LIM,
                             unit = 'adu')
                             
combined_bias.meta['combined'] = True
combined_bias.write(calibrated_path / 'combined_bias.fit', overwrite = True)


# Dark
dark_files = ccdp.ImageFileCollection(Path(DARK_DIR))
dark_exp_time = set(dark_files.summary['exptime'])

for exp_time in dark_exp_time:
    darks_list = dark_files.files_filtered(imagetyp = 'dark',
                                           exptime = exp_time,
                                           include_path = True)

    combined_dark = ccdp.combine(darks_list,
                                 method = 'average',
                                 sigma_clip = True,
                                 sigma_clip_low_thresh = 5,
                                 sigma_clip_high_thresh = 5,
                                 sigma_clip_func = np.ma.median,
                                 signma_clip_dev_func = mad_std,
                                 mem_limit = MEM_LIM,
                                 unit = 'adu')

    combined_dark.meta['combined'] = True               
    combined_dark.write(calibrated_path / 'combined_dark_{}.fit'.format(exp_time), overwrite = True)


# Flat
if os.path.isdir(FLAT_DIR) is True:
    Path('.', RDCED_FLATS_DIR).mkdir(exist_ok = True)
    flats_path = Path(RDCED_FLATS_DIR)

    flat_files = ccdp.ImageFileCollection(Path(FLAT_DIR))
    a_flat = CCDData.read(flat_files.files_filtered(imagetyp='flat', include_path=True)[0], unit='adu')

    flat_filter = set(flat_files.summary['filter'])
    flat_exp_time = set(flat_files.summary['exptime'])

    def inv_median(a):
        return 1 / np.median(a)

    for flt in flat_filter:
        flat_list = flat_files.files_filtered(imagetyp = 'flat',
                                            filter = flt,
                                            include_path = True)
        combined_flat = ccdp.combine(flat_list,
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
        combined_flat.write(calibrated_path / 'combined_flat_{}_{}.fit'.format(combined_flat.header['exptime'], flt), overwrite = True)


    # Dark Flat
    dark_flat_files = ccdp.ImageFileCollection(Path(DARK_FLAT_DIR))

    dark_falt_exp_time = set(dark_flat_files.summary['exptime'])

    for exp_time in dark_falt_exp_time:
        darks_flat_list = dark_flat_files.files_filtered(imagetyp = 'dark',
                                                        exptime = exp_time,
                                                        include_path = True)

        combined_dark_flat = ccdp.combine(darks_flat_list,
                                        method = 'average',
                                        sigma_clip = True,
                                        sigma_clip_low_thresh = 5,
                                        sigma_clip_high_thresh = 5,
                                        sigma_clip_func = np.ma.median,
                                        signma_clip_dev_func = mad_std,
                                        mem_limit = MEM_LIM,
                                        unit = 'adu')
                                        
        combined_dark_flat.meta['combined'] = True
        combined_dark_flat.write(calibrated_path / 'combined_dark_flat_{}.fit'.format(exp_time), overwrite = True)


    #Flat Reduction
    def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=60):

        dark_exposures = np.array(list(dark_exposure_times))
        idx = np.argmin(np.abs(dark_exposures - image.header['exptime']))
        closest_dark_exposure = dark_exposures[idx]

        if (tolerance is not None and 
            np.abs(image.header['exptime'] - closest_dark_exposure) > tolerance):
            
            raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                            'time {}.'.format(closest_dark_exposure, image.header['exptime']))
            
        
        return closest_dark_exposure

    calibrate_files = ccdp.ImageFileCollection(calibrated_path)

    combined_darks = {ccd.header['exptime']: ccd for ccd in calibrate_files.ccds(imagetyp='dark', combined=True)}
    combined_flats = {ccd.header['filter']: ccd for ccd in calibrate_files.ccds(imagetyp='flat', combined=True)}
    combined_dark_flats = {ccd.header['exptime']: ccd for ccd in calibrate_files.ccds(imagetyp='dark', combined=True)}

    for flats, flats_fname in calibrate_files.ccds(imagetyp='flat',
                                                return_fname=True,
                                                ccd_kwargs=dict(unit='adu')):

        reduced_flat = ccdp.subtract_bias(flats, combined_bias)

        closest_dark = find_nearest_dark_exposure(reduced_flat, combined_darks.keys())
        reduced_flat = ccdp.subtract_dark(reduced_flat,
                                        combined_darks[closest_dark], 
                                        exposure_time='exptime',
                                        exposure_unit=units.second,
                                        scale=True)

        reduced_flat.write(flats_path / 'reducted_flat_{}.fits'.format(reduced_flat.header['filter']), overwrite = True)

else:
    pass


# Light (Single)
light_files = ccdp.ImageFileCollection(Path(LIGHT_DIR))
light_set = set(light_files.files_filtered(imagetyp = 'light', include_path = True))
light_name = set(light_files.summary['file'])
light_exp_time = set(light_files.summary['exptime'])
light_filter = set(light_files.summary['filter'])


# Reduction
def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=0.5):

    dark_exposures = np.array(list(dark_exposure_times))
    idx = np.argmin(np.abs(dark_exposures - image.header['exptime']))
    closest_dark_exposure = dark_exposures[idx]

    if (tolerance is not None and 
        np.abs(image.header['exptime'] - closest_dark_exposure) > tolerance):
        
        raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           'time {}.'.format(closest_dark_exposure, a_flat.header['exptime']))
        
    
    return closest_dark_exposure

for light, file_name in light_files.ccds(imagetyp='light',
                                         return_fname=True,
                                         ccd_kwargs=dict(unit='adu')):

    reduced = ccdp.subtract_bias(light, combined_bias)

    closest_dark = find_nearest_dark_exposure(reduced, combined_darks.keys())
    reduced = ccdp.subtract_dark(reduced,
                                 combined_darks[closest_dark],
                                 exposure_time='exptime',
                                 exposure_unit=units.second)
    
    if os.path.isdir(FLAT_DIR) is True:
        reduced_flat_files = ccdp.ImageFileCollection(flats_path)
        reduced_flats = {ccd.header['filter']: ccd for ccd in reduced_flat_files.ccds(imagetyp='flat', combined=True)}
        
        good_flat = reduced_flats[reduced.header['filter']]
        reduced = ccdp.flat_correct(reduced, good_flat)
    else:
        pass
    
    reduced.write(lights_path / file_name, overwrite = True)


#Info
print('light_exp:    ', light_exp_time)
print('light_filter: ', light_filter)
print('dark_exp:     ', dark_exp_time)
print('flat_filter:  ', flat_filter)
print('flat_exp:     ', flat_exp_time)
print('dark_flat_exp:', dark_falt_exp_time)
# %%
