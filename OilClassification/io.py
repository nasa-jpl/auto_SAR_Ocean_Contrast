import numpy as np
from os import path, makedirs
from os.path import abspath, basename
from netCDF4 import Dataset
from osgeo import gdal
import spectral

__all__ = ['get_most_recent_git_tag',
           'parse_prodname_UAVSAR',
           'writedata']

import subprocess
def get_most_recent_git_tag():
    try:
        git_tag = str(
            subprocess.check_output(['git', 'describe', '--abbrev=0', '--tags'], stderr=subprocess.STDOUT)
        ).strip('\'b\\n')
    except subprocess.CalledProcessError as exc_info:
        raise Exception(str(exc_info.output))
    return git_tag

def parse_prodname_UAVSAR(filename):
    """parse_prodname_UAVSAR(filename)

    Given a filename using the UAVSAR product naming convention, returns different
    parameters if found.  If not found returns [] for that field.  As written, this
    assumes products are for 90Deg beam steering.

    The UAVSAR filename convention is:
    {site}_{lineID}_{fltID}_{datatake}_{fltdate}_L090([HV]{4}|.?*)_XX_{int}.*
    """
    ProdPolarizations = ['HHHH', 'HVHV', 'VVVV', 'HHVV', 'HVVV', 'HHHV']

    # Get the UAVSAR product information (dataName and product name)

    #prod = [join(dirname(fn), basename(fn).split('.')[0]) for fn in filename] # full path to first '.'
    filename = basename(filename).split('.')[0]

    assert('L090' in filename) # Assumes products with 90 deg beam steering
    pol = filename.split('L090')[1][0:4] # first 4 chars after 'L090'
    if pol in ProdPolarizations:
        dataName = filename.replace(thispol,'')
    else:
        dataName = filename
        pol      = None

    _, lineID, fltID, datatake, fltdate, *rest = dataName.split('_')

    return [dataName, filename, pol, lineID, fltID, datatake, fltdate]

def writedata(datadict, infile, outdir, cfg):
    """Write data from processing infile to outfile."""
    (infile_basename, suffix) = path.splitext(basename(infile))
    ID     = cfg['ID']
    format = cfg['format'].upper()

    if format == 'NETCDF':
        outfile = abspath(path.join(outdir, f'{infile_basename}_{ID}_contrast_ratio.nc'))
        with Dataset(infile,'r') as ds, Dataset(outfile,'w') as of:
            # global attributes
            of.setncatts(ds.__dict__)
            for attr in cfg['rm_attr']:
                try:
                    of.delncattr(attr)
                except:
                    print(f'attribute {attr}: nothing to remove')
            for attr in cfg['global_attr']:
                of.setncattr(attr, cfg['global_attr'][attr])
            # dimensions
            for dim in ds.dimensions.values():
                size = len(dim) * (not dim.isunlimited())
                of.createDimension(dim.name, size)
            # variables
            for var in cfg['vars']:
                if var in ds.variables:
                    of.createVariable(var,ds[var].dtype,ds[var].dimensions,'zlib')
                    of[var].setncatts(ds[var].__dict__)
                    of[var][:] = ds[var][:]
                else:
                    print(f'variable {var} not found in input file')
            # main results: dampratio
            for key in datadict:
                of.createVariable(key,'f4',('y','x'),'zlib',fill_value=np.nan)
                of[key].setncatts(ds['sigma'].__dict__)
                of[key][:] = datadict[key]['data']
                of[key].long_name = datadict[key]['long_name']


            ### NOAA special request: x,y,nx,ny,lat,lon etc. ###
            of.createVariable('northernmost_latitude','f4')
            of['northernmost_latitude'][:] = np.amax(ds['latitude'])
            of.createVariable('southernmost_latitude','f4')
            of['southernmost_latitude'][:] = np.amin(ds['latitude'])
            of.createVariable('easternmost_longitude','f4')
            of['easternmost_longitude'][:] = np.amax(ds['longitude'])
            of.createVariable('westernmost_longitude','f4')
            of['westernmost_longitude'][:] = np.amin(ds['longitude'])

            ny, nx = ds['sigma'].shape

            of.createVariable('nx','i4',fill_value=1-2**31)
            of['nx'].setncatts({'units': 'Dimensionless', 'long_name': 'Image dimensions'})
            of['nx'][:] = nx

            of.createVariable('ny','i4',fill_value=1-2**31)
            of['ny'].setncatts({'units': 'Dimensionless', 'long_name': 'Image dimensions'})
            of['ny'][:] = ny

            if 'x' not in of.variables:
                of.createVariable('x','i4',('x',),fill_value=1-2**31)
                of['x'].setncatts({'units': '1',
                                   'standard_name': 'projection_x_coordinate',
                                   'long_name': 'Row Index of Pixel Centers',
                                   'axis': 'X'})

            if 'y' not in of.variables:
                of.createVariable('y','i4',('y',),fill_value=1-2**31)
                of['y'].setncatts({'units': '1',
                                   'standard_name': 'projection_y_coordinate',
                                   'long_name': 'Row Index of Pixel Centers',
                                   'axis': 'Y'})

            # arrays for which we do not have an input or do not generate and output.
            garbage_long_name = {'sigma2'          : 'Normalized Radar Cross Section 2nd array',
                                 'coil'            : 'Surfactant classification',
                                 'clowwind'        : 'Low wind classification',
                                 'thickness_class' : 'Surfactant thickness classification' }

            for var in garbage_long_name:
                if var not in of.variables:
                    of.createVariable(var,'f4',('y','x'),fill_value=1e20)
                    of[var].setncatts(ds['sigma'].__dict__)
                    of[var].long_name = garbage_long_name[var]

            # cr_levels_str -- added as ATTRIBUTE, not VARIABLE
            cr_cum_levels = (0, .5, .80, .90, .95, .98, 1.0)
            cr_levels = []
            CU = datadict['cumulative']['data']
            CR = datadict['contrast_ratio']['data']
            for frac in cr_cum_levels:
                idx0,idx1 = np.unravel_index(np.argmin(np.abs(CU - frac)),CU.shape)
                cr_levels.append(CR[idx0,idx1])

            of.setncattr_string('cr_cum_levels_str',' '.join('%.2f'%x for x in cr_cum_levels))
            of.setncattr_string('cr_levels_str'    ,' '.join('%.2f'%x for x in cr_levels))


    elif format == 'GEOTIFF':
        ## https://gis.stackexchange.com/questions/164853/reading-modifying-and-writing-a-geotiff-with-gdal-in-python
        ## https://gdal.org/drivers/raster/gtiff.html#raster-gtiff
        ## https://gdal.org/tutorials/raster_api_tut.html

        # get metadata from input file
        ds           = gdal.Open(infile)
        GeoTransform = ds.GetGeoTransform()
        Projection   = ds.GetProjection()
        ds           = None

        driver = gdal.GetDriverByName("GTiff")

        for key in datadict:
            outfile = abspath(path.join(outdir, f'{infile_basename}_{ID}_{key}.tif'))
            rows, cols = datadict[key]['data'].shape

            outarr = datadict[key]['data'].data.copy()
            outarr[datadict[key]['data'].mask] = np.nan

            outdata = driver.Create(outfile, cols, rows, 1, gdal.GDT_Float32)
            outdata.SetGeoTransform(GeoTransform)          ##sets same geotransform as input
            outdata.SetProjection(Projection)              ##sets same projection as input

            outdata.GetRasterBand(1).WriteArray(outarr)
            outdata.FlushCache()                           ##saves to disk
            outdata = None

    elif format == 'ENVI' or True:
        sfx = suffix if (suffix == '.grd') else '.dat'
        for key in datadict:
            outfile = abspath(path.join(outdir, f'{infile_basename}_{ID}_{key+sfx}'))
            spectral.envi.save_image(outfile + '.hdr',
                                     datadict[key]['data'],
                                     force      = True,
                                     dtype      = np.float32,
                                     ext        = None,
                                     interleave = cfg['interleave'],
                                     metadata   = cfg['metadata'])
