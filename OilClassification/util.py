import numpy as np
from matplotlib.pyplot import *
from scipy import signal
from numpy.polynomial.polynomial import polyfit
from numpy.polynomial.polynomial import polyval
from os import path, makedirs
from os.path import abspath, basename
from netCDF4 import Dataset
from osgeo import gdal
import yaml, spectral

__all__ = ['id_pix_as_cleanwater',
           'fit_xy_from_AOImask']

def AOI_histogram(data, edges='scott', maskfile=None, mask=False, percentclip=[5,95]):
    """Histogram the area-of-interest.

    DATA is a masked array (add functionality for regular array).
    EDGES is a numeric list or a string indicating the binning algorithm:
          {auto|fd|doane|scott|stone|rice|sturges|sqrt}. Empirically, 'stone' is slow.
    MASKFILE specifies a path to a mask header file.
    MASK is a mask array (reject True, keep False).
    PERCENTCLIP is a 2-element list with the low and high rejection percentiles.

    Returns:
    0: Probability distribution histogram
    1: Array of bin edges
    2: Bin map (digitized data in masked array)"""

    # Read mask and incorporate it into DATA.mask
    if maskfile is not None:
        print("NOT IMPLEMENTED: read the mask file %s"%maskfile)
        print("and OR it with optional input MASK")
    if not isinstance(data, np.ma.core.MaskedArray):
        data = np.ma.MaskedArray(data)
    data.mask |= mask

    if isinstance(edges, str):
        clip = np.percentile(data.data[~data.mask], percentclip)
        # data = np.ma.masked_outside(data, clip[0], clip[1])

        # matlab code has mechanisms to set #bins between ~20 and 200. (scale by 1.5 if #edges < 20)
        edges = np.histogram_bin_edges(data.data[~data.mask], edges, range=clip)

    (histo, edges) = np.histogram(data.data[~data.mask],edges)
    binmap         = np.ma.masked_array(np.digitize(data,edges))
    binmap.mask    = data.mask
    prob           = histo/np.sum(~data.mask)
    return (prob, edges, binmap)

def smooth(data, nsamples):
    """Boxcar smoothing of DATA over NSAMPLES

    Linearly extrapolates NSAMPLES of DATA on either side."""
    assert(len(data) > nsamples) # Can't smooth with window bigger than data

    idx = np.r_[0:nsamples]
    # linear fits to the ends
    pp_front = polyfit(idx, data[0:nsamples], 1)
    pp_back  = polyfit(idx, data[-nsamples:], 1)
    # linear extrapolations
    front = polyval(np.r_[0:nsamples] - nsamples, pp_front)
    back  = polyval(np.r_[0:nsamples] + nsamples, pp_back)

    # cat the extrapolations before smoothing
    sdata = np.convolve(np.r_[front,data,back], np.ones(nsamples)/nsamples, 'same')

    # chop off the ends upon return.
    return sdata[nsamples:-nsamples]

def id_pix_as_cleanwater(data, tgtval=None, edges=None):
    """Given a PDF, return index values of pixels around the center of the peak
    that is most likely to be clean water.

    Input:
    DATA is a masked array of limited angular extent in which clean water pixels will be identified
    TGTVAL is the bootstrap target intensity of DATA. If TGTVAL is None, then bootstrapping steps are skipped.
    EDGES is the binning algorithm name or a list of the AOI histogram edge values.
    INFO is an optional string to be used as the plot title.

    Output:
    CLEAN is a logical array of shape BINMAP with TRUE on clean pixels (or a scalar FALSE)
    FAILS = TRUE if can't find a peak
    LASTVAL is the TGTVAL for the next iteration
    """
    # histogram the (masked) data, using supplied edges or algorithm if specified
    if edges is None:
        prob, edges, binmap = AOI_histogram(data) # use default edges
    else:
        prob, edges, binmap = AOI_histogram(data, edges)
    if np.amin(binmap) == len(edges):
        print("No data in this area of interest", end=" ")
        return False, True, None
    sigma = (edges[:-1]+edges[1:])/2
    binwidth = np.mean(np.diff(edges))
    prob = prob / binwidth # convert probability to probability density

    # smoothing and peakfinding parameters
    maxsmooth = 10
    peakfrac  =  0.5
    closelim  =  0.25 # selected binval must be within this fraction of the full range of TGTVAL
    minDeltaSigma = max(round(0.5/binwidth), 1)  # [bins] = [sigma]/[binwidth]
    minProminence = 0.01 # probability density
    closeTol  = None # defined if bootstrapping happens

    nbins = len(prob)

    # matlab's SMOOTH handles endpoints differently from convolution and filters
    # by shortening the length of the filter near the endpoints.
    boxlen = min(round(len(prob)/10), maxsmooth)
    sprob = smooth(prob, boxlen)

    peakidx , _properties = signal.find_peaks(sprob,
                                              distance=minDeltaSigma,
                                              prominence=minProminence)
    prominences, _leftbase, _rightbase = signal.peak_prominences(sprob, peakidx)
    _promNorm  = prominences/(max(sprob) - min(sprob))
    sigmaPeaks = sigma[peakidx]

    figure()
    plot(sigma, sprob)
    plot(sigma, prob, '.', markersize=2, alpha=.5)
    xlabel('$\sigma$ [dB]')
    ylabel('probability density')
    plot(sigma[peakidx], sprob[peakidx], 'r*', label="candidates")

    # Assumes that cleanwater is the upper prominent peak
    # indx0 is the index into the peakidx|prom* arrays
    if len(peakidx) == 0:
        print("No clean water peak in this area of interest", end=" ")
        clean   = False
        fails   = True
        lastval = None
    else:
        if len(peakidx) == 1:
            indx0 = -1
        else: # the matlab code takes the brighter peak if there are only 2
            # provisionally pick the most prominent peak, but disregard the dimmest peak
            indx0 = np.argmax(np.r_[0, prominences[1:]])
            if tgtval is None:
                # No bootstrap :: use brightest peak if it has higher peak probability than the most prominent one.
                if sprob[peakidx[-1]] > sprob[peakidx[indx0]]:
                    indx0 = -1
            else:
                # identify peaks that are within CLOSELIM tolerance of TGTVAL
                closeTol = closelim * (edges[-1] - edges[0])
                closeIdx = np.where(abs(sigmaPeaks - tgtval) <= closeTol)[0]
                if len(closeIdx) > 0:
                    if indx0 in closeIdx:
                        pass # keep indx0 as is
                    else: # pick the most prominent of the close peaks
                        indx0 = np.argmin(abs(prominences - max(prominences[closeIdx])))
                else: # no peaks in tolerance range; use max peak
                    indx0 = -1 # matlab code does not follow comments -- just take the last one for now.

        CleanH2O_bin  = peakidx[indx0]
        CleanH2O_prob = sprob[CleanH2O_bin]
        limval  = min(sprob) + (1 - peakfrac) * (CleanH2O_prob - min(sprob)) # corrects for the background values

        plot(sigma[CleanH2O_bin], CleanH2O_prob, 'go', markersize=8, alpha=.67, label='clean water peak')

        # find the first below limval points to either side of the peak
        tmp = np.where(sprob < limval)[0]

        if len(tmp) == 0 or tmp[0] >= CleanH2O_bin:
            binlo = 0
        else:
            binlo = tmp[tmp < CleanH2O_bin][-1]

        if len(tmp) == 0 or tmp[-1] <= CleanH2O_bin:
            binhi = nbins - 1
        else:
            binhi = tmp[tmp > CleanH2O_bin][0]

        # "Added to keep values near peak (not lopsided)"
        if tgtval is not None:
            deltalo = abs(CleanH2O_bin - binlo);
            deltahi = abs(CleanH2O_bin - binhi);
            if ((deltalo > 1.5 * deltahi and deltahi != nbins) or
                (deltahi > 1.5 * deltalo and deltalo != 1)):
                delta = min(deltalo, deltahi);
                binlo = max(CleanH2O_bin - delta, 0);
                binhi = min(CleanH2O_bin + delta, nbins - 1);

        plot(sigma[[binlo, binhi]], sprob[[binlo,binhi]],'bx', label='"near-peak" range')
        if closeTol is not None:
            yrange = ylim()
            plot(np.ones(2)*(tgtval - closeTol), yrange,'k:',label='bootstrap window')
            plot(np.ones(2)*(tgtval + closeTol), yrange,'k:')
        # find indices of cleanwater in image -- consider outputting this as a mask
        clean = np.logical_and(binlo <= binmap, binmap <= binhi)
        clean = np.logical_and(clean, ~binmap.mask)
        fails = False
        lastval = sigmaPeaks[indx0]

    legend()
    return clean, fails, lastval

def fit_xy_from_AOImask(xdata, ydata, cleanmask, Xlabel='', Ylabel='', Title=''):
    """Linear (or polynomial) fit to clean y(x) within an AOI.
    XDATA: 2D array of x values
    YDATA: 2D array of y values
    CLEANMASK: logical mask (clean = TRUE, dirty = FALSE)

    return fit results for y(x)
    """
    pfitDeg =  5 # degree of polynomial fit
    nbins   = 50

    if xdata.shape != ydata.shape or xdata.shape != cleanmask.shape:
        print('Error: Input arrays are not the same size.');
        status = False;
        return status, None, str(pfitDeg), nbins

    # exclude unclean water and invalid data in either X or Y.
    datamask = ~cleanmask | xdata.mask | ydata.mask

    xvals = xdata.data[~datamask]
    yvals = ydata.data[~datamask]

    xEdges = np.linspace(np.amin(xvals), np.amax(xvals), nbins+1)
    bincnt = np.histogram(xvals, xEdges)
    indx   =  np.digitize(xvals, xEdges) - 1

    xBins   = (xEdges[0:-1] + xEdges[1:])/2
    ybinned = np.array([np.mean(yvals[indx==kk]) for kk in range(nbins)])

    idxOK = np.isfinite(xBins) & np.isfinite(ybinned)
    pfitCoeff = polyfit(xBins[idxOK], ybinned[idxOK], pfitDeg)
    yfit      = polyval(xEdges, pfitCoeff)

    figure()
    subplot(2,1,1)
    hist(xvals,xEdges);
    ylabel('# pixels')
    title(Title)

    subplot(2,1,2)
    plot(xBins , ybinned, 'o')
    plot(xEdges, yfit   , '-')
    xlabel(Xlabel)
    ylabel(Ylabel)

    status = True

    return status, pfitCoeff, str(pfitDeg), nbins

def writedata(datadict, infile, outdir, cfg):
    """Write data from processing infile to outfile."""
    (infile_basename, suffix) = path.splitext(basename(infile))
    ID     = cfg['ID']
    format = cfg['format'].upper()

    if format == 'NETCDF':
        outfile = abspath(path.join(outdir, f'{infile_basename}_{ID}_contrastratio.nc'))
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
            outfile = abspath(path.join(outdir, f'{infile_basename}-{ID}-{key}.tif'))
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
            outfile = abspath(path.join(outdir, f'{infile_basename}-{ID}-{key+sfx}'))
            spectral.envi.save_image(outfile + '.hdr',
                                     datadict[key]['data'],
                                     force      = True,
                                     dtype      = np.float32,
                                     ext        = None,
                                     interleave = cfg['interleave'],
                                     metadata   = cfg['metadata'])


def dumb_sigma(data, ratio=None, nsig=None):
    """Get a rough estimate for \sigma from the major peak of a distribution.

    Log-bins the data into 1000 bins and finds the low-side HWHM.  Use the
    Gaussian relationship HWHM = \sigma \sqrt{2 ln 2} to estimate \sigma.

    Accepts masked arrays as-is.

    If optional argument `ratio` is given, then return the data value for
    which the value-frequency exceeds the Gaussian-model by that factor.
    RATIO takes precedence over NSIG.

    If optional argument `nsig` is given, then return the data value that is
    NSIG higher than the peak value.

    """
    if np.ma.is_masked(data):
        thisdata = data[~data.mask]
    else:
        thisdata = data.reshape(-1)

    if any(thisdata <= 0):
        print("Data contains negative values.  Truncating...")
        thisdata[thisdata <= 0] = np.nan

    logmin = np.log10(np.nanmin(thisdata))
    logmax = np.log10(np.nanmax(thisdata))

    bins   = np.logspace(logmin, logmax, 1000)
    xx     = np.sqrt(bins[:-1] * bins[1:])
    yy,_   = np.histogram(thisdata, bins)

    imax   = np.argmax(yy)
    hmax   = np.argmin(abs(yy[:imax] - max(yy)/2))

    HWHM   = xx[imax] - xx[hmax]
    sigma  = HWHM / np.sqrt(2 * np.log(2))

    if nsig is None and ratio is None:
        return sigma
    elif nsig is not None and ratio is None:
        return xx[imax] + nsig * sigma
    else:
        dx = np.diff(bins)
        gauss = (yy[imax]/dx[imax]) * dx[imax:] * np.exp(-0.5 * ((xx[imax:] - xx[imax])/sigma)**2)
        ithreshold = np.argwhere(yy[imax:]/gauss > ratio)[0][0] + imax
        return xx[ithreshold]

def ma_expand(ma, dist=1):
    """Expand the valid segments of a masked array (MA) by DIST."""
    from skimage import segmentation as seg

    data = ma.data.copy()
    mask = ma.mask.copy()

    dmin = np.amin(data[~mask])
    data += 1 - dmin # offset data above zero by 1
    data[mask] = 0   # set masked values to zero for expansion into these positions

    data = seg.expand_labels(data, dist)
    mask[data > 0] = False # update the mask
    data += dmin - 1 # reset data values

    return np.ma.array(data, mask=mask)

def ma_remove_small(img, area=1, connectivity=2):
    """Remove isolated regiong with area less than AREA.

    default AREA is 1 pixel."""
    from skimage import morphology

    # connectivity:: 1 for edges, 2 for corners
    mask = morphology.remove_small_objects(~img.mask, area+.01, connectivity)

    return np.ma.array(img.data, mask=~mask)
