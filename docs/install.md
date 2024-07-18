# auto_SAR_Ocean_Contrast \(autonomous Ocean Contrast Estimation for SAR Images\)

# Installation

## 1. Requirements
The GeoTIFF interface requires *libgdal 3.8.4* or later.  There are many ways 
to install this on your system, e.g., macports (*port install gdal*).

## 2. Installation without using Conda environments)

To install the libraries for this package, after 'git pull ...', run this from the OilClassification directory:
  
    pip install -e . 

## 3. Installation in new conda environment

If you use conda, then do this:
  
    conda create --name NEW-ENV python=<version no.>
    conda activate NEW-ENV
    conda install gdal
    cd <OilClassification>
    pip install -e .

