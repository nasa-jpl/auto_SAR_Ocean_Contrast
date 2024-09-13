# auto_SAR_Ocean_Contrast \(autonomous Ocean Contrast Estimation for SAR Images\)

[![Language](https://img.shields.io/badge/Python-3.8-blue)](https://www.python.org/)
[![Latest version](https://img.shields.io/badge/latest%20version-v1.3-yellowgreen.svg)]()
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](http://www.apache.org/licenses/LICENSE-2.0)

### Update Notes:

```diff
+ 
```

**A Python module for estimating the contrast in a SAR image of the ocean surface, relative to clean water pixels.  The code is primarily intended to identify oil slicks, but can be used to identify any radar-dark feature in a scene that is not entirely radar-dark. The contrast ratio is often referred to as the damping ratio in the scientific literature concerning mineral oil slicks.**

**The contrast is defined as $\sigma^{clean}(\theta)\over\sigma(\theta)$, where $\sigma(\theta)$ is the calibrated Normalized Radar Cross Section (NRCS), which depends upon the incidence angle, $\theta$. The algorithm identifies high confidence clean sea pixels based on the statistics of the SAR intensity to calculate the ratio, as described in \[Jones, 2023\].**

**auto\_SAR\_Ocean\_Contrast is a Python module that has been tested with SAR images acquired in L\-, S\-, C\-, and X\-band and in polarizations modes VV, HH, HV, and VH.  It works on images in radar coordinates or on a georeferenced grid, and requires a map of the incidence angle for each pixel in the scene.  If land is in the scene, the user can mask it out in the provided NRCS data or optionally can specify that a user-provided land mask file be used.**

**The outputs are 1\) the contrast ratio for all ocean pixels in the scene, and 2\) the cumulative distribution function's value for those pixels that are identified as being radar-dark, with 0 corresponding the the radar-dark pixel with the lower contrast ratio and 1 corresponding to the radar-dark pixel with the highest contrast ratio. Radar-dark pixels are those identified as having a contrast ratio significantly higher than that of the clean water peak, which is centered about a contrast ratio value of 1.  The threshold is set adaptively based upon the statistics of the contrast ratio for the scene.  Various png figures are also output.**

**Sample input and output files for each input format option are available at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12702460.svg)](https://doi.org/10.5281/zenodo.12702460)**

Copyright (c) 2023-24 California Institute of Technology (“Caltech”). U.S. Government
sponsorship acknowledged.  All rights reserved.  
Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.  
You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0

Link: https://github.com/nasa-jpl/auto_SAR_Ocean_Contrast

Citation:  Mao, P., & Jones, C. E. (2024). Autonomous Ocean Contrast Estimation for SAR Images (auto_SAR_Ocean_Contrast), [doi.org/10.48588/h6h7-8y85](https://doi.org/10.48588/h6h7-8y85).


## 1. Authors

Cathleen E. Jones (JPL/Caltech; cathleen.e.jones@jpl.nasa.gov; ORCID: 0000-0002-2739-1545) developed the algorithm, which is described in (Jones, 2013), developed the first version in MATLAB, and tested and refined it for SAR data from UAVSAR, FSAR, ISRO/ASAR, Sentinel-1, TerraSAR-X, and Radarsat-2.

Peter Mao (JPL/Caltech; peter.mao@jpl.nasa.gov; ORCID: 0009-0002-5543-153X) translated the MATLAB code to Python, and further optimized the algorithm, and is developing its related module auto_OilSpill_ID for oil spill classification.

**Reference:** 
[Jones, C. E. (2023). An automated algorithm for calculating the ocean contrast in support of oil spill response. Marine Pollution Bulletin, 191, 114952.](https://www.sciencedirect.com/science/article/pii/S0025326X23003843)

## 2. [Installation](/docs/install.md)

## 3. [Instructions](/docs/instructions.md)

## 4. [Algorithm](/docs/algorithm.md)

## 5. Related Open Source Code

NOTE: We intend to add open source code to classify the radar-dark pixels in the scene based upon damping ratio, to aid in identifying the thickest oil in a slick for emergency response.

## 6. License
Copyright (c) 2023-24 California Institute of Technology (“Caltech”). U.S. Government
sponsorship acknowledged.
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided
that the following conditions are met:  
• Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.  
• Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.  
• Neither the name of Caltech nor its operating division, the Jet Propulsion Laboratory, nor the
names of its contributors may be used to endorse or promote products derived from this software
without specific prior written permission.  

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Open Source License Approved by Caltech/JPL
APACHE LICENSE, VERSION 2.0
• Text version: https://www.apache.org/licenses/LICENSE-2.0.txt
• SPDX short identifier: Apache-2.0
• OSI Approved License: https://opensource.org/licenses/Apache-2.0

## 7. Acknowledgement

This effort was funded by the NASA 2018 *A.37 Earth Science Applications: Disaster Risk Reduction and Response* call's  Marine Oil Thickness (MOST) project, P.I. Francis Monaldo (NOAA).

## 8. MOST Project
[NASA-NOAA tech will aid marine oil spill response](https://phys.org/news/2021-12-nasa-noaa-tech-aid-marine-oil.html)

[Hear about the MOST Project on YouTube](https://youtu.be/6brucBsqR-g)

[Detecting oil thickness to aid oil spill response](https://appliedsciences.nasa.gov/our-impact/news/detecting-oil-thickness-aid-oil-spill-response)

[A profile of MOST PI, Frank Monaldo](https://appliedsciences.nasa.gov/our-impact/people/frank-monaldo-making-most-technology-detect-oil-spills)



