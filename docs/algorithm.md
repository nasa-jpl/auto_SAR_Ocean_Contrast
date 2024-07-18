# auto_SAR_Ocean_Contrast \(autonomous Ocean Contrast Estimation for SAR Images\)

## Algorithm Description

**The algorithm is described in [Jones 2023](https://www.sciencedirect.com/science/article/pii/S0025326X23003843).**

From [Jones 2023](https://www.sciencedirect.com/science/article/pii/S0025326X23003843):
*In remote sensing of the ocean, contrast in the measured intensity between clean water and other features is used to identify different objects on the ocean surface either directly or indirectly via alteration of the ocean wave spectrum.  The damping ratio, a measure of contrast, is increasingly used for operational oil spill monitoring as an aid or alternative to visual inspections by trained personnel, and can in some cases identify thicker oil in a slick.  A method is proposed for automatically calculating the contrast based upon the statistical properties of the measured intensity signals from the ocean surface, and shown to work well even for complex slick geometries.  The algorithm is demonstrated using synthetic aperture radar (SAR) data from UAVSAR and Sentinel-1 to show that it can handle multi-frequency and medium-to-high resolution data.  The algorithm's flexibility and computational simplicity makes it suitable for real-time processing to support oil spill response.*

*The damping ratio is a measure of the ocean contrast, which for synthetic aperture radar is calculated as the ratio of the backscatter intensity for unslicked water,* $\sigma^{clean}$*,  to the measured intensity,* $\sigma$,

$DR = \left(\sigma^{clean}(\theta)\over\sigma(\theta)\right)$.

*The calibrated intensity is referred to as the normalized radar cross section, NRCS, and is known to vary significantly with the incidence angle, θ.  (It is important to know that the values of the NRCS used in this equation are in linear units, not dB).
For oil spill applications in particular, calculating the contrast is often the first step in mapping the slick.  Values of the DR near unity indicate unslicked water and increasing values correlate with thicker oil layers or emulsified oil.*

*The algorithm identifies high-confidence clean ocean pixels within the scene based upon statistical analysis, uses those pixels to estimate* $σ^{clean}$, *which are then used to calculate the damping ratio for the entire scene.  This is done as a function of incidence angle.*

[Schematic](algorithm_schematic.png)

