# VeloPro

VeloPro (Velocity Processing) package is the collection of functions to return the associated spectrum and subsequent processing from raw velocities.

![Examplefig300dpi](https://user-images.githubusercontent.com/108192400/175784368-b4b668f1-b05f-45ba-876c-2642cadfe47d.png)

It offers two Fourier Transform methods, the fast Fourier transform (FFT) and the Date-Compensated Discrete Fourier Transform (DCDFT) proposed by Sylvio Ferraz-Mello:

Ferraz-Mello, S. (1981). Estimation of Periods from Unequally Spaced Observations. The Astronomical Journal, 86, 619. https://doi.org/10.1086/112924

Additional functionality includes the fitting of the spectrum to a hyperbolic tangent function to automatically extract breakpoints in the spectrum for logistic regression used to study the slope of the spectrum. The use of this is described in:

Hocking, W. K., Dempsey, S., Wright, M., Taylor, P., & Fabry, F. (2021). Studies of relative contributions of internal gravity waves and 2‐D turbulence to tropospheric and lower‐stratospheric temporal wind spectra measured by a network of VHF windprofiler radars using a decade‐long data set in Canada. Quarterly Journal of the Royal Meteorological Society, 147(740), 3735–3758. https://doi.org/10.1002/qj.4152

for processing of velocities measured by VHF doppler radar.

Please cite the above if using this code for research purposes.
