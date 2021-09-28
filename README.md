# RayDec
RayDec is a matlab code to estimate the ellipticity of Rayleigh waves from 3-component single-station recordings.

You can find the details about the code in the original publication:
M. Hobiger, P.-Y. Bard, C. Cornou, and N. Le Bihan (2009). Single station determination of Rayleigh wave ellipticity by using the random decrement technique (RayDec), Geophys. Res. Lett. 36, L14303, doi: 10.1029/2009GL038863.
The attached file "RayDec-Description.pdf" is an edited version of this paper and was extracted from my PhD thesis ("Polarization of surface waves : characterization, inversion and application to seismic hazard assessment") at the Université de Grenoble. You can download the full text of this thesis at https://hal.univ-grenoble-alpes.fr/tel-00577887.

Example of usage:
[flx, elx] = raydec1station2021(vert, north, east, time, .1, 50, 100, 10, .1, 10)
analyzes the given signals between 0.1 and 50 Hz with 100 frequency steps, using 10 cycles, .1 as dfpar, and cutting the signal in 10 time windows.
The output of the given RayDec code are the two matrices flx and elx, both of size "frequency steps" x "time windows" and each column corresponds to one time window. All columns of flx are the same, so the frequency list is fl=flx(:,1).
To obtain an average RayDec ellipticity of all time windows, you can use
ellipticity_mean = exp(mean(log(elx), 2)),
ellipticity_error = exp(std(log(elx), 0, 2)).
This gives the ellipticity in a logarithmic scale and the error factor. 
ellipticity_mean .* ellipticity_error and ellipticity_mean ./ ellipticity_error are then the +/- 1 standard deviation ellipticity values. 
It is best to save a file with [fl, ellipticity_mean, ellipticity_error].

Usage in dinver
The ellipticity value and error factor obtained in the way presented above can be loaded in dinver (included in the geopsy package, http://www.geopsy.org/) and be inverted for the velocity profile of the underground. To read the file saved above in a correct way, use "Ellipticity (H/V)" for the second column and "log(H/V) stddev (approx.)" for the third column. The factor has to be set to "-1" for the second and third columns to attribute them to retrograde particle motion and stay "+1" for prograde particle motion. The fundamental mode of Rayleigh waves is always retrograde at low and high frequencies, but can be prograde in the intermediate range if a strong velocity contrast exists in the subsurface.

Be aware that an ellipticity curve alone is not sufficient to retrieve the velocity profile without additional constraints (e.g. Love/Rayleigh wave dispersion curves, interface depths, constraints on the velocities of certain layers)!

For a study on which parts of the ellipticity curve to use, please read:
M. Hobiger, C. Cornou, M. Wathelet, G. Di Giulio, B. Knapmeyer-Endrun, F. Renalier, P.-Y. Bard, A. Savvaidis, S. Hailemikael, N. Le Bihan, et al. (2013). Ground structure imaging by inversions of Rayleigh wave ellipticity: Sensitivity analysis and application to European strong motion sites, Geophys. J. Int. 192, 201–229.

For a more recent example of how to use the Rayleigh wave ellipticity, please read:
M. Hobiger, P. Bergamo, W. Imperatori, F. Panzera, A. Marrios Lontsi, V. Perron, C. Michel, J. Burjánek, and D. Fäh (2021). Site Characterization of Swiss Strong-Motion Stations: The Benefit of Advanced Processing Algorithms, Bull. Seismol. Soc. Am., doi: 10.1785/0120200316

https://zenodo.org/badge/378194540.svg
