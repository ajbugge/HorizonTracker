# HorizonTracker<br/>
<br/>
Code for the paper
### Automatic extraction of dislocated horizons from 3D seismic data using non-local trace matching <br/>
*Aina Juell Bugge, Jan Erik Lie, Andreas Kjelsrud Evensen, Jan Inge Faleide, and Stuart Clark*<br/>
Geophysics, 2019. Contact: aina.juell.bugge@gmail.com<br/>
<br/>
This repository contains python code and jupyter notebook for data-driven seismic horizon tracking using non-local dynamic time warping and unwrapped instantaneous phase. The jupyter notebook includes a tutorial with a small seismic survey.<br/>


![Tutorial results](Tutorial_results.png)<br/>

Required python packages: <br/>
-- numpy<br/>
-- scipy<br/>
-- matplotlib<br/>
-- skimage<br/>
-- time<br/>
-- tslearn<br/>
-- IPython.display<br/>
-- sys<br/>
<br/>

Extracting key horizons from seismic images is an important element of the seismic interpretation workflow. Although numerous computer-assisted horizon extraction methods exist, they are typically sensitive to structural and stratigraphic discontinuities. As a result, these computer-assisted methods have difficulties extracting non-coherent dislocated horizons. We present a new data-driven method to correlate, track and extract horizons from seismic volumes with complex geological structures. The proposed method correlates seismic horizons across discontinuities and does not require user input in the form of seed points or prior identification of faults. Furthermore, the method is robust towards amplitude changes along a seismic horizon and does not jump from peak to trough or vice versa. We use a large sliding window and match full-length seismic traces using non-local dynamic time warping to extract grids of correlated points for our target horizons. Through computed accuracy measurements, we discard non-accurate correlations before interpolating complete seismic horizons. Because the proposed method does not require manually picked seed points or prior structural restoration, it does not rely on interpretive experience or geological knowledge. <br/>
![Horizon correlation and interpolation](Figure11.png)<br/>
Figure from Bugge et al., 2019
