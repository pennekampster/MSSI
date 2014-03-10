Multiscale straightness index (MSSI)
====================================

This R script implements two functions to calculate and plot the MSSI measure for animal trajectory analysis proposed 
by Postlethwaite, C.M., P. Brown, und T.E. Dennis: „A New Multi-Scale Measure for Analysing Animal Movement Data“. 
Journal of Theoretical Biology 317 (2013): 175–185. doi:10.1016/j.jtbi.2012.10.007.

The calculate_MSSI function requires a unique ID, the time and X- and Y-coordinates. Different window sizes and time 
resolutions (i.e. granulosity) can be specified.

The plot_MSSI function provides a visualization of the MSSI across window sizes and time for a given trajectory.
