# HW 3

## Description
This folder contains the functions to convert the ephemerides matrix for a
given satellite at a given time (i.e., 3 input arguments). The function  'get_satpos.m' returns 
corresponding X, Y, and Z coordinates in the ECEF frame, using time (t), satellite number (sv), and ephemerides matrix (eph) for instance: 
[X, Y,Z] = get_satpos(t,sv,eph).
The file 'test.m' is an activation file that uses the attached functions to implement 2.2.3 assignment 1 in the book (including the answer to question 3).  
'residuals_plot.png' shows the vector difference (norm) between sp3 and broadcast XYZ 
position as a function of time for satellite 31.