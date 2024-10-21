# SIOG230_GNS_Geodesy
Homework repository for SIOG239, Fall 2024

hw03:   
(Due 2024-10-17 2pm 50%; 2024-10-18 10am 50%)  
====================================================  
Write a program to convert the ephemerides in a RINEX navigation file into an ECEF cartesian  coordinate system (Assignment 2.2.3 Geodesy Lab Book page 67) 
The folder products contains the files you will need. Please use the files for the days 2019-10-19 and 2019-10-20.  
Functions in matlab_utils are available for a head start to read the files.  
Note: the read_rinexn.m function has difficulty for some reason with the files from 2012-01-16.  

hw02:  
Write a function xyz2wgs.m to convert from xyz to wgs84 following the equations in the geodesy_lab book section 1.1.2 p. 3.  
Please write the function so that the input is a N x 4 array of time(sec), x(m), y(m), z(m) and the output is an array of time, long(deg), lat(deg), height(m).  
Set the time = 1 in your input data. 
  
To check your code, verify that:  
Lat = 45.8791, Lon = 4.6766, Ele = 432.4222 m  
corresponds to:  
X = 4433469.9438 m, Y = 362672.7267 m, Z = 4556211.6409 m  
  
To check that your code works for a list of N positions at N different time points repeat the x, y, z values above N times and set the time values to be 1:N

hw01:  
do a matlab tutorial such as matlab onramp:    
    https://matlabacademy.mathworks.com/details/matlab-onramp/gettingstarted    
familiarize yourself with github and create a repository branch  
    https://gcapes.github.io/git-course/02-local/index.html  
    Here is a cheatsheet:  
    https://training.github.com/downloads/github-git-cheat-sheet.pdf  
    http://ndpsoftware.com/git-cheatsheet.html#loc=remote_repo  
