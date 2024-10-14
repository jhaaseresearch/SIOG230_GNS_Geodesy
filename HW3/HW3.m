clear all; close all; clc;

% Create function to read the broadcast ephemerides, file brdc0160.12n
em_fn = "brdc0160.12n"
em_file = "data/" + em_fn;

% read single line
fileID = fopen(em_file,'r');

% File header information is up to line 8, pass up to line 8
for h = 1:8
    fgetl(fileID);
end

% read 22 rows of the data
for id = 1:22
    % read in each line
    for ln = 1:
end