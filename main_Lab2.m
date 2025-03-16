clear;
close all;
clc;
addpath(genpath('scripts'))
addpath(genpath("functions"))
addpath(genpath("datas"))

%% Room dim, Source and receiver positions
room =[3.14, 4.38, 3.27]; 

rS_1 = [0.16, 0.155, 0.155];
rS_2 = [0.16, 0.155, 0.155];
rS_3 = [0.16, 0.155, 0.155];
rS_4 = [0.16, 0.155, 0.155];
rS_5 = [0.035, 2.19, 0.155];
rS_6 = [0.13, 0.13, 1.64];

rM_1 = [0.03, room(2) - 0.03, 0.05];
rM_2 = [0.025, room(2)/2, 0.05];
rM_3 = [0.015, room(2)/2, 1.64];
rM_4 = [1.58, room(2)/2, 1.64];
rM_5 = [0.16, 0.155, 0.155];
rM_6 = [0.015, 1.46, 1.64];
%% Compute green function from measurement : Example
% This script describe the way to compute the gren function for a room from meas
eg_gf_room_from_meas

%% plot_room_config
% This script is plotting the room and an example of measurement configuration, to have a visual representation
plot_room_config()

%% plot cross section with green function : Example
% This part serve just as an exemple, for understanding how the Gf simulation work
% for plotting use template !
eg_compute_gf_cross_sec

%% plot green function for a specific pair of source receiver : Example
% This part serve just as an exemple, for understanding how the Gf simulation work
% for plotting use template !
eg_gf_simulation_1pos

%% plot the green function reciprocity between sources and measurement
plot_gf_reciprocity

%% plot the green function for source and receiver in the corner
plot_gf_src_rec_corner

%% plot the green function for source in corner and receiver in [0,ly/2,0]
plot_gf_src_corner_rec_midY
%% plot the green function for source in corner and receiver in [lx/2,ly/2,lz/2]
plot_gf_src_corner_rec_Fmid

%% Contourplot
plot_gf_cross_sec

%% plot 8.4 of the book
plot_gf_through_room