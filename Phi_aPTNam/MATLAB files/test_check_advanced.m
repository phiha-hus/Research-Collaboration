% Checking the non-advancedness of the index 2 DDAE
clear all; close all; clc

%% Example 1. Phi/P.T.Nam21
E = [1 0;0 0]
A = [0 1;1 0]
ep1 = 1 
ep2 = 0
Ad = [0 ep1;0 ep2]

[is_non_advanced,is_impulse_free] = check_advanced(E,A,Ad)

%% Example 2.

