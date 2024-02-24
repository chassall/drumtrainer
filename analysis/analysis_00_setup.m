%Code for Yan, Y., Hunt, L. T., & Hassall, C. D. (2023). 
%Reward positivity biases interval production in a continuous timing task. bioRxiv, 2023-07.

%run this code before all other matlab code.
currentDir = fileparts(which("analysis_00_setup.m"));
cd(currentDir);
addpath("./functions");
%below are some MATLAB add-ons. Change this to your local address.
matlabDir = "/Users/rh/Documents/MATLAB/";
addpath(matlabDir+"unfold");
addpath(matlabDir+"cbrewer");
addpath(matlabDir+"makefigure");
addpath(matlabDir+"subtightplot");

%EEGlab
addpath("/Users/rh/Documents/eeglab");

close all; clear all; init_unfold();