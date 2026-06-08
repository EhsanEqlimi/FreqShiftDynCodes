% Main script to run FreqShiftDyn Project
clc;clear;close all;%clearvars;
% Ehsan Eqlimi, columbia univerity
% June 2026
%% Set directories
% Add fieldtrip toolbox
addpath('/Users/ehsaneqlimi/fieldtrip-20260518');
ft_defaults; % -> really important!
% Add eeglab toolbox
addpath(genpath('/Users/ehsaneqlimi/eeglab2026.0.0'));
eeglab
addpath("MATLAB scripts_Julio_WM_PeerJ2024/");
EEGPath='/Volumes/HaegensLab/DATA/project_categorization/EEG/'; %lab shared folder
Dir=dir(fullfile(EEGPath,'*.vhdr'));
for i=1:length(length(Dir))
    % Load the file
    %EEG=pop_loadbv(EEGPath,Dir(i).name);
    cgf=[];
    cfg.dataset=[EEGPath Dir(i).name];
    EEGDataFT=ft_preprocessing(cfg);
    Event=ft_read_event(cfg.dataset);

end