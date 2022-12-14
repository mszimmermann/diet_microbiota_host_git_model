% File containing global variables used by most scripts
%
addpath(genpath('.\'));

rawdataFolder = '.\InputData\';
inputFolder = '.\ProcessedData\metabolomics\';
inputFolderSeq = '.\ProcessedData\sequencing\';
resultsFolder = '.\ProcessedData\example_output\'; 
processeddataFolder = '.\ProcessedData\';
outputFolder = '.\ProcessedData\output\';
figureFolder = '.\Figures\';
stringtieFolderRNA = '.\InputData\sequencing_data\ballGown_RNA\eggNOGann\';   
stringtieFolderDNA = '.\InputData\sequencing_data\ballGown_DNA\eggNOGann\';   

% metabolomics
massThreshold = 0.001;
RTthreshold = 0.15;

intensityNoise = 5000;

testing_mode_flag = 1;

% define colors to display speceies
species_colors = [70 170 150;...
            0.8*[70 170 150];...
            222 45 38;...
            0.9*[222 45 38];...
            0.8*[222 45 38];...
            0.7*[222 45 38];...
            0.6*[222 45 38];...
            0.5*[222 45 38];...
            8 81 156;...
            0.9*[8 81 156];...
            0.8*[8 81 156];...
            0.7*[8 81 156];...
            0.6*[8 81 156];...
            205 145 60]/255;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create matlab object with HMDB pathways from tables prepared in Python

% call script defining file dependenciesand global variables
%addpath(genpath('.\'));
%add_global_and_file_dependencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data requirements: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
% Files:
%
% Figures:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

