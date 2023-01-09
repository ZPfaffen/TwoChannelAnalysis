%% Two Channel Wrapper
% NAME:
%               
%
% AUTHOR:
%               Originally written by Zechariah Pfaffenberger
%               Last updated: 07-22-22
%
% PURPOSE:
%               This is a wrapper script to do everything needed for two
%               channel emission analysis. 
% 
% CATEGORY:
%               Image Analysis
%
%
% DEPENDENCIES:
%               SMALL-LABS-master:
%
%               fQuickAnalysis:
%               
%               NP_localization: 
%               
%               TwoChannelEmissionAnalysis:
%
%             
% INPUTS:

%
% OUTPUTS:
%            

%% User input

%Fullpath to either a directory, a movie file, or a .txt list of movies
%just as the typical input to SMALL-LABS

file_or_directory_name = '';

%% Run overlay scripts to put the two channels together and save the movies
selectNanoparticlesForOverlay(file_or_directory_name);

%% Set up name to overlay movie
[dlocs,dnames,exts] = importFile(file_or_directory_name);
for ii_=1:size(dlocs,1)
    spltname = strsplit(dnames{ml},'_');
    %Check if the final part of the name is overlay. It should be if you
    %are inputting things correctly.
    if strcmp(spltname{end},'overlay')
        ch1fname = strjoin([spltname(1:end-1),{'center'}],'_');
        ch2fname = strjoin([spltname(1:end-1),{'offcenter'}],'_');
        movFnames = {ch1fname,ch2fname};
    else
        error('Make sure overlay movie was run in SMALL-LABS and selected.')
    end
    
end

overlay_file_or_directory = finalOutput;
%% SMALL-LABS and Nanoparticle finding
%Run typical SMALL-LABS fitting on on the overlay movies (based on
%fQuickAnlaysis)
fQuickAnalysis(overlay_file_or_directory);

%% Append the nanoparticle positions and find their off frames

appendNPs();


%% Find each channel intensity based on the guess list generated for the
%overlay movie
subtractChannelFromFits():


