function [off_frames,moloffwin_np] = findNP_mol_off_win(fname_guess,npBoxSz,moloffwin_np)
% NAME:
%               findGNR_mol_off_win
% AUTHOR:
%               Originally written by Zechariah Pfaffenberger
%               Last updated: 07-22-22
% PURPOSE:
%               This script will determine which frames are off-frames 
%               when there is no molecule localized near a nanoparticle.
%
% CATEGORY:
%               Image Analysis
%
% CALLING SEQUENCE:
%                off_frames = findNP_mol_off_win(fname_guess,npBoxSz,moloffwin_np)
%
% DEPENDENCIES:
%               appendNPs: this script should be called from inside
%               appendNPs which is the overarching analysis script to take
%               the localized nanoparticles, add them to the guess list,
%               and determine the off-frames for the nanoparticles
%
%
% INPUTS:
%              fname_guess: the name of the .mat file with the guesses,
%              full path
%
%              npBoxSz: is the  size of a diffraction limited spot for particles in pixels. It's the
%              nominal diameter, NOT the FWHM or something similar. Integer please!
%
%              moloffwin_np: the number of frames around the current frame to use for the BGSUB. For
%              example if moloffwin=50 then you would subtract 25 frames before and after the
%              current frame. Even number please!
%
%
% OUTPUTS:
%             off_frames: cell array, list of frames for each guess that
%             had no additional molecules localized in a region. This gets
%             saved to a newly created file in the same directory as the
%             overlay movie you gave to appendNPs.m with the filename of
%             that movie plus "_Mol_off_frames_NP.mat"

%% Import files

% %Pull out the file name and path to be used later for saving
[~,fname,~] = fileparts(fname_guess);
%load in the guesses & the movie size
load(fname_guess,'guesses','movsz');


%% Set up inputs

% Erroring out if dfrlmsz isn't an integer, because it's a "strong"
% parameter. You should really input what you mean.
if npBoxSz~=round(npBoxSz);error('dfrlmsz must be an integer');end
%rounding moloffwin
if moloffwin_np~=(ceil(moloffwin_np/2)*2)
    moloffwin_np=(ceil(moloffwin_np/2)*2);
    warning(['moloffwin must be an even integer. It has been reset to avgwin = ',num2str(moloffwin_np)])
end

%% Loop through all NP guesses and determine the off frames

%Note: This code is copied directly from Mol_off_frames in SMALL-LABS.
%Credit to Ben Isachoff as original author. 

%cell array of off frames vectors, for each localization (each row of
%guesses) a list of frames to include
off_frames=cell(size(guesses,1),1);
tic;
h1=waitbar(0);
set(findall(h1,'type','text'),'Interpreter','none');
waitbar(0,h1,['Creating off frames list for ',fname]);
for ii=1:movsz(3)
    try; waitbar(ii/movsz(3),h1); end
    %number of molecules in the current frame
    frmrows=find(guesses(:,1)==ii);
    nummol=length(frmrows);
    
    %skip it if there aren't any guesses in the current frame
    if nummol~=0
        %determine the frame list of frames to check for the current frame
        if ii<=(moloffwin_np/2)%the first group of frames
            frmlst=1:(moloffwin_np+1);
            frmlst=frmlst(frmlst~=ii);
        elseif ii>=(movsz(3)-moloffwin_np/2)%the last group of frames
            frmlst=movsz(3)+((-moloffwin_np):0);
            frmlst=frmlst(frmlst~=ii);
        else %all the frames in the middle
            frmlst=ii+[(-moloffwin_np/2):-1,1:(moloffwin_np/2)];
        end
        
        %find the rows of the guesses that are in the current frame list,
        %then pull out the row & column #'s
        mols2frms=find(ismembc(guesses(:,1),frmlst));
        allr=guesses(mols2frms,2);
        allc=guesses(mols2frms,3);
        
        for jj=frmrows(1):frmrows(end)
            %current molecule's position
            molr=guesses(jj,2);
            molc=guesses(jj,3);
            dists=(abs(allr-molr)<2*npBoxSz) & (abs(allc-molc)<2*npBoxSz);
            
            %save the lists of frames in which the current localization is off.
            off_frames{jj}=frmlst(~ismembc(frmlst,guesses((mols2frms(dists,1)),1)));
            
            % NOTE : when I first wrote this I thought that there needed to
            % be a unique as shown below. In reviewing this though, I don't
            % think it's needed and there is a big speed boost from leaving
            % it out. Please let me know if you find that this is needed.
            %off_frames{jj}=frmlst(~ismembc(frmlst,guesses(unique(mols2frms(dists,1)),1)));
            
            if length(off_frames{jj})<(0.05*moloffwin_np)
                warning(['Guess ',num2str(jj),' only has ',...
                    num2str(length(off_frames{jj})),' off frames. Consider increasing moloffwin'])
            end
        end
    end
end
try;close(h1);end

tictoc=toc;%the time to run the entire program
end

