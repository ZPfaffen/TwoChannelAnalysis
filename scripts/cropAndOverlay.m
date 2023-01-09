function [overlayMov,centMov,regOffCentMov,tformSim,RA,cropLims,ssimVals] = cropAndOverlay(ROIimg,cent_prt_position,offcent_prt_position,mov,params)
% NAME:
%               cropAndOverlay
% AUTHOR:
%               Originally written by Zechariah Pfaffenberger
%               Last updated: 01-08-23
% PURPOSE:
%               This script will crop an area out surrouding off-center
%               and center channel nanoparticles that is the same size,
%               perform image registration on them, and then create an
%               overlay movie. The center channel is considered fixed and
%               the off-center channel is transformed for the registration.
%
% CATEGORY:
%               Image Analysis
%
% CALLING SEQUENCE:
%                cropAndOverlay(cent_prt_position,offcent_prt_position)
%
% DEPENDENCIES:
%
%
% INPUTS:
%              ROIimg: summed stack image from movie created in
%              selectNanoparticlesForOverlay.m
%
%              cent_prt_position: fit position of a particle in the center
%              channel
%
%              offcent_prt_position: fit position of a particle in the offcenter
%              channel
%              
%              mov: 3D movie array variable for the raw data taken on the
%              camera converted into matlab format
%             
%              
% OUTPUTS:
%             overlayMov- 3D array, uint16, sum of the registered off and on center
%             channel movies
%
%             centMov- 3D array, uint16, center channel movie alone
%
%             regOffCentMov- 3D array, uint16, off-center channel movie
%             alone
%
%             tformSim- image transformation object, data structure

%% Parameter set up

%Size of box around selected nanoparticle to crop out. Should be around the
%FWHM of the gaussian beam profile. Units in pixels.
params.ROIsize = 155;
params.checkROIs = 0;

%Positioning of the crop box. A few different options here.
%'centered' = particle will be the exact center of the box (e.g. plus or
%minus ROIsize/2 from that center particle position in both x and y
%dimensions)
%'upperleft' = particles you selected in selectNanoparticlesForOverlay
%will be considered the upper left of your crop box
%'bottomleft' = particles you selected in selectNanoparticlesForOverlay
%will be considered the bottom left of your crop box
%'bottommidleft' = particles you selected in selectNanoparticlesForOverlay
%will be considered the bottom middle left particle (basically 75% is in
%positive y and 25% negative y)
%'topmidleft' =  particles you selected in selectNanoparticlesForOverlay
%will be considered the top middle left particle (basically 
params.cropType = 'topmidleft';

%Type of image registration. Should be a string. Currently implemented
%Control points = 'cp'
%
%Optimize affine = 'oa'
%
%Phase correlation followed by optimization = 'phCorr'
params.registrationType = 'cp';

%Parameters for control point registration
params.manualCP = 0;
params.updateCp = 1;
params.saveTForm = 1;

%% Section 1: Find region for center and off-center channels
% Center channel ROI
roiCenter = squeeze(cent_prt_position);
roiOffCenter = squeeze(offcent_prt_position);

%Displacement vector is defined a center minus off center
dispVector = squeeze(cent_prt_position) - squeeze(offcent_prt_position);
%calculate the center ROI by rounding
roiOffCenter = roiCenter - dispVector;
roiCenter = round(roiCenter);
roiOffCenter = round(roiOffCenter);


switch params.cropType
    case 'centered'
        
        
        centROIpos = [(roiCenter(2) - params.ROIsize/2) (roiCenter(1) - params.ROIsize/2) params.ROIsize, params.ROIsize];
        offcentROIpos = [(roiOffCenter(2) - params.ROIsize/2) (roiOffCenter(1) - params.ROIsize/2) params.ROIsize params.ROIsize];
        
    case 'upperleft'
        
        centROIpos = [(roiCenter(2) - 0.125*(params.ROIsize)) (roiCenter(1) - 0.15*(params.ROIsize)) params.ROIsize, params.ROIsize];
        offcentROIpos = [(roiOffCenter(2) - 0.125*(params.ROIsize)) (roiOffCenter(1) - 0.15*(params.ROIsize)) params.ROIsize, params.ROIsize];
    case 'bottomleft'
        %We're specifying the "bottom" corner of the rectangle but the y
        %axis is different for our image (basically positive y is going
        %down) So you have to subtract to get a corner that is more towards the 
        %top in the image than the specified point
        centROIpos = [(roiCenter(2) - 0.125*(params.ROIsize)) (roiCenter(1) - 0.85*(params.ROIsize)) params.ROIsize, params.ROIsize];
        offcentROIpos = [(roiOffCenter(2) - 0.125*(params.ROIsize)) (roiOffCenter(1) - 0.85*(params.ROIsize)) params.ROIsize, params.ROIsize];
    case 'bottommidleft'
        centROIpos = [(roiCenter(2) - 0.15*(params.ROIsize)) (roiCenter(1) - 0.60*(params.ROIsize)) params.ROIsize, params.ROIsize];
        offcentROIpos = [(roiOffCenter(2) - 0.15*(params.ROIsize)) (roiOffCenter(1) - 0.60*(params.ROIsize)) params.ROIsize, params.ROIsize];
    case 'topmidleft'
        centROIpos = [(roiCenter(2) - 0.3*(params.ROIsize)) (roiCenter(1)  - 0.1*(params.ROIsize)) params.ROIsize, params.ROIsize];
        offcentROIpos = [(roiOffCenter(2) - 0.3*(params.ROIsize)) (roiOffCenter(1) - 0.1*(params.ROIsize)) params.ROIsize, params.ROIsize];
    case 'topmidright'
        centROIpos = [(roiCenter(2) - 0.65*(params.ROIsize)) (roiCenter(1)  - 0.1*(params.ROIsize)) params.ROIsize, params.ROIsize];
        offcentROIpos = [(roiOffCenter(2) - 0.65*(params.ROIsize)) (roiOffCenter(1) - 0.1*(params.ROIsize)) params.ROIsize, params.ROIsize];
        
        
end

%place the ROI positions in a variable called cropLims
cropLims = [centROIpos; offcentROIpos];

if params.checkROIs
    ROItxtin = 110;
    while (ROItxtin ~= 121 || isempty(ROItxtin))
        figure(11);
        imshow(ROIimg,params.zlims);
        hold on
        %Specify pos as a four-element vector of the form [x y w h] in data units.
        %x and y are bottom lower left corner of rectangle
        cTitle = text(roiCenter(2),(roiCenter(1) - params.ROIsize/2)-5,'Center');
        cTitle.Color = 'm';
        rectangle('Position',centROIpos,'EdgeColor','m')
        offTitle = text(roiOffCenter(2),(roiOffCenter(1) - params.ROIsize/2)-5,'Off-Center');
        offTitle.Color = 'b';
        rectangle('Position',offcentROIpos,'EdgeColor','b')
        hold off
        
        ROItxtin=input('Do these ROI selections look OK to you (y/n)?  ','s');
        
    end
end

%% Section 2: Crop out center and off-center regions
tic;
centMov = mov(centROIpos(2):centROIpos(2)+params.ROIsize-1,centROIpos(1):centROIpos(1)+params.ROIsize-1,:);
offCentMov =  mov(offcentROIpos(2):offcentROIpos(2)+params.ROIsize-1,offcentROIpos(1):offcentROIpos(1)+params.ROIsize-1,:);
zstackOffCent = sum(offCentMov,3);
zstackCent = sum(centMov,3);

centImLims = [min(min(zstackCent)) 0.9*max(max(zstackCent))];
offcentImLims = [min(min(zstackOffCent)) 0.9*max(max(zstackOffCent))];
centImg = mat2gray(zstackCent,centImLims);
offCentImg = mat2gray(zstackOffCent,offcentImLims);

%% Section 3: Register the z-stacks and apply transformation to every frame
figure;
imshowpair(offCentImg,centImg)
title('Unregistered')

ssimVals.ssimvalOg = ssim(centImg,offCentImg);

switch params.registrationType
    case 'oa'
        [optimizer,metric] = imregconfig('multimodal');
        
        regDefault = imregister(offCentImg,centImg,'affine',optimizer,metric);
        
        figure;
        imshowpair(regDefault,centImg);
        title('A: Default Registration')
        
        optimizer.InitialRadius = 0.0001;
        optimizer.MaximumIterations = 1000;
        optimizer.Epsilon = 1e-3;
        regUpdate = imregister(offCentImg,centImg,'affine',optimizer,metric,'PyramidLevels',1);
        
        figure;
        imshowpair(regUpdate,centImg);
        title('A: Improved Registration')
        tformSim = imregtform(offCentImg,centImg,'affine',optimizer,metric);
        RA = imref2d(size(centImg));
        ssimVals.ssimvalRg = ssim(centImg,regUpdate);
    case 'phCorr'
        if params.saveTForm
            %Run initial alignment with phase correlation
            tformEstimate = imregcorr(offCentImg,centImg);
            RA = imref2d(size(centImg));
            regPhase = imwarp(offCentImg,tformEstimate,'OutputView',RA);
            figure;
            imshowpair(regPhase,centImg,'falsecolor');
            title('A: Phase Correlation Register');
            
            %Run finer alignment with monomodal optimization based correlation
            [optimizer, metric] = imregconfig('multimodal');
            optimizer.InitialRadius = 1e-6;
            regUpdate = imregister(offCentImg, centImg,...
                'affine', optimizer, metric,'InitialTransformation',tformEstimate);
            
            figure;
            imshowpair(regUpdate,centImg);
            title('B: Optimized Registration')
            
            tformSim = imregtform(offCentImg,centImg,'affine',optimizer,metric);
            ssimVals.ssimvalRg = ssim(centImg,regUpdate);
            save('TransformationMatrix_AllMovies.mat','tformSim');
        else
            RA = imref2d(size(centImg)); 
            load('TransformationMatrix_AllMovies.mat','tformSim');
            regUpdate = imwarp(offCentImg,tformSim,'OutputView',RA);
            figure;
            imshowpair(regUpdate,centImg);
            title('B: Optimized Registration')
        end
    case 'cp'
        if params.manualCP
            [mp,fp] = cpselect(offCentImg,centImg,'Wait',true);
            tformSim = fitgeotrans(mp,fp,'projective');
            if params.saveTForm
                save('TransformationMatrix_AllMovies.mat','tformSim','mp','fp');
            end
        else
            load('TransformationMatrix_AllMovies.mat','tformSim','mp','fp');
        end
        RA = imref2d(size(centImg)); 
        regUpdate = imwarp(offCentImg,tformSim,'OutputView',RA);
        ssimVals.ssimvalRg = ssim(centImg,regUpdate);
        
        figure;
        imshowpair(regUpdate,centImg);
        title('A: Control Point Registration')
        
        if params.updateCp
            %Run cross correlation to adjust loaded in moving and fixed
            %points from last movie
            mpAdj = cpcorr(mp,fp,offCentImg,centImg);
            tformSimAdj = fitgeotrans(mpAdj,fp,'projective');
            regUpdateAdj = imwarp(offCentImg,tformSimAdj,'OutputView',RA);
            
            figure;
            imshowpair(regUpdateAdj,centImg); title('B: Control Point Registration w/ Corr')
            ssimValAdj = ssim(centImg,regUpdateAdj);
            if params.cpTest
                if isfile('CPTestSSimScores.mat')
                    load('CPTestSSimScores.mat','cpScores');
                    cpScores.ssimPChanges = [cpScores.ssimPChanges; (ssimValAdj/ssimVals.ssimvalOg-1)...
                        *100, (ssimVals.ssimvalRg/ssimVals.ssimvalOg-1)*100];
                    save('CPTestSSimScores.mat','cpScores');
                else
                    cpScores.ssimPChanges = [(ssimValAdj/ssimVals.ssimvalOg-1)...
                        *100, (ssimVals.ssimvalRg/ssimVals.ssimvalOg-1)*100];
                    save('CPTestSSimScores.mat','cpScores');
                end
            end
            %For update purposes, we will always save the t-matrix and
            %positions that results in the higher SSIM.
            if ssimValAdj > ssimVals.ssimvalRg
                ssimVals.ssimvalRg = ssimValAdj;
                regUpdate = regUpdateAdj;
               
                tformSim = tformSimAdj;
                mp = mpAdj;
                if params.saveTForm
                    save('TransformationMatrix_AllMovies.mat','tformSim','mp','fp');
                end
            end
        end
end

% Additional useful commands to pull out transformation matrix and spatial
% referencing object


%RregistrationEstimator(centImg,offCentImg);
% %Debugging
% figure;
% imshow(sum(centMov,3),params.zlims);


regOffCentMov = zeros(size(offCentMov));
for ii_=1:size(centMov,3)
    %Select current movie frame from the center channel
    currFrCent = centMov(:,:,ii_);
    %Select current movie frame from off center channel
    currFrOff = offCentMov(:,:,ii_);
    
    if params.debug
        figure;
        imshowpair(currFrOff,currFrCent)
        title('Unregistered')
    end
    currFrOff_ref = RA;
    %Use the imwarp command to apply the geometric transform object to the
    %current frame with the spatial referncing object.
    [testFr,currFrOff_ref] = imwarp(currFrOff,tformSim,'OutputView',currFrOff_ref);
    regOffCentMov(:,:,ii_) = testFr;
    
    if params.debug
        figure;
        imshowpair(regOffCentMov(:,:,ii_),currFrOff);
        title('Register v Unregistered')
    end
    
end
regOffCentMov = uint16(regOffCentMov);
overlayMov = regOffCentMov + centMov;
close all
end

