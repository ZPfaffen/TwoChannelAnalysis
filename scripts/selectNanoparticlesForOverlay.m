function selectNanoparticlesForOverlay(file_or_directory_name,varargin)
% NAME:
%               selectNanoparticlesForOverlay
%
% AUTHOR:
%               Originally written by Zechariah Pfaffenberger
%               Last updated: 07-22-22
%
% PURPOSE:
%               This script allows the user to select a nanoparticle
%
% CATEGORY:
%               Image Analysis
%
% CALLING SEQUENCE:
%                selectNanoparticlesForOverlay(file_or_directory_name,varargin)
%
% DEPENDENCIES:
%
%
% INPUTS:
%              file_or_directory_name: directory, filepath, or text file listing
%              filepaths to movie saved in matlab v7.3 format (should be
%              created if you run SMALL-LABS)
%
%              varargin: variable arguments
%
% OUTPUTS:
%             Creates three matlab files that contain the 3D array (LxMxN,
%             L is row, M is column, N is frame)
%             variable mov which are the movies in the overlay, center, and
%             offcenter channels. 

%% Create filename list to loop through

[dlocs,dnames,exts] = importFile(file_or_directory_name);

%% Set up and parameters
%Set up constants


%pad size
params.pdsz=50;

%Turn on Autofind
params.autofind = 0;

%Turn on and off previousGuesses
params.loadGuessPos = 1;

%Turn on and off saving guess positions;
params.saveGuessPos = 0;

%Turn on whether to use the fit positions from last movie to guess in the
%next sequential moving. This could help account for drift.
params.updateSaveGuesses = 1;

%Toggle for checking nanoparticle guesses or not
params.checkGuesses = 0;

%Set up constants for box size

%pixel size in nm/pixel
params.pxsz = 102;

%numerical aperture of fluorescence objective
params.NA = 1.3;

%max wavelength of fluorescence emission in nm
params.fl_wl = 700;

%threshold for automatic finding particles
params.bpthrsh = 96;

%This is diffraction limited size in fluorescence (for fluorophore)
% params.fl_dfrlmsz = round(1.4 * (params.fl_wl/55)/(2*params.NA));
params.fl_dfrlmsz = 9;

%This is box size for fluorescence
params.flBoxSize = 4*params.fl_dfrlmsz;

%This is whether to use gaussian filter to find particles
params.gaussFilt = 0;

params.debug = 0;

%This is control variable whether testing the control point updating
params.cpTest = 1;

egdesz = params.fl_dfrlmsz;
mask_fname=[];


if ~mod(params.flBoxSize,2)
    params.flBoxSize = params.flBoxSize - 1;
end

paramsnames=fieldnames(params);
% if any parameters are included as inputs, change the parameter mentioned
if nargin>1
    for ii=1:2:nargin-2
        whichField = strcmp(paramsnames,varargin{ii});
        try
            eval(['params.' paramsnames{whichField} ' = varargin{ii+1};'])
        catch
            error([varargin{ii}, '  is not an input parameter. Check the spelling.'])
        end
    end
end

for ii_=1:numel(dlocs)
    %sublength of movie (should always be 1)
    params.sublen = 1;
    %% Section 1: Make Z Stack of movie and prepare for display
    [filepath,filename,ext] = fileparts([dlocs{ii_},filesep,dnames{ii_},exts{ii_}]);
    disp(['Running overlay on ',filename]);
    load([filepath filesep filename,'.mat'],'mov');
    params.movsz=size(mov);
    
    if params.sublen==1
        params.sublen = params.movsz(3);
    else
        params.sublen = params.sublen;
    end
    %This is box size for fluorescence
    params.flBoxSize = 3*params.fl_dfrlmsz;
    
    %making the phasemask logical map
    if ~isempty(mask_fname)
        if ischar(mask_fname)
            %the step is to get rid of the avgsub, note that this shouldn't
            %do anything if bgsub=0
            try
                load([pathstr,filesep,mask_fname,'.mat'],'PhaseMask')
            catch
                [datalist,dataloc,~]=uigetfile([pathstr,filesep,'*.*']);
                if ~iscell(datalist); datalist={datalist}; end
                datalist=[dataloc datalist];
                [dlocs,dnames,~]=cellfun(@fileparts,datalist,'uniformoutput',false);
                load([dlocs{1,1},filesep,dnames{1,2}],'PhaseMask')
            end
        else
            load([pathstr,filesep,strrep(fname,'_avgsub',[]),'_PhaseMask.mat'],'PhaseMask')
            
        end
        PhaseMasklg=PhaseMask;
        PhaseMask=PhaseMask;
        PhaseMasklg(PhaseMasklg~=0)=1;
        PhaseMasklg=logical(PhaseMasklg);
    else
        PhaseMasklg=true(params.movsz([1,2]));
        PhaseMask=true(params.movsz([1,2]));
    end
    
    %Create the z stack to display and pad it
    ROIshow=sum(double(mov(:,:,:)),3);
%     params.zlims = [0.001*max(max(ROIshow)),0.2*max(max(ROIshow))];
        params.zlims = [];
    ROIpad=padarray(ROIshow,[params.pdsz,params.pdsz],'symmetric');%pad it
    PhaseMasklgPad=padarray(PhaseMasklg,[params.pdsz,params.pdsz],'symmetric');%Pad phase mask logical
    PhaseMaskPad=padarray(PhaseMask,[params.pdsz,params.pdsz],'symmetric');%Pad phase mask
    guesses=zeros(1,2);
    
    %% Section 2: Select Nanoparticles in center and off center
    if ~params.loadGuessPos
        % Choose center channel nanoparticle
        ROItxtin = 110; %110 is the ASCII code for the letter n
        while (ROItxtin ~= 121 || isempty(ROItxtin))
            
            figure(1);
            imshow(ROIpad,params.zlims)
            title('Click on a nanoparticle in the center channel');
            rectangle('Position',[2,2,size(ROIpad,2)-2,size(ROIpad,1)-2],'EdgeColor','magenta',...
                'LineWidth',3)
            % Click and choose fiduciaries ginput returns
            % [col_1,row_1;col_2,row_2;...;col_n,row_n];
            % this is important because rows correspond to y points and columns to
            % x
            centGuesses = round(ginput);
            numFids = size(centGuesses, 1);
            hold on
            for jj=1:numFids
                pxb_x=centGuesses(jj,1)+[-(params.flBoxSize-1)/2,(params.flBoxSize-1)/2];
                pxb_y=centGuesses(jj,2)+[-(params.flBoxSize-1)/2,(params.flBoxSize-1)/2];
                
                rectangle('Position',[pxb_x(1) pxb_y(1) params.flBoxSize params.flBoxSize],'EdgeColor','r','LineWidth',2);
                
                text(centGuesses(jj,1)-(params.flBoxSize+3),centGuesses(jj,2)+(params.flBoxSize+3),num2str(jj),'color','red','fontsize',13,'fontweight','bold');
            end
            hold off
            ROItxtin=input('Do these guesses for center channel nanoparticle look OK to you (y/n)?  ','s');
        end
        centGuesses = fliplr(centGuesses);
        close(figure(1));
   
        
        %Choose off-center channel nanoparticle
        ROItxtin = 110; %110 is the ASCII code for the letter n
        while (ROItxtin ~= 121 || isempty(ROItxtin))
            %    imshow(ROIshow,prctile(ROIshow(:),[.1,99.8]))
            figure(2);
            imshow(ROIpad,params.zlims)
            title('Click on the same nanoparticle in the off-center channel');
            rectangle('Position',[2,2,size(ROIpad,2)-2,size(ROIpad,1)-2],'EdgeColor','blue',...
                'LineWidth',3)
            % Click and choose fiduciaries ginput returns
            % [col_1,row_1;col_2,row_2;...;col_n,row_n];
            % this is important because rows correspond to y points and columns to
            % x
            offcentGuesses = round(ginput);
            numFids = size(offcentGuesses, 1);
            hold on
            for jj=1:numFids
                pxb_x=offcentGuesses(jj,1)+[-(params.flBoxSize-1)/2,(params.flBoxSize-1)/2];
                pxb_y=offcentGuesses(jj,2)+[-(params.flBoxSize-1)/2,(params.flBoxSize-1)/2];
                
                rectangle('Position',[pxb_x(1) pxb_y(1) params.flBoxSize params.flBoxSize],'EdgeColor','r','LineWidth',2);
                
                text(offcentGuesses(jj,1)-(params.flBoxSize+3),offcentGuesses(jj,2)+(params.flBoxSize+3),num2str(jj),'color','red','fontsize',13,'fontweight','bold');
            end
            hold off
            ROItxtin=input('Do these guesses for off-center channel nanoparticle look OK to you (y/n)?  ','s');
        end
        offcentGuesses = fliplr(offcentGuesses);
        
        close(figure(2));
        if params.saveGuessPos
            savedGuesses.centRow = centGuesses(:,1);
            savedGuesses.centCol = centGuesses(:,2);
            savedGuesses.offcentRow = offcentGuesses(:,1);
            savedGuesses.offcentCol = offcentGuesses(:,2);
            save([filepath,filesep,'Saved_NPGuesses_ForOverlay.mat'],'savedGuesses');
        end
    elseif params.loadGuessPos
        if exist([filepath filesep 'Saved_NPguesses_ForOverlay.mat'], 'file') == 2
            load([filepath filesep 'Saved_NPguesses_ForOverlay.mat'],'savedGuesses');
        else
            error('Please create a file with the saved guesses at nanoparticle positions.')
            return
        end
        
        centguesses = [savedGuesses.centCol,savedGuesses.centRow];
        offcentguesses = [savedGuesses.offcentCol,savedGuesses.offcentRow];
        allGuesses = [centguesses;offcentguesses];
        numROIs = size(allGuesses, 1);
        if params.checkGuesses
            figure;
            imshow(ROIpad,params.zlims);
            rectangle('Position',[2,2,size(ROIpad,2)-2,size(ROIpad,1)-2],'EdgeColor',[0.9290 0.6940 0.1250],...
                'LineWidth',3)
            ROItxtin = 110; %110 is the ASCII code for the letter n
            while (ROItxtin ~= 121 || isempty(ROItxtin))
                hold on
                for jj=1:numROIs
                    pxb_x=allGuesses(jj,1)+[-(params.flBoxSize-1)/2,(params.flBoxSize-1)/2];
                    pxb_y=allGuesses(jj,2)+[-(params.flBoxSize-1)/2,(params.flBoxSize-1)/2];
                    
                    rectangle('Position',[pxb_x(1) pxb_y(1) params.flBoxSize params.flBoxSize],'EdgeColor','r','LineWidth',2);
                    
                    text(allGuesses(jj,1)-(params.flBoxSize+3),allGuesses(jj,2)+(params.flBoxSize+3),num2str(jj),'color','red','fontsize',13,'fontweight','bold');
                end
                hold off
                ROItxtin=input('Does this ROI selection for center (1) and offcenter (2) look OK to you (y/n)?  ','s');
            end
        end
        close gcf
        params.flBoxSize = params.flBoxSize*ones(numROIs,1);
        centGuesses = fliplr(centguesses);
        offcentGuesses =  fliplr(offcentguesses);
    end
    %Fit center channel positions
    [cent_prt_crds,cent_conf95log,cent_zstacks]=fitPartPosition(mov,centGuesses,0,params);
    
    %Fit enhancement particle positions
    [offcent_prt_crds,offcent_conf95log,offcent_zstacks]=fitPartPosition(mov,offcentGuesses,0,params);
    
    if params.updateSaveGuesses
        %We round the position guessed for the particle to the be the
        %nearest integer so they can be used as the guesses for the next
        %chronological movie.
        savedGuesses.centRow = round(squeeze(cent_prt_crds(:,:,1)))+params.pdsz;
        savedGuesses.centCol = round(squeeze(cent_prt_crds(:,:,2)))+params.pdsz;
        savedGuesses.offcentRow = round(squeeze(offcent_prt_crds(:,:,1)))+params.pdsz;
        savedGuesses.offcentCol = round(squeeze(offcent_prt_crds(:,:,2)))+params.pdsz;
        save([filepath,filesep,'Saved_NPGuesses_ForOverlay.mat'],'savedGuesses');
    end
    %% Section 3: Run cropAndOverlay
    [overlayMov,centMov,regOffCentMov,tformSim,RA,cropLims,ssimVals] = cropAndOverlay(ROIshow,cent_prt_crds,offcent_prt_crds,mov,params);
    
    %% Section 4: Save overlay and individual channel movies
    mov = overlayMov;
    save([filepath filesep filename,'_overlay.mat'],'mov','-v7.3');
    mov = centMov;
    save([filepath filesep filename,'_center.mat'],'mov','-v7.3');
    mov = regOffCentMov;
    save([filepath filesep filename,'_offcenter.mat'],'mov','-v7.3');
    save([filepath filesep filename,'_overlayinfo.mat'],'tformSim','RA','cropLims','ssimVals');
    if params.cpTest
        load('CPTestSSimScores.mat','cpScores');
        if isfield(cpScores,'movNum')
            cpScores.movNum = [cpScores.movNum; ii_];
        else
            cpScores.movNum = ii_;
        end
        save('CPTestSSimScores.mat','cpScores');
    end
    close all
end
end

