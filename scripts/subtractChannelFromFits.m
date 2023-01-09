function subtractChannelFromFits(file_or_directory_name,varargin)
% NAME:
%               subtractChannelFromFits
% AUTHOR:
%               Originally written by Zechariah Pfaffenberger
%               Last updated: 08-03-22
% PURPOSE:
%              This will be the final step in two channel emission
%              analysis. It will calculate the summed intensity of
%              molecules in the center and off center channels and the
%              intensity of the nanoparticles themselves if needed.
%
% CATEGORY:
%               Image Analysis
%
% CALLING SEQUENCE:
%                subtractChannelFromFits(file_or_directory_name)
%
% DEPENDENCIES:
%              This script assumes you've run selectNanoparticlesForOverlay
%              and cropAndOverlay scripts to get an overlay movie and
%              center and offcenter movies. Additionaly, if you are running
%              background subtraction you need to have run appendNPs and
%              findNP_mol_off_win to ensure you have off frames lists
%
%
% INPUTS:
%              file_or_directory_name: the names of the overlay movies
%
%              Changeable parameters
%              bgsub: 1 or 0, MAKE SURE YOU KNOW WHETHER YOU HAVE DONE
%              BACKGROUND SUBTRACTION IN SMALL-LABS AND SET THIS
%              ACCORDINGLY. 1 is do background subtraction, 0 is don't.
% OUTPUTS:
%
%% Create filename list to loop through

[dlocs,dnames,exts] = importFile(file_or_directory_name);

%% Set up and parameters
%Set up control variables and constants
%Tell whether or not you used background subtraction on the movie
params.bgsub=1;

%Tell whether to do a mean or median for the subtraction
params.do_avgsub = 0;

%Tell whether or not to include the NPs in the calculation of the true
%intensity in each channel
params.includeNPs = 1;

%sublength of movie (should always be 1)
params.sublen = 1;


%Bigger box size for subtraction
params.NP_biggerbox=10;

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



%% Loop through the files
for ml=1:1:numel(dlocs)
    %Load in the name of the movie. Should be the overlay movie
    [filepath,filename,~] = fileparts([dlocs{ml},filesep,dnames{ml},exts{ml}]);
    
    %Find and load the guessfile. Different filename depending on whether
    %you have done background subtraction or not.
    if params.bgsub
        fname_guess = [filepath filesep filename,'_avgsub_guesses.mat'];
        fname_fits = [filepath filesep filename,'_AccBGSUB_fits.mat'];
        np_off_frames_fname=fullfile([filepath filesep filename '_Mol_off_frames_NP.mat']);
        
        load(fname_guess,'guesses','mol_np_log');
        load(np_off_frames_fname,'off_frames','moloffwin_np');
        fitdata = load(fname_fits);
        num_NPs= size(fitdata.directOutput.particle_crds,2);
        dfrlmsz=fitdata.dfrlmsz; %load the radius of the box size used for intensity calculation
        movsz=fitdata.movsz; %load the size of the movie
        moloffwin=fitdata.moloffwin;%load the moloffwin of the molecules
    else
        fname_guess = [filepath filesep filename,'_guesses.mat'];
        fname_fits = [filepath filesep filename,'_fits.mat'];
        load(fname_guess,'guesses','mol_np_log');
        fitdata = load(fname_fits);
        dfrlmsz=fitdata.dfrlmsz; %load the radius of the box size used for intensity calculation
        movsz=fitdata.movsz; %load the size of the movie
    end
    
    
    
    molr=guesses(:,2);
    molc=guesses(:,3);
    framelist=guesses(:,1);
    
    %Separate the filename and find what should be the center and off
    %center movie names along with the fits file
    % Split the filename string at underscores
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
    
    
    
    %% Loop through the two channels
    %We will loop through the same procedure for both on and off center
    %channels. 1 will be on center, 2 off center
    for ch=1:2
        fname_mov = [filepath filesep movFnames{ch},'.mat'];
        load(fname_mov,'mov');
        
        %% Initialize fits structure for saving
        %initializing the fits structure
        fits.frame=NaN(size(guesses,1),1);%frame numbers
        fits.sum=NaN(size(guesses,1),1);%sum of pixels in ROI around guess
        fits.bg_sum=NaN(size(guesses,1),1);%sum of pixels in ROI around guess
        
        h1=waitbar(0);
        set(findall(h1,'type','text'),'Interpreter','none');
        waitbar(0,h1,['Subtracting ',movFnames{ch}]);
        disp(['Subtracting ',movFnames{ch}]);
        NP_biggerbox = params.NP_biggerbox;
        if params.bgsub
            fits.bg_sum=NaN(size(guesses,1),1);%sum of pixels in ROI around guess
            %Get the filename for the Nanoparticles off frames
    
            for ii=1:size(guesses,1)
                try; waitbar(ii/size(guesses,1),h1); end
                curfrmnum=framelist(ii);
                curmolr=molr(ii);
                curmolc=molc(ii);
                frmlst=off_frames{ii};
                if isnan(curmolr)
                    continue
                end
                if params.do_avgsub
                    mean_mov=mean(single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),frmlst)),3);
                else
                    mean_mov=median(single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),frmlst)),3);
                end
                
                
                if mol_np_log(ii)==0
                    %                     moloffwin_whichone=moloffwin_GNR;
                    
                    
                    mean_mov_bigger=mean(mov(curmolr+(-NP_biggerbox:NP_biggerbox),curmolc+(-NP_biggerbox:NP_biggerbox),frmlst),3);
                    GNR_biggersum=sum(mean_mov_bigger(:));
                    GNR_sum=sum(mean_mov(:));
                    %avg intensity per pixel of the background times the pixel size
                    %of GNR dfrlmsz.so it is the sum background intensity of GNR
                    %dfrlmsz
                    avg_bg=(dfrlmsz*2+1)^2*(GNR_biggersum-GNR_sum)/((NP_biggerbox*2+1)^2-(dfrlmsz*2+1)^2); %avg bg intensity times the size of diffraction limit spot
                    GNR_insty=GNR_sum-avg_bg;
                    
                    
                    fits.bg_sum(ii)=GNR_insty;
                    fits.sum(ii)=0; %fit.sum is for molecule intensity but because the way the code has been written. the nanorod intensity is also recognized
                    %as molecule. So set them to zeros. it
                    
                else
                    
                    %the molecule image
                    molim=single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),curfrmnum));
                    %the subtracted image
                    data=molim-mean_mov;
                    
                    %putting the fit results into the fits structure
                    fits.sum(ii)=sum(data(:));%sum of pixels in ROI around guess
                end
            end
        else
            
            for ii=1:size(guesses,1)
                try; waitbar(ii/size(guesses,1),h1); end
                curfrmnum=framelist(ii);
                curmolr=molr(ii);
                curmolc=molc(ii);
                
                
                
                if mol_np_log(ii)==0
                    %the molecule image
                    try
                        data = single(mean(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),:),3));
                        mean_mov_bigger=mean(mov(curmolr+(-NP_biggerbox:NP_biggerbox),curmolc+(-NP_biggerbox:NP_biggerbox),:),3);
                        GNR_biggersum=sum(mean_mov_bigger(:));
                        GNR_sum=sum(data(:));
                        %avg intensity per pixel of the background times the pixel size
                        %of GNR dfrlmsz. so it is the sum background intensity of GNR
                        %dfrlmsz
                        %You need to ensure that you are just getting
                        %photon counts that are only background counts to subtract away later, so you
                        %subtract away GNR_sum (counts in the diffraction
                        %limited spot size) from GNR_biggersum which is the
                        %background counts plus the GNR signal counts.
                        %However, these counts are not just from 1 pixel,
                        %they are from an area of pixels that's equal to
                        %the additional area you added onto the original
                        %diffraction limit box size so you divide by this
                        %area to get an average counts per pixel and finally
                        %multiply by the number of pixels in original
                        %diffraction limit spot size to get the total
                        %amount of background in the spot
                        avg_bg=(dfrlmsz*2+1)^2*(GNR_biggersum-GNR_sum)/((NP_biggerbox*2+1)^2-(dfrlmsz*2+1)^2); %avg bg intensity times the size of diffraction limit spot
                        GNR_insty=GNR_sum-avg_bg;
                        
                        
                        fits.bg_sum(ii)=avg_bg;
                        fits.sum(ii)=GNR_insty; %fit.sum is for molecule intensity but because the way the code has been written. the nanorod intensity is also recognized
                        %as molecule. So set them to zeros. it
                        
                    catch
                        warning('One of the nanoparticles had a NaN for localization');
                        continue
                    end
                else
                    data = single(mov(curmolr+(-dfrlmsz:dfrlmsz),curmolc+(-dfrlmsz:dfrlmsz),curfrmnum));
                    fits.sum(ii)=sum(data(:));
                end
            end
        end
        %% Save the summed intensities for each channel
        
        overlayfname=fullfile(fname_fits);
        overlay=matfile(overlayfname,'Writable',true);
        overlay_fits=overlay.fits;
        if ch==1
            %save the mole intenisty in the center channel
            overlay_fits.mol_intensity_center=fits.sum;
            overlay_fits.NP_center=fits.bg_sum;
        else
            %save the mole intensity in the offcenter channel
            overlay_fits.mol_intensity_offcenter=fits.sum;
            overlay_fits.NP_offcenter=fits.bg_sum;
        end
        
        if params.bgsub
            overlay_fits.NP_bigboxsize=NP_biggerbox;
            overlay_fits.NP_boxsize=dfrlmsz;
        end
        overlay.fits=overlay_fits;
        try
            close(h1);
        catch
        end
    end
    
end

% Work in progress function to check molecule and nanoparticle positions to
% make sure they aren't too close to an edge so that the program will not
% be able to crop out the box size needed to take the mean of the off
% frames
%
%     function mean_mov = make_mean_mov(boxSz,edgeDist)
%                 %checking that it's not outside the frame and that off_frames for this
%             %guess isn't empty
%             
%             if (curmolc>edgeDist && curmolc<(movsz(2)-edgeDist) && curmolr>edgeDist && curmolr<(movsz(1)-edgeDist)) && ...
%                     (curmolc>edgeDist && curmolc<(movsz(2)-edgeDist) && curmolr>edgeDist && curmolr<(movsz(1)-edgeDist))&& ...
%                     ~isempty(offFrames)
%                 %the average (or median) frame
%                 % mean_mov is made only from off frames. frmlst first index tells you
%                 % what window you are looking at in the current movie so
%                 % you can adjust the off frames to be relative to that
%                 try
%                     %the average (or median) frame
%                     if params.do_avgsub
%                         mean_mov=mean(single(mov(molRow+(-boxSz:boxSz),molCol+(-boxSz:boxSz),frmlst)),3);
%                     else
%                         mean_mov=median(single(mov(molRow+(-boxSz:boxSz),molCol+(-boxSz:boxSz),frmlst)),3);
%                     end
%                 catch
%                     disp([fname, ' Index exceeds matrix dimensions at mean_mov'])
%                 end
%             else
%             end
%     end

end

