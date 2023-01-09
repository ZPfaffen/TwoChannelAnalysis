function appendNPs(file_or_directory_name,varargin)
% NAME:
%               appendNPs
% AUTHOR:
%               Originally written by Zechariah Pfaffenberger
%               Last updated: 08-02-22
% PURPOSE:
%               This script will append nanoparticle localizations from
%               findNanoparticle script to the end of the guess list
%               output from SMALL-LABS and designate which are
%               nanoparticles or molecules.
%
% CATEGORY:
%               Image Analysis
%
% CALLING SEQUENCE:
%                appendNPs(file_or_directory_name)
%
% DEPENDENCIES:
%
%
% INPUTS:
%              file_or_directory_name: the name of the file for the
%              overlay movie of two channels
%
%
%             Optional Parameters:
%
%              moloffwin_np: integer, the number of frames around the current frame to use for the BGSUB. For
%              example if moloffwin=50 then you would subtract 25 frames before and after the
%              current frame. Even number please!
%
%              npBoxSz: integer, size of the box around a nanoparticle that
%              is considered to find frames where there are no guesses
%
% OUTPUTS:
%               guesses: within your SMALL_LABS guess output (in a file
%               that is your movie filename plus "_avgsub_guesses.mat")
%               This script will append on the nanoparticle guess positions
%               to the "guesses" variable which is a Nx3 array (where N is
%               the number of localizations), Column 1 is the frame, column
%               2 is row position, column 3 is column position. Note you
%               can tell which rows in this array represent nanoparticles
%               because they will always be at the very end and have their
%               frame number set to 1
%
%               mol_np_log: additionally, in the same file as the guesses,
%               this script creates an Nx1 double variable (where N is the
%               number of guessed molecules plus the number of particles)
%               which is a logical variable telling you which guesses are
%               molecules and which are particles. 1 represents molecules,
%               0 represents particles.
%


%% Create filename list to loop through

[dlocs,dnames,exts] = importFile(file_or_directory_name);

%% Set up and parameters
%Set up control variables and constants
%Tell whether or not you used background subtraction on the movie
params.bgsub=1;

%Tell whether or not to find the off window for the nanoparticles
params.findNP_off_win = 1;

%sublength of movie (should always be 1)
params.sublen = 1;

%%% Parameters for finding the nanoparticle off frames
%The number of frames to use for the off frame window
params.moloffwin_np=2000;

%The size of the box that should have no guess in it to be considered an
%off frame
params.npBoxSz = 7;

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

%% Load in guess list
for ml=1:1:numel(dlocs)
    %Load in the name of the movie. Should be the overlay movie
    [filepath,filename,~] = fileparts([dlocs{ml},filesep,dnames{ml},exts{ml}]);
    
    if params.bgsub
        fitdata = load([filepath filesep filename,'_AccBGSUB_fits.mat']);
        fname_guess = [filepath filesep filename,'_avgsub_guesses.mat'];
        guessdata = matfile(fname_guess,'Writable',true); %make it writable
    else
        fitdata = load([filepath filesep filename,'_fits.mat']);
        fname_guess = [filepath filesep filename,'_guesses.mat'];
        guessdata = matfile(fname_guess,'Writable',true); %make it writable
    end
    
    numNPs = size(fitdata.directOutput.particle_crds,2);
    %format the nanoparticle locations along with a dummy frame number
    np_guesses = round([ones(1,size(fitdata.directOutput.particle_crds,2));...
        fitdata.directOutput.particle_crds(1,:,1);fitdata.directOutput.particle_crds(1,:,2)]');
    
    %Append nanoparticle locations to guess list
    if isprop(guessdata,'mol_np_log')
        guessdata.guesses(end-numNPs+1:end,:)=np_guesses;
%         fitdata.fits.goodfit(end-numNPs+1:end,:) = true(numNPs,1);
    else
        
        % append the NP location to the guess list like a molecule, also
        % create a new vector in the guess file called mol_np_log which
        % is the same length as the guess list
        % with 1 indicating it's a nanoparticle and 0 indicating a molecule
        guessdata.mol_np_log = [ones(size(guessdata.guesses,1),1);zeros(numNPs,1)];
        guessdata.guesses= [guessdata.guesses;np_guesses];
%         fitdata.fits.goodfit = [fitdata.fits.goodfit;true(numNPs,1)];
    end
    %% Run off window for nanoparticles
    if params.findNP_off_win
        [off_frames,moloffwin_np] = findNP_mol_off_win(fname_guess,params.npBoxSz,params.moloffwin_np);
        
        %Save the off window for the nanoparticles
        %the five lines below combines the GNR moloffwin with molecules moloffwin
        %for an example, 2900 moloffwin for GNR and 50 moloffwin for other
        %molecules
        moloffwin_fname=fullfile(filepath,[filename,'_avgsub_guesses_Mol_off_frames.mat']);
        off_frame_file_orginal=load(moloffwin_fname);
        off_frames_orginal=off_frame_file_orginal.off_frames;
        off_frames(1:end-numNPs)=off_frames_orginal;
        
        %save the data
        save([filepath filesep filename '_Mol_off_frames_NP.mat'],'off_frames','moloffwin_np','params')
    end
end


end




