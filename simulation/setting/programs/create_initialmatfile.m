function [] = create_initialmatfile(id,newT,reinforcer,outInd,muscscale,yoke,plotOn)
% BABBLE_DASPNET_RESERVOIR Neural network model of the development of reduplicated canonical babbling in human infancy.
%
%   Modification of Izhikevich's (2007 Cerebral Cortex) daspnet.m and of a previous model described in Warlaumont (2012, 2013 ICDL-EpiRob).
%
%   Estimates of auditory salience calucated from a modified version of
%   Coath et. al. (2009) auditory salience algorithms.
%
%   For further reading:
%
%       Warlaumont A.S. (2012) Salience-based reinforcement of a spiking neural network leads to increased syllable production.
%       Proceedings of the 2013 IEEE Third Joint International Conference on Development and Learning and Epigenetic Robotics (ICDL). doi: 10.1109/DevLrn.2013.6652547
%
%       Izhikevich E.M. (2007) Solving the Distal Reward Problem through Linkage of STDP and Dopamine Signaling. Cerebral Cortex. doi: 10.1093/cercor/bhl152
%
%       The auditory salience estimator and auditory event detector used in salience-reinforced versions:
%           http://emcap.iua.upf.edu/downloads/content_final/auditory_saliency_model.html
%           http://www.mcg.uva.nl/papers/Coath-et-al-2007.pdf
%
%   Description of Input Arguments:
%       id          % Unique identifier for this simulation. Must not contain white space.
%       newT        % Time experiment is to run in seconds. Can specify new times (longer or shorter)
%                      for experimental runs by changing this value when a simulation is restarted.
%       reinforcer  % Type of reinforcement. Can be 'human', 'relhisal', or 'range'.
%       outInd      % Index of reservoir neurons that project to motor neurons. Length of this vector must be even. Recommended vector is 1:100
%       muscscale   % Scales activation sent to Praat. Recommended value is 4
%       yoke        % Indicates whether to run an experiment or yoked control simulation. Set to 'false' to run a regular simulation.
%                      Set to 'true' to run a yoked control. There must alread have been a simulation of the same name run with its
%                      data on the MATLAB path for the simulation to yoke to.
%       plotOn      % Enables plots of several simulation parameters. Set to 0 to disable plots, and 1 to enable.
%
%   Example of Use:
%       babble_daspnet_reservoir('Mortimer',7200,'relhisal',1:100,4,'false',0);
%
% Authors: Anne S. Warlaumont and Megan K. Finnegan
% Cognitive and Information Sciences
% University of California, Merced
% email: awarlaumont2@ucmerced.edu or anne.warlaumont@gmail.com or
% Website: http://www.annewarlaumont.org/lab/
% December 2014
% For updates, see https://github.com/AnneSWarlaumont/BabbleNN

%INITIALIZATIONS AND LOADING%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng shuffle;

% Initialization.
salthresh = 4.5;            % Initial salience value for reward (used in 'relhisal' reinforcment).
DAinc = 1;                  % Amount of dopamine given during reward.
sm = 4;                     % Maximum synaptic weight.
testint = 1;                % Number of seconds between vocalizations.

% Directory names for data.
wavdir = [id, '_Wave'];
firingsdir = [id, '_Firings'];
workspacedir = [id, '_Workspace'];
yokeworkspacedir = [id, '_YokedWorkspace'];

% Error Checking.
if(any(isspace(id)))
    disp('Please choose an id without spaces.');
    return
end

% Creating data directories.
if ~exist(wavdir, 'dir')
    mkdir(wavdir);
else
    addpath(wavdir);
end
if ~exist(firingsdir, 'dir')
    mkdir(firingsdir);
else
    addpath(firingsdir);
end
if ~exist(workspacedir, 'dir')
    mkdir(workspacedir);
else
    addpath(workspacedir);
end
if strcmp(yoke, 'true')
    if ~exist(yokeworkspacedir, 'dir')
        mkdir(yokeworkspacedir);
    else
        addpath(yokeworkspacedir);
    end
end

% Creating workspace names.
if strcmp(yoke,'false')
    workspaceFilename=[workspacedir,'/babble_daspnet_reservoir_',id,'.mat'];
elseif strcmp(yoke,'true')
    % Where to put the yoked control simulation data.
    workspaceFilename=[yokeworkspacedir,'/babble_daspnet_reservoir_',id,'_yoke.mat'];
    % Where to find the original simulation data.
    yokeSourceFilename=[workspacedir,'/babble_daspnet_reservoir_',id,'.mat'];
end

% Directory for Coath et. al. Saliency Detector.
addpath('auditorysaliencymodel');


if exist(workspaceFilename) > 0
    load(workspaceFilename);
else

    M=100;                 % number of synapses per neuron
    D=1;                   % maximal conduction delay
    % excitatory neurons   % inhibitory neurons      % total number
    Ne=800;                Ni=200;                   N=Ne+Ni;
    Nout = length(outInd);
    Nmot=Nout; % Number of motor neurons that the output neurons in the reservoir connect to.
    a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];     % Sets time scales of membrane recovery variable.
    d=[   8*ones(Ne,1);    2*ones(Ni,1)];       % Membrane recovery variable after-spike shift.
    a_mot=.02*ones(Nmot,1);
    d_mot=8*ones(Nmot,1);
    post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); % Assign the postsynaptic neurons for each neuron's synapse in the reservoir.
    post_mot=repmat(1:Nmot,Nout,1); % All output neurons connect to all motor neurons.
    
    post(:,8:end)=0; %amount of synap

    
    
    
    s=[rand(Ne,M);-rand(Ni,M)];         % Synaptic weights in the reservoir.
    sout=rand(Nout,Nmot); % Synaptic weights from the reservoir output neurons to the motor neurons.

    
    %for STDP
    for i=1:N
           pre{i}=find(post==i&s>0); %find pre excitatory neurons
           aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
    end;
    
    % Normalizing the synaptic weights.
    sout=sout./(mean(mean(sout)));

    sd=zeros(Nout,Nmot); % The change to be made to sout.

    for i=1:N
        delays{i,1}=1:M;
    end
    for i=1:Nout
        delays_mot{i,1}=1:Nmot;
    end
    STDP = zeros(Nout,1001+D);
    v = -65*ones(N,1);          % Membrane potentials.
    v_mot = -65*ones(Nmot,1);
    u = 0.2.*v;                 % Membrane recovery variable.
    u_mot = 0.2.*v_mot;
    firings=[-D 0];     % All reservoir neuron firings for the current second.
    outFirings=[-D 0];  % Output neuron spike timings.
    motFirings=[-D 0];  % Motor neuron spike timings.

    DA=0; % Level of dopamine above the baseline.

    muscsmooth=100; % Spike train data sent to Praat is smoothed by doing a 100 ms moving average.
    sec=0;

    rewcount=0;
    rew=[];

    % History variables.
    v_mot_hist={};
    sout_hist={};

    % Initializing reward policy variables.
    if strcmp(reinforcer,'relhisal')
        temprewhist=zeros(1,10); % Keeps track of rewards given at a threshold value for up to 10 previous sounds.
    end

    if strcmp(reinforcer, 'range')
        temprewhist=zeros(1,10);
        range_hist = NaN(newT, 1); % Keeps a record of the range over the entire simulation run time.
        rangethresh = 0.75; % Starting threshold for reward.
        Ranginc = .05; % How much to increase the threshold by.
    end

    % Variables for saving data.
    vlstsec = 0; % Record of when v_mot_hist was last saved.
    switch reinforcer  % Sets how often to save data.
        case 'human'
            SAVINTV = 10;
        otherwise
            SAVINTV = 100;
    end

end

T=newT;
clearvars newT;

% Absolute path where Praat can be found.
praatPathpc = 'c:\users\Takimto\Praat\Praat';
praatPathmac = '/Applications/Praat.app/Contents/MacOS/Praat';


% Special initializations for a yoked control.
if strcmp(yoke,'true')
    %load(yokeSourceFilename,'rew', 'yokedruntime');
    %if(T > yokedruntime)
     %   T = yokedruntime;
    %end
end


    % Every so often, save the workspace in case the simulation is interupted all data is not lost.
    if mod(sec,SAVINTV*testint)==0 || sec==T
        display('Data Saving..... Do not exit program.');
        save(workspaceFilename, '-regexp', '^(?!(v_mot_hist)$).');
    end


    % Writing motor neuron membrane potentials to a single large text file.
    if mod(sec,SAVINTV*testint)==0 || sec==T
        if strcmp(yoke, 'false')
            vmhist_fid = fopen([workspacedir,'/v_mot_hist_',id,'.txt'],'a');
        elseif strcmp(yoke, 'true')
            vmhist_fid = fopen([yokeworkspacedir,'/v_mot_hist_',id,'_yoke.txt'],'a');
        end
        if(vlstsec == 0)
            fprintf(vmhist_fid, 'History of Motor Neuron Membrane Potentials:\n');
        end
        for sindx = (vlstsec+1):sec % Going through all the seconds needed data saved for.
            % Information formated to make data more human readable.
            fprintf(vmhist_fid, '\nSecond: %i\n',sindx);
            fprintf(vmhist_fid, 'Millisecond:\n\t\t');
            fprintf(vmhist_fid, '%f\t\t%f\t\t', 1:1000);
            fprintf(vmhist_fid, '\n');
            fprintf(vmhist_fid, 'Neuron: \n');
            % Appending new voltage data.
            for nrow = 1:size(v_mot_hist{sindx},1)
               fprintf(vmhist_fid, '%i\t\t', nrow);
               fprintf(vmhist_fid, '%f\t\t%f\t\t', v_mot_hist{sindx}(nrow, :));
               fprintf(vmhist_fid, '\n');
            end
            fprintf(vmhist_fid, '\n');
        end

        fclose(vmhist_fid);
        vlstsec = sec; % Latest second of saving.
        save(workspaceFilename, 'vlstsec', '-append'); % Saving the value of the last written second in case the simulation is terminated and restarted.
        display('Data Saved.');
    end


end
