function [] = babbling(ID,newT,reinforce,outInd,muscscale,yoke,plotOn,feedbacktime,learningratio,speinplate,STDP,debug,IP,separatephase,Network,reward,feedbacktype)
% BABBLE_DASPNET_RESERVOIR Neural network model of the development of reduplicated canonical babbling in human infancy.
%
%   Modification of Izhikevich's (2007 Cerebral Cortex) daspnet.m and of a previous model described in Warlaumont (2012, 2013 ICDL-EpiRob).
%
%   Estimates of auditory salience calucated from a modified version of
%   Coath et. al. (2009) auditory salience algorithms.
%
%   For further reading:
%
%       Warlaumont A.S. (2012) Salience-based reinforcement of a spiking neural network BabbleNN-remakeleads to increased syllable production.
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
%       feedbacktime % Determine frequency of feedback.
%
%   Example of Use:
%       babbling('170123_1000_reinforce_100_4_NY_1_5_0.001_2_STDP',1000,'reinforce',1:100,4,'NY',1,5,0.001,2,'STDP',1,'IP_or_Tonic');
%
% Authors: Anne S. Warlaumont and Megan K. Finnegan
% Cognitive and Information Sciences
% University of California, Merced
% email: awarlaumont2@ucmerced.edu or anne.warlaumont@gmail.com or
% Website: http://www.annewarlaumont.org/lab/
% December 2014
% For updates, see https://github.com/AnneSWarlaumont/BabbleNN

%INITIALIZATIONS AND LOADING%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
rng shuffle;
warning('off','all');

global id itenumber iterate

% Initialization.
salthresh = 4.5;            % Initial salience value for reward (used in 'relhisal' reinforcment).
DAinc = 1;                  % Amount of dopamine given during reward.
sm = 4;                     % Maximum synaptic weight.
smr=10;                     % Maximum synaptic weight for reservoir
testint = 1;                % Number of seconds between vocalizations.
fftsize=2048;             %fft????????????????????????????????????
lpcsize=8;                %LPC size
inpInd=outInd;
tau=20;                   %decay parameter
SAVINTV=200;
LST_hist=[1:1000];
muscle_number=0;
negativereward = 0;
nega_rate=1/5;
mot_thre = 0;
constant_inplate = 0.01;
LTD = 0.00525;      %default 1.5  Li defo ->0.00525
LTP = 0.005;        %default 1.0  LI defo ->0.005
LTD_m = 0.00525;        %default 1.5 Li defo ->0.00525
LTP_m = 0.005;          %default 1.0 Li defo ->0.005
muscle_his = zeros(1000,newT);  %muscle_history for Sc, verticle=time horizontal=sec
yokemode = 'time';     %timespace or time 

%
sparse_mot = 1;   %1 or 0 
sparse_degree = 40; % number of synapse from out to motor Nout - sparse_degree default=40


%for low pass filter
cutoff_frequency = 10; % Design a 70th order lowpass FIR filter with cutoff frequency of 75 Hz.
Fs = 1000;                    % sample rate in Hz
Fnorm = cutoff_frequency/(Fs/2);           % Normalized frequency
df = designfilt('lowpassfir','FilterOrder',70,'CutoffFrequency',Fnorm);


%Create data file
createddata_dir = ['../created_data/'];
if ~exist([createddata_dir,id],'dir')
 mkdir([createddata_dir,id]);
end

mkdir([createddata_dir,id,'/data']);
mkdir([createddata_dir,id,'/csv']);

datadir = [createddata_dir,id,'/data/'];


% Directory names for data.
wavdir = [datadir,ID, '_Wave'];
firingsdir = [datadir,ID, '_Firings'];
workspacedir = [datadir,ID, '_Workspace'];
setdir = ['setting'];




% Error Checking.
if(any(isspace(ID)))
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

if ~exist(setdir, 'dir')
    mkdir(setdir);
else
    addpath(setdir);
end

% Creating workspace names.

    workspaceFilename=[workspacedir,'/babble_daspnet_reservoir_',ID,'.mat'];


    % Creating import names.

    %importFilename=[setdir,'/initial.mat'];
    importFilename=[setdir,'/babble_daspnet_reservoir_randominitial.mat']; %sparse random data



% Directory for Coath et. al. Saliency Detector.
addpath('auditorysaliencymodel');


%Import initial value
table=importdata([setdir,'/table_notnormalization_from0.5to1.5KHz.csv']);

if strcmp(separatephase,'nseparate')
    
    %load(importFilename,'s','sout','post','post_spe','post_mot','pre','aux');
    
    
elseif strcmp(separatephase,'separate')
    
    if strcmp(yoke,'No')
        load([setdir,'/babble_daspnet_reservoir_1_171201_3000_reinforce_100_4_No_1_1_0.03_0.3_NSTD_IP_separatephase_lattice_negativereward.mat'],'s','sout','post','post_spe','post_mot','pre','aux','TE','threshold','TEI','NeuronID_Position','InputneuronID','OutputneuronID','Post_position','Post','size_neuron','sesum');
    elseif strcmp(yoke,'Sc')
        load([setdir,'/babble_daspnet_reservoir_1_171201_3000_reinforce_100_4_Sc_1_1_0.03_0.3_NSTD_IP_separatephase_lattice_negativereward.mat'],'s','sout','post','post_spe','post_mot','pre','aux','TE','threshold','TEI','NeuronID_Position','InputneuronID','OutputneuronID','Post_position','Post','size_neuron','sesum');
    end
    
elseif strcmp(separatephase,'randSc') && strcmp(yoke,'Sc')
    NoID = [num2str(itenumber),'_',id,'_',num2str(newT),'_reinforce_100_4_No_',num2str(plotOn),'_',num2str(feedbacktime),'_',num2str(learningratio),'_',num2str(speinplate),'_',STDP,'_',IP,'_',separatephase,'_',Network,'_',reward,'_',feedbacktype];
    load([datadir,NoID, '_Workspace/',num2str(iterate),'_babble_daspnet_reservoir_',NoID,'.mat'],'muscle_his','salhist');
    randtime=randperm(newT);
    Nomuscle_his=muscle_his(:,randtime);
    muscle_his = 0;
end


if exist(workspaceFilename) > 0
    load(workspaceFilename);
else





    M=100;                 % number of synapses per neuron
    D=1;                   % maximal conduction delay
    % excitatory neurons   % inhibitory neurons      % total number
    Ne=800;                Ni=200;                   N=Ne+Ni;
    Nout = length(outInd);
    Nmot=Nout; % Number of motor neurons that the output neurons in the reservoir connect to.
    Ninp=Nout; % Number of input neurons
    a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];     % Sets time scales of membrane recovery variable.
    d=[   8*ones(Ne,1);    2*ones(Ni,1)];       % Membrane recovery variable after-spike shift.
    a_mot=.02*ones(Nmot,1);
    d_mot=8*ones(Nmot,1); %
    s=[6*rand(Ne,M);-5*rand(Ni,M)];         % synaptic weights (default of walaumont is s=[rand(Ne,M);-rand(Ni,M)]
    sout=rand(Nout,Nmot); % Synaptic weights from the reservoir output neurons to the motor neurons.
    sd=zeros(N,M);                          % their derivatives

    post_mot=repmat(1:Nmot,Nout,1);         % All output neurons connect to all motor neurons.

    % Normalizing the synaptic weights.
    %sout=sout./(mean(mean(sout)));
    %sinp=sinp./(mean(mean(sinp)));
    %s(1:Ne,:)=s(1:Ne,:)./(mean(mean(s(1:Ne,:)))); %seikika

    sd_mot=zeros(Nout,Nmot); % The change to be made to sout. %STDP
    
    if strcmp(IP,'LiIP')
        D_noise=0.5;  %gausian noise coefficient
        b_max=0.2;
        b_min=0.12;
        b=[b_min*ones(Ne,1)+(b_max-b_min)*rand(Ne,1);b_min*ones(Ni,1)+(b_max-b_min)*rand(Ni,1)];       %default b=[b_min*ones(Ne,1)+(b_max-b_min)*rand(Ne,1);0.2*ones(Ni,1)]
        Te_max=110;     %default 110 ,lognormal 100000
        Te_min=90;      %default 90,lognormal 25
        Ti_max=25;%25   lognormal 100000
        Ti_min=15;%7     lognormal 25
        h=[0.012*ones(Ne,1);0.012*ones(Ni,1)];      %larning rate of z
        z=zeros(N,1);       %larning rate of fai and b
        all_fired=[];
        fired_time=zeros(N,1);
        ISI=zeros(N,1);
        
    elseif strcmp(IP,'NoIP')
        D_noise=0.5;  %gausian noise coefficient
        b_max=0.2;
        b_min=0.12;
        b=[b_min*ones(Ne,1)+(b_max-b_min)*rand(Ne,1);b_min*ones(Ni,1)+(b_max-b_min)*rand(Ni,1)];       %default b=[b_min*ones(Ne,1)+(b_max-b_min)*rand(Ne,1);0.2*ones(Ni,1)]

    elseif strcmp(IP,'threIP')
        TE.max = -15;             %for IP
        TI.max = -15.2;             %for IP
        etaIP = 0.005;          %for IP 0.01
        threshold = -55.1;          %for IP  -55.1
        sumweight =5000;          %for SN(amount of wight of output-motor)
        %HIP = 1/200  ;  %target firing rate (100 = numer of input neuron) defalt 2*input/Ne
        %%%%%%%%%%%%%%%%%%%%
        m_IP = ones(100,1)/200;

        %{
        %normal distribution
        variance_IP = ones(1000,1)/200;
        m_IP = ones(100,1)/200;
        V_HIP = normrnd(variance_IP,0.001);
        m_HIP = normrnd(m_IP,0.001);
        %}
  
        %lognormal distribution
        m = 1/100;
        v = 0.0001;
        mu = log((m^2)/sqrt(v+m^2));
        sigma = sqrt(log(v/(m^2)+1));
        V_HIP = lognrnd(mu,sigma,1000,1);
        %m_HIP = lognrnd(mu,sigma,100,1);      %lonnomal of motor neuron
        m_HIP = normrnd(m_IP,0.001);
        %for debag
        %debag_HIP = V_HIP*Fs;
        %edge = logspace(0, 10, 300);
        %h=histogram(debag_HIP,edge);set(gca,'Xscale','log');xlim([0 100]);
        %%%%%%%%%%%%%%%%%%%%%%

    end


    if strcmp(STDP,'STDP')
    sd_reservor=zeros(N,M); %inter reservor STDP
    end

    for i=1:N
        delays{i,1}=1:M;
    end
    for i=1:Nout
        delays_mot{i,1}=1:Nmot;
    end

    STDP_out = zeros(Nout,1001+D); % out motor STDP matrix
    STDP_mot = zeros(Nout,1001+D); % spe inp STDP matrix


    if strcmp(STDP,'STDP')
     for i=1:N
           % pre{i}=find(post==i&s>0); %find pre excitatory neurons
           % aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
     end;

     STDP_reservor = zeros(N,1001+D); % inter reservor STDP matrix
    end


    if strcmp(IP,'threIP')
        
        TE.r = TE.max*(rand(Ne,1));
        TI = TI.max*(rand(Ni,1));
        TEI=[TE.r;TI]; 
        
        TE.m = TE.max*(rand(Nout,1));
        
        %%%%%%% SN
        sesum = sum(sum(s(1:Ne,:)));
        %smsum = sum(sum(sout));
%        s = s/sum(sum(s));
         %sout = sumweight*sout/sum(sum(sout));
         
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%% create random network
    if strcmp(Network,'random')
        %シナプスのつなぎ方を決定する
        post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
        for i=1:Ne
            while length(find(post(i,1:M)==i)) >=1
                numbering=find(post(i,1:M)==i);
                for j=1:length(find(post(i,1:M)==i))
                    post(i,numbering(j))=N*ceil(rand);
                end
            end
        end
        
        
        for i=1:N
            if i<=Ne
                for j=1:D
                    delays{i,j}=M/D*(j-1)+(1:M/D);
                end
            else
                delays{i,1}=1:M;
            end
            pre{i}=find(post==i&s>0);             % pre excitatory neurons
            aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%% create lattice network
    if strcmp(Network,'lattice') && strcmp(separatephase,'nseparate')
        %initialization
        Position_xy = zeros(1000,2);        %xy date of neuron position
        Post_position = zeros(1000,100);        %matrix of post neuron (data is position not neuron ID)
        Post.lattice = zeros(1000,100);      %matrix of post neuron(vertical= neuron ID , horizontal = post neuron ID)

        in_nrows = 100;
        in_ncols = 8;
        nrows = 100;
        ncols = 2;
        internal_neuron = reshape(randperm(in_ncols*in_nrows), [in_nrows in_ncols])+200; % matrix of random neuron Position 100*10
        inputoutput = reshape(randperm(ncols*nrows), [nrows ncols]);
        InputneuronID = inputoutput(:,1);
        OutputneuronID = inputoutput(:,end);
        
        NeuronID_Position = [InputneuronID,internal_neuron,OutputneuronID];         % matrix of random neuron Position 100*10
        
        InputneuronID2 = NeuronID_Position(:,2);

        r = 1*sqrt(2);      % range of connection (for example r = 1, connect next to a neuron ) 


        %create Position information
        k = 1;
        for i = 1:10
            for j = 1:100

                Position_xy(k,:) = [i,j];
                k = k + 1;
            end
        end

        %caliculate the distance and nearly neuron ID
        [idx, dist] = rangesearch(Position_xy,Position_xy,r);

        %caliculate Post neuron (position data)
        for i = 1: 1000
            size_idx = size(idx{i});
            Post_position(i,1:size_idx(1,2)) = idx{i};
        end
        %delete number of themselve
        Post_position(:,1) = 0;

        %caliculate Post neuron (neuron ID data)
        for i = 1:1000
        size_neuron{i} = 2:max(find(Post_position(find(NeuronID_Position==i),:)~=0));
        Post.lattice(i,size_neuron{i}) = NeuronID_Position(Post_position(find(NeuronID_Position==i),size_neuron{i})); %this sentense create new post matrix(neuron ID)
        end
        
        
        %create pre and aux for STDP
        
        for i=1:N
            
            pre{i}=find(Post.lattice==i&s>0);             % pre excitatory neuron
            aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%% end of create lattice network
    
    
    if strcmp(Network,'random') && ~strcmp(separatephase,'separate')
        
        random.a = size(find(post(1,:)~=0));
        
        
        
        
    end
    
    
    
    
    
    v = -65*ones(N,1);          % Membrane potentials.
    v_mot = -65*ones(Nmot,1);   %
    u = 0.2.*v;                 % Membrane recovery variable.
    u_mot = 0.2.*v_mot;
    firings=[-D 0];       % All reservoir neuron firings for the current second.
    outFirings=[-D 0];  % Output neuron spike timings.
    motFirings=[-D 0];  % Motor neuron spike timings.
    inpFirings=[-D 0];  % Input neuron spike timings.

    DA=0; % Level of dopamine above the baseline.

    muscsmooth=100; % Spike train data sent to Praat is smoothed by doing a 100 ms moving average.
    sec=0;

    rewcount=0;
    rew=[];

    % History variables.
    v_mot_hist={};
    sout_hist={};
    v_spe_hist={};
    sinp_hist={};
    feedbackhist=zeros(100,1);
    
    
    % sparse out-motor synapse
    if sparse_mot == 1
        
        for i = 1:Nout
            zero_mot{i} = randperm(100,sparse_degree);
        end

    end
    

    % Initializing reward policy variables.
        %strcmp
        temprewhist=zeros(1,10); % Keeps track of rewards given at a threshold value for up to 10 previous sounds.

end


T=newT;
clearvars newT;

% Absolute path where Praat can be found.
praatPathpc = 'c:\users\Takimto\Praat\Praat';
praatPathmac = '/Applications/Praat.app/Contents/MacOS/Praat';
praatPathlinux = '/usr/bin/praat';


datahistsize=((1000-muscsmooth)/feedbacktime);



%RUNNING THE SIMULATION%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sec=(sec+1):T % T is the duration of the simulation in seconds.
          tic;
    display('********************************************');
    display(['iterate=',num2str(itenumber),'yoke=',yoke,'_phase=',STDP,'_',IP,'_Second ',num2str(sec),' of ',num2str(T)]);

    v_mot_hist{sec}=[]; % Record of all the membrane voltages of the motor neurons.


    %Initialize decaysmooth
    decaysmoothneg=zeros(1,1000);
    decaysmoothpos=zeros(1,1000);
    datatime=1;

    for t=1:1000                         % Millisecond timesteps

        %Random Thalamic Input or No random input (IP)
        if strcmp(IP,'Tonic')
            I=13*(rand(N,1)-0.5);
            I_mot=13*(rand(Nmot,1)-0.5);
        elseif strcmp(IP,'threIP')
            I=zeros(N,1);
            I_mot=zeros(Nmot,1);
        elseif strcmp(IP,'afterIP')
            I=zeros(N,1);
            I_mot=zeros(Nmot,1);
        elseif strcmp(IP,'LiIP')
            I=[6*ones(Ne,1);0*ones(Ni,1)];      %default I=[6*ones(Ne,1);0*ones(Ni,1)]
            I_mot=zeros(Nmot,1);
            %I_mot=12*(rand(Nmot,1)-0.5);
        elseif strcmp(IP,'NoIP')
            I=[6*ones(Ne,1);0*ones(Ni,1)];
            I_mot=zeros(Nmot,1);
        end
        
        %for ip test
        %I(1:Ninp)=13*(rand(Nmot,1)-0.5);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%   feedback every time
        if strcmp(Network,'random')
           if t-1>muscsmooth&mod(t-1,feedbacktime)==0
               
               

            
            if strcmp(yoke,'No') && strcmp(feedbacktype,'fft')
                feedback1=speinplate*table(:,muscle_number);
                I(1:Ninp)=I(1:Ninp)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                I((Ninp*2+1):(Ninp*2+Ninp))=I((Ninp*2+1):(Ninp*2+Ninp))+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,feedback1];
            elseif strcmp(yoke,'Sc') && strcmp(feedbacktype,'fft') && ~strcmp(separatephase,'randSc')
                feedback1=speinplate*table(:,muscle_number);
                randid=randperm(100);
                yokedfeedback=feedback1(randid);
                I(1:Ninp)=I(1:Ninp)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                I((Ninp*2+1):(Ninp*2+Ninp))=I((Ninp*2+1):(Ninp*2+Ninp))+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,yokedfeedback];
                
                
            elseif strcmp(yoke,'Sc') && strcmp(feedbacktype,'fft') && strcmp(separatephase,'randSc')
                feedback1=speinplate*table(:,Nomuscle_his(t,sec));
                if strcmp(yokemode,'timespace')
                    randid=randperm(100);
                    yokedfeedback=feedback1(randid);
                else
                    yokedfeedback = feedback1;
                end
                
                I(1:Ninp)=I(1:Ninp)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                I((Ninp*2+1):(Ninp*2+Ninp))=I((Ninp*2+1):(Ninp*2+Ninp))+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,yokedfeedback];
                
                
            elseif strcmp(yoke,'No') && strcmp(feedbacktype,'consonant')
                
                if muscle_number>16000
                    feedback1 = constant_inplate*ones(100,1);
                else
                    feedback1 = zeros(100,1);
                end
                
                I(1:Ninp)=I(1:Ninp)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                I((Ninp+1):(Ninp+Ninp))=I((Ninp+1):(Ninp+Ninp))+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,feedback1];
                
            elseif strcmp(yoke,'Sc') && strcmp(feedbacktype,'consonant')
                              
                if rand>0.9
                    yokedfeedback = constant_inplate*ones(100,1);
                else
                    yokedfeedback = zeros(100,1);
                end
                I(1:Ninp)=I(1:Ninp)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                I((Ninp+1):(Ninp+Ninp))=I((Ninp+1):(Ninp+Ninp))+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,yokedfeedback];
            elseif strcmp(feedbacktype,'none')
            end

           end
        elseif strcmp(Network,'lattice')
            
           if t-1>muscsmooth&mod(t-1,feedbacktime)==0

            if strcmp(yoke,'No') && strcmp(feedbacktype,'fft')
                feedback1=speinplate*table(:,muscle_number);
                I(InputneuronID)=I(InputneuronID)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                I(InputneuronID2)=I(InputneuronID2)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,feedback1];
            elseif strcmp(yoke,'Sc') && strcmp(feedbacktype,'fft')
                feedback1=speinplate*table(:,muscle_number);
                randid=randperm(100);
                yokedfeedback=feedback1(randid);
                I(InputneuronID)=I(InputneuronID)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                I(InputneuronID2)=I(InputneuronID2)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,yokedfeedback];
                
            elseif strcmp(yoke,'No') && strcmp(feedbacktype,'consonant')
                if muscle_number>16000
                    feedback1 = constant_inplate*ones(100,1);
                else
                    feedback1 = zeros(100,1);
                end
                I(InputneuronID)=I(InputneuronID)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                I(InputneuronID2)=I(InputneuronID2)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,feedback1];
                
            elseif strcmp(yoke,'Sc') && strcmp(feedbacktype,'consonant')
                              
                if rand>0.9
                    yokedfeedback = constant_inplate*ones(100,1);
                else
                    yokedfeedback = zeros(100,1);
                end
                I(InputneuronID)=I(InputneuronID)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                I(InputneuronID2)=I(InputneuronID2)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,yokedfeedback];
            elseif strcmp(feedbacktype,'none')
                

            end
            

           end
        end
            

        %%%%%%%%%%%%%%%%%%%%%

        if strcmp(IP,'Tonic')
            
            fired = find(v>=30);                % Indices of fired neurons
            fired_out = find(v(Ninp+outInd)>=30);    %100~200
            fired_mot = find(v_mot>=30);        %motorneurons
            fired_inp = find(v(inpInd)>=30);
        elseif strcmp(IP,'threIP')

            fired = find((v-TEI)>=threshold);
            fired_out = find((v(Ninp + outInd)-TEI(Ninp + outInd))>=threshold);
            fired_inp = find((v(inpInd) - TEI(inpInd))>=threshold);
            
            %fired_mot = find((v_mot-TE.m)>=threshold);     %IP for motor
            fired_mot = find((v_mot)>=mot_thre);      %not IP for motor neuron
            
            
            %%%%%%%%% SN
            s_i = s(801:1000,:);
            s_e = sesum*s(1:Ne,:)/sum(sum(s(1:Ne,:)));
            %s_e = s(1:Ne,:)/sum(sum(s(1:Ne,:)));%実験
            s = [s_e ; s_i];
            %sout = smsum*sout/sum(sum(sout));
            %sout = sumweight*sout/sum(sum(sout));%実験
            
            
        elseif strcmp(IP,'afterIP')
            
            fired = find((v-TEI)>=threshold);
            fired_out = find((v(Ninp + outInd)-TEI(Ninp + outInd))>=threshold);
            fired_inp = find((v(inpInd) - TEI(inpInd))>=threshold);
            
            %fired_mot = find((v_mot-TE.m)>=threshold);
            fired_mot = find((v_mot)>=mot_thre);      %not IP for motor neuron
            
            %%%%%%%%% SN
            s_i = s(801:1000,:);
            s_e = sesum*s(1:Ne,:)/sum(sum(s(1:Ne,:)));
            %s_e = s(1:Ne,:)/sum(sum(s(1:Ne,:)));%実験
            s = [s_e ; s_i];
            %sout = smsum*sout/sum(sum(sout));
            %sout = sumweight*sout/sum(sum(sout));%実験
            
        elseif strcmp(IP,'LiIP')
            fired = find(v>=30);                % Indices of fired neurons
            fired_out = find(v(Ninp+outInd)>=30);    %100~200
            fired_mot = find(v_mot>=30);        %motorneurons
            fired_inp = find(v(inpInd)>=30);
            
          
            ISI(fired)=t+(sec-1)*1000-fired_time(fired);
            fired_time(fired)=t+(sec-1)*1000;
            
        elseif strcmp(IP,'NoIP')
            fired = find(v>=30);                % Indices of fired neurons
            fired_out = find(v(Ninp+outInd)>=30);    %100~200
            fired_mot = find(v_mot>=30);        %motorneurons
            fired_inp = find(v(inpInd)>=30);
         
        end
        
        


        v(fired)=-65;                       % Reset the voltages for those neurons that fired
        v_mot(fired_mot)=-65;
        u(fired)=u(fired)+d(fired);         % Individual neuronal dynamics
        u_mot(fired_mot)=u_mot(fired_mot)+d_mot(fired_mot);

        % Spike-timing dependent plasticity computations:
        STDP_out(fired_out,t+D)=0.1; % Keep a record of when the output neurons spiked.
        STDP_mot(fired_mot,t+D)=0.1; % Keep a record of when the spectrum neurons spiked.
        if strcmp(STDP,'STDP')
        STDP_reservor(fired,t+D)=0.1; %inter reservor STDPmatrix record


         for k=1:length(fired)
                sd_reservor(pre{fired(k)})=sd_reservor(pre{fired(k)})+LTP*STDP_reservor(N*t+aux{fired(k)}); %LTP A plus 1
                % pre{fired(1)=54}=52,78,210,277,350,372,73
               
                
         end;


        end


        for k=1:length(fired_mot)
                sd_mot(:,fired_mot(k))=sd_mot(:,fired_mot(k))+LTP_m*STDP_out(:,t); % Adjusting sd for synapses eligible for potentiation.LTP
        end
   

        firings=[firings;t*ones(length(fired),1),fired];                % Update the record of when neuronal firings occurred.
        outFirings=[outFirings;t*ones(length(fired_out),1),fired_out];  %
        motFirings=[motFirings;t*ones(length(fired_mot),1),fired_mot];  %
        inpFirings=[inpFirings;t*ones(length(fired_inp),1),fired_inp];  %


        % For any presynaptic neuron that just fired, calculate the current to add
        % as proportional to the synaptic strengths from its postsynaptic neurons.
        if strcmp(Network,'random')
            k=size(firings,1);
            while firings(k,1)>t-D
            %del=delays{firings(k,2),t-firings(k,1)+1}; %1:M synaptic data
            del = 1:random.a(1,2);
            ind = post(firings(k,2),del);          % post neurons of fired neron id [23,45,13,14,53,??????]
                I(ind)=I(ind)+s(firings(k,2), del)';%'
             if strcmp(STDP,'STDP')

              sd_reservor(firings(k,2),del)=sd_reservor(firings(k,2),del)-LTD*STDP_reservor(ind,t+D)';%'LTD A-1.5
             end

                k=k-1;
            end;
        end
        
        if strcmp(Network,'lattice')
            k=size(firings,1);
            while firings(k,1)>t-D
            %del=delays{firings(k,2),t-firings(k,1)+1}; %1:M  synaptic data
            %ind = post(firings(k,2),del);          % post neurons of fired neron id [23,45,13,14,53,??????]
            
            del = size_neuron{firings(k,2)};
            ind = Post.lattice(firings(k,2),del);
        
                I(ind)=I(ind)+s(firings(k,2), del)';%'
             if strcmp(STDP,'STDP')

              sd_reservor(firings(k,2),del)=sd_reservor(firings(k,2),del)-LTD*STDP_reservor(ind,t+D)';
              %'LTD A-1.5,firings(k,2) = pre firing,ind = post neuron,STDP_reservor(ind,t+D)=firing history of post neuron(pre at t) at t+Dlater 
             end

                k=k-1;
            end
        end
        
        
        
        
        % Calculating currents to add for motor neurons.
        k=size(outFirings,1);
        while outFirings(k,1)>t-D       %D = 1
            del_mot=delays_mot{outFirings(k,2),t-outFirings(k,1)+1};
            ind_mot = post_mot(outFirings(k,2),del_mot);   %del_mot=1:Nmot ind_mot=
            I_mot(ind_mot)=I_mot(ind_mot)+2*sout(outFirings(k,2), del_mot)';%'
            
            sd_mot(outFirings(k,2),:)=sd_mot(outFirings(k,2),:)-LTD_m*STDP_mot(:,t+D)'; % LTD of motor neuron 
            
            k=k-1;
        end;
        
        
        %caliculate Li et al 17 (7)
        if strcmp(IP,'LiIP') && mod(t,50)==0
           for i=1:Ne
            if ISI(i)<Te_min
                z(i)=-h(i)*exp((Te_min-ISI(i))/Te_min);
            elseif ISI(i)>Te_max
                z(i)=h(i)*exp((ISI(i)-Te_max)/Te_max);
            elseif ISI(i)>=Te_min && ISI(i)<=Te_max
                z(i)=0;
            end
           end
           
           for i=Ne+1:N
                if ISI(i)<Ti_min
                    z(i)=-h(i)*exp((Ti_min-ISI(i))/Ti_min);
                elseif ISI(i)>Ti_max
                    z(i)=h(i)*exp((ISI(i)-Ti_max)/Ti_max);
                elseif ISI(i)>=Ti_min && ISI(i)<=Ti_max
                    z(i)=0;
               end
           end
           
        b=b+b_max*z;
        for i=1:N
            b(i)=max(0.12,b(i));
            b(i)=min(0.2,b(i));
        end
        
        end
        


        

%%%%%%%%%%%%%%%
        %caliculate LST
        LAP=mean(mean(v));
        LAP_hist(t)=LAP;
%%%%%%%%%%%%%%%


        % Individual neuronal dynamics computations:
        v=v+0.5*((0.04*v+5).*v+140-u+I);                            % for numerical I...randam+from spe neurons(1*100)
        v=v+0.5*((0.04*v+5).*v+140-u+I);                            % stability time
        v_mot=v_mot+0.5*((0.04*v_mot+5).*v_mot+140-u_mot+I_mot);    % step is 0.5 ms I_mot...randam+from out neurons
        v_mot=v_mot+0.5*((0.04*v_mot+5).*v_mot+140-u_mot+I_mot);
        if strcmp(IP,'LiIP') || strcmp(IP,'NoIP')
            u=u+a.*(b.*v-u)+D_noise*randn(N,1); % stability,gausian noise
        else
            u=u+a.*(0.2*v-u);
        end
        u_mot=u_mot+a_mot.*(0.2*v_mot-u_mot);



                % Exponential decay of the traces of presynaptic neuron firing
        STDP_out(:,t+D+1)=0.95*STDP_out(:,t+D);                             % tau = 20 ms
        STDP_mot(:,t+D+1)=0.95*STDP_mot(:,t+D);                             % tau = 20 ms

        if strcmp(STDP,'STDP')

          STDP_reservor(:,t+D+1)=0.95*STDP_reservor(:,t+D);

        end

                % Exponential decay of the dopamine concentration over time.
        DA=DA*0.995;

        
      
        if sparse_mot == 1
            
            for i = 1:Nout
                sd_mot(i,zero_mot{i}) = 0;
                sout(i,zero_mot{i}) = 0;
            end
            
        end
  
        
        % Modify synaptic weights.
        if (mod(t,10)==0)
            if strcmp(reinforce,'reinforce')
                %sout=0.5*sout./(mean(mean(sout)));              % Normalizing the synaptic weights. out-mot
                if sec>50
                sout=max(0,min(sm,sout+DA*sd_mot));
                end
                
            end

            %sinp=max(0,min(sm,sinp+sd_spe));%
            %sinp=sinp./((mean(mean(sinp)))*speinplate);            % Normalizing the synaptic weights.  spe-inp

            if strcmp(STDP,'STDP')
             s(1:Ne,:)=max(0,min(smr,s(1:Ne,:)+(learningratio)*sd_reservor(1:Ne,:)));%reservor
             
             %s(1:Ne,:)=s(1:Ne,:)./(mean(mean(s(1:Ne,:)))); %seikika
             sd_reservor=0.99*sd_reservor;
            end



            sd_mot=0.99*sd_mot; % Decrease in synapse's' eligibility to change over time.
    %       sd_spe=0.09*sd_spe;

        end;   %
        
        
        % IP lerning
        if strcmp(IP,'threIP')
            
            fire_all = zeros(N,1);
            fire_m = zeros(Nmot,1);
            fire_all(fired) = 1;     %fired neuron of RSN
            fire_m(fired_mot) = 1;
            
            %IP only excitatory neurons
            TEI(:,1) = TEI(:,1) + etaIP*(fire_all(:,1) - V_HIP(:,1));
            %TEI(1:Ne,1) = TEI(1:Ne,1) + etaIP*(fire_all(1:Ne,1) - V_HIP(1:Ne,1));
            %TEI(1:Ne,1) = TEI(1:Ne,1) + etaIP*(fire_all(1:Ne,1) - HIP);
            %IP all neurons
            %TEI(:,1) = TEI(:,1) + etaIP*(fire_all(:,1) - HIP);
            
            TE.m = TE.m + etaIP*(fire_m - m_HIP);
            %TE.m = TE.m + etaIP*(fire_m - HIP);  %IP
           
        elseif strcmp(IP,'LiIP')
            
            
            
        end
        
        
        
        
        % Every testint seconds, use the motor neuron spikes to generate a sound.
        if (mod(sec,testint)==0)

            if strcmp(IP,'Tonic')
                firedmusc1pos=find(v_mot(1:Nmot/2)>=30); % Find out which of the jaw/lip motor neurons fired.
                firedmusc1neg=find(v_mot(Nmot/2+1:end)>=30);
                
            elseif strcmp(IP,'threIP')
                %firedmusc1pos=find(v_mot(1:Nmot/2)-TE.m(1:Nmot/2)>=threshold); % Find out which of the jaw/lip motor neurons fired.
                %firedmusc1neg=find(v_mot(Nmot/2+1:end)-TE.m(Nmot/2+1:end)>=threshold);
                
                firedmusc1pos=find(v_mot(1:Nmot/2)>=mot_thre); % Find out which of the jaw/lip motor neurons fired.
                firedmusc1neg=find(v_mot(Nmot/2+1:end)>=mot_thre);
                
            elseif strcmp(IP,'afterIP')
                firedmusc1pos=find(v_mot(1:Nmot/2)-TE.m(1:Nmot/2)>=threshold); % Find out which of the jaw/lip motor neurons fired.
                firedmusc1neg=find(v_mot(Nmot/2+1:end)-TE.m(Nmot/2+1:end)>=threshold);
                
            elseif strcmp(IP,'LiIP') || strcmp(IP,'NoIP')
                firedmusc1pos=find(v_mot(1:Nmot/2)>=30); % Find out which of the jaw/lip motor neurons fired.
                firedmusc1neg=find(v_mot(Nmot/2+1:end)>=30);
                
            end
            
            summusc1posspikes(t)=size(firedmusc1pos,1); % Sum the spikes at each timestep across the set of motor neurons.
            summusc1negspikes(t)=size(firedmusc1neg,1);

            %create decaysmoothmusc
            %muscsmooth is 100
            if(t>muscsmooth)
                for i=1:muscsmooth
                    decaysmoothneg(t)=decaysmoothneg(t)+summusc1negspikes(t-muscsmooth+i)*(1-exp(-i/tau));
                    decaysmoothpos(t)=decaysmoothpos(t)+summusc1posspikes(t-muscsmooth+i)*(1-exp(-i/tau));
                end
                decaysmoothneg(t)=decaysmoothneg(t)/(muscsmooth*(1-tau/100));
                decaysmoothpos(t)=decaysmoothpos(t)/(muscsmooth*(1-tau/100));

              smoothmusc(t)=muscscale*(decaysmoothpos(t)-decaysmoothneg(t));  %positive-negative
              if smoothmusc(t)>1
                  smoothmusc(t)=1;
              elseif smoothmusc(t)<-0.9999
                  smoothmusc(t)=-0.9999;
              end


              muscle_number=round((smoothmusc(t)+1)*10000);
              
              muscle_his(t,sec) = muscle_number;
              
              

            end
            
            

            


            if t==1000 % Based on the 1 s timeseries of smoothed summed motor neuron spikes, generate a sound.


                if ~strcmp(reinforce,'range')
                    
                    %low pass filter for negative reward
                    if strcmp(reward,'nega')
                        
                    filter_smoothmusc = filter(df,smoothmusc);      %filter smoothmusc
                  
                    for i=2:1000        %caliculate muscle cost
                        
                        negativereward = negativereward + sqrt((filter_smoothmusc(i)-filter_smoothmusc(i-1))^2);
                    end
                    
                    end
                    
                    
                    
                  
                    
                    % Write the Praat script:

                    fid = fopen([wavdir,'/ressynth_',ID,num2str(sec,'%d'),'.praat'],'w');
                    fprintf(fid,'Create Speaker... speaker Female 2\n');
                    fprintf(fid,['Create Artword... babble ' num2str((1000-muscsmooth)/1000,'%.3f') '\n']);
                    fprintf(fid,'select Artword babble\n');
                    fprintf(fid,['Set target... ' num2str((0)/1000,'%.3f') ' ' num2str(0.1,'%.3f') ' Lungs\n']);
                    fprintf(fid,['Set target... ' num2str(0.02,'%.3f') ' ' num2str(0.1,'%.3f') ' Lungs\n']);
                    fprintf(fid,['Set target... ' num2str(0.05,'%.3f') ' ' num2str(0,'%.3f') ' Lungs\n']);
                    fprintf(fid,['Set target... ' num2str((1000-muscsmooth)/1000,'%.3f') ' ' num2str(0,'%.3f') ' Lungs\n']);
                    fprintf(fid,['Set target... ' num2str((0)/1000,'%.3f') ' ' num2str(0.5,'%.3f') ' Interarytenoid\n']);
                    fprintf(fid,['Set target... ' num2str((1000-muscsmooth)/1000,'%.3f') ' ' num2str(0.5,'%.3f') ' Interarytenoid\n']);
                    fprintf(fid,['Set target... ' num2str((0)/1000,'%.3f') ' ' num2str(0.4,'%.3f') ' ' ' Hyoglossus\n']);
                    fprintf(fid,['Set target... ' num2str((1000-muscsmooth)/1000,'%.3f') ' ' num2str(0.4,'%.3f') ' Hyoglossus\n']);
                    for praatt=0:(1000-muscsmooth)
                        fprintf(fid,['Set target... ' num2str((praatt)/1000,'%.3f') ' ' num2str(smoothmusc(praatt+muscsmooth),'%.3f') ' Masseter\n']);
                        fprintf(fid,['Set target... ' num2str((praatt)/1000,'%.3f') ' ' num2str(smoothmusc(praatt+muscsmooth),'%.3f') ' OrbicularisOris\n']);
                    end
                    fprintf(fid,'select Speaker speaker\n');
                    fprintf(fid,'plus Artword babble\n');
                    fprintf(fid,'To Sound... 22050 25 0 0 0 0 0 0 0 0 0\n');
                    fprintf(fid,'\tselect Sound babble_speaker\n');

                        fprintf(fid,['\tWrite to WAV file... synth_',ID,'_',num2str(sec),'.wav\n']);


                    fclose(fid);

                   % Execute the Praat script -- produces a wave file:
                    % Housekeeping
                    if ismac
                            system([praatPathmac, ' ', wavdir,'/ressynth_',ID,num2str(sec,'%d'),'.praat']);
                            delete([wavdir,'/ressynth_',ID,num2str(sec,'%d'),'.praat']);  %delete praat script
                    elseif ispc
                            system([praatPathpc, ' --run ', wavdir,'\ressynth_',ID,num2str(sec,'%d'),'.praat']);
                            delete([wavdir,'\ressynth_',ID,num2str(sec,'%d'),'.praat']);  %delete praat script
                            
                    elseif isunix
                        
                        system([praatPathlinux, ' --run ', wavdir,'/ressynth_',ID,num2str(sec,'%d'),'.praat']);
                        delete([wavdir,'/ressynth_',ID,num2str(sec,'%d'),'.praat']);  %delete praat script

                    end



                    % Find the auditory salience of the sound:
                        salienceResults = auditorySalience([wavdir,'/synth_',ID,'_',num2str(sec),'.wav'],0);

                    salience = sum(abs(salienceResults.saliency(31:180))); % Summing over salience trace to produce a single value.
                                                                           % abs(X)= return absolute value    neglect first 30  saliency=[1,180]

                    salhist(sec,1) = salience; % History of salience over entire simulation.
                    display(['salience =',num2str(salience)]);

                    if mod(sec,1)~=0
                      delete([wavdir,'/synth_',ID,'_',num2str(sec),'.wav']);%
                    end

                end

                % Assign Reward.

%%%%%%%%%%%%%%%%%%%%%%%%


                    display(['salthresh: ',num2str(salthresh)]);
                    temprewhist(1:9)=temprewhist(2:10);
                    % Reward if the salience of the sound is above
                    % threshold value.
                    if salience>salthresh
                        display('rewarded');
                        rew=[rew,sec*1000+t];
                        rewcount=rewcount+1;
                        temprewhist(10)=1;
                        % If at least 3 of the last 10 sounds were above
                        % threshold, raise the threshold value and reset the count.
                        if sum(temprewhist==1)>=3
                            salthresh=salthresh+.1;
                            temprewhist=zeros(1,10);
                        end
                    elseif salience<salthresh
                        display('not rewarded');
                        temprewhist(10)=-1;

                        if sum(temprewhist==-1)>=10
                            salthresh=salthresh-0.1;
                            temprewhist=zeros(1,10);
                        end
                    end
                    display(['temprewhist: ',num2str(temprewhist)]);
                    display(['mean(temprewhist): ',num2str(mean(temprewhist))]);

%%%%%%%%%%%%%%%%%%%
            end
            
        end

        % If the human listener decided to reinforce (or if the yoked control
        % schedule says to reinforce), increase the dopamine concentration.
        if strcmp(reward,'normal')
            nega(sec) = 0;
            if any(rew==sec*1000+t)   % what mean?
                DA=DA+DAinc;
            end
            DA_history(sec) = DA;
        end
        
        if strcmp(reward,'nega')
            nega(sec) = nega_rate*negativereward;
            
            if any(rew == sec*1000+t)
                DA = DA + DAinc - nega(sec);       %negativereward
            end
            DA_history(sec) = DA;
            negativereward = 0;     %initialize negativereward
        end
        
       
    end
    
    
    








    % Writing reservoir neuron firings for this second to a text file.     %recoring firings matrix
    if mod(sec,SAVINTV*testint)==0 || sec==T

            firings_fid = fopen([firingsdir,'/R_firings_',ID,'_',num2str(sec),'.txt'],'w');

        for firingsrow = 1:size(firings,1)
            fprintf(firings_fid,'%i\t',sec);
            fprintf(firings_fid,'%i\t%i',firings(firingsrow,:));
            fprintf(firings_fid,'\n');
        end
        fclose(firings_fid);
        
            motfirings_fid = fopen([firingsdir,'/Mot_firings_',ID,'_',num2str(sec),'.txt'],'w');

        for motfiringsrow = 1:size(motFirings,1)
            fprintf(motfirings_fid,'%i\t',sec);
            fprintf(motfirings_fid,'%i\t%i',motFirings(motfiringsrow,:));
            fprintf(motfirings_fid,'\n');
        end
        fclose(motfirings_fid);

        % ---- plot for gif ------
        
        if sec == T && strcmp(Network,'lattice')
            mkdir([firingsdir,'/gif_sec=',num2str(sec)]);
            figure(23);
            for i = 1:1000
                colormap(gray);
                %firing_position = zeros(100,10);
                A = firings(find(firings(:,1)==i),2);
                row = zeros(1);
                col = zeros(1);
                for j=1:size(A)
                    [row(j), col(j)] = find(NeuronID_Position==A(j));
                end


%{
                %firing_position_reverse = imcomplement(firing_position);
                    fig15 = plot(col(1:end),row(1:end),'.'); % Plot the output neurons'' spikes
                    title(['Neuron Firings',num2str(i)], 'fontweight','bold');
                    axis([0 10 0 100]);
                    saveas(fig15,[id, '_Firings/gif_sec=',num2str(sec),'/time=',num2str(i),'.png']);
                    clearvars row col;
%}
            end
        end
        
        
        %HIImageConvert2GIF([id, '_Firings/gif_sec=',num2str(sec),'/*.png'], [id, '_Firings/gif_sec=',num2str(sec),'/',num2str(sec),'.gif'], 0.1) ;
        % ----- end for gif ------
    end
    
    
    
    % Make axis labels and titles for plots that are being kept.
    % ---- plot -------
    if plotOn
        hNeural = figure(103);
        set(hNeural, 'name', ['Neural Spiking for Second: ', num2str(sec)], 'numbertitle','off');
        subplot(4,1,1)
        plot(firings(:,1),firings(:,2),'.'); % Plot all the neurons'' spikes
        title('Reservoir Firings', 'fontweight','bold');
        axis([0 1000 0 N]);
        subplot(4,1,2)
        plot(outFirings(:,1),outFirings(:,2),'.'); % Plot the output neurons'' spikes
        title('Output Neuron Firings', 'fontweight','bold');
        axis([0 1000 0 Nout]);
     %   subplot(5,1,3)
     %   plot(speFirings(:,1),speFirings(:,2),'.'); % Plot the output neurons'' spikes
     %   title('Input Neuron Firings', 'fontweight','bold');
     %   axis([0 1000 0 Ninp]);
        subplot(4,1,3)
        plot(motFirings(:,1),motFirings(:,2),'.'); % Plot the motor neurons'' spikes
        title('Motor Neuron Firings', 'fontweight','bold');
        axis([0 1000 0 Nmot]);
        if strcmp(reward,'nega')
            subplot(4,1,4);
            plot(filter_smoothmusc(muscsmooth:1000)); ylim([-.5,.5]); xlim([-100,900]); % Plot the smoothed sum of motor neuron spikes 1 s timeseries
            title('Sum of Agonist/Antagonist Motor Neuron Activation', 'fontweight','bold');
        else
            subplot(4,1,4);
            plot(smoothmusc(muscsmooth:1000)); ylim([-.5,.5]); xlim([-100,900]); % Plot the smoothed sum of motor neuron spikes 1 s timeseries
            title('Sum of Agonist/Antagonist Motor Neuron Activation', 'fontweight','bold');
        end
        

        drawnow;

        hSyn = figure(113);%
        set(hSyn, 'name', ['Synaptic Strengths for Second: ', num2str(sec)], 'numbertitle','off');%
        subplot(2,1,1);

        imagesc(sout)
        set(gca,'YDir','normal')
        colorbar;
        title('Synapse Strength between Output Neurons and Motor Neurons', 'fontweight','bold');
        xlabel('postsynaptic motor neuron index', 'fontweight','bold');
        ylabel('presynaptic output neuron index', 'fontweight','bold');

        subplot(2,1,2);

        imagesc(s)
        set(gca,'YDir','normal')
        colorbar;
        title('Synapse Strength in reservoir', 'fontweight','bold');
        xlabel('Synaptic index', 'fontweight','bold');
        ylabel('Neurons index', 'fontweight','bold');
        
      if debug==1  
        fd=figure(114);
        set(fd, 'name', ['feedback ', num2str(sec)], 'numbertitle','off');%
        imagesc(feedbackhist);
        colorbar;
        drawnow;
      end
      %  subplot(3,1,3);
      %  imagesc(sinp)
      % set(gca,'YDir','normal')
      %  colorbar;
      %  title('Synapse Strength between spectrum and input', 'fontweight','bold');
      %  xlabel('postSynaptic input neuron index', 'fontweight','bold');
      %  ylabel('presynaptic spectrum neurons index', 'fontweight','bold');

    end
    
    
    if debug
        hNeural = figure(103);
        set(hNeural, 'name', ['Neural Spiking for Second: ', num2str(sec)], 'numbertitle','off');
        subplot(4,1,1)
        plot(firings(:,1),firings(:,2),'.'); % Plot all the neurons'' spikes
        title('Reservoir Firings', 'fontweight','bold');
        axis([0 1000 0 N]);
        subplot(4,1,2)
        plot(outFirings(:,1),outFirings(:,2),'.'); % Plot the output neurons'' spikes
        title('Output Neuron Firings', 'fontweight','bold');
        axis([0 1000 0 Nout]);
     %   subplot(5,1,3)
     %   plot(speFirings(:,1),speFirings(:,2),'.'); % Plot the output neurons'' spikes
     %   title('Input Neuron Firings', 'fontweight','bold');
     %   axis([0 1000 0 Ninp]);
        subplot(4,1,3)
        plot(motFirings(:,1),motFirings(:,2),'.'); % Plot the motor neurons'' spikes
        title('Motor Neuron Firings', 'fontweight','bold');
        axis([0 1000 0 Nmot]);
        if strcmp(reward,'nega')
            subplot(4,1,4);
            plot(filter_smoothmusc(muscsmooth:1000)); ylim([-.5,.5]); xlim([-100,900]); % Plot the smoothed sum of motor neuron spikes 1 s timeseries
            title('Sum of Agonist/Antagonist Motor Neuron Activation', 'fontweight','bold');
        else
            subplot(4,1,4);
            plot(smoothmusc(muscsmooth:1000)); ylim([-.5,.5]); xlim([-100,900]); % Plot the smoothed sum of motor neuron spikes 1 s timeseries
            title('Sum of Agonist/Antagonist Motor Neuron Activation', 'fontweight','bold');
        end
        

        saveas(figure(103),[firingsdir,'/Neural_Spiking_',num2str(sec),'.png']);

      %  subplot(3,1,3);
      %  imagesc(sinp)
      % set(gca,'YDir','normal')
      %  colorbar;
      %  title('Synapse Strength between spectrum and input', 'fontweight','bold');
      %  xlabel('postSynaptic input neuron index', 'fontweight','bold');
      %  ylabel('presynaptic spectrum neurons index', 'fontweight','bold');

    end
    
    % ---- end plot ------
    

    
    
    
    %for IP, caliculate average firing rate
    f_rate = size(find(firings(:,2)==500));
    firing_rate_500 = f_rate(1,1);
    %f_rate = size(firings);
    %firing_rate = (f_rate(1,1)-1)/1000;
    display(['firing rate of 500= ',num2str(firing_rate_500)]);
    %display(['neuron500 threshold = ',num2str(TEI(500,1))]);

    sout_hist{sec}=sout;
    
    %for negative reward
    %display(['nega= ',num2str(nega(sec))]);
    %display(['DA= ',num2str(DA)]);

    % Preparing STDP and firings for the following 1000 ms.
    STDP_out(:,1:D+1)=STDP_out(:,1001:1001+D);
    STDP_mot(:,1:D+1)=STDP_mot(:,1001:1001+D);
    ind = find(firings(:,1) > 1001-D);
    firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
    ind_out = find(outFirings(:,1) > 1001-D);
    outFirings=[-D 0;outFirings(ind_out,1)-1000,outFirings(ind_out,2)];
    ind_mot = find(motFirings(:,1) > 1001-D);
    motFirings=[-D 0;motFirings(ind_mot,1)-1000,motFirings(ind_mot,2)];%
    ind_inp = find(inpFirings(:,1) > 1001-D);
    inpFirings=[-D 0;inpFirings(ind_inp,1)-1000,inpFirings(ind_inp,2)];
    feedbackhist=zeros(100,1);

    if strcmp(STDP,'STDP')
       STDP_reservor(:,1:D+1)=STDP_reservor(:,1001:1001+D);%
    end



    % Every so often, save the workspace in case the simulation is interupted all data is not lost.
    if mod(sec,SAVINTV*testint)==0 || sec==T || sec==1
        display('Data Saving..... Do not exit program.');
        save([workspacedir,'/',num2str(sec),'_babble_daspnet_reservoir_',ID,'.mat'], '-regexp', '^(?!(v_mot_hist)$).');

        if plotOn==1
        %wright LAP
        figure(6);
        fig1=plot(LAP_hist);ylim([-80,-60]);
        saveas(fig1,[workspacedir,'/LAP_',num2str(sec),'.png']);
        m = length(LAP_hist);          % Window length
        y = fft(LAP_hist);           % DFT
        power = abs(y);   % Power of the DFT
        figure(7);
        fig2=plot(power);xlim([0,100]);ylim([0,100]);
        saveas(fig2,[workspacedir,'/LAP_spectrum',num2str(sec),'.png']);

        % wright weight
        csvwrite([workspacedir,'/sout_',num2str(sec),'.csv'],sout);
        csvwrite([workspacedir,'/s',num2str(sec),'.csv'],s);
        %csvwrite([workspacedir,'/sinp',sec,'.csv'],sinp);
        end
        


    end


    toc;

        if sec==T
            number=[1:T].'; %' 
            salhistdate=[salhist,number,nega.',DA_history.'];
            csvwrite([workspacedir,'/',ID,'.csv'],salhistdate);     %write salhist data and negativereward to csv
        end
        

    

        
        
    end

end




