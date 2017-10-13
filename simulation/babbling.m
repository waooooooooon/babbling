function [] = babbling(id,newT,reinforce,outInd,muscscale,yoke,plotOn,feedbacktime,learningratio,speinplate,STDP,debug,IP)
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

% Initialization.
salthresh = 4.5;            % Initial salience value for reward (used in 'relhisal' reinforcment).
DAinc = 1;                  % Amount of dopamine given during reward.
sm = 4;                     % Maximum synaptic weight.
smr=4;                     % Maximum synaptic weight for reservoir
testint = 1;                % Number of seconds between vocalizations.
fftsize=2048;             %fft????????????????????????????????????
lpcsize=8;                %LPC size
inpInd=outInd;
tau=20;                   %decay parameter
SAVINTV=100;
LST_hist=[1:1000];
muscle_number=0;
TE.max = -2;             %for IP
TI.max = 0;             %for IP
etaIP = 0.001;          %for IP
threshold = -70;          %for IP




% Directory names for data.
wavdir = [id, '_Wave'];
firingsdir = [id, '_Firings'];
workspacedir = [id, '_Workspace'];
setdir = ['setting'];


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

if ~exist(setdir, 'dir')
    mkdir(setdir);
else
    addpath(setdir);
end

% Creating workspace names.

    workspaceFilename=[workspacedir,'/babble_daspnet_reservoir_',id,'.mat'];

% Creating import names.

    importFilename=[setdir,'/initial.mat'];


% Directory for Coath et. al. Saliency Detector.
addpath('auditorysaliencymodel');


%Import initial value
table=importdata([setdir,'/table.csv']);
load(importFilename,'s','sout','post','post_spe','post_mot','pre','aux');







    M=100;                 % number of synapses per neuron
    D=1;                   % maximal conduction delay
    % excitatory neurons   % inhibitory neurons      % total number
    Ne=800;                Ni=200;                   N=Ne+Ni;
    Nout = length(outInd);
    Nmot=Nout; % Number of motor neurons that the output neurons in the reservoir connect to.
%    Nspe=Nout; % Number of input neurons for feedback
    Ninp=Nout; % Number of input neurons
    a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];     % Sets time scales of membrane recovery variable.
    d=[   8*ones(Ne,1);    2*ones(Ni,1)];       % Membrane recovery variable after-spike shift.
    a_mot=.02*ones(Nmot,1);
    d_mot=8*ones(Nmot,1); %
%   a_spe=.02*ones(Ninp,1);%
%   d_spe=8*ones(Ninp,1); %

    %post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); % Assign the postsynaptic neurons for each neuron''s synapse in the reservoir.

    %post_spe=repmat(1:Ninp,Nspe,1); %input to Neuron ID 1:Nspe
    %post_mot=repmat(1:Nmot,Nout,1); % All output neurons connect to all motor neurons. repmat(A,1,4)


    %sout=rand(Nout,Nmot); % Synaptic weights from the reservoir output neurons to the motor neurons.
    %sinp=rand(Nspe,Ninp); % Synaptic weights from the spectrum neurons to the input neurons.

    % Normalizing the synaptic weights.
    %sout=sout./(mean(mean(sout)));
    %sinp=sinp./(mean(mean(sinp)));
    %s(1:Ne,:)=s(1:Ne,:)./(mean(mean(s(1:Ne,:)))); %seikika

    sd_mot=zeros(Nout,Nmot); % The change to be made to sout. %STDP
 %   sd_spe=zeros(Ninp,Nspe);

    if strcmp(STDP,'STDP')
    sd_reservor=zeros(N,M); %inter reservor STDP
    end

    for i=1:N
        delays{i,1}=1:M;
    end
    for i=1:Nout
        delays_mot{i,1}=1:Nmot;
    end
  %  for i=1:Ninp
  %      delays_inp{i,1}=1:Nspe;
  % end
    STDP_mot = zeros(Nout,1001+D); % out motor STDP matrix
 %  STDP_spe = zeros(Nspe,1001+D); % spe inp STDP matrix


    if strcmp(STDP,'STDP')
     for i=1:N
           % pre{i}=find(post==i&s>0); %find pre excitatory neurons
           % aux{i}=N*(D-1-ceil(ceil(pre{i}/N)/(M/D)))+1+mod(pre{i}-1,N);
     end;

     STDP_reservor = zeros(N,1001+D); % inter reservor STDP matrix
    end


    if strcmp(IP,'IP')
        HIP = 1/100  ;  %target firing rate (100 = numer of input neuron) defalt 2*input/Ne
        
        TE.r = TE.max*(rand(Ne,1));
        TI = TI.max*(rand(Ni,1));
        TEI=[TE.r;TI]; 
        
        TE.m = TE.max*(rand(Nout,1));
        
        %%%%%%% SN
        s(1:Ne,:) = s(1:Ne,:)/sum(sum(s(1:Ne,:)));
        sout = sout/sum(sum(sout));
    end
    
    v = -65*ones(N,1);          % Membrane potentials.
    v_mot = -65*ones(Nmot,1);   %
 %  v_spe = -65*ones(Nspe,1);   %
    u = 0.2.*v;                 % Membrane recovery variable.
    u_mot = 0.2.*v_mot;
 %  u_spe = 0.2.*v_spe;
    firings=[-D 0];       % All reservoir neuron firings for the current second.
    outFirings=[-D 0];  % Output neuron spike timings.
    motFirings=[-D 0];  % Motor neuron spike timings.
    inpFirings=[-D 0];  % Input neuron spike timings.
 %  speFirings=[-D 0];  % Spectrum neuron spike timings.

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

    % Initializing reward policy variables.
        %strcmp
        temprewhist=zeros(1,10); % Keeps track of rewards given at a threshold value for up to 10 previous sounds.




T=newT;
clearvars newT;

% Absolute path where Praat can be found.
praatPathpc = 'c:\users\Takimto\Praat\Praat';
praatPathmac = '/Applications/Praat.app/Contents/MacOS/Praat';

datahistsize=((1000-muscsmooth)/feedbacktime);

%RUNNING THE SIMULATION%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sec=(sec+1):T % T is the duration of the simulation in seconds.
          tic;
    display('********************************************');
    display(['Second ',num2str(sec),' of ',num2str(T)]);

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
        elseif strcmp(IP,'IP')
        I=zeros(N,1);
        I_mot=zeros(Nmot,1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%   feedback every time
       if t-1>muscsmooth&mod(t-1,feedbacktime)==0

        feedback1=speinplate*table(:,muscle_number);
        if strcmp(yoke,'NY')
            I(1:Ninp)=I(1:Ninp)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
            feedbackhist=[feedbackhist,feedback1];
        elseif strcmp(yoke,'Yoked')
            randid=randperm(100);
            yokedfeedback=feedback1(randid);
            I(1:Ninp)=I(1:Ninp)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
            feedbackhist=[feedbackhist,yokedfeedback];
        end

        end

        %%%%%%%%%%%%%%%%%%%%%


        if strcmp(IP,'Tonic')
            
            fired = find(v>=30);                % Indices of fired neurons
            fired_out = find(v(Ninp+outInd)>=30);    %100~200
            fired_mot = find(v_mot>=30);        %motorneurons
   %        fired_spe = find(v_spe>=30);
            fired_inp = find(v(inpInd)>=30);
        elseif strcmp(IP,'IP')

            fired = find((v-TEI)>=threshold);
            fired_out = find((v(Ninp + outInd)-TEI(Ninp + outInd))>=threshold);
            fired_inp = find((v(inpInd) - TEI(inpInd))>=threshold);
            
            fired_mot = find((v_mot-TE.m)>=threshold);
            
            
            %%%%%%%%% SN
            s(1:Ne,:) = s(1:Ne,:)/sum(sum(s(1:Ne,:)));
            sout = sout/sum(sum(sout));
            
        end
        
        


        v(fired)=-65;                       % Reset the voltages for those neurons that fired
        v_mot(fired_mot)=-65;
   %    v_spe(fired_spe)=-65;
        u(fired)=u(fired)+d(fired);         % Individual neuronal dynamics
        u_mot(fired_mot)=u_mot(fired_mot)+d_mot(fired_mot);
   %    u_spe(fired_spe)=u_spe(fired_spe)+d_spe(fired_spe);

        % Spike-timing dependent plasticity computations:
        STDP_mot(fired_out,t+D)=0.1; % Keep a record of when the output neurons spiked.
   %    STDP_spe(fired_spe,t+D)=0.1; % Keep a record of when the spectrum neurons spiked.
        if strcmp(STDP,'STDP')
        STDP_reservor(fired,t+D)=0.1; %inter reservor STDPmatrix record


         for k=1:length(fired)
                sd_reservor(pre{fired(k)})=sd_reservor(pre{fired(k)})+STDP_reservor(N*t+aux{fired(k)}); %LTP A plus 1
                % pre{fired(1)=54}=52,78,210,277,350,372,734
         end;


        end


        for k=1:length(fired_mot)
                sd_mot(:,fired_mot(k))=sd_mot(:,fired_mot(k))+STDP_mot(:,t); % Adjusting sd for synapses eligible for potentiation.LTP
        end
   %    for k=1:length(fired_inp)
   %            sd_spe(:,fired_inp(k))=sd_spe(:,fired_inp(k))+STDP_spe(:,t); % Adjusting sd for synapses eligible for potentiation.LTP
   %    end


        firings=[firings;t*ones(length(fired),1),fired];                % Update the record of when neuronal firings occurred.
        outFirings=[outFirings;t*ones(length(fired_out),1),fired_out];  %
        motFirings=[motFirings;t*ones(length(fired_mot),1),fired_mot];  %
        inpFirings=[inpFirings;t*ones(length(fired_inp),1),fired_inp];  %
   %    speFirings=[speFirings;t*ones(length(fired_spe),1),fired_spe];  %


        % For any presynaptic neuron that just fired, calculate the current to add
        % as proportional to the synaptic strengths from its postsynaptic neurons.
        k=size(firings,1);
        while firings(k,1)>t-D
        del=delays{firings(k,2),t-firings(k,1)+1}; %1:M
        ind = post(firings(k,2),del);          % post neurons of fired neron id [23,45,13,14,53,??????]
            I(ind)=I(ind)+s(firings(k,2), del)';%'
         if strcmp(STDP,'STDP')

          sd_reservor(firings(k,2),del)=sd_reservor(firings(k,2),del)-1.5*STDP_reservor(ind,t+D)';%'LTD A-1.5
         end

            k=k-1;
        end;

        % Calculating currents to add for motor neurons.
        k=size(outFirings,1);
        while outFirings(k,1)>t-D
            del_mot=delays_mot{outFirings(k,2),t-outFirings(k,1)+1};
            ind_mot = post_mot(outFirings(k,2),del_mot);   %del_mot=1:Nmot ind_mot=
            I_mot(ind_mot)=I_mot(ind_mot)+2*sout(outFirings(k,2), del_mot)';%'
            k=k-1;
        end;

        % Calculating currents to add for input neurons.
   %    k=size(speFirings,1);
   %    while speFirings(k,1)>t-D
   %        del_inp=delays_inp{speFirings(k,2),t-speFirings(k,1)+1};
   %        ind_inp = post_spe(speFirings(k,2),del_inp);%del_inp=1:Ninp ind_inp=
   %        I(ind_inp)=I(ind_inp)+2*sinp(speFirings(k,2), del_inp)';%'
   %        k=k-1;
   
  %     end;


%%%%%%%%%%%%%%%
        %caliculate LST
        LAP=mean(mean(v));
        LAP_hist(t)=LAP;
%%%%%%%%%%%%%%%%%

        % Individual neuronal dynamics computations:
        v=v+0.5*((0.04*v+5).*v+140-u+I);                            % for numerical I...randam+from spe neurons(1*100)
        v=v+0.5*((0.04*v+5).*v+140-u+I);                            % stability time
        v_mot=v_mot+0.5*((0.04*v_mot+5).*v_mot+140-u_mot+I_mot);    % step is 0.5 ms I_mot...randam+from out neurons
        v_mot=v_mot+0.5*((0.04*v_mot+5).*v_mot+140-u_mot+I_mot);
  %     v_spe=v_spe+0.5*((0.04*v_spe+5).*v_spe+140-u_spe+I_spe);    % step is 0.5 ms I_spe...randam+from feedback
  %     v_spe=v_spe+0.5*((0.04*v_spe+5).*v_spe+140-u_spe+I_spe);
        u=u+a.*(0.2*v-u);
        u_mot=u_mot+a_mot.*(0.2*v_mot-u_mot);
  %     u_spe=u_spe+a_spe.*(0.2*v_spe-u_spe);



                % Exponential decay of the traces of presynaptic neuron firing
        STDP_mot(:,t+D+1)=0.95*STDP_mot(:,t+D);                             % tau = 20 ms
  %     STDP_spe(:,t+D+1)=0.95*STDP_spe(:,t+D);                             % tau = 20 ms

        if strcmp(STDP,'STDP')

          STDP_reservor(:,t+D+1)=0.95*STDP_reservor(:,t+D);

        end

                % Exponential decay of the dopamine concentration over time.
        DA=DA*0.995;

        % Modify synaptic weights.
        if (mod(t,10)==0)
            if strcmp(reinforce,'reinforce')
                sout=max(0,min(sm,sout+DA*sd_mot));
                sout=sout./(mean(mean(sout)));              % Normalizing the synaptic weights. out-mot
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
        if strcmp(IP,'Tonic')
            
           TEI(1:Ne,1) = TEI(1:Ne,1) + etaIP*(fired(1:Ne,1) - HIP);  
           TE.m = TE.m + etaIP*(fired_mot - HIP);
         
        end
        




        % Every testint seconds, use the motor neuron spikes to generate a sound.
        if (mod(sec,testint)==0)

            firedmusc1pos=find(v_mot(1:Nmot/2)>=30); % Find out which of the jaw/lip motor neurons fired.
            firedmusc1neg=find(v_mot(Nmot/2+1:end)>=30);
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



            end



            if t==1000 % Based on the 1 s timeseries of smoothed summed motor neuron spikes, generate a sound.


                if ~strcmp(reinforce,'range')
                    % Write the Praat script:

                    fid = fopen([wavdir,'/ressynth_',id,num2str(sec,'%d'),'.praat'],'w');
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

                        fprintf(fid,['\tWrite to WAV file... synth_',id,'_',num2str(sec),'.wav\n']);


                    fclose(fid);

                   % Execute the Praat script -- produces a wave file:
                    % Housekeeping
                    if ismac
                            system([praatPathmac, ' ', wavdir,'/ressynth_',id,num2str(sec,'%d'),'.praat']);
                            delete([wavdir,'/ressynth_',id,num2str(sec,'%d'),'.praat']);  %delete praat script
                    elseif ispc
                            system([praatPathpc, ' --run ', wavdir,'\ressynth_',id,num2str(sec,'%d'),'.praat']);
                            delete([wavdir,'\ressynth_',id,num2str(sec,'%d'),'.praat']);  %delete praat script

                    end



                    % Find the auditory salience of the sound:
                        salienceResults = auditorySalience([wavdir,'/synth_',id,'_',num2str(sec),'.wav'],0);

                    salience = sum(abs(salienceResults.saliency(31:180))); % Summing over salience trace to produce a single value.
                                                                           % abs(X)= return absolute value    neglect first 30  saliency=[1,180]

                    salhist(sec,1) = salience; % History of salience over entire simulation.
                    display(['salience =',num2str(salience)]);

                    if mod(sec,1)~=0
                      delete([wavdir,'/synth_',id,'_',num2str(sec),'.wav']);%
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
        if any(rew==sec*1000+t)   % what mean?
            DA=DA+DAinc;
        end
    end








    % Writing reservoir neuron firings for this second to a text file.     %recoring firings matrix
    if mod(sec,SAVINTV*testint)==0 || sec==T

            firings_fid = fopen([firingsdir,'/babble_daspnet_firings_',id,'_',num2str(sec),'.txt'],'w');

        for firingsrow = 1:size(firings,1)
            fprintf(firings_fid,'%i\t',sec);
            fprintf(firings_fid,'%i\t%i',firings(firingsrow,:));
            fprintf(firings_fid,'\n');
        end
        fclose(firings_fid);
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
        subplot(4,1,4);
        plot(smoothmusc(muscsmooth:1000)); ylim([-.5,.5]); xlim([-100,900]); % Plot the smoothed sum of motor neuron spikes 1 s timeseries
        title('Sum of Agonist/Antagonist Motor Neuron Activation', 'fontweight','bold');

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
    % ---- end plot ------

    sout_hist{sec}=sout;

    % Preparing STDP and firings for the following 1000 ms.
    STDP_mot(:,1:D+1)=STDP_mot(:,1001:1001+D);
%    STDP_spe(:,1:D+1)=STDP_spe(:,1001:1001+D);
    ind = find(firings(:,1) > 1001-D);
    firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
    ind_out = find(outFirings(:,1) > 1001-D);
    outFirings=[-D 0;outFirings(ind_out,1)-1000,outFirings(ind_out,2)];
    ind_mot = find(motFirings(:,1) > 1001-D);
    motFirings=[-D 0;motFirings(ind_mot,1)-1000,motFirings(ind_mot,2)];%
    ind_inp = find(inpFirings(:,1) > 1001-D);
    inpFirings=[-D 0;inpFirings(ind_inp,1)-1000,inpFirings(ind_inp,2)];
    feedbackhist=zeros(100,1);
 %   ind_spe = find(speFirings(:,1) > 1001-D);
 %   speFirings=[-D 0;speFirings(ind_spe,1)-1000,speFirings(ind_spe,2)];%

    if strcmp(STDP,'STDP')
       STDP_reservor(:,1:D+1)=STDP_reservor(:,1001:1001+D);%
    end



    % Every so often, save the workspace in case the simulation is interupted all data is not lost.
    if mod(sec,SAVINTV*testint)==0 || sec==T || sec==1
        display('Data Saving..... Do not exit program.');
        save(workspaceFilename, '-regexp', '^(?!(v_mot_hist)$).');

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
        csvwrite([workspacedir,'/sout_',sec,'.csv'],sout);
        csvwrite([workspacedir,'/s',sec,'.csv'],s);
 %       csvwrite([workspacedir,'/sinp',sec,'.csv'],sinp);
    end


    toc;

end

%create datefiles

if sec==T
  number=[1:T].'; %'
  salhistdate=[salhist,number];
  csvwrite([workspacedir,'/',id,'.csv'],salhistdate);
end