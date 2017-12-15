function [] = babbling(id,newT,reinforce,outInd,muscscale,yoke,plotOn,feedbacktime,learningratio,speinplate,STDP,debug,IP,separatephase,Network,reward)
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


global simutime tag



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
SAVINTV=1000;
LST_hist=[1:1000];
muscle_number=0;
negativereward = 0;
nega_rate=1/10;
mot_thre = 0;


%
sparse_mot = 1;   %1 or 0 
sparse_degree = 70; % number of synapse from out to motor Nout - sparse_degree


%for low pass filter
cutoff_frequency = 10; % Design a 70th order lowpass FIR filter with cutoff frequency of 75 Hz.
Fs = 1000;                    % sample rate in Hz
Fnorm = cutoff_frequency/(Fs/2);           % Normalized frequency
df = designfilt('lowpassfir','FilterOrder',70,'CutoffFrequency',Fnorm);


if strcmp(IP,'IP')
 TE.max = -15;             %for IP
 TI.max = -15.2;             %for IP
 etaIP = 0.005;          %for IP 0.01
 threshold = -55.1;          %for IP  -55.1
 sumweight =5000;          %for SN(amount of wight of output-motor)
 %HIP = 1/200  ;  %target firing rate (100 = numer of input neuron) defalt 2*input/Ne
 %%%%%%%%%%%%%%%%%%%%
 %normal distribution
 %variance_IP = ones(1000,1)/200;
 %m_IP = ones(100,1)/200;
 %V_HIP = normrnd(variance_IP,0.001);
 %m_HIP = normrnd(m_IP,0.001);
 
 %lognormal distribution
 m = 1/100;
 v = 0.0001;
 mu = log((m^2)/sqrt(v+m^2));
 sigma = sqrt(log(v/(m^2)+1));
 V_HIP = lognrnd(mu,sigma,1000,1);
 %m_HIP = lognrnd(mu,sigma,100,1);      %lonnomal of motor neuron
 m_HIP = normrnd(ones(100,1)/200,0.0000001);        %norml distribution
 %for debag
 %debag_HIP = V_HIP*Fs;
 %edge = logspace(0, 10, 300);
 %h=histogram(debag_HIP,edge);set(gca,'Xscale','log');xlim([0 100]);
 %%%%%%%%%%%%%%%%%%%%%%

end




% Directory names for data.
wavdir = [id, '_Wave'];
firingsdir = [id, '_Firings'];
workspacedir = [id, '_Workspace'];
setdir = ['setting'];
createddata_dir = ['~/babbling/created_data/'];     %data dir
id_dir = [tag,'/'];
outdir = [createddata_dir,id_dir,'network_analysis'];
firingdir = [outdir,'/firing_data'];


% Error Checking.
if(any(isspace(id)))
    disp('Please choose an id without spaces.');
    return
end

% Creating data directories.

if ~exist(outdir, 'dir')
    mkdir(outdir);
else
    addpath(outdir);
end

if ~exist(firingdir, 'dir')
    mkdir(firingdir);
else
    addpath(firingdir);
end


% Creating workspace names.
workspaceFilename=[workspacedir,'/babble_daspnet_reservoir_',id,'.mat'];



% Directory for Coath et. al. Saliency Detector.
addpath('auditorysaliencymodel');


%Import initial value
table=importdata([setdir,'/table_notnormalization_from0.5to1.5KHz.csv']);

%Import matfile
load([createddata_dir,id_dir,'data/',workspaceFilename]);

%initialization
sec = 0;


T=newT;
clearvars newT;

% Absolute path where Praat can be found.
praatPathpc = 'c:\users\Takimto\Praat\Praat';
praatPathmac = '/Applications/Praat.app/Contents/MacOS/Praat';
praatPathlinux = '/usr/bin/praat';


datahistsize=((1000-muscsmooth)/feedbacktime);



%RUNNING THE SIMULATION%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sec=1:2 % T is the duration of the simulation in seconds.
          tic;
    display('********************************************');
    display(['yoke=',yoke,'_phase=',STDP,'_',IP,'_Second ',num2str(sec),' of ',num2str(2)]);

    v_mot_hist{sec}=[]; % Record of all the membrane voltages of the motor neurons.


    %Initialize decaysmooth
    decaysmoothneg=zeros(1,simutime);
    decaysmoothpos=zeros(1,simutime);
    datatime=1;

    for t=1:simutime                         % Millisecond timesteps

        %Random Thalamic Input or No random input (IP)
        if strcmp(IP,'Tonic')
            I=13*(rand(N,1)-0.5);
            I_mot=13*(rand(Nmot,1)-0.5);
        elseif strcmp(IP,'IP')
            I=zeros(N,1);
            I_mot=zeros(Nmot,1);
        elseif strcmp(IP,'afterIP')
            I=zeros(N,1);
            I_mot=zeros(Nmot,1);
        end
        
        %for ip test
        %I(1:Ninp)=13*(rand(Nmot,1)-0.5);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%   feedback every time
        if strcmp(Network,'random')
           if t-1>muscsmooth&mod(t-1,feedbacktime)==0

            feedback1=speinplate*table(:,muscle_number);
            if strcmp(yoke,'No')
                I(1:Ninp)=I(1:Ninp)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,feedback1];
            elseif strcmp(yoke,'Sc')
                randid=randperm(100);
                yokedfeedback=feedback1(randid);
                I(1:Ninp)=I(1:Ninp)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,yokedfeedback];
            end

           end
        elseif strcmp(Network,'lattice')
            
           if t-1>muscsmooth&mod(t-1,feedbacktime)==0

            feedback1=speinplate*table(:,muscle_number);
            if strcmp(yoke,'No')
                I(InputneuronID)=I(InputneuronID)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,feedback1];
            elseif strcmp(yoke,'Sc')
                randid=randperm(100);
                yokedfeedback=feedback1(randid);
                I(InputneuronID)=I(InputneuronID)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,yokedfeedback];
            end

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
            
            %fired_mot = find((v_mot-TE.m)>=threshold);     %IP for motor
            fired_mot = find((v_mot)>=mot_thre);      %not IP for motor neuron
            
            
            %%%%%%%%% SN
            s_i = s(801:1000,:);
            s_e = sesum*s(1:Ne,:)/sum(sum(s(1:Ne,:)));
            %s_e = s(1:Ne,:)/sum(sum(s(1:Ne,:)));%ŽÀŒ±
            s = [s_e ; s_i];
            %sout = smsum*sout/sum(sum(sout));
            %sout = sumweight*sout/sum(sum(sout));%ŽÀŒ±
            
            
        elseif strcmp(IP,'afterIP')
            
            fired = find((v-TEI)>=threshold);
            fired_out = find((v(Ninp + outInd)-TEI(Ninp + outInd))>=threshold);
            fired_inp = find((v(inpInd) - TEI(inpInd))>=threshold);
            
            fired_mot = find((v_mot-TE.m)>=threshold);
            
            %%%%%%%%% SN
            s_i = s(801:1000,:);
            s_e = sesum*s(1:Ne,:)/sum(sum(s(1:Ne,:)));
            %s_e = s(1:Ne,:)/sum(sum(s(1:Ne,:)));%ŽÀŒ±
            s = [s_e ; s_i];
            %sout = smsum*sout/sum(sum(sout));
            %sout = sumweight*sout/sum(sum(sout));%ŽÀŒ±
        end
        
        


        v(fired)=-65;                       % Reset the voltages for those neurons that fired
        v_mot(fired_mot)=-65;
   %    v_spe(fired_spe)=-65;
        u(fired)=u(fired)+d(fired);         % Individual neuronal dynamics
        u_mot(fired_mot)=u_mot(fired_mot)+d_mot(fired_mot);
   %    u_spe(fired_spe)=u_spe(fired_spe)+d_spe(fired_spe);

        % Spike-timing dependent plasticity computations:
        STDP_out(fired_out,t+D)=0.1; % Keep a record of when the output neurons spiked.
        STDP_mot(fired_mot,t+D)=0.1; % Keep a record of when the spectrum neurons spiked.
        if strcmp(STDP,'STDP')
        STDP_reservor(fired,t+D)=0.1; %inter reservor STDPmatrix record


         for k=1:length(fired)
                sd_reservor(pre{fired(k)})=sd_reservor(pre{fired(k)})+STDP_reservor(N*t+aux{fired(k)}); %LTP A plus 1
                % pre{fired(1)=54}=52,78,210,277,350,372,73
               
                
         end;


        end


        for k=1:length(fired_mot)
                sd_mot(:,fired_mot(k))=sd_mot(:,fired_mot(k))+STDP_out(:,t); % Adjusting sd for synapses eligible for potentiation.LTP
        end
   

        firings=[firings;t*ones(length(fired),1),fired];                % Update the record of when neuronal firings occurred.
        outFirings=[outFirings;t*ones(length(fired_out),1),fired_out];  %
        motFirings=[motFirings;t*ones(length(fired_mot),1),fired_mot];  %
        inpFirings=[inpFirings;t*ones(length(fired_inp),1),fired_inp];  %
   %    speFirings=[speFirings;t*ones(length(fired_spe),1),fired_spe];  %


        % For any presynaptic neuron that just fired, calculate the current to add
        % as proportional to the synaptic strengths from its postsynaptic neurons.
        if strcmp(Network,'random')
            k=size(firings,1);
            while firings(k,1)>t-D
            del=delays{firings(k,2),t-firings(k,1)+1}; %1:M synaptic data
            ind = post(firings(k,2),del);          % post neurons of fired neron id [23,45,13,14,53,??????]
                I(ind)=I(ind)+s(firings(k,2), del)';%'
             if strcmp(STDP,'STDP')

              sd_reservor(firings(k,2),del)=sd_reservor(firings(k,2),del)-1.5*STDP_reservor(ind,t+D)';%'LTD A-1.5
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

              sd_reservor(firings(k,2),del)=sd_reservor(firings(k,2),del)-1.5*STDP_reservor(ind,t+D)';
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
            
            sd_mot(outFirings(k,2),:)=sd_mot(outFirings(k,2),:)-1.5*STDP_mot(:,t+D)'; % LTD of motor neuron 
            
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
  
        
        
        
        
        
        
        % Every testint seconds, use the motor neuron spikes to generate a sound.


            if strcmp(IP,'Tonic')
                firedmusc1pos=find(v_mot(1:Nmot/2)>=30); % Find out which of the jaw/lip motor neurons fired.
                firedmusc1neg=find(v_mot(Nmot/2+1:end)>=30);
                
            elseif strcmp(IP,'IP')
                %firedmusc1pos=find(v_mot(1:Nmot/2)-TE.m(1:Nmot/2)>=threshold); % Find out which of the jaw/lip motor neurons fired.
                %firedmusc1neg=find(v_mot(Nmot/2+1:end)-TE.m(Nmot/2+1:end)>=threshold);
                
                firedmusc1pos=find(v_mot(1:Nmot/2)>=mot_thre); % Find out which of the jaw/lip motor neurons fired.
                firedmusc1neg=find(v_mot(Nmot/2+1:end)>=mot_thre);
                
            elseif strcmp(IP,'afterIP')
                firedmusc1pos=find(v_mot(1:Nmot/2)-TE.m(1:Nmot/2)>=threshold); % Find out which of the jaw/lip motor neurons fired.
                firedmusc1neg=find(v_mot(Nmot/2+1:end)-TE.m(Nmot/2+1:end)>=threshold);
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
              
              

            end
    end         %end fo firings
        
            
            

            
    
    
    
    % Make axis labels and titles for plots that are being kept.
    % ---- plot -------
    if 0
        hNeural = figure(103);
        set(hNeural, 'name', ['Neural Spiking for Second: ', num2str(sec)], 'numbertitle','off');
        subplot(4,1,1)
        plot(firings(:,1),firings(:,2),'.'); % Plot all the neurons'' spikes
        title('Reservoir Firings', 'fontweight','bold');
        axis([0 simutime 0 N]);
        subplot(4,1,2)
        plot(outFirings(:,1),outFirings(:,2),'.'); % Plot the output neurons'' spikes
        title('Output Neuron Firings', 'fontweight','bold');
        axis([0 simutime 0 Nout]);
     %   subplot(5,1,3)
     %   plot(speFirings(:,1),speFirings(:,2),'.'); % Plot the output neurons'' spikes
     %   title('Input Neuron Firings', 'fontweight','bold');
     %   axis([0 1000 0 Ninp]);
        subplot(4,1,3)
        plot(motFirings(:,1),motFirings(:,2),'.'); % Plot the motor neurons'' spikes
        title('Motor Neuron Firings', 'fontweight','bold');
        axis([0 simutime 0 Nmot]);
        if strcmp(reward,'negativereward')
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
        if strcmp(reward,'negativereward')
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
    

     
    if sec==2
        % Writing reservoir neuron firings for this second to a text file.     %recoring firings matrix

            firings_fid = fopen([firingdir,'/firing_',id,'_',num2str(simutime),'.txt'],'w');

        for firingsrow = 1:size(firings,1)
            fprintf(firings_fid,'%i\t',sec);
            fprintf(firings_fid,'%i\t%i',firings(firingsrow,:));
            fprintf(firings_fid,'\n');
        end
        fclose(firings_fid);
        
        
        motorcommand = fopen([firingdir,'/motorcommand_',id,'_',num2str(simutime),'.txt'],'w');

        size_musc = size(smoothmusc);
        for i = 1:size_musc(1,2)
            fprintf(firings_fid,'%i\t',i);
            fprintf(firings_fid,'%i\t%i',smoothmusc(i));
            fprintf(firings_fid,'\n');
        end
        fclose(motorcommand);       
        

    end
    
    
    %for IP, caliculate average firing rate
    f_rate = size(find(firings(:,2)==500));
    firing_rate_500 = f_rate(1,1);
    %f_rate = size(firings);
    %firing_rate = (f_rate(1,1)-1)/1000;
    display(['firing rate of 500= ',num2str(firing_rate_500)]);
    display(['neuron500 threshold = ',num2str(TEI(500,1))]);

    sout_hist{sec}=sout;
    
    %for negative reward
    display(['nega= ',num2str(nega(sec))]);
    display(['DA= ',num2str(DA)]);


    
    % Preparing STDP and firings for the following 1000 ms.
    STDP_out(:,1:D+1)=STDP_out(:,simutime+1:simutime+1+D);
    STDP_mot(:,1:D+1)=STDP_mot(:,simutime+1:simutime+1+D);
    ind = find(firings(:,1) > simutime+1-D);
    firings=[-D 0;firings(ind,1)-simutime,firings(ind,2)];
    ind_out = find(outFirings(:,1) > simutime+1-D);
    outFirings=[-D 0;outFirings(ind_out,1)-simutime,outFirings(ind_out,2)];
    ind_mot = find(motFirings(:,1) > simutime+1-D);
    motFirings=[-D 0;motFirings(ind_mot,1)-simutime,motFirings(ind_mot,2)];%
    ind_inp = find(inpFirings(:,1) > simutime+1-D);
    inpFirings=[-D 0;inpFirings(ind_inp,1)-simutime,inpFirings(ind_inp,2)];
    feedbackhist=zeros(100,1);
 %   ind_spe = find(speFirings(:,1) > 1001-D);
 %   speFirings=[-D 0;speFirings(ind_spe,1)-1000,speFirings(ind_spe,2)];%

    if strcmp(STDP,'STDP')
       STDP_reservor(:,1:D+1)=STDP_reservor(:,simutime+1:simutime+1+D);%
    end



    end


    toc;
  
end


