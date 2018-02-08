function [] = babbling(ID,newT,reinforce,outInd,muscscale,yoke,plotOn,feedbacktime,learningratio,speinplate,STDP,debug,IP,separatephase,Network,reward,feedbacktype,NoID)
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


global simutime tag firingdir motortype

smoothmusc = zeros(1,simutime);
summusc1posspikes = zeros(1,simutime);
summusc1negspikes = zeros(1,simutime);


% Directory names for data.
wavdir = [ID, '_Wave'];
firingsdir = [ID, '_Firings'];
workspacedir = [ID, '_Workspace'];
setdir = ['~/babbling/simulation/setting'];
createddata_dir = ['~/babbling/created_data/'];     %data dir
id_dir = [tag,'/'];
%id_dir = ['180120_SctimeNone/'];
outdir = [createddata_dir,id_dir,'network_analysis'];
outfiringdir = [outdir,'/onfeedback_firing_data'];





% Error Checking.
if(any(isspace(ID)))
    disp('Please choose an id without spaces.');
    return
end

% Creating data directories.

if ~exist(outdir, 'dir')
    mkdir(outdir);
else
    addpath(outdir);
end

if ~exist(outfiringdir, 'dir')
    mkdir(outfiringdir);
else
    addpath(outfiringdir);
end


% Creating workspace names.
workspaceFilename=[workspacedir,'/babble_daspnet_reservoir_',ID,'.mat'];



% Directory for Coath et. al. Saliency Detector.
addpath('auditorysaliencymodel');



%Import initial value
table=importdata([setdir,'/table_notnormalization_from0.5to1.5KHz.csv']);


ID_backup = ID;
%Import matfile
load([createddata_dir,id_dir,'data/',workspaceFilename]);
ID = ID_backup;

%initialization
sec = 0;


T=newT;
clearvars newT;

% Absolute path where Praat can be found.
praatPathpc = 'c:\users\Takimto\Praat\Praat';
praatPathmac = '/Applications/Praat.app/Contents/MacOS/Praat';
praatPathlinux = '/usr/bin/praat';


datahistsize=((1000-muscsmooth)/feedbacktime);




%%%%%%%% create feedback
if strcmp(motortype,'sin')
muscle_sin = sin(linspace(pi,-pi,500))*0.9;
muscle_sin = repmat(muscle_sin,1,simutime/500);
end


%%%%%%% import No feedback
if strcmp(yoke,'Sc') && strcmp(feedbacktype,'fft') && strcmp(separatephase,'randSc') && ~strcmp(motortype,'sin')
motcommanddata=importdata([firingdir,'/motorcommand_onfeedback_',motortype,'_',NoID,'_',num2str(simutime),'.txt']);
end


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
                if strcmp(motortype,'sin')
                    feedback1=speinplate*table(:,round((muscle_sin(t)+1)*10000));
                elseif strcmp(motortype,'feedback')
                    feedback1=speinplate*table(:,muscle_number);
                end
                
                I(1:Ninp)=I(1:Ninp)+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                I((Ninp*2+1):(Ninp*2+Ninp))=I((Ninp*2+1):(Ninp*2+Ninp))+feedback1; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,feedback1];
            elseif strcmp(yoke,'Sc') && strcmp(feedbacktype,'fft') && ~strcmp(separatephase,'randSc')
                if strcmp(motortype,'sin')
                    feedback1=speinplate*table(:,muscle_sin(t));
                elseif strcmp(motortype,'feedback')
                    feedback1=speinplate*table(:,muscle_number);
                end
                randid=randperm(100);
                yokedfeedback=feedback1(randid);
                I(1:Ninp)=I(1:Ninp)+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                I((Ninp*2+1):(Ninp*2+Ninp))=I((Ninp*2+1):(Ninp*2+Ninp))+yokedfeedback; %refrect feedback to spectrum neurons 0~2000hz/20
                feedbackhist=[feedbackhist,yokedfeedback];
                
                
            elseif strcmp(yoke,'Sc') && strcmp(feedbacktype,'fft') && strcmp(separatephase,'randSc')    
                if strcmp(motortype,'sin')
                    feedback1=speinplate*table(:,round((muscle_sin(t)+1)*10000));
                elseif strcmp(motortype,'feedback')
                    feedback1=speinplate*table(:,round((motcommanddata(t,2)+1)*10000));   %‚Ç‚¤‚·‚é‚©H
                end
                
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
            %s_e = s(1:Ne,:)/sum(sum(s(1:Ne,:)));%ŽÀŒ±
            s = [s_e ; s_i];
            %sout = smsum*sout/sum(sum(sout));
            %sout = sumweight*sout/sum(sum(sout));%ŽÀŒ±
            
            
        elseif strcmp(IP,'afterIP')
            
            fired = find((v-TEI)>=threshold);
            fired_out = find((v(Ninp + outInd)-TEI(Ninp + outInd))>=threshold);
            fired_inp = find((v(inpInd) - TEI(inpInd))>=threshold);
            
            %fired_mot = find((v_mot-TE.m)>=threshold);
            fired_mot = find((v_mot)>=mot_thre);      %not IP for motor neuron
            
            %%%%%%%%% SN
            s_i = s(801:1000,:);
            s_e = sesum*s(1:Ne,:)/sum(sum(s(1:Ne,:)));
            %s_e = s(1:Ne,:)/sum(sum(s(1:Ne,:)));%ŽÀŒ±
            s = [s_e ; s_i];
            %sout = smsum*sout/sum(sum(sout));
            %sout = sumweight*sout/sum(sum(sout));%ŽÀŒ±
            
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
  
        
        
        
        
        % Every testint seconds, use the motor neuron spikes to generate a sound.


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
              
              %muscle_his(t,sec) = muscle_number;
              
              

            end
            
            
        end         %end fo firings

        
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
    

    

     
    if sec==2
        % Writing reservoir neuron firings for this second to a text file.     %recoring firings matrix

            firings_fid = fopen([outfiringdir,'/firing_onfeedback_',motortype,'_',ID,'_',num2str(simutime),'.txt'],'w');

        for firingsrow = 1:size(firings,1)
            fprintf(firings_fid,'%i\t',sec);
            fprintf(firings_fid,'%i\t%i',firings(firingsrow,:));
            fprintf(firings_fid,'\n');
        end
        fclose(firings_fid);
        
         
        firings_mot = fopen([outfiringdir,'/motfiring_onfeedback_',motortype,'_',ID,'_',num2str(simutime),'.txt'],'w');

        for firingsrow = 1:size(motFirings,1)
            fprintf(firings_mot,'%i\t',sec);
            fprintf(firings_mot,'%i\t%i',motFirings(firingsrow,:));
            fprintf(firings_mot,'\n');
        end
        fclose(firings_mot);
        
        
        motorcommand = fopen([outfiringdir,'/motorcommand_onfeedback_',motortype,'_',ID,'_',num2str(simutime),'.txt'],'w');

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

    sout_hist{sec}=sout;
    
    %for negative reward

    
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
    
    
    
    toc;
    end


    
  
end


