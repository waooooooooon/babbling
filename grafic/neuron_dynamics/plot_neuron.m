rng shuffle;

% Initialization.
salthresh = 4.5;            % Initial salience value for reward (used in 'relhisal' reinforcment).
DAinc = 1;                  % Amount of dopamine given during reward.
sm = 4;                     % Maximum synaptic weight.
testint = 1;                % Number of seconds between vocalizations.


% Directory for Coath et. al. Saliency Detector.
addpath('auditorysaliencymodel');



    M=100;                 % number of synapses per neuron
    D=1;                   % maximal conduction delay
    % excitatory neurons   % inhibitory neurons      % total number
    Ne=800;                Ni=200;                   N=Ne+Ni;
    
    Nmot = 100;
    Nout = 100;
    a=[0.02*ones(Ne,1);    0.1*ones(Ni,1)];     % Sets time scales of membrane recovery variable.
    d=[   8*ones(Ne,1);    2*ones(Ni,1)];       % Membrane recovery variable after-spike shift.
    a_mot=.02*ones(Nmot,1);
    d_mot=8*ones(Nmot,1);
    post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); % Assign the postsynaptic neurons for each neuron's synapse in the reservoir.
    post_mot=repmat(1:Nmot,Nout,1); % All output neurons connect to all motor neurons.

    s=[rand(Ne,M);-rand(Ni,M)];         % Synaptic weights in the reservoir.
    sout=rand(Nout,Nmot); % Synaptic weights from the reservoir output neurons to the motor neurons.

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


sec =1;
time = 300;
v_hist= zeros(1000,time);
I = 0;
I_mot = 0;


%RUNNING THE SIMULATION%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for t=1:time                          % Millisecond timesteps

        %Random Thalamic Input.
        if t>50
        I=20*(ones(N,1));
        I_mot=13*(rand(Nmot,1)-0.5);
        end
        I_hist(t) =I(1,1);

        fired = find(v>=30);                % Indices of fired neurons
        fired_mot = find(v_mot>=30);
        v(fired)=-65;                       % Reset the voltages for those neurons that fired
        v_mot(fired_mot)=-65;
        u(fired)=u(fired)+d(fired);         % Individual neuronal dynamics
        u_mot(fired_mot)=u_mot(fired_mot)+d_mot(fired_mot);

        % Spike-timing dependent plasticity computations:

        firings=[firings;t*ones(length(fired),1),fired];                % Update the record of when neuronal firings occurred.
        motFirings=[motFirings;t*ones(length(fired_mot),1),fired_mot];
        % For any presynaptic neuron that just fired, calculate the current to add
        % as proportional to the synaptic strengths from its postsynaptic neurons.
        k=size(firings,1);


        % Individual neuronal dynamics computations:
        v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
        v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
        %v=v+0.5*((0.04*v+5).*v+140-u+I);                            % for numerical
        %v=v+0.5*((0.04*v+5).*v+140-u+I);                            % stability time
        v_mot=v_mot+0.5*((0.04*v_mot+5).*v_mot+140-u_mot+I_mot);    % step is 0.5 ms
        v_mot=v_mot+0.5*((0.04*v_mot+5).*v_mot+140-u_mot+I_mot);
        u=u+a.*(0.2*v-u);
        u_mot=u_mot+a_mot.*(0.2*v_mot-u_mot);
        
        
        v_hist(:,t) = v(:,1);

    end
    

    v_hist(find(v_hist>30))=30;
    % Writing reservoir neuron firings for this second to a text file.


    % Make axis labels and titles for plots that are being kept.
    % ---- plot -------

        hNeural = figure(103);
        hNeural = plot(v_hist(500,:),'k-','LineWidth',2); % Plot all the neurons' spike
        
        axis([0 time -90 40]);
        set(gca,'FontSize',20);
        saveas(hNeural,['./RS.png']);
        
        
        
        fig111 = plot(v_hist(850,:),'k','LineWidth',2); % Plot all the neurons' spike
        
        axis([0 time -90 40]);
        set(gca,'FontSize',20);
        saveas(fig111,['./FS.png']);
        
        fig222 = plot(I_hist,'r','LineWidth',2); % Plot all the neurons' spike
        
        axis([0 time 0 40]);
        set(gca,'FontSize',20);
        saveas(fig222,['./Input.png']);     
        
