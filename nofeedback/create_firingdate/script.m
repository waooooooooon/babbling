%initialization
!rm -r ../firing_data

d=5;    %Iterate number of simulation
yoked=['NY'];   %Yoked or NY
stdp=['STDP'];
id=['170222'];
feedbacktime=1;
iterate=2000;     %network iterate number (usually 2000)
speinplate=0.5;
debug=1;
simutime=10000;  %simutime of plotfiringslong

%Create mean file
%outdir=['output'];
%if ~exist(outdir,'dir')
% mkdir(outdir);
%end


%Babbling&mean caliculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  %mkdir([outdir,'/p=',num2str(p),'_',yoked,'_',stdp]);

   for i=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoked,'_',stdp,'/babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',num2str(simutime),'.txt']);

           y=plot_firingslong([num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp],iterate,'reinforce',1:100,4,yoked,1,feedbacktime,p,speinplate,stdp,debug);

           %csvwrite([outdir,'/p=',num2str(p),'_',yoked,'_',stdp,'/',num2str(i),'_p=',num2str(p),'_',yoked,'_',stdp,'.csv'],y);
   end


 clearvars y;
end


yoked=['NY'];   %Yoked or NY
stdp=['NSTDP'];

%Babbling&mean caliculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

%  mkdir([outdir,'/p=',num2str(p),'_',yoked,'_',stdp]);

   for i=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoked,'_',stdp,'/babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',num2str(simutime),'.txt']);

           y=plot_firingslong([num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp],iterate,'reinforce',1:100,4,yoked,1,feedbacktime,p,speinplate,stdp,debug);

           %csvwrite([outdir,'/p=',num2str(p),'_',yoked,'_',stdp,'/',num2str(i),'_p=',num2str(p),'_',yoked,'_',stdp,'.csv'],y);
   end


 clearvars y;
end


yoked=['Yoked'];   %Yoked or NY
stdp=['STDP'];

%Babbling&mean caliculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

%  mkdir([outdir,'/p=',num2str(p),'_',yoked,'_',stdp]);

   for i=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoked,'_',stdp,'/babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',num2str(simutime),'.txt']);

           y=plot_firingslong([num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp],iterate,'reinforce',1:100,4,yoked,1,feedbacktime,p,speinplate,stdp,debug);

           %csvwrite([outdir,'/p=',num2str(p),'_',yoked,'_',stdp,'/',num2str(i),'_p=',num2str(p),'_',yoked,'_',stdp,'.csv'],y);
   end


 clearvars y;
end


yoked=['Yoked'];   %Yoked or NY
stdp=['NSTDP'];

%Babbling&mean caliculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

 % mkdir([outdir,'/p=',num2str(p),'_',yoked,'_',stdp]);

   for i=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoked,'_',stdp,'/babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',num2str(simutime),'.txt']);

           y=plot_firingslong([num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp],iterate,'reinforce',1:100,4,yoked,1,feedbacktime,p,speinplate,stdp,debug);

           %csvwrite([outdir,'/p=',num2str(p),'_',yoked,'_',stdp,'/',num2str(i),'_p=',num2str(p),'_',yoked,'_',stdp,'.csv'],y);
   end


 clearvars y;
end


