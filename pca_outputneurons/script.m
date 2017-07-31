%initialization
d=10;    %Iterate number
yoked=['NY'];   %Yoked or NY
stdp=['STDP'];
id=['170222'];
feedbacktime=1;
iterate=2000;
speinplate=0.5;
debug=1;

%Create mean file
outdir=['output'];
if ~exist(outdir,'dir')
 mkdir(outdir);
end


%Babbling&mean caliculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.1,0.075,0.025,0.01]

  mkdir([outdir,'/p=',num2str(p),'_',yoked,'_',stdp]);

   for i=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoked,'_',stdp,'/babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_2000.txt']);

           y=plot_firings([num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp],iterate,'reinforce',1:100,4,yoked,1,feedbacktime,p,speinplate,stdp,debug);

           %csvwrite([outdir,'/p=',num2str(p),'_',yoked,'_',stdp,'/',num2str(i),'_p=',num2str(p),'_',yoked,'_',stdp,'.csv'],y);
   end


 clearvars y;
end