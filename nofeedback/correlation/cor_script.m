%initialization
global d k yoke STDP outdir p dim simutime

!rm -r output

d=5;    %Iterate number
yoke=['NY'];   %Yoked or NY
STDP=['STDP'];
id=['170222'];
feedbacktime=1;
iterate=2000;     %iterate of created network
speinplate=0.5;
debug=1;
simutime=10000;   %iterate of simulation when create plot

%Create mean file
outdir=['output'];
if ~exist(outdir,'dir')
 mkdir(outdir);
end



%caliculate correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
   for k=1:d
        
        firingid = ['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt'];
        motorid  = ['../firing_data/p=',num2str(p),'_',yoke,'_',STDP,'/motorcommand_p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt'];
        display(['correlation_firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);
        correlation(motorid,firingid);   
   end
   
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



yoke=['Yoked'];   %Yoked or NY
STDP=['STDP'];


%caliculate correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
   for k=1:d
        
        firingid = ['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt'];
        motorid  = ['../firing_data/p=',num2str(p),'_',yoke,'_',STDP,'/motorcommand_p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt'];
        display(['correlation_firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);
        correlation(motorid,firingid);   
   end
   
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


yoke=['NY'];   %Yoked or NY
STDP=['NSTDP'];


%caliculate correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
   for k=1:d
        
        firingid = ['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt'];
        motorid  = ['../firing_data/p=',num2str(p),'_',yoke,'_',STDP,'/motorcommand_p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt'];
        display(['correlation_firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);
        correlation(motorid,firingid);   
   end
   
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yoke=['Yoked'];   %Yoked or NY
STDP=['NSTDP'];


%caliculate correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
   for k=1:d
        
        firingid = ['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt'];
        motorid  = ['../firing_data/p=',num2str(p),'_',yoke,'_',STDP,'/motorcommand_p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt'];
        display(['correlation_firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);
        correlation(motorid,firingid);   
   end
   
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!mkdir output/png
!find ./output/ -maxdepth 4 -name '*.png' | xargs -J % cp % ./output/png/
