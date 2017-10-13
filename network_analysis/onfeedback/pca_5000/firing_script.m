%initialization
global d k yoke STDP outdir p dim simutime

!rm -r output

d=2;    %Iterate number
yoke=['NY'];   %Yoked or NY
STDP=['STDP'];
id=['170222'];
feedbacktime=1;
iterate=2000;     %iterate of created network
speinplate=0.5;
debug=1;
simutime=100000;   %iterate of simulation when create plot

%Create mean file
outdir=['output'];
if ~exist(outdir,'dir')
 mkdir(outdir);
end

ma=zeros(3,6);


%Babbling&mean caliculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  transfer1=[];
  dim=[];
  mean_ruiseki=[];
  mean_kiyo=[];
   for k=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);

           [y,ruiseki,kiyo,transfer_score]=dimention(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);

           dim=[dim;y];
           if k==1
           mean_ruiseki=ruiseki;
           mean_kiyo=kiyo;
           else
               mean_ruiseki=mean_ruiseki+ruiseki;
               mean_kiyo=mean_kiyo+kiyo;
           end
           
           transfer1=[transfer1;transfer_score];  %save transfer score
           
   end
%p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',id,'_.txt';
%p=0.05_NY_STDPbabble_daspnet_firings_6_170222_2000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_

% analysis of dim
m=mean(dim)
s=std(dim)
rui=mean_ruiseki/d;
ki=mean_kiyo/d;
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/dim_p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);

% analysis of transfer
m2=mean(transfer1)
s2=std(transfer1)
transfer1=[m2;s2;0;transfer1];
%csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/transfer_p=',num2str(p),'_',yoke,'_',STDP,'.csv'],transfer1);




%%whitebg('white')
%fig1=plot(rui,'b');
%hold on
%fig1=bar(ki);

%saveas(fig1,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.png']);

%clf('reset');
 clearvars y dim m s rui ki ruiseki kiyo;
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



yoke=['Yoked'];   %Yoked or NY
STDP=['STDP'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  transfer2=[];
  dim=[];
  mean_ruiseki=[];
  mean_kiyo=[];
   for k=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);

           [y,ruiseki,kiyo,transfer_score]=dimention(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);

           dim=[dim;y];
           if k==1
           mean_ruiseki=ruiseki;
           mean_kiyo=kiyo;
           else
               mean_ruiseki=mean_ruiseki+ruiseki;
               mean_kiyo=mean_kiyo+kiyo;
           end
           
           transfer2=[transfer2;transfer_score];  %save transfer score
           
   end
%p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',id,'_.txt';
%p=0.05_NY_STDPbabble_daspnet_firings_6_170222_2000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_

% analysis of dim
m=mean(dim)
s=std(dim)
rui=mean_ruiseki/d;
ki=mean_kiyo/d;
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/dim_p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);

% analysis of transfer
m2=mean(transfer2)
s2=std(transfer2)
transfer2=[m2;s2;0;transfer2];
%csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/transfer_p=',num2str(p),'_',yoke,'_',STDP,'.csv'],transfer2);




%%whitebg('white')
%fig1=plot(rui,'b');
%hold on
%fig1=bar(ki);

%saveas(fig1,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.png']);

%clf('reset');
 clearvars y dim m s rui ki ruiseki kiyo;
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


yoke=['NY'];   %Yoked or NY
STDP=['NSTDP'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  transfer3=[];
  dim=[];
  mean_ruiseki=[];
  mean_kiyo=[];
   for k=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);

           [y,ruiseki,kiyo,transfer_score]=dimention(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);

           dim=[dim;y];
           if k==1
           mean_ruiseki=ruiseki;
           mean_kiyo=kiyo;
           else
               mean_ruiseki=mean_ruiseki+ruiseki;
               mean_kiyo=mean_kiyo+kiyo;
           end
           
           transfer3=[transfer3;transfer_score];  %save transfer score
           
   end
%p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',id,'_.txt';
%p=0.05_NY_STDPbabble_daspnet_firings_6_170222_2000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_

% analysis of dim
m=mean(dim)
s=std(dim)
rui=mean_ruiseki/d;
ki=mean_kiyo/d;
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/dim_p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);

% analysis of transfer
m2=mean(transfer3)
s2=std(transfer3)
transfer3=[m2;s2;0;transfer3];
%csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/transfer_p=',num2str(p),'_',yoke,'_',STDP,'.csv'],transfer3);




%ylim([0 1]);
%whitebg('white')
%fig1=plot(rui,'b');
%hold on
%fig1=bar(ki);

%saveas(fig1,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.png']);

%clf('reset');
 clearvars y dim m s rui ki ruiseki kiyo;
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yoke=['Yoked'];   %Yoked or NY
STDP=['NSTDP'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  transfer4=[];
  dim=[];
  mean_ruiseki=[];
  mean_kiyo=[];
   for k=1:d

        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);

           [y,ruiseki,kiyo,transfer_score]=dimention(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_',num2str(simutime),'.txt']);

           dim=[dim;y];
           if k==1
           mean_ruiseki=ruiseki;
           mean_kiyo=kiyo;
           else
               mean_ruiseki=mean_ruiseki+ruiseki;
               mean_kiyo=mean_kiyo+kiyo;
           end
           
           transfer4=[transfer4;transfer_score];  %save transfer score
           
   end
%p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',id,'_.txt';
%p=0.05_NY_STDPbabble_daspnet_firings_6_170222_2000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_

% analysis of dim
m=mean(dim)
s=std(dim)
rui=mean_ruiseki/d;
ki=mean_kiyo/d;
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/dim_p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);

% analysis of transfer
m2=mean(transfer4)
s2=std(transfer4)
transfer4=[m2;s2;0;transfer4];
%csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/transfer_p=',num2str(p),'_',yoke,'_',STDP,'.csv'],transfer4);




%ylim([0 1]);
%whitebg('white')
%fig1=plot(rui,'b');
%hold on
%fig1=bar(ki);

%saveas(fig1,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.png']);

%clf('reset');
 clearvars y dim m s rui ki ruiseki kiyo;
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% write transfer.csv
transfer=[transfer1,transfer2,transfer3,transfer4];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/transfer_onfeedback_p=',num2str(p),'.csv'],transfer);
%%%%%%%%%%%%%%%%

!mkdir output/png
!mkdir output/csv
!find ./output/ -maxdepth 4 -name '*.png' | xargs -J % cp % ./output/png/
!find ./output/ -maxdepth 4 -name '*.csv' | xargs -J % cp % ./output/csv/



display('END');
%plot FIg with R
