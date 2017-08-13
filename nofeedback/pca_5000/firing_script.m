%initialization
global d k yoke STDP outdir p dim 



d=5;    %Iterate number
yoke=['NY'];   %Yoked or NY
STDP=['STDP'];
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

ma=zeros(3,6);


%Babbling&mean caliculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  dim=[];
  mean_ruiseki=[];
  mean_kiyo=[];
   for k=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           [y,ruiseki,kiyo]=dimention(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           dim=[dim;y];
           if k==1
           mean_ruiseki=ruiseki;
           mean_kiyo=kiyo;
           else
               mean_ruiseki=mean_ruiseki+ruiseki;
               mean_kiyo=mean_kiyo+kiyo;
           end

   end
%p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',id,'_.txt';
%p=0.05_NY_STDPbabble_daspnet_firings_6_170222_2000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_

m=mean(dim)
s=std(dim)
rui=mean_ruiseki/d;
ki=mean_kiyo/d;
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);
ylim([0 1]);
whitebg('white')
fig1=plot(rui,'b');
hold on
%fig1=bar(ki);

saveas(fig1,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.png']);

%clf('reset');
 clearvars y dim m s rui ki ruiseki kiyo;
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



yoke=['Yoked'];   %Yoked or NY
STDP=['STDP'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  dim=[];
  mean_ruiseki=[];
  mean_kiyo=[];
   for k=1:d
        %conduct babbling
           display(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           [y,ruiseki,kiyo]=dimention(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           dim=[dim;y];
           if k==1
           mean_ruiseki=ruiseki;
           mean_kiyo=kiyo;
           else
               mean_ruiseki=mean_ruiseki+ruiseki;
               mean_kiyo=mean_kiyo+kiyo;
           end

   end
%p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',id,'_.txt';
%p=0.05_NY_STDPbabble_daspnet_firings_6_170222_2000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_

m=mean(dim)
s=std(dim)
rui=mean_ruiseki/d;
ki=mean_kiyo/d;
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);
ylim([0 1]);
fig2=plot(rui,'r');
hold on
%fig2=bar(ki);
saveas(fig2,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.png']);
%clf('reset');
 clearvars y dim m s rui ki ruiseki kiyo;
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



yoke=['NY'];   %Yoked or NY
STDP=['NSTDP'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  dim=[];
  mean_ruiseki=[];
  mean_kiyo=[];
   for k=1:d
        %conduct babbling
           display(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           [y,ruiseki,kiyo]=dimention(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           dim=[dim;y];
           if k==1
           mean_ruiseki=ruiseki;
           mean_kiyo=kiyo;
           else
               mean_ruiseki=mean_ruiseki+ruiseki;
               mean_kiyo=mean_kiyo+kiyo;
           end

   end
%p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',id,'_.txt';
%p=0.05_NY_STDPbabble_daspnet_firings_6_170222_2000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_

m=mean(dim)
s=std(dim)
rui=mean_ruiseki/d;
ki=mean_kiyo/d;
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);
ylim([0 1]);
fig3=plot(rui,'g');
hold on
%fig3=bar(ki);
saveas(fig3,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.png']);
%clf('reset');
 clearvars y dim m s rui ki ruiseki kiyo;
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


yoke=['Yoked'];   %Yoked or NY
STDP=['NSTDP'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  dim=[];
  mean_ruiseki=[];
  mean_kiyo=[];
   for k=1:d
        %conduct babbling
           display(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           [y,ruiseki,kiyo]=dimention(['../firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(k),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           dim=[dim;y];
           if k==1
           mean_ruiseki=ruiseki;
           mean_kiyo=kiyo;
           else
               mean_ruiseki=mean_ruiseki+ruiseki;
               mean_kiyo=mean_kiyo+kiyo;
           end

   end
%p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',id,'_.txt';
%p=0.05_NY_STDPbabble_daspnet_firings_6_170222_2000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_

m=mean(dim)
s=std(dim)
rui=mean_ruiseki/d;
ki=mean_kiyo/d;
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);
ylim([0 1]);
fig4=plot(rui,'k');
hold on
%fig4=bar(ki);
legend('No-STDP','Sc-STDP','Np-NSTDP','Sc-NSTDP')
saveas(fig4,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.png']);
%clf('reset');
 clearvars y dim m s rui ki ruiseki kiyo;
 close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!mkdir output/png
!find ./output/ -maxdepth 4 -name '*.png' | xargs -J % cp % ./output/png/


%{
yoke=['Yoked'];   %Yoked or NY
STDP=['STDP'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  dim=[];
   for i=1:d
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           y=dimention(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           y
           dim=[dim;y];
   end


m=mean(dim)
s=std(dim)
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);

 clearvars y dim m s;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yoke=['NY'];   %Yoked or NY
STDP=['NSTDP'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  dim=[];
   for i=1:5
        %conduct babbling
            display(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           y=dimention(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           y
           dim=[dim;y];
   end


m=mean(dim)
s=std(dim)
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);

 clearvars y dim m s;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yoke=['Yoked'];   %Yoked or NY
STDP=['NSTDP'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]

  mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP]);
  dim=[];
   for i=1:5
        %conduct babbling
           display(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           y=dimention(['firing_data/','p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'babble_daspnet_firings_',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoke,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',STDP,'_.txt']);

           y
           dim=[dim;y];
   end


m=mean(dim)
s=std(dim)
dim=[m;s;0;dim];
csvwrite([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/p=',num2str(p),'_',yoke,'_',STDP,'.csv'],dim);

 clearvars y dim m s;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
display('END');
%plot FIg with R
