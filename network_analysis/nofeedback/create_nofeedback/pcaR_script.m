%importdata
global tag simutime createddata_dir id_dir outdir firingdir pca_dir

d=1;    %Iterate number
YOKED=['No';'Sc'];   %Sc or No
ploton = 0; % 1 or 0
STDP=['STDP';'NSTD'];
tag=['180120_Sctime'];
IP =['LiIP'];        %threIP or Tonic or NoIP or LiIP
separatephase = ['randSc'];      %separate or nseparate or randSc
Network = ['random'];      %lattice or random
reward = ['normal'];        %nega or normal
feedbacktype = ['fft'];        %consonant or fft or none
feedbacktime=1;
iterate=2000;
speinplate=0.3;
debug=0;
p = 0.03;
created_data = ['../created_data/'];



%for firing
simutime=10000;        %simutime of plotfiringslong


createddata_dir = ['~/babbling/created_data/'];     %data dir
id_dir = [tag,'/'];
outdir = [createddata_dir,id_dir,'network_analysis'];
firingdir = [outdir,'/firing_data'];
pca_dir = [outdir,'/PCA_reservoir'];

if ~exist(pca_dir, 'dir')
    mkdir(pca_dir);
else
    addpath(pca_dir);
end

dime = zeros(2,4);
for j = 1:1
    yoked = YOKED(j,:);
    for i =1:2
        stdp = STDP(i,:);
        %Babbling&mean caliculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ID = [tag,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward,'_',feedbacktype];

           for k=1:d
                %conduct create_firing
                   display([num2str(k),ID]);
                   eighty_dime(k) = dimention([num2str(k),'_',ID]);

           end
         
           dime(1,(j-1)*2+i)=mean(eighty_dime);
           dime(2,(j-1)*2+i)=std(eighty_dime);


      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end


%create dimention file
Name = {[YOKED(1,:),STDP(1,:)];[YOKED(1,:),STDP(2,:)];[YOKED(2,:),STDP(1,:)];[YOKED(2,:),STDP(2,:)];};
ave_dime = [dime(1,1);dime(1,2);dime(1,3);dime(1,4);];
std_dime = [dime(2,1);dime(2,2);dime(2,3);dime(2,4);];
T = table(ave_dime,std_dime,'RowNames',Name);
writetable(T,[pca_dir,'/dimention_',tag,'.csv'],'WriteRowNames',true)



display('END');
%plot FIg with R
