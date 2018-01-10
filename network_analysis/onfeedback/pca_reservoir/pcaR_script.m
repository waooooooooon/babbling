%importdata
global tag simutime createddata_dir id_dir outdir firingdir pca_dir

d=5;    %Iterate number
YOKED=['No';'Sc'];   %Sc or No
ploton = 0; % 1 or 0
STDP=['STDP';'NSTD'];
tag=['180102_p0.01'];
IP =['IP'];        %IP or Tonic or afterIP
separatephase = ['notseparate'];      %separatephase or notseparate
Network = ['lattice'];      %lattice or random
reward = ['negativereward'];
feedbacktime=1;
iterate=3000;
speinplate=0.3;
debug=0;
simutime=50000;        %simutime of plotfiringslong


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
for j = 1:2
    yoked = YOKED(j,:);
    for i =1:2
        stdp = STDP(i,:);
        %Babbling&mean caliculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for p=[0.01]
        ID = ['_',tag,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward];
        

           for k=1:d
                %conduct create_firing
                   display([num2str(k),ID]);
                   eighty_dime(k) = dimention([num2str(k),ID]);

           end
         
           dime(1,(j-1)*2+i)=mean(eighty_dime);
           dime(2,(j-1)*2+i)=std(eighty_dime);


        end
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
