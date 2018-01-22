%importdata
global tag simutime itenumber createddata_dir id_dir outdir firingdir pca_dir nofeedbackdir

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
firingdir = [outdir,'/nofeedback_firing_data'];
pca_dir = [outdir,'/PCA_reservoir'];
nofeedbackdir = [pca_dir,'/nofeedback'];

if ~exist(nofeedbackdir, 'dir')
    mkdir(nofeedbackdir);
else
    addpath(nofeedbackdir);
end






for j = 1:2
    yoked = YOKED(j,:);
    for i =1:2
        stdp = STDP(i,:);
        %create firings
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ID = [tag,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward,'_',feedbacktype];


           for k=1:d
                %conduct create_firing
                   display([num2str(k),ID]);
                   %create_firingnofeedback([num2str(k),'_',ID],iterate,'reinforce',1:100,4,yoked,ploton,feedbacktime,p,speinplate,stdp,debug,IP,separatephase,Network,reward,feedbacktype);
                   

           end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end



dime = zeros(2,4);
for j = 1:2
    yoked = YOKED(j,:);
    for i =1:2
        stdp = STDP(i,:);
        %caliculate pca etc..
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
writetable(T,[nofeedbackdir,'/dimention_',tag,'.csv'],'WriteRowNames',true)



display('END');
%plot FIg with R