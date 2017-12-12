%initialization
d=1;    %Iterate number
YOKED=['No';'Sc'];   %Sc or No
ploton = 0; % 1 or 0
STDP=['STDP';'NSTD'];
id=['1701212_long'];
IP =['IP'];        %IP or Tonic or afterIP
separatephase = ['notseparate'];      %separatephase or notseparate
Network = ['lattice'];
reward = ['negativereward'];
feedbacktime=1;
iterate=6000;
speinplate=0.3;
debug=0;

saliencedata_dir = ['../salience_analysis/saliency_data/'];

%Create mean file
meandir=['mean'];
if ~exist(meandir,'dir')
 mkdir(meandir);
end
mkdir([meandir,'/',meandir],'dir');

for j = 1:2
    yoked = YOKED(j,:);
    for i =1:2
        stdp = STDP(i,:);
        %Babbling&mean caliculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for p=[0.03]
        ID = ['_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward];
          mkdir([meandir,'/p=',num2str(p)]);

           for i=1:d
                %conduct babbling
                   display([num2str(i),ID]);

                   babbling([num2str(i),ID],iterate,'reinforce',1:100,4,yoked,ploton,feedbacktime,p,speinplate,stdp,debug,IP,separatephase,Network,reward);

                   copyfile([num2str(i),ID,'_Workspace/',num2str(i),ID,'.csv'],[meandir,'/p=',num2str(p),'/']);
           end

         %caliculate mean
         %cd [meandir,'/p=',num2str(p),'/'];
        display(['Calculating mean of p=',num2str(p)]);
        lenght=size(importdata([meandir,'/p=',num2str(p),'/',num2str(1),ID,'.csv']));
        datahist=zeros(lenght);
%{
          for i=1:d
                  data(:,:,i)=importdata([meandir,'/p=',num2str(p),'/',num2str(i),ID,'.csv']);
                  datahist(:,:)=datahist(:,:)+data(:,:,i);
          end

         meanout=datahist/d;
         output=[meanout(:,1),data(:,2,1)];
         csvwrite([meandir,'/p=',num2str(p),'/mean_p=',num2str(p),'.csv'],output);


         copyfile([meandir,'/p=',num2str(p),'/mean_p=',num2str(p),'.csv'],[meandir,'/',meandir,'/']);
        %}

         clearvars meanout data datahist output;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

mkdir([saliencedata_dir,id]);




display('END');
%plot FIg with R
