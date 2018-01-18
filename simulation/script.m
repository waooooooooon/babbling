global id

%initialization
d=4;    %Iterate number
YOKED=['No';'Sc'];   %Sc or No
ploton = 0; % 1 or 0
STDP=['STDP';'NSTD'];
id=['180118_Li=0DA0.001'];
IP =['LiIP'];        %threIP or Tonic or NoIP or LiIP
separatephase = ['nseparate'];      %separate or nseparate
Network = ['random'];      %lattice or random
reward = ['normal'];        %nega or normal
feedbacktype = ['fft'];        %consonant or fft or none
feedbacktime=1;
iterate=2000;
speinplate=0.3;
debug=0;
p = 1;
created_data = ['../created_data/'];


%caliculate babbling
for j = 1:1
    yoked = YOKED(j,:);
    for i =1:1
        stdp = STDP(i,:);
        %Babbling&mean caliculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ID = ['_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward,'_',feedbacktype];


           for i=1:d
                %conduct babbling
                   display([num2str(i),'_',ID]);

                   babbling([num2str(i),ID],iterate,'reinforce',1:100,4,yoked,ploton,feedbacktime,p,speinplate,stdp,debug,IP,separatephase,Network,reward,feedbacktype);
           end

         clearvars meanout data datahist output;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%caliculate mean
for j = 1:2
    yoked = YOKED(j,:);
    for i =2:2
        stdp = STDP(i,:);
        %Babbling&mean caliculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ID = [id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward,'_',feedbacktype];

           for k=1:d
                %conduct babbling
                   display([num2str(k),'_',ID]);
                   workspacedir = [num2str(k),'_',ID,'_Workspace/'];
                   
                   if ~strcmp(feedbacktype,'none')
                   salhist{k}=csvread([created_data,id,'/data/',workspacedir,num2str(k),'_',ID,'.csv']);
                   end
                   
                   plotfirings([num2str(k),'_',ID],iterate,'reinforce',1:100,4,yoked,ploton,feedbacktime,p,speinplate,stdp,debug,IP,separatephase,Network,reward,feedbacktype);
                   
           end
           
           if ~strcmp(feedbacktype,'none')
               for k=1:d
                   if k==1
                       mean_sal = salhist{k};
                   else
                       mean_sal = mean_sal + salhist{k};
                   end
               end
               mean_sal = mean_sal/d;
               csvwrite([created_data,id,'/csv/',ID,'.csv'],mean_sal); 
           end
           
         

         clearvars salhist mean_sal;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end


display('END');
%plot FIg with R
