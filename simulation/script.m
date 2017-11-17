%initialization
d=1;    %Iterate number
yoked=['Yoked'];   %Yoked or NY
stdp=['STDP'];
id=['170222'];
IP =['IP'];
separatephase = ['separatephase'];
Network = ['lattice'];
feedbacktime=1;
iterate=3000;
speinplate=0.1;
debug=1;


%Create mean file
meandir=['mean'];
if ~exist(meandir,'dir')
 mkdir(meandir);
end
mkdir([meandir,'/',meandir],'dir');

%Babbling&mean caliculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p=[0.05]
ID = ['_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network];
  mkdir([meandir,'/p=',num2str(p)]);

   for i=1:d
        %conduct babbling
           display([num2str(i),ID]);

           babbling([num2str(i),ID],iterate,'reinforce',1:100,4,yoked,1,feedbacktime,p,speinplate,stdp,debug,IP,separatephase,Network);

           copyfile([num2str(i),ID,'_Workspace/',num2str(i),ID,'.csv'],[meandir,'/p=',num2str(p),'/']);
   end

 %caliculate mean
 %cd [meandir,'/p=',num2str(p),'/'];
display(['Calculating mean of p=',num2str(p)]);
lenght=size(importdata([meandir,'/p=',num2str(p),'/',num2str(1),ID,'.csv']));
datahist=zeros(lenght);

  for i=1:d
          data(:,:,i)=importdata([meandir,'/p=',num2str(p),'/',num2str(i),ID,'.csv']);
          datahist(:,:)=datahist(:,:)+data(:,:,i);
  end

 meanout=datahist/d;
 output=[meanout(:,1),data(:,2,1)];
 csvwrite([meandir,'/p=',num2str(p),'/mean_p=',num2str(p),'.csv'],output);


 copyfile([meandir,'/p=',num2str(p),'/mean_p=',num2str(p),'.csv'],[meandir,'/',meandir,'/']);

 clearvars meanout data datahist output;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('END');
%plot FIg with R
