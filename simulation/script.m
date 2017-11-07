%initialization
d=5;    %Iterate number
yoked=['NY'];   %Yoked or NY
stdp=['STDP'];
id=['170222'];
IP =['IP'];
separatephase = ['separatephase'];
feedbacktime=1;
iterate=1000;
speinplate=0.5;
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

  mkdir([meandir,'/p=',num2str(p)]);

   for i=1:d
        %conduct babbling
           display([num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase]);

           babbling([num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,IP,separatephase],iterate,'reinforce',1:100,4,yoked,1,feedbacktime,p,speinplate,stdp,debug,IP,separatephase);

           copyfile([num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_Workspace/',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'.csv'],[meandir,'/p=',num2str(p),'/']);
   end

 %caliculate mean
 %cd [meandir,'/p=',num2str(p),'/'];
display(['Calculating mean of p=',num2str(p)]);
lenght=size(importdata([meandir,'/p=',num2str(p),'/',num2str(1),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'.csv']));
datahist=zeros(lenght);

  for i=1:d
          data(:,:,i)=importdata([meandir,'/p=',num2str(p),'/',num2str(i),'_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_1_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'.csv']);
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
