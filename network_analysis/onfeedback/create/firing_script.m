%importdata
global tag simutime 

d=5;    %Iterate number
YOKED=['No';'Sc'];   %Sc or No
ploton = 0; % 1 or 0
STDP=['STDP';'NSTD'];
tag=['1701231_p0.07'];
IP =['IP'];        %IP or Tonic or afterIP
separatephase = ['notseparate'];      %separatephase or notseparate
Network = ['lattice'];      %lattice or random
reward = ['negativereward'];
feedbacktime=1;
iterate=3000;
speinplate=0.3;
debug=0;

%for firing
simutime=50000;        %simutime of plotfiringslong



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
        for p=[0.07]
        ID = ['_',tag,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward];
        

           for i=1:d
                %conduct create_firing
                   display([num2str(i),ID]);
                   create_firing([num2str(i),ID],iterate,'reinforce',1:100,4,yoked,ploton,feedbacktime,p,speinplate,stdp,debug,IP,separatephase,Network,reward);

           end


        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end



display('END');
%plot FIg with R