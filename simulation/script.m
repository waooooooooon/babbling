global id

%initialization
d=5;    %Iterate number
YOKED=['No';'Sc'];   %Sc or No
ploton = 0; % 1 or 0
STDP=['STDP';'NSTD'];
id=['test'];
IP =['IP'];        %IP or Tonic or afterIP
separatephase = ['notseparate'];      %separatephase or notseparate
Network = ['lattice'];      %lattice or random
reward = ['negativereward'];
feedbacktype = ['consonant']        %consonant or fft
feedbacktime=1;
iterate=3000;
speinplate=0.3;
debug=0;

saliencedata_dir = ['../salience_analysis/saliency_data/'];



for j = 1:1
    yoked = YOKED(j,:);
    for i =1:1
        stdp = STDP(i,:);
        %Babbling&mean caliculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for p=[0.03]
        ID = ['_',id,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward,'_',feedbacktype];


           for i=1:d
                %conduct babbling
                   display([num2str(i),ID]);

                   babbling([num2str(i),ID],iterate,'reinforce',1:100,4,yoked,ploton,feedbacktime,p,speinplate,stdp,debug,IP,separatephase,Network,reward,feedbacktype);
           end

         clearvars meanout data datahist output;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end



display('END');
%plot FIg with R
