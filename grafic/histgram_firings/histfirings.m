global id itenumber

%initialization
d=1;    %Iterate number
YOKED=['No';'Sc'];   %Sc or No
ploton = 0; % 1 or 0
STDP=['STDP';'NSTD'];
id=['180120_Sctime'];
IP =['LiIP'];        %threIP or Tonic or NoIP or LiIP
separatephase = ['randSc'];      %separate or nseparate or randSc
Network = ['random'];      %lattice or random
reward = ['normal'];       %nega or normal
feedbacktype = ['fft'];        %consonant or fft or none
feedbacktime=1;
iterate=2000;
speinplate=0.3;
debug=0;
p = 0.03;
created_data = ['../created_data/'];

createddata_dir = ['~/babbling/created_data/'];
datadir = [createddata_dir,id,'/data/'];



color = ['r';'k'];
feedback = ['Sc']; %聴覚フィードバック: No サロゲート入力: Sc です．

%caliculate mean
for j = 1:2
    stdp = STDP(j,:);
    col = color(j,:);

 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        NoID = [id,'_',num2str(iterate),'_reinforce_100_4_',feedback,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward,'_',feedbacktype];

        
        Nofiringsdir = [datadir,'1_',NoID, '_Firings'];
        
        
        No=dlmread([Nofiringsdir,'/R_firings_','1_',NoID,'_',num2str(iterate),'.txt']);
        

        %histogram of firings

        for j =1:800
            NoBr = size(find(No(:,3)==j));
            Nofiringrate(j) = NoBr(1,1);
            
        end

        edge = logspace(0, 10, 300);
        hold on;
        h=histogram(Nofiringrate,edge);set(gca,'Xscale','log');xlim([0 200]);
         
        h.FaceColor = col;
        h.EdgeColor = 'k';
        saveas(h,['./Hist.png']);
        

end



