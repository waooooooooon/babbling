global tag simutime createddata_dir id_dir outdir firingdir pca_dir onfeedbackdir motortype yoked stdp ploton feedbacktime p speinplate IP separatephase Network reward feedbacktype d




outputdir = [onfeedbackdir,'/through_simulation/ago_antago'];
if ~exist(outputdir, 'dir')
    mkdir(outputdir);
else
    addpath(outputdir);
end


for k = 1:d

    ID = [tag,'_',num2str(iterate),'_reinforce_100_4_',yoked,'_',num2str(ploton),'_',num2str(feedbacktime),'_',num2str(p),'_',num2str(speinplate),'_',stdp,'_',IP,'_',separatephase,'_',Network,'_',reward,'_',feedbacktype];
    motfiringsdata=importdata([firingdir,'/motfiring_onfeedback_',motortype,'_',num2str(k),'_',ID,'_',num2str(simutime),'.txt']);

    display([num2str(k),ID]);
    %%%%%%%%%%%%%%caliculate the Firings
    firings=zeros(100,simutime);
    for i=1:simutime

        I=motfiringsdata(find(motfiringsdata(:,2)==i),3); %t=i???????????j???[????id??I??
        C=size(I);  %I???T?C?Y
        for j=1:C(1,1)
        firings(I(j,1),i)=1;
        end


    end
    %%%%%%%%%%%%%


    pos = sum(firings(1:50,:));
    nega = sum(firings(51:end,:));


    poshist(1,k) = 1000*sum(pos)/simutime;
    poshist(2,k) = var(pos);
    negahist(1,k) =1000*sum(nega)/simutime;
    negahist(2,k) = var(nega);
end



meanpos= mean(poshist,2);
meannega= mean(negahist,2);
std_pos(1,1) = std(poshist(1,:));
std_pos(2,1) = std(poshist(2,:));
std_nega(1,1) = std(negahist(1,:));
std_nega(2,1) = std(negahist(2,:));



Name = {['posmean'];['posvariance'];['negamean'];['negavariance']};
mean = [meanpos(1,1);meanpos(2,1);meannega(1,1);meannega(2,1);];
std = [std_pos(1,1);std_pos(2,1);std_nega(1,1);std_nega(2,1);];

T = table(mean,std,'RowNames',Name);
writetable(T,[outputdir,'/',motortype,'_Difference_',ID,'_',num2str(simutime),'.csv'],'WriteRowNames',true);

A = table(poshist,negahist);
writetable(A,[outputdir,'/allhist',motortype,'_Difference_',ID,'_',num2str(simutime),'.csv']);





