function [] = plotfirings(ID,newT,reinforce,outInd,muscscale,yoke,plotOn,feedbacktime,learningratio,speinplate,STDP,debug,IP,separatephase,Network,reward,feedbacktype)
global id

createddata_dir = ['~/babbling/created_data/'];
datadir = [createddata_dir,id,'/data/'];
wavdir = [datadir,ID, '_Wave'];
firingsdir = [datadir,ID, '_Firings'];
workspacedir = [datadir,ID, '_Workspace'];
setdir = ['setting'];



    for i = 1000:1000:newT

        A=dlmread([firingsdir,'/R_firings_',ID,'_',num2str(i),'.txt']);
        B=dlmread([firingsdir,'/Mot_firings_',ID,'_',num2str(i),'.txt']);


        %%%%%%%%%%%%%%% plot Resevoir
        fig1 = scatter(A(:,2),A(:,3),'.','k'); % Plot all the neurons'' spikes
        title('Reservoir Firings', 'fontweight','bold');
        axis([0 1000 0 1000]);
        xlabel('Millisecond');
        ylabel('Neuron Number');


        saveas(fig1,[firingsdir,'/Firings_',ID,'_',num2str(i),'.png']);

        %histogram of firings

        for i =1:1000
            Br = size(find(A(:,3)==i));
            firingrate(i) = Br(1,1);

        end

        edge = logspace(0, 10, 300);
        %h=histogram(firingrate,edge);set(gca,'Xscale','log');xlim([0 100]);
        h=histogram(firingrate);xlim([0 200]);


        saveas(h,[firingsdir,'/Hist_',ID,'_',num2str(i),'.png']);
        %%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%% plot Motor
        fig1 = scatter(B(:,2),B(:,3),'.','k'); % Plot all the neurons'' spikes
        title('Motor Firings', 'fontweight','bold');
        axis([0 1000 0 100]);
        xlabel('Millisecond');
        ylabel('Neuron Number');


        saveas(fig1,[firingsdir,'/MotFirings_',ID,'_',num2str(i),'.png']);

        %histogram of firings

        for i =1:1000
            Bm = size(find(B(:,3)==i));
            firingratemot(i) = Bm(1,1);

        end

        edge = logspace(0, 10, 300);
        %h=histogram(firingratemot,edge);set(gca,'Xscale','log');xlim([0 100]);
        h=histogram(firingratemot);xlim([0 100]);
        saveas(h,[firingsdir,'/MotHist_',ID,'_',num2str(i),'.png']);
        %%%%%%%%%%%%%%%%%%%%%%
    end

end

