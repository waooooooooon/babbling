function [] = plotfirings(ID,newT,reinforce,outInd,muscscale,yoke,plotOn,feedbacktime,learningratio,speinplate,STDP,debug,IP,separatephase,Network,reward,feedbacktype)
global id

createddata_dir = ['~/babbling/created_data/'];
datadir = [createddata_dir,id,'/data/'];
wavdir = [datadir,ID, '_Wave'];
firingsdir = [datadir,ID, '_Firings'];
workspacedir = [datadir,ID, '_Workspace'];
setdir = ['setting'];


if newT>1000
    for i = 1000:1000:newT

        A=dlmread([firingsdir,'/R_firings_',ID,'_',num2str(i),'.txt']);
        B=dlmread([firingsdir,'/Mot_firings_',ID,'_',num2str(i),'.txt']);


        %%%%%%%%%%%%%%% plot Resevoir
        nega=find(A(:,3)>800);
        input = find(A(:,3)<201);
        output = find(A(:,3)>200 & A(:,3)<301);
        exi = find(A(:,3)>300 & A(:,3)<801);
        
        hold on;
        fig1 = scatter(A(nega,2),A(nega,3),'.','b'); % Plot all the neurons'' spikes
        fig1 = scatter(A(input,2),A(input,3),'.','k'); % Plot all the neurons'' spikes
        fig1 = scatter(A(output,2),A(output,3),'.','g'); % Plot all the neurons'' spikes
        fig1 = scatter(A(exi,2),A(exi,3),'.','r'); % Plot all the neurons'' spikes
        hold off;
        title('Reservoir Firings', 'fontweight','bold');
        axis([0 1000 0 1000]);
        xlabel('Millisecond');
        ylabel('Neuron Number');
        hold off;
        %legend('inhibitory','input','output','exitatory');


        saveas(fig1,[firingsdir,'/Firings_',ID,'_',num2str(i),'.png']);
        close all;

        %histogram of firings

        for j=1:800
            Br = size(find(A(:,3)==j));
            firingrate(j) = Br(1,1);

        end

        edge = logspace(0, 10, 300);
        hold on;
        h=histogram(firingrate,edge);set(gca,'Xscale','log');xlim([0 200]);
        %h=histogram(firingrate);xlim([0 200]);
        hold off;
        


        saveas(h,[firingsdir,'/Histlog_',ID,'_',num2str(i),'.png']);
        close all;
        h=histogram(firingrate);xlim([0 200]);
        saveas(h,[firingsdir,'/Hist_',ID,'_',num2str(i),'.png']);
        close all;
        %%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%% plot Motor
        negamot=find(B(:,3)>50);
        eximot = find(B(:,3)<51);
        

        hold on;
        fig2 = scatter(B(negamot,2),B(negamot,3),'.','b'); % Plot all the neurons'' spikes
        fig2 = scatter(B(eximot,2),B(eximot,3),'.','r'); % Plot all the neurons'' spikes

        hold off;
        %legend('exitatory','inhibitory');
        %fig1 = scatter(B(:,2),B(:,3),'.','k'); % Plot all the neurons'' spikes
        title('Motor Firings', 'fontweight','bold');
        axis([0 1000 0 100]);
        xlabel('Millisecond');
        ylabel('Neuron Number');


        saveas(fig2,[firingsdir,'/MotFirings_',ID,'_',num2str(i),'.png']);
        close all;

        %histogram of firings

        for j =1:1000
            Bm = size(find(B(:,3)==j));
            firingratemot(j) = Bm(1,1);

        end

        edge = logspace(0, 10, 300);
        %h=histogram(firingratemot,edge);set(gca,'Xscale','log');xlim([0 100]);
        h=histogram(firingratemot);xlim([0 100]);
        saveas(h,[firingsdir,'/MotHist_',ID,'_',num2str(i),'.png']);
        %%%%%%%%%%%%%%%%%%%%%%
    end
else
    i=newT;
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

        for j =1:800
            Br = size(find(A(:,3)==j));
            firingrate(j) = Br(1,1);

        end

        edge = logspace(0, 10, 300);
        h=histogram(firingrate,edge);set(gca,'Xscale','log');xlim([0 200]);
        %h=histogram(firingrate);xlim([0 200]);


        saveas(h,[firingsdir,'/Histlog_',ID,'_',num2str(i),'.png']);
        
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

        for j =1:1000
            Bm = size(find(B(:,3)==j));
            firingratemot(j) = Bm(1,1);

        end

        edge = logspace(0, 10, 300);
        %h=histogram(firingratemot,edge);set(gca,'Xscale','log');xlim([0 100]);
        h=histogram(firingratemot);xlim([0 100]);
        saveas(h,[firingsdir,'/MotHist_',ID,'_',num2str(i),'.png']);
        %%%%%%%%%%%%%%%%%%%%%%
end

end

