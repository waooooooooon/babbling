function [] = plotfirings(ID,newT,reinforce,outInd,muscscale,yoke,plotOn,feedbacktime,learningratio,speinplate,STDP,debug,IP,separatephase,Network,reward,feedbacktype)
global id

createddata_dir = ['~/babbling/created_data/'];
datadir = [createddata_dir,id,'/data/'];
wavdir = [datadir,ID, '_Wave'];
firingsdir = [datadir,ID, '_Firings'];
workspacedir = [datadir,ID, '_Workspace'];
setdir = ['setting'];



A=dlmread([firingsdir,'/babble_daspnet_firings_',ID,'_',num2str(newT),'.txt']);

fig1 = scatter(A(:,2),A(:,3),'.','k'); % Plot all the neurons'' spikes
title('Reservoir Firings', 'fontweight','bold');
axis([0 1000 0 1000]);
xlabel('Millisecond');
ylabel('Neuron Number');


saveas(fig1,[firingsdir,'/Firings_',ID,'_',num2str(newT),'.png']);

%histogram of firings

for i =1:1000
    B = size(find(A(:,3)==i));
    firingrate(i) = B(1,1);
    
end

edge = logspace(0, 10, 300);
h=histogram(firingrate,edge);set(gca,'Xscale','log');xlim([0 100]);

saveas(h,[firingsdir,'/Hist_',ID,'_',num2str(newT),'.png']);





