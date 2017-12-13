id = ['1_1701212_long_6000_reinforce_100_4_No_0_1_0.03_0.3_STDP_IP_notseparate_lattice_negativereward_Firings/'];
dir = ['firingdata/',id];
filename = ['babble_daspnet_firings_1_1701212_long_6000_reinforce_100_4_No_0_1_0.03_0.3_STDP_IP_notseparate_lattice_negativereward_6000.txt']


A=dlmread([dir,filename]);

fig1 = scatter(A(:,2),A(:,3),'.','k'); % Plot all the neurons'' spikes
title('Reservoir Firings', 'fontweight','bold');
axis([0 1000 0 1000]);
xlabel('Millisecond');
ylabel('Neuron Number');


saveas(fig1,[dir,'/',filename,'.png']);

%histogram of firings

for i =1:1000
    B = size(find(A(:,3)==i));
    firingrate(i) = B(1,1);
    
end

edge = logspace(0, 10, 300);
h=histogram(firingrate,edge);set(gca,'Xscale','log');xlim([0 100]);

saveas(h,[dir,'/',filename,'_histogram.png']);





