function [y,ruiseki,kiyo,transfer_score,collelation]= dimention_correlation(id,motorid)
global k yoke STDP outdir p dim simutime
mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d']);
y=0;
consonant = 0.4;   %threshold of consonant
threshold = 0.01;   %threshold of correlation




%%%%%%%%%%%%%% import firings
firings=importdata(id);
Firings=zeros(1000,simutime);


%%%%%%%%%%%%%%%import motor command
motor = importdata(motorid);


%%%%%%%%%%%%%%caliculate the Firings
for i=1:simutime
  
    I=firings(find(firings(:,2)==i),3); %t=i???????????j???[????id??I??
    C=size(I);  %I???T?C?Y
    for j=1:C(1,1)
    Firings(I(j,1),i)=1;
    end
    
end

NeFirings=Firings(1:800,:);  %NeFirings=NeFirings.'; %invert to caliculate pca
NeFirings_conv=NeFirings.';  %invert NeFirings
%Firings=Firings.';

%Firings=zscore(Firings);






%%%%%%%%%%%%%caliculate the PCA
[COEFF,SCORE,latent] = pca(NeFirings_conv);   %conduct pca

%NeFirings 800,5000  need to convert
%pca(obserbed_data,variable)



% initialization of transfer score
sizescore=size(SCORE);
size_score=sizescore(1,1);    %initialization ofsize of SCORE
%transfer_score=zeros(size_score);
transfer_score=0;   %initialization of transfer_score
transfer_score_A=0;   %initialization of transfer_score_A
transfer_history=zeros(size_score,2);
transfer_history_A=zeros(size_score,2);


%%%%%%%%%%%%%%%%% initialization and caliculate the contribution rate of pca
latentsize=size(latent);
latent(950:latentsize(1,1),:)=[];
kiyo=latent/sum(latent);
ruiseki=zeros(size(kiyo));

for i=1:size(kiyo)
    if i==1
    ruiseki(1)=kiyo(i);
    
    else
        ruiseki(i)=ruiseki(i-1)+kiyo(i);
    end
    
end



%%%%%%%%%%%%%%%%%%%% caliculate the correlation
motor(find(motor(:,2)>consonant),2)=1;
motor(find(motor(:,2)<=consonant),2)=0;

% initialization
neuronsize = size(Firings);
neuron_size = neuronsize(1,1);
neuron = zeros(1,neuron_size);

motor=motor(:,2)';

for i=1:neuron_size  % all neurons
neuron(1,i) = corr2(Firings(i,:),motor);   %���֕���
end



numberofconsonantneuron = size(find(neuron > threshold));

collelation = numberofconsonantneuron(1,2);

%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%   %caliculate transfer_score in the principal component space
for i=1:(size_score-1)  
    temp_score=((SCORE(i,1)-SCORE(i+1,1))^2+(SCORE(i,2)-SCORE(i+1,2))^2+(SCORE(i,3)-SCORE(i+1,3))^2)^0.5;
    transfer_score = transfer_score + temp_score;
    transfer_history(i,1)=i; transfer_history(i,2)=temp_score;
    
end
%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%   %caliculate transfer_score of firings  
for i=1:simutime-1 
    temp_score_A=((NeFirings(:,i)-NeFirings(:,i+1))'*(NeFirings(:,i)-NeFirings(:,i+1)))^0.5;
    transfer_score_A = transfer_score_A + temp_score_A;
    transfer_history_A(i,1)=i; transfer_history_A(i,2)=temp_score_A;
    
end
%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%     plot pca
y=min(find(ruiseki>0.8));
%PCA3d = figure(103);
fig103=plot3(SCORE(:,1),SCORE(:,2),SCORE(:,3));
saveas(fig103,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d/Ne_p=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'simutime=',num2str(simutime),'transferscore=',num2str(transfer_score),'.png']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%% plot transfer_history
fig204=plot(transfer_history(1:size_score-1,1),transfer_history(1:size_score-1,2));
ylim([0,2]);
saveas(fig204,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d/transferhis_p=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'simutime=',num2str(simutime),'transferscore=',num2str(transfer_score),'.png']);
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% plot transfer_history_A
fig205=plot(transfer_history_A(1:simutime-1,1),transfer_history_A(1:simutime-1,2));
ylim([0,5]);
saveas(fig205,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d/Eucliddistancd_p=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'simutime=',num2str(simutime),'transferscore=',num2str(transfer_score),'.png']);
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% fft of transfer_history
[f,p1]= fft_script(transfer_history(1:size_score-1,2),0);

fig203=plot(f,p1);
xlabel('Frequency');
ylabel('Power');
xlim([0,10]);
saveas(fig203,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d/fft_of_transfer_p=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'simutime=',num2str(simutime),'transferscore=',num2str(transfer_score),'.png']);
%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% fft of transfer_history_A
[f,p1]= fft_script(transfer_history_A(1:simutime-1,2),0);

fig201=plot(f,p1);
xlabel('Frequency');
ylabel('Power');
xlim([0,10]);
saveas(fig201,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d/fft_of_transfer_p_A=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'simutime=',num2str(simutime),'transferscore=',num2str(transfer_score),'.png']);
%%%%%%%%%%%%%%%%%%%%%%%%
%size(kiyo)
%sum(kiyo)
