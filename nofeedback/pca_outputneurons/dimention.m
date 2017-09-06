function [trans,y,ruiseki,kiyo]= dimention(id)
global k yoke STDP outdir p dim 
mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d']);
y=0;


firings=importdata(id);

Firings=zeros(1000,5000);

for i=1:5000
  
    I=firings(find(firings(:,2)==i),3); %t=i???????????j???[????id??I??
    C=size(I);  %I???T?C?Y
    for j=1:C(1,1)
    Firings(I(j,1),i)=1;
    end
    
end

outFirings=Firings(201:300,:);  %
outFirings=outFirings.'; %invert to caliculate pca

%Firings=Firings.';

%Firings=zscore(Firings);

[COEFF,SCORE,latent] = pca(outFirings);


sizescore=size(SCORE);
size_score=sizescore(1,1);    %initialization ofsize of SCORE
%transfer_score=zeros(size_score);
transfer_score=0;   %initialization of transfer_score


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

y=min(find(ruiseki>0.8));
%PCA3d = figure(103);
fig103=plot3(SCORE(:,1),SCORE(:,2),SCORE(:,3));

for i=1:(size_score-1)     %caliculate transfer_score
    
    transfer_score = transfer_score + ((SCORE(i,1)-SCORE(i+1,1))^2+(SCORE(i,2)-SCORE(i+1,2))^2+(SCORE(i,3)-SCORE(i+1,3))^2)^0.5;
end

trans=transfer_score;

saveas(fig103,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d/Out_p=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'transferscore=',num2str(transfer_score),'.png']);

%size(kiyo)
%sum(kiyo)

