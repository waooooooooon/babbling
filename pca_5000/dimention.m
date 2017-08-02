function [y,ruiseki,kiyo]= dimention(id)
global k yoke STDP outdir p dim 
mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d']);
y=0;


firings=importdata(id);

Firings=zeros(1000,1000);

for i=1:1000
  
    I=firings(find(firings(:,2)==i),3); %t=iに発火したニューロンidをIに
    C=size(I);  %Iのサイズ
    for j=1:C(1,1)
    Firings(I(j,1),i)=1;
    end
    
end

NeFirings=Firings(1:800,:);  %NeFirings=NeFirings.'; %invert to caliculate pca

%Firings=Firings.';

%Firings=zscore(Firings);

[COEFF,SCORE,latent] = pca(NeFirings);

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
saveas(fig103,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/PCA3d/Ne_p=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'.png']);


%size(kiyo)
%sum(kiyo)

