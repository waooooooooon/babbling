function [all_sum,input_sum,output_sum,inh_sum,other_sum]=separate_firingrate(id,Noid)
global tag simutime createddata_dir id_dir outdir firingdir pca_dir onfeedbackdir motortype

display('importing data');

    firingsdata=importdata([firingdir,'/firing_onfeedback_',motortype,'_',id,'_',num2str(simutime),'.txt']);
    %firingsdata=importdata([firingdir,'/firing_onfeedback_',id,'_',num2str(simutime),'.txt']);
    outputdir = [onfeedbackdir,'/through_simulation/significant_difference'];
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    else
        addpath(outputdir);
    end
    
  
    
    %%%%%%% import No feedback
    if strcmp(motortype,'feedback')
        motcommanddata=importdata([firingdir,'/motorcommand_onfeedback_',Noid,'_',num2str(simutime),'.txt']);
    %consonant = find(motcommanddata(:,2)>0.5);
    %vowel = find(motcommanddata(:,2)<=0.5);
    elseif strcmp(motortype,'sin')
        muscle_sin = sin(linspace(pi,-pi,1000))*0.9;
        muscle_sin = repmat(muscle_sin,1,simutime/1000);
        motcommanddata =[[1:simutime]',muscle_sin'];
        
    end
    
    
    
    



    
    
display('caliculating the Firings');
%%%%%%%%%%%%%%caliculate the Firings
firings=zeros(1000,simutime);
for i=1:simutime
  
    I=firingsdata(find(firingsdata(:,2)==i),3); %t=i???????????j???[????id??I??
    C=size(I);  %I???T?C?Y
    for j=1:C(1,1)
    firings(I(j,1),i)=1;
    end
    
    
end
%%%%%%%%%%%%%

sizefiring = size(firings);
firings = [firings;motcommanddata(:,2)'];

%%%% 100ŒÂ‚Ìƒf[ƒ^‚É•ªŠ„
Firings = zeros(1001,sizefiring(1,2)/100,100);
consofiringrate = zeros(1000,100);
vowelfiringrate = zeros(1000,100);

display('caliculating consonants and vowels');

for i =1:100
    Firings(:,:,i) = firings(:,1+(i-1)*sizefiring(1,2)/100:(i)*sizefiring(1,2)/100);
end

for i =1:100
    
    sepaconsonant=find(Firings(1001,:,i)>0.5);
    sepavowel=find(Firings(1001,:,i)<0.5);
    consofiringrate(:,i) = mean((Firings(1:1000,sepaconsonant,i)),2);
    vowelfiringrate(:,i) = mean((Firings(1:1000,sepavowel,i)),2); 
end

% firing rate of each state
%consofiringrate = mean((firings(:,consonant)),2);
%vowelfiringrate = mean((firings(:,vowel)),2);



%T
%h = ttest(x,y)
display('caliculating ttest');


for i = 1:1000
     [significant_difference(i,1),significant_difference(i,2)] = ttest2(consofiringrate(i,:),vowelfiringrate(i,:));
end


    
        
%%%%%%%%% 0~200:input neuron 200~300:output neuron 800~1000:inhibitory
%%%%%%%%% neuron
all_sum = sum(significant_difference(:,1));
input_sum = sum(significant_difference(1:200,1));
output_sum = sum(significant_difference(201:300,1));
inh_sum = sum(significant_difference(801:1000,1));
other_sum = sum(significant_difference(301:800,1));


data = [all_sum;input_sum;output_sum;inh_sum;other_sum];

csvwrite([outputdir,'/',motortype,'_significant_difference_',id,'_',num2str(simutime),'.txt'],data);
    

end




