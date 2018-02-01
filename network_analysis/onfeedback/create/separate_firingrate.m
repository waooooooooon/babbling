function [] = separate_firingrate(id,Noid)
global tag simutime createddata_dir id_dir outdir firingdir pca_dir onfeedbackdir motortype

display('importing data');

    %firingsdata=importdata([firingdir,'/firing_onfeedback_',motortype,'_',id,'_',num2str(simutime),'.txt']);
    firingsdata=importdata([firingdir,'/firing_onfeedback_feedback_',id,'_',num2str(simutime),'.txt']);
    %firingsdata=importdata([firingdir,'/firing_onfeedback_',id,'_',num2str(simutime),'.txt']);
    outputdir = [onfeedbackdir,'/through_simulation/significant_difference'];
    datadir = [outputdir,'/matdata'];
    workspacefilename = [datadir,'/',id,'.mat'];
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    else
        addpath(outputdir);
    end
    if ~exist(datadir, 'dir')
        mkdir(datadir);
    else
        addpath(datadir);
    end
    
  
    
    %%%%%%% import No feedback
    if strcmp(motortype,'feedback')
        motcommanddata=importdata([firingdir,'/motorcommand_onfeedback_',motortype,'_',Noid,'_',num2str(simutime),'.txt']);
        %motcommanddata=importdata([firingdir,'/motorcommand_onfeedback_',motortype,'_',id,'_',num2str(simutime),'.txt']);%inportminemotorcommand
        
    %consonant = find(motcommanddata(:,2)>0.5);
    %vowel = find(motcommanddata(:,2)<=0.5);
    elseif strcmp(motortype,'sin')
        muscle_sin = sin(linspace(pi,-pi,500))*0.9;
        muscle_sin = repmat(muscle_sin,1,simutime/500);
        motcommanddata =[[1:simutime]',muscle_sin'];
    elseif strcmp(motortype,'output')
        %motcommanddata=importdata([firingdir,'/motorcommand_onfeedback_',motortype,'_',id,'_',num2str(simutime),'.txt']);
        motcommanddata=importdata([firingdir,'/motorcommand_onfeedback_feedback_',id,'_',num2str(simutime),'.txt']);
        
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
openfiringrate = zeros(1000,100);
closefiringrate = zeros(1000,100);


display('caliculating consonants and vowels');

for i =1:100
    Firings(:,:,i) = firings(:,1+(i-1)*sizefiring(1,2)/100:(i)*sizefiring(1,2)/100);
end


%firingrate of consonant and vowel 
for i =1:100
    
    sepaconsonant=find(Firings(1001,:,i)>0.5);
    sepavowel=find(Firings(1001,:,i)<0.5);
    consofiringrate(:,i) = mean((Firings(1:1000,sepaconsonant,i)),2);
    vowelfiringrate(:,i) = mean((Firings(1:1000,sepavowel,i)),2); 
end


for i =1:100
    
    openmouth=find(diff(Firings(1001,:,i))>0);
    closemouth=find(diff(Firings(1001,:,i))<0);
    openfiringrate(:,i) = mean((Firings(1:1000,openmouth,i)),2);
    closefiringrate(:,i) = mean((Firings(1:1000,closemouth,i)),2); 
end

for i =1:100
    
    lowmotor=find(Firings(1001,:,i)<-0.5);
    otherhighmotor=find(Firings(1001,:,i)>-0.5);
    lowmotorfiringrate(:,i) = mean((Firings(1:1000,lowmotor,i)),2);
    otherhighmotorfiringrate(:,i) = mean((Firings(1:1000,otherhighmotor,i)),2); 
end



for i =1:100
    
    middlemotor = find(Firings(1001,:,i)>-0.5);
    middlemotor = find(Firings(1001,middlemotor,i)<0.5);
    other1=find(Firings(1001,:,i)>0.5);
    other2 = find(Firings(1001,:,i)<-0.5);
    others = [other1,other2];
    middlefiringrate(:,i) = mean((Firings(1:1000,middlemotor,i)),2);
    othersfiringrate(:,i) = mean((Firings(1:1000,others,i)),2); 
end





% firing rate of each state
%consofiringrate = mean((firings(:,consonant)),2);
%vowelfiringrate = mean((firings(:,vowel)),2);



%T
%h = ttest(x,y)
display('caliculating ttest');


for i = 1:1000
     [significant_difference_conso(i,1),significant_difference_conso(i,2)] = ttest2(consofiringrate(i,:),vowelfiringrate(i,:));
end

for i = 1:1000
     [significant_difference_move(i,1),significant_difference_move(i,2)] = ttest2(openfiringrate(i,:),closefiringrate(i,:));
end

for i = 1:1000
     [significant_difference_low(i,1),significant_difference_low(i,2)] = ttest2(lowmotorfiringrate(i,:),otherhighmotorfiringrate(i,:));
end

for i = 1:1000
     [significant_difference_middle(i,1),significant_difference_middle(i,2)] = ttest2(middlefiringrate(i,:),othersfiringrate(i,:));
end

        
%%%%%%%%% 0~200:input neuron 200~300:output neuron 800~1000:inhibitory
%%%%%%%%% neuron

%%%%%%%%%%%%conso
A =[mean(consofiringrate,2)-mean(vowelfiringrate,2),significant_difference_conso];
A_sub = size(find(A(find(A(:,2)==1),1)>0));

all_sum1 = sum(significant_difference_conso(:,1));
input_sum1 = sum(significant_difference_conso(1:200,1));
output_sum1 = sum(significant_difference_conso(201:300,1));
inh_sum1 = sum(significant_difference_conso(801:1000,1));
other_sum1 = sum(significant_difference_conso(301:800,1));
data1 = [all_sum1;input_sum1;output_sum1;inh_sum1;other_sum1;A_sub(1,1)];


csvwrite([outputdir,'/',motortype,'_difference_conso_',id,'_',num2str(simutime),'.txt'],data1);
csvwrite([outputdir,'/',motortype,'_A_conso_',id,'_',num2str(simutime),'.txt'],A);


%%%%%%%%%%%%  move
B =[mean(closefiringrate,2)-mean(openfiringrate,2),significant_difference_move];
B_sub = size(find(B(find(B(:,2)==1),1)>0));

all_sum2 = sum(significant_difference_move(:,1));
input_sum2 = sum(significant_difference_move(1:200,1));
output_sum2 = sum(significant_difference_move(201:300,1));
inh_sum2 = sum(significant_difference_move(801:1000,1));
other_sum2 = sum(significant_difference_move(301:800,1));
data2 = [all_sum2;input_sum2;output_sum2;inh_sum2;other_sum2;B_sub(1,1)];

csvwrite([outputdir,'/',motortype,'_difference_move_',id,'_',num2str(simutime),'.txt'],data2);
csvwrite([outputdir,'/',motortype,'_B_move_',id,'_',num2str(simutime),'.txt'],B);


%%%%%%%%%%%% low motor
D =[mean(lowmotorfiringrate,2)-mean(otherhighmotorfiringrate,2),significant_difference_low];
D_sub = size(find(D(find(D(:,2)==1),1)>0));


all_sum3 = sum(significant_difference_low(:,1));
input_sum3 = sum(significant_difference_low(1:200,1));
output_sum3 = sum(significant_difference_low(201:300,1));
inh_sum3 = sum(significant_difference_low(801:1000,1));
other_sum3 = sum(significant_difference_low(301:800,1));
data3 = [all_sum3;input_sum3;output_sum3;inh_sum3;other_sum3;D_sub(1,1)];

csvwrite([outputdir,'/',motortype,'_difference_lowmotor_',id,'_',num2str(simutime),'.txt'],data3);
csvwrite([outputdir,'/',motortype,'_D_low_',id,'_',num2str(simutime),'.txt'],D);



%%%%%%%%%%%%  middle
E =[mean(middlefiringrate,2)-mean(othersfiringrate,2),significant_difference_middle];
E_sub = size(find(B(find(E(:,2)==1),1)>0));

all_sum4 = sum(significant_difference_middle(:,1));
input_sum4 = sum(significant_difference_middle(1:200,1));
output_sum4 = sum(significant_difference_middle(201:300,1));
inh_sum4 = sum(significant_difference_middle(801:1000,1));
other_sum4 = sum(significant_difference_middle(301:800,1));
data4 = [all_sum4;input_sum4;output_sum4;inh_sum4;other_sum4;E_sub(1,1)];

csvwrite([outputdir,'/',motortype,'_difference_middle_',id,'_',num2str(simutime),'.txt'],data4);
csvwrite([outputdir,'/',motortype,'_E_middle_',id,'_',num2str(simutime),'.txt'],E);


moveandconso=sum(A(find(B(:,2)==1),2));



Name = {['consonant'];['move(up)'];['low motor(vowel)'];['middle'];};
all_neuron = [all_sum1;all_sum2;all_sum3;all_sum4;];
input = [input_sum1;input_sum2;input_sum3;input_sum4;];
output = [output_sum1;output_sum2;output_sum3;output_sum4;];
inhibitory = [inh_sum1;inh_sum2;inh_sum3;inh_sum4;];
others = [other_sum1;other_sum2;other_sum3;other_sum4;];
which_mode = [A_sub(1,1);B_sub(1,1);E_sub(1,1);D_sub(1,1);];
T = table(all_neuron,input,output,inhibitory,others,which_mode,'RowNames',Name);
writetable(T,[outputdir,'/',motortype,'_Difference_',id,'_',num2str(simutime),'.csv'],'WriteRowNames',true);



save(workspacefilename);
    
    

end




