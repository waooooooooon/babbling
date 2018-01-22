function [historycor] = correlation(id)
global tag simutime createddata_dir id_dir outdir firingdir pca_dir onfeedbackdir

    firingsdata=importdata([firingdir,'/firing_onfeedback_',id,'_',num2str(simutime),'.txt']);
    outputdir = [onfeedbackdir,'/through_simulation/correlation'];
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    else
        addpath(outputdir);
    end
  
    setdir = ['~/babbling/simulation/setting'];
    table=importdata([setdir,'/table_notnormalization_from0.5to1.5KHz.csv']);
%%%%%%%% create feedback
muscle_his = 10000*sin(linspace(pi,-pi,1000))*0.8+10000;
muscle_his = repmat(muscle_his,1,10);
muscle_his = round(muscle_his);



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





%%%%%%%%%%%%%% initialization
NeFirings=firings(1:800,:);  %Exitatory neuron activity
threshold = 17990;             %threshold of consonant.


%%%%%%%%%%%%%% deta decision
Firings = firings;           %decide the data to be used
lengthdata = size(Firings);   %length of data
%%%%%%%%%%%%%%

%%%%
avemot=zeros(simutime,1);
history=zeros(lengthdata(1,1),1);
historycor=zeros(lengthdata(1,1),1);

%%%consonant vector initializaton
historyvec=zeros(1000,201);
consonant = find(muscle_his > threshold);
sizeconsonant = size(consonant);
consonantvec = consonant;
size(find(Firings>0.1))



muscle_conso = muscle_his;
muscle_conso = muscle_conso > threshold;

%%%%%%%%%%%%% caliculate correration
for i=1:lengthdata(1,1)
    historycor(i,1) = abs(corr2(Firings(i,:),muscle_conso));
end
%highcor = find(historycor(:,1)>0.02 || historycor(:,1)<-0.02 )
%%%%%%%%%%%%%

%%%%%%%% plot hist for pwp(normal)
fig455=histogram(historycor);
fig455.NumBins = 30;
fig455.BinEdges = [-0.06:0.001:0.06];
axis([-0.03 0.06 0 150]);
saveas(fig455,[outputdir,'/correlation_',id,'_',num2str(simutime),'.png']);


end




