function [highcor] = correlation(commandid,firingid)
global k yoke STDP outdir p dim simutime
mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/correlation']);




%%%%%%%%%%%%%% import data
motorcommand=importdata(commandid);
firingsdata=importdata(firingid);
%%%%%%%%%%%%%%


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
threshold = 0.5;             %threshold of consonant.


%%%%%%%%%%%%%% deta decision
Firings = firings;           %decide the data to be used
lengthdata = size(Firings);   %length of data
%%%%%%%%%%%%%%

%%%%
avemot=zeros(simutime,1);
history=zeros(lengthdata(1,1),1);


%%%consonant vector initializaton
historyvec=zeros(1000,201);
consonant = find(motorcommand(:,2) > threshold);
sizeconsonant = size(consonant);
consonantvec = consonant;
size(find(Firings>0.1))





%%%%%%%%%%%%% create consonant vectar to caliculate correlation
for i=1:100
    consonantvec=[consonantvec consonant+i];
end
for i=1:100
    consonantvec=[consonant-i consonantvec];
end
%%%%%%%%%%%%% size(consonantvec)=[?,201]


%%%%%%%%%%%%% 100ms average
windowSize = 100; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
average_motor = filter(b,a,motorcommand(:,2));
avemot(1:(simutime-windowSize/2+1))=average_motor(windowSize/2:simutime);
%%%%%%%%%%%%%



%%%%%%%%%%%%% arrange data
sabun = diff(avemot);
sabun(find(sabun<0))=0;
sabun(find(sabun>0))=1;
%%%%%%%%%%%%%


sabun1 = find(sabun==1);
size1=size(sabun1);
sabun0 = find(sabun==0);
size0=size(sabun0);




for i=1:lengthdata(1,1)
    for a=1:size1(1,1)
    history(i,1)=history(i,1) + Firings(i,sabun1(a));    
    end
    history(i,1) = history(i,1)/size1(1,1);
end


highcor = find(history(:,1)>0.0016);

%%%%%%%%%%%%%%%%%%%%%%%%% plot histogram
fig555=histogram(history);
fig555.NumBins = 30;
fig555.BinEdges = [0:0.0002:0.003];
fig555.FaceColor = 'k';
fig555.EdgeColor = 'k';
axis([0 0.003 0 400]);
set(gca,'FontSize',16);
ytickformat('%.2f');
saveas(fig555,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/correlation/correlation_p=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'simutime=',num2str(simutime),'.png']);
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    




end




