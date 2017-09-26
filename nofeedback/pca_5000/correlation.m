function correlation(commandid,Firings)
global k yoke STDP outdir p dim simutime
mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/correlation']);




%%%%%%%%%%%%%% import command
motorcommand=importdata(commandid);



%%%%%%%%%%%%%% initialization
NeFirings=Firings(1:800,:);  %Exitatory neuron activity
threshold = 0.5;             %threshold of consonant.




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



for a=1:consonant(1,1)
    for i=1:1000
    Firings(i,consonantvec(a,:));
    
    
    
    
    end
end

    

    
    
    




end




