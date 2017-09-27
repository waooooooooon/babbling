N = 1000;
Nmot = 100;
Nout = 100;
simutime = 100000;

I_randominput = 13*(rand(N,simutime)-0.5);
Imot_randominput=13*(rand(Nmot,simutime)-0.5);
sout = rand(Nout,Nmot); 

   
dlmwrite('I.txt',I_randominput);
dlmwrite('Imot_randominput.txt',Imot_randominput);
dlmwrite('sout.txt',sout);