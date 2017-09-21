
% initializaton
!rm -r result
!mkdir result



%%%%%%%%% files nemes

file1='../../nofeedback/pca_5000/output/csv/dim_p=0.05_NY_NSTDP.csv';
file2='../../nofeedback/pca_5000/output/csv/dim_p=0.05_Yoked_NSTDP.csv';
file3='../../nofeedback/pca_5000/output/csv/dim_p=0.05_NY_STDP.csv';
file4='../../nofeedback/pca_5000/output/csv/dim_p=0.05_Yoked_STDP.csv';

file11='../../onfeedback/pca_5000/output/csv/dim_p=0.05_NY_NSTDP.csv';
file12='../../onfeedback/pca_5000/output/csv/dim_p=0.05_Yoked_NSTDP.csv';
file13='../../onfeedback/pca_5000/output/csv/dim_p=0.05_NY_STDP.csv';
file14='../../onfeedback/pca_5000/output/csv/dim_p=0.05_Yoked_STDP.csv';

file21='../../nofeedback/pca_outputneurons/output/csv/dim_p=0.05_NY_NSTDP.csv';
file22='../../nofeedback/pca_outputneurons/output/csv/dim_p=0.05_Yoked_NSTDP.csv';
file23='../../nofeedback/pca_outputneurons/output/csv/dim_p=0.05_NY_STDP.csv';
file24='../../nofeedback/pca_outputneurons/output/csv/dim_p=0.05_Yoked_STDP.csv';

file211='../../onfeedback/pca_outputneurons/output/csv/dim_p=0.05_NY_NSTDP.csv';
file212='../../onfeedback/pca_outputneurons/output/csv/dim_p=0.05_Yoked_NSTDP.csv';
file213='../../onfeedback/pca_outputneurons/output/csv/dim_p=0.05_NY_STDP.csv';
file214='../../onfeedback/pca_outputneurons/output/csv/dim_p=0.05_Yoked_STDP.csv';
%%%%%%%%

%%%%%%%% read scv files
Ex.nofeedback_Nomal_NSTDP = csvread(file1);
Ex.nofeedback_Scramble_NSTDP = csvread(file2);
Ex.nofeedback_Nomal_STDP = csvread(file3);
Ex.nofeedback_Scramble_STDP = csvread(file4);

Ex.onfeedback_Nomal_NSTDP = csvread(file11);
Ex.onfeedback_Scramble_NSTDP = csvread(file12);
Ex.onfeedback_Nomal_STDP = csvread(file13);
Ex.onfeedback_Scramble_STDP = csvread(file14);

Ou.nofeedback_Nomal_NSTDP = csvread(file21);
Ou.nofeedback_Scramble_NSTDP = csvread(file22);
Ou.nofeedback_Nomal_STDP = csvread(file23);
Ou.nofeedback_Scramble_STDP = csvread(file24);

Ou.onfeedback_Nomal_NSTDP = csvread(file211);
Ou.onfeedback_Scramble_NSTDP = csvread(file212);
Ou.onfeedback_Nomal_STDP = csvread(file213);
Ou.onfeedback_Scramble_STDP = csvread(file214);

%%%%%%%%


%%%%%%%% arrange date
%nofeedback_Nomal_NSTDP= {'stdp';num2str(nofeedback_Nomal_NSTDP(1,1));num2str(nofeedback_Nomal_NSTDP(2,1)'.')};
%nofeedback_Nomal_NSTDP = mat2dataset(nofeedback_Nomal_NSTDP,'VarNames',{'Nofeedback_Nomal_NSTDP'})
Nomal_NSTDP = [Ex.nofeedback_Nomal_NSTDP(1,1);Ex.onfeedback_Nomal_NSTDP(1,1);Ou.nofeedback_Nomal_NSTDP(1,1);Ou.onfeedback_Nomal_NSTDP(1,1)];
Scramble_NSTDP = [Ex.nofeedback_Scramble_NSTDP(1,1);Ex.onfeedback_Scramble_NSTDP(1,1);Ou.nofeedback_Scramble_NSTDP(1,1);Ou.onfeedback_Scramble_NSTDP(1,1)];
Nomal_STDP = [Ex.nofeedback_Nomal_STDP(1,1);Ex.onfeedback_Nomal_STDP(1,1);Ou.nofeedback_Nomal_STDP(1,1);Ou.onfeedback_Nomal_STDP(1,1)];
Scramble_STDP = [Ex.nofeedback_Scramble_STDP(1,1);Ex.onfeedback_Scramble_STDP(1,1);Ou.nofeedback_Scramble_STDP(1,1);Ou.onfeedback_Scramble_STDP(1,1)];
%%%%%%%%

Name = {'All.Ex.Nofeedback';'All.Ex.Onfeedback';'Out.Ex.Nofeedback';'Out.Ex.Onfeedback';};
T = table(Name,Nomal_NSTDP,Scramble_NSTDP,Nomal_STDP,Scramble_STDP);
writetable(T,'./result/contributed_dimentions.csv');



nofeedback_Scramble_NSTDP = csvread(file2);
nofeedback_Nomal_STDP = csvread(file3);
nofeedback_Scramble_STDP = csvread(file4);
