
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

%%%%%%%%

%%%%%%%% read scv files
nofeedback_Nomal_NSTDP = csvread(file1);
nofeedback_Scramble_NSTDP = csvread(file2);
nofeedback_Nomal_STDP = csvread(file3);
nofeedback_Scramble_STDP = csvread(file4);

onfeedback_Nomal_NSTDP = csvread(file11);
onfeedback_Scramble_NSTDP = csvread(file12);
onfeedback_Nomal_STDP = csvread(file13);
onfeedback_Scramble_STDP = csvread(file14);
%%%%%%%%


%%%%%%%% arrange date
%nofeedback_Nomal_NSTDP= {'stdp';num2str(nofeedback_Nomal_NSTDP(1,1));num2str(nofeedback_Nomal_NSTDP(2,1)'.')};
%nofeedback_Nomal_NSTDP = mat2dataset(nofeedback_Nomal_NSTDP,'VarNames',{'Nofeedback_Nomal_NSTDP'})
Nomal_NSTDP = [nofeedback_Nomal_NSTDP(1,1);onfeedback_Nomal_NSTDP(1,1)];
Scramble_NSTDP = [nofeedback_Scramble_NSTDP(1,1);onfeedback_Scramble_NSTDP(1,1)];
Nomal_STDP = [nofeedback_Nomal_STDP(1,1);onfeedback_Nomal_STDP(1,1)];
Scramble_STDP = [nofeedback_Scramble_STDP(1,1);onfeedback_Scramble_STDP(1,1)];
%%%%%%%%

Name = {'Nofeedback';'Onfeedback';};
T = table(Name,Nomal_NSTDP,Scramble_NSTDP,Nomal_STDP,Scramble_STDP);
writetable(T,'./result/contributed_dimentions.csv');



nofeedback_Scramble_NSTDP = csvread(file2);
nofeedback_Nomal_STDP = csvread(file3);
nofeedback_Scramble_STDP = csvread(file4);
