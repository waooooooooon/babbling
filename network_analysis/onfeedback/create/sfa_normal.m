function sfa_normal(id)
% see Wiskott, L. and Sejnowski, T.J. (2002), "Slow Feature Analysis:
% Unsupervised Learning of Invariances", Neural Computation, 14(4):715-770,
% Figure 2

global tag simutime createddata_dir id_dir outdir firingdir pca_dir onfeedbackdir

    firingsdata=importdata([firingdir,'/firing_onfeedback_',id,'_',num2str(simutime),'.txt']);
    motcommanddata=importdata([firingdir,'/motorcommand_onfeedback_',id,'_',num2str(simutime),'.txt']);
    outputdir = [onfeedbackdir,'/through_simulation/sfa'];
    if ~exist(outputdir, 'dir')
        mkdir(outputdir);
    else
        addpath(outputdir);
    end
  
    
msg = ['This demo reproduces an example from Wiskott, L. and Sejnowski,' ...
       ' T.J. (2002), "Slow Feature Analysis: Unsupervised Learning of' ...
       ' Invariances", Neural Computation, 14(4):715-770, Figure 2\n\n'];

fprintf(msg);


% init graphics
figure; clf; set(gcf, 'Position', [89 477 866 338]);

% create the input signal
T = 5000;
t = linspace(0, 1, T);

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


Firings=firings.'; %invert to caliculate



% plot the input signal

plot(firingsdata(:,2),firingsdata(:,3),'.'); % Plot all the neurons'' spikes
title('Reservoir Firings', 'fontweight','bold');
axis([0 simutime 0 1000]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Slow Feature Analysis
%
% create a SFA object
hdl = sfa2_create(2, xp_dim(2), 'PCA');
% perform the preprocessing step
sfa_step(hdl, Firings, 'preprocessing');
% perform the expansion step
sfa_step(hdl, Firings, 'expansion');
% close the algorithm
sfa_step(hdl, [], 'sfa');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% compute the output signal
y = sfa_execute(hdl, Firings);

% it would have been quicker (but less instructive) to write:
% [y, hdl] = sfa2(x);

subplot(4,1,1);
% plot the output of the slowest varying function
fig121 = plot(y(:,1));
xlim([0 simutime])
ylim([-5 5])
%set(gca, 'Xlim', [0,2*pi], 'PlotBoxAspectRatio', [1,1,1])
title('output of the slowest varying function');
xlabel('t'); ylabel('y_1(t)');
%saveas(fig101,[outputdir,'/sfa_',id,'_',num2str(simutime),'.png']);
subplot(4,1,2);
fig131 = plot(y(:,2));
xlim([0 simutime])
ylim([-5 5])
%set(gca, 'Xlim', [0,2*pi], 'PlotBoxAspectRatio', [1,1,1])
title('output of the 2nd slowest varying function');
xlabel('t'); ylabel('y_2(t)');
subplot(4,1,3);
fig141 = plot(y(:,3));
xlim([0 simutime])
ylim([-10 10])
%set(gca, 'Xlim', [0,2*pi], 'PlotBoxAspectRatio', [1,1,1])
title('output of the 3nd slowest varying function');
xlabel('t'); ylabel('y_2(t)');
subplot(4,1,4);
fig237 = plot(motcommanddata(:,1),motcommanddata(:,2));
saveas(fig121,[outputdir,'/sfa_',id,'_',num2str(simutime),'.png']);
close all;


% plot the contours of the slowest varying function
%{
n = 100;
[X1,X2] = meshgrid(linspace(-1.2,2.2,n),linspace(-1.2,1.2,n));
testpts = [X1(:) X2(:)];

out = sfa_execute(hdl, testpts);
out = reshape(out(:,1),n,n);

subplot(1,3,3); contourf(X2,X1,out,10);
title('contours of the slowest varying function g_1(x)');
xlabel('x_2'); ylabel('x_1');
%}