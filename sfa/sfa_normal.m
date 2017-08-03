function [y]= sfa_normal(id)
% see Wiskott, L. and Sejnowski, T.J. (2002), "Slow Feature Analysis:
% Unsupervised Learning of Invariances", Neural Computation, 14(4):715-770,
% Figure 2

msg = ['This demo reproduces an example from Wiskott, L. and Sejnowski,' ...
       ' T.J. (2002), "Slow Feature Analysis: Unsupervised Learning of' ...
       ' Invariances", Neural Computation, 14(4):715-770, Figure 2\n\n'];

fprintf(msg);

global k yoke STDP outdir p dim 
mkdir([outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/output_of_the_slowest_varying_function']);


% init graphics
figure; clf; set(gcf, 'Position', [89 477 866 338]);

% create the input signal
T = 5000;
t = linspace(0, 1, T);

%x1 = sin(t)+cos(11*t).^2;
%x2 = cos(11*t);
%x = [x1; x2]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%firing input
%%%%%%%%%%%%%%%%%%%%%%%%%%%
firings=importdata(id);
Firings=zeros(1000,5000);
for i=1:5000
  
    I=firings(find(firings(:,2)==i),3); %t=iに発火したニューロンidをIに
    C=size(I);  %Iのサイズ
    for j=1:C(1,1)
    Firings(I(j,1),i)=1;
    end
    
end

Firings=Firings.'; %invert to caliculate

% plot the input signal

subplot(1,3,1)
plot(firings(:,2),firings(:,3),'.'); % Plot all the neurons'' spikes
title('Reservoir Firings', 'fontweight','bold');
axis([0 5000 0 1000]);


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

% plot the output of the slowest varying function
subplot(1,3,2); 
fig101=plot(y(:,1));
xlim([0 5000])
%set(gca, 'Xlim', [0,2*pi], 'PlotBoxAspectRatio', [1,1,1])
title('output of the slowest varying function');
xlabel('t'); ylabel('y_1(t)');
saveas(fig101,[outdir,'/p=',num2str(p),'_',yoke,'_',STDP,'/output_of_the_slowest_varying_function/sfa_p=',num2str(p),'_',yoke,'_',STDP,'d=',num2str(k),'.png']);
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