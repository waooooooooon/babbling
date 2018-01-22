function [f,P1]= fft_script(input,ploton)
global k yoke STDP outdir p dim simutime
%{
fs=1000; %sample frequency (Hz)
y=fft(input); %fft
n=length(input); %number of samples
f=(0:n-1)*(fs/n); %frequency range
power = abs(y).^2/n; %power of the DFT


y0 = fftshift(y);         % shift y values
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
power0 = abs(y0).^2/n;    % 0-centered power

if plot==1
plot(f0,power0);
xlabel('Frequency');
ylabel('Power');
end
%}

Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = simutime;             % Length of signal
time = (0:L-1)*T;        % Time vector

%{  
S = 0.7*sin(2*pi*50*time) + sin(2*pi*120*time);    %for debug deta
X = S + 2*randn(size(time));

plot(1000*time(1:50),X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')
%}

Y = fft(input);


P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

%{
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
xlim([0,1]);
ylabel('|P1(f)|')
%}