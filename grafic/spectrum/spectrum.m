function [salience] = spectrum2(id)
%analsys wav file
%
%
%  Example of Use:
%     spectrum('1_0.5')
%
% Authors: Tomohiro Takimoto
% Cognitive and Information Sciences
% University of Osaka, Merced
% email: tomohiro.takimoto@ams.eng.osaka-u.ac.jp
%
% data:180120Sctime/1_No_STDP
%low salience ->12
%high salience ->1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lpcsize=8;
fftsize=2048;
% Directory for Coath et. al. Saliency Detector.
addpath('auditorysaliencymodel');


wavname=['wavdata/synth_1_180120_Sctime_2000_reinforce_100_4_No_0_1_0.03_0.3_STDP_LiIP_randSc_random_normal_fft_1989.wav'];
        [data,Fs]=audioread(wavname); % dft=fft(data,fftsize); pdft=abs(dft).^2; %power spectrum


Fs = 44100; %ÉTÉìÉvÉäÉìÉOé¸îgêî

%subplot(2,1,1);
%plot(data);
%xlabel('éûä‘[sample]');
%xlim([1 length(data)]);
%ylabel('êUïù');
subplot(1,1,1);
spectrogram(data, hamming(64), 32, 256, Fs, 'yaxis');
xlabel('Time(s)');
ylabel('Frequency(kHz)');
set(gca,'FontSize',20);
saveas(gca,['./wavdata/1989.png']);

        %{
        center = fix(length(data)/2);
        cuttime=0.04; %
        wavdata = data(center-fix(cuttime/2*Fs) : center+fix(cuttime/2*Fs));%
        %
        han_window = 0.5 - 0.5 * cos(2 * pi * [0 : 1/length(wavdata) : 1]);%
        wavdata = han_window(1:length(wavdata))' .* wavdata;%'

        [P,f]=pyulear(wavdata,lpcsize,fftsize,Fs); %lpc
        AP=abs(P)/fftsize;
        fscale=linspace(0,Fs,fftsize);  %

        plot(fscale(1:fftsize/4),AP(1:fftsize/4));ylim([0,0.000000003]);xlim([0,5000]) %[Hz]
       %dft = fft(wavdata, fftsize);                     %
       %Adft = abs(dft) / fftsize;
%}
