
%%%%%%%%%%importdata
i=916;
id = ['synth_1_170222_1000_reinforce_100_4_NY_1_1_0.05_0.5_STDP_',num2str(i)];
dataname = ['p=0.05_stdp_wave/',id,'.wav'];
[data,Fs]=audioread(dataname);
%%%%%%%%%%

%%%% create directry
feedbackdir = [id, '_forgif'];
mkdir(feedbackdir);
%%%%

%%%%%% initialization

fftsize=2048;             %fft????????????????????????????????????
lpcsize=8;                %LPC size
 frameSize = 0.05;               % �t���[�����F0.025�b�i25ms�j
 frameShift = 0.010;              % �t���[���V�t�g���F0.010�b�i10ms�j


 %data = filter([1 -0.97],1,data);     %�v���G���t�@�V�X

 frameSizeSample = fix( Fs * frameSize );        % �t���[�����F�T���v�����Z
 frameShiftSample = fix( Fs * frameShift );      % �t���[���V�t�g���F�T���v�����Z
 disp( strcat('�T�E���h�f�[�^�̃T���v���� ',int2str(length(data)),'[sample], ',...
     '�t���[���� ',int2str(frameSizeSample),'[sample], ',...
     '�t���[���V�t�g�� ',int2str(frameShiftSample),'[sample]') );

 %�t���[���̑��������߂�
 maxFrame = ...
  fix( ( length(data) - (frameSizeSample - frameShiftSample) ) / frameShiftSample ) - 1;

 %�t���[���̌������������J��Ԃ�
 startThisFrame = 1;                                  % �t���[���̊J�n�T���v���ԍ�
 endThisFrame = startThisFrame + frameSizeSample - 1; % �t���[���̏I���T���v���ԍ�
 for countFrame = 1 : 1 : maxFrame
     disp(strcat('�t���[���ԍ� ',int2str(countFrame),' / ',int2str(maxFrame)));
     disp(strcat('���̃t���[���̊J�n�T���v���ԍ� ',num2str(startThisFrame)));
     disp(strcat('���̃t���[���̏I���T���v���ԍ� ',num2str(endThisFrame)));
     thisData = data(startThisFrame : endThisFrame);

     % �����ł��̃t���[���ɑ΂���v�Z���s��
%     subplot(10, 9, countFrame); plot(thisData); ylim([-3,3]);
%    axis off; %�ڐ��������
    
            %%%%%%%% create png file of spectrum of created sound
        center = fix(length(thisData)/2);
        cuttime=0.04; %
        wavdata = thisData(center-fix(cuttime/2*Fs) : center+fix(cuttime/2*Fs));%
        %
        han_window = 0.5 - 0.5 * cos(2 * pi * [0 : 1/length(wavdata) : 1]);%
        wavdata = han_window(1:length(wavdata))' .* wavdata;%'

        [P,f]=pyulear(wavdata,lpcsize,fftsize,Fs); %lpc
        AP=abs(P)/fftsize;
        fscale=linspace(0,Fs,fftsize);  %
        

        
        figure(1);
        %plot(fscale(1:fftsize/4),AP(1:fftsize/4),'k','LineWidth',2);ylim([0,0.000000003]);xlim([0,2000]) %[Hz]
        set(gca,'FontSize',20);
        
        %feedback=13*AP(1:2:((fftsize-48)/4)-300)/max(AP(1:2:((fftsize-48)/4)-300)); %normalization 1~2KHz
        feedback=13*AP(51:1:((fftsize-48)/4)-400+50); %does not normalization  % 500~1.5KHz
        
        fig146 = figure('visible', 'off');
        %auditory feedback
        
        fig146 = plot(feedback,'sk','LineWidth',2);
        set(gca,'FontSize',20);
        ylim([0 1*10^-7]);
        
        %yoked feedback
        %randid=randperm(100);
        %yokedfeedback=feedback(randid);
        
        %fig146 = plot(yokedfeedback,'sk','LineWidth',2);
        %set(gca,'FontSize',20);
        
        saveas(fig146,[feedbackdir,'/',num2str(countFrame),'.png']);
        
        close all;
        %%%%%%%% end create png file of spectrum of created sound
        
        

     startThisFrame = startThisFrame + frameShiftSample;
    endThisFrame = startThisFrame + frameSizeSample - 1;
 end