
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
 frameSize = 0.05;               % フレーム長：0.025秒（25ms）
 frameShift = 0.010;              % フレームシフト長：0.010秒（10ms）


 %data = filter([1 -0.97],1,data);     %プリエンファシス

 frameSizeSample = fix( Fs * frameSize );        % フレーム長：サンプル換算
 frameShiftSample = fix( Fs * frameShift );      % フレームシフト長：サンプル換算
 disp( strcat('サウンドデータのサンプル数 ',int2str(length(data)),'[sample], ',...
     'フレーム長 ',int2str(frameSizeSample),'[sample], ',...
     'フレームシフト長 ',int2str(frameShiftSample),'[sample]') );

 %フレームの総数を求める
 maxFrame = ...
  fix( ( length(data) - (frameSizeSample - frameShiftSample) ) / frameShiftSample ) - 1;

 %フレームの個数だけ処理を繰り返す
 startThisFrame = 1;                                  % フレームの開始サンプル番号
 endThisFrame = startThisFrame + frameSizeSample - 1; % フレームの終了サンプル番号
 for countFrame = 1 : 1 : maxFrame
     disp(strcat('フレーム番号 ',int2str(countFrame),' / ',int2str(maxFrame)));
     disp(strcat('このフレームの開始サンプル番号 ',num2str(startThisFrame)));
     disp(strcat('このフレームの終了サンプル番号 ',num2str(endThisFrame)));
     thisData = data(startThisFrame : endThisFrame);

     % ここでこのフレームに対する計算を行う
%     subplot(10, 9, countFrame); plot(thisData); ylim([-3,3]);
%    axis off; %目盛りを消す
    
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