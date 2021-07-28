%% Simulating neuronal responses of a core and a belt auditory region
% For details on the model, see Zulfiqar et al. 2020
% This script simulates neuronal responses with varying
% spectro-temporal properties in the Belt region. This is acheieved by
% manipulating the: 
%    Time constant, tau (for temporal dynamics)
%    The connectivity kernel (controlling connectivity between the simulated core and belt region), and
%    Sigmas (spatial spread of activation)
% Overall, these parameters control the dynamics of the belt region.

% Input signal is passed throught the model to mimic the stimulus presented
% in the fMRI experiment (Santoro et al., 2017) in a fast event related
% design where stimulus is presented in silent interval after acquisition.

% Uses gammatone filterbank implementation by Ma et al. (2007) available at
% http://staffwww.dcs.shef.ac.uk/people/N.Ma/resources/gammatone/

% =====================
% Author: Isma Zulfiqar
% =====================

%%
clear all
close all
clc

models = 1; % multiple models can also be simulated

% connectivity kernel for the model; connectivity between core and belt
connectivity = [0.5 1    1  1  1  1  1  1   0.5];

% varying parameters of the simulated Belt
taus = [3;1]; % time constant for the model, reduces along tonotopic axis
all_sigmas = [200; 300];% sigma controlling the spread of the activity                
all_Qs = 3; % Q signifies the quality factor, i.e., the measure of sharpness of tuning. Low value indicates broad tuning.

for model = 1:length(models)
    
    sigmaEE4 = all_sigmas(1,model);
    sigmaEI4 = all_sigmas(2,model);
    
    DT_high = taus(1,1);
    DT_low = taus(2,1);
    
    Q = all_Qs(model);
    
    % =========================== setting up model
    
    Size = 100;  %Spatial (Frequency) size of array - number of units
    
    % change these for different length stimuli
    lim1_frqaxis = 50;
    lim2_frqaxis = 8000;
    F = MakeErbCFs(lim1_frqaxis,lim2_frqaxis,Size);  % ERB freq axis
    X = F;
    
    % the TR and TA are specific to the fMRI experiment
    TR = 2.6; % s,
    TA = 1.2; % s
    jitter = 0;
    
    nSounds = 1;
    % generating a sample sound (1kHz tone) to be processed by the model
    Fs = 16000; % sampling rate
    f = 1000; % Hz
    stim = sin(2*pi*f*(1/Fs:1/Fs:1));
    
    for cursound = 1:nSounds
        
        s = zeros(1, TR*Fs);
        s(1600:1600+Fs-1) = stim(cursound,:);
        
        tic
        %% Peripheral Processing Stage
        % passing sound through gammatone filterbank
        bm_upsampled = gammatoneFast(s',F,Fs);
        
        y2_h = bm_upsampled(:,Size);
        y3_h = 0; y4 = []; y5 =[]; v5 = [];
        for ch = Size-1:-1:1
            y2 = bm_upsampled(:,ch);
            
            % ===== Lateral inhibition (masked by higher
            % (frequency) spatial response)
            
            y3 = y2 - y2_h;
            y2_h = y2;
            
            % ==== half wave rectification --- > y4
            y4(:,ch) = max(y3,0);
            
            
            % temporal integration window ---> y5
            y5 = filter(1, [1 -0.9845],y4(:,ch));
            v5(:,ch) = decimate(y5,64);
            
            
        end
        F_new = F(2:Size-1);
        bm = y4(:,2:Size-1);
              
        Last = floor(1000 * (length(s)/Fs));  %last time in computation in ms
        
        DelT = 1000/Fs; %0.05; % dt - time resolution of the simulations in ms
        
        time = 0:DelT:Last;
        step = DelT/1000; % ms
        timei = 0:step:Last/1000;  % time axis with DelT for simulations (in s.)
        FS = 1/step;
        bm = bm(ceil(linspace(1,length(s),length(timei))),:);
        bm = 1E3*bm;
        clear bm_upsampled;
        
        
        Stim = bm';
        color='k'; figureSpec(timei,F_new,Stim,1,color)
        
        %% Cortical Processing Stage
        
        EE1 = zeros(Size-2, length(time)); IN1 = EE1; EEresp1 = EE1; INresp1 = EE1;
        EE4 = zeros(Size-2, length(time)); IN4 = EE4; EEresp4 = EE4; INresp4 = EE4;
        
        % Stimulus
        P = Stim(:,1:length(time)); Q_in=0; % input to the IN population is set to 0
        DX = 20; 
        Xsyn1 = DX*(-5:5);Xsyn2 = DX*(-5:5); Xsyn3 = DX*(-2:2); % filter width
        Xsyn4 = DX*(-2:2);
        
        % naka rushton constants
        
        % A1 (core)
        thExc1 = 60;
        thIn1 = 80;
        max1 = 100;
        
        % belt
        thExc4 = 60;%20;
        thIn4 = 80;%40;
        max4 = 100;%12.5;
        
        % maximum synaptic strength for active transient response
        EEgain = 1.5;
        EIgain = 1.3;
        IIgain = 1.5;
        IEgain = EIgain;
        
        % tau - time constant of the network in ms;       
        DT1 = 10; % for A1
        DT4 = linspace(DT_high,DT_low,Size-2)'; % for Belt region
        
        % spatial spread
        % A1
        sigmaEE1 = 40;%20;%40;
        sigmaEI1 = 160;%80;%160;
        sigmaIE1 = sigmaEI1;
        sigmaII1 = 10;
        disp('For A1:'); checkSigmas(sigmaEE1, sigmaEI1, sigmaIE1, sigmaII1, thExc1, thIn1, EEgain, EIgain, IEgain, IIgain , max1);
        
        % simulated belt area        
        sigmaIE4 = sigmaEI4;
        sigmaII4 = 30;
        disp('For simulated area:'); checkSigmas(sigmaEE4, sigmaEI4, sigmaIE4, sigmaII4, thExc4, thIn4, EEgain, EIgain, IEgain, IIgain , max4);
             
        % synaptic weights
        % A1
        synEE1 = EEgain*exp(-abs(Xsyn1)./sigmaEE1);
        synEI1 = EIgain*exp(-abs(Xsyn1)./sigmaEI1);
        synII1 = IIgain*exp(-abs(Xsyn1)./sigmaII1);
        % Belt
        synEE4 = EEgain*exp(-abs(Xsyn4)./sigmaEE4);
        synEI4 = EIgain*exp(-abs(Xsyn4)./sigmaEI4);
        synII4 = IIgain*exp(-abs(Xsyn4)./sigmaII4);
        
        % % % ===================== Simulating responses
        
        % ================= A1
        smoothing_kernel = [0.5 1 0.5];
        
        for T = 2:length(time)  %Loop in ms, Euler solution method
            
            input1 = conv(P(:,T),smoothing_kernel,'same')./sum(smoothing_kernel);%(((100-CC_inRatio)/100) * P(:,T)) + ((CC_inRatio/100) * CC_in(:,T));
            
            EEresp1(:,T) = NeuralConv(synEE1, EE1(:,T-1)) - NeuralConv(synEI1, IN1(:,T-1)) + input1;
            EEresp1(:,T) = (EEresp1(:,T).*(EEresp1(:,T) > 0)).^2;
            INresp1(:,T) = NeuralConv(synEI1, EE1(:,T-1)) - NeuralConv(synII1, IN1(:,T-1)) + Q_in;
            INresp1(:,T) = (INresp1(:,T).*(INresp1(:,T) > 0)).^2;
            
            EE1(:,T) = EE1(:,T-1) + (DelT/DT1)*(-EE1(:,T-1) + (max1)*EEresp1(:,T)./(thExc1^2 + EEresp1(:,T)));
            IN1(:,T) = IN1(:,T-1) + (DelT/DT1)*(-IN1(:,T-1) + (max1)*INresp1(:,T)./(thIn1^2 + INresp1(:,T)));
        end
        clear EEresp1 INresp1;
        
        figure(2); 
        color='b'; figureSpec(timei,F_new,EE1,2,color); title('A1');
        
        % =============== Belt       
        
        smoothing_kernel_new = connectivity(model,:);
        
        for T = 2:length(time)
            
            % adding a smoothing kernel (simple 1D)
            input4 = conv(EE1(:,T),smoothing_kernel_new,'same')./sum(smoothing_kernel_new);
            
            EEresp4(:,T) = NeuralConv(synEE4, EE4(:,T-1)) - NeuralConv(synEI4, IN4(:,T-1)) + input4;
            EEresp4(:,T) = (EEresp4(:,T).*(EEresp4(:,T) > 0)).^2;
            INresp4(:,T) = NeuralConv(synEI4, EE4(:,T-1)) - NeuralConv(synII4, IN4(:,T-1)) + Q_in;
            INresp4(:,T) = (INresp4(:,T).*(INresp4(:,T) > 0)).^2;
            
            EE4(:,T) = EE4(:,T-1) + (DelT./DT4).*(-EE4(:,T-1) + (max4)*EEresp4(:,T)./(thExc4^2 + EEresp4(:,T)));
            IN4(:,T) = IN4(:,T-1) + (DelT./DT4).*(-IN4(:,T-1) + (max4)*INresp4(:,T)./(thIn4^2 + INresp4(:,T)));
            
            
        end
        clear EEresp4 INresp4;
        
        figure(3);
        color='b'; figureSpec(timei,F_new,EE4,3,color); title('Belt');
        
        % ======== downsampling the data
        index = 1; n = 10; % space
        temp1 = []; temp4 = [];
        for i = 1:n:Size-n % downsampling over space
            temp1(index,:) = mean(EE1(i:i+n-1,:));
            temp4(index,:) = mean(EE4(i:i+n-1,:));
            
            index = index+1;
        end
        temp1(index,:) =  mean(EE1(91:end,:));
        temp4(index,:) =  mean(EE4(91:end,:));
        
        n = 1600; % time
        index = 1;
        for i = 1:n:length(temp1)-n % downsampling over time
            simulated_area(cursound,:, index) = mean(temp4(:,i:i+n-1),2);
            
            index = index+1;
        end
    end
    toc
    disp('Finished');
end
