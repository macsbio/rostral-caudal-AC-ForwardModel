% This script simulates BOLD responses from Neuronal Responses generated
% via models with varying spectral and temporal processing in the auditory
% belt. The script also performs comparison with measured responses using
% VB-GLM to predict the best neuronal model for the measured response.
% 
% Uses fMRI data from Santoro et al., 2017 (Experiment 2)
% =======================
% Author: Martin Havlicek
% =======================

clear all; close all;
set(0,'DefaultAxesFontSize', 14, ...
    'DefaultAxesFontAngle','normal', ... % Not sure the difference here
    'DefaultAxesFontWeight','normal', ... % Not sure the difference here
    'defaultLineLineWidth', 2, ...
    'defaultLineMarkerSize',5,...
    'DefaultAxesFontSize', 14, ...
    'DefaultAxesTitleFontWeight','normal', ...
    'DefaultAxesTitleFontSizeMultiplier', 1);


dir = cd;

% the neuronal model space consisted of 28 models sepcifying temporal
% dynamics (using tau) and spectral specificity (Q) of simulated belt
tau = [ 3 1; 3 1; 3 1; 3 1;  20 18; 20 18; 20 18; 20 18; 50 48; 50 48; 50 48; 50 48; 100 70; 100 70; 100 70; 100 70; 200 170; 200 170; 200 170; 200 170; 300 270; 300 270; 300 270; 300 270; 400 370;400 370; 400 370; 400 370];
Qs = repmat(fliplr([3 6 9 12]),1,7);

tic
runs = 12; % from the fMRI experiment
models = length(tau);
channels = 10; % the downsampled spectral axis in the neuronal simulations
subject = 'subj1';
hemisphere = 'RH';

colors = colormap(jet(models));
load([pwd '\stimOrder\' subject '_stimOrder.mat']); % this mat file contains order in which sound stimuli were presented during the fMRI experiment

neu_resp_simulated_area=  zeros(40,10,12);
% Get deconvolved (average) neuronal responses
% using the same technique as for deconvolving HRF function

for m = 1:models
    disp(['****Model ',num2str(m),' ****']);
  
    for i = 1:runs
        % load fast  Updated the path here!!!
        load([dir,'\neuronal simulations\' subject '\Q_', num2str(Qs(m)),'_tau_',num2str(tau(m,1)), '_' num2str(tau(m,2)), '\model_out_run' num2str(i) '.mat']);
        
        Xon   = zeros(size(simulated_area,2),40);
        S{i} = simulated_area';
        
        for k = 1:40,
            Xon(find(eval(['stimOrder.run',num2str(i),'>0']))*26 - 20+k-1,k) = 1;
        end
        XonB   = zeros(size(simulated_area,2),200);
        for k = 1:200,
            XonB(find(eval(['stimOrder.run',num2str(i),'>0']))*26 - 20+k-1,k)  = 1;
        end
        check = resample(simulated_area(1,:)',1,26);
        XonB_TR   = zeros(size(check,1),8);
        for k = 1:8,
            XonB_TR(find(eval(['stimOrder.run',num2str(i),'>0'])) - 1+k-1,k)  = 1;
        end
        for ch = 1:10
            neu_resp_simulated_area(:,ch,i,m)  = Xon(1:size(simulated_area',1),:)\simulated_area(ch,:)';
            
            [BOLD_tc_simulated_area{i}(:,ch) Neu_tc_simulated_area{i}(:,ch) Flow_tc_simulated_area{i}(:,ch)] = hemodynamic_model(simulated_area(ch,:)',0.1,2.5);
            BOLD_tc_simulated_area_TR{i}(:,ch) = interp1([1:1:size(BOLD_tc_simulated_area{i}(:,ch),1)],BOLD_tc_simulated_area{i}(:,ch),[6:26:size(BOLD_tc_simulated_area{i}(:,ch),1)],'linear');
            Neu_tc_simulated_area_TR{i}(:,ch) = interp1([1:1:size(Neu_tc_simulated_area{i}(:,ch),1)],Neu_tc_simulated_area{i}(:,ch),[12:26:size(Neu_tc_simulated_area{i}(:,ch),1)],'linear');
            
            
            bold_resp_simulated_area(:,ch,i) = XonB(1:size(simulated_area',1),:)\BOLD_tc_simulated_area{i}(1:size(simulated_area',1),ch);
            Neu_resp_simulated_area(:,ch,i) = XonB(1:size(simulated_area',1),:)\Neu_tc_simulated_area{i}(1:size(simulated_area',1),ch);
            Flow_resp_simulated_area(:,ch,i) = XonB(1:size(simulated_area',1),:)\Flow_tc_simulated_area{i}(1:size(simulated_area',1),ch);
            bold_TR_simulated_area(:,ch,i) = XonB_TR(1:size(check,1),:)\BOLD_tc_simulated_area_TR{i}(1:size(check,1),ch);
%             
%             figure(ch); subplot(1,2,1);
%             plot(neu_resp_simulated_area(:,ch,i,m),'color','b','linewidth',5);  x = gca; x.XLim = [-2 45]; x.YLim = [-1 8];
%             
%             subplot(1,2,2);
%             plot(bold_TR_simulated_area(:,ch,i),'color','r','linewidth',5); x = gca; x.XLim = [0 9]; x.YLim = [-1 1];
%             
        end
        
    end
    
    close all;
    
    %load stimOrder;
    % Specific for subjects!!!
    load([pwd '\masked timecourses\' subject '_masked_timeCourseData_' hemisphere]);
    
    % build a GLM model
    region = {'fast','slow'};
    for reg = 1:length(region)
        ROI = eval(region{reg});
        TR = 2.6;
        cut_off = 128;
        model1 = eval(['BOLD_tc_simulated_area_TR']);

        d = 2;
        
        Xx1 = [];
        Xc1 = [];
        
        Vx1 = [];
        Xx1 = [];
        Y  = [];
        Ym = [];
        CM = [];
        for run = 1:12,
            
            Xr1 = model1{run}(:,:);
            
            
            % do it also each run separatetly.... for comparison...
            Xx1 = [Xx1;Xr1];
            Xc1 = blkdiag(Xc1,ones(size(Xr1,1),1));

            xY  = pca_clean(eval(['ROI.tc.run',num2str(run)]),TR,cut_off);
            
            
            
            Y = [Y;double(xY.y(1+d:size(Xr1,1)+d,:))];
            Ym = [Ym;double(xY.y(1+d:size(Xr1,1)+d,:))-repmat(mean(double(xY.y(1+d:size(Xr1,1)+d,:)),1),size(Xr1,1),1)];
            % calculate the mean HRF
            
            block = 7;
            stim = eval(['stimOrder.run',num2str(run)]);
            ConvMat = zeros(size(xY.y,1),block);
            pos = find(stim>0);
            for k = 0:block
                ConvMat(pos+k-1,k+1) = 1;
            end
            ConvMat = ConvMat(1:size(Xr1,1),:);
            CM      = [CM;ConvMat]; % if not z-scored
            
            
        end
        m1 = repmat(mean(Xx1,1),size(Xx1,1),1);
        Xx1 = Xx1 - m1;
               
        % test the effect of standardization!!!!
        Xx1 = Xx1./repmat(std(Xx1,[],1),size(Xx1,1),1);
         
        
        % model fitting
        
        % PCA:
        C1 = Xx1'*Xx1;
                
        [V1, D1] = eig(C1);
        
        [D1,sind1] =sort(diag(D1),'descend');
        
        V1 = V1(:,sind1);
        
        % using first three eigenvectors
        V1r = V1(:,1:3);
        PC1 = Xx1*V1r;
        F1  = PC1'*Xx1;        
        X1 = [Xx1,Xc1];      
        
        for vox = 1:size(Y,2)
                      
            temp_hrf = [CM,Xc1]\Y(:,vox);
            HRF_est_data{reg}(:,vox) = temp_hrf(1:8);

            %
            % create basis and transformation matrix
            [betasPC{reg,1}(:,vox,m),~,FEpc{reg,1}(vox,m)] = vb_glm(Y(:,vox),[PC1,Xc1],0);
            Recon1 = PC1*betasPC{reg,1}(1:3,vox,m);
            F1_test = PC1'*Recon1;
            
            
            for i = 1:channels
                [~, ~, F_channel{reg,1}(i,vox,m)] = vb_glm(F1_test,F1(:,i),0);
                
            end
            
            temp_Fch1 = F_channel{reg,1}(:,vox,m) - squeeze(min(F_channel{reg,1}(:,vox,m),[],1));
            p_channelF{reg,1}(:,vox,m) = exp(temp_Fch1)./sum(exp(temp_Fch1));
                       
            
            HRF_predPC{reg,1}(:,vox,m) = [CM]\(PC1(:,1:3)*betasPC{reg,1}(1:3,vox,m));            
            
            res1pc = Y(:,vox)-[PC1,Xc1]*betasPC{reg,1}(:,vox,m);
            SSE1pc = res1pc'*res1pc;
            SSR1pc = (PC1*betasPC{reg,1}(1:3,vox,m))'*(PC1*betasPC{reg,1}(1:3,vox,m));
            SST  =  Ym(:,vox)'*Ym(:,vox);
            r2pc{reg,1}(vox,m)  = 1 - SSE1pc./SST;
            
            
            z = zeros(1,3);
            for ch =1:3
                z(ch) = 1;
                conX = z;
                conI = zeros(1,12);
                con = [conX,conI];
                tsPC{reg,1}(vox,ch,m) = con*betasPC{reg,1}(:,vox,m)./sqrt(SSE1pc./(size(Y,1)-(size([PC1,Xc1],2)+1))*con*inv([PC1,Xc1]'*[PC1,Xc1])*con');
            end
            % add F test!!!
            C = [eye(3);zeros(12,3)];
            X1pc = [PC1,Xc1];
            
            C0   = eye(size(X1pc,2))-C*pinv(C);
            X10  = X1pc*C0;
            R1   = eye(size(X1pc,1)) - X1pc*pinv(X1pc);
            R10  = eye(size(X1pc,1)) - X10*pinv(X10);
            
            M1  = R10 - R1;
            dfr1 = rank(X1pc)-rank(X10);
            dfe1 = size(X1pc,1)-rank(X1pc);
            
            F_value{reg,1}(vox,m) = (betasPC{reg,1}(:,vox,m)'*X1pc'*M1*X1pc*betasPC{reg,1}(:,vox,m))./(Y(:,vox)'*R1*Y(:,vox)).*dfe1./dfr1;
            
            p_value{reg,1}(vox,m) = 1 - cdf('F',F_value{reg,1}(vox,m),dfr1,dfe1);
            
        end
    end
end
toc


% IF you can run it and save all the results here
save([subject '_results_28models_' hemisphere '.mat']);

