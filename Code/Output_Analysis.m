% Generating Plots from the analysis. This script only plots the results
% using data from one single partipant and one hemishphere and is intended
% to show how the output from model fitting can be analyzed.

% Uses output of Model_fitting_comparison.m

% ========================================
% Authors: Isma Zulfiqar, Martin Havlicek
% ========================================

clear
clc
close all

subjects = {'subj1'};
hemispheres = {'RH'};

% model space
Qs = [12 9 6 3];
taus = [3 20 50 100 200 300 400];

count_model_slow = zeros(4,7);
count_model_fast = zeros(4,7);

counter_T_slow = []; counter_T_fast =[];
counter_Q_slow = []; counter_Q_fast = [];

for sub = 1:length(subjects)
    for hemi = 1:length(hemispheres)
        
        best_model = [];
        
        % this file is generated
        load([subjects{sub} '_results_28models_' hemispheres{hemi} '.mat']);
        
        % Free energy for the two labelled regions
        FE1  = [FEpc{1}'];
        dFE1 = FE1 - repmat(min(FE1,[],1),size(FE1,1),1);
        p10  =  exp(dFE1)./repmat(sum(exp(dFE1),1),size(FE1,1),1);
        
        FE2  = [FEpc{2}'];
        dFE2 = FE2 - repmat(min(FE2,[],1),size(FE2,1),1);
        p20  =  exp(dFE2)./repmat(sum(exp(dFE2),1),size(FE2,1),1);
        
        XYZ = [];
        Label = [];
        
        for reg = 1:2
            ROI     = eval(region{reg});
            XYZtemp = ROI.vertices_coordinates;
            XYZ     = [XYZ;XYZtemp];
            Label   = [Label;ones(size(XYZtemp,1),1)*reg];
        end;
        mXYZ = mean(XYZ,1);
        
        [V, D] = eig(cov(XYZ-repmat(mXYZ,size(XYZ,1),1)));
        
        e1=[mXYZ',mXYZ'+V(:,1)*D(1,1)]';
        e2=[mXYZ',mXYZ'+V(:,2)*D(1,1)]';
        e3=[mXYZ',mXYZ'+V(:,3)*D(1,1)]';
        
        v1 = V(:,1);
        v2 = V(:,2);
        v3 = V(:,3);
        XYZc =XYZ-repmat(mXYZ,size(XYZ,1),1);
        for i = 1:size(XYZ,1)
            pXYZ(i,:) = (v2'*XYZc(i,:)'/(v2'*v2)*v2 + v3'*XYZ(i,:)'/(v3'*v3)*v3);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        i1 = 0;
        i2 = 0;
        for i = 1:length(Label),
            
            if Label(i) == 1
                m_sel = find(p10(:,i)==max(p10(:,i)));
                i1 =  i1+1;
                
                b_channel1(i1,:) = p_channelF{1,1}(:,i1,m_sel);
            elseif Label(i) == 2
                m_sel = find(p20(:,i-sum(Label==1))==max(p20(:,i-sum(Label==1))));
                
                i2 =  i2+1;
                
                b_channel2(i2,:) = p_channelF{2,1}(:,i2,m_sel);
            end
            
            best_model(i,:) = m_sel;
            
            
        end;
        
        % removing any common vertices from the RoIs
        
        [~,~,ib] = intersect(a1.vertices, fast.vertices);
        [~,~,id] = intersect(r.vertices, slow.vertices);
        
        [~,ix,~] = intersect(fast.vertices, slow.vertices);
        
        id = id + length(dFE1);
        
        best_model(ib) = nan; % removing vertices of primary areas
        best_model(id) = nan;
        
        best_model(ix) = nan; % putting common vertices in slow area
        
        nvcs_fast = length(dFE1);
        nvcs_slow = length(dFE2);
        
        add_factor_fast = 1/(nvcs_fast-length(ib)-length(ix));
        add_factor_slow = 1/(nvcs_slow-length(id));
        
        % get the total number of voxels best represented by each model in
        % the model space for both Slow and Fast regions
        
        for i = 1:nvcs_fast
            if ~isnan(best_model(i))
                [x1, y1] = ind2sub([4 7],best_model(i));
                count_model_fast(x1, y1) = count_model_fast(x1, y1) + add_factor_fast;
            end
        end
        
        
        for i = nvcs_fast+1:length(best_model)
            if ~isnan(best_model(i))
                [x1, y1] = ind2sub([4 7],best_model(i));
                count_model_slow(x1, y1) = count_model_slow(x1, y1) + add_factor_slow;
            end
        end
        
        % get the tau (T) values from best fit model for each voxel in both labelled regions
        
        for i = 1:length(Label)
            if ~isnan(best_model(i))
                if Label(i) == 1
                    m_sel = find(p10(:,i)==max(p10(:,i)));
                    [x1,y1] = ind2sub([4 7],m_sel);
                    
                    counter_T_fast(sub,hemi,i,:) = taus(y1);
                elseif Label(i) == 2
                    m_sel = find(p20(:,i-sum(Label==1))==max(p20(:,i-sum(Label==1))));
                    [x1,y1] = ind2sub([4 7],m_sel);
                    
                    counter_T_slow(sub,hemi,i,:) = taus(y1);
                end
                
                
            end
        end
        
        % get the Q values from best fit model for each voxel in both labelled regions
        
        for i = 1:length(Label)
            if ~isnan(best_model(i))
                if Label(i) == 1
                    m_sel = find(p10(:,i)==max(p10(:,i)));
                    [x1,y1] = ind2sub([4 7],m_sel);
                    
                    counter_Q_fast(sub,hemi,i,:) = Qs(x1);
                elseif Label(i) == 2
                    m_sel = find(p20(:,i-sum(Label==1))==max(p20(:,i-sum(Label==1))));
                    [x1,y1] = ind2sub([4 7],m_sel);
                    
                    counter_Q_slow(sub,hemi,i,:) = Qs(x1);
                    
                end
            end
            
        end
    end
    
    disp(subjects{sub});
end

%% ============= Figure 6: Model Predictions

figure(6); subplot(1,4,1)
contourf(flipud(reshape(p10(:,4), [4 7]))); colormap(flipud(bone(12)));
x = gca; x.YLim = [1 4]; x.YTick = 1:4; x.YTickLabel = {'3','6','9','12'}; ylabel('Quality factor')
x.XLim = [1 7]; x.XTick = 1:7; x.XTickLabel = {'3','20','50','100', '200', '300', '400'}; xlabel('Time constant')
x.FontSize = 14; x.Box = 'off'; x.TickLength = [0.015 1]; colorbar
grid on; caxis([0 0.5]);
hold on;

figure(6); subplot(1,4,2)
contourf(flipud(reshape(p20(:,27), [4 7]))); colormap(flipud(bone(12)));
x = gca; x.YLim = [1 4]; x.YTick = 1:4; x.YTickLabel = {'3','6','9','12'}; ylabel('Quality factor')
x.XLim = [1 7]; x.XTick = 1:7; x.XTickLabel = {'3','20','50','100', '200', '300', '400'}; xlabel('Time constant')
x.FontSize = 14; caxis([0 0.5]);
x.FontSize = 14; x.Box = 'off'; x.TickLength = [0.015 1]; colorbar
grid on


figure(6); subplot(1,4,3)
xx = HRF_est_data{1, 1}(:,4);
plot(xx./xx(3),'b-','linewidth',1.5);
hold on
[~,yy] = max(p10(:,4));
xx = squeeze(HRF_predPC{1}(:,4,yy));
plot(xx./xx(3),'b:', 'LineWidth', 1.5);

xx = HRF_est_data{1, 2}(:,27);
plot(xx./xx(3),'r','linewidth',1.5);
[~,yy] = max(p20(:,27));
xx = squeeze(HRF_predPC{2}(:,27,yy));
plot(xx./xx(3),'r:', 'LineWidth', 1.5);

x = gca; x.YLim = [-0.5 1.5]; x.YTick = [-0.5:0.5:1.5];ylabel('Normalized BOLD response')
x.XLim = [0.1 8.2]; x.XTick = [0.38*3:0.38*5:8]; x.XTickLabel = {'1','5','10','15'}; xlabel('Time (s)');
x.TickLength = [0.015 1]; x.FontSize = 14; x.Box = 'off';
legend('Measured', 'Predicted')


figure(6); subplot(1,4,4)
[x1,y1] = max(b_channel1(4,:));
xx = mean(neu_resp_simulated_area(:,y1,:,11),3);
plot(xx./mean(xx(15:17)),'b-.','Linewidth',1.5);
hold on;
[x1,y1] = max(b_channel2(27,:));
xx = mean(neu_resp_simulated_area(:,y1,:,27),3);
plot(xx./mean(xx(15:17)),'r-.','Linewidth',1.5);
x = gca; x.YLim = [0 2.5]; x.YTick = 0:1:2; x.YTickLabel = {'0','1','2'};ylabel('Normalized Firing Rate')
x.XLim = [0 41]; x.XTick = [6 16 26 36]; x.XTickLabel = {'0','1','2', '3'}; xlabel('Time (s)');
x.FontSize = 14; x.Box = 'off'; x.TickLength = [0.015 1];

%% ========== Figure 9: Distribution of voxels best represented by each of the models in the model space
figure(9);
subplot(1,3,1)
contourf(flipud(count_model_slow*100/length(subjects)),'LineWidth', 2);
x = gca; x.YTick = [1:4]; x.XTick = [1:7];  x.YTickLabel = {'3','6','9','12'}; x.XTickLabel = {'3','20','50','100','200','300','400'}; x.TickLength = [0;0]; x.FontSize = 14;x.Box = 'off';
xlabel('Time constant'); ylabel('Quality Factor')
colorbar; caxis([0 20]);
subplot(1,3,2)
contourf(flipud(count_model_fast*100/length(subjects)),'LineWidth', 2);
x = gca; x.YTick = [1:4]; x.XTick = [1:7];  x.YTickLabel = {'3','6','9','12'}; x.XTickLabel = {'3','20','50','100','200','300','400'}; x.TickLength = [0;0]; x.FontSize = 14;x.Box = 'off';
xlabel('Time constant'); ylabel('Quality Factor')
colormap jet; colorbar; caxis([0 20])

subplot(1,3,3);
contourf(flipud((count_model_slow*100/length(subjects)) - (count_model_fast*100/length(subjects))),'LineWidth', 2);
x = gca; x.YTick = [1:4]; x.XTick = [1:7];  x.YTickLabel = {'3','6','9','12'}; x.XTickLabel = {'3','20','50','100','200','300','400'}; x.TickLength = [0;0]; x.FontSize = 14; x.Box = 'off';
xlabel('Time constant'); ylabel('Quality Factor')
colormap jet; colorbar;
