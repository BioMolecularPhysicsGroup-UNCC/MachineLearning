% clean; rng(now);
% Do SPLOC between the 

initializeSPLOC("gtype","png")                          % Initialize SPLOC
% This line initializes the log file and global variables used in sploc
% functions!

                                                       % Get training data
                                                       
[~,~,FnameF_tr] = readFileNameList('trainingdata','sType','F');  % Tr Func
[~,~,FnameN_tr] = readFileNameList('trainingdata','sType','N'); % Tr NFunc

[~,~,FnameF_ts] = readFileNameList('testingdata','sType','F');   % Ts Func
[~,~,FnameN_ts] = readFileNameList('testingdata','sType','N');  % Ts NFunc
[~,~,FnameU_ts] = readFileNameList('testingdata','sType','U');  % Ts Uknow

% Input data
mFormat = setDataMatrixFormat('notVectored',1);     % Sets the data format

                                                        % Read in matrices
                                                                % Training
[AmatrixF_tr,tableF_tr] = readDataMatrices('functional',FnameF_tr,mFormat);
[AmatrixN_tr,tableN_tr] = readDataMatrices('nonfunctional',FnameN_tr,mFormat);
                                                                 % Testing
[AmatrixF_ts,tableF_ts] = readDataMatrices('functional',FnameF_ts,mFormat);
[AmatrixN_ts,tableN_ts] = readDataMatrices('nonfunctional',FnameN_ts,mFormat);
[AmatrixU_ts,tableU_ts] = readDataMatrices('undetermined',FnameU_ts,mFormat);
% This outputs a struct which hold the data and is accessible by other
% functions in the SPLOC toolkit

                                           % Calculate training statistics
CovF_tr = getMultivariateStats4sploc(AmatrixF_tr,samp_siz,1,'cov');
CovN_tr = getMultivariateStats4sploc(AmatrixN_tr,samp_siz,1,'cov');

CovF_ts = getMultivariateStats4sploc(AmatrixF_ts,samp_siz,1,'cov');
CovN_ts = getMultivariateStats4sploc(AmatrixN_ts,samp_siz,1,'cov');
CovU_ts = getMultivariateStats4sploc(AmatrixU_ts,samp_siz,1,'cov');
% This extracts means and covariance from all data packets for SPLOC to use
% in the learning process

% SPLOC it!
splocResults = sploc(0, 1, 'trial1', CovF_tr, CovN_tr, 0); % Perform SPLOC

plotCongruencySpectrum(splocResults,3);        % Show final SPLOC solution

U = getDiscriminantSBV(splocResults);   % Extract Disc Modes from solution
if(do_save == 1); dlmwrite('DiscriminateModes.txt', U); end

% {
% Classify it!
if(size(U,1) > 0)      % Ensure you only enter this code if U is populated
                                                          % Classification
    [classifyResults_F,Table_F] = classifyDataStream('classify_F', ...
        AmatrixF_ts, AmatrixF_tr, AmatrixN_tr, U, 1);
    if(do_save == 1)
    writetable(Table_F, 'ClassificationRanks_F.txt','Delimiter',',');
    save('ClassificationResults_F.mat', 'classifyResults_F');
    movefile('ClassificationRanks_F.txt', 'classification');
    movefile('ClassificationResults_F.mat', 'classification');
    end
    
    [classifyResults_N,Table_N] = classifyDataStream('classify_N', ...
        AmatrixN_ts, AmatrixF_tr, AmatrixN_tr, U, 1);
    if(do_save == 1)
    writetable(Table_N, 'ClassificationRanks_N.txt','Delimiter',',');
    save('ClassificationResults_N.mat', 'classifyResults_N');
    movefile('ClassificationRanks_N.txt', 'classification');
    movefile('ClassificationResults_N.mat', 'classification');
    end
    
    [classifyResults_U,Table_U] = classifyDataStream('classify_U', ...
        AmatrixU_ts, AmatrixF_tr, AmatrixN_tr, U, 1);
    if(do_save == 1)
    writetable(Table_U, 'ClassificationRanks_U.txt','Delimiter',',');
    save('ClassificationResults_U.mat', 'classifyResults_U');
    movefile('ClassificationRanks_U.txt', 'classification');
    movefile('ClassificationResults_U.mat', 'classification');
    end
    % The outputs are structure containing the Classifiaction Results
    
                                                              % Clustering
                              % Use the following functions to extract MFS
    featF_tr = getFeatureVectors(CovF_tr, U);       % Training: Functional
    featN_tr = getFeatureVectors(CovN_tr, U);   % Training: Non-Functional
    
    featF_ts = getFeatureVectors(CovF_ts, U);        % Testing: Functional
    featN_ts = getFeatureVectors(CovN_ts, U);    % Testing: Non-Functional
    featU_ts = getFeatureVectors(CovU_ts, U);           % Testing: Unknown
    
                                            % Plot the MFSP for each Mode!
    for mod = 1:size(U,2)
        cluster = figure; hold on;
        
        % Testing Data
        scatter(featF_ts.Fmatrix(2*(mod-1)+1,:),featF_ts.Fmatrix(2*mod,:),'o','r');
        scatter(featN_ts.Fmatrix(2*(mod-1)+1,:),featN_ts.Fmatrix(2*mod,:),'o','b');
        scatter(featU_ts.Fmatrix(2*(mod-1)+1,:),featU_ts.Fmatrix(2*mod,:),'o','g');
        
        % Training Data
        scatter(featF_tr.Fmatrix(2*(mod-1)+1,:),featF_tr.Fmatrix(2*mod,:),'.','r');
        scatter(featN_tr.Fmatrix(2*(mod-1)+1,:),featN_tr.Fmatrix(2*mod,:),'.','b');
        hold off; 
        
        legend(leg(1), leg(2), "Unknown: "+leg(3));
        xlabel("\mu: Mode "+mod); ylabel("\sigma: Mode "+mod);
        title("Feature Vector "+mod);
        
        if(do_save == 1)
        saveas(cluster, "Clustering_Vec"+mod, 'png');
        movefile("Clustering_Vec"+mod+".png", 'classification');
        end
    end
end
%}

% {
% Project it!
if (size(U,2) > 0)
                  % These functions project the data into projection space
                                                                % Training
    Project_F = dataStreamProjection('proj_F',AmatrixF_ts,U);
    Project_N = dataStreamProjection('proj_N',AmatrixN_ts,U);
    Project_U = dataStreamProjection('proj_U',AmatrixU_ts,U);
                                                                 % Testing
    Project_F_tr = dataStreamProjection('proj_F_tr',AmatrixF_tr,U);
    Project_N_tr = dataStreamProjection('proj_N_tr',AmatrixN_tr,U);
    
                                                         % Plot the figure
          % Cases for handing U spaces with 1, 2, or more modes is handled                                 
    figure(100);
    if (size(U,2) == 1)
        hold on;
        for k = 1:N_bs
            plot(Project_F.A{k},'r');
            plot(Project_N.A{k},'b');
            plot(Project_U.A{k},'g');
        end
        hold off; legend(leg);
        title("Projection Plots"); xlabel("Sample Index"); ylabel("SPLOC Mode 1");
    elseif (size(U,2) > 1)
        if (size(U,2) > 2)
            hold on;
            for k = 1:N_bs
                scatter3(Project_F.A{k}(1,:), Project_F.A{k}(2,:), Project_F.A{k}(3,:), 'r', 'd', 'filled')
                scatter3(Project_N.A{k}(1,:), Project_N.A{k}(2,:), Project_N.A{k}(3,:), 'b', 's', 'filled')
                scatter3(Project_U.A{k}(1,:), Project_U.A{k}(2,:), Project_U.A{k}(3,:), 'g', 'filled')
                
                scatter3(Project_F_tr.A{k}(1,:), Project_F_tr.A{k}(2,:), Project_F_tr.A{k}(3,:), 'r', 'd')
                scatter3(Project_N_tr.A{k}(1,:), Project_N_tr.A{k}(2,:), Project_N_tr.A{k}(3,:), 'b', 's')
            end
            hold off; legend(leg);
            title("Projection Plots"); xlabel("SPLOC Mode 1"); ylabel("SPLOC Mode 2");    
        else
            hold on;
            for k = 1:N_bs
                scatter(Project_F.A{k}(1,:), Project_F.A{k}(2,:), 'd', 'r', 'filled')
                scatter(Project_N.A{k}(1,:), Project_N.A{k}(2,:), 's', 'b', 'filled')
                scatter(Project_U.A{k}(1,:), Project_U.A{k}(2,:), 'g', 'filled')
                
                scatter(Project_F_tr.A{k}(1,:), Project_F_tr.A{k}(2,:), 'd', 'r')
                scatter(Project_N_tr.A{k}(1,:), Project_N_tr.A{k}(2,:), 's', 'b')
            end
            hold off; legend(leg);
            title("Projection Plots");
            xlabel("SPLOC Mode 1"); ylabel("SPLOC Mode 2"); zlabel("SPLOC Mode 3");
        end
    end
end
%}
 
% {
% Plot The Loadings For Each Mode
if (size(U,2) > 0)
    lege = [];
    modes = figure; hold on;
    for mode = 1:size(U, 2)
        lege = [lege, "Mode "+mode];
        plot(1:size(U,1), U(:,mode).^2);
    end; hold off;
    legend(lege); title("Squared Modes");
    xlabel("Original Feature Index"); ylabel("Squared Loading");
    if(do_save == 1)
        saveas(gcf, "SqModes", 'png');
        movefile("SqModes.png", 'analysis');
    end
end
%}

% Save Data (Optional)
if(do_save == 1)
    save('splocResults.mat', 'splocResults');
    saveas(figure(100), 'ProjectionPlot', 'fig');
end
