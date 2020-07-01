% Run file for SPLOC IRIS Analysis
% Author: Chris Avery
% To run this code, SPLOC ToolSet must be added to the matlab path
% 1) Download SPLOC ToolSet from: https://github.com/BioMolecularPhysicsGroup-UNCC
% 2) Place SPLOC code with the rest of your matlab functions, this can be
%    anywhere accessible on your computer
% 3) Set matlab path using addpath( genpath( path/to/SPLOCToolSet/ ) )
%    (if you place this line in your startup.m file in matlabs root
%     directory then this will be done automatically upon startup!)
% 4) Run this script from the same directory as sploc_runner.m
%
% If do_save = 1 => All data will be saved and cataluged in a diretory
% called 'data/'

clear all; clc; clf; close all;                 % Start with a clean slate 
trials = 6;                     % Select number of trials to do (trials=6)
rngnum = 1:trials;                             % Make seeds to control RNG
do_save = 0;                     % (1, 0) -> (save data, do not save data)

                                          % Select the compairison to make
                                    % [functional, nonfucntional, unknown]
chooseClasses = [1,2,3];                           % (Setosa vs Virginica)
%chooseClasses = [1,3,2];                         % (Setosa vs Versicolor)
%chooseClasses = [2,3,1];                      % (Virginica vs Versicolor)

train_set = 1:25;                  % Let 25 samples make the training pool
test_set = max(train_set)+1:50;      % Let the rest be in the testing pool

[iris, Y] = iris_dataset;                 % Load the Iris Data from matlab

                     % Extract the classes (these are known ahead of time)
setosa = iris(:,1:50);
versic = iris(:,51:100);
vergin = iris(:,101:150);

                                                  % Do Each of the Trials!
for trial = 1:trials
    %% Set Training/Testing Data
    rng(rngnum(trial));              % Set Random Seed for Reproducability
    
    %disp(rand(1,3)) % Check seed is working
    
    leg = ["Setosa", "Virginica", "Versicolor"]; % Base classes (no touch)
    tr_ts = randperm(50);                                % Shuffle Classes
    
                                               % Extract the training data
    setosa_train = setosa(:,tr_ts(train_set));
    versic_train = versic(:,tr_ts(train_set));
    vergin_train = vergin(:,tr_ts(train_set));
    N_tr = length(tr_ts(train_set));

                                                % Extract the testing data
    setosa_test = setosa(:,tr_ts(test_set));
    versic_test = versic(:,tr_ts(test_set));
    vergin_test = vergin(:,tr_ts(test_set));
    N_ts = length(tr_ts(test_set));

    data_train = {setosa_train, versic_train, vergin_train};
    data_test = {setosa_test, versic_test, vergin_test};

    %% Create Bootstrap Samples
    num_bootstrap = 30;                   % Number of Bootstrapped Samples
    samp_siz = 15;                % Select number of samples per bootstrap
    n_k = floor(50/samp_siz);
    start = 1:samp_siz:50;
    ends = start+samp_siz-1;
    n_samples = num_bootstrap;

                                        % Create the Bootstrapped Samples!
                                        % Input directory is created here
    if(exist('input','dir')); rmdir('input','s'); end 
    mkdir('input');
    for j = 1:3
        for boot = 1:num_bootstrap
            number = boot;
            name1 = leg(j)+"_train_"+number+".txt";
            name2 = leg(j)+"_test_"+number+".txt";

            tmp1 = data_train{j}(:,randperm(N_tr, samp_siz)); % training
            dlmwrite(name1, tmp1, '\t');
            movefile(name1, 'input');

            tmp2 = data_test{j}(:,randperm(N_ts, samp_siz)); % testing
            dlmwrite(name2, tmp2, '\t');
            movefile(name2, 'input');
        end
    end

    %% Create input file
    leg = leg(chooseClasses);              % Permute  to correct class IDs

    fileID_tr = fopen('trainingdata', 'w');    % This is used for training
    fileID_ts = fopen('testingdata', 'w');      % This is used for testing
    %for flower = leg
        for i = 1:boot
            fprintf(fileID_tr, "F "+leg(1)+"_train_"+i+".txt\n");
            fprintf(fileID_tr, "N "+leg(2)+"_train_"+i+".txt\n");
            fprintf(fileID_tr, "U "+leg(3)+"_train_"+i+".txt\n");

            fprintf(fileID_ts, "F "+leg(1)+"_test_"+i+".txt\n");
            fprintf(fileID_ts, "N "+leg(2)+"_test_"+i+".txt\n");
            fprintf(fileID_ts, "U "+leg(3)+"_test_"+i+".txt\n");
        end
        %}
    %end
    fclose(fileID_tr);
    fclose(fileID_ts);
    movefile('trainingdata', 'input');
    movefile('testingdata', 'input');

    %% Train and Test SPLOC
    sploc_runner;    % This is the driver that defines the sploc workflow

    %% Clean Up
    % Move Files
    if(do_save == 1)
        fclose('all');
        outfile = "run"+trial;
        movefile('analysis', outfile);
        movefile('training', outfile);
        if(size(U,2) > 0); movefile('DiscriminateModes.txt', outfile);
        movefile('ProjectionPlot.fig', outfile);
        movefile('ProjectionPlot.png', outfile);
        movefile('splocResults.mat', outfile);
        movefile('input', outfile);
        movefile('classification', outfile); end
    rmdir('splocLog','s');
    outoutfile = "data\"+leg(1)+"_"+leg(2)+"\";
    movefile("run"+trial, outoutfile)
    end
end

fclose all;                               % Close any files left unclosed
