% Run file for SPLOC WINE Analysis
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
num_trial = 6;                  % Select number of trials to do (trials=6)
rngnum = 1:num_trial;                          % Make seeds to control RNG
do_save = 0;                     % (1, 0) -> (save data, do not save data)

                                          % Select the compairison to make
                                    % [functional, nonfucntional, unknown]
chooseClasses = [1,2,3];                           % (Wine 1 vs Wine 2)
%chooseClasses = [1,3,2];                          % (Wine 1 vs Wine 3)
%chooseClasses = [2,3,1];                          % (Wine 2 vs Wine 3)

variables = ["Alcohol","Malic acid","Ash","Alcalinity of ash","Magnesium",...
    "Total phenols","Flavanoids","Nonflavanoid phenols","Proanthocyanins",...
    "Color intensity","Hue","OD280/OD315","Proline"];

leg = ["Wine1", "Wine2", "Wine3"];
train = 25;                            % Number of samples in training set

[x_data, labels] = wine_dataset;          % Load the Wine data from matlab
N_class = size(labels, 1);                     % Extract Number of Classes 
N_tot = size(labels, 2);                 % Extract Total Number of Samples
N_feat = size(x_data, 1);                     % Extract Number of Features

leg = leg(chooseClasses);         % Permute classes to correct compairison
for do_trial = 1:num_trial                                 % Do each trial
    rng(rngnum(do_trial));           % Set random seed for reproducability
    
    %% Extract Classes From Data
    wine1 = x_data(:, find(labels(1,:)==1));              % Extract Wine 1
    wine2 = x_data(:, find(labels(2,:)==1));              % Extract Wine 2
    wine3 = x_data(:, find(labels(3,:)==1));              % Extract Wine 3

    data = {wine1, wine2, wine3};
    data_tr = cell(1,length(data));         % This will hold training data
    data_ts = cell(1,length(data));          % This will hold testing data
    N = [];
    for i = 1:N_class
        tmp = data{i};
        N(i) = size(tmp, 2);
        inx = randperm(N(i));
        data_tr{i} = tmp(:, inx(1:train));
        data_ts{i} = tmp(:, inx(train+1:N(i)));
    end
    %clear inx tmp wine1 wine2 wine3 % Time to clean up!
    
                                    % Permute to correct class compairison
    data_tr = data_tr(chooseClasses);
    data_ts = data_ts(chooseClasses);

    %% Create Bootstrapped Samples
    if(exist('input','dir')); rmdir('input','s'); end    %Create input dir
    mkdir('input');

    N_bs = 30;                    % Number of bootstrapped samples to make
    samp_siz = 15;                       % Number of samples per bootstrap
    
    for i = 1:N_class           % Create Training and Testing data packets
        tmp_tr = data_tr{i};
        tmp_ts = data_ts{i};
        for j = 1:N_bs
            inxtr = randperm(size(tmp_tr,2), samp_siz);
            inxts = randperm(size(tmp_ts,2), samp_siz);
            dlmwrite("train_"+leg(i)+"_sample"+j, tmp_tr(:,inxtr), '\t');
            dlmwrite("test_"+leg(i)+"_sample"+j, tmp_ts(:,inxts), '\t');
            movefile("train_"+leg(i)+"_sample"+j, 'input');
            movefile("test_"+leg(i)+"_sample"+j, 'input');
        end
    end

    fileID_tr = fopen('trainingdata', 'w');   % The file for training data
    fileID_ts = fopen('testingdata', 'w');     % The file for testing data
    for i = 1:N_bs
            fprintf(fileID_tr, "F "+"train_"+leg(1)+"_sample"+i+"\n");
            fprintf(fileID_tr, "N "+"train_"+leg(2)+"_sample"+i+"\n");
            fprintf(fileID_tr, "U "+"train_"+leg(3)+"_sample"+i+"\n");

            fprintf(fileID_ts, "F "+"test_"+leg(1)+"_sample"+i+"\n");
            fprintf(fileID_ts, "N "+"test_"+leg(2)+"_sample"+i+"\n");
            fprintf(fileID_ts, "U "+"test_"+leg(3)+"_sample"+i+"\n");
    end
    fclose(fileID_tr);
    fclose(fileID_ts);
    movefile('trainingdata', 'input');       % Move to the input directory
    movefile('testingdata', 'input');        %                  for SPLOC!

    %
    %% Train and Test SPLOC
    sploc_runner;        % This driver file runs the basic SPLOC workflow!

    %% Clean Up
    % Move Files
    if(do_save == 1)
        fclose('all');
        outfile = "data\"+leg(1)+"_"+leg(2)+"\run"+do_trial+"\";
        movefile('analysis', outfile);
        movefile('training', outfile);
        movefile('DiscriminateModes.txt', outfile);
        movefile('ProjectionPlot.fig', outfile);
        movefile('splocResults.mat', outfile);
        movefile('input', outfile);
        if(size(U,2) > 0); movefile('classification', outfile); end
        rmdir('splocLog','s');
        close all;
    end
end
fclose all;
%}