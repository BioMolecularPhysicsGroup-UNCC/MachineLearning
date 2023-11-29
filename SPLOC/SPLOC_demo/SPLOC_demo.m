%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File:    SPLOC_demo.m                         %%%                                   %%%
%%%  Purpose: Preprocess ECG Data from Ethelyn     %%%    BioMolecular Physics Group     %%%
%%%           and perform SPLOC for classification %%%   University of North Carolina    %%%
%%%  Created: 07-14-2023                           %%%           at Charlotte            %%%
%%%  Authors: Tyler J Grear and Ethelyn Ofei       %%%                                   %%%
%%----------------------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------| Load ECG Data |-------------------------------------%
mcpStart = cputime;       % Get current CPU time
initializeDemo;           % Initialize DSP dirs and parameters
if preprocess == 1
    preprocessECG4SPLOC;  % Run preprocessing algorithm to correct for baseline wander
else
    prepDemo;               % Partition ECG recording sessions without preprocessing
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------| Parse ECG Traces into Data Packets |--------------------------%
n = 5000;                   % Set desired number of samples for each SPLOC data packet
[n_0,p] = size(A0{1,1}');   % This assumes all raw data matrices are the same shape
nF = size(A0,2);            % Get number of functional traces
nN = size(A1,2);            % Get number of nonfunctional traces
nPackets_0 = floor(n_0/n);  % Get # of continuous samples for 1 data packet
nPackets = nF*nPackets_0;   % Final number of data packets per class
start = 1:n:n_0;                        % Index for starting samples of data packets
start = start(1,1:(end-1));             % Remove last element of index
stop = n:n:n_0;                         % Index for stopping samples of data packets
class0_0 = cell(nF,nPackets_0);         % Initialize cell to populate with F data packets
class0_tags = cell(1,(nF*nPackets_0));  % Initialize cell to populate with F trace (packet) tags
class1_0 = cell(nN,nPackets_0);         % Initialize cell to populate with N data packets
class1_tags = cell(1,(nN*nPackets_0));  % Initialize cell to populate with N trace (packet) tags

for i = 1:nF                                    % Loop over # of F recording sessions (raw data)
    X_i = A0{1,i}';                             % Get i-th F recording session data
    for j = 1:nPackets_0                        % Loop over # of desired F packets per session
        X_j = X_i((start(1,j):(stop(1,j))),:);  % Extract j-th F data packet from i-th session
        class0_0{i,j} = X_j;                    % Store j-th F data packet
    end
    if rando == 1
        rando_F = randperm(nPackets_0,nPackets_0);
        class0_0(i,:) = class0_0(i,rando_F);
    end
    if i == nF                                   % If last recording session, do:
        class0_0 = reshape(class0_0,1,[]);       % Squish into singleton cell array
        for k = 1:(nF*nPackets_0)                % Loop over total number of packets
            class0_tags{1,k} = "tr"+k;           % Add trace tag to cell array
        end
        class0 = vertcat(class0_0,class0_tags);  % Stack F traces and tags for final F systems
    end
end

for i = 1:nN                                    % Loop over # of N recording sessions (raw data)
    X_i = A1{1,i}';                             % Get i-th N recording session data
    for j = 1:nPackets_0                        % Loop over # of desired N packets per session
        X_j = X_i((start(1,j):(stop(1,j))),:);  % Extract j-th N data packet from i-th session
        class1_0{i,j} = X_j;                    % Store j-th N data packet
    end
    if rando == 1
        rando_N = randperm(nPackets_0,nPackets_0);
        class1_0(i,:) = class1_0(i,rando_N);
    end
    if i == nN                                   % If last recording session, do:
        class1_0 = reshape(class1_0,1,[]);       % Squish into singleton cell array
        for k = 1:(nN*nPackets_0)                % Loop over total number of packets
            class1_tags{1,k} = "tr"+k;           % Add trace tag to cell array
        end
        class1 = vertcat(class1_0,class1_tags);  % Stack N traces and tags for final N systems
    end
end

if rando == 1
rando_F = randperm(size(class0,2),size(class0,2));  % Get index for random F packet selection
rando_N = randperm(size(class1,2),size(class1,2));  % Get index for random N packet selection
class0 = class0(:,rando_F);
class1 = class1(:,rando_N);
end
         
nRetain_F = 25;  % !!!! =-> TG: Taking subset just so I can run for debugging, remove eventually
nRetain_N = 25;  % !!!! =-> TG: Taking subset just so I can run for debugging, remove eventually
class0 = class0(:,1:nRetain_F);
class1 = class1(:,1:nRetain_N);

disp("Samples-to-dimensions ratio (n/p): "+(n)/p);
disp("Number of class 0 data packets: "+size(class0,2))
disp("Number of class 1 data packets: "+size(class1,2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------| Supervised Projective Learning for Orthogonal Completeness |--------------%
runSPLOC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcpTime = cputime - mcpStart;
disp("SPLOC demo CPU time: "+mcpTime+" seconds ("+round(mcpTime/60,2)+" minutes)"); disp(" ")
BMPG_Quote_Generator;  % Print random BioMolecular Physics Group quote
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%