%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File:    prepEO.m                             %%%                                   %%%
%%%  Purpose: Preprocess ECG Data from Ethyln      %%%    BioMolecular Physics Group     %%%
%%%           and perform SPLOC for classification %%%   University of North Carolina    %%%
%%%  Created: 07-14-2023                           %%%           at Charlotte            %%%
%%%  Authors: Tyler J Grear and Ethelyn Ofei       %%%                                   %%%
%%----------------------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------| Load ECG Data |-------------------------------------%
F_dir = "Functional";     % Set dir to pull functional recording sessions
N_dir = "Hypertrophy";  % Set dir to pull nonfunctional recording sessions
F_files = struct2cell(dir("Data"+fsep+"Labeled ECG Data"+fsep+F_dir+fsep+"*.mat"));  % Get F info
N_files = struct2cell(dir("Data"+fsep+"Labeled ECG Data"+fsep+N_dir+fsep+"*.mat"));  % Get N info
nF = size(F_files,2);  % Number of functional traces
nN = size(N_files,2);  % Number of nonfunctional traces
A0 = cell(1,nF);        % Initialize cell to populate with initial F data matrices
A1 = cell(1,nN);        % Initialize cell to populate with initial N data matrices
 
for i = 1:nF  % Loop over the number of functional recording sessions (raw data)
    A_i = struct2cell(load("Data"+fsep+"Labeled ECG Data"+fsep+F_dir+fsep+F_files{1,i}));  % Load i-th F tr
    A0{1,i} = (cell2mat(A_i));                                     % Store i-th F tr
end

for j = 1:nN  % Loop over the number of nonfunctional recording sessions (raw data)
    B_i = struct2cell(load("Data"+fsep+"Labeled ECG Data"+fsep+N_dir+fsep+N_files{1,j}));  % Load j-th N tr
    A1{1,j} = (cell2mat(B_i));                                     % Store i-th N tr
end

%A0(40) = [];
for i = 1:size(A0,2)
    session_i = A0{1,i};
    n_i = size(session_i,2);
    if n_i > 115200
        A0{1,i} = session_i(:,1:115200);
    end
end
for j = 1:size(A1,2)
    session_j = A1{1,j};
    n_j = size(session_j,2);
    if n_j > 115200
        A1{1,j} = session_j(:,1:115200);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%