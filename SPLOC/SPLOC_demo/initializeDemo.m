%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File:    initializeDemo.m                    %%%                                   %%%
%%%  Purpose: Initialize demo parameters & dirs   %%%    BioMolecular Physics Group     %%%
%%%           then perform miscellaneous tasks    %%%   University of North Carolina    %%%
%%%  Created: 07-16-2023                          %%%           at Charlotte            %%%
%%%  Author:  Tyler J Grear and Ethelyn Ofei      %%%                                   %%%
%%---------------------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------| Initialize MCP Parameters |------------------------------%
                       %                                                                  %
preprocess = 0;        %  Activate correction for baseline wandering noise of ECG.        %
                       %                                                                  %
%-----------------------------------------------------------------------------------------%
                       %                                                                  %
res = 150;             %  Set the resolution of generated figures in dots-per-inch (DPI). %
                       %  The default resoution is 150 DPI, 300 DPI is print quality.     %
                       %                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------| SPLOC |----------------------------------------%
                       %                                                                  %
pdgm = 0;              %  Sets the SPLOC operational paradigm for DSP protocols. The      %
                       %  available options are:                                          %
                       %                                                                  %
                       %            2: Maximize d-modes aggressively.                     %
                       %            1: Maximize d-modes conservatively.                   %
                       %            0: Maximize d-modes & i-modes simultaneously.         %
                       %           -1: Maximize i-modes conservatively.                   %
                       %           -2: Maximize i-modes aggressively.                     %
                       %                                                                  %
%-----------------------------------------------------------------------------------------%
                       %                                                                  %
qType = 1;             %  Set the SPLOC decision boundary method for clustering quality   %
                       %  calculation. The available options are:                         %
                       %                                                                  %
                       %            1: Hard linear decision boundary.                     %
                       %            2: Soft geometry decision boundary.                   %
                       %            3: Hybrid, combination of 1 & 2.                      %
                       %                                                                  %
%-----------------------------------------------------------------------------------------%
                       %                                                                  %
MFSP = 'Y';            %  A setting of 'Y' (on) displays the mode-feature-space plane     %
                       %  (MFSP) output of SPLOC. A setting of 'N' (off) will disable     %
                       %  this feature.                                                   %
                       %                                                                  %
%-----------------------------------------------------------------------------------------%
                       %                                                                  %
mkSize = 8;            %  Set the marker size of scatter plot data. It is recommended     %
                       %  that mkSize is decreased as the # of data packets increases.    %
                       %                                                                  %
%-----------------------------------------------------------------------------------------%
                       %                                                                  %
rando = 0;             %  Set the marker size of scatter plot data. It is recommended     %
                       %  that mkSize is decreased as the # of data packets increases.    %
                       %                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------| DEMO Workflow |------------------------------------%
%                                                                                         %
%                                                                                         %
%                                                                                         %
%                _______      _______     _______     ______     __________               %
%               |       |    |       |   |       |   |      |   |          |              %
%               | Input |--> | Split |-->| SPLOC |-->| MFSP |-->| Analysis |              %
%               |_______|    |_______|   |_______|   |______|   |__________|              %
%                                                                                         %
%                                                                                         %
%                                                                                         %
%                                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------| Miscellaneous Initialization Tasks |--------------------------%
fsep = filesep;                      % Get operating system file separator
rtime = string(datetime('now'));     % Get current time
OUTDIR = "Output_"+rtime;
OUTDIR = regexprep(OUTDIR,' ','_');
OUTDIR = regexprep(OUTDIR,':','-');
if ~exist(OUTDIR,'dir')              % Check for output directory
    mkdir(OUTDIR);                   % Create directory if output dir does not exist
end
launchDir = string(pwd);                  % Get full path for launch directory
outPath = strcat(launchDir,fsep+OUTDIR);  % Get full path for output directory
paths.launchDir = launchDir; paths.outPath = outPath;  % Store paths for later use
addpath(launchDir); addpath(outPath);                               % Add search paths
addpath("Data"); addpath("Data"+fsep+"scripts"); addpath("SPLOC");  % Add search paths
%-----------------------------------------------------------------------------------------%
mFormat = struct; mFormat.pType = 'notVectored'; mFormat.dim = 1;
if qType == 1; qType = 'hardLinear'; end
if qType == 2; qType = 'softGeometry'; end
if qType == 3; qType = 'hybrid'; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%