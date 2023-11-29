%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File:    runSPLOC.m                           %%%                                   %%%
%%%  Purpose: Prepare training set for SPLOC       %%%    BioMolecular Physics Group     %%%
%%%           training then generate MFSPs         %%%   University of North Carolina    %%%
%%%  Created: 06-03-2023                           %%%           at Charlotte            %%%
%%%  Author:  Tyler J Grear                        %%%                                   %%%
%%----------------------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd SPLOC; disp("SPLOC: writing training data..."); splocStart = cputime;
n0 = size(class0,2); n1 = size(class1,2);  % Get number of class 0 & 1 data packets
if n0 == 1; class0 = horzcat(class0,class0); n0 = size(class0,2); end   % Duct tape solution
if ~exist("input",'dir'); mkdir("input"); end                           % Get OS file sep.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%================================| Prepare Data for SPLOC |================================%
fileID = fopen("input"+fsep+"DEMO_trainingSet",'w');  % Open training set txt
for i = 1:n0  % Write class 0 (functional) files for SPLOC input
    writematrix((class0{1,i})',"input"+fsep+"class0_"+class0{2,i}+".txt",'Delimiter','tab');
    fprintf(fileID,"F class0_"+class0{2,i}+".txt\n",'Delimeter','tab');
end
for j = 1:n1  % Write class 1 (nonfunctional) files for SPLOC input
    writematrix((class1{1,j})',"input"+fsep+"class1_"+class1{2,j}+".txt",'Delimiter','tab');
    fprintf(fileID,"N class1_"+class1{2,j}+".txt\n",'Delimeter','tab');
end
fclose('all'); disp("SPLOC: writing complete.");
disp("SPLOC: reading input data...");
initializeSPLOC("gtype","png","qType",qType);
[~,~,FnameF] = readFileNameList('DEMO_trainingSet','sType','F');  % Read class 0
[~,~,FnameN] = readFileNameList('DEMO_trainingSet','sType','N');  % Read class 1
[AmatF,~] = readDataMatrices('functional',FnameF,mFormat);      % Read F data
[AmatN,~] = readDataMatrices('nonfunctional',FnameN,mFormat);   % Read N data
CovF = getMultivariateStats4sploc(AmatF,n,2,'cov');  % Class 0 (F) training stats
CovN = getMultivariateStats4sploc(AmatN,n,2,'cov');  % Class 1 (N) training stats
disp("SPLOC: performing projection pursuit machine learning...");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==============| Supervised Projective Learning for Orthogonal Completeness |==============%
splocResults = sploc(pdgm,0,'DEMO',CovF,CovN,2);  % Run SPLOC
splocTime = cputime - splocStart; disp("SPLOC complete.");
disp("SPLOC CPU time: "+splocTime+" seconds ("+round(splocTime/60,2)+" minutes)");
plotCongruencySpectrum(splocResults,2,'basisSpectra');
movefile("basisSpectra",paths.outPath);
movefile("training"+fsep+"DEMO_sploc.log",paths.outPath);
close all
%==========================================================================================%
D = getDiscriminantSBV(splocResults);  % Get SPLOC discriminant basis vectors
U = getUndeterminedSBV(splocResults);  % Get SPLOC undetermined basis vectors
I = getIndifferenceSBV(splocResults);  % Get SPLOC indifference basis vectors
if MFSP == "Y"                         % If MFSPs are set to be generated
plotMFSP;                              % Plot mode-feature-space planes
end; cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%