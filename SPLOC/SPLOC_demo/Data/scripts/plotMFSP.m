%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  File:    plotMFSP.m                           %%%                                   %%%
%%%  Purpose: Generate mode-feature-space planes   %%%    BioMolecular Physics Group     %%%
%%%           from SPLOC basis                     %%%   University of North Carolina    %%%
%%%  Created: 06-05-2023                           %%%           at Charlotte            %%%
%%%  Author:  Tyler J Grear                        %%%                                   %%%
%%----------------------------------------------------------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = horzcat(D,U,I); disp(" ");        % Direct sum of subspace D and I
if isempty(V) == 0; disp("Constructing MFSPs..."); end
modeNames = []; dim = size(V,2);
disp("Number of d-modes: "+size(D,2))
dnames = repmat({'d-mode'},1,size(D,2));
disp("Number of u-modes: "+size(U,2))
unames = repmat({'u-mode'},1,size(U,2));
disp("Number of i-modes: "+size(I,2))
inames = repmat({'i-mode'},1,size(I,2));
modeNames = horzcat(dnames,unames,inames);
featF = getFeatureVectors(CovF,V);  % Extract class 0 (functional) covariance matrix
featN = getFeatureVectors(CovN,V);  % Extract class 1 (nonfunctional) covariance matrix
if ~exist(outPath+fsep+"MFSPs",'dir')
    mkdir(outPath+fsep+"MFSPs");
end
cSpectrum = imread(outPath+fsep+"basisSpectra"+fsep+"DEMO_congruencySpectrum.png");
%==========================================================================================%
for mode = 1:size(V,2)  % Loop over number of column vectors in V and plot MFSPs
    figure(42)
    scatter(featF.Fmatrix(2*(mode-1)+1,:),featF.Fmatrix(2*mode,:),mkSize,'o', ...
    'MarkerFaceColor','red','MarkerEdgeColor','red'); hold on
    scatter(featN.Fmatrix(2*(mode-1)+1,:),featN.Fmatrix(2*mode,:),mkSize,'o', ...
    'MarkerFaceColor','blue','MarkerEdgeColor','blue');
    legend('Class 0','Class 1','Location','southoutside','Orientation','horizontal');
    xlabel("Mode "+mode+" \mu"); ylabel("Mode "+mode+" \sigma");
    xtickformat('%1d'); ytickformat('%05.2f'); fontsize(gcf,8,"points")
    title("Feature vector "+mode+": "+modeNames{1,mode}); axis square
    exportgraphics(gcf,"temp"+mode+".png","Resolution",res);
    mfsp = imread("temp"+mode+".png");
    nrows = size(mfsp,1);
    cSpectrum = imresize(cSpectrum,[nrows NaN]);
    figure(43)
    montage({mfsp,cSpectrum},'Size',[1 2],'ThumbnailSize',[nrows NaN],'BackgroundColor','white')
    set(gcf,'Units','Normalized','OuterPosition',[0 0 1 1]); delete("temp"+mode+".png");
    exportgraphics(gcf,outPath+fsep+"MFSPs"+fsep+"MFSP"+mode+".png","Resolution",res);
    hold off; close all;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%