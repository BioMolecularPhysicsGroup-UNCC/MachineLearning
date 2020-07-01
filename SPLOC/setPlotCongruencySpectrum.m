function setPlotCongruencySpectrum(splitSpectrum,showGMean)
% allows gvSPLOC.splocSplit and gvSPLOC.showGeometricMean to be reset
%
% INPUT
% splitSpectrum = (-1,1) => (split, not split)
%                 -1 => spectrum is shown as positive and negative
%                 1  => spectrum is shown only as positive (SPLOC only)
% showGMean     = (0,1) = (do not, do) show the geometrical mean line
%                 The geometric mean of the discriminant modes provides a 
%                 merit of figure for the discriminant selection power. 
%
global gvSPLOC
%%                                               check number of arguments
% ------------------------------------------------------------ three cases
flagChange = -1;                               % assumes no change is made
   switch nargin
     case 0                                       % reset back to defaults
       if( gvSPLOC.splocSplit ~= -1)
       flagChange = 1;
       end
       if( gvSPLOC.showGeometricMean ~= 0)
       flagChange = 1;
       end
     gvSPLOC.splocSplit = -1;   % default is to use split format for SPLOC
     gvSPLOC.showGeometricMean = 0;    % default is not to show Gmean line
     case 1              % changes splocSplit, preserves showGeometricMean
       if( gvSPLOC.splocSplit ~= splitSpectrum)
       flagChange = 1;
       end
       if( splitSpectrum == 1 )
       gvSPLOC.splocSplit = 1;      % now SPLOC spectrum will not be split
       else
       gvSPLOC.splocSplit = -1; % default is to use split format for SPLOC
       end
     case 2
       if( gvSPLOC.splocSplit ~= splitSpectrum)
       flagChange = 1;
       end
       if( splitSpectrum == 1 )
       gvSPLOC.splocSplit = 1;      % now SPLOC spectrum will not be split
       else
       gvSPLOC.splocSplit = -1; % default is to use split format for SPLOC
       end
       if( gvSPLOC.showGeometricMean ~= showGMean)
       flagChange = 1;
       end
       if( showGMean == 1 )
       gvSPLOC.showGeometricMean = 1;  % show geometrical mean for d-modes
       else
       gvSPLOC.showGeometricMean = 0;    % default: no show geometric mean
       end
   end   
%%                         reset setPlotCongruencySpectrum characteristics
   if(  flagChange > 0 )     %=> properties of congruency spectrum changed
   splocLogFile = gvSPLOC.splocLogFile;
   fid = fopen(splocLogFile,'a');                 % append new information
   fprintf(fid,'%s \n','  ');
% ---------------------------- write info on output file type for spectrum
   fprintf(fid,'%s \n',[mfilename,'()']);
   msg = dividerLine('congruency spectrum toggles');
   fprintf(fid,'%s \n',msg);
   msg = ['split SPLOC spectrum = ',num2str(gvSPLOC.splocSplit)];
   fprintf(fid,'%s \n',msg);
   msg = [' show geometric mean = ',num2str(gvSPLOC.showGeometricMean)];
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n',dividerLine);
   fclose(fid);  
   % else                             no notice of change needs to be made
   end
end
