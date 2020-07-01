function writeSPLOCresults(fileName,splocResults)
% Write sploc results (a data structure) to the specified file: fileName
% INPUT
% fileName = name to output barebones information (may include full path)
%
%         splocResults <-- data structure
% ========================================================================
% splocResults.sType       = type of spectrum: MCSPLOC
% splocResults.pursuitType = 1,0,-1  => d, d&i, i  set DEFAULT VALUE = 0
% splocResults.pType       = packing format describing the vector space
% splocResults.dim         = # of components in a local vector
% splocResults.baseFname   = base file name
% splocResults.SBV         = selection basis vectors
% splocResults.vT          = voting threshold to establish consensus
% splocResults.efficacy    = quantifies the ability for SBV to cluster
% splocResults.Dd          = # of discriminant only modes
% splocResults.Ddi         = # of (discriminant & indifference) modes
% splocResults.Di          = # of indifference only modes
% splocResults.Du          = # of undetermined modes
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% splocResults.EEVd = efficacy eigenvalues for discriminant congruence
% splocResults.EEVi = efficacy eigenvalues for indifference congruence 
% splocResults.USVd = conditional upper similar discriminant-eigenvalues
% splocResults.USVi = conditional upper similar indifference-eigenvalues
% splocResults.LSVd = conditional lower similar discriminant-eigenvalues
% splocResults.LSVi = conditional lower similar indifference-eigenvalues
% splocResults.QEVd  = conditonal quality discriminant-eigenvalues
% splocResults.QEVi  = conditonal quality indifference-eigenvalues
% splocResults.SEVd  = conditional selection discriminant-eigenvalues
% splocResults.SEVi  = conditional selection indifference-eigenvalues 
% splocResults.CEVd  = conditional consensus discriminant-eigenvalues 
% splocResults.CEVi  = conditional consensus indifference-eigenvalues 
% splocResults.Cind  = congruency indicator (2,1,0,-1)
%                    2 => discriminant and indifference congruences
%                    1 => projections for discriminant congruences
%                    0 => projections that are undetermined
%                   -1 => projections for indifference congruences
%%                                                prepare to record action
global gvSPLOC
splocLogFile = gvSPLOC.splocLogFile;
%%                                                simple consistency check
Dd = splocResults.Dd;
Ddi= splocResults.Ddi;
Du = splocResults.Du;
Di = splocResults.Di;
nV = Dd + Ddi + Du + Di;
SBV = splocResults.SBV;                          % selection basis vectors
[m1,m2] = size(SBV);
mCorrect = 2;
   if( m1 ~= nV )
   mCorrect = mCorrect - 1;
   end
     if( m2 ~= nV )
     mCorrect = mCorrect - 1;
     end
   if( mCorrect == 1 )
   error('SBV is not a square matrix');
   elseif( mCorrect == 0 )
   error(['SBV size is incompatiable with nV = ',num2str(nV)]);
   end
%%                              construct matrix for barebones information
EEVd = splocResults.EEVd;
USVd = splocResults.USVd;
LSVd = splocResults.LSVd;
QEVd = splocResults.QEVd;
SEVd = splocResults.SEVd;
CEVd = splocResults.CEVd;
Cind = splocResults.Cind;
EEVi = splocResults.EEVi;
USVi = splocResults.USVi;
LSVi = splocResults.LSVi;
QEVi = splocResults.QEVi;
SEVi = splocResults.SEVi;
CEVi = splocResults.CEVi;
barebonesInfo = [EEVd; USVd; LSVd; QEVd; SEVd; CEVd; Cind; ...
                 EEVi; USVi; LSVi; QEVi; SEVi; CEVi; SBV];
%%                               write information in splocResults to file
   if( contains(fileName,'_splocResults.dlm') )
   fOUT = fileName;
   else
   fOUT = [fileName,'_splocResults.dlm'];
   end
% splocResults.sType       = type of spectrum: MCSPLOC
% splocResults.pursuitType = 1,0,-1  => d, d&i, i  set DEFAULT VALUE = 0
% splocResults.pType       = packing format describing the vector space
% splocResults.dim         = # of components in a local vector
% splocResults.baseFname   = base file name
% splocResults.SBV         = selection basis vectors
% splocResults.vT          = voting threshold to establish consensus
% splocResults.efficacy    = quantifies the ability for SBV to cluster
% splocResults.Dd          = # of discriminant only modes
% splocResults.Ddi         = # of (discriminant & indifference) modes
% splocResults.Di          = # of indifference only modes
% splocResults.Du          = # of undetermined modes
%
fileID = fopen(fOUT,'w');
fprintf(fileID,'%s\n',splocResults.sType);                             % 1
fprintf(fileID,'%d\n',splocResults.pursuitType);                       % 2
fprintf(fileID,'%s\n',splocResults.pType);                             % 3
fprintf(fileID,'%d\n',splocResults.dim);                               % 4
fprintf(fileID,'%s\n',splocResults.baseFname);                         % 5
fprintf(fileID,'%13.11f\n',splocResults.vT);                           % 6
fprintf(fileID,'%19.5f\n',splocResults.efficacy);                     % 7
fprintf(fileID,'%d\n',Dd);                                             % 8
fprintf(fileID,'%d\n',Ddi);                                            % 9
fprintf(fileID,'%d\n',Du);                                            % 10
fprintf(fileID,'%d\n',Di);                                            % 11
fclose(fileID);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dlmwrite(fOUT,barebonesInfo,'precision',10,'-append');
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
fprintf(fid,'%s \n', dividerLine('writing current splocResults') );
msg = ['output file = ',fileName];
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);  
end
