function splocResults = readSPLOCresults(fname,fullPath)
% reads sploc results (as data structure) from the file specified by fname
%
% INPUT
% fname = name of file without standard subdirectory or ending. fullPath=0
%      OR complete name using relative path or absolute path.   fullPath=1
% fullPath: 0 => look for inputFname file within the training directory
%           1 => assume inputFname contains the full path of interest.
%
% OUTPUT
% barebones splocResults datastructure that is populated
% 
% NOTE: Some information about source of data is missing and cannot be
%       recovered from the retrieval of barebones information
%
%         splocResults <-- data structure
% ========================================================================
% splocResults.sType       = type of spectrum: MCSPLOC
% splocResults.pursuitType = 1,0,-1  => d, d&i, i  set DEFAULT VALUE = 0
% splocResults.pType       = packing format describing the vector space
% splocResults.baseFname   = base file name
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^-> lost above info
% splocResults.dim         = # of components in a local vector
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
%%                                                   standardize file name
   if( contains(fname,'_splocResults.dlm') )
   fIN = fname;                                   %=> have standard ending
   else
   fIN = [fname,'_splocResults.dlm'];             % append standard ending
   end
%%                                                               get input
   if( nargin == 1 ) 
   fullPath = 0;                                           % default value
   end
   if( fullPath == 0 )                 % process file naming the sploc way
   targetFolder = 'training';
   fName = getInputFileName(targetFolder,fIN);         % add standard path
   else
   fName = fIN;
   end
%%                                                               read data
fileID = fopen(fName,'r');
sType = fgetl(fileID); 
pursuitType = str2num( fgetl(fileID) );
pType = fgetl(fileID);
dim = str2num( fgetl(fileID) );
baseFname = fgetl(fileID);
vT = str2double( fgetl(fileID) );
efficacy = str2double( fgetl(fileID) );
Dd = str2num( fgetl(fileID) );
Ddi = str2num( fgetl(fileID) );
Du = str2num( fgetl(fileID) );
Di = str2num( fgetl(fileID) );
fclose(fileID);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
barebonesInfo = dlmread(fName,',',11,0);
% size(barebonesInfo)
% barebonesInfo = [EEVd; USVd; LSVd; QEVd; SEVd; CEVd; Cind; ...
%                  EEVi; USVi; LSVi; QEVi; SEVi; CEVi; SBV];
EEVd = barebonesInfo( 1,:);
USVd = barebonesInfo( 2,:);
LSVd = barebonesInfo( 3,:);
QEVd = barebonesInfo( 4,:);
SEVd = barebonesInfo( 5,:);
CEVd = barebonesInfo( 6,:);
Cind = barebonesInfo( 7,:);
EEVi = barebonesInfo( 8,:);
USVi = barebonesInfo( 9,:);
LSVi = barebonesInfo(10,:);
QEVi = barebonesInfo(11,:);
SEVi = barebonesInfo(12,:);
CEVi = barebonesInfo(13,:);
SBV  = barebonesInfo(14:end,:);
%%                              build barebones splocResults datastructure
splocResults = struct;
splocResults.sType = sType;   
splocResults.pursuitType = pursuitType;
splocResults.pType = pType;
splocResults.dim = dim;
splocResults.baseFname = baseFname;
splocResults.SBV = SBV;                          % selection basis vectors
splocResults.vT = vT;
splocResults.efficacy = efficacy;
splocResults.Dd = Dd;
splocResults.Ddi= Ddi;
splocResults.Du = Du;
splocResults.Di = Di;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
splocResults.EEVd = EEVd;
splocResults.EEVi = EEVi;
splocResults.USVd = USVd;
splocResults.USVi = USVi;
splocResults.LSVd = LSVd;
splocResults.LSVi = LSVi;
splocResults.QEVd = QEVd;
splocResults.QEVi = QEVi;
splocResults.SEVd = SEVd;
splocResults.SEVi = SEVi;
splocResults.CEVd = CEVd;
splocResults.CEVi = CEVi;
splocResults.Cind = Cind;


% fprintf(fileID,'%s\n',splocResults.sType);                             % 1
% fprintf(fileID,'%d\n',);                       % 2
% fprintf(fileID,'%s\n',splocResults.pType);                             % 3
% fprintf(fileID,'%d\n',splocResults.dim);                               % 4
% fprintf(fileID,'%s\n',splocResults.baseFname);                         % 5
% fprintf(fileID,'%13.11f\n',splocResults.vT);                           % 6
% fprintf(fileID,'%19.10f\n',splocResults.efficacy);                     % 7
% fprintf(fileID,'%d\n',Dd);                                             % 8
% fprintf(fileID,'%d\n',Ddi);                                            % 9
% fprintf(fileID,'%d\n',Du);                                            % 10
% fprintf(fileID,'%d\n',Di);                                            % 11


% ------------------------------------------------------------ error check
[~,nV] = size(barebonesInfo);            % nV = number of variables = nDOF
% ----------------------------------------------- find congruent subspaces
L3 = ( Cind > 0 );                                   % captures Dd and Ddi
L2 = ( Cind == 2 );                  % Ddi: dual-purpose (d & i) subspaces
Ld = ( Cind == 1 );                                  % Dd: d subspace only
Lu = ( Cind == 0 );                                  % Du: u subspace only
Li = ( Cind == -1);                                  % Di: i subspace only
dd = sum(Ld);                               % # of discriminant only modes
d2 = sum(L2);      % # of dual-purpose (discriminant & indifference) modes
du = sum(Lu);                                    % # of undetermined modes
di = sum(Li);                               % # of indifference only modes
nDOF = dd + d2 + du + di;
   if( nDOF ~= nV )
   error('Mismatch between vector space dimension and DOF');
   end
efficacy = sum( EEVd(L3) ) + sum( EEVi(L2) ) + sum( EEVi(Li) ); 
   if( abs(efficacy - splocResults.efficacy)  > 0.001 )
   disp(['from       Cind: efficacy = ',num2str(efficacy)]);
   disp(['from splcResults.efficacy = ',num2str(splocResults.efficacy)]);
   error('efficacy inconsistency');
   end
   if( Dd ~= dd )
   error('Dd is inconsistent with Cind');
   end
     if( Ddi ~= d2 )
     error('Ddi is inconsistent with Cind');
     end
   if( Du ~= du )
   error('Du is inconsistent with Cind');
   end
     if( Di ~= di )
     error('Di is inconsistent with Cind');
     end
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
fprintf(fid,'%s \n', dividerLine('basic summary') );
msg = ['input file = ',fname];
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);  
end