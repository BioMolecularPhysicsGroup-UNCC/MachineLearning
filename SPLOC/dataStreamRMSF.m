function RMSFinfo = dataStreamRMSF(Ainfo)
% calculates RMSF on the DOF of the system based on inputted A-datamatrix
%
% DEFINITIONS
% dataMatrixInfo. <-- data structure
% dataRefName    = reference name for A matrix data with similar traits
% dataMatrixName = cell array for file names that store the A matrix data 
%     nVariables = number of variables defining size of vector space
%          pType = packing type (format) applied to the vector space
%            dim = # of components in local vector (e.g. x,y,z for 1 atom)
%              n = number of data matrices
%           A{:} = cell array for the data matrices in the rowVar format
%    nSamples(:) = array for the number of samples in each data matrix
%
% -------------------------------------------------------- local variables 
% INPUT
% fname = file name to specify origin of data and/or identify the output.
% Ainfo = dataMatrixInfo containing one or more systems of interest. It is
%             particularly convenient to group all systems together that 
%             have a similar property.   
%
% PROCESS
% calculate the root mean square fluctuation based on the pType packing
% of the A-matrix. 
%
% OUTPUT  RMSFinfo. <-- data structure
% dataRefName    = inherits name of the A matrix used to create the RMSF
% dataMatrixName = cell array for file names that store the A matrix data 
%     nVariables = number of variables defining size of vector space
%              n = number of RMSF descriptions
%        rmsf{:} = cell array for the RMSF function & length(RMSF) = #DOF.
% REMARK:          the index for DOF automatically ranges from 1 to end
%
%%                                                    set sploc parameters
global gvSPLOC
% --------------------------------------------------- for recording action
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                 unpackage data matrices
N = Ainfo.n;                                        % # of unknown systems
A = Ainfo.A;                                    % cell array of A-matrices
aMname = Ainfo.dataMatrixName;           % cell array of data matrix names
pType = Ainfo.pType;
[nDOF,~] = size(A{1});   % nDOF = # of degrees of freedom = # of variables
dim = Ainfo.dim;
%%                                                          calculate RMSF
rmsf = cell(1,N);
   switch pType
   case 'notVectored'
      m = nDOF;
         for k=1:N
         rmsf{k} = std( A{k}, 0,2)';
         end
   case 'xyz-xyz-xyz'
      msf = cell(1,N);
         for k=1:N
         msf{k} = var( A{k}, 0,2)';
         end
      m = round(nDOF/dim);
        for k=1:N
        s2 = zeros(1,m);
           for j=1:dim
           index = j:dim:nDOF;
           s2 = s2 + msf{k}(index);
           end
        rmsf{k} = sqrt(s2);
        end
   case 'xxx-yyy-zzz'
      msf = cell(1,N);
         for k=1:N
         msf{k} = var( A{k}, 0,2)';
         end
      m = round(nDOF/dim);
        for k=1:N
        s2 = zeros(1,m);
        index0 = 1:m;
        shift = 0;
           for j=1:dim
           index = index0 + shift;
           s2 = s2 + msf{k}(index);
           shift = shift + m;
           end
        rmsf{k} = sqrt(s2);
        end
   otherwise
   error('pType is unknown');
   end
%%                                                      package dataOutput
RMSFinfo = struct;
RMSFinfo.dataRefName = Ainfo.dataRefName;
RMSFinfo.dataMatrixName = aMname;
RMSFinfo.nVariables = Ainfo.nVariables;
RMSFinfo.n = Ainfo.n;
RMSFinfo.rmsf = rmsf;
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
str = Ainfo.dataRefName;
fprintf(fid,'%s \n',['       input dataMatrixInfo name = ',str]);
fprintf(fid,'%s \n',['            output RMSFinfo name = ',str]);
fprintf(fid,'%s \n',['      dimension of local vectors = ',num2str(dim)]);
fprintf(fid,'%s \n',['   dimension of input A-matrices = ',num2str(nDOF)]);
fprintf(fid,'%s \n',['dimension of output RMSF vectors = ',num2str(m)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end

