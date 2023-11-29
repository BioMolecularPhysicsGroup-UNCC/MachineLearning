function dataOutput = dataStreamFilter(fname,dataInput,U)
% Filters out information from the A-matrix that is outside the U subspace
% A projection operator is created as U*U' which is then applied to A. The
% output data matrices can be used in the same what that the original A 
% matrices can be used, except the parts of information orthogonal to the
% properties of interest are projected out. Notice that no shift is done.
% A user may want to shift the origin of coordinates in various ways, so
% this projection function makes no assumptions, and outputs non shifted
% data, which is true to the inputted data.
%
% Afiltered = (U*U')*Ainput
% ------------------------------------------------------------------------
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
% dataInput = dataMatrixInfo for one or more systems of interest. It would
%             be particularly convenient to group all systems together 
%             that have a similar property. 
% U = basis vectors employed to project data into subspace.   
%
% PROCESS
% create a projection matrix to project all vectors into the subspace that
% is defined by the U basis vectors. Then apply the projection matrix on
% the A matrix to arrive at the filtered A-matrix, which removes the 
% information orthogonal to the subspace defined by U.
%
% OUTPUT  dataMatrixInfo. <-- data structure
% The same type of data structure is outputed, with virtually no change in 
% meta-data descriptions. The A matrices represent the filtered versions.
%
%%                                                    set sploc parameters
global gvSPLOC
% --------------------------------------------------- for recording action
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                 unpackage data matrices
N = dataInput.n;                                    % # of unknown systems
A = dataInput.A;                                % cell array of A-matrices
aMname = dataInput.dataMatrixName;       % cell array of data matrix names
%%                                                      simple error check
[nDOF,~] = size(A{1});   % nDOF = # of degrees of freedom = # of variables
[nV,nModes] = size(U);         % nModes = # of basis vectors to be checked
   if( nDOF ~= nV )
   disp([' A:   nDOF = ',num2str(nDOF)]);
   disp([' U:   nDOF = ',num2str(nV)]);
   disp([' U: nmodes = ',num2str(nModes)]);
   error('number of variables differ between A and basis set dimension');
   end
%%                                                        do the filtering
P = U*U';                     %  P = projection matrix   = 1 when complete
   for k=1:N
   A{k} = P*A{k};    % all motions orthogonal to the U-subspace is removed
   end
%%                                                      package dataOutput
dataOutput = struct;
dataOutput.dataRefName = fname;
dataOutput.dataMatrixName = aMname;
dataOutput.nVariables = dataInput.nVariables;
dataOutput.pType = dataInput.pType;             % format type is preserved
dataOutput.dim = dataInput.dim;
dataOutput.n = dataInput.n;
   for k=1:N
   dataOutput.A{k} = A{k};
   end
dataOutput.nSamples = dataInput.nSamples;
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
str = dataInput.dataRefName;
fprintf(fid,'%s \n',['     input dataMatrixInfo name = ',str]);
fprintf(fid,'%s \n',['    output dataMatrixInfo name = ',fname]);
fprintf(fid,'%s \n',[' dimension of input A-matrices = ',nDOF]);
fprintf(fid,'%s \n',['dimension of output A-matrices = ',nModes]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end

