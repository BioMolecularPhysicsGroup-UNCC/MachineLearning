function dataOutput = dataStreamProjection(fname,dataInput,U)
% projects A-matrix data into a subspace defined by basis vectors in U.
% Given N degrees of freedom and D-observations, the NxD data matrix will
% be projected into a subspace comprising M degrees of freedom associated
% with a generalized coordinate per mode (or basis vector). If the U
% matrix was constructed from PCA modes, then these projections would be
% essentially principal components. The output data matrices can be used
% to create (1) fuzzballs by plotting one row of data against another, or 
% (2) time sequence plots by taking any row to represent a time sequence,
% or (3) different rows of the output data matrix can be used to calculate
% PDFs of the generalized coordinates. Notice that no shifting is done.
% A user may want to shift the origin of coordinates in various ways, so
% this projection function makes no assumptions, and outputs non shifted
% data, which is true to the inputted data.
%
% Aprojected = U'*Ainput
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
% Project a N-dimensional vector into M-modes. The M-modes define a basis
% set of vectors. Each mode has N-components. After the projection is made
% the k-th mode yields a generalized coordinate; qk(t) = <mode-k|A(t)>. 
% Here |A(t)> is one column vector from the A-data-matrix, where it has 
% N-rows and D columns such that t=1 to D. As a matrix, size(A) = NxD, and
% size(q) = MxD, where M is .LE. N. U could be complete, in which case 
% N = M, but in general U could define a small number of basis vectors,
% and this subspace is of interest.
%
% OUTPUT  dataMatrixInfo. <-- data structure
% The same type of data structure is outputed, with appropriate changes 
% made in the meta-data description. For example, pType is modified to 
% indicate the generalized coordinates represent notVectored data.
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
%%                                                       do the projection
   for k=1:N
   A{k} = U'*A{k};       % expresses subspace in terms of U-basis vectors
   end
%%                                                      package dataOutput
dataOutput = struct;
dataOutput.dataRefName = fname;
dataOutput.dataMatrixName = aMname;
dataOutput.nVariables = dataInput.nVariables;
dataOutput.pType = 'notVectored';      % info about xyz-components is lost
dataOutput.dim = 1;
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

