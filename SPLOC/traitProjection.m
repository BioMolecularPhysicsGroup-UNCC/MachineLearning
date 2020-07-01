function traitOut = traitProjection(fname,traitIn,U)
% projects statistical C-matrix into subspace defined by U basis vectors.
% Given a DxD C-matrix of some type, it will be projected into a subspace
% of Ds dimension. Note that D= total # of degrees of freedom (DOF) in the
% system, and Ds = # of DOF in the subspace of interest defined by basis
% vectors contained in the U matrix. size(U) = DxDs. 
% Then, traitOut will be a DsxDs matrix where the variables correspond to 
% the collective variables associated with the basis vectors.
% Css = C-matrix in the subspace 
% Css = U' * C * U     because the basis vectors must form an orthornormal
% set of vectors, and as such U'*U = identity matrix in a DsxDs subspace.
% Note that U' is the inverse of U for full rank matrices.  
% ------------------------------------------------------------------------
%
% DEFINITIONS
% traitData.  <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%  nVariables = number of variables defining size of vector space
%       pType = packing type (format) that is applied to the vector space
%         dim = # of components in a local vector (e.g. x,y,z for 1 atom)
%           n = number of independent subsamples in this collection
%          mu = cell array for mean vectors representing the system
%          cM = cell array for cMatrix (covariance, correlation, etc) 
%  sampleSize = sample size for partitioning data matrices (a fixed value)
%     nDtotal = total # of data samples = nis*sampleSize
%
% -------------------------------------------------------- local variables 
% INPUT
% fname = file name to specify origin of data and/or identify the output.
% traitIn = traitData for one or more systems of interest. It would
%             be particularly convenient to group all systems together 
%             that have a similar property. 
% U = basis vectors employed to project data into subspace.  
%
% PROCESS
% Project a DxD matrix into a DsxDs matrix describing couplings between
% the basis vectors that define the subspace.
%
% OUTPUT  traitData.  <-- data structure
% The same type of data structure is outputed, with appropriate changes 
% made in the meta-data description. For example, pType is modified to 
% indicate the generalized coordinates represent notVectored data.
%
%%                                                    set sploc parameters
global gvSPLOC
% --------------------------------------------------- for recording action
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                 unpackage data matrices
N = traitIn.n;                                      % # of unknown systems
cM = traitIn.cM;                                % cell array of c-matrices
mu = traitIn.mu;       % cell array of mean coordinates that cM relates to
mMatrixName = traitIn.mMatrixName;       % cell array of mean vector names
cMatrixName = traitIn.cMatrixName;          % cell array of c-matrix names
%%                                                      simple error check
[nDOF,temp] = size(cM{1});      % # of degrees of freedom = # of variables
[nV,nModes] = size(U);         % nModes = # of basis vectors to be checked
   if( nDOF ~= temp )
   error([' cM{1} matrix is not square: size = [', ... 
          num2str([nDOF,temp]),']']);
   end
   if( nDOF ~= nV )
   disp([' cM:  nDOF = ',num2str(nDOF)]);
   disp([' U:   nDOF = ',num2str(nV)]);
   disp([' U: nmodes = ',num2str(nModes)]);
   error('number of DOF differ between cM and U ');
   end
%%                                                       do the projection
   for k=1:N
   mu{k} = U'*mu{k};      % expresses subspace in terms of U-basis vectors
   cM{k} = U'*cM{k}*U;
   end
%%                                                      package traitOutput
traitOut = struct;
traitOut.dataRefName = fname;
traitOut.mMatrixName = mMatrixName;
traitOut.cMatrixName = cMatrixName;
traitOut.nVariables = nModes;
traitOut.pType = 'notVectored';      % info about xyz-components is lost
traitOut.dim = 1;
traitOut.n = N;
   for k=1:N
   traitOut.mu{k} = mu{k};
   traitOut.cM{k} = cM{k};
   end
traitOut.sampleSize = traitIn.sampleSize;
traitOut.nDtotal = traitIn.nDtotal;
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
str = traitIn.dataRefName;
fprintf(fid,'%s \n',['     input trait name = ',str]);
fprintf(fid,'%s \n',['    output trait name = ',num2str(fname)]);
fprintf(fid,'%s \n',[' dimension of input C-matrices = ',num2str(nDOF)]);
fprintf(fid,'%s \n',['dimension of output C-matrices = ',num2str(nModes)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end

