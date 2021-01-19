function traitOut = traitFilter(fname,traitIn,U,reduce)
% filters out from the C-matrix statistical data outside the U-subspace.
% Given a DxD C-matrix of some type, it will be projected into a subspace
% of Ds dimension. However, the basis vectors that are used to define the
% U-subspace are expressed in terms of the original coordinates, so the 
% final C-matrix will be expressed as a DxD matrix. Mathematically, any
% statistical information related to variables outside the U-subpace are
% filtered out.  
% C_filtered = (U*U') * C * (U*U') 
% Note that if U was complete, no filtering would occur, because U*U' = 1. 
% Because the basis vectors form an orthornormal set of vectors, but span
% a subspace that is not complete, U*U' is a projection operator that can
% be obtained by   sum_k  |k><k|    for the set of k of interest. 
% Note that U' is the inverse of U for full rank matrices. 
% 
% when reduce = 'reduce' the matrices will be reduced from having xyz 
% components in vectorized form to a form that is not vectorized. 
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
% when reduce = 'reduce'
% Reduces a (dimxN) by (dimxN) matrix down to a NxN matrix by grouping 
% corresponding vectorized components. For instance, the packing in the 
% input Cmatrix can be of the form xyz-xyz-xyz or xxx-yyy-zzz. At the end,
% the reduced matrix is really of the form: xx + yy + zz per element. This 
% function generalizes the number of components in a vector from 3 to dim.
%
% OUTPUT  traitData.  <-- data structure
% The same type of data structure is outputed, with appropriate changes 
% made in the meta-data description. For example, pType & dim are modified
% to indicate the generalized coordinates represent notVectored data if
% vectorized input is reduced. 
%
%%                                                             parse input
   switch nargin
     case 3
     reduce = 'no';                  % default is not to reduce the matrix
     otherwise
        if( nargin < 3 )
        error('minimally must specify: {fname, traitIn, U}');
        end
   end
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
%%                                                        do the filtering
   for k=1:N
   P = U*U';                                    % DxD projection matrix
   mu{k} = P*mu{k};  % filter out parts of mu-vector outside of U-subspace
   cM{k} = P*cM{k}*P; % filter out parts of c-matrix outside of U-subspace
   end
%%                                                    reduce vectorization
pType = traitIn.pType;
dim = traitIn.dim;
% ------------------------------------------------------------ error check
   if( strcmp(pType,'notVectored') )
   reduce = 'no';                                          % non-reducible
   end
   if( strcmp(reduce,'reduce') )
   nDOFout = round(nDOF/dim);
   pType = 'notVectored';              % info about xyz-components is lost
   dim = 1;
      for k=1:N
      cM{k} = reduceCmatrix(cM{k},pType,dim);
      mu{k} = [];                                % must set to null vector
%     REMARK: It is not possible to reduce the mu vectors. This operation
% is not defined mathematically because reducing a C-matrix is based on
% inner products, and mu is a vector as well as its reduced form. As such,
% to prevent misuse of a meaningless output, the mu is set to null.
      nV = nDOFout;
      end
   else
   nDOFout = nDOF;                                             % no change
   nV = traitIn.nVariables;
   end
%%                                                     package traitOutput
traitOut = struct;
traitOut.dataRefName = fname;
traitOut.mMatrixName = mMatrixName;
traitOut.cMatrixName = cMatrixName;
traitOut.nVariables = nV;
traitOut.pType = pType; 
traitOut.dim = dim;
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
fprintf(fid,'%s \n',['              input trait name = ',str]);
fprintf(fid,'%s \n',['             output trait name = ',fname]);
fprintf(fid,'%s \n',[' dimension of input C-matrices = ',num2str(nDOF)]);
fprintf(fid,'%s \n',['number of modes used in filter = ', ...
                     num2str(nModes)]);
    if( strcmp(reduce,'reduce') )
    fprintf(fid,'%s \n','         reduced output matrix = true');
    else
    fprintf(fid,'%s \n','         reduced output matrix = false');
    end
fprintf(fid,'%s \n',['dimension of output C-matrices = ', ...
                     num2str(nDOFout)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end

