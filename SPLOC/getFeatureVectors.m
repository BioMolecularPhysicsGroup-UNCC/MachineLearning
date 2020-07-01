function featureX = getFeatureVectors(traitX,U)
% calculates two features per basis vector given in U per input system X
% The features are the mean and standard deviation per selection basis 
% vector. A system of M dimensions gets mapped into feature space of 2M
% dimensions. However, the dimension of the discriminant subspace, within 
% feature space will generally be much smaller than number of DOF (nDOF)
% in the system. 
%
% DEFINITIONS
% ------------------------------------------------------- trait definition
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
% -------------------------------------- define local language for sploc()
% binary states: 1 => "on" => Function    and    0 => "off" => Nonfunction
% M => mean column vector
% Q => symmetric matrix, such as a covariance matrix.
%
% -------------------------------------------------------- local variables
% INPUT
% traitX <---- data structure (see above for definition)
% U = basis vectors that are to be used for the classification process.
%     Note that U can be catenated from an ensemble of viable sets of U.
%     In the case that U represents multiple U sets, the vectors contained
%     within U (ie the modes) will not satisfy all-pairwise orthogonality.
%     This redundancy can help reduce noise in the classification process.
% ------------------------------------------------------------------------
%
% PROCESS
%+++ Per basis vector (or per mode):
% 1) calculate projected variance and mean
% 2) pack mean and then variance in back to back components
% 3) cycle through all modes listed in U
% 4) output a feature vector having 2xD dimenions where D is the number
%    of modes contained in U.
% 5) A feature vector is created per system considered.
%
% OUTPUT
% featureX. <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%   nXsystems = number of systems being projected into feature space
%               X represents any system (0,1,unknown) that is ranked
%                 1-system => from training set labeled as functional
%                 0-system => from training set labeled as nonfunctional
%      nModes = # of discriminant modes contained in U.
%   nFeatures = 2xnModes = number of distinct features
%     Fmatrix = nFeatures x nXsystems  
% ---------------------------------------------------------- basis vectors
%            U = nV x nModes matrix containing nModes used to discriminate
%                data, where each mode lives in a nV dimensional space. 
%%                                                    set sploc parameters
global gvSPLOC
% -------------------------------------------- scoring function parameters
add0 = gvSPLOC.add0;                 % prevents variance ratios to diverge
% --------------------------------------------------- for recording action
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                               check input specification
   if( nargin ~= 2 )
   error('input format/data is required to be: traitX,U');
   end
%%                                                        unpackage traits
Nx = traitX.n;
Mx = traitX.mu;
Qx = traitX.cM;
dataRefName = traitX.dataRefName; 
mMname = traitX.mMatrixName;
cMname = traitX.cMatrixName;                
%%                                                      simple error check
[nDOF,~] = size(Qx{1});  % nDOF = # of degrees of freedom = # of variables
[nV,nModes] = size(U);         % nModes = # of basis vectors to be checked
   if( nDOF ~= nV )
   disp(['Qx:   nDOF = ',num2str(nDOF)]);
   disp([' U:   nDOF = ',num2str(nV)]);
   disp([' U: nmodes = ',num2str(nModes)]);
   error('number of variables differ between Qx and basis set dimension');
   end
%%                                      set verbosity output level details
% if( verbosity > 0 )
 subFolder = 'classification';
 fName = [dataRefName,'_feature.log']; 
 %               ^^^^^^^^^^^--------> user controls uniqueness in basename
 iLogFileName = getOutputFileName(subFolder,fName);
%  disp(iLogFileName);                                     % for debugging
   if( isfile(iLogFileName) )
   fid = fopen(iLogFileName,'a');     % => do not need header: append data
   else                                % => virgin file: needs header info
   fid = fopen(iLogFileName,'w');          % => will overwrite old results
   msg ='layout for how feature vectors are packed for classification';
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n','   ');
   msg ='     Uk = basis vector in original space';
   fprintf(fid,'%s \n',msg);
   msg ='    MUk = average component along the Uk projection.';
   fprintf(fid,'%s \n',msg);
   msg ='    SDk = standard deviation along the Uk projection.';
   fprintf(fid,'%s \n',msg);
   msg ='   Fvec = [MU1,SD1, MU2,SD2, MU3,SD3, ... MUn,SDn]';
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n',dividerLine);
   end
% FIX ME write some useful information here
% end
%%                                         pre-calculate projected moments
aveX = zeros(Nx,nModes);
varX = zeros(Nx,nModes);
   for k=1:nModes                    % consider each k-th vector-direction
   vec = U(:,k);                                     % the selected vector
% ------------------------------------ project into unknown system moments
      for jx=1:Nx                            % all unknown X-state systems
      aveX(jx,k) = vec'*Mx{jx};                       % projected  average
      varX(jx,k) = max(vec'*Qx{jx}*vec,add0);         % projected variance
      end
   end
%%                               pack this information as an output vector
nComponents = 2*nModes;
Fmatrix = zeros(nComponents,Nx);
   for jx=1:Nx
   jj = 0;
      for k=1:nModes
      jj = jj + 1;
      Fmatrix(jj,jx) = aveX(jx,k);                                  % mean
      jj = jj + 1;        
      Fmatrix(jj,jx) = sqrt( varX(jx,k) );                     % std. dev.
      end
   end
% --------------------------------- record information into data structure
featureX = struct;
featureX.dataRefName = dataRefName;
featureX.mMname = mMname;
featureX.cMname = cMname;
featureX.nXsystems = Nx;
featureX.nModes = nModes;
featureX.nFeatures = 2*nModes;
featureX.Fmatrix = Fmatrix;
%%                                         record action in sploc log file
nM = nModes;
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['           feature extraction: ']);
fprintf(fid,'%s \n',['  (from traitX) reference name = ',dataRefName]);
fprintf(fid,'%s \n',['number of X-systems classified = ',num2str(Nx)]);
fprintf(fid,'%s \n',['  number of discriminant modes = ',num2str(nM)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end
