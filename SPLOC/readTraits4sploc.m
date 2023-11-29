function [traitData,T] = readTraits4sploc(dataRefName, ...
                                       cFname,mFname,nSamples,mFormat)
% read list of file names for cMatrices & mVectors, then read the matrices
%
% INPUT
% dataRefName   = name for the bundle of data as a required reference.
% cFname        = file name cell array for symmetric correlation matrices
% mFname        = file name cell array for the corresponding mean vectors
% nSamples      = number of samples used to calculate all matrices/vectors
% mFormat       = data structure:  pType and dim
% mFormat.pType = {xyz-xyz-xyz,xxx-yyy-zzz,notFormated}
% mFormat.dim   = number of components defining local vector property
% 
% USAGE
% readTraits4sploc(dataRefName,cFname,mFname,nSamples)
% readTraits4sploc(dataRefName,cFname,mFname,nSamples,mFormat)
%
% Remark: The format for the correlation matrices and mean vectors cannot
% be determined from the data. The user must know the input format.
% 
% PROCESS
% Determine if input is a single string or an cell array of strings. 
% define a cell array of 1 element or use N elements
% read the data matrices using the assumed format or specified format. 
%
% OUTPUT: traitData <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%  nVariables = number of variables defining size of vector space
%       pType = packing type (format) that is applied to the vector space
%         dim = # of components in a local vector (e.g. x,y,z for 1 atom)
%         nis = number of independent subsamples in this collection
%          mu = cell array for mean vectors representing the system
%          cM = cell array for cMatrix (covariance, correlation, etc) 
%  sampleSize = sample size for partitioning data matrices (a fixed value)
%     nDtotal = total # of data samples = nis*sampleSize
%%                                            associate with SPLOC toolset
global gvSPLOC
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                           process input
   if( nargin < 4 )
   error('missing required information: See Usage');
   elseif( nargin == 4 )
   mFormat = setDataMatrixFormat;                 % => launches user query
   elseif( nargin == 5 )
   mFormat = setDataMatrixFormat(mFormat);      % => checks mFormat validy
   end 
%%                                                         make cell array
   if( iscell(cFname) )                   % => multiple files will be read
   nis = length(cFname);
   else
   nis = 1;
   temp = cFname;
   cFname = cell(1,nis);
   cFname{1} = temp;
   temp = mFname;
   mFname = cell(1,nis);
   mFname{1} = temp;
   end
% REMARK: It is assumed that cFname and mFname come in associated pairs.
%%                               read data matrices within input directory
mu = cell(1,nis);
cM = cell(1,nis);
ref_j = -1;
   for j=1:nis
   cFn = getInputFileName('input',cFname{j}); 
   mFn = getInputFileName('input',mFname{j});
   cM{j} = importdata(cFn);
   mu{j} = importdata(mFn);
   [a,b] = size(cM{j});
   [c,d] = size(mu{j});
   c = max(c,d);
      if( a ~= b )
      error('correlation matrix is not square');
      end
      if( a ~= c )
      error('mean vector is not associated with correlation matrix');
      end
      if( ref_j < 0 )
      ref_j = a;
      else
          if( a ~= ref_j )
          error(['all matrices/vectors do not share same ',
                 'vector space dimension']);
          end
      end
   end
%%                                                       consistency check
jj = floor(ref_j/mFormat.dim)*mFormat.dim;
   if( ref_j ~= jj )
   error(['specified matrix format for packing is ', ...
          'inconsistent with matrix size']);
   end
%%                                      create/build output data structure
nDtotal = nis*nSamples;
traitData = struct;
traitData.dataRefName = dataRefName;               % file name as a string
traitData.mMatrixName = mFname;                    % cell array of strings
traitData.cMatrixName = cFname;                    % cell array of strings
traitData.nVariables = ref_j;
traitData.pType = mFormat.pType;                   % packing type (format)
traitData.dim = mFormat.dim;             % # of components in local vector
traitData.n = nis;                          % number of cM and mu matrices
traitData.mu = mu;                                % cell array: mu => mean
traitData.cM = cM;              % cell array: covariance, correlation, etc
traitData.sampleSize = nSamples;                       % applies uniformly
traitData.nDtotal = nDtotal;          % total sample size = nis*sampleSize
%%                                                            create table
rSamples = nSamples*ones(1,nis);
nVars = ref_j*ones(1,nis);
k = 1:nis;
T = table(k',nVars',rSamples',cFname',mFname','VariableNames', ...
          {dataRefName,'nVariables','nSamples', ...
          'cMatrixName','mVectorName'});
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['    reference name = ',dataRefName]);
fprintf(fid,'%s \n',['matrix sample size = ',num2str(nSamples)]);
fprintf(fid,'%s \n',['matrix format type = ',traitData.pType]);
  if( traitData.dim > 1 )
  fprintf(fid,'%s \n',['   # of components = ',num2str(traitData.dim)]);
  end
fprintf(fid,'%s \n',['  # of cM matrices = ',num2str(nis)]);
fprintf(fid,'%s \n',[' total sample size = ',num2str(nDtotal)]);
fprintf(fid,'%s \n',['    # of variables = ',num2str(ref_j)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end

