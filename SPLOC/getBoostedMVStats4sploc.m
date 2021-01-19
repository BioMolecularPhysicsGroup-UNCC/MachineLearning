function traitData = getBoostedMVStats4sploc(dataMatrixInfo, ...
                     nPairs,verbosity,cType)
% convert data matrices to various multivariate statistical moments
%
% INPUT
% datastructure: dataMatrixInfo
%                dataMatrixInfo.n = number of data matrices
%                dataMatrixInfo.A{:} = data matrix indexed from 1 to n
%                dataMatrixInfo.dataMatrixName = cell array of file names
%                dataMatrixInfo.nSamples = # of samples per data matrix
%                dataMatrixInfo.pType = string defines packing type
%                dataMatrixInfo.dim = # of components in local vector
% nPairs = # of pairs of half-samples. 1 sample = 2*(half-samples)
%          input sample is divided up into nPairs # of random half-samples
% verbosity: Specifies amount of intermediate processing steps to report.
% default: 0 => minimal checklist report in sploc log file.
% process: 1 => same as 0, plus detailed report in specialized log file.
% summary: 2 => same as 1, plus writing figures to show key relationships.
% display: 3 => same as 2, with figures displayed on the screen, sometimes
%               with a pause, and/or with additional output printed to the 
%               command window. Useful to understand how the code works!
%       ---> All printed figures have the same file type (fig, png, etc).
%       ---> Specialized log files are written in separate directories.
% ------------------------------------------------------------------------
% cType = {'cov','cor'} for covariance or correlation matrix (default=cov)
%
% USAGE:
%       traitCovN = getBoostedMVStats4sploc(AmatrixInfoN,1);
%       traitCovF = getBoostedMVStats4sploc(AmatrixInfoF,2,2);
%       traitCovN = getBoostedMVStats4sploc(AmatrixInfoN,5,0,'cov');
%       traitCorN = getBoostedMVStats4sploc(AmatrixInfoN,5,1,'cor');
%
% REMARK: expect   size(data matrix) = (# variables) x (# of data samples)
%    nV = # of variables
%    nD = # of data samples
%
% PROCESS
% split data matrices into independent windows specified by sampleSize
% As such AF{k} --> nis independent sample sets.
% => nF --> nis*nF0    and    nN --> nis*nN0
% General consistency and error checks will be made.
%
% OUTPUT: traitData <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%  nVariables = number of variables defining size of vector space
%       pType = packing type (format) that is applied to the vector space
%         dim = # of components in a local vector (e.g. x,y,z for 1 atom)
%         nis = number of instituted subsamples in this collection
%          mu = cell array for mean vectors representing the system
%          cM = cell array for cMatrix (covariance, correlation, etc) 
%  sampleSize = sample size for partitioning data matrices (a fixed value)
%     nDtotal = total # of data samples = nis*sampleSize
%%                                            associate with SPLOC toolset
global gvSPLOC                         % shared information across toolset
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                           process input
dataRefName = dataMatrixInfo.dataRefName;  % reference name stays the same
 n = dataMatrixInfo.n;                      % # of different data matrices
nD = dataMatrixInfo.nSamples;          % # of data samples per data matrix
nV = dataMatrixInfo.nVariables;    % # of variables is the same for 1 to n
% ------------------------------------------ data requirement error checks
  if( nV < 2 )
  error('nV < 2  Must have at least two variables');
  end
nDmax = max(nD);
nDmin = min(nD);
   if( nDmax ~= nDmin )
   msg1 = ['nDmax = ',num2str(nDmax),'  nDmin = ',num2str(nDmin)];
   msg2 = 'nDmax must equal nDmin to employ boosting';
   % This condition could be relaxed some, but the programming is easier
   % if all systems have the same sample size. It is best practice to make
   % this condition true, so this restriction is put here. If a user wants
   % to override sound rationale protocols, he/she can simply take care of
   % partitioning uneven batches of data into equal size batches of data
   % in their own programming for preprocessing the data streams. 
   disp(msg1);
   error(msg2);
   end
sampleSize = nDmax;
% ------------------------
   if( sampleSize < 2*nV )
   msg1 = ['# samples = ',num2str(sampleSize), ...
           '  # variables = ',num2str(nV)];
   msg2 = '# of samples must equal at least twice # of variables';
%  This is done so that after splitting the data into two halfs, each half
%  is a full rank matrix. While this is not essential, this is a minimum
%  requirement. OPV >> 1 is best. OPV = 1 is bad, which is set here as a
%  minimum requirement. If a user wants to override rationale protocols, 
%  he/she can simply take care of this by duplicating samples outside of
%  this code. This is highly not advisable though. 
   disp(msg1);
   error(msg2);
   end
% -------------------------------- check for input errors and set defaults
   if( nargin < 1 )
   msg = ['Must specify dataMatrixInfo data-structure ', ...
          'from readDataMatrices().'];
   error(msg);
   elseif( nargin == 1 )
   msg = 'Must specify number of pairs of instituded half-samples';
   error(msg);
   elseif( nargin > 1 )
      if( nPairs < 1 )
      error('nPairs must be 1 or greater');
      end
% --------------------------- work with other variables besides sampleSize
      if( nargin == 2 )
      verbosity = 0;                                             % default
      cType = 'cov';           % default value specifies covariance matrix
      elseif( nargin == 3 )
      verbosity = setVerbosity(verbosity);  % specifies options: {0,1,2,3}
      cType = 'cov';           % default value specifies covariance matrix
      else
      verbosity = setVerbosity(verbosity);  % specifies options: {0,1,2,3}
      cType = lower(cType);
         switch cType
            case {'cov','cor'}
            % do nothing
             otherwise
            error('unrecognized second moment matrix: Try cov or cor.');
         end
      end
   end
A = dataMatrixInfo.A;                    % cell array for each data matrix
%%                                                            process data
sLast = round(sampleSize/2);          % sample size of 1 covariance matrix
prefix_cM = [cType,num2str(sLast),'_'];
prefix_mu = ['ave',num2str(sLast),'_'];
dataMatrixName = dataMatrixInfo.dataMatrixName;    % cell array of strings
   for k=1:n
      if( dataMatrixName{k}(end-4) == '.' )
      dataMatrixName{k} = dataMatrixName{k}(1:end-4);
      end
   end
nis = n*(2*nPairs);
nDtotal = nPairs*sampleSize;
cMatrixName = cell(1,nis);
mMatrixName = cell(1,nis);
mu = cell(1,nis);
cM = cell(1,nis);
sStart = sLast + 1;
ns1 = sLast;                             % # of samples in 1st half sample
ns2 = sampleSize - ns1;                  % # of samples in 2nd half sample
padd = '%03i';
ndigits = min(1 + floor( log10(0.0001 + nPairs) ),9);
padd(3) = num2str(ndigits);
m = 0;                                     % counter ranging from 1 to nis
   for k=1:n
   [mV,mD] = size(A{k});              % as already checked mD = sampleSize
      if( mV ~= nV )
      error('number of variables across data matrices are not equal');
      end
      if( mD ~= sampleSize )
      error('number of samples across data matrices are not equal');
      end
   B = A{k};
   j = 0;                             % counter ranging from 1 to 2*nPairs
      for i=1:nPairs
      B = B(:, randperm(sampleSize) );                   % shuffle samples
% ------------------------------------------------------- 1st half dataset
      j = j + 1;
      B1 = B(:,1:sLast);
      m = m + 1;                    % odd values of m => 1st half datasets
      mu{m} = mean(B1,2);
      B1 = B1 - mu{m};                        % row center the data matrix
      cM{m} = (B1*B1')/ns1;                            % covariance matrix
      fname = [prefix_cM,num2str(j,padd),'_',dataMatrixName{k}];
      cMatrixName{m} = fname;
      fname = [prefix_mu,num2str(j,padd),'_',dataMatrixName{k}];
      mMatrixName{m} = fname;
% ------------------------------------------------------- 2nd half dataset
      j = j + 1;
      B2 = B(:,sStart:sampleSize);
      m = m + 1;                   % even values of m => 2nd half datasets
      mu{m} = mean(B2,2);
      B2 = B2 - mu{m};                        % row center the data matrix
      cM{m} = (B2*B2')/ns2;                            % covariance matrix
      fname = [prefix_cM,num2str(j,padd),'_',dataMatrixName{k}];
      cMatrixName{m} = fname;
      fname = [prefix_mu,num2str(j,padd),'_',dataMatrixName{k}];
      mMatrixName{m} = fname;
      end  
   end
%%                                convert covariance to correlation matrix
   if( strcmp(cType,'cor') == 1 )
      for m=1:nis
      cM{m} = corrcov(cM{m});             % use MATLAB conversion function
      % mu{m} remains the same.
      end
   end
%%                                             create/build data-structure
traitData = struct;
traitData.dataRefName = dataRefName;               % file name as a string
traitData.mMatrixName = mMatrixName;               % cell array of strings
traitData.cMatrixName = cMatrixName;               % cell array of strings
traitData.nVariables = nV;
traitData.pType = dataMatrixInfo.pType;            % packing type (format)
traitData.dim = dataMatrixInfo.dim;      % # of components in local vector
traitData.n = nis;                          % number of cM and mu matrices
traitData.mu = mu;                                % cell array: mu => mean
traitData.cM = cM;                 % cell array: covariance or correlation
traitData.sampleSize = sLast;                          % applies uniformly
traitData.nDtotal = nDtotal;          % total sample size = nis*sampleSize
%%                              write summary information to spoc log file
fid = fopen(splocLogFile,'a');
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('calculation summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['    reference name = ',dataRefName]);
fprintf(fid,'%s \n',['matrix sample size = ',num2str(sampleSize)]);
fprintf(fid,'%s \n',['matrix format type = ',traitData.pType]);
  if( traitData.dim > 1 )
  fprintf(fid,'%s \n',['   # of components = ',num2str(traitData.dim)]);
  end
fprintf(fid,'%s \n',['  # of cM matrices = ',num2str(nis)]);
fprintf(fid,'%s \n',[' total sample size = ',num2str(nDtotal)]);
fprintf(fid,'%s \n',['    # of variables = ',num2str(nV)]);
fprintf(fid,'%s \n',['           cM type = ',cType]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
%%                     write additional information depending on verbosity
gFileType = gvSPLOC.gFileType;
   if( verbosity > 0 ) 
   subFolder = 'input';
   cMfname = cell(1,nis);
   mMfname = cell(1,nis);
      for m=1:nis
      cMfname{m} = getOutputFileName(subFolder,cMatrixName{m});
      mMfname{m} = getOutputFileName(subFolder,mMatrixName{m});
      cMresults = [cMfname{m},'.dlm']; 
      mMresults = [mMfname{m},'.dlm']; 
%       disp([cMatrixName{m},'   ',cMresults]);               % checks out
%       disp([mMatrixName{m},'   ',mMresults]);               % checks out
%       disp('  ');
      dlmwrite(cMresults,cM{m},'precision',10);   % need at least 9 digits
      dlmwrite(mMresults,mu{m},'precision',10);   % need at least 9 digits
      end
% ---------------------------------------- plot mean vectors on same graph
      if( verbosity == 3 )
      figure;
      figure_number = get(gcf,'Number');
      figMean = figure(figure_number);
      else
      figMean = figure('visible','off');
      end
   clf
   nComponents = 1:nV;
   plot(nComponents',mu{1});
   hold on;
      for k=2:n
      plot(nComponents',mu{k});
      end
   xlabel('number of components');
   ylabel('mean row centered vector');
   title('collection of mean row centered vectors');
   mMf = [prefix_mu,dataRefName];
   mMf = getOutputFileName(subFolder,mMf);
   saveas(figMean,mMf,gvSPLOC.gFileType);
   end
% ------------------------------------- consider writing graphics to files
   if( verbosity == 3 )
   colorMatrixTool(cM,0,'fName',cMfname,'fType',gFileType, ...
                   'xLabel','coordinate index', ...
                   'yLabel','coordinate index', ...
                   'cType','b0r','cutupper',95,'cutlower',95, ...
                   'morezero',5);
   elseif( verbosity == 2 )
   colorMatrixTool(cM,0,'fName',cMfname,'fType',gFileType, ...
                   'xLabel','coordinate index', ...
                   'yLabel','coordinate index', ...
                   'cType','b0r','cutupper',95,'cutlower',95, ...
                   'morezero',5,'fShow','no');
   end
end
