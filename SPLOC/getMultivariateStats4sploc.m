function traitData = getMultivariateStats4sploc(dataMatrixInfo, ...
                     sampleSize,verbosity,cType)
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
% sampleSize = # of samples used to calculate multivariate statistics
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
% traitCovF = getMultivariateStats4sploc(AmatrixInfoF);
% traitCovN = getMultivariateStats4sploc(AmatrixInfoN,2000);
% traitCovF = getMultivariateStats4sploc(AmatrixInfoF,2000,2);
% traitCovN = getMultivariateStats4sploc(AmatrixInfoN,5000,0,'cov');
% traitCorN = getMultivariateStats4sploc(AmatrixInfoN,5000,1,'cor');
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
%         nis = number of independent subsamples in this collection
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
% -------------------------------- check for input errors and set defaults
  if( nV < 2 )
  error('nV < 2  Must have at least two variables');
  end
   if( nargin < 1 )
   msg = ['Must specify dataMatrixInfo data-structure ', ...
          'from readDataMatrices().'];
   error(msg);
   elseif( nargin == 1 )
   nDmax = max(nD);
   sampleSize = nDmax;             % pick the best sampling space possible
   verbosity = 0;                                                % default
   cType = 'cov';              % default value specifies covariance matrix 
   elseif( nargin > 1 )
   nDmax = max(nD); 
      if( sampleSize < (nV + 1) )
      error('Input sampleSize must be at least (1 + # of variables).');
      elseif( sampleSize > nDmax )
      msg = ['Input sampleSize cannot exceed maximum number ', ...
             'of data samples.'];
      error(msg);
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
%%                                     filter out data that cannot be used
nis = 0;
keep = -ones(1,n);                                        % => no not keep
%nKeep = 0;
   for k=1:n
   i = floor( 0.0000000001 + nD(k)/sampleSize );
      if( i > 0 )
      nis = nis + i;
      keep(k) = i;     % => # of independent segments for this data matrix
%     nKeep = nKeep + 1;
      end
   end
A = dataMatrixInfo.A;                    % cell array for each data matrix
%%                                                            process data
prefix_cM = [cType,num2str(sampleSize),'_'];
prefix_mu = ['ave',num2str(sampleSize),'_'];
dataMatrixName = dataMatrixInfo.dataMatrixName;    % cell array of strings
   for k=1:n
      if( dataMatrixName{k}(end-4) == '.' )
      dataMatrixName{k} = dataMatrixName{k}(1:end-4);
      end
   end
cMatrixName = cell(1,nis);
mMatrixName = cell(1,nis);
jDmin = zeros(1,nis);
jDmax = zeros(1,nis);
mu = cell(1,nis);
cM = cell(1,nis);
nDtotal = 0;
m = 0;
   for k=1:n
      if( keep(k) > 0 )
      [mV,mD] = size(A{k});
         if( mV ~= nV )
         error('number of variables across data matrices are not equal');
         end
      i = floor( 0.0000000001 + nD(k)/sampleSize );
      nDtotal = nDtotal + mD;
      jDmin(i) = mD - sampleSize + 1;
      jDmax(i) = mD;
         for j=i-1:-1:1
         jDmax(j) = jDmin(j+1) - 1;
         jDmin(j) = jDmax(j) - sampleSize + 1;
         end
% ----------------------------------------------------- work on file names
      padd = '%03i';
      ndigits = min(1 + floor( log10(0.0001 + i) ),9);
      padd(3) = num2str(ndigits);
         for j=1:i
         m = m + 1;
         s1 = jDmin(j);
         s2 = jDmax(j);
         B = A{k}(:,s1:s2);
         mu{m} = mean(B,2);
         B = B - mu{m};                       % row center the data matrix
         cM{m} = (B*B')/(s2 - s1);                     % covariance matrix
         %disp( size(cM{m}) );
         %disp( (s2 - s1 + 1) );
         fname = [prefix_cM,num2str(j,padd),'_',dataMatrixName{k}];
         cMatrixName{m} = fname;
         fname = [prefix_mu,num2str(j,padd),'_',dataMatrixName{k}];
         mMatrixName{m} = fname;
         end 
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
traitData.sampleSize = sampleSize;                     % applies uniformly
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
