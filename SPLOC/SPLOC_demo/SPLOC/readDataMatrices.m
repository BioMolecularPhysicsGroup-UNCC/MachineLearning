function [dataMatrixInfo,T] = readDataMatrices(dataRefName, ...
                                        dataMatrixName,mFormat,fType)
% read list of file names for data matrices, then read those data matrices
%
% INPUT
% dataRefName    = name for the bundle of data as a required reference.
% dataMatrixName = file name for the data matrix to be read, or a cell
%                  array for list of file names to read each data matrix.
% mFormat        = data structure:  pType and dim
% mFormat.pType  = {xyz-xyz-xyz,xxx-yyy-zzz,notFormated}
% mFormat.dim    = number of components defining local vector property
% fType          = optional input with two allowed values (rowVar,colVar).
%         rowVar => rows => variables and columns => samples     (default)
%         colVar => columns => variables and rows => samples
% 
% USAGE
% readDataMatrices(dataRefName,dataMatrixName)
% readDataMatrices(dataRefName,dataMatrixName,mFormat)
% readDataMatrices(dataRefName,dataMatrixName,mFormat,'rowVar')  (default)
% readDataMatrices(dataRefName,dataMatrixName,mFormat,'colVar')
%
% Remark: The output A matrix has rows as variables. No actual check can 
% be made. The user must know the input format.
% 
% PROCESS
% Determine if input is a single string or an cell array of strings. 
% define a cell array of 1 element or use N elements
% read the data matrix using the assumed format or specified format. 
% convert data matrix format to the expected rowVar format
%
% OUTPUT dataMatrixInfo <-- data structure
% dataRefName    = reference name for A matrix data with similar traits
% dataMatrixName = cell array for file names that store the A matrix data 
%     nVariables = number of variables defining size of vector space
%          pType = packing type (format) applied to the vector space
%            dim = # of components in local vector (e.g. x,y,z for 1 atom)
%              n = number of data matrices
%           A{:} = cell array for the data matrices in the rowVar format
%    nSamples(:) = array for the number of samples in each data matrix
%%                                            associate with SPLOC toolset
global gvSPLOC
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                           process input
   if( nargin < 1 )
   error('missing data reference name and all data matrix file names!');
   elseif( nargin == 1 )
   error('missing all data matrix file names!');
   elseif( nargin == 2 )
   mFormat = setDataMatrixFormat;                 % => launches user query
   fType = 'rowVar';
   elseif( nargin == 3 )
   mFormat = setDataMatrixFormat(mFormat);      % => checks mFormat validy
   fType = 'rowVar';
   end
%%                                                      simple error check
flagError = 1;                                  % assume there is an error
   if( strcmp(fType,'rowVar') == 0 )
   flagError = -1;
   elseif( strcmp(fType,'colVar') == 0 )
   flagError = -1;
   end
   if( flagError > 0 )                      % => stop if there is an error
   error('2nd argument for format type can only be: rowVar  or  colVar');
   end  
%%                                                         make cell array
   if( iscell(dataMatrixName) )          % => multiple files will be read
   nM = length(dataMatrixName);
   else
   nM = 1;
   temp = dataMatrixName;
   dataMatrixName = cell(1,nM);
   dataMatrixName{1} = temp;
   end
% --------------------------------------------------- check for null input
   if( nM == 0 )
   error('dataMatrixName. List of unknown data matrices is empty!');
   end   
%%                               read data matrices within input directory
nSamples = zeros(1,nM);
nVar = zeros(1,nM);
A = cell(1,nM);
   if( strcmp(fType,'rowVar') == 1 )
      for j=1:nM
      fName = getInputFileName('input',dataMatrixName{j});
      A{j} = importdata(fName);
      [mV,mS] = size(A{j});
% ---------------------------------------- verify statistical significance
         if( mS < (mV+1) )
         error('# of samples is less than (1 + # of variables)');
         else
         nVar(j) = mV;
         nSamples(j) = mS;
         end
% ---------------------------------------- verify uniformality in # of DOF
         if( j == 1 )
         ref_j = nVar(j);
         else
            if( ref_j ~= nVar(j) )
            error('# of variables is not the same for all data matrices');
            end
         end
      end
   else                     % format of type: colVar, need transpose of it
      for j=1:nM
      fName = getInputFileName('input',dataMatrixName{j});
      B = importdata(fName);
      A{j} = B';
      [mV,mS] = size(A{j});
% ---------------------------------------- verify statistical significance
         if( mS < (mV+1) )
         error('# of samples is less than (1 + # of variables)');
         else
         nVar(j) = mV;
         nSamples(j) = mS;
         end
% ---------------------------------------- verify uniformality in # of DOF           
         if( j == 1 )
         ref_j = nVar(j);
         else
            if( ref_j ~= nVar(j) )
            error('# of variables is not the same for all data matrices');
            end
         end
      end
   end
%%                                                       consistency check
jj = floor(nVar(j)/mFormat.dim)*mFormat.dim;          % ref_j <--- nVar(j)
   if( nVar(j) ~= jj )
   error(['specified matrix format for packing is ', ...
          'inconsistent with matrix size']);
   end
%%                                      create/build output data structure
dataMatrixInfo = struct;
dataMatrixInfo.dataRefName = dataRefName;                         % string
dataMatrixInfo.dataMatrixName = dataMatrixName;    % cell array of strings
dataMatrixInfo.nVariables = ref_j;                    % numerical variable
dataMatrixInfo.pType = mFormat.pType;        % string defines packing type
dataMatrixInfo.dim = mFormat.dim;        % # of components in local vector
dataMatrixInfo.n = nM;                                % number of matrices
dataMatrixInfo.A = A;                                   % numerical matrix
dataMatrixInfo.nSamples = nSamples;                      % numerical array
%%                                                            create table
k = 1:nM;
nVars = ref_j*ones(1,nM); 
T = table(k',nVars',nSamples',dataMatrixName','VariableNames', ...
          {dataRefName,'nVariables','nSamples','dataMatrixName'});
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['    reference name = ',dataRefName]);
fprintf(fid,'%s \n',['matrix format type = ',mFormat.pType]);
  if( mFormat.dim > 1 )
  fprintf(fid,'%s \n',['   # of components = ',num2str(mFormat.dim)]);
  end
fprintf(fid,'%s \n',['    number of rows = ',num2str(ref_j)]);
m = min(nSamples);
fprintf(fid,'%s \n',['  min # of columns = ',num2str(m)]);
m = max(nSamples);
fprintf(fid,'%s \n',['  max # of columns = ',num2str(m)]);
fprintf(fid,'%s \n',['   # of A matrices = ',num2str(nM)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end

