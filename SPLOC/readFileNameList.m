function [nFiles,varargout] = readFileNameList(inputFname,varargin)
% read an input file specifying a list of file names with identifiers
%
% INPUT
% inputFname = name of input file
%
%   readInputFnameSet(inputFname)
%   readInputFnameSet(inputFname,verbosity)
%   readInputFnameSet(inputFname,verbosity,fullPath)
%   readInputFnameSet(inputFname,'command','value',...)
%   readInputFnameSet(inputFname,verbosity,'command','value',...)
%   readInputFnameSet(inputFname,verbosity,fullPath,'command','value',...)
% 
% commands      purpose/options
% 'IDcol'       to specify the column used for the classification ID
%               {0,1,2,...}
%               0 => no column is used as a classifier
%               1 => default value when # of columns > 1 in input file.
%               0 => default value when # of columns = 0 in input file.
% 'cType'       overrides classification type on the specified ID column.
%               when IDcol = 0, a new classification column is augmented. 
%               {F,U,N} are the allowed classifications. 
%               F => functional, U => unknown, N => nonfunctional
% 'sType'       selects output to only comprise of the class ID specified.
%               {F,U,N} are the allowed classifications.
%               
% verbosity: 0 => no write to screen
%            1,2 => write output to screnn
%            3 => same as 1,2 plus some additional 
%  fullPath: 0 => look for inputFname file within the input directory
%            1 => assume inputFname contains the full path of interest.
%
% Remark: cType options can be capital or small case
% 
% Allowed file formats:
% ----------------------------------------------------------- unclassified
% <fileName1>  <fileName2>  <fileName3>  ...
% ---------------------------------------------------- expected classified
% ID <fileName1>  <fileName2>  <fileName3>  ... 
% --------------------------------------- example of unexpected classified
% <fileName1>  <fileName2>  ID  <fileName3>  ...      **=> not recommended
% -------------------------------------- most common intended applications
% ID <AmatrixFileName>                            classify by ID = {F,U,N}
% ID <QmatrixFileName> <MvectorFileName>          classify by ID = {F,U,N}
% <AmatrixFileName>                                      (file classifies)
% <QmatrixFileName> <MvectorFileName>                    (file classifies)
%
% PROCESS
% A file with 1 column  is assumed by default to be unclassified.
% A file with 2 columns is assumed by default to be classified datastream.
% A file with 3 columns is assumed by default to be classified (Q and M).
%
% Sometimes it is useful to override the classification that is specified
% in a file, and make all filenames to be functional or nonfunctional or 
% unknown. If a file has various lines classified as F or U or N, and the 
% user specifies cType = 'U' then all files read in will be treated as 
% unclassified. On the other hand, one might want to list all files as 
% either functional or nonfunctional or unknown in a single file. In this
% case, there is no point in specifying the classifer per line. The user
% simply reads the file with the classifer column missing, and specifies
% the classification. 
%
% Error checking. For the column spacified as the classification ID, a 
% check is made that the ID is either F, U or N. 
%
% OUTPUT (different number of cell arrays are possible)
% [nFiles, classifer, {verbatim columns}]
% cell array for each column read in represents {verbatim columns} 
% The order of information in the input file is preserved. The user takes
% responsibility for manipulating the column order as it is required.
%%                                            associate with SPLOC toolset
global gvSPLOC
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                  process input commands
p = inputParser;
% -------------------------------------------------------------- verbosity
validVerbosity = @(x) isnumeric(x) && isscalar(x) && (x > -1) && (x<4);
addOptional(p,'verbosity',0,validVerbosity)   % default prints out nothing
% REMARK: {0,1,2} => no write  {3} => write out information
% --------------------------------------------------------------- fullPath
validPath = @(x) isnumeric(x) && isscalar(x) && (x > -0.000000001) ...
                                             && (x <  1.000000001);
addOptional(p,'fullPath',0,validPath) % look for file in expected location 
% ------------------------------------------- override classification type
% cType = classification type (F,U,N).     % => allowed classification IDs
expectedCtype = {'F','f','U','u','N','n'};        % remove case dependence
addParameter(p,'cType','-', ...             % default is no classification 
             @(x) any(validatestring(x,expectedCtype)));
% ------------------------------------ make selection of classification ID
% sType = classification ID type (F,U,N).  % => allowed classification IDs
expectedCtype = {'F','f','U','u','N','n'};        % remove case dependence
addParameter(p,'sType','-', ...                  % default is no selection 
             @(x) any(validatestring(x,expectedCtype)));
% -------------------------------------------------- check ID class column
% IDcol = integer from 0, 1, 2, 3, ... 
validIDcol = @(x) isnumeric(x) && isscalar(x) && (x > -0.00000000001);
addParameter(p,'IDcol',-1,validIDcol)      % usual default is first column
%%                                                  parse input parameters
parse(p,varargin{:});
% ------------------------------------------------- create local variables
cType = p.Results.cType;
sType = p.Results.sType;
cIDnum = round(p.Results.IDcol);
fullPath = round(p.Results.fullPath);
verbosity = round(p.Results.verbosity);
%%                                                   read contents of file
   if( fullPath == 0 )
   inputFname = getInputFileName('input',inputFname);
   end
   if( verbosity ~= 0 )
   disp('    ');
      if( fullPath == 0 )
      disp( 'looking within directory: input');
      disp(['Searching for input file: ',inputFname,' ...']);
      else                      % => user must supply full path + filename
      disp(['Searching for input file: ',inputFname,' ...']);
      end
   disp('    ');
   end
fid = fopen(inputFname);
   if( fid < 0 )
   error(['file: ',inputFname,' not found!']);
   end
txtFile = textscan(fid,'%s','delimiter','\n');
fclose(fid);
tempStr = split( txtFile{1}(:) );
[nLines,nCol] = size(tempStr);
%%                                         set defaults and error checking
   if( (nLines == 0) || (nCol == 0) )
   error('No lines are detected in the input file!');
   end
   if( cIDnum > nCol )                             % => not enough columns
   error('classification column does not exist!');
   end
   if( nCol == 1 )                         % => nCol = 1 is a special case
      if( cIDnum == -1 )   % special flag => default operation is expected
      cIDnum = 0;  % new default value: no column gives the classification      
      elseif( cIDnum == 1 )              % => there is nothing to classify
      error('must classify something');     % => 1st column gives filename
      end
   else
   cIDnum = abs(cIDnum);      % resets default -1 to 1 as the real default
%                            % resets any other value to itself, such as 0
   end
cType = upper(cType);               % force cType to be in capital letters
sType = upper(sType);               % force sType to be in capital letters
   if( cType ~= '-' )
       if( sType ~= '-' )                      % => selection is specified
          if( sType ~= cType )
          error('empty set will result: sType and cType are not equal');
          end
       end
   end
nOutCol = nargout - 1;      % reserve one extra output variable for nFiles
% ----------------------------------------------------- for debugging only
   if( verbosity == 3 )
   disp(['   nCol = ',num2str(nCol)]);
   disp([' cIDnum = ',num2str(cIDnum)]);
   disp(['nOutCol = ',num2str(nOutCol)]);
   disp(['  cType = ',upper(cType)]);
   end
% ------------------------------------------ check for augmented classifer
   if( (cType ~= '-') && (cIDnum == 0) )       % => add a classifer column
   nCol = nCol + 1;                               % will augment column #1
      if( nOutCol ~= nCol )
      disp('   ');
      disp('classifer ID column is augmented to # of input columns');
      error('# of output columns is not equal to 1 + # of input columns');
      end
   classiferID = cell(nLines,1);             % create the augmented column
      for j=1:nLines
      classiferID{j} = cType;            % populate with degenerate values
      end
   outList = cell(nCol,nLines);
   outList{1} = classiferID;                        % the augmented column
      for j=2:nCol                            % the original input columns
      outList{j} = tempStr(:,j-1);           % output columns shifted by 1
      end
   cIDnum = 1;                         % the 1st column is a classifer now
   else
      if( nOutCol ~= nCol ) 
      error('# of input columns is not equal with # of output columns');
      end
   outList = cell(nCol,nLines);
      for j=1:nCol                            % the original input columns
      outList{j} = tempStr(:,j); 
      end
      if( cType ~= '-' )            % must override all classification IDs
      classiferID = cell(nLines,1);                % create the new column
         for j=1:nLines
         classiferID{j} = cType;         % populate with degenerate values
         end   
      outList{cIDnum} = classiferID;     % replace with requested override
      else
      end
   end
%%                                  print results to screen when requested
   if( verbosity ~= 0 )
   A = [outList{:}]; 
   T = cell2table(A);
   disp('   ');
   disp(T);
   end   
% --------------------------------- final check: Finite classification IDs
   if( cIDnum ~= 0 )            % => classification into F,U,N is required
      for j=1:nLines
      flag = -1;
          if( outList{cIDnum}{j} == 'F' )
          flag = 1;
          elseif( outList{cIDnum}{j} == 'U' )
          flag = 1;
          elseif( outList{cIDnum}{j} == 'N' )
          flag = 1;
          end
         if( flag < 0 )
         error('classification IDs are not all F, U or N');
         end
      end
   end
% ------------------------------------- apply final selection if requested
nFiles = nLines;
   if( sType ~= '-' )
   nFiles = 0;
      for j=1:nLines
         if( outList{cIDnum}{j} == sType )
         nFiles = nFiles + 1;                                       % keep
            for i=1:nCol
            outList{i}{nFiles} = outList{i}{j}; 
            end
         end
      end
   end
   if( nFiles == 0 )
   error(['Empty set found: No Files are detected of type ',sType]);
   end
% ------------------------------------------------------ copy to varargout
varargout = cell(nCol,nFiles);
   for j=1:nCol
       for k=1:nFiles
       varargout{j}{k} = outList{j}{k};
       end
   end
% ----------------------------------------------- record in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['       input file = ',inputFname]);
   if( strcmp(sType,'-') == 0 )   
   fprintf(fid,'%s \n',['        selection = ',sType]);
   else
   fprintf(fid,'%s \n', '        selection = FUN');
   end
fprintf(fid,'%s \n',['number of columns = ',num2str(nCol)]);
fprintf(fid,'%s \n',['  number of files = ',num2str(nFiles)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end

