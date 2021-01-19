function [nFiles,pID,fileList] = readPropertySpecs(inputFname,varargin)
% process input file specifying list of file names with property IDs
%
% INPUT
% inputFname = name of input file
%
% generic example input file format having three distinct properties
% Each property is assigned F, U or N. 
% F => property of interest is observed (functioning).
% N => property of interest is not observed (nonfunctioning).
% U => property of interest is not measured or undetermined (unknown).
% P1 => property 1.    P2 => property 2.    P3 => property 3.
%
% expected format: 
% S    filename1 filename2 ...   => skip system altogether 
% UUU  filename1 filename2 ...   => unknown properties across the board
% FUU  filename1 filename2 ...   => P1 is functional, P2,P3 are unknown
% FNU  filename1 filename2 ...   => P1 is functional, P2 is nonfunctional
%                                    and P3 is unknown
% FFN  filename1 filename2 ...   => P1,P2 are functional, P3 nonfunctional
% NNF  filename1 filename2 ...   => P1,P2 are nonfunctional, P3 functional
% FUN  filename1 filename2 ...   => P1 is functional, P2 is unknown and P3
%                                   is nonfunctional
%%
% USAGE:
%   readPropertySpecs(inputFname)
%   readPropertySpecs(inputFname,verbosity)
%   readPropertySpecs(inputFname,verbosity,fullPath)
%               
% verbosity: 0 => no write to screen
%            1,2 => write output to screnn
%            3 => same as 1,2 plus some additional 
%  fullPath: 0 => look for inputFname file within the input directory
%            1 => assume inputFname contains the full path of interest.
%% 
% PROCESS
% A list of filenames are on the end of each row.
% The 1st column defines a string with length equal to # of properties.
% Each character position is given a character F, U or N.
% F => the property is functional. 
% N => the property is nonfunctional. 
% U => whether the property is functional or nonfunctional is unknown.
% all strings have the same length except if it equals S. 
% S is a special case that allows a user to skip data.
% This function reads this information and skips what needs to be skipped.
%%
% OUTPUT (different number of cell arrays are possible)
% [nFiles, pID, {verbatim filename columns}]
% nFiles = number of files that are to be processed
% pID = property identification string that specifies known functional
%       state of each property (F, N or U for unknown).
% fileList = cell array for all columns read after the pID for filenames
%
% REMARK: The order of information in the input file is preserved.
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
%%                                                  parse input parameters
parse(p,varargin{:});
% ------------------------------------------------- create local variables
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
   if( nCol == 1 )                         % => nCol = 1 is a special case
   error('only 1 column found: must classify something'); 
   end
%%                                         remove all systems labeled as S
nFiles = 0;
strLength = -1;
   for k=1:nLines
   str = upper(tempStr{k,1});               % force all capital letter IDs
      if( strcmp(str,'S') == 0 )                          % => do not skip
      nFiles = nFiles + 1;
% ---------------------------------------------------------- error check 1
        if( strLength < 0 )
        strLength = length(str);
        else
        itest = length(str);
           if( itest ~= strLength )
           error('property ID string length is not same for all systems'); 
           end
        end
% ---------------------------------------------------------- error check 2
        for j=1:strLength 
        cstr = str(j);
           if( cstr ~= 'F' && cstr ~= 'U' && cstr ~= 'N' )
           error('Property ID characters are not restricted to FUN');
           end
        end
      end
   tempStr{k,1} = str;
   end
% ------------------------------------------------------------ error check
   if( nFiles < 2 )
   error('need at least two files for splocing');
   end
%%                                          create corresponding file list
jj = nCol - 1;                % remove first column for string property ID
fileList = cell(nFiles,jj);
pID = cell(nFiles,1);
nFiles = 0;
   for k=1:nLines
      if( strcmp(tempStr{k,1},'S') == 0 )                 % => do not skip
      nFiles = nFiles + 1;
         for j=2:nCol
         jj = j - 1;
         fileList{nFiles,jj} = tempStr{k,j};
         end
      pID{nFiles} = tempStr{k,1}; 
      end
   end
%%                                  print results to screen when requested
   if( verbosity ~= 0 )
   T = table(pID,fileList);
   disp('   ');
   disp(T);
   end    
nFiles = nLines;
fileList = fileList';
pID = pID';
% ----------------------------------------------- record in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['       input file = ',inputFname]);
fprintf(fid,'%s \n',['number of columns = ',num2str(nCol)]);
fprintf(fid,'%s \n',['  number of files = ',num2str(nFiles)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end

