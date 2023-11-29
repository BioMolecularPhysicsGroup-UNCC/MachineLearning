function inFileName = getInputFileName(targetFolder,fName)
% looks in targetfolder by making targetfolder a prefix to input filename
% use as a simple tool for SPLOC toolset
   if( isfolder(targetFolder) == 0 )         % targetfolder does not exist 
   error(['target folder not found. Target: ',targetFolder]);
   end
% ----------------------------------- path format dependent on system type 
   if( ispc )              
   inFileName = [targetFolder,'\',fName];            % Windows/DOS systems
   else
   inFileName = [targetFolder,'/',fName];             % Linux/UNIX systems 
   end
   if( isfile(inFileName) == 0 )                     % file does not exist
   disp('   ');
   disp('target folder exist, but requested file is not present');
   error(['target file: ',inFileName]);
   end
end
