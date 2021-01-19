function longFileName = getOutputFileName(subFolder,fName)
% creates subfolder and/or augments folder as a prefix to output file name
% use as a simple tool for SPLOC toolset
   if( isfolder(subFolder) == 0 )                       % create subfolder
   mkdir(subFolder);
   end
% ----------------------------------- path format dependent on system type 
   if( ispc )              
   longFileName = [subFolder,'\',fName];             % Windows/DOS systems
   else
   longFileName = [subFolder,'/',fName];              % Linux/UNIX systems 
   end
end
