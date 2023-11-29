function setGraphicFileType(gFileType)
% allows global gvSPLOC.gFileType to be reset
% INPUT
% gFileType
%%                                                         reset gFileType
global gvSPLOC
splocLogFile = gvSPLOC.splocLogFile;
strFileTypes = 'fig_eps_png_jpg_pdf_bmp_tif';
   if( length(gFileType) ~= 3 )
   error('File type not recognized or currently supported');
   elseif( contains(strFileTypes,gFileType) == 0 )
   error('File type not recognized or currently supported');
   end
% --------------------------------------------------- change graphics type
disp(gvSPLOC.gFileType);

   if(  strcmp(gFileType,gvSPLOC.gFileType) == 0 )
% ------------------------------------------ record change in splocLogFile
   gvSPLOC.gFileType = gFileType;
   fid = fopen(splocLogFile,'a');                 % append new information
   fprintf(fid,'%s \n','  ');
% ---------------------------- write info on output file type for graphics
   fprintf(fid,'%s \n',[mfilename,'()']);
   msg = dividerLine('file type for output graphics');
   fprintf(fid,'%s \n',msg);
   msg = ['graphics file type = ',gvSPLOC.gFileType];
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n',dividerLine);
   fclose(fid);  
   % else                                       no change needs to be made
   end
end
