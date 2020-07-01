function verbosity = setVerbosity(verbosity)
% enforces verbosity input by user to have the admissible values {0,1,2,3} 
% used in many functions that are part of the SPLOC toolbox
% verbosity: Specifies amount of intermediate processing steps to report.
%            Functions exist to output certain data, while verbosity has
%            no connection to that type of output and how it is reported.
%       ---> verbosity extracts hidden information within SPLOC functions.
%       ---> verbosity controls optional intermediate results ONLY:
% default: 0 => no output except what might get dumped into main log file.
% process: 1 => same as 0, writing key processing data in other log files.
% summary: 2 => same as 1, writing figures to files showing key relations.
% display: 3 => same as 2, with figures displayed on the screen, sometimes
%               with a pause, and/or with additional output printed to the 
%               command window. Useful to understand how the code works!
%       ---> All printed figures have the same file type (fig, png, etc).
%       ---> Other log files are written in separate folders/directories.
% ------------------------------------------------------------------------
   if( isnumeric(verbosity) )
   verbosity = round(verbosity);
      if( verbosity < 0 )
      verbosity = 0;                              % => no output to screen
      elseif( verbosity > 3 )
      verbosity = 3;
      end
   else
   error('verbosity must be {0, 1, 2, 3}');
   end
end
