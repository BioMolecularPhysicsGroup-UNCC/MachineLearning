function lineName = dividerLine(str)
% creates a string variable used to print divider lines in various files
% INPUT
% str = string to label the line
%
% OUTPUT
% lineName = string variable that is a line with name appended at the end
%%                                                     create divider line
line80 = ['----------------------------------------', ...
          '----------------------------------------'];     % 80 characters
   if( nargin < 1 )
   lineName = line80;
   else
   m = length(str);
      if( m > 74 )
      lineName = ['----- ',str];
      else
      n = 79 - m;
      lineName = [line80(1:n),' ',str];
      end
   end
end

