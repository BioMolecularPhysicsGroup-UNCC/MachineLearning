function fDate = getDateString()
% makes a string out of the current date
% use as a simple tool for SPLOC toolset
str = datestr(datetime);
L = isspace(str);
L = or(L, (str == ':') );
L = or(L, (str == '-') );
str(L) = '_';
fDate = str;
end
