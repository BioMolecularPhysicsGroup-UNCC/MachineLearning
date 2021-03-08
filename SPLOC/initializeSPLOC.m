function initializeSPLOC(varargin) 
% set all global parameters and functions shared across the SPLOC toolset
% 
% INPUT (all optional to avoid default values)
%
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
% commands      purpose/options
% 'fName'       file name of the sploc log file.     
%               'sploc' is default name
% 'gType'       graphics file type:  'read' is the default value
%               read => command line query to user
%               allowed types: {'fig','eps','png','jpg','pdf','bmp','tif'}
% 'wType'       {'fresh','append'} => start new file or append to previous
% 'dType'       density estimation type for ClassifyDataMoments
%               {'kde', 'pdfest'} -> MATLAB kernel density estimator OR 
%                                    Jenny Farmer's PDF Estimator
%                     Note: To use PDF Estimator it must be installed which
%                     can be done easily through MATLABs add-on browser
%                     Details: https://github.com/jennyfarmer/PDFAnalyze
% USAGE: initializeSPLOC                              % all defaults apply
%        initializeSPLOC(3)   % set verbosity = 3 all other defaults apply
%        initializeSPLOC('gType','jpg') 
%        initializeSPLOC('wType','append','gType','png') 
%        initializeSPLOC('fName','test','wType','append','gType','tif')
%        initializeSPLOC(3,'fName','test','wType','append','gType','pdf')
% ------------------------------------------------------------------------
% PROCESS
% clear all figures on screen
% sets all parameters that are dependent on user-defined parameters 
% creates all look-up tables used in other SPLOC functions
% sets voting threshold either by default calculation or user-specified
%
% OUTPUT 
% constant global parameters shared across SPLOC functions
%%                                                                get date
strDate = getDateString();
%%                                                              read input
p = inputParser;
% -------------------------------------------------------------- verbosity
validVerbosity = @(x) isnumeric(x) && isscalar(x) && (x > -1) && (x<4);
addOptional(p,'verbosity',0,validVerbosity)   % default prints out nothing
% -------------------------------------------------------------- file name
check4char = @(x) ischar(x); 
addParameter(p,'fName','sploc',check4char);          % default name: sploc
% ---------------------------------------------------------- graphics type
expectedFtype = {'eps','fig','png','jpg','pdf','bmp','tif'};
addParameter(p,'gType','read', ...
             @(x) any(validatestring(x,expectedFtype)));
% ------------------------------------------------------------- write type
expectedWtype = {'fresh','append'};
addParameter(p,'wType','fresh', ...
             @(x) any(validatestring(x,expectedWtype)));
% ----------------------------------------------------------- density type
expectedDtype = {'kde', 'pdfest'};
addParameter(p,'dType','kde', ...
             @(x) any(validatestring(x,expectedDtype)));
% ------------------------------------------------------------------------
parse(p,varargin{:});
%%                                                  define local variables
verbosity = round(p.Results.verbosity);
gFileType = p.Results.gType;
wType = p.Results.wType;
fName = p.Results.fName;
dType = p.Results.dType;
% -------------------------------------------------- create sploc log file
splocLogFile = [fName,'.log'];
subFolder = 'splocLog';
splocLogFile = getOutputFileName(subFolder,splocLogFile);
% ----------------------------------------------- default verbosity levels
   if( strcmp(wType,'fresh') )    % write basic information on top of file
      if( verbosity == 0 )        % => basic header information is lacking
      verbosity = 1;              % this enforces basic header information 
      end
   end
%%                                                            error checks
   if( strcmp(wType,'append') )   % write basic information on top of file
      if( isfile(splocLogFile) == 0 )
      error(['file: ',splocLogFile,' not found!']);
      end
   end
% ------------------------------------------------------------- user input
%%                                      query for graphic output file type
   if( strcmp(gFileType,'read') )
   disp('  ');
   disp('File types: {fig,eps,png,jpg,pdf,bmp,tif}');
   disp('  <ENTER> => fig');
   gFileType = input('   Enter file type for all graphic files: ','s');
   mmm = length(gFileType);
      if( mmm == 0 )
      gFileType = 'fig';                                   % default value
      else
         if( length(gFileType) ~= 3 )
         error('File type not recognized or currently supported')
         end
      strFileTypes = 'fig_eps_png_jpg_pdf_bmp_tif';
         if( contains(strFileTypes,gFileType) == 0 )
         error('File type not recognized or currently supported')
         end
      end
   end
%%                                         set default graphics parameters
% REMARKS: These set functions are to make the plots ~publication quality.
%          Using the default set functions for axes, font size and figure
%          color are for convenience only.
%fig.InvertHardcopy = 'off';         
set(0,'DefaultFigureColor','white')
set(0,'defaultAxesFontSize',14); 
set(0,'defaultLineLineWidth',1.4);   
set(0,'defaultLineMarkerSize',7); 
set(0,'defaultAxesLineWidth',1.4);
%set(0,'defaultGeographicAxes','on');
%%                                       randomize seed for random numbers
rng('shuffle');
%%                                             clear all figures on screen
close all    
%%                                          create global shared variables
global gvSPLOC
gvSPLOC = struct;
gvSPLOC.strDate = strDate;                         % may be used elsewhere
gvSPLOC.gFileType = gFileType;                    % used in saveas command
gvSPLOC.splocLogFile = splocLogFile;
gvSPLOC.dType = dType;                    % This if for classifyDataStream
%%                                                verbosity control levels
%
% FIX ME  Not implemented yet. This will be an upgrade from 4 to 7 levels.
%               
% Universal definition of verbosity levels will be encoded as:
%                          default=2                     start truth table
%   verbosity = 0       1       2       3       4       5       6     
% splocLogInfo: false   true    true    true    true    true    true  
% detailed_log: false   false   true    true    true    true    true   
% print2Screen: false   false   false   true    true    false   true 
% prntOut2File: false   false   false   false   true    true    false   
% debug2Screen: false   false   false   false   false   false   true
% -------------------------------------------------------- end truth table
% REMARKS:
%               0 => no screen/write records   for debugging <= 6
%     only sploc.log <= 1    quitely write to output <= 5   
%      detailed summary only <= 2       3 => (2) + on screen checking    
%  on screen checking & write to output file <= 4
splocLogInfo = [false, true,  true,  true,  true,  true,  true ];
detailed_log = [false, false, true,  true,  true,  true,  true ];
print2Screen = [false, false, false, true,  true,  false, true ];
prntOut2File = [false, false, false, false, true,  true,  false];
debug2Screen = [false, false, false, false, false, false, true ];
% ------------------------------------------------------------------------
gvSPLOC.splocLogInfo = splocLogInfo;   % =>  log sploc summary information
gvSPLOC.detailed_log = detailed_log;   % => write specific log information
gvSPLOC.print2Screen = print2Screen;   % => print MATLAB figures on screen
gvSPLOC.prntOut2File = prntOut2File;   % =>  print information into a file
gvSPLOC.debug2Screen = debug2Screen;   % => write to screen debugging info
% ---------------------------- set properties for congruency spectrum plot
gvSPLOC.splocSplit = -1;         % default is to use split-spectrum format
gvSPLOC.showGeometricMean = 0;     % default is not to show geometric mean
%%                             set adjustable parameters to default values
% -------------------------------------------- scoring function parameters
gvSPLOC.minScore1 = 2.0;        % score > minScore1 => discriminant effect
gvSPLOC.maxScore0 = 1.3;        % score < maxScore0 => indifference effect
gvSPLOC.add0 = 1.0e-12;      % needed for: variance --> max(variance,add0)
qT = 1.0;     % probably 0.5 is fine, but 1 is better if results are found
gvSPLOC.qualityThreshold = qT;   % quality threshold for an effective mode
% ---------------------------------------------------- efficacy parameters
std_sVS = 0.1; % standard reduction in vecScore for i-modes w.r.t. d-modes
gvSPLOC.std_sVS = std_sVS;                        % standard value for sVS
gvSPLOC.sVS = std_sVS;              % initialize sVS to the standard value
% ----------------------------------------------- random search parameters
nss = 5;                % # of successive successes to inflate scale by 2x
gvSPLOC.nss = nss; 
nsf = nss*5;            % # of successive failures to  shrink scale by x/2
gvSPLOC.nsf = nsf;    
maxFail = 10*nsf;             % shrinks scale by ~ 1/1000 before giving up
gvSPLOC.maxFail = maxFail; 
gvSPLOC.initialScale = 0.02;           % controls step size of random walk
% - - - - - - - - - - - - - - - - - - - - - - - - - - convergence criteria
% REMARK:    slowPC = 5.0 with minPC = 1.0  works well while cuts out fast
% REMARK:    slowPC = 1.0 with minPC = 0.2  works well but a bit too slow?
gvSPLOC.slowPC = 2.0; % percentChange < slowPC => stop reseting nfail to 0
gvSPLOC.minPC = 0.4;                   % finish when percentChange < minPC
%%                                       check score thresholds for sanity
maxScore0 = gvSPLOC.maxScore0;     % need this to calculate overlap pivots
minScore1 = gvSPLOC.minScore1;     % need this to calculate overlap pivots
lnMinScore1 = log(minScore1);
lnMaxScore0 = log(maxScore0);
   if( ( lnMinScore1 - log(1.49999999999) ) < lnMaxScore0 )
   error('MinScore1/MaxScore0 must be at least 1.5 --- preferably 2');
   end
   if( maxScore0 < 1.04999999999 )
   error('naxScore0 must be at least 1.05 --- preferably 1.3');
   end
%%                               check if Gaussian overlap integrals exist
fM = 'overlapMatrix4SPLOC.txt';
fX = 'overlapXvalue4SPLOC.txt';
fY = 'overlapYvalue4SPLOC.txt';
specialFolder = [userpath,'/splocLibrary'];
fileM = getOutputFileName(specialFolder,fM);
fileX = getOutputFileName(specialFolder,fX);
fileY = getOutputFileName(specialFolder,fY);
filesFound = true;
% ---------------------------------------------------
   if( isfile(fileM) == false )
   filesFound = false;
   disp('   ');
   disp('FILE MISSING: overlapMatrix4SPLOC.txt');     
   end
% ---------------------------------------------------
   if( isfile(fileX) == false )
      if( filesFound == true )
      disp('  ');
      disp('FILE MISSING: overlapXvalue4SPLOC.txt');
      filesFound = false;
      else
      disp('              overlapXvalue4SPLOC.txt');
      end
   end
% ---------------------------------------------------
   if( isfile(fileY) == false )
      if( filesFound == true )
      disp('  ');
      disp('FILE MISSING: overlapYvalue4SPLOC.txt');
      filesFound = false;
      else
      disp('              overlapYvalue4SPLOC.txt');
      end
   end
%%                      populate table for interpolating overlap integrals
   if( filesFound )                               % => read the data files
   V0_overlap = dlmread(fileM);
   X0_signal  = dlmread(fileX);
   Y0_ln_r    = dlmread(fileY);
   [nY0,nX0]  = size(V0_overlap);
   else                                    % warn users if data is missing
   disp('    ');
   beep
   disp(['****   WARNING: Place the missing files into ',specialFolder]);
   disp('     OTHERWISE: These data files will be generated/stored next');
   disp('                Estimated time: ~5 CPU minutes to re-calculate');
   disp('  ');
   nextStep = input('enter 0/1 to (stop->place->restart) or calculate: ');
      if( nextStep == 0 )
      disp('   ');
      error(['Next Step: Move missing files to ',specialFolder]);
      end
% -------------------------------- calculate Gaussian overlap lookup table
   t0 = cputime();
   ln_r0 = 0:0.002:0.2;
   ln_r1 = 0.205:0.005:0.5;
   ln_r2 = 0.51:0.01:1;
   ln_r3 = 1.02:0.02:2;
   ln_r4 = 2.04:0.04:4;
   ln_r5 = 4.05:0.05:5;
   Y0_ln_r = -[ln_r0, ln_r1, ln_r2, ln_r3, ln_r4, ln_r5];
   Y0_ln_r = sort(Y0_ln_r);
   nY0 = length(Y0_ln_r);
% --------------------------------------
   signal0 = 0:0.002:0.2;
   signal1 = 0.205:0.005:0.5;
   signal2 = 0.51:0.01:2;
   signal3 = 2.02:0.02:4;
   signal4 = 4.05:0.05:8;
   X0_signal = [signal0, signal1, signal2, signal3, signal4];
   nX0 = length(X0_signal);
   V0_overlap = zeros(nY0,nX0);       % rows => y-axis   columns => x-axis
      if( verbosity == 3 )                % <= this level is for debugging
      test = 1.0e50;
      Xtarget = [0, 0.5, 1, 1.5, 2, 3, 4, 6, 8];
      mX = length(Xtarget);
      iX = zeros(1,mX);
      i = 1;
         for j=1:nX0
         temp = abs( Xtarget(i) - X0_signal(j) ); 
            if( temp < test )
            test = temp;
            iX(i) = j;
            else
            i = i + 1;
               if( i > mX )
               break;
               end
            test = 1.0e50;                                 % reset process
            end
         end
%       disp( num2str(Xtarget) );            % verification all is correct
%       disp( num2str(iX) );
%       disp( num2str(X0_signal(iX)) );
% ------------------------------------------------------------------------
      test = 1.0e50;
      Ytarget = [0, -0.5, -1.0, -2.0, -3.25, -4.5];
      mY = length(Ytarget);
      jY = zeros(1,mY);
      i = 1;
         for j=nY0:-1:1
         temp = abs( Ytarget(i) - Y0_ln_r(j) ); 
            if( temp < test )
            test = temp;
            jY(i) = j;
            else
            i = i + 1;
               if( i > mY )
               break;
               end
            test = 1.0e50;                                 % reset process
            end
         end
%       disp( num2str(Ytarget) );            % verification all is correct
%       disp( num2str(jY) );
%       disp( num2str(Y0_ln_r(jY)) ); 
      end
% ------------------------------------------------------------------------
      for i=1:nX0                                      % loop over columns
      signal = X0_signal(i);
         for j=1:nY0                                      % loop over rows
         r = exp(Y0_ln_r(j));
         dx = r/100;
         xmin = -5;
         xmax = max(5,signal + 5*r); 
         x = xmin:dx:xmax;
         pdf1 = normpdf(x); 
         pdf0 = normpdf(x,signal,r);
% ------------------------------------------ begin check: previously done!
%       norm1 = dx*sum(pdf1);
%       norm0 = dx*sum(pdf0);
%          if( abs(norm1 - 1) > 0.0005 )
%          error('norm1 is not 1');
%          end
%          if( abs(norm0 - 1) > 0.0005 )
%          error('norm0 is not 1');
%          end
% -------------------------------------------------------------- end check
         temp = dx*sum( abs(pdf1 - pdf0) );
         overlap = (2 - temp)/2;
         V0_overlap(j,i) = overlap;
            if( verbosity == 3 )          % <= this level is for debugging
            s1 = sum( (i == iX) );
            s2 = sum( (j == jY) );
               if( s1*s2 ~= 0 )       % plot only a select number of cases
% ---------------------------------------------- plot intermediate results
               figure(1)
               clf
               plot(x,pdf0,'r',x,pdf1,'b');
               xlabel('projected value');
               ylabel('Gaussian probability density');
               title(['sig1 = 1  sig0 = ',num2str(r), ...
                     ' ave1 = 0  ave0 = ',num2str(signal), ...
                     ' overlap = ',num2str(overlap)]); 
               pause
% ----------------------------------------------
               end
            end
         end
      end
   dlmwrite(fileM,V0_overlap);
   dlmwrite(fileX,X0_signal);
   dlmwrite(fileY,Y0_ln_r);
   dtMakeTable = cputime() - t0;  
   end
gvSPLOC.V0_overlap = V0_overlap;
gvSPLOC.X0_signal = X0_signal;
gvSPLOC.Y0_ln_r = Y0_ln_r;
gvSPLOC.ln_rMIN = min(Y0_ln_r) + 0.01;
gvSPLOC.signalMAX = max(X0_signal) - 0.01; 
   if( verbosity == 3 )
   close all
   end
%%                                               generate overlap heat map
   if( verbosity > 1 )
   xmin = 0;
   xmax = 1.5;
   ymin = -log(10^0.601);
   ymax = 0;      
% ------------------------------ must make points on image uniform density
   nLinear = 500;
   dy = (ymax-ymin)/nLinear;  
   dx = (xmax-xmin)/nLinear;
   nLinear1 = nLinear + 1;
   V0uniform = zeros(nLinear1);
   X0uniform = zeros(1,nLinear1);
   Y0uniform = zeros(1,nLinear1);
      for i=1:nLinear1
      ym = ymin + dy*(i - 1);
      xm = xmin + dx*(i - 1);
      X0uniform(i) = xm;
      Y0uniform(i) = ym;
      end
% -------------------------------------
      for i=1:nLinear1
      ym = Y0uniform(i);      
         for j=1:nLinear1
         xm = X0uniform(j);
% ========================================================================
% use explicit linear interpolation method: interp2() is ~75 times slower!
         [~,indx] = sort( abs(X0_signal - xm) );
         iLa = indx(1);
         iLb = indx(2);
         [~,indx] = sort( abs(Y0_ln_r - ym) );
         jLa = indx(1);
         jLb = indx(2);
         %-----------------> y   x
           Vaa = V0_overlap(jLa,iLa);
           Vba = V0_overlap(jLa,iLb);       % Vba =>  b => x    a => y
           Vbb = V0_overlap(jLb,iLb);
           Vab = V0_overlap(jLb,iLa);       % Vab =>  a => x    b => y
         %-----------------> y   x
         xa = X0_signal(iLa);
         xb = X0_signal(iLb);
         ya = Y0_ln_r(jLa);
         yb = Y0_ln_r(jLb);
         yba = yb - ya;
         yma = ym - ya;
         Vam = Vaa + yma*(Vab - Vaa)/yba;
         Vbm = Vba + yma*(Vbb - Vba)/yba;
         V0uniform(i,j) = Vam + (xm - xa)*(Vbm - Vam)/(xb - xa); 
         end
      end
   xmin0 = min(X0_signal);
   xmax0 = max(X0_signal);
   ymin0 = min(Y0_ln_r);
   ymax0 = max(Y0_ln_r);
% ------------------------------------------------------- define file name
   fileOM = getOutputFileName(subFolder,'overlapMap');
      if( verbosity == 3 )
      colorMatrixTool(V0_overlap,1,'fName',fileOM,'fType',gFileType, ...
                     'xLabel','signal','xLimit',[xmin0,xmax0], ...
                     'yLabel','ln(ratio)','yLimit',[ymin0,ymax0], ...
                     'cType','0b');
      else
      colorMatrixTool(V0_overlap,1,'fName',fileOM,'fType',gFileType, ...
                     'xLabel','signal','xLimit',[xmin0,xmax0], ...
                     'yLabel','ln(ratio)','yLimit',[ymin0,ymax0], ...
                     'cType','0b','fShow','no');
      end
   end
%%                                                 generate overlap pivots
% NGOT-check -------------------------------------------------- NGOT-check
Lgotit = (gvSPLOC.minScore1 == 2.0) & ...
         (gvSPLOC.maxScore0 == 1.3);
%
% These thresholds yield the values: aveOverlap0, aveOverlap1
%                                    maxOverlap0, minOverlap0
% NGOT-check -------------------------------------------------- NGOT-check

   if(Lgotit)       % true Lgotit => previously determined characteristics
% ------------------------ previously calculated values with >384k samples
   aveOverlap0 = 0.85305;
   maxOverlap0 = 0.88113;
   aveOverlap1 = 0.56108;
   minOverlap1 = 0.47078;
   count = 37808135;
   dtOverlapPivots = 16.9;
   else       % run simulation to characterize overlap activation function
   t0 = cputime();
      if( verbosity == 3 )
      disp('   ');
      disp(['running simulation to characterize overlap ', ...
            'activation functions']);
      end
   ln_rMIN = gvSPLOC.ln_rMIN;
   signalMAX = gvSPLOC.signalMAX;
% ------------------------------------- set range of reasonable parameters
   max_ave = 10;
   min_ave = -10;
   max_sig = 10.0;
   min_sig = 0.1;
% ---------------------------------
   ln_max_sig = log(max_sig);
   ln_min_sig = log(min_sig);
   d_ln_sig = ln_max_sig - ln_min_sig;
   d_ave = max_ave - min_ave;
% ------------------------------ set band of scores around critical values
   min_lnScore0 = log(maxScore0 - 0.02);
   max_lnScore0 = log(maxScore0 + 0.02);
   min_lnScore1 = log(minScore1 - 0.02);
   max_lnScore1 = log(minScore1 + 0.02);
% ------------------------------------------------------ randomize samples
   Nsamples = 200000;
   overlap = zeros(1,Nsamples);
   lnScore = zeros(1,Nsamples);
   add0 = 1.0e-12;
   count = 0;
   s = 0;
      while( s < Nsamples )
      count = count + 1;
      ln_sig0 = ln_min_sig + rand()*d_ln_sig;
      ln_sig1 = ln_min_sig + rand()*d_ln_sig;
      sig0 = exp(ln_sig0);
      sig1 = exp(ln_sig1);
      ave0 = min_ave + rand()*d_ave;
      ave1 = min_ave + rand()*d_ave; 
      var0 = max(sig0*sig0,add0); 
      var1 = max(sig1*sig1,add0);
      r1 = sqrt(var0/var1);
      r2 = 1/r1;
      r = max(r1,r2);
      noise = sqrt(var1 + var0); 
      signal = abs(ave1 - ave0);
      s2nr = signal/noise;                         % signal to noise ratio
      r = r - 1;
      test_lnScore = log(sqrt(r*r + s2nr*s2nr) + 1);   % add the 1 back in
      L0 = and( (test_lnScore > min_lnScore0) , ... 
                (test_lnScore < max_lnScore0) );
      L1 = and( (test_lnScore > min_lnScore1) , ...
                (test_lnScore < max_lnScore1) );
         if( L0 || L1 )
         s = s + 1;
         else
         continue
         end
      lnScore(s) = test_lnScore;
% ------------------------------------------------- now let r = min(r1,r2)
      r = min(r1,r2);
      max_sig = max(sig0,sig1);
      xm = abs(ave1 - ave0)/max_sig;                         % xm = signal
      ym = log(r);                                             % ym = ln_r
         if( (ym < ln_rMIN) || (xm > signalMAX) )
         overlap(s) = 0;                            % drop overlap to zero
         else
         %overlap(s) = interp2(X0,Y0,V0_overlap,signal,ln_r);
% ========================================================================
% Remark:  The linear interpolation method: interp2() is ~75 times slower!
         [~,indx] = sort( abs(X0_signal - xm) );
         iLa = indx(1);
         iLb = indx(2);
         [~,indx] = sort( abs(Y0_ln_r - ym) );
         jLa = indx(1);
         jLb = indx(2);
         %-----------------> y   x
           Vaa = V0_overlap(jLa,iLa);
           Vba = V0_overlap(jLa,iLb);       % Vba =>  b => x    a => y
           Vbb = V0_overlap(jLb,iLb);
           Vab = V0_overlap(jLb,iLa);       % Vab =>  a => x    b => y
         %-----------------> y   x
         xa = X0_signal(iLa);
         xb = X0_signal(iLb);
         ya = Y0_ln_r(jLa);
         yb = Y0_ln_r(jLb);
         yba = yb - ya;
         yma = ym - ya;
         Vam = Vaa + yma*(Vab - Vaa)/yba;
         Vbm = Vba + yma*(Vbb - Vba)/yba;
         overlap(s) = Vam + (xm - xa)*(Vbm - Vam)/(xb - xa); 
% ========================================================================
         end
      end
% ---------------------------------------- find magic 50/50 overlap values
   L0 = and( (lnScore > min_lnScore0) , (lnScore < max_lnScore0) );
      if( sum(L0) == 0 )
      error('estimate for overlap pivot for maxScore0 failed');
      end
   aveOverlap0 = mean( overlap(L0) );
   maxOverlap0 = max( overlap(L0) );
   L1 = and( (lnScore > min_lnScore1) , (lnScore < max_lnScore1) );
      if( sum(L1) == 0 )
      error('estimate for overlap pivot for minScore1 failed');
      end
   aveOverlap1 = mean( overlap(L1) );
   minOverlap1 = min( overlap(L1) );
   dtOverlapPivots = cputime() - t0;
      if( verbosity > 1 )
         if( verbosity == 3 )
         figure;
         figure_number = get(gcf,'Number');
         figMagicOverlapVals = figure(figure_number);
         else
         figMagicOverlapVals = figure('visible','off');
         end
      clf
% ---------------------------------------------------------- graph results
      max_lnScore = max(lnScore);
      [xx,indx] = sort(overlap);
      yy = lnScore(indx);
      plot(xx,yy,'bo');
% ------------------------------------ horizontal line
      x1 = 0:0.01:1;
      y1 = log(1.3)*ones( size(x1) );
      hold on
      plot(x1,y1,'r');
% ------------------------------------ horizontal line
      x2 = x1;
      y2 = log(2)*ones( size(x1) );
      plot(x2,y2,'r');
% ------------------------------------ vertical line
      x3 = [aveOverlap0,aveOverlap0];
      y3 = [0,max_lnScore];
      plot(x3,y3,'k');
% ------------------------------------ vertical line
      x4 = [aveOverlap1,aveOverlap1];
      y4 = [0,max_lnScore];
      plot(x4,y4,'k');
      xlabel('overlap');
      ylabel('ln[score]');
      figNameMagic = getOutputFileName(subFolder, ...
                                  'overlapThresholds');
      saveas(figMagicOverlapVals,figNameMagic,gvSPLOC.gFileType);
      end
   end
%%                               initialize ln(score) activation functions
% --------------------------------------- with respect to scoring function
x = 0:0.001:1.502;        % dx = 0.001 must stay a constant DO NOT CHANGE!
xMax = 1.5;                                   % defined here as a constant
% ---------------------------------------------------- discriminant weight
c1 = log(9)/16;              % characterizes spread in activation function
dWt = 1./(1 + exp( -16*(x - lnMinScore1) ) );        % 16 is a good number
pw = max(1, (2 - x/lnMinScore1) );
dWt = dWt.^pw;
% ---------------------------------------------------- indifference weight
c0 = log(9)/16;              % characterizes spread in activation function
iWt = 1 - 1./(1 + exp( -16*(x - lnMaxScore0) ) );    % 16 is a good number
pw = min(1, (x/lnMaxScore0).^2 );
iWt = iWt.^pw;
% ---------------------------------------------------- undetermined weight
uWt = max(dWt,iWt);
% -------------------------------- record results in global SPLOC varibles
gvSPLOC.dWtS = dWt;               % discriminant probability weight factor
gvSPLOC.iWtS = iWt;               % indifference probability weight factor
gvSPLOC.uWtS = uWt;               % undetermined probability weight factor
gvSPLOC.xS = x;              % defines range of scores that are considered
gvSPLOC.xMax = xMax;                     % maximum value of score required
% ---------------------------------------------- record based on verbosity
   if( verbosity > 1 )
      if( verbosity == 3 )
      figure;
      figure_number = get(gcf,'Number');
      figVoteWeight = figure(figure_number);
      else
      figVoteWeight = figure('visible','off');
      end
   clf
   Ld = ( x > lnMinScore1 );
   Lu = and( (x > lnMaxScore0) , (x < lnMinScore1) );
   Li = ( x < lnMaxScore0 );
   plot(x(Ld),dWt(Ld),'b');
   hold on;
   plot(x(Li),iWt(Li),'k');
   plot(x(Lu),uWt(Lu),'r');
   x0 = zeros(1,2);
   x1 = zeros(1,2);
   y2 = zeros(1,2);
   x0(1) = lnMaxScore0;
   x0(2) = lnMaxScore0;
   x1(1) = lnMinScore1;
   x1(2) = lnMinScore1;
   y2(1) = 0;
   y2(2) = 1;
   plot(x0,y2,'m--');
   plot(x1,y2,'g--');
   xlabel('ln(score)');
   ylabel('voter fraction');
   legend('discriminant','indifference','undetermined', ...
          'lower-band','upper-band','location','best');
   figNameVote = getOutputFileName(subFolder, ...
                                  'lnScoreActivationFunctions');
   saveas(figVoteWeight,figNameVote,gvSPLOC.gFileType);
   end
%%                                 initialize overlap activation functions
% ------------------------------------------------ with respect to overlap
dx = 0.001;
x = 0:dx:1.0;                                           % => overlap value
% ---------------------------------------------------- discriminant weight
b1 = aveOverlap1 - minOverlap1;
a1 = log(9)/b1;                       % => at x=minOverlap1  set dWt = 0.9
dWt = 1./(1 + exp( a1*(x - aveOverlap1) ) );
pw = min(sqrt(x/aveOverlap1),1);
dWt = dWt.^pw;
% ---------------------------------------------------- indifference weight
b0 = maxOverlap0 - aveOverlap0;
a0 = log(9)/b0;                       % => at x=maxOverlap0  set iWt = 0.9
iWt = 1 - 1./(1 + exp( a0*(x - aveOverlap0) ) );
pw = min(sqrt(1 - x)/sqrt(1 - aveOverlap0) , 1); 
iWt = iWt.^pw;
% -------------------------------- record results in global SPLOC varibles
gvSPLOC.dWtOVL = dWt;             % discriminant probability weight factor
gvSPLOC.iWtOVL = iWt;             % indifference probability weight factor
gvSPLOC.dxOVL = dx;
% ----------------------------------------------------- check and document
   if( verbosity > 1 ) 
   uWt = max(dWt,iWt);                               % undetermined weight
   Ld = ( x < aveOverlap1 ); 
   Lu = and( (x > aveOverlap1) , (x < aveOverlap0) );
   Li = ( x > aveOverlap0 );
      if( verbosity == 3 )
      figure;
      figure_number = get(gcf,'Number');
      figVoteWeight = figure(figure_number);
      else
      figVoteWeight = figure('visible','off');
      end
   clf
   plot(x(Ld),dWt(Ld),'b');
   hold on;
   plot(x(Li),iWt(Li),'k');
   plot(x(Lu),uWt(Lu),'r');
   x0 = zeros(1,2);
   x1 = zeros(1,2);
   y2 = zeros(1,2);
   x0(1) = aveOverlap0;
   x0(2) = aveOverlap0;
   x1(1) = aveOverlap1;
   x1(2) = aveOverlap1;
   y2(1) = 0;
   y2(2) = 1;
   plot(x0,y2,'m--');
   plot(x1,y2,'g--');
   xlabel('overlap');
   ylabel('voter fraction');
   legend('discriminant','indifference','undetermined', ...
          'lower-band','upper-band','location','best');
   figNameVote = getOutputFileName(subFolder, ...
                                  'overlapActivationFunctions');
   saveas(figVoteWeight,figNameVote,gvSPLOC.gFileType);
% ----------------------------------------- plot voteweight on 2D heat map
      if( verbosity == 3 )
      gridData = zeros(size(V0uniform));             % temporary 2D matrix
      for j=1:nLinear1           % loop over columns in V0uniform (x-axis)
         for i=1:nLinear1           % loop over rows in V0uniform (y-axis)
         x = V0uniform(i,j);
         k = ceil(x/dx + 0.1);
         yd = dWt(k);
         yi = iWt(k);
            if( yd > yi )
               if( yd > 0.01 )
               gridData(i,j) = (1 - yd)/2;
               else
               gridData(i,j) = 0.5;
               end
            else
               if( yi > 0.01 )
               gridData(i,j) = (1 + yi)/2;
               else
               gridData(i,j) = 0.5;
               end
            end
         end
      end  
      bLstr = ['discriminant (~0) to undetermined (~0.5) to ', ...
              'indifference (~1)'];
% ------------------------------------------------------- define file name
      fileO = getOutputFileName(subFolder,'opinionMap');
         if( verbosity == 3 )
         colorMatrixTool(gridData,1,'fName',fileO,'fType',gFileType, ...
                      'xLabel','scaled signal','xLimit',[xmin0,xmax0], ...
                      'yLabel','log10(ratio)', 'yLimit',[ymin0,ymax0], ...
                      'bLabel',bLstr,'tLabel','opinion map', ...
                      'cType','bwr');
         else
         colorMatrixTool(gridData,1,'fName',fileO,'fType',gFileType, ...
                        'xLabel','signal','xLimit',[xmin0,xmax0], ...
                        'yLabel','ln(ratio)','yLimit',[ymin0,ymax0], ...
                        'bLabel',bLstr,'tLabel','opinion map', ...
                        'cType','bwr','fShow','no');
         end
      end
% -------------------------------------------- plot voter threshold curves
      if( verbosity == 3 )
      figure;
      figure_number = get(gcf,'Number');
      figVT = figure(figure_number);
      else
      figVT = figure('visible','off');
      end
   hold on;
   nVref = 100;                                       % set as a reference
   nSmin = logspace(2,6,100); 
   x = nSmin/nVref;
   ratio = 1;
      for i=1:12
      % CAUTION: when function is changed here, make corresponding change
      % to the function called:  getDefaultVoteThreshold(nS0,nS1,nDOF)
      % REMARKS:
      % nD1 = n1*samplesize       and   nD0 = n0*sampleSize
      % nS0 = nD0/sqrt(n0);     % nS0 = sqrt(n0)*sampleSize = nD0/sqrt(n0)
      % nS1 = nD1/sqrt(n1);     % nS1 = sqrt(n1)*sampleSize = nD1/sqrt(n1)
      % nSmin = min(nS1,nS0);
      % nSmax = max(nS1,nS0);
      % Neff = 0.5*( nSmin + nSmin* ( 2*nSmax/(nSmin + nSmax) )^2 );      
      Neff = nSmin*0.5*( 1 + ( 2/(1 + 1/ratio) )^2 );
      pw = max(0,0.5 - 0.49999*sqrt(nVref./Neff) );
      voteThreshold = min(0.95,0.5 + 0.7*(nVref./Neff).^pw);
      plot(log10(x),voteThreshold)
      ratio = 3*ratio;
      end
   xlabel('log10(nS/nV)');
   ylabel('vote threshold');
   title('1 < nSmax/nSmin < 177200'); 
   figNameVote = getOutputFileName(subFolder, ...
                                  'voteThreshold');
   saveas(figVT,figNameVote,gvSPLOC.gFileType);
   end
%error('stop here for now');
%%                                           write information to log file
 if( verbosity == 0 )                     % => append to existing log file
fid = fopen(splocLogFile,'a');          % continues in previous named file
fprintf(fid,'%s \n',' ');
msg='Supervised Projective Learning with Orthogonal Completeness (SPLOC)';
fprintf(fid,'%s \n',msg);
version = splocToolsetVersion();
msg = ['MATLAB SPLOC toolset ',version];
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',[mfilename,'()']);
msg = ['date = ',strDate];
fprintf(fid,'%s \n',msg);
% ---------------------------- write info on output file type for graphics
msg = dividerLine('file type for output graphics');
fprintf(fid,'%s \n',msg);
msg = ['graphics file type = ',gFileType];
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',dividerLine);
 else                                         % => starting a new log file
    if( strcmp(wType,'fresh') )
    fid = fopen(splocLogFile,'w');                 % spawns a new log file
    else
    fid = fopen(splocLogFile,'a');      % continues in previous named file
    fprintf(fid,'%s \n',' ');
    end
msg ='Supervised Projective Learning for Orthogonal Congruences (SPLOC)';
fprintf(fid,'%s \n',msg);
version = splocToolsetVersion();
msg = ['MATLAB SPLOC toolset ',version];
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',[mfilename,'()']);
msg = ['date = ',strDate];
fprintf(fid,'%s \n',msg);
% ---------------------------- write info on output file type for graphics
msg = dividerLine('file type for output graphics');
fprintf(fid,'%s \n',msg);
msg = ['graphics file type = ',gFileType];
fprintf(fid,'%s \n',msg);
% -------------------------- write info on congruency spectrum preferences
msg = dividerLine('congruency spectrum toggles');
fprintf(fid,'%s \n',msg);
msg = ['split SPLOC spectrum = ',num2str(gvSPLOC.splocSplit)];
fprintf(fid,'%s \n',msg);
msg = [' show geometric mean = ',num2str(gvSPLOC.showGeometricMean)];
fprintf(fid,'%s \n',msg);
msg = [' i-mode/d-mode scale = ',num2str(std_sVS)];
fprintf(fid,'%s \n',msg);
% ---------------------------------------------------- employed thresholds
msg = dividerLine('Employed thresholds');
fprintf(fid,'%s \n',msg);
msg = ['      maximum allowed indifference score = ',num2str(maxScore0)];
fprintf(fid,'%s \n',msg);
msg = ['      minimum allowed discriminant score = ',num2str(minScore1)];
fprintf(fid,'%s \n',msg);
msg = ['minimum allowed quality for d- or i-mode = ',num2str(qT)];
fprintf(fid,'%s \n',msg);
% --------------------------- write info on ln(score) activation functions
msg = dividerLine('ln(score) activation functions');
fprintf(fid,'%s \n',msg);
msg = ['indifference ln(score) pivot = ',num2str(lnMaxScore0)];
fprintf(fid,'%s \n',msg);
msg = ['indifference ln(score) range = ',num2str(c0)];
fprintf(fid,'%s \n',msg);
msg = ['discriminant ln(score) pivot = ',num2str(lnMinScore1)];
fprintf(fid,'%s \n',msg);
msg = ['discriminant ln(score) range = ',num2str(c1)];
fprintf(fid,'%s \n',msg);
% ------------------------------- write info on Gaussian overlap integrals
msg = dividerLine('Gaussian overlap integral tables');
fprintf(fid,'%s \n',msg);
   if( filesFound )
   msg = 'input files: overlapMatrix4SPLOC.txt';
   fprintf(fid,'%s \n',msg);
   msg = '             overlapXvalue4SPLOC.txt';
   fprintf(fid,'%s \n',msg);
   msg = '             overlapYvalue4SPLOC.txt';
   fprintf(fid,'%s \n',msg);
   else
   msg = 'output files: overlapMatrix4SPLOC.txt';
   fprintf(fid,'%s \n',msg);
   msg = '              overlapXvalue4SPLOC.txt';
   fprintf(fid,'%s \n',msg);
   msg = '              overlapYvalue4SPLOC.txt';
   fprintf(fid,'%s \n',msg);
   msg = ['file generation: total CPU seconds = ',num2str(dtMakeTable)];
   fprintf(fid,'%s \n',msg); 
   end
msg = ['   length of X0_signal = ',num2str(nX0)];
fprintf(fid,'%s \n',msg);
msg = ['   length of Y0_ln_r   = ',num2str(nY0)];
fprintf(fid,'%s \n',msg);
msg = ['total # of grid points = ',num2str(nX0*nY0)];
fprintf(fid,'%s \n',msg);
   if( Lgotit == false )
   msg = 'Caution: New Gaussian overlap table!';
   fprintf(fid,'%s \n',msg);
   msg = 'Accuracy/speed character will change';
   fprintf(fid,'%s \n',msg);
   msg = 'Modify the program at NGOT-check!';
   fprintf(fid,'%s \n',msg);
   msg = 'This warning implies potential error!';
   fprintf(fid,'%s \n',msg);
   msg = 'Once benchmarked, warning goes away';
   fprintf(fid,'%s \n',msg);
   end
% ----------------------------- write info on overlap activation functions
msg = dividerLine('overlap activation functions');
fprintf(fid,'%s \n',msg);
msg = ['indifference overlap pivot = ',num2str(aveOverlap0)];
fprintf(fid,'%s \n',msg);
msg = ['indifference overlap range = ',num2str(b0)];
fprintf(fid,'%s \n',msg);
msg = ['discriminant overlap pivot = ',num2str(aveOverlap1)];
fprintf(fid,'%s \n',msg);
msg = ['discriminant overlap range = ',num2str(b1)];
fprintf(fid,'%s \n',msg);
msg = ['          number of trials = ',num2str(count)];
fprintf(fid,'%s \n',msg);
msg = [' net sampling: cpu seconds = ',num2str(dtOverlapPivots)];
fprintf(fid,'%s \n',msg);
% ---------------------- write info on heristic formula for vote threshold
msg = dividerLine('vote threshold heristic formula');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n','   nV = # of variables for all systems');
fprintf(fid,'%s \n','  nD0 = total # of data observations for systems-0');
fprintf(fid,'%s \n','  nD1 = total # of data observations for systems-1');
fprintf(fid,'%s \n','nDmin = MIN(nD0,nD1)');
fprintf(fid,'%s \n','nDmax = MAX(nD0,nD1)');
fprintf(fid,'%s \n',[' Neff = 0.5*( nDmin + nDmin*', ...
                     '( 2*nDmax/(nDmin + nDmax) )^2 )']);
fprintf(fid,'%s \n','   pw = max(0,0.5 - 0.49999*sqrt(nV/Neff) )');
fprintf(fid,'%s \n','   vT = min(0.95, 0.5 + 0.9*(nV/Neff).^pw)');
% -------------------------------- write info on consistency score formula
msg = dividerLine('consistency rank score formula');
fprintf(fid,'%s \n',msg);
msg = '      N = # of sploc solutions under identical conditions';
fprintf(fid,'%s \n',msg);
msg = '   MSIP = mean squared inner product';
fprintf(fid,'%s \n',msg);
msg = 'msip_dd = NxN matrix for discriminant-discriminant pairwise MSIPs';
fprintf(fid,'%s \n',msg);
msg = 'msip_ii = NxN matrix for indifference-indifference pairwise MSIPs';
fprintf(fid,'%s \n',msg);
msg = 'msip_di = NxN matrix for discriminant-indifference pairwise MSIPs';
fprintf(fid,'%s \n',msg);
msg = '  score = 50*[sum(msip_dd) + sum(msip_ii) - 2*sum(msip_di)]/N';
fprintf(fid,'%s \n',msg);
msg = '          maximum possible score = 100, used to rank consistency';
fprintf(fid,'%s \n',msg);
msg = '          a greater score indicates a greater consistency';
fprintf(fid,'%s \n',msg);
msg = '          a smaller score indicates high likelihood of outliers';
fprintf(fid,'%s \n',msg);
% new section goes here
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
 end
end
