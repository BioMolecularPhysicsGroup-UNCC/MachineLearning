% automate experiment B many times for the PNAS paper
% COMMENTS:
% This is a script with a function in it.
% The design is made to launch this function as a job array.
% Recommendation for HPC:
% Set nExp = 15.
% run 1 bucket and 1 particular pair at a time. 
% request 1 core per job.
% This strategy will require 432 cores, and you would have the request
% inside a job array. Estimated time to complete 15 jobs is 24 hours.
% Hence, one could set wall time for 50 hours for each request.
%%                                                  output data file names
%   xyz_abc_jkSout.dlm
%             ^------------> S or B for small or big system size
%            ^-------------> minGap10 = 0,1,2, ... 99  (two digits)
%           ^--------------> ic       = 0,1    
%         ^----------------> c        = F, L, T
%        ^-----------------> b        = F, L, T, S
%       ^------------------> a        = E, F
%     ^--------------------> z        = F, L, T
%    ^---------------------> y        = F, L, T, S
%   ^----------------------> x        = E, F
%
% examples: EST_EFF_120Sout.dlm
%           FLF_EFT_025Bout.dlm
%           ETT_FFF_130Sout.dlm
%           FTL_FSL_015Bout.dlm
% ------------------------------------------------------- clear out memory
close all
clearvars
%%                                        example ways to run the function
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                             EXAMPLE 1
% Play with the parameters to understand the capabilities of the program
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%%                                      case 1: Do what the user specified
%{
%%                                                 user defined parameters
jL_bucket = 4;      % options are 1,2,3,4 for xLz,xSz,xTz,xFz respectively
startingPair = 1;    % allows one to fill in last set of unfinished runs
%                      % typical values will be 1 or (1 + file count)
endingPair = 15;            % allows one to stop at any point of interest
%                               % typical values will be file count or 108
nExp = 2;         % = # of experiments per  xyzF & abcN  initial condition
%                 % typical range will be from 1 to 10      do not use >25
% -------------------------------------------------------- select data set
nSamples = 20000;    % (500,20000)         using pregenerated synthetic data
%                  % 500   => trainingSet1
%                  % 20000 => trainingSet1B
% ------------------------------------------ less frequent flags to adjust
% DEFAULT SETTINGS:  deletePrevious = outputScopeOnly = showT = false   
deletePrevious = false;  % (true,false) => delete previous files or append
outputScopeOnly = false;                    % false => run the experiments
%               % true => prints to the screen the pairs that will be done
%                         No calculations will be performed. 
%                         The true flag is to verify you get what you want
showT = false;  % (true,false) => (show, no show) classification on screen
show_vb = 1;   % verbosity level => vb = 1,2,3
               % 1 => minimal output from sploc
               % 2 => iterative output from sploc
               % 3 => full fledge navigation windows from sploc
% ----------------------------------------- least frequent flags to adjust
% DEFAULT SETTINGS:  minGap = 2   ic = 0
minGap = 2;      % minimum gap in log probability between F/N molecules
%                % 0 => no gap criteria is used. 
%          REMARK: small minGap => will likely miss functional molecules
%                  large minGap => cost of failed experiments will be high
%     
ic = 0;      % ic = (0,1) => (best PCA, identity)-basis initial conditions
%              ic = 0 should be the option to use.
M = true;    % (true,false) to use memory during iterative learning
BS = true;   % bootstrap = (true,false)
% ------------------------------------------------------------ error check
   if( (nSamples ~= 500) && (nSamples ~= 20000) )
   error('incompatible # of samples: nSamples can be 500 or 20000');
   end
   if( (show_vb < 1) || (show_vb > 3) )
   error('show_vb is out of range: allowed values are: 1,2,3');
   end
runAutomateExpB(nSamples,jL_bucket,startingPair,endingPair,nExp,M,BS, ...
                deletePrevious,outputScopeOnly,showT,show_vb,minGap,ic);
%}

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                             EXAMPLE 2
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%%                                Translate this to a job array on the HPC
% {
nExp = 3;         % = # of experiments per  xyzF & abcN  initial condition
%                 % typical range will be from 1 to 10      do not use >25
% ------------------------------------------ less frequent flags to adjust
% DEFAULT SETTINGS:  deletePrevious = outputScopeOnly = showT = false   
deletePrevious = false;  % (true,false) => delete previous files or append
outputScopeOnly = false;                    % false => run the experiments
%               % true => prints to the screen the pairs that will be done
%                         No calculations will be performed. 
%                         The true flag is to verify you get what you want
showT = false;  % (true,false) => (show, no show) classification on screen
show_vb = 1;   % verbosity level => vb = 1,2,3
               % 1 => minimal output from sploc
               % 2 => iterative output from sploc
               % 3 => full fledge navigation windows from sploc
% -------------------------------------------------------- select data set
nSamples = 20000;  % (500,20000)
%                  % 500   => trainingSet1
%                  % 20000 => trainingSet1B
% ----------------------------------------- least frequent flags to adjust
% DEFAULT SETTINGS:  minGap = 2   ic = 0
minGap = 2;      % minimum gap in log probability between F/N molecules
%          REMARK: small minGap => will likely miss functional molecules
%                  large minGap => cost of failed experiments will be high
ic = 0;      % ic = (0,1) => (best PCA, identity)-basis initial conditions
M = true;    % (true,false) to use memory during iterative learning
BS = true;   % bootstrap = (true,false)
% ------------------------------------------------------------ error check
   if( (nSamples ~= 500) && (nSamples ~= 20000) )
   error('incompatible # of samples: nSamples can be 500 or 20000');
   end
% ----------------------------------------------------------- do job array
   for jL_bucket=1:4
      for i=1:108
      startingPair = i;   % DO THIS ASSIGNMENT:  Never cross memory!!!!!!!
      endingPair = i;     % DO THIS ASSIGNMENT:  Never cross memory!!!!!!!
      runAutomateExpB(nSamples,jL_bucket,startingPair,endingPair,nExp, ...
             M,BS,deletePrevious,outputScopeOnly,showT,show_vb,minGap,ic);
      end
   end 
%}
%%                                                            the function
function runAutomateExpB(nS,jL_bucket,startingPair,endingPair,nExp,M, ...
                BS,deletePrevious,outputScopeOnly,showT,show_vb,minGap,ic)
% INPUT 
% nS               number of samples in simulation data
% jL_bucket        options are 1,2,3,4 for xLz,xSz,xTz,xFz respectively
% startingPair     specifies starting case from 1 to 108 
% endingPair       specifies ending case from startingPair to 108 
% nExp             # of experiments per  xyzF & abcN  initial condition
% M                (true,false) => (use memory, start from scratch)
% BS               (true,false) => (bootstrap, split by division)
% deletePrevious   (true,false) => delete previous file or append to file
% outputScopeOnly  (true,false) => show the plan  OR  run the experiments
% showT            (true,false) => (show, no show) classification results
% show_vb          1,2,3                  => different levels of verbosity
% minGap           minimum gap in log probability between F/N molecules
% ic               (best PCA, identity)-basis initial conditions
%
% INTERNAL
% minGap = the gap between largest -log10(probability of nonfunctional) 
%                        - smallest -log10(probability of functional)
% minGap = minGap10/10;
%                                   
% 
% OUTPUT
%                                           many output files of the form:
%   xyz_abc_jkSout.dlm
%             ^------------> S or B for small or big system size
%            ^-------------> minGap10 = 0,1,2, ... 99  (two digits)
%           ^--------------> ic       = 0,1    
%         ^----------------> c        = F, L, T
%        ^-----------------> b        = F, L, T, S
%       ^------------------> a        = E, F
%     ^--------------------> z        = F, L, T
%    ^---------------------> y        = F, L, T, S
%   ^----------------------> x        = E, F
%
% examples: EST_EFF_120Sout.dlm
%           FLF_EFT_025Bout.dlm
%           ETT_FFF_130Sout.dlm
%           FTL_FSL_015Bout.dlm
%
% PROCESS overview:
% There are 24 molecules of the form xyz
% x = E or F
% y = L or S or T or F
% z = L or T or F
% ------------------------------------------------------ summary of script
%   for y = L, S, T, F
%      for x = E or F
%         for z = L or T or F
%         DEFINE  xyz                               => functional molecule
% --------------------------------------- get first nonfunctional molecule
%            for b = L, S, T, F
%               if( b ~= y )
%                  for a = E or F
%                     for c = L or T or F
%                     DEFINE abc                 => nonfunctional molecule
%                     USE xyz and abc as INITIAL CONDITION
%                     Initialize RECORD 
%                        for trials=1:25 
%                        automate SPLOC
%                        UPDATE RECORD by tracking path
%                        break out after 6 functional are found
%                        end
%                     RECORD all data for 26 trials
%                     end
%                  end
%               end
%            end
%         end
%      end
%   end
% ----------------------------------------------------------------- counts
% [(4*2*3)* (3*2*3)]* 25 = 432*25 = 10800 runs
%                          432*3 = 1296       takes 25 days on my computer
% Each run is saved 
% summarize all runs, use to report statistics and summary graph.
%%                                              initialize data structures
xyzLabels = cell(4, 2, 3 );
numLabels = cell(4, 2, 3 );
%                kL,kM,kR
%                xyz =>  kM -- kL -- kR
% ------------------------------------------------------
xyzLabels{1,1,1} = {'FLF'};   numLabels{1,1,1} = {'1'};
xyzLabels{1,1,2} = {'FLL'};   numLabels{1,1,2} = {'2'};
xyzLabels{1,1,3} = {'FLT'};   numLabels{1,1,3} = {'3'};
% ------------------------------------------------------
xyzLabels{1,2,1} = {'ELF'};   numLabels{1,2,1} = {'4'};
xyzLabels{1,2,2} = {'ELL'};   numLabels{1,2,2} = {'5'};
xyzLabels{1,2,3} = {'ELT'};   numLabels{1,2,3} = {'6'};
% =======================================================
xyzLabels{2,1,1} = {'FSF'};   numLabels{2,1,1} = {'7'};
xyzLabels{2,1,2} = {'FSL'};   numLabels{2,1,2} = {'8'};
xyzLabels{2,1,3} = {'FST'};   numLabels{2,1,3} = {'9'};
% ------------------------------------------------------
xyzLabels{2,2,1} = {'ESF'};   numLabels{2,2,1} = {'10'};
xyzLabels{2,2,2} = {'ESL'};   numLabels{2,2,2} = {'11'};
xyzLabels{2,2,3} = {'EST'};   numLabels{2,2,3} = {'12'};
% =======================================================
xyzLabels{3,1,1} = {'FTF'};   numLabels{3,1,1} = {'13'};
xyzLabels{3,1,2} = {'FTL'};   numLabels{3,1,2} = {'14'};
xyzLabels{3,1,3} = {'FTT'};   numLabels{3,1,3} = {'15'};
% ------------------------------------------------------
xyzLabels{3,2,1} = {'ETF'};   numLabels{3,2,1} = {'16'};
xyzLabels{3,2,2} = {'ETL'};   numLabels{3,2,2} = {'17'};
xyzLabels{3,2,3} = {'ETT'};   numLabels{3,2,3} = {'18'};
% =======================================================
xyzLabels{4,1,1} = {'FFF'};   numLabels{4,1,1} = {'19'};
xyzLabels{4,1,2} = {'FFL'};   numLabels{4,1,2} = {'20'};
xyzLabels{4,1,3} = {'FFT'};   numLabels{4,1,3} = {'21'};
% ------------------------------------------------------
xyzLabels(4,2,1) = {'EFF'};   numLabels{4,2,1} = {'22'};
xyzLabels(4,2,2) = {'EFL'};   numLabels{4,2,2} = {'23'};
xyzLabels(4,2,3) = {'EFT'};   numLabels{4,2,3} = {'24'};
% ------------------------------------------------------------- # --> name
num2Name = cell(1,24);
triplet2num = zeros(4,2,3);
m = 0;
   for jL=1:4
      for jM=1:2
         for jR=1:3
         msg1 = char( xyzLabels{jL,jM,jR} );
         num1 = str2num( char( numLabels{jL,jM,jR} ) );
         m = m + 1;
         num2Name{m} = msg1;
         triplet2num(jL,jM,jR) = m;
         %disp([num2str(m),'  ',num2str(num1),'  ',msg1]);
         end
      end
   end
%%                                                 % functional ID numbers
FID = cell(1,4);                               % molecules of the form xyz
FID{1} = [1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];     % => y = L
FID{2} = [0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];     % => y = S
FID{3} = [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0];     % => y = T
FID{4} = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1];     % => y = F
masterList = zeros(1,24); %=> -1, 0, 1 => nonfunctional,unknown,functional
% ---------------------------------------------------------- sanity checks
   if( (jL_bucket < 1) || (jL_bucket > 4) )
   error('jL_bucket is out of range: can use {1, 2, 3, 4}');
   end
   if( startingPair > 108 )
   error('startingPair is out of range: Max value = 108');
   end
   if( startingPair < 1 )
   error('startingPair is out of range: Min value = 1');
   end
   if( endingPair < startingPair)
   disp('   ');
   disp('   ');
   disp(['startingPair = ',num2str(startingPair)]);
   error('endingPair must be equal to or greater than startingPair');
   end
   if( endingPair > 108 )
   error('endingPair is out of range: Max value = 108');
   end
   if( nExp < 1 )
   nExp = 1;
   end
   if( nExp > 25 )
   error('Do in batches, nExp should not be > 25');
   end
nExp = round(nExp);
% -------------------------------------------------------- restrict minGap
minGap10 = round(10*minGap);              % round to nearest decimal point
   if( minGap10 < 0 )
   minGap10 = 0;
   end
   if( minGap10 > 99 )
   minGap10 = 99;
   end
minGap = minGap10/10;
% ------------------------------------------------------------ restrict ic
ic = round(ic);
   if( (ic < 0) || (ic > 1) )
   error('ic is out of range: ic = {0,1}');
   end
% ---------------------------------------------------------- set constants
WRONG = false;
CORRECT = true;
%%                                        initialize SPLOC data collectors
% REMARK: Requires having the training data prepared. This data comes with
%         the sploc toolbox. Similar to workflow01 and workflow02
%
% ------------------------------------------------- determine training set
logFile = ['all_',mfilename];
initializeSPLOC(1,'fName',logFile,'gType','png');
mFormat = setDataMatrixFormat('xxx-yyy-zzz',2);  % by-pass screen input
mType = 'cov';
  if( nS == 500 )
  [~,~,FnameAdataU] = readFileNameList('trainingSet1',0,'cType','U');
  else
  [~,~,FnameAdataU] = readFileNameList('trainingSet1B',0,'cType','U');
  end
FnameAdataX = FnameAdataU;
% --- put file names in sync with the numbering system used in this script
   for k=1:24
   xyz = FnameAdataU{k}(end-2:end);
      switch xyz(1)
          case 'F'
          kM = 1;
          case 'E'
          kM = 2;
          otherwise
          error('unknown first letter in the xyz format');
      end
% ------------------------------------------------------------------------
      switch xyz(2)
          case 'L'
          kL = 1;
          case 'S'
          kL = 2;
          case 'T'
          kL = 3;
          case 'F'
          kL = 4;
          otherwise
          error('unknown first letter in the xyz format');
      end
% ------------------------------------------------------------------------
      switch xyz(3)
          case 'F'
          kR = 1;
          case 'L'
          kR = 2;
          case 'T'
          kR = 3;
          otherwise
          error('unknown first letter in the xyz format');
      end
% ------------------------------------------------------------------------
   j = triplet2num(kL,kM,kR);
      if( strcmp( num2Name{j} , xyz ) )
      FnameAdataX{j} = FnameAdataU{k};
      end   
   end
% ------------------------------------------------------------------ check
%    for j=1:24
%    disp( [num2str(j),'  ',FnameAdataX{j}] );
%    end
%%                                                  create reference trait
refNameX = 'refX';
[AmatrixInfoX,tableX] = readDataMatrices(refNameX,FnameAdataX,mFormat);
nSamples = AmatrixInfoX.nSamples(1);
ns2 = round(nSamples/2);
%disp(tableX);                                                % checks out
% ------------------------------------------------------------------ check
% traitX2 = getBoostedMVStats4sploc(AmatrixInfoX,1,0,mType);
% traitX2
% disp('   ');
%    for j=1:24
%    j2b = 2*j;
%    j2a = j2b - 1; 
%    disp( [num2str(j),'  ',char(traitX2.mMatrixName{j2a}), ...
%                      '  ',char(traitX2.mMatrixName{j2b})] );
%    end
% % --------------
% disp('   ');
%    for j=1:24
%    j2b = 2*j;
%    j2a = j2b - 1; 
%    disp( [num2str(j),'  ',char(traitX2.cMatrixName{j2a}), ...
%                      '  ',char(traitX2.cMatrixName{j2b})] );
%    end
% --------------------------------- initialize containers (to be modified)
AmatrixInfoF = AmatrixInfoX;
AmatrixInfoN = AmatrixInfoX;
traitX = getMultivariateStats4sploc(AmatrixInfoX,nSamples,0,mType);
%        all samples are being used, so use getMultivariateStats4sploc
%traitX             
idCode = [num2str(ic),num2str(minGap10)];
   if( showT )
   vb = show_vb;
   else
   vb = 0;
   end
%%                                                    initialize variables
% ------------------------------------------------------- output container
nRows = 1 + 7 + 4 + 1 + 24;
%       |   |   |   |   ^^--> # of molecules
%       |   |   |   ^-------> selected molecule (1 to 24, canonical order)
%       |   |   ^-----------> counts for true/false -- positives/negatives
%       |   ^---------------> train, Dd,Di, minPF, maxPN, gap, #remaining
%       |-------------------> label for model prediction #
% ------------------------------------------------------------------------
nCols = 1 + 22;
%       |   ^^-------------> maximum # of experiments pass the initial two
%       ^------------------> trial # (simulated series of experiments)
saveWrite = zeros(nRows,nCols);
sIndex = zeros(1,24);
mType = 'cov';
nPair = 0;               % nPair = counter for distinct starting F/N pairs
%%                                         turn 6 nexted loops into 1 loop
save_splocOut = cell(1,108);
save_prefix   = cell(1,108);
save_xyzF     = cell(1,108);
save_numF     = cell(1,108);
save_abcN     = cell(1,108);
save_numN     = cell(1,108);
save_nPair    = ones(1,108);
save_j        = ones(1,108);
save_k        = ones(1,108);
% -------------------- push through all 6 loops without sploc calculations
kCases = 0;                                   % counts actual cases to run
   for jL=jL_bucket:jL_bucket                   % start of MAIN OUTER LOOP
% ------------------------------------------------------- initialize SPLOC
      switch jL
          case 1
          logFile = ['xLz_',mfilename];
          case 2
          logFile = ['xSz_',mfilename];
          case 3
          logFile = ['xTz_',mfilename];
          case 4
          logFile = ['xFz_',mfilename];
      end
% ------------------------------------------------------- initialize SPLOC
   initializeSPLOC(1,'fName',logFile,'gType','png');
   %disp('   ');
   %disp( dividerLine(logFile) );
% ------------------------------------------------------------------------
      for jM=1:2
         for jR=1:3
         xyzF = char( xyzLabels{jL,jM,jR} );
         numF = str2num( char( numLabels{jL,jM,jR} ) );
         j = triplet2num(jL,jM,jR);
% ---------------------------------
            for kL=1:4
               for kM=1:2
                 if( kL ~= jL )
                   for kR=1:3
                   nPair = nPair + 1;
                     if( (nPair >= startingPair) && ...
                         (nPair <= endingPair) )
                     k = triplet2num(kL,kM,kR);
                     abcN = char( xyzLabels{kL,kM,kR} );
                     numN = str2num( char( numLabels{kL,kM,kR} ) );
                     
                     %disp('   ');
                     %disp(['--------------------------------------', ...
                     %      '--- nPair = ',num2str(nPair), ...
                     %      '  (',num2str(numF),',',num2str(numN), ...
                     %      ')  => F:',xyzF,'  N:',abcN]);
                     %disp( num2Name{j} );
                     %disp( num2Name{k} );
                     %disp( [j,k] );
                     
% ----------------------------------------------- perform nExp experiments
                        if( nS == 500 )
                        prefix = [xyzF,'_',abcN,'_',idCode,'S'];
                        else
                        prefix = [xyzF,'_',abcN,'_',idCode,'B'];
                        end
                     splocOut = [prefix,'out.dlm'];
                     
                     disp( splocOut );
                     
                        if( deletePrevious )
                           if( isfile(splocOut) )
                           delete(splocOut);
                           end
                        end
                        
                     kCases = kCases + 1;
                     save_splocOut{kCases} = splocOut;
                     save_prefix{kCases}   = prefix;
                     save_abcN{kCases}     = abcN;
                     save_numN{kCases}     = numN;
                     save_numF{kCases}     = numF;
                     save_xyzF{kCases}     = xyzF;
                     save_nPair(kCases)    = nPair;  
                     save_j(kCases)        = j;
                     save_k(kCases)        = k;
% --------------------------------------------------- run experiment check
                        for r=1:nExp                       % r => repeat #
                        prefixExp = [prefix,'Exp',num2str(r,'%02i'),'_'];
                        %disp('  ');
                        %msg = ['trial ',num2str(r,'%02i')];
                        %disp( dividerLine(msg) );
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         for t=2:23                                % t-th training session
         splocName = [prefixExp,num2str(t,'%02i')]; 
%          disp( splocName )
%          strF = char( traitF2.mMatrixName );
%          strN = char( traitN2.mMatrixName );
%          disp('   ');
%          disp( dividerLine('Functional') );
%          disp( strF );
%          disp( dividerLine('Nonfunctional') );
%          disp( strN );
         splocNameLog = [splocName,'_sploc.log'];
         %disp( splocNameLog );
         fileLog2Remove = getOutputFileName('training',splocNameLog);
            if( isfile(fileLog2Remove) )
            delete(fileLog2Remove);
            end
         end                      
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        end
                     end
                   end
                 end
               end
            end
         end
      end
   end
   
% % % %   if( jL > 0 ) 
% % % %   error('stop here for now');
% % % %   end

% ---------------------------------------- run 6 loops using a single loop
jL = jL_bucket;                     % this variable is needed for a bucket
   if( outputScopeOnly )        % use random predictions in place of SPLOC
      for kC=1:kCases
      xyzF = save_xyzF{kC};
      numF = save_numF{kC};
      splocOut = save_splocOut{kC};
      prefix = save_prefix{kC};
      abcN = save_abcN{kC};
      numN = save_numN{kC};              
      nPair = save_nPair(kC);
      j = save_j(kC);
      k = save_k(kC);
% ---------------------------------- at this point, emulate origional code
                     r1Tag = 1;
                        for r=1:nExp                       % r => repeat #
                        prefixExp = [prefix,'Exp',num2str(r,'%02i'),'_'];
                        %disp('  ');
                        %msg = ['trial ',num2str(r,'%02i')];
                        %disp( dividerLine(msg) );
                        masterList = 0*masterList;     % => all is unknown
                        % ------------------------------------------------
                        masterList(j) = 1;             %=> j is functional
                        nF = 1;
                        AmatrixInfoF.n = nF;
                        AmatrixInfoF.dataMatrixName{nF} ...
                                       = AmatrixInfoX.dataMatrixName{j};
                        AmatrixInfoF.A{nF} = AmatrixInfoX.A{j};
                           if( BS )
                           traitF2 = getBoostedMVStats4sploc( ...
                                     AmatrixInfoF,1,0,mType);
                           else
                           traitF2 = getMultivariateStats4sploc( ...
                                     AmatrixInfoF,ns2,0,mType);
                           end
                        % ------------------------------------------------
                        masterList(k) = -1;         %=> k is nonfunctional
                        nN = 1;
                        AmatrixInfoN.n = nN;
                        AmatrixInfoN.dataMatrixName{nN} ...
                                       = AmatrixInfoX.dataMatrixName{k};
                        AmatrixInfoN.A{nN} = AmatrixInfoX.A{k};
                           if( BS )
                           traitN2 = getBoostedMVStats4sploc( ...
                                     AmatrixInfoN,1,0,mType);
                           else
                        traitN2 = getMultivariateStats4sploc( ...
                                  AmatrixInfoN,ns2,0,mType);
                           end
% -------------------------------------------------------- run experiments
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         for t=2:23                                % t-th training session
         splocName = [prefixExp,num2str(t,'%02i')]; 
%          disp( splocName )
         strF = char( traitF2.mMatrixName );
         strN = char( traitN2.mMatrixName );
%          disp('   ');
%          disp( dividerLine('Functional') );
%          disp( strF );
%          disp( dividerLine('Nonfunctional') );
%          disp( strN );
         splocNameLog = [splocName,'_sploc.log'];
%          disp( splocNameLog )
%          pause
         prb = rand(1,24);                                      % unsorted
         log10invprb = -log10( max(1.0e-20,prb) );
         gap = 0;
         %disp('>>>>>random generated probabilities');           
         fp = prb.*( 1 - 2*abs(masterList) );
         [~,jNext] = max(fp);
            if( r == r1Tag )
            disp('  ');
            disp(['------------- minGap= ',num2str(minGap), ...
                               ' ic= ',num2str(ic), ...
                               ' nExp= ',num2str(nExp) ...
                               ' bucket= ',num2str(jL_bucket), ...
                               ' nPair= ',num2str(nPair), ...
                               ' (',num2str(numF),',',num2str(numN), ...
                               ')  => F:',xyzF,'  N:',abcN]);
            disp('   ');
            r1Tag = 0;
            end
% ------------------------------------ decide F or N for next training set
           if( FID{jL}(jNext) == 1 )                       % => functional
           %disp('  ');
           %disp(['F: ', char( num2Name{jNext} ) ]);
           % -------------------------------------------------------------
           nF = nF + 1;
           masterList(jNext) = 1;                  %=> jNext is functional
           AmatrixInfoF.n = nF;
           AmatrixInfoF.dataMatrixName{nF} ...
                          = AmatrixInfoX.dataMatrixName{jNext};
           AmatrixInfoF.A{nF} = AmatrixInfoX.A{jNext};
              if( BS )
              traitF2 = getBoostedMVStats4sploc(AmatrixInfoF,1,0,mType);
              else
              traitF2 = ... 
                   getMultivariateStats4sploc(AmatrixInfoF,ns2,0,mType);
              end
           else                                         % => nonfunctional
           %disp('  ');
           %disp(['N: ', char( num2Name{jNext} ) ]);
           % -------------------------------------------------------------
           nN = nN + 1;
           masterList(jNext) = -1;              %=> jNext is nonfunctional
           AmatrixInfoN.n = nN;
           AmatrixInfoN.dataMatrixName{nN} ...
                          = AmatrixInfoX.dataMatrixName{jNext};
           AmatrixInfoN.A{nN} = AmatrixInfoX.A{jNext};
              if( BS )
              traitN2 = getBoostedMVStats4sploc(AmatrixInfoN,1,0,mType);
              else
              traitN2 = ...
                   getMultivariateStats4sploc(AmatrixInfoN,ns2,0,mType);
              end
           end
% --------------------- break out early when nF = 6 and nN > 2 and gap > 3
           %disp('   ');
           %disp(['nF = ',num2str(nF),'  nN = ',num2str(nN), ...
           %      '  FNgap = ',num2str(gap)]);
           if( (nF == 6) && (gap > minGap) )              % => break early
           break;
           end
         end                           % end of t = 2:23 training sessions
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        end          % end of r from 1 to nExp experiments
      
      
      end
   else                                % run 1 loop for SPLOC calculations
   useSBV = M;
      for kC=1:kCases
      %xyzF = save_xyzF{kC};
      %numF = save_numF{kC};
      splocOut = save_splocOut{kC};
      prefix = save_prefix{kC};
      %abcN = save_abcN{kC};
      %numN = save_numN{kC};             
      %nPair = save_nPair(kC);
      j = save_j(kC);
      k = save_k(kC);
% ---------------------------------- at this point, emulate origional code                  
                     %r1Tag = 1;
                        for r=1:nExp                       % r => repeat #
                        prefixExp = [prefix,'Exp',num2str(r,'%02i'),'_'];
                        %disp('  ');
                        %msg = ['trial ',num2str(r,'%02i')];
                        %disp( dividerLine(msg) );
                        masterList = 0*masterList;     % => all is unknown
                        % ------------------------------------------------
                        masterList(j) = 1;             %=> j is functional
                        nF = 1;
                        AmatrixInfoF.n = nF;
                        AmatrixInfoF.dataMatrixName{nF} ...
                                       = AmatrixInfoX.dataMatrixName{j};
                        AmatrixInfoF.A{nF} = AmatrixInfoX.A{j};
                           if( BS )
                           traitF2 = getBoostedMVStats4sploc( ...
                                     AmatrixInfoF,1,0,mType);
                           else
                           traitF2 = getMultivariateStats4sploc( ...
                                     AmatrixInfoF,ns2,0,mType);
                           end
                        % ------------------------------------------------
                        masterList(k) = -1;         %=> k is nonfunctional
                        nN = 1;
                        AmatrixInfoN.n = nN;
                        AmatrixInfoN.dataMatrixName{nN} ...
                                       = AmatrixInfoX.dataMatrixName{k};
                        AmatrixInfoN.A{nN} = AmatrixInfoX.A{k};
                           if( BS )
                           traitN2 = getBoostedMVStats4sploc( ...
                                     AmatrixInfoN,1,0,mType);
                           else
                           traitN2 = getMultivariateStats4sploc( ...
                                     AmatrixInfoN,ns2,0,mType);
                           end
% -------------------------------------------------------- run experiments
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     saveWrite = 0*saveWrite;
     enthusiasm = true;       % always start a new project with enthusiasm
     contextF = traitF2;
     contextN = traitN2;
     prediction = false;            % initialize no prediction => as false
     gap = 1000;
     countTP = 0;     
     countTN = 0;     
     countFP = 0;
     countFN = 0;
     tLast = 15;
        for t=2:tLast                              % t-th training session
        splocName = [prefixExp,num2str(t,'%02i')]; 
        disp( splocName )
        strF = char( traitF2.mMatrixName );
        strN = char( traitN2.mMatrixName );
        disp('   ');
        disp( dividerLine('Functional') );
        disp( strF );
        disp( dividerLine('Nonfunctional') );
        disp( strN );
        %splocNameLog = [splocName,'_sploc.log'];
        %pause
% ------------------------------ update context upon false prediction only
% ---------------------------------------------------- make new prediction
           if( t == 2 )                                         % SPLOC it
           splocResult = sploc(0,ic,splocName,traitF2,traitN2,vb); 
           SBV = splocResult.SBV;
           Di = splocResult.Di;
           Dd = splocResult.Dd;
           % -------------------------------------- set initial perception
              if( Dd == 0 )
              Up = splocResult.SBV(:,1);   
              else
              Up = getDiscriminantSBV(splocResult); 
              end
           elseif( enthusiasm && (prediction == WRONG) ) % => do an update
           %clear splocResult                % nothing seems to fix MATLAB
              if( useSBV )
                if( gap > 0.5 )
                %disp('SPLOC using SBV memory');
                splocResult = sploc(0,SBV,splocName,traitF2,traitN2,vb);
                                    % ^^^------------> using prior memory
                else
                %disp('SPLOC using no memory');
                splocResult = sploc(0,0,splocName,traitF2,traitN2,vb);
                end                 % ^-> erase memory: start from scratch
              else
              %disp('SPLOC using no memory');
              splocResult = sploc(0,0,splocName,traitF2,traitN2,vb);
              end                 % ^---> erase memory: start from scratch
           Di = splocResult.Di;
           Dd = splocResult.Dd;
% -------------------------------------- check if new perception is viable
             if( Dd > 0 )                              % update perception
             contextF = traitF2;                             % new context
             contextN = traitN2;                             % new context
             Up = getDiscriminantSBV(splocResult);        % new perception
             SBV = splocResult.SBV;                    % <= retains memory
             else                            % failure after three strikes 
             enthusiasm = false;   % no hope without discriminant subspace 
             end
           %else   previous prediction is correct: keep current perception
           end
% ------------------- make predictions with current perception and context
        rankName = [prefix,'rankX',num2str(t,'%02i')];
        %[rankX,T] = classifyDataStream(rankName,AmatrixInfoX, ... 
        %                              AmatrixInfoF,AmatrixInfoN,Ud,0);
        [~,T] = classifyMoments(rankName,traitX,contextF,contextN,Up,0);
% --- put file names in sync with the numbering system used in this script
           for kk=1:24
           uvw = T.(2){kk}(end-2:end);
              switch uvw(1)
              case 'F'
              kkM = 1;
              case 'E'
              kkM = 2;
              otherwise
              error('unknown first letter in the xyz format');
              end
% ------------------------------------------------------------------------
              switch uvw(2)
              case 'L'
              kkL = 1;
              case 'S'
              kkL = 2;
              case 'T'
              kkL = 3;
              case 'F'
              kkL = 4;
              otherwise
              error('unknown first letter in the xyz format');
              end
% ------------------------------------------------------------------------
              switch uvw(3)
              case 'F'
              kkR = 1;
              case 'L'
              kkR = 2;
              case 'T'
              kkR = 3;
              otherwise
              error('unknown first letter in the xyz format');
              end
% ------------------------------------------------------------------------
           jj = triplet2num(kkL,kkM,kkR);
              if( strcmp( num2Name{jj} , uvw ) )
              sIndex(jj) = kk;
              end   
           end
        prbSorted = T.(3)';                                    
        log10invprbSorted = T.(5)';                           
        prb = prbSorted(sIndex);                         % cononical order
        log10invprb = log10invprbSorted(sIndex);         % cononical order
       %^^^^^^^^^^^------------> = negative of log10(prob) = log10(1/prob)
% ----------------------------------------------------- verification check
        %disp(T);
%               for jj=1:24
%               kk = sIndex(jj);
%               msg = [char( T.(2){kk} ),'  ',num2str( prb(jj) ), ...
%                      '  ',num2str( log10invprb(jj) )];
%               disp(msg);
%               end
        %pause
% --------------------------------------------------------------- analysis 
        LF = ( masterList > 0  );
        LN = ( masterList < -0.5 );
        gap = min( log10invprb(LN) ) - max( log10invprb(LF) ); 
        minPF = min( prb(LF) );            % minimum prob to be functional
        maxPN = min( prb(LN) );         % maximum prob to be nonfunctional
% --------------------------------------------------
        LU = ( abs(masterList) < 0.5 );
        pCut = minPF/10^minGap;      % threshold is based on minGap
        countRemaining = sum( prb(LU) > pCut ); 
        fp = prb.*( 1 - 2*abs(masterList) );
        [~,jNext] = max(fp);
        pSelect = prb(jNext);
% ------------------------------------ decide F or N for next training set
           if( FID{jL}(jNext) == 1 )                       % => functional
           %disp('  ');
           %disp(['F: ', char( num2Name{jNext} ) ]);
           % -------------------------------------------------------------
           nF = nF + 1;
           masterList(jNext) = 1;                  %=> jNext is functional
           AmatrixInfoF.n = nF;
           AmatrixInfoF.dataMatrixName{nF} ...
                          = AmatrixInfoX.dataMatrixName{jNext};
           AmatrixInfoF.A{nF} = AmatrixInfoX.A{jNext};
              if( BS )
              traitF2 = getBoostedMVStats4sploc(AmatrixInfoF,1,0,mType);
              else
              traitF2 = ...
                   getMultivariateStats4sploc(AmatrixInfoF,ns2,0,mType);
              end
% ------------------- apply conditions for true positive or false negative
              if( pSelect > minPF )                        % => functional
              prediction = CORRECT;                     % => true positive
              countTP = countTP + 1;
              elseif( pSelect < maxPN )
              prediction = WRONG;                      % => false negative
              countFN = countFN + 1;
              elseif( pSelect/minPF > 1/10^minGap )        % => functional
              prediction = CORRECT;                     % => true positive
              countTP = countTP + 1;
              else                      % the model predicts nonfunctional
              prediction = WRONG;                      % => false negative
              countFN = countFN + 1;
              end
           else                                         % => nonfunctional
           %disp('  ');
           %disp(['N: ', char( num2Name{jNext} ) ]);
           % -------------------------------------------------------------
           nN = nN + 1;
           masterList(jNext) = -1;              %=> jNext is nonfunctional
           AmatrixInfoN.n = nN;
           AmatrixInfoN.dataMatrixName{nN} ...
                          = AmatrixInfoX.dataMatrixName{jNext};
           AmatrixInfoN.A{nN} = AmatrixInfoX.A{jNext};
              if( BS )
              traitN2 = getBoostedMVStats4sploc(AmatrixInfoN,1,0,mType);
              else
              traitN2 = ...
                   getMultivariateStats4sploc(AmatrixInfoN,ns2,0,mType);
              end
% ------------------- apply conditions for true negative or false positive
              if( pSelect > minPF )                        % => functional
              prediction = WRONG;                      % => false positive
              countFP = countFP + 1;
              elseif( pSelect < maxPN )
              prediction = CORRECT;                     % => true negative
              countTN = countTN + 1;
              elseif( pSelect/minPF > 1/10^minGap )        % => functional
              prediction = WRONG;                      % => false positive
              countFP = countFP + 1;
              else                      % the model predicts nonfunctional
              prediction = CORRECT;                     % => true negative
              countTN = countTN + 1;
              end 
           end
% --------------------------- display results of experimental verification
           if( showT )
           disp('  ');
           disp(T);
           disp(['minPF= ',num2str( round(1000*minPF)/1000 ), ...
                 '  maxPN= ',num2str(maxPN), ...
                 '  gap= ',num2str(gap), ...
                 '  #remaining = ',num2str(countRemaining), ...
                 '  TP= ',num2str(countTP), ...
                 '  TN= ',num2str(countTN), ...
                 '  FP= ',num2str(countFP), ...
                 '  FN= ',num2str(countFN)]);
           disp(dividerLine);
           disp('  ');
           end
% ------------------------------ save results of experimental verification
% output container definition
% nRows = 1 + 7 + 4 + 1 + 24;
% %       |   |   |   |   ^^> # of molecules
% %       |   |   |   ^-----> selected molecule (1 to 24, canonical order)
% %       |   |   ^----------> counts for true/false - positives/negatives
% %       |   ^--------------> train, Dd,Di, minPF, maxPN, gap, #remaining
% %       |------------------> label for model prediction #
% % ----------------------------------------------------------------------
% nCols = 1 + 22;
% %       |   ^^-----------> maximum # of experiments pass the initial two
% %       ^----------------> trial # (simulated series of experiments) 
% ------------------------------------------------------------------------
         nTrial = r*ones(nRows,1);
         nPred = t - 1; 
         saveWrite(:,1) = nTrial;                     % all rows  column 1
% ------------------------------------------          1 row,   the 1st row
         saveWrite(1,t) = nPred;                     
% ------------------------------------------          7 rows:  2 through 8
            if( prediction )
            xtrain = 0;                            % => no need to retrain
            else
               if( enthusiasm )
               xtrain = 1;            % => retrain after a mistake is made
               else
               xtrain = -1;                   % should retrain, but do not
               end
            
            end
         saveWrite(2,t) = xtrain;
         saveWrite(3,t) = Dd;
         saveWrite(4,t) = Di;
         saveWrite(5,t) = minPF;
         saveWrite(6,t) = maxPN;
         saveWrite(7,t) = gap;                       
         saveWrite(8,t) = countRemaining;
% ------------------------------------------          4 rows: 9 through 12
         saveWrite(9,t) = countTP;
         saveWrite(10,t)= countTN;
         saveWrite(11,t)= countFP;
         saveWrite(12,t)= countFN;
% ------------------------------------------           1 row: the 13th
         saveWrite(13,t)= jNext;
% ------------------------------------------          24 rows: 14th to end
         saveWrite(14:end,t) = log10invprb';  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Go or no go decision:             criteria to stop prediction/experiment
% If go => make a new prediction followed by experimental verification
% REMARK: To get the data we need for the paper, instead of breaking 
% prematurely, the break is prevented until all 6 functional molecules are
% identified. However, in practice an early break would mean functional
% molecules will be missed. 
           if( countRemaining == 0 )      % must check all active unknowns
              if( nF == 6 )  
              break;
              end
           end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %pause                             % for slow manual observations
         end                           % end of t = 2:23 training sessions
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                        fileID = fopen(splocOut,'a');
                        outputFormat = ['%d ', ...
                            '%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f ', ...
                            '%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f ', ...
                            '%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f ', ...
                            '%6.3f %6.3f %6.3f %6.3f \n'];
                           for iL=1:nRows
                           fprintf(fileID,outputFormat,saveWrite(iL,:)); 
                           end
                        fclose(fileID);
                        end          % end of r from 1 to nExp experiments
      end
   end
end