% test driver for classification functions within the SPLOC toolset
% ------------------------------------------------------------------------
% SPLOC = Supervised Projective Learning with Orthogonal Completeness
%
% This test driver checks sploc functions related to classification that
% employs the discriminant selection basis vectors obtained by splocing.  
% 
% In addition to testing these functions, the purpose of this test driver
% is to serve as an example for how to use these sploc functions. 
%
% Toolset: I/O directory structure:
% 
% current directory --> sub-directories: input
%                                        splocLog
%                                        training
%                                        basisComparison
%                                        classification
%                                        analysis
%                 
% INPUT: (general description) 
% Many types of files can be read by a variety of sploc functions. Each
% such file must be formated according to specifications. In many cases, 
% the output files generated by some sploc functions will serve as input
% files for other sploc functions, which will usually have dependence. For
% example, under supervised training on various datastreams that are known
% to be functional or nonfunctional, a basis set of vectors is optimized
% to successfully classify an unknown datastream. At a later time, new
% data that comes in can be processed simply by reading in an input file
% that contains the relevant basis vectors for discrimination, without 
% redoing the training. Hence, the output from training becomes the input
% for classification along with more data. The read functions are designed
% to have flexibility so that no specific rules need to be followed.
% 
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
% 
% PROCESS: 
% SPLOC is comprised of many MATLAB functions that define the SPLOCtoolset
% See splocToolsetVersion() for detail comments about the SPLOC toolset.
%
% Each function that is tested here is listed in the menu.
%%                                                     start from sctratch
clc 
disp('  ');
disp(mfilename);
disp( dividerLine('apply new setup') );
iii = input('enter 1 to create new setup: ');
   if( iii == 1 )               % cleans out memory usage and starts fresh
   clear all
   iii = 1;
   else
   iii = 0;
   end
close all
%%                                                           perform setup
if( iii == 1 )
% ------------------------------------------------- determine training set
   disp('  ');
   disp( dividerLine('select training set # (1 through 6)') );
   setNumber = input('Enter set #: ');
      if( setNumber < 1 )
      error('set # must not be less than 1: Available range 1-6');
      end
      if( setNumber > 6 )
      error('set # must not be greater than 6: Available range 1-6');
      end
      switch setNumber
       case 1
       fName1a = 'trainingSet1';
       fName1b = 'trainingSet1B';
       case 2
       fName1a = 'trainingSet2';
       fName1b = 'trainingSet2B';
       case 3
       fName1a = 'trainingSet3';
       fName1b = 'trainingSet3B';
       case 4
       fName1a = 'trainingSet4';
       fName1b = 'trainingSet4B';
       case 5
       fName1a = 'trainingSet5';
       fName1b = 'trainingSet5B';
       case 6
       fName1a = 'trainingSet6';
       fName1b = 'trainingSet6B';
       otherwise
       error('unknown option');
      end
   fName1 = fName1a;
% ----------------------------- determine covariance or correlation matrix
   disp('  ');
   disp( dividerLine('select covariance or correlation matrix') );
   disp('1. use covariance  matrix')
   disp('2. use correlation matrix')
   iMtype = input('   Enter option: ');
   switch iMtype
       case 1
       mType = 'cov';
       case 2
       mType = 'cor';
       otherwise
       error('unknown option');
   end
end
%%                                                     build function menu
test_classifyMoments = false;
test_classifyDataStream = false;
test_plotRankDecomposition = false;
test_classifyMoments_switch = false;
test_classifyDataStream_switch = false;

disp('  ');
disp( dividerLine('select function to test') );
disp('1. classifyMoments() ')
disp('2. classifyDataStream() ')
disp('3. plotRankDecomposition');
disp('4. classifyMoments() with switch: F -> N and N -> F')
disp('5. classifyDataStream() with switch: F -> N and N -> F')
disp('  ');
nMenu = input('   Enter option: ');
   switch nMenu
       case 1
       test_classifyMoments = true;
       case 2
       test_classifyDataStream = true;
       case 3
       test_plotRankDecomposition = true;
       case 4
       test_classifyMoments_switch = true;
       case 5
       test_classifyDataStream_switch = true;
       otherwise
       error('unknown option');
   end
prefix = ['opt',num2str(nMenu,'%02i'),mType,'_'];
disp('  ');
disp( dividerLine(['setup = ',prefix]) );
%%                                                        initialize SPLOC
   if( iii == 1 )
   logFile = [prefix,mfilename];
   initializeSPLOC(1,'fName',logFile,'gType','png');
   mFormat = setDataMatrixFormat('xxx-yyy-zzz',2);  % by-pass screen input
% =========================================================== get traits1a
   disp('   ');
   disp( dividerLine(['constructing traits from: ',fName1a]) );
% ------------------------------------------------------------------------
   [~,~,FnameAdata1aF] = readFileNameList(fName1a,0,'sType','F');
   [~,~,FnameAdata1aN] = readFileNameList(fName1a,0,'sType','N');
   [~,~,FnameAdata1aX] = readFileNameList(fName1a,0);           % X => FUN
% ------------------------------------------------------------------------
   refName1aF = [fName1a,'_F'];           % make this a suitable file name
   [AmatrixInfo1aF,table1aF] = readDataMatrices(refName1aF, ... 
                                                FnameAdata1aF,mFormat);
   disp('  ');
   disp(table1aF);
   disp('   ');
% -----------------------------
   refName1aN = [fName1a,'_N'];           % make this a suitable file name
   [AmatrixInfo1aN,table1aN] = readDataMatrices(refName1aN, ... 
                                                FnameAdata1aN,mFormat);
   disp('  ');
   disp(table1aN);
   disp('  ');
% ------------------------------------------------------------------------
   ns = 125;  %=> split data into 4 groups per system
   trait1aF4 = getMultivariateStats4sploc(AmatrixInfo1aF,ns,0,mType);
   trait1aN4 = getMultivariateStats4sploc(AmatrixInfo1aN,ns,0,mType);
   ns = 500;  %=> do not split the data: only 1 group per system
   trait1aF1 = getMultivariateStats4sploc(AmatrixInfo1aF,ns,0,mType);
   trait1aN1 = getMultivariateStats4sploc(AmatrixInfo1aN,ns,0,mType);
   disp('   ');
% =========================================================== get traits1b
   disp( dividerLine(['constructing traits from: ',fName1b]) );
% ------------------------------------------------------------------------
   [~,~,FnameAdata1bF] = readFileNameList(fName1b,0,'sType','F');
   [~,~,FnameAdata1bN] = readFileNameList(fName1b,0,'sType','N');
   [~,~,FnameAdata1bX] = readFileNameList(fName1b,0);           % X => FUN
% ------------------------------------------------------------------------
   refName1bF = [fName1b,'_F'];           % make this a suitable file name
   [AmatrixInfo1bF,table1bF] = readDataMatrices(refName1bF, ... 
                                                FnameAdata1bF,mFormat);
   disp('  ');
   disp(table1bF);
   disp('   ');
% -----------------------------
   refName1bN = [fName1b,'_N'];           % make this a suitable file name
   [AmatrixInfo1bN,table1bN] = readDataMatrices(refName1bN, ... 
                                                FnameAdata1bN,mFormat);
   disp('  ');
   disp(table1bN);
   disp('   ');
% ------------------------------------------------------------------------
   ns = 5000; %=> split data into 4 groups per system
   trait1bF4 = getMultivariateStats4sploc(AmatrixInfo1bF,ns,0,mType);
   trait1bN4 = getMultivariateStats4sploc(AmatrixInfo1bN,ns,0,mType);
   ns= 20000; %=> do not split the data: only 1 group per system
   trait1bF1 = getMultivariateStats4sploc(AmatrixInfo1bF,ns,0,mType);
   trait1bN1 = getMultivariateStats4sploc(AmatrixInfo1bN,ns,0,mType);
   disp('   ');
% ============================================= get discriminant basis set
% ------------------------------------------------------------------ sploc
   pName1a = [prefix,fName1a];
   testFile = getOutputFileName('training',[pName1a,'_splocResults.dlm']);
      if( isfile(testFile) )
      disp(['SBV file = ',testFile]);
      disp( dividerLine('reading discriminant basis set') );
      splocResults1a = readSPLOCresults(testFile,1);
      else
      disp(['sploc base file name = ',pName1a]);
      disp( dividerLine('generating sploc discriminant basis set') );
      splocResults1a = sploc(0,0,pName1a,trait1aF4,trait1aN4,0);
     %splocResults1a = sploc(0,0,pName1a,trait1aF4,trait1aN4,0,0.7);
      end
   U1a = getDiscriminantSBV(splocResults1a);
% ------------------------------------------------- get another basis: U1b 
   pName1b = [prefix,fName1b];
   testFile = getOutputFileName('training',[pName1b,'_splocResults.dlm']);
      if( isfile(testFile) )
      disp(['SBV file = ',testFile]);
      disp( dividerLine('reading discriminant basis set') );
      splocResults1b = readSPLOCresults(testFile,1);
      else
      disp(['sploc base file name = ',pName1b]);
      disp( dividerLine('generating sploc discriminant basis set') );
      splocResults1b = sploc(0,0,pName1b,trait1bF4,trait1bN4,0);
      end  
   U1b = getDiscriminantSBV(splocResults1b);
% ------------------------------------------------------------ error check
   [~,j1a] = size(U1a);
      if( j1a == 0 )
      disp('U1a does not exist, fiddle with vT');
      end
   [~,j1b] = size(U1b);
      if( j1b == 0 )
      disp('U1b does not exist, fiddle with vT');
      end
      if( (j1a == 0) || (j1b == 0) )
      error('manually change script, by seting vT lower to get a SBV');
      end
% ========================================================= classification
   disp('   ');
   msg = 'Prepare systems to be classified with 500 samples';   % X => FUN
   disp( dividerLine(msg) );   
   [AmatrixInfo1aX,table1aX] = readDataMatrices(fName1a, ... 
                                                FnameAdata1aX,mFormat);
   disp('  ');
   disp(table1aX);
   ns = 500; %=> do not split the data: only 1 group per system
   trait1aX = getMultivariateStats4sploc(AmatrixInfo1aX,ns,0,mType);
   disp('   ');
   msg = 'Prepare systems to be classified with 20000 samples'; 
   disp( dividerLine(msg) );
   [AmatrixInfo1bX,table1bX] = readDataMatrices(fName1b, ... 
                                                FnameAdata1bX,mFormat);
   disp('  ');
   disp(table1bX);
   ns=20000; %=> do not split the data: only 1 group per system
   trait1bX = getMultivariateStats4sploc(AmatrixInfo1bX,ns,0,mType);
   end
%%                                                  test classifyMoments()
   if( test_classifyMoments )
% =============================================================== batch 1a
   disp('   ');
   msg = 'key input: trait1aX, trait1aF4, trait1aN4, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_a1_a4_a4_U1a'];
   [rank_a1_a4_a4_U1a,T] = classifyMoments(cName,trait1aX, ... 
                                   trait1aF4,trait1aN4,U1a,1);
   disp('   ');
   disp(T);
   disp('   ');
   disp( dividerLine('visualize mode contributions') );
      for m=1:rank_a1_a4_a4_U1a.nXsystems
      plotRankDecomposition(rank_a1_a4_a4_U1a,0,m);
      pause
      end
   pause
%% ============================================================== batch 2a
   disp('   ');
   msg = 'key input: trait1bX, trait1aF4, trait1aN4, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_b1_a4_a4_U1a'];
   [~,T] = classifyMoments(cName,trait1bX,trait1aF4, ...
                                   trait1aN4,U1a,1);
   disp('   ');
   disp(T);
   pause
%% ============================================================== batch 3a
   disp('   ');
   msg = 'key input: trait1aX, trait1aF1, trait1aN1, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_a1_a1_a1_U1a'];
   [rank_a1_a1_a1_U1a,T] = classifyMoments(cName,trait1aX, ...
                                   trait1aF1,trait1aN1,U1a,1);
   disp('   ');
   disp(T);
   pause
%% ============================================================== batch 4a
   disp('   ');
   msg = 'key input: trait1bX, trait1aF1, trait1aN1, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'rank_b1_a1_a1_U1a'];
   [rank_b1_a1_a1_U1a,T] = classifyMoments(cName,trait1bX,trait1aF1, ...
                                   trait1aN1,U1a,1);
   disp('   ');
   disp(T);
%    disp('   ');
%    disp( dividerLine('visualize mode contributions') );
%       for m=1:rank_b1_a1_a1_U1a.nXsystems
%       plotRankDecomposition(rank_b1_a1_a1_U1a,0,m);
%       pause
%       end
   pause
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeat the same analysis as above, but with more statistics on the 
% training set.  However, keep the basis vectors from low statistics.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============================================================== batch 5a
   disp('   ');
   msg = 'key input: trait1aX, trait1bF4, trait1bN4, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_a1_b4_b4_U1a'];
   [~,T] = classifyMoments(cName,trait1aX,trait1bF4, ...
                                   trait1bN4,U1a,1);
   disp('   ');
   disp(T);
   pause
%% ============================================================== batch 6a
   disp('   ');
   msg = 'key input: trait1bX, trait1bF4, trait1bN4, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_b1_b4_b4_U1a'];
   [~,T] = classifyMoments(cName,trait1bX,trait1bF4, ...
                                   trait1bN4,U1a,1);
   disp('   ');
   disp(T);
   pause
%% ============================================================== batch 7a
   disp('   ');
   msg = 'key input: trait1aX, trait1bF1, trait1bN1, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_a1_b1_b1_U1a'];
   [~,T] = classifyMoments(cName,trait1aX,trait1bF1, ...
                                   trait1bN1,U1a,1);
   disp('   ');
   disp(T);
   pause
%% ============================================================== batch 8a
   disp('   ');
   msg = 'key input: trait1bX, trait1bF1, trait1bN1, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_b1_b1_b1_U1a'];
   [rank_b1_b1_b1_U1a,T] = classifyMoments(cName,trait1bX, ...
                                   trait1aF4,trait1aN4,U1a,1);
   disp('   ');
   disp(T);
   pause
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare classification quality based on sploc results using maximum 
% statistics and the best estimate for the discriminant basis vectors. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ============================================================== batch 1b
   disp('   ');
   msg = 'key input: trait1bX, trait1bF1, trait1bN1, U1b';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_b1_b1_b1_U1b'];
   [rank_b1_b1_b1_U1b,T] = classifyMoments(cName,trait1bX, ...
                                   trait1bF1,trait1bN1,U1b,1);
   disp('   ');
   disp(T);
   end
%%                                               test classifyDataStream()
   if( test_classifyDataStream )
%% =============================================================== batch 1
   disp('   ');
   msg = 'key input: AmatrixInfo1aX, AmatrixInfo1aF, AmatrixInfo1aN, U1a';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_aaa_U1a'];
   [rank_aaa_U1a,T] = classifyDataStream(cName,AmatrixInfo1aX, ... 
                                   AmatrixInfo1aF,AmatrixInfo1aN,U1a,1);
   disp('   ');
   disp(T);
   disp('   ');
   disp( dividerLine('visualize mode contributions') );
      for m=1:rank_aaa_U1a.nXsystems
      plotRankDecomposition(rank_aaa_U1a,0,m);
      pause
      end
   pause
%% =============================================================== batch 2
   disp('   ');
   msg = 'key input: AmatrixInfo1bX, AmatrixInfo1bF, AmatrixInfo1bN, U1b';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_bbb_U1b'];
   [rank_bbb_U1b,T] = classifyDataStream(cName,AmatrixInfo1bX, ... 
                                   AmatrixInfo1bF,AmatrixInfo1bN,U1b,1);
   disp('   ');
   disp(T);
   disp('   ');
   disp( dividerLine('visualize mode contributions') );
      for m=1:rank_bbb_U1b.nXsystems
      plotRankDecomposition(rank_bbb_U1b,0,m);
      pause
      end
   end
%%                                            test plotRankDecomposition()
   if( test_plotRankDecomposition )
   disp('verbosity: {0,1,2,3} --> for plotting functions:');
   disp('            0,1 => print figures only to the screen.');
   disp('            2   => print figures to file only.');
   disp('            3   => print figures to file and screen.');
   disp('  ');
   v = input('Enter the verbosity level: ');   
% -------------------------------------------- work with rank_b1_b1_b1_U1b
   disp('   ');
   msg = 'key input: trait1bX, trait1bF1, trait1bN1, U1b';
   disp( dividerLine(msg) );
   cName = [fName1,'_rank_b1_b1_b1_U1b'];
   [rank_b1_b1_b1_U1b,T] = classifyMoments(cName,trait1bX, ...
                                   trait1bF1,trait1bN1,U1b,0);
   disp('   ');
   disp(T);
   disp('   ');
   %disp( dividerLine('rank_b1_b1_b1_U1b') );
      for m=1:rank_b1_b1_b1_U1b.nXsystems
      plotRankDecomposition(rank_b1_b1_b1_U1b,v,m);
      pause
      end
   end
%%                                      test classifyMoments() with switch
   if( test_classifyMoments_switch )
%% =============================================================== batch 1
   disp('   ');
   msg = 'key input: trait1bX, trait1bN1, trait1bF1, U1b';
   disp( dividerLine(msg) );
   cName = [fName1,'_sw_rank_b1_b1_b1_U1b'];
   [rank_b1_b1_b1_U1b,T] = classifyMoments(cName,trait1bX, ...
                                   trait1bN1,trait1bF1,U1b,1);
   disp('   ');
   disp(T);
   end
%%                                   test classifyDataStream() with switch
   if( test_classifyDataStream_switch )
%% =============================================================== batch 1
   disp('   ');
   msg = 'key input: AmatrixInfo1bX, AmatrixInfo1bN, AmatrixInfo1bF, U1b';
   disp( dividerLine(msg) );
   cName = [fName1,'_sw_rank_b1_b1_b1_U1b'];
   [rank_sw_bbb_U1b,T] = classifyDataStream(cName,AmatrixInfo1bX, ... 
                                   AmatrixInfo1bN,AmatrixInfo1bF,U1b,1);
   disp('   ');
   disp(T);
   end 