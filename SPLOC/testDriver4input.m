% input test driver for the SPLOC toolset
% ------------------------------------------------------------------------
% SPLOC = Supervised Projective Learning with Orthogonal Completeness
% 
% This test driver checks sploc functions that substantially involve the
% reading data from and/or writing data to the input directory.
% 
% In addition to testing these functions, the purpose of this test driver
% is to serve as an example for how to use these sploc functions. 
%
% Toolset: I/O directory structure:
% 
% current directory --> sub-directories: splocLibrary
%                                        input
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
clear all
clc
%%                                                     build function menu
test_readFileNameList = false;
test_readDataMatrices = false;
test_getMultivariateStats4sploc = false;
test_getBoostedMVStats4sploc = false;
test_readTraits4sploc = false;
test_readSPLOCresults = false;
test_readPropertySpecs = false;
disp('  ');
disp( dividerLine('select function to test') )
disp('1. readFileNameList()')
disp('2. readDataMatrices()')
disp('3. getMultivariateStats4sploc()')
disp('4. getBoostedMVStats4sploc()')
disp('5. readTraits4sploc()')
disp('6. readSPLOCresults()   and  writeSPLOCresults()')
disp('7. readPropertySpecs()')
disp('  ');
nMenu = input('   Enter option: ');
   switch nMenu
       case 1
       test_readFileNameList = true;
       case 2
       test_readDataMatrices = true;
       case 3
       test_getMultivariateStats4sploc = true;
       case 4
       test_getBoostedMVStats4sploc = true;
       case 5
       test_readTraits4sploc = true;
       case 6
       test_readSPLOCresults = true;
       case 7
       test_readPropertySpecs = true;
       otherwise
       error('unknown option');
   end
prefix = ['opt',num2str(nMenu,'%02i')];
logFile = [prefix,mfilename];
initializeSPLOC('fName',logFile,'gType','png');
%%                                                 test readFileNameList()
% test files        format                         Remarks
% ex_1col.txt       fileNames                      no classification
% ex_1colDS.txt     datamatrix                     file name classifies
% ex_2colQM.txt     Qmatrix aveVector              file name classifies
% ex_FUN2colDS.txt  classifer datamatrix           1st column classifies
% ex_FN2colDS.txt   classifer datamatrix           1st column classifies
% ex_FUN3colQM.txt  classifer Qmatrix aveVector    1st column classifies
% ex_4columns       generic list                   no classification
% ------------------------------------------------------------------------
if( test_readFileNameList )
[~,a] = readFileNameList('ex_1col.txt'); 
disp('  ');
disp(a);
disp('  ');
pause
[~,a] = readFileNameList('ex_1col.txt',3);
pause
[~,classifer,a] = readFileNameList('ex_1col.txt','cType','F');
[~,classifer,a] = readFileNameList('input/ex_1col.txt',3,1,'cType','F');
pause
[~,classifer,a] = readFileNameList('input/ex_1colDS.txt',3,1,'cType','N');
pause
[~,classifer,a] = readFileNameList('ex_FUN2colDS.txt',0);
[~,classifer,a] = readFileNameList('input/ex_FUN2colDS.txt',3,1);
pause
[~,classifer,a,b] = readFileNameList('ex_FUN3colQM.txt',0,0);
[~,classifer,a,b] = readFileNameList('input/ex_FUN3colQM.txt',3,1);
pause
[~,classID,a,b] = readFileNameList('input/ex_FUN3colQM.txt',3,1, ...
                                   'cType','U');
% REMARKS: The full or relative path can be specified. An override of the
% classification in the file can be enforced. It is useful to reclassify
% function and nonfunction systems to unknown systems to check classifer. 

pause
[~,a,b,c,d] = readFileNameList('ex_4columns',3,'IDcol',0);
pause
[~,a,b,classID,d] = readFileNameList('input/ex_4columns',3,1, ...
                                     'IDcol',3,'cType','U');
pause
% ========================================================== use selection
disp('  ');
disp( dividerLine('examples using selection: sType') );
[m,classifer,a] = readFileNameList('input/ex_FUN2colDS.txt',3,1, ...
                                   'sType','f');
disp(m)
disp([classifer',a'])
pause
[m,classifer,a] = readFileNameList('ex_FUN2colDS.txt',3,0, ...
                                   'sType','U');
disp(m)
disp([classifer',a'])
pause

[m,classifer,a] = readFileNameList('ex_FN2colDS.txt',3,'sType','F');
disp(m)
disp([classifer',a'])
pause

% --------------------------------------------------- example application1
disp('  ');
disp( dividerLine('example 1: reading in trainingSet1') );
[nLines,classifer,dataMatrixFname] = readFileNameList('trainingSet1',3);
disp(['total number of files = ',num2str(nLines)]);
disp('   ');
% -------------------------------------- extract file names per categories
nF0 = sum( strcmp( classifer,'F' ) );
nN0 = sum( strcmp( classifer,'N' ) );
nU0 = sum( strcmp( classifer,'U' ) );
FnameAdataF = cell(1,nF0);
FnameAdataN = cell(1,nN0);
FnameAdataU = cell(1,nU0);
kF = 0;
kN = 0;
kU = 0;
   for k=1:nLines
      if( strcmp( classifer{k},'F' ) )
      kF = kF + 1;
      FnameAdataF{kF} = dataMatrixFname{k};
      elseif( strcmp( classifer{k},'N' ) )
      kN = kN + 1;
      FnameAdataN{kN} = dataMatrixFname{k};
      else
      kU = kU + 1;
      FnameAdataU{kU} = dataMatrixFname{k};
      end
   end
pause
% --------------------------------------------------- example application2
disp('  ');
disp( dividerLine('example 2: sorting FUN within trainingSet1') );
% ----------------- read in input data file names and their classification
[nF0,~,FnameAdataF] = readFileNameList('trainingSet1','sType','F');
[nN0,~,FnameAdataN] = readFileNameList('trainingSet1','sType','N');
[nU0,~,FnameAdataU] = readFileNameList('trainingSet1','sType','U');
% ---------------------------------------------------------- confirm input
disp('  ');
disp('List of file names for data matrices that are functional');
disp('  ');
   for k=1:nF0
   disp( FnameAdataF{k} );
   end
% ---------------------
disp('  ');
disp('List of file names for data matrices that are nonfunctional');
disp('  ');
   for k=1:nN0
   disp( FnameAdataN{k} );
   end
% ---------------------
disp('  ');
disp('List of file names for data matrices that are unclassified');
disp('  ');
   for k=1:nU0
   disp( FnameAdataU{k} );
   end

disp('  ');
disp( dividerLine('EXPECTED ERROR') );
msg = ['copy from script and paste to command window ', ...
       'to create more errors'];
disp( dividerLine(msg) );
disp('  ');
   fprintf(2,'Example of ERRORS:  You are suppose to get an error!\n');
   fprintf(2,'investigate this script to see why error occurred. \n');
   disp('  ');
disp('  ');
% ----------------------------------------------------- examples of errors
[~,classifer,a] = readFileNameList('ex_FUN2colDS.txt',0,1);
[~,classifer,a] = readFileNameList('input/ex_FUN2colDS.txt');
[m,classifer,a] = readFileNameList('input/ex_FN2colDS.txt',3, ...
                                   'sType','F','cType','u');
[m,classifer,a] = readFileNameList('input/ex_FN2colDS.txt',3, ...
                                   'sType','U');
[~,a,b,c,d] = readFileNameList('input/ex_4columns',3,'IDcol',3);
[~,classID,a,b] = readFileNameList('input/ex_4columns',3,'IDcol',0); 
[~,a,b,c,d] = readFileNameList('input/ex_4columns',3,'IDcol',8);
[~,classifer,a,b,c] = readFileNameList('input/ex_FUN3colQM.txt');
[~,classifer,a] = readFileNameList('input/ex_FUN3colQM.txt'); % too few
[~,a] = readFileNameList('input/ex_1col.txt',3,'cType','F');
end
%%                                                   read data matrix data
   if( test_readDataMatrices )
   [~,~,FnameAdataF] = readFileNameList('trainingSet1','sType','F');
   [~,~,FnameAdataN] = readFileNameList('trainingSet1','sType','N');
   [~,~,FnameAdataU] = readFileNameList('trainingSet1','sType','U');
%  mFormat = setDataMatrixFormat;            % note that mFormat is needed
% ========================================================== start testing
   mFormat = setDataMatrixFormat('xxx-yyy-zzz',2); 
   [AmatrixInfoF,tableF] = readDataMatrices('functional', ...
                                         FnameAdataF,mFormat);
   disp('   ');
   disp(dividerLine('Table summary of functional data matrices read'));
   disp(tableF);
   disp( dividerLine('contents of data structure: AmatrixInfoF') );
   disp('   ');
   disp( AmatrixInfoF);
   pause
% ------------------------------------------------------------------------
   %[AmatrixInfoN,~] = readDataMatrices('nonfunctional',FnameAdataN);
   % -------------^----> The table is for checking: do not need to use it.
   [AmatrixInfoN,tableN] = readDataMatrices('nonfunctional', ...
                                            FnameAdataN,mFormat);
   disp('   ');
   disp(dividerLine('Table summary of nonfunctional data matrices read'));
   disp(tableN);
   disp( dividerLine('contents of data structure: AmatrixInfoN') );
   disp('   ');
   disp(AmatrixInfoN);
   pause
% ------------------------------------------------------------------------
   [AmatrixInfoU,tableU] = readDataMatrices('unknown',FnameAdataU);
   disp('   ');
   disp('Can read in format of data matrices on the fly if needed');
   disp(dividerLine('Table summary of unknown data matrices read'));
   disp(tableU);
   disp( dividerLine('contents of data structure: AmatrixInfoU') );
   disp('   ');
   disp(AmatrixInfoU);
   pause
disp('  ');
disp( dividerLine('EXPECTED ERROR') );
msg = ['copy from script and paste to command window ', ...
       'to create more errors'];
disp( dividerLine(msg) );
disp('  ');
   fprintf(2,'Example of ERRORS:  You are suppose to get an error!\n');
   fprintf(2,'investigate this script to see why error occurred. \n');
   disp('  ');
disp('  ');
% ----------------------------------------------------- examples of errors
   mFormat = setDataMatrixFormat('xyz-xyz-xyz',3);  % error: 3 should be 2
   [AmatrixInfoF,tableF] = readDataMatrices('functional', ...
                                         FnameAdataF,mFormat);                                       
   [AmatrixInfo1,table1] = readDataMatrices('functional', ...
                                            FnameAdataF,mFormat,'colVar');
   end
%%                  test conversion from A matrices to covariance matrices
   if( test_getMultivariateStats4sploc )
   [~,~,FnameAdataF] = readFileNameList('trainingSet1','sType','F');
   [~,~,FnameAdataN] = readFileNameList('trainingSet1','sType','N');
   mFormat = setDataMatrixFormat('xyz-xyz-xyz',2);
   [AmatrixInfoF,~] = readDataMatrices('functional',FnameAdataF,mFormat);
   [AmatrixInfoN,~] = readDataMatrices('nonfunctional', ...
                                       FnameAdataN,mFormat);
   traitCovF = getMultivariateStats4sploc(AmatrixInfoF,125)
   pause
   traitCovN = getMultivariateStats4sploc(AmatrixInfoN,125,0)
   pause
   traitCorF = getMultivariateStats4sploc(AmatrixInfoF,125,1,'cor')
   pause
   disp( dividerLine('verbosity = 3 => show the matrices for checking') );
   disp('  ');
   traitCorN = getMultivariateStats4sploc(AmatrixInfoN,125,3,'cor')
   end
%%                  test conversion from A matrices to covariance matrices
   if( test_getBoostedMVStats4sploc )
   [~,~,FnameAdataF] = readFileNameList('trainingSet1','sType','F');
   [~,~,FnameAdataN] = readFileNameList('trainingSet1','sType','N');
   mFormat = setDataMatrixFormat('xyz-xyz-xyz',2);
   [AmatrixInfoF,~] = readDataMatrices('functional',FnameAdataF,mFormat);
   [AmatrixInfoN,~] = readDataMatrices('nonfunctional', ...
                                       FnameAdataN,mFormat);
   traitCovF = getBoostedMVStats4sploc(AmatrixInfoF,2)
   pause
   traitCovN = getBoostedMVStats4sploc(AmatrixInfoN,2,0)
   pause
   traitCorF = getBoostedMVStats4sploc(AmatrixInfoF,2,1,'cor')
   pause
   disp( dividerLine('verbosity = 3 => show the matrices for checking') );
   disp('  ');
   traitCorN = getBoostedMVStats4sploc(AmatrixInfoN,2,3,'cor')
   end
%%                                                   read traits for sploc
   if( test_readTraits4sploc )
   [n1,~,cFname1,mFname1,str_ns1] = readFileNameList('trainingMVS1', ...
                                                     'sType','F');
   ns1 = zeros(1,n1);
      for k=1:n1
      ns1(k) = str2num(str_ns1{k});
      end
   nSamples = min(ns1);
   b = max(ns1);
      if( nSamples ~= b )
      error('all data does not have the same statistical significance');
      end
   mFormat = setDataMatrixFormat('xxx-yyy-zzz',2); 
   [traitDataF,T] = readTraits4sploc('testF', ...
                                       cFname1,mFname1,nSamples,mFormat);
   disp('  ');
   msg = 'Table summary of functional correlation matrix data';
   disp(dividerLine(msg));
   disp('  ');
   disp(T);
   disp( dividerLine('contents of data structure: traitDataF') );
   disp('   ');
   disp(traitDataF);
   pause
   [n0,~,cFname0,mFname0,str_ns0] = readFileNameList('trainingMVS1', ...
                                                     'sType','N');
   ns0 = zeros(1,n0);
      for k=1:n0
      ns0(k) = str2num(str_ns0{k});
      end
   nSamples = min(ns0);
   b = max(ns0);
      if( nSamples ~= b )
      error('all data does not have the same statistical significance');
      end
   [traitDataN,T] = readTraits4sploc('testN', ...
                                       cFname0,mFname0,nSamples,mFormat);
   disp('  ');
   msg = 'Table summary of nonfunctional correlation matrix data';
   disp(dividerLine(msg));
   disp('  ');
   disp(T);
   disp( dividerLine('contents of data structure: traitDataN') );
   disp('   ');
   disp(traitDataN);
   end
%%                                                   test readSPLOCresults
   if( test_readSPLOCresults )
   [~,~,FnameAdataF] = readFileNameList('trainingSet1','sType','F');
   [~,~,FnameAdataN] = readFileNameList('trainingSet1','sType','N');
   mFormat = setDataMatrixFormat('xxx-yyy-zzz',2);
   [AmatrixInfoF,~] = readDataMatrices('testSet1_F',FnameAdataF,mFormat);
   [AmatrixInfoN,~] = readDataMatrices('testSet1_N',FnameAdataN,mFormat);
   traitCovF = getMultivariateStats4sploc(AmatrixInfoF,125,0);
   traitCovN = getMultivariateStats4sploc(AmatrixInfoN,125,0);
   disp('   ');
   disp('   ');
   fprintf(2,'Please wait: Running sploc() on training data!\n');
   disp('............  ');
   splocResults1 = sploc(0,0,'testSet1',traitCovF,traitCovN,0);
%                  ^^^^^---> uses writeSPLOCresults()
%%
% ------------------------------------------------------------------------
   disp('  ');
   disp( dividerLine('require numerical error < 1.0E-09') );
   disp('  ');
   fname = 'testSet1'; 
   splocResults2 = readSPLOCresults(fname);             % use default path
  %splocResults2 = readSPLOCresults(fname,0);           % use default path
  %splocResults2 = readSPLOCresults(fname,1);   % fname includes full path
% ---------------------------------------------------- do some comparisons
   sbv1 = splocResults1.SBV;
   sbv2 = splocResults2.SBV;
   temp = max( max( abs(sbv2 - sbv1) ) );
   disp(['comparison of SBV1 and SBV2: max( |diff| ) = ',num2str(temp)]);
  %
   SEVd1 = splocResults1.SEVd;
   SEVd2 = splocResults2.SEVd;
   temp = max( abs(SEVd2 - SEVd1) );
   disp(['comparison of SEVd2 & SEVd1: max( |diff| ) = ',num2str(temp)]);
  %
   SEVi1 = splocResults1.SEVi;
   SEVi2 = splocResults2.SEVi;
   temp = max( abs(SEVi2 - SEVi1) );
   disp(['comparison of SEVi2 & SEVi1: max( |diff| ) = ',num2str(temp)]);
  %
   CEVd1 = splocResults1.CEVd;
   CEVd2 = splocResults2.CEVd;
   temp = max( abs(CEVd2 - CEVd1) );
   disp(['comparison of CEVd2 & CEVd1: max( |diff| ) = ',num2str(temp)]);
  %
   CEVi1 = splocResults1.CEVi;
   CEVi2 = splocResults2.CEVi;
   temp = max( abs(CEVi2 - CEVi1) );
   disp(['comparison of CEVi2 & CEVi1: max( |diff| ) = ',num2str(temp)]);
  %
   QEVd1 = splocResults1.QEVd;
   QEVd2 = splocResults2.QEVd;
   temp = max( abs(QEVd2 - QEVd1) );
   disp(['comparison of QEVd2 & QEVd1: max( |diff| ) = ',num2str(temp)]);
  %
   QEVi1 = splocResults1.QEVi;
   QEVi2 = splocResults2.QEVi;
   temp = max( abs(QEVi2 - QEVi1) );
   disp(['comparison of QEVi2 & QEVi1: max( |diff| ) = ',num2str(temp)]);
  %
%   EEVd1 = splocResults1.EEVd;
%   EEVd2 = splocResults2.EEVd;
%   temp = max( abs(EEVd2 - EEVd1) );
%   disp(['comparison of EEVd2 & EEVd1: max( |diff| ) = ',num2str(temp)]);
  %
   EEVi1 = splocResults1.EEVi;
   EEVi2 = splocResults2.EEVi;
   temp = max( abs(EEVi2 - EEVi1) );
   disp(['comparison of EEVi2 & EEVi1: max( |diff| ) = ',num2str(temp)]);
  %
   Cind1 = splocResults1.Cind;
   Cind2 = splocResults2.Cind;
   temp = max( abs(Cind2 - Cind1) );
   disp(['comparison of Cind2 & Cind1: max( |diff| ) = ',num2str(temp)]);
  % 
   vT1 = splocResults1.vT;
   vT2 = splocResults2.vT;
   temp = max( abs(vT2 - vT1) );
   disp(['comparison of vT2  and  vT1: max( |diff| ) = ',num2str(temp)]);
   % 
%   efficacy1 = splocResults1.efficacy;
%   efficacy2 = splocResults2.efficacy;
%   temp = max( abs(efficacy2 - efficacy1) );
%   disp(['comp: efficacy2 & efficacy1: max( |diff| ) = ',num2str(temp)]);
  % 
   Dd1 = splocResults1.Dd;
   Dd2 = splocResults2.Dd;
   temp = max( abs(Dd2 - Dd1) );
   disp(['comparison of Dd2  and  Dd1: max( |diff| ) = ',num2str(temp)]);
   disp('   ');
   pause
% ------------------------------------------------------------------------
   fname = 'training/testSet1_splocResults.dlm';
   splocResults3 = readSPLOCresults(fname,1);   % fname includes full path
% ---------------------------------------------------------- do comparison
   sbv3 = splocResults3.SBV;
   temp = max( max( abs(sbv3 - sbv1) ) );
   disp(['comparison of SBV3 and SBV1: max( |diff| ) = ',num2str(temp)]);
  %
   Cind3 = splocResults3.Cind;
   temp = max( abs(Cind3 - Cind1) );
   disp(['comparison of Cind3 & Cind1: max( |diff| ) = ',num2str(temp)]);
  % 
   Dd3 = splocResults3.Dd;
   temp = max( abs(Dd3 - Dd1) );
   disp(['comparison of Dd3  and  Dd1: max( |diff| ) = ',num2str(temp)]);
   disp('   ');
   pause
   disp( dividerLine('works if you forget the .dlm extension') );
   fname = 'training/testSet1';
   disp( fname );
   disp('  ');
   splocResults4 = readSPLOCresults(fname,1);   % should include full path
% ---------------------------------------------------------- do comparison
   sbv4 = splocResults4.SBV;
   temp = max( max( abs(sbv4 - sbv1) ) );
   disp(['comparison of SBV4 and SBV1: max( |diff| ) = ',num2str(temp)]);
  %
   Cind4 = splocResults4.Cind;
   temp = max( abs(Cind4 - Cind1) );
   disp(['comparison of Cind4 & Cind1: max( |diff| ) = ',num2str(temp)]);
  % 
   Dd4 = splocResults4.Dd;
   temp = max( abs(Dd4 - Dd1) );
   disp(['comparison of Dd4  and  Dd1: max( |diff| ) = ',num2str(temp)]);
   disp('   ');
   pause
   disp('  ');
   disp( dividerLine('EXPECTED ERROR') );
   msg = ['copy from script and paste to command window ', ...
          'to create more errors'];
   disp( dividerLine(msg) );
   disp('  ');
   fprintf(2,'Example of ERRORS:  You are suppose to get an error!\n');
   fprintf(2,'investigate this script to see why error occurred. \n');
   disp('  ');
% ----------------------------------------------------- examples of errors
   fname = 'testSet1';
   splocResults5 = readSPLOCresults(fname,1);   % should include full path
   end
%%                                                test readPropertySpecs()    
   if( test_readPropertySpecs )
   [nFiles1,pID1,fileList1] = readPropertySpecs('propertySet1',3);
   pause
   disp('  ');
   disp( dividerLine('propertySet1B') );
   [nFiles1B,pID1B,fileList1B] = readPropertySpecs('propertySet1B',0);
   T = table(pID1B',fileList1B');
   disp('   ');
   disp(T);
   fullPath1 = 'input/propertySet1'; 
   [nFiles1f,pID1f,fileList1f] = readPropertySpecs(fullPath1,3,1);
%  REMARK: The full or relative path can be specified.
   pause
   disp('  ');
   disp( dividerLine('EXPECTED ERROR') );
   msg = ['copy from script and paste to command window ', ...
          'to create more errors'];
   disp( dividerLine(msg) );
   disp('  ');
   fprintf(2,'Example of ERRORS:  You are suppose to get an error!\n');
   fprintf(2,'investigate this script to see why error occurred. \n');
   disp('  ');
   disp('  ');
% ----------------------------------------------------- examples of errors
   [nFiles1f,pID1f,fileList1f] = readPropertySpecs('fullPath1',3,1);
   end

   
