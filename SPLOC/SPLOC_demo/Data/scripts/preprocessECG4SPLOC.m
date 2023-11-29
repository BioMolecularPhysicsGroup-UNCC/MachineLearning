% This script is used to preprocess ECG data. It first identifies the QRS
% complex by marking the R-peaks and it also removes drift. There are 
% many flags that can separate the data into training and testing based 
% on user defined percentage and allows to user to select which type of 
% files are functional (healthy/control) and non-functional (various 
% diseases). Note that the flags can be generalized, but are designed for
% a particular dataset obtained from PTB PhysioNet-Cardiovascular data.
% For all channels 1 through 15, except, 7,8,9 the R-peak is considered 
% to be positive, and these channels are called V1, V2, V3 respectively.
% The program will allow the user to select the channels that should 
% process the R-peak as a negative or positive value. In addition, the 
% user can select the number of channels to consider. 
%
% INPUT 
% A directory structure (folders) that is identical to the PTB PhysioNet-
% cardiovascular data with the files placed appropriately within the 
% various folders. Beyond this raw data, all other input is user-defined
% parameters. 
%
% percentTest = allows a user to randomly select the data and partition it
% into training and testing data. Testing is referred to as checking here.
% 
% seed = way to control different randomizations
% channelList = is a way to select which channels to consider
% negChannelList = is a way to select channels that use negative R peaks
%
%
% checkReadDataXY => visually check X=F/N and Y=train/test (four cases)
%
% flag_fileNames => select files within certain folders to process or not
%
% OUTPUT: Main cell arrays. There is more: See OUTPUT ARRAYS
% level_Ftrain{1:nFtrain} = cell array for nFtrain leveled data matrices
% level_Fcheck{1:nFcheck} = cell array for nFcheck leveled data matrices
% level_Ntrain{1:nNtrain} = cell array for nNtrain leveled data matrices
% level_Ncheck{1:nNcheck} = cell array for nNcheck leveled data matrices
% ------------------------------------------------------------------------
%clear all                                 % never not do this (sometimes)
%%                                           user input control parameters
percentTest = 60;
checkReadDataFtrain = false;   % (true,false) 
checkReadDataFcheck = false;   % (true,false)    % Note check => test data
checkReadDataNtrain = false;   % (true,false)
checkReadDataNcheck = false;   % (true,false)    % Note check => test data
seed = 1209;
nCmax = 15;                                   % maximum number of channels
% ------------------------------------------------------ channel selection
channelList = 1:15;  %[6,7,8,9,10];
% ----------------------------------------- negative QRS complex selection
negChannels = [7,8,9];
% --------------------------------- get training/testing datasets & labels
%                                          REMARK: not training => testing
%                                          REMARK:     testing => checking
% ---------------------------------------- choose functional training data
% -1 => nonfunctional data to be used
%  0 => do not use this data
%  1 => functional data to be used
% ------------------------------------------------------- folder selection
flag_HealthyControls = 1;             % likely this will never be modified
% -----------------------------
flag_BundelBranchBlock    = 0; 
flag_Dysrhythmia          = 0;
flag_Hypertrophy          = -1;
flag_Myocarditis          = 0;
flag_Cardiomyopathy       = 0;
flag_MyocardialInfarction = 0;
flag_ValvularHeartDisease = 0;
% ========================================================================
folderNameList = {'Bundle Branch Block', ... 
                  'Dysrhythmia', ...
                  'Hypertrophy', ...
                  'Myocarditis', ...
                  'Cardiomyopathy', ...
                  'Healthy Controls', ...
                  'Myocardial Infarction', ...
                  'Valvular Heart Disease'};
folderFlagList = [flag_BundelBranchBlock, ...
                  flag_Dysrhythmia, ...
                  flag_Hypertrophy, ...
                  flag_Myocarditis, ...
                  flag_Cardiomyopathy, ...
                  flag_HealthyControls, ...
                  flag_MyocardialInfarction, ...
                  flag_ValvularHeartDisease];
nFlaggedFolders = length(folderFlagList);
% ----------------------------------------------------- silly sanity check
  if( nFlaggedFolders ~= length(folderNameList) )
  error(['Number of file names in folderNameList is ', ... 
         'not equal to nFlaggedFolders']);
  end
testSum = sum( (channelList > nCmax) );
  if( testSum > 0 )
  error('channelList must not have any channel exceeding 15');
  end
testSum = sum( (channelList < 1) );
  if( testSum > 0 )
  error('channelList must not have any channel under 1');
  end
testSum = sum( (negChannels > nCmax) );
  if( testSum > 0 )
  error('negChannels must not have any channel exceeding 15');
  end
testSum = sum( (negChannels < 1) );
  if( testSum > 0 )
  error('negChannels must not have any channel under 1');
  end
% ----------------------------------------------------------- define signs
peakSign = ones(1,nCmax);                             % pSign => peak sign
peakSign(negChannels) = -1;
%%                                                           read the data
rng(seed);                       % for separation between training/testing
dirList = dir('**/*.mat');
fileList = struct2cell(dirList);
% --------------------------------------------- fetch the selected folders
[~,nFolders0] = size(fileList);
countTraceTypes = zeros(1,nFlaggedFolders);
nF = 0;
nN = 0;
  for j=1:nFolders0
    for k=1:nFlaggedFolders
      if( folderFlagList(k) ~= 0 )
      checkTrue = contains( fileList{2,j} ,folderNameList{k} );
      else
      checkTrue = false; % data is skipped
      end
      if( checkTrue )          % count only the training and testing files
      countTraceTypes(k) = countTraceTypes(k) + 1;
        if( folderFlagList(k) > 0 )
        nF = nF + 1;
        else
        nN = nN + 1;
        end
      end
    end
  end
% ---------------------------------------------------- simple saniety test
  if( (percentTest > 80) || (percentTest < 20) )
  error('percentTest must be 20 to 80 percent!')
  end
% ------------------ split data into functional/nonfunctional & train/test
nFtrain = 0;
nNtrain = 0;
nFcheck = 0;
nNcheck = 0;
flagTrain = cell(1,nFlaggedFolders);
  for k=1:nFlaggedFolders
  n = countTraceTypes(k);
  flagTrain{k} = ones(1,n);
  m = ceil( percentTest*countTraceTypes(k)/100 );
  flagTrain{k}(1:m) = -1;
  index = randperm(n);
  flagTrain{k} = flagTrain{k}(index);
    if( folderFlagList(k) > 0 )
    nFtrain = nFtrain + (n - m);
    nFcheck = nFcheck + m;
    else
    nNtrain = nNtrain + (n - m);
    nNcheck = nNcheck + m;
    end
  end
% ----------------------------------------------------- silly sanity check
  if( nF ~= (nFtrain + nFcheck) )
  error('sum rule error: nF is not equal to nFtrain + nFcheck');
  end
  if( nN ~= (nNtrain + nNcheck) )
  error('sum rule error: nN is not equal to nNtrain + nNcheck');
  end
% ------------------------------------------------- collect names of files
inputFileNameFtrain =  cell(1,nFtrain);
inputFileTypeFtrain = zeros(1,nFtrain);
subjectLabelFtrain  =  cell(1,nFtrain);
% -------------------------------------
inputFileNameFcheck =  cell(1,nFcheck);
inputFileTypeFcheck = zeros(1,nFcheck);
subjectLabelFcheck  =  cell(1,nFcheck);
% -------------------------------------
inputFileNameNtrain =  cell(1,nNtrain);
inputFileTypeNtrain = zeros(1,nNtrain);
subjectLabelNtrain  =  cell(1,nNtrain);
% -------------------------------------
inputFileNameNcheck =  cell(1,nNcheck);
inputFileTypeNcheck = zeros(1,nNcheck);
subjectLabelNcheck  =  cell(1,nNcheck);
% ------------------------------------------- fetch selected folders again
countTraceTypes = 0*countTraceTypes;
nFtrain = 0;
nNtrain = 0;
nFcheck = 0;
nNcheck = 0;
  for j=1:nFolders0
    for k=1:nFlaggedFolders
      if( folderFlagList(k) ~= 0 )
      checkTrue = contains( fileList{2,j} ,folderNameList{k} );
      else
      checkTrue = false;                                 % data is skipped
      end
      if( checkTrue )          % count only the training and testing files
        if( folderFlagList(k) > 0 )                        % => functional
        countTraceTypes(k) = countTraceTypes(k) + 1;
        n = countTraceTypes(k);          % F & N have mutually exclusive k
          if( flagTrain{k}(n) > 0 )                          % => training
          nFtrain = nFtrain + 1;
            if( isunix )
            inputFileNameFtrain{nFtrain} = ...
                         [fileList{2,j},'/',fileList{1,j}];
            else
            inputFileNameFtrain{nFtrain} = ...
                         [fileList{2,j},'\',fileList{1,j}];
            end
          inputFileTypeFtrain(1,nFtrain) = k;
          subjectLabelFtrain{1,nFtrain} = fileList{1,j}(1:5);
          else                                               % => checking
          nFcheck = nFcheck + 1;
            if( isunix )
            inputFileNameFcheck{nFcheck} = ...
                         [fileList{2,j},'/',fileList{1,j}];
            else
            inputFileNameFcheck{nFcheck} = ...
                         [fileList{2,j},'\',fileList{1,j}];
            end
          inputFileTypeFcheck(1,nFcheck) = k;
          subjectLabelFcheck{1,nFcheck} = fileList{1,j}(1:5);
          end
        else                                            % => nonfunctional
        countTraceTypes(k) = countTraceTypes(k) + 1;
        n = countTraceTypes(k);          % F & N have mutually exclusive k
          if( flagTrain{k}(n) > 0 )                          % => training
          nNtrain = nNtrain + 1;
            if( isunix )
            inputFileNameNtrain{nNtrain} = ...
                         [fileList{2,j},'/',fileList{1,j}];
            else
            inputFileNameNtrain{nNtrain} = ...
                         [fileList{2,j},'\',fileList{1,j}];
            end
          inputFileTypeNtrain(1,nNtrain) = k;
          subjectLabelNtrain{1,nNtrain} = fileList{1,j}(1:5);
          else                                               % => checking
          nNcheck = nNcheck + 1;
            if( isunix )
            inputFileNameNcheck{nNcheck} = ...
                         [fileList{2,j},'/',fileList{1,j}];
            else
            inputFileNameNcheck{nNcheck} = ...
                         [fileList{2,j},'\',fileList{1,j}];
            end
          inputFileTypeNcheck(1,nNcheck) = k;
          subjectLabelNcheck{1,nNcheck} = fileList{1,j}(1:5);
          end
        end
      end
    end
  end
%%                                   read raw data matrices (all channels)
allRawFtrain = cell(1,nFtrain);                 % functional training data
allRawFcheck = cell(1,nFcheck);                 % functional testing  data
allRawNtrain = cell(1,nFtrain);              % nonfunctional training data
allRawNcheck = cell(1,nFcheck);              % nonfunctional testing  data
  for j=1:nFtrain
  allRawFtrain{j} = load(inputFileNameFtrain{j});
  end
  for j=1:nNtrain
  allRawNtrain{j} = load(inputFileNameNtrain{j});
  end
  for j=1:nFcheck
  allRawFcheck{j} = load(inputFileNameFcheck{j});
  end
  for j=1:nNcheck
  allRawNcheck{j} = load(inputFileNameNcheck{j});
  end
%%                                                parse data into channels
rawFtrain = cell(1,nFtrain);                    % functional training data
rawFcheck = cell(1,nFcheck);                    % functional testing  data
rawNtrain = cell(1,nNtrain);                 % nonfunctional training data
rawNcheck = cell(1,nNcheck);                 % nonfunctional testing  data
nC = length(channelList);
imin = 1000000000;
imax = -imin;
jmin = imin;
jmax = -jmin;
  for j=1:nFtrain
  [ii,jj] = size(allRawFtrain{j}.val);
  imin = min(imin,ii);
  imax = max(imax,ii);
  jmin = min(jmin,jj);
  jmax = max(jmax,jj);
  rawFtrain{j} = allRawFtrain{j}.val(channelList,:);
  end
  for j=1:nNtrain
  [ii,jj] = size(allRawNtrain{j}.val);
  imin = min(imin,ii);
  imax = max(imax,ii);
  jmin = min(jmin,jj);
  jmax = max(jmax,jj);
  rawNtrain{j} = allRawNtrain{j}.val(channelList,:);
  end
  for j=1:nFcheck
  [ii,jj] = size(allRawFcheck{j}.val);
  imin = min(imin,ii);
  imax = max(imax,ii);
  jmin = min(jmin,jj);
  jmax = max(jmax,jj);
  rawFcheck{j} = allRawFcheck{j}.val(channelList,:);
  end
  for j=1:nNcheck
  [ii,jj] = size(allRawNcheck{j}.val);
  imin = min(imin,ii);
  imax = max(imax,ii);
  jmin = min(jmin,jj);
  jmax = max(jmax,jj);
  rawNcheck{j} = allRawNcheck{j}.val(channelList,:);
  end
% ---------------------------------------------------------- write summary
disp(['minimum number of channels = ',num2str(imin)]);
disp(['maximum number of channels = ',num2str(imax)]);
disp(['minimum  # of observations = ',num2str(jmin)]);
disp(['maximum  # of observations = ',num2str(jmax)]);
  if( imin ~= imax )
  error('stop: all subjects do not have the same number of channels')
  end
% different jmin is allowed. We will deal with this in subsampling
%%                                        apply level tranform to the data
% REMARK: After the R-peaks are identified, the averaging for each 
%         window is apply to obtain yLevel. The output is therefore, 
%         transformed y => yLevel, nPeaks, tPeaks and hPeaks.  
% OUTPUT arrays
% ----------------------------------------------- functional training data
level_Ftrain = cell(nC,nFtrain);   
nPeaksFtrain =zeros(nC,nFtrain);
tPeaksFtrain = cell(nC,nFtrain);
hPeaksFtrain = cell(nC,nFtrain);
% ------------------------------------------------ functional testing data
level_Fcheck = cell(nC,nFcheck);
nPeaksFcheck =zeros(nC,nFcheck);
tPeaksFcheck = cell(nC,nFcheck);
hPeaksFcheck = cell(nC,nFcheck);
% -------------------------------------------- nonfunctional training data
level_Ntrain = cell(nC,nNtrain);  
nPeaksNtrain =zeros(nC,nNtrain);
tPeaksNtrain = cell(nC,nNtrain);
hPeaksNtrain = cell(nC,nNtrain);
% --------------------------------------------- nonfunctional testing data
level_Ncheck = cell(nC,nNcheck);
nPeaksNcheck =zeros(nC,nNcheck);
tPeaksNcheck = cell(nC,nNcheck); 
hPeaksNcheck = cell(nC,nNcheck);
% -------------------------------------------------
  for j=1:nFtrain
    for ic=1:nC
    c = channelList(ic);             % c is the channel that user selected
    ps = peakSign(c);
    y = rawFtrain{j}(ic,:);
    [yLevel,nPeaks,tPeaks,hPeaks] = levelData(y,ps,checkReadDataFtrain);
    level_Ftrain{ic,j} = yLevel;
    nPeaksFtrain(ic,j) = nPeaks;
    tPeaksFtrain{ic,j} = tPeaks;
    hPeaksFtrain{ic,j} = hPeaks;
      if( checkReadDataFtrain )
      figure(1);
      clf;
      hold on;
      t = 1:length(y);
      plot(t,y,'b');
      plot(tPeaks,y(tPeaks),'or');
      xlabel('sample points');
      ylabel('measured signal');
      title([subjectLabelFtrain{j},':  channel= ',num2str(c), ... 
            '  sign= ',num2str(ps),' nFtrain= ',num2str(j)]);
      figure(2);
      clf;
      hold on;
      t = tPeaks(1):tPeaks(end);
      plot(t,yLevel,'b');
      plot(tPeaks,hPeaks,'or');
      xlabel('sample points');
      ylabel('leveled signal');
      aveT = round( 10*mean( diff(tPeaks) ) )/10;
      stdT = round( 10*std( diff(tPeaks) ) )/10;
      aveH = round( 10*mean(hPeaks) )/10;
      stdH = round( 10*std(hPeaks) )/10;
      title(['nPeaks= ',num2str(nPeaks), ...
             '     aveT= ',num2str(aveT), ...
             ' stdT= ',num2str(stdT), ...
             '     aveH= ',num2str(aveH), ...
             ' stdH= ',num2str(stdH)]);
      pause(0.8);
      end
    end
  end
% -------------------------------------------------
  for j=1:nFcheck
    for ic=1:nC
    c = channelList(ic);             % c is the channel that user selected
    ps = peakSign(c);
    y = rawFcheck{j}(ic,:);
    [yLevel,nPeaks,tPeaks,hPeaks] = levelData(y,ps,checkReadDataFcheck);
    level_Fcheck{ic,j} = yLevel;
    nPeaksFcheck(ic,j) = nPeaks;
    tPeaksFcheck{ic,j} = tPeaks;
    hPeaksFcheck{ic,j} = hPeaks;
      if( checkReadDataFcheck )
      figure(1);
      clf;
      hold on;
      t = 1:length(y);
      plot(t,y,'b');
      plot(tPeaks,y(tPeaks),'or');
      xlabel('sample points');
      ylabel('measured signal');
      title([subjectLabelFcheck{j},':  channel= ',num2str(c), ... 
            '  sign= ',num2str(ps),' nFcheck= ',num2str(j)]);
      figure(2);
      clf;
      hold on;
      t = tPeaks(1):tPeaks(end);
      plot(t,yLevel,'b');
      plot(tPeaks,hPeaks,'or');
      xlabel('sample points');
      ylabel('leveled signal');
      aveT = round( 10*mean( diff(tPeaks) ) )/10;
      stdT = round( 10*std( diff(tPeaks) ) )/10;
      aveH = round( 10*mean(hPeaks) )/10;
      stdH = round( 10*std(hPeaks) )/10;
      title(['nPeaks= ',num2str(nPeaks), ...
             '     aveT= ',num2str(aveT), ...
             ' stdT= ',num2str(stdT), ...
             '     aveH= ',num2str(aveH), ...
             ' stdH= ',num2str(stdH)]);
      pause(0.8);
      end
    end
  end
% -------------------------------------------------
  for j=1:nNtrain
    for ic=1:nC
    c = channelList(ic);             % c is the channel that user selected
    ps = peakSign(c);
    y = rawNtrain{j}(ic,:);
    [yLevel,nPeaks,tPeaks,hPeaks] = levelData(y,ps,checkReadDataNtrain);
    level_Ntrain{ic,j} = yLevel;
    nPeaksNtrain(ic,j) = nPeaks;
    tPeaksNtrain{ic,j} = tPeaks;
    hPeaksNtrain{ic,j} = hPeaks;
      if( checkReadDataNtrain )
      figure(1);
      clf;
      hold on;
      t = 1:length(y);
      plot(t,y,'b');
      plot(tPeaks,y(tPeaks),'or');
      xlabel('sample points');
      ylabel('measured signal');
      title([subjectLabelNtrain{j},':  channel= ',num2str(c), ... 
            '  sign= ',num2str(ps),' nNtrain= ',num2str(j)]);
      figure(2);
      clf;
      hold on;
      t = tPeaks(1):tPeaks(end);
      plot(t,yLevel,'b');
      plot(tPeaks,hPeaks,'or');
      xlabel('sample points');
      ylabel('leveled signal');
      aveT = round( 10*mean( diff(tPeaks) ) )/10;
      stdT = round( 10*std( diff(tPeaks) ) )/10;
      aveH = round( 10*mean(hPeaks) )/10;
      stdH = round( 10*std(hPeaks) )/10;
      title(['nPeaks= ',num2str(nPeaks), ...
             '     aveT= ',num2str(aveT), ...
             ' stdT= ',num2str(stdT), ...
             '     aveH= ',num2str(aveH), ...
             ' stdH= ',num2str(stdH)]);
      pause(0.8);
      end
    end
  end
% -------------------------------------------------
  for j=1:nNcheck
    for ic=1:nC
    c = channelList(ic);             % c is the channel that user selected
    ps = peakSign(c);
    y = rawNcheck{j}(ic,:);
    [yLevel,nPeaks,tPeaks,hPeaks] = levelData(y,ps,checkReadDataNcheck);
    level_Ncheck{ic,j} = yLevel;
    nPeaksNcheck(ic,j) = nPeaks;
    tPeaksNcheck{ic,j} = tPeaks;
    hPeaksNcheck{ic,j} = hPeaks;
      if( checkReadDataNcheck )
      figure(1);
      clf;
      hold on;
      t = 1:length(y);
      plot(t,y,'b');
      plot(tPeaks,y(tPeaks),'or');
      xlabel('sample points');
      ylabel('measured signal');
      title([subjectLabelNcheck{j},':  channel= ',num2str(c), ... 
            '  sign= ',num2str(ps),'  nNtrain= ',num2str(j)]);
      figure(2);
      clf;
      hold on;
      t = tPeaks(1):tPeaks(end);
      plot(t,yLevel,'b');
      plot(tPeaks,hPeaks,'or');
      xlabel('sample points');
      ylabel('leveled signal');
      aveT = round( 10*mean( diff(tPeaks) ) )/10;
      stdT = round( 10*std( diff(tPeaks) ) )/10;
      aveH = round( 10*mean(hPeaks) )/10;
      stdH = round( 10*std(hPeaks) )/10;
      title(['nPeaks= ',num2str(nPeaks), ...
             '     aveT= ',num2str(aveT), ...
             ' stdT= ',num2str(stdT), ...
             '     aveH= ',num2str(aveH), ...
             ' stdH= ',num2str(stdH)]);
      pause(0.8);
      end
    end
  end

% % OUTPUT arrays to this point (handshake requires some more processing)
% % % ------------------------------------------- functional training data
% % level_Ftrain = cell(nC,nFtrain);   
% % nPeaksFtrain =zeros(nC,nFtrain);
% % tPeaksFtrain = cell(nC,nFtrain);
% % hPeaksFtrain = cell(nC,nFtrain);
% % % -------------------------------------------- functional testing data
% % level_Fcheck = cell(nC,nFcheck);
% % nPeaksFcheck =zeros(nC,nFcheck);
% % tPeaksFcheck = cell(nC,nFcheck);
% % hPeaksFcheck = cell(nC,nFcheck);
% % % ---------------------------------------- nonfunctional training data
% % level_Ntrain = cell(nC,nNtrain);  
% % nPeaksNtrain =zeros(nC,nNtrain);
% % tPeaksNtrain = cell(nC,nNtrain);
% % hPeaksNtrain = cell(nC,nNtrain);
% % % ----------------------------------------- nonfunctional testing data
% % level_Ncheck = cell(nC,nNcheck);
% % nPeaksNcheck =zeros(nC,nNcheck);
% % tPeaksNcheck = cell(nC,nNcheck); 
% % hPeaksNcheck = cell(nC,nNcheck);
%%                                                convert to standard form
% Synchronize the time across all channels
% Start and stop times will be the same across all channels
% Throw out any subject-data with little data (supposed to be the same)
% Each subject is associated with a standard data matrix: p x n
% We end with a cell array of subjects: each subject is a data matrix
% ------------------------------------------------------- get minimum time
max_tMin = -1;
max_tMax = -1;
ave_tMax = 0;
min_tMax = 1000000000;
k = 0;
  for s=1:nFtrain
    for ic=1:nC
    max_tMin = max(max_tMin,tPeaksFtrain{ic,s}(1)); 
    k = k + 1;
    tMax = max(tPeaksFtrain{ic,s}(end));
    ave_tMax = ave_tMax + tMax;
    max_tMax = max(tMax,max_tMax);
    min_tMax = min(tMax,min_tMax);
    end
  end
% ------------------
  for s=1:nFcheck
    for ic=1:nC
    max_tMin = max(max_tMin,tPeaksFcheck{ic,s}(1)); 
    k = k + 1;
    tMax = max(tPeaksFcheck{ic,s}(end));
    ave_tMax = ave_tMax + tMax;
    max_tMax = max(tMax,max_tMax);
    min_tMax = min(tMax,min_tMax);
    end
  end
% ------------------
  for s=1:nNtrain
    for ic=1:nC
    max_tMin = max(max_tMin,tPeaksNtrain{ic,s}(1)); 
    k = k + 1;
    tMax = max(tPeaksNtrain{ic,s}(end));
    ave_tMax = ave_tMax + tMax;
    max_tMax = max(tMax,max_tMax);
    min_tMax = min(tMax,min_tMax);
    end
  end
% ------------------
  for s=1:nNcheck
    for ic=1:nC
    max_tMin = max(max_tMin,tPeaksNcheck{ic,s}(1)); 
    k = k + 1;
    tMax = max(tPeaksNcheck{ic,s}(end));
    ave_tMax = ave_tMax + tMax;
    max_tMax = max(tMax,max_tMax);
    min_tMax = min(tMax,min_tMax);
    end
  end
ave_tMax = floor(ave_tMax/k);
% disp(tMin)
% disp(min_tMax)
% disp(ave_tMax)
% disp(max_tMax)
t1 = max_tMin;
t2 = ave_tMax - (max_tMax - ave_tMax);
disp(['time range: t1= ',num2str(t1),'  t2= ',num2str(t2)]);
% -------------------------------------- throw out data that is incomplete
keepFtrain = -ones(1,nFtrain);
keepNtrain = -ones(1,nNtrain);
keepFcheck = -ones(1,nFcheck);
keepNcheck = -ones(1,nNcheck);
% ----------------------------------------- flag the data that can be kept
  for s=1:nFtrain
  k = 0;
    for ic=1:nC
      if( tPeaksFtrain{ic,s}(end) >= t2 )
      k = k + 1;
      end
    end
    if( k == nC )   % this is a confirmation
    keepFtrain(s) = 1;
    end
  end
% -------------------------------------------
  for s=1:nFcheck
  k = 0;
    for ic=1:nC
      if( tPeaksFcheck{ic,s}(end) >= t2 )
      k = k + 1;
      end
    end
    if( k == nC )   % this is a confirmation
    keepFcheck(s) = 1;
    end
  end
% -------------------------------------------
  for s=1:nNtrain
  k = 0;
    for ic=1:nC
      if( tPeaksNtrain{ic,s}(end) >= t2 )
      k = k + 1;
      end
    end
    if( k == nC )   % this is a confirmation
    keepNtrain(s) = 1;
    end
  end
% -------------------------------------------
  for s=1:nNcheck
  k = 0;
    for ic=1:nC
      if( tPeaksNcheck{ic,s}(end) >= t2 )
      k = k + 1;
      end
    end
    if( k == nC )   % this is a confirmation
    keepNcheck(s) = 1;
    end
  end
kFtrain = sum( keepFtrain > 0 );
kFcheck = sum( keepFcheck > 0 );
kNtrain = sum( keepNtrain > 0 );
kNcheck = sum( keepNcheck > 0 );
disp(['before: nFtrain= ',num2str(nFtrain), ...
      '  nFcheck= ',num2str(nFcheck), ...
      '  nNtrain= ',num2str(nNtrain), ...
      '  nNcheck= ',num2str(nNcheck)]);
disp(['after:  nFtrain= ',num2str(kFtrain), ...
      '  nFcheck= ',num2str(kFcheck), ...
      '  nNtrain= ',num2str(kNtrain), ...
      '  nNcheck= ',num2str(kNcheck)]);
% ---------------------------------------------- preallocate data matrices
A0 = cell(1,kFtrain);              % Functional data matrices for training
A1 = cell(1,kNtrain);           % Nonfunctional data matrices for training
B0 = cell(1,kFcheck);               % Functional data matrices for testing
B1 = cell(1,kNcheck);            % Nonfunctional data matrices for testing
% -------------------------------------------------------- do the transfer
nSamples = t2 - t1 + 1;
% ---------------------------------------
k = 0;
  for s=1:nFtrain
    if( keepFtrain(s) > 0 )
    k = k + 1;
    A0{k} = zeros(nC,nSamples);
      for ic=1:nC
      tArray = tPeaksFtrain{ic,s}(1):tPeaksFtrain{ic,s}(end);
      L1 = ( tArray >= t1 );
      L2 = ( tArray <= t2 );
      L = and(L1,L2);
      A0{k}(ic,:) = level_Ftrain{ic,s}(L);
      end
    end
  end
nFtrain = k;
% ---------------------------------------
k = 0;
  for s=1:nFcheck
    if( keepFcheck(s) > 0 )
    k = k + 1;
    B0{k} = zeros(nC,nSamples);
      for ic=1:nC
      tArray = tPeaksFcheck{ic,s}(1):tPeaksFcheck{ic,s}(end);
      L1 = ( tArray >= t1 );
      L2 = ( tArray <= t2 );
      L = and(L1,L2);
      B0{k}(ic,:) = level_Fcheck{ic,s}(L);
      end
    end
  end
nFcheck = k;
% ---------------------------------------
k = 0;
  for s=1:nNtrain
    if( keepNtrain(s) > 0 )
    k = k + 1;
    A1{k} = zeros(nC,nSamples);
      for ic=1:nC
      tArray = tPeaksNtrain{ic,s}(1):tPeaksNtrain{ic,s}(end);
      L1 = ( tArray >= t1 );
      L2 = ( tArray <= t2 );
      L = and(L1,L2);
      A1{k}(ic,:) = level_Ntrain{ic,s}(L);
      end
    end
  end
nNtrain = k;
% ---------------------------------------
k = 0;
  for s=1:nNcheck
    if( keepNcheck(s) > 0 )
    k = k + 1;
    B1{k} = zeros(nC,nSamples);
      for ic=1:nC
      tArray = tPeaksNcheck{ic,s}(1):tPeaksNcheck{ic,s}(end);
      L1 = ( tArray >= t1 );
      L2 = ( tArray <= t2 );
      L = and(L1,L2);
      B1{k}(ic,:) = level_Ncheck{ic,s}(L);
      end
    end
  end
nNcheck = k;
tArray = t1:t2;
% ------------------------------------- handshake to Tyler's SPLOC machine
% OUTPUT
% A0 = cell(1,kFtrain);            % Functional data matrices for training
% A1 = cell(1,kNtrain);         % Nonfunctional data matrices for training
% B0 = cell(1,kFcheck);             % Functional data matrices for testing
% B1 = cell(1,kNcheck);          % Nonfunctional data matrices for testing
% tArray same for all data matrices: specifies time of measurement

%%                                                               functions
function [yLevel,nPeaks,tPeaks,hPeaks] = levelData(y,peakSign,debug)
% find R-peaks or peaks that appear to be R-peaks
%
% INPUT 
% y gives a datastreams
% peakSign = +/- one   (+ => normal, - => flip data before and after
% 
% PROCESS
% between two R-peaks, set the peaks to approximately the same height 
% subtract out a linear drift from the data between the two peaks.
%
% OUTPUT
% yLevel = transformed datastream after being leveled
% nPeaks = number of presumed R-peaks
% tPeaks = indices where the R-peak occurs in the datastream
% hPeaks = heights of the R-peaks in datastream after the trace is leveled
% ------------------------------ do not debug unless user request to do so
  if( nargin == 2 )
  debug = false;
  end
% ------------------------------------------------------------ switch sign
  if( peakSign < 0 )
  y = -y;
  end
% ------------------------------------------------------------- get hScale
w = 15;
jMin = w + 1;
jMax = length(y) - w;
hScale = zeros(size(y));
kk = 0;
  for j=jMin:jMax
  yWindow = y(j-w:j+w);
  kk = kk + 1;
  hScale(kk) = max(yWindow) - min(yWindow);
  end
kMax = kk;
hScale = sort(hScale(1:kMax),'descend');
% disp(hScale(1:100))
% pause 
identifyQRScomplex = true;
% ------------------------------------------- initialize hPeaks and tPeaks
isOkay = -1;
t = 1:length(y);
hPeaks = y + 0.001*rand(size(y));                % add noise to break ties
tPeaks = t;
  if( identifyQRScomplex )   %=> keep y that could belong to a QRS complex
    for kk=1:100:kMax
    hQRSmin = 0.4*hScale(kk);
    keepY = -ones(size(y));                         % => do not keep any y
      for j=jMin:jMax
      tempL = y(j) - y(j-w);
      tempR = y(j) - y(j+w);
        if( (tempL > hQRSmin) && (tempR > hQRSmin) )
        keepY(j) = 1;
        end
      end
    L = ( keepY > 0 );
      if( sum(L) > 1500 )
      %disp(num2str([sum(L),kk,hQRSmin]));
      hPeaks = hPeaks(L);
      tPeaks = tPeaks(L);
      isOkay = 1;
      break;
      end
    end
% else general case means do nothing: already start with everything
  end
% --------------------------------------------------------- verify or stop
  if( isOkay < 0 )
  disp('failed to obtain 1500 potential initial points to check')
  figure(10);
  clf;
  plot(t,y);
  pause
  end
% ------------------------------------------------------------------------
                      %  |||||------> input must set debug = true to debug
debugVisual = and(false, debug);  
debugList   = and(true , debug);
nPow = 0:8;     % ^^^^^---------------------> user can select how to debug
cutDtArray = [2.^nPow,275:25:350,360:10:480,485:5:510,511:525];
maxAveT = 800;
% -------------------------------------------- visualize the pluck process
        if( debugVisual )
        figure(10);
        clf; 
        plot(t,y,'b');
        hold on;
        plot(tPeaks,hPeaks,'or');
        xlabel('sample points');
        ylabel('measured signal');
        title('original signal with drift');
        pause;
        end
% --------------------------- remove all peaks that are not well separated
cut = 0;
    for cutDt=cutDtArray
    cut = cut + 1;
    j2 = 2;
    jE = length(tPeaks) - 1;
% ------------------------------------------------- remove irregular peaks
    flag2remove = -ones(size(tPeaks));          % initialize to no removal
    k = 0;                              % counter for points to be removed
      for j=j2:jE
      jL = j - 1;
      jR = j + 1;
      dtL = tPeaks(j) - tPeaks(jL);
      dtR = tPeaks(jR) - tPeaks(j);
        if( dtL < cutDt )              
        dh = hPeaks(jL) - hPeaks(j);                     % must avoid ties
          if( dh > 0.001 )
          flag2remove(j) = 1;
          elseif( dh < -0.001 )
          flag2remove(jL) = 1;
          %else skip ties
          end
        k = k + 1;
        end
        if( dtR < cutDt )
        dh = hPeaks(jR) - hPeaks(j);                     % must avoid ties
          if( dh > 0.001 )
          flag2remove(j) = 1;
          elseif( dh < -0.001 )
          flag2remove(jR) = 1;
          %else skip ties
          end
        k = k + 1;
        end
      end
% ----------------------------------------------- remove points when k > 0
      if( k > 0 ) 
      L = ( flag2remove < 0 );                          % keep these peaks
      hPeaks = hPeaks(L);
      tPeaks = tPeaks(L);
      hPeaks = hPeaks + 0.01*rand(size(tPeaks)); % add noise to break ties       
% -------------------------------------------- visualize the pluck process
        if( debugVisual )
        figure(10);
        clf; 
        plot(t,y,'b');
        hold on;
        plot(tPeaks,hPeaks,'or');
        xlabel('sample points');
        ylabel('measured signal');
        title('original signal with drift');
        pause;
        end
      end
% --------------------------------------------- check for a required break
    dtPeaks = sort(diff(tPeaks),'descend');
    m = min(90,length(dtPeaks));
    temp = 0.75*mean(dtPeaks(1:m));
      if( length(tPeaks) > 180 )       % => period is approximately 2xtemp
      s = (length(tPeaks) - 180)/60;
      temp = (1.2 + 0.8*s)*temp;
      end
    maxCutDt = max(temp,425);
    aveT = mean( diff(tPeaks(1:m) ) );
% --------------------------------------------- summarize list of numerics
      if( debugList && (k > 0) )
      nPts = length(tPeaks);
      disp(['cutDt= ',num2str(cutDt), ...
            '  pts-removed= ',num2str(k), ...
            '  maxCutDt= ',num2str(maxCutDt), ...
            '  aveT= ',num2str(aveT), ...
            '  maxAveT= ',num2str(maxAveT), ...
            '  nPts= ',num2str(nPts)]); 
      end
% -------------------------------------------------- check for early break
        if( ~identifyQRScomplex )        % => prone to greater # of errors
          if( (cutDt > maxCutDt) || (aveT > maxAveT) )
          break;       % better to identifyQRScomplex for greater accuracy
          end
        end
    end 
% -------------------------------------------------------- window the data
bMax = length(tPeaks) - 1;
tMinBlock = zeros(1,bMax);
tMaxBlock = zeros(1,bMax);
blockLpeak = zeros(1,bMax);
blockRpeak = zeros(1,bMax);
blockMean = zeros(1,bMax);
  for b=1:bMax
  tMinBlock(b) = tPeaks(b);
  tMaxBlock(b) = tPeaks(b+1);
  blockLpeak(b) = hPeaks(b);
  blockRpeak(b) = hPeaks(b+1);
  end
% ------------------------------------------------------- level each block
% figure(100)
% clf;
% hold on;
yShift = zeros(size(y));
  for b=1:bMax
  L1 = ( t > tMinBlock(b) );
  L2 = ( t < tMaxBlock(b) );
  L = and(L1,L2);
  yBlock = y(L);
  blockMean(b) = mean(yBlock);
  yShift(L) = yBlock - blockMean(b);
  % plot(t(L),y(L),'b');
  % plot(t(L),yShift(L),'r');
  % plot(tPeaks(b),hPeaks(b),'g*');
  % plot(tPeaks(b+1),hPeaks(b+1),'g*');
  % pause
  end
% ------------------------------------------------ level each peak divider
hPeaks = zeros( size(tPeaks) );
t1 = tPeaks(1);
yShift(t1) = y(t1) - blockMean(1);
hPeaks(1) = yShift(t1);
tE = tPeaks(end);
yShift(tE) = y(tE) - blockMean(end);
hPeaks(end) = yShift(tE);
b2 = 2;
  for b=b2:bMax
  yAve = 0.5*( blockRpeak(b-1) + blockLpeak(b) ) - blockMean(b);
  tt = tMinBlock(b);
  yShift(tt) = yAve;
  hPeaks(b) = yAve;
  end 
% --------------- reposition peaks due to uncertainty in left/right shifts
  for b=2:bMax-1
  t0 = tPeaks(b);
  [hMax,t1] = max(yShift(t0-40:t0+40));
  tMax = t0 - 40 + (t1 - 1);
  hPeaks(b) = hMax;
  tPeaks(b) = tMax;
  end
% ------------------------------------------- cut out first and last peaks
hPeaks = hPeaks(2:end-1);
tPeaks = tPeaks(2:end-1);
% -------------------------------------------- remove outlier peak heights
hPeaksSorted = sort(hPeaks,'descend');
m = ceil(0.9*length(hPeaksSorted));
aveH = mean(hPeaksSorted(1:m));
cutH = 0.5*aveH;
L = ( hPeaks > cutH );
hPeaks = hPeaks(L);
tPeaks = tPeaks(L);
% ------------------------------------------- finalize peak identification
nPeaks = length(tPeaks);
t1 = tPeaks(1);
tE = tPeaks(end);
yLevel = yShift(t1:tE);
  if( debugList )
  disp('  ');
  end
% ------------------------------------------------------------ switch sign
  if( peakSign < 0 )
  yLevel = -yLevel;
  hPeaks = -hPeaks;
  end
end
