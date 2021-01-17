% test driver for workflow using the SPLOC toolset
% ------------------------------------------------------------------------
% SPLOC = Supervised Projective Learning with Orthogonal Completeness
%
% This test driver checks sploc functions in a way that represents a 
% possible workflow. This is called workflow 1, as various workflows are
% possible. These workflows can serve as starter-scripts to solve actual
% problems related to protein dynamics. The various workflow scripts could
% even be the same as the ones used in drug/protein/solvent design.   
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
% INPUT: synthetic 2D molecules
%
% PROCESS: Various routes of interest. If interesting things diverge too
% much, make new script with different workflow. 
%
% OUTPUT: stuff
%%                                                      start from scratch
close all
clear all
clc 
disp('  ');
disp(mfilename);
disp( dividerLine('select training set') );
disp('1. trainingSet1 ')
disp('2. trainingSet1B')
disp('  ');
nMenu = input('   Enter option: ');
   switch nMenu
       case 1
       fName = 'trainingSet1';
       Ns = 500;
       case 2
       fName = 'trainingSet1B';
       Ns = 20000;
       otherwise
       error('unknown option');
   end
%%                                                           perform setup
disp('  ');
prefix = input('enter name of scenario: ','s');
   if( length(prefix) < 2 )
   error('name of scenario is too short');
   else
   prefix = [prefix,'_'];
   end
logFile = [prefix,mfilename];
disp('   ');
disp(logFile);
initializeSPLOC(1,'fName',logFile,'gType','png');
mFormat = setDataMatrixFormat('xxx-yyy-zzz',2);  % by-pass screen input
mType = 'cov';
splocResults = cell(1,24);
% ------------------------------------------------- determine training set
[~,cLabel,FnameAdataU] = readFileNameList(fName,0,'cType','U');   
FnameAdataX = FnameAdataU;
%%                                           prepare for iterative process
nF = 0;
nN = 0;
nU = 24;
% ------------------------------------- summarize initial unclassifed list
disp('   ');
disp(  dividerLine('Unclassified simulation data') );
g = 1:nU;
T = table(g',FnameAdataU','variableNames',{'count' 'Usystem'});
disp(T);
% -------------
SBV = 0;                                              %  initial basis set
flagF = -ones(1,24);
flagN = -ones(1,24);
flagFINISH = -1;
iter = 0;
%%                                              perform iterative workflow
 while( nU > 0 )
 iter = iter + 1;
%% --------------------------------------------------- modify training set
 flag_finish = 1;
   while( flag_finish ~= 0 )
   disp('  ');
      if( flag_finish == 1 )
      disp('Experiments performed:');
      disp( dividerLine('molecules exhibiting FUNCTION') );
      disp('{-1,0} => {FINISH process, STOP adding molecules}');
      disp('   ');
      else
      fprintf(2,'Remove mistakes!\n')
      %disp('Remove mistakes!');
      disp( dividerLine('molecules exhibiting FUNCTION') );
      disp('0 => STOP removing molecules');
      end
      while(1)
      moleculeName = input('Enter F-molecule name: ','s');
         if( strcmp(moleculeName,'-1') == 1)
         moleculeName = '0';
            if( nN > 0 )
            flagFINISH = 1;
            end
         end
         if( strcmp(moleculeName,'0') == 1 )
             if( nF > 0 )
             break;
             end
         else
            if( length(moleculeName) ~= 3 )
            moleculeName = 'x';
            end
            if( flag_finish == 1 )
               for ju=1:24
                  if( strcmp(cLabel{ju},'U') == 1 )
                     if( contains(FnameAdataX{ju},moleculeName) )
                     cLabel{ju} = 'F';
                     flagF(ju) = 1;
                     nU = nU - 1;
                     nF = nF + 1;
                     flagSkipped = -1;
                     break;
                     else
                     flagSkipped = 1;
                     end
                  end
               end
            else
               for jF=1:24
                  if( strcmp(cLabel{jF},'F') == 1 )
                     if( contains(FnameAdataX{jF},moleculeName) )
                     cLabel{jF} = 'U';
                     flagF(jF) = -1;
                     nU = nU + 1;
                     nF = nF - 1;
                     flagSkipped = -1;
                     break;
                     else
                     flagSkipped = 1;
                     end
                  end
               end
            end
            if( flagSkipped > 0 )
            disp('choice skipped');
            end
         end
      end
      if( flagFINISH > 0 )
      break;
      end
% ------------------------------------------------------------------------
   disp('  ');
   disp( dividerLine('molecules exhibiting NONFUNCTION') );
   disp('0 => stop');
      while(1)
      moleculeName = input('Enter N-molecule name: ','s');
         if( strcmp(moleculeName,'0') == 1 )
             if( nN > 0 )
             break;
             end
         else
            if( length(moleculeName) ~= 3 )
            moleculeName = 'x';
            end
            if( flag_finish == 1 )
               for ju=1:24
                  if( strcmp(cLabel{ju},'U') == 1 )
                     if( contains(FnameAdataX{ju},moleculeName) )
                     cLabel{ju} = 'N';
                     flagN(ju) = 1;
                     nU = nU - 1;
                     nN = nN + 1;
                     flagSkipped = -1;
                     break;
                     else
                     flagSkipped = 1;
                     end
                  end
               end
            else 
               for jN=1:24
                  if( strcmp(cLabel{jN},'N') == 1 )
                     if( contains(FnameAdataX{jN},moleculeName) )
                     cLabel{jN} = 'U';
                     flagN(jN) = -1;
                     nU = nU + 1;
                     nN = nN - 1;
                     flagSkipped = -1;
                     break;
                     else
                     flagSkipped = 1;
                     end
                  end
               end
            end
            if( flagSkipped > 0 )
            disp('choice skipped');
            end
         end
      end
% ---------------------------------------------- get F-molecule file names
   FnameAdataF = cell(1,nF);
   k = 0;
      for j=1:24
         if( flagF(j) > 0 )
         k = k + 1;
         FnameAdataF{k} = FnameAdataX{j};  
         end
      end
   if( k == 0 )
   error('No molecules that function are present');
   end
% ---------------------------------------------- get N-molecule file names
   FnameAdataN = cell(1,nN);
   k = 0;
      for j=1:24
         if( flagN(j) > 0 )
         k = k + 1;
         FnameAdataN{k} = FnameAdataX{j};  
         end
      end
   if( k == 0 )
   error('No molecules that do not function are present');
   end
% ---------------------------------------------- get U-molecule file names
   FnameAdataU = cell(1,nU);
   k = 0;
      for j=1:24
         if( (flagN(j) < 0) && (flagF(j) < 0) )
         k = k + 1;
         FnameAdataU{k} = FnameAdataX{j};
         end
      end
   disp( dividerLine('molecules with FUNCTION') );
   g = 1:nF;
   T = table(g',FnameAdataF','variableNames',{'count'  'Fsystem'});
   disp(T);
   disp( dividerLine('molecules with NONFUNCTION') );
   g = 1:nN;
   T = table(g',FnameAdataN','variableNames',{'count'  'Nsystem'});
   disp(T);
      while(1)
      msg = 'Enter -100, -1, 0, 1 to (FLIP, remove, accept, add): ';
      flag_finish = input(msg);
% ------------------------------------------------------------------- FLIP
         if( flag_finish == -100 )
         FnameAdataN = cell(1,nF);
         FnameAdataF = cell(1,nN);
         SBV = 0;                       % restart from scratch with a flip
         temp = flagN;
         flagN = flagF;
         flagF = temp;
         nF = 0;
         nN = 0;
            for j=1:24
               if( flagN(j) > 0 )
               nN = nN + 1;
               FnameAdataN{nN} = FnameAdataX{j};
               elseif( flagF(j) > 0 )
               nF = nF + 1;
               FnameAdataF{nF} = FnameAdataX{j};
               end
            end
         disp('  ');
         disp( dividerLine('molecules with FUNCTION') );
         g = 1:nF;
         T = table(g',FnameAdataF','variableNames',{'count'  'Fsystem'});
         disp(T);
         disp( dividerLine('molecules with NONFUNCTION') );
         g = 1:nN;
         T = table(g',FnameAdataN','variableNames',{'count'  'Nsystem'});
         disp(T);
         end
      flag_finish = floor(flag_finish);
         if( abs(flag_finish) < 2 )
         break;
         end
      end
   end
   if( flagFINISH > 0 )
   break;
   else
%%                                       collect simulation data for SPLOC
   refNameF = [prefix,'F',num2str(iter,'%02i')];
   [AmatrixInfoF,tableF] = readDataMatrices(refNameF,FnameAdataF,mFormat);                                                
   disp('  ');
   disp(tableF);
   disp('   ');
   refNameN = [prefix,'N',num2str(iter,'%02i')];
   [AmatrixInfoN,tableN] = readDataMatrices(refNameN,FnameAdataN,mFormat);                                                
   disp('  ');
   disp(tableN);
      if( nU > 0 )
      disp('   ');
      refNameU = [prefix,'U',num2str(iter,'%02i')];
      [AmatrixInfoU,tableU]= readDataMatrices(refNameU, ...
                                              FnameAdataU,mFormat);
      disp('  ');
      disp(tableU);
      end
   disp('   ');
   refNameX = [prefix,'X',num2str(iter,'%02i')];
   [AmatrixInfoX,~] = readDataMatrices(refNameX,FnameAdataX,mFormat);
   disp('  ');
%    disp(tableX);
%    disp('   ');
   pause(2);
%%                                                           create traits
   ns2 = floor(Ns/2);              %=> split data into 2 groups per system
   traitF2 = getMultivariateStats4sploc(AmatrixInfoF,ns2,0,mType);
   traitN2 = getMultivariateStats4sploc(AmatrixInfoN,ns2,0,mType);
   %traitU2 = getMultivariateStats4sploc(AmatrixInfoU,ns2,0,mType);
   close all
%%                                                                SPLOC it
   splocName = [prefix,'sploc',num2str(iter,'%02i')];
   [ma,mb] = size(SBV);
      if( (ma == 1) && (mb == 1) )                   % => no guess for SBV
      splocResult = sploc(0,0,splocName,traitF2,traitN2,1);
                                                      % 1
      else                                          % => have previous SBV
      splocResult = sploc(0,SBV,splocName,traitF2,traitN2,1);
      %splocResult = sploc(-1,0,splocName,traitF2,traitN2,1);
      end
   
   
   splocResults{iter} = splocResult;
   SBV = splocResult.SBV;
   plotCongruencySpectrum(splocResult,3);
   Ud = getDiscriminantSBV(splocResult);         % discriminant modes only
   [~,jjj] = size(Ud);
      if( jjj < 1 )
      disp('NO discriminant modes: Cannot classify!!!!');
      end
   pause(2);
%%                                          classify U-molecules as F or N
      if( jjj > 0 )
      disp('   ');
      msg = 'classify most updated training data';
      disp( dividerLine(msg) );
      rankName = [prefix,'rankX',num2str(iter,'%02i')];
      [rankX,T] = classifyDataStream(rankName,AmatrixInfoX, ... 
                                     AmatrixInfoF,AmatrixInfoN,Ud,0);
                                         % plot visual results => 3
                               % skip visual results on screen => 2
%       traitX2 = getMultivariateStats4sploc(AmatrixInfoX,ns2,0,mType);
%       [rankX,T] = classifyMoments(rankName,traitX2, ...
%                                   traitF2,traitN2,Ud,3);
      disp('   ');
      disp(T);
      disp('   ');
% {
      disp( dividerLine('visualize mode contributions') );
      m1 = max(1,nF - 2);
      m2 = min(24,nF + 2);
         for m=m1:m2
         plotRankDecomposition(rankX,0,m);
         pause(5)
         end
%}
         if( nU > 0 )
         rankName = [prefix,'rankU',num2str(iter,'%02i')];
         [~,T] = classifyDataStream(rankName,AmatrixInfoU, ... 
                                    AmatrixInfoF,AmatrixInfoN,Ud,1);
%          traitU2 = getMultivariateStats4sploc(AmatrixInfoU,ns2,0,mType);
%          [rankU,T] = classifyMoments(rankName,traitU2, ...
%                                   traitF2,traitN2,Ud,3);


         disp('   ');
         disp(T);
         end

%%                                                 plot feature clustering
% {
      traitX2 = getMultivariateStats4sploc(AmatrixInfoX,ns2,0,mType);
      featureF2 = getFeatureVectors(traitF2,Ud);
      featureN2 = getFeatureVectors(traitN2,Ud);
      featureX2 = getFeatureVectors(traitX2,Ud);
% featureX. <-- data structure
% dataRefName = reference name for data with similar traits for sploc
% mMatrixName = cell array for file names that store the mMatrix data
% cMatrixName = cell array for file names that store the cMatrix data
%   nXsystems = number of systems being projected into feature space
%               X represents any system (0,1,unknown) that is ranked
%                 1-system => from training set labeled as functional
%                 0-system => from training set labeled as nonfunctional
%      nModes = # of discriminant modes contained in U.
%   nFeatures = 2xnModes = number of distinct features
%     Fmatrix = nFeatures x nXsystems 
      FmatrixF = featureF2.Fmatrix;
      FmatrixN = featureN2.Fmatrix;
      FmatrixX = featureX2.Fmatrix;
      jStDv = 0;
         for k=1:featureX2.nModes
         figure(100 + k);
         clf;
         hold on;
         jMean = jStDv + 1;
         jStDv = jMean + 1;
         % -----------------------------------
         x1 = FmatrixX(jMean,:);
         x2 = FmatrixX(jStDv,:);
         scatter(x1,x2,25,'k');
         % -----------------------------------
         x1 = FmatrixF(jMean,:);
         x2 = FmatrixF(jStDv,:);
         scatter(x1,x2,75,'filled','b');
         % -----------------------------------
         x1 = FmatrixN(jMean,:);
         x2 = FmatrixN(jStDv,:);
         scatter(x1,x2,75,'filled','r');
         xlabel('mean');
         ylabel('standard deviation');
         title(['Feature space for mode: ',num2str(k)]);
         pause
         end
%}
% ========================================================================
      end
   end
 end
