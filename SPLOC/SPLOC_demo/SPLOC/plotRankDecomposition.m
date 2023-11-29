function plotRankDecomposition(systemRanking,verbosity,mPick)
% plot weight factor per basis vector used within the classification set
%
% INPUT
% systemRanking. <-- data structure
%     rankType = ('dataStream','moments')
%    baseFname = base file name useful for ID and spawning more file names
%    nXsystems = number of systems assigned a classification rank
%                X represents any system (0,1,unknown) that is ranked
%                1-system => from training set labeled as functional
%                0-system => from training set labeled as nonfunctional
% dataRefNameX = dataMatrixInfo.dataRefName    for rankType = 'dataStream'
%              = traitData.dataRefName         for rankType =    'moments'
%                X represents any system (0,1,unknown) that is ranked
% dataRefNameF = dataMatrixInfo.dataRefName    for rankType = 'dataStream'
%              = traitData.dataRefName         for rankType =    'moments'
%                Applies to the 1-systems that were used in a training set
%                F represents 1-systems that are regarded as functional
% dataRefNameN = dataMatrixInfo.dataRefName    for rankType = 'dataStream'
%              = traitData.dataRefName         for rankType =    'moments'
%                Applies to the 0-systems that were used in a training set
%                N represents 0-systems that are regarded as nonfunctional
%       nModes = # of discriminant modes used to classify unknown systems 
% ------------------------------------------------ list of cell arrays {:}
% systemXname1 = traitData.mMatrixName         for rankType =    'moments'
% systemXname2 = traitData.cMatrixName         for rankType =    'moments'
%                X represents any system (0,1,unknown) that is ranked
% -------------- list of arrays that reveal classification characteristics
%     pFave(j) = estimated Likelihood the j-th X-system is functional. 
%     pFstd(j) = standard deviation in the estimated likelihood.
%        wm(k) = average weight for k-th mode taken over all training data
% ---------------------------------------------------------- double arrays
%    pXF(j)(k) = MODE likelihood the j-th X-system is functional in the
%                k-th mode without comparing to known N-systems.
%    qXN(j)(k) = MODE likelihood the j-th X-system is not nonfunctional in
%                the k-th mode without comparing to known F-systems.
%    pqX(j)(k) = MODE likelihood the j-th X-systm is functional based on 
%                being similar to {F}-systems & dissimilar to {N}-systems
% ------------------------------------------------ discriminant SBV matrix
%            U = nV x nModes matrix containing nModes used to discriminate
%                data, where each mode lives in a nV dimensional space. 
% verbosity: {0,1,2,3} --> for plotting functions:
%             0,1 => print figures only to the screen.
%             2   => print figures to file only.
%             3   => print figures to file and screen. 
% ------------------------------------------------------------------------
% mPick = index for the X-system to plot the weights for. Without mPick,
%         the default is to plot the results for only the first X-system
%         in the list. Note that mPick corresponds to the rank order. 
%
% USAGE:                           
%                                   %  defaults: verbosity = 0   mPick = 1
%   plotRankDecomposition(systemRanking)          
%   plotRankDecomposition(systemRanking,verbosity)
%   plotRankDecomposition(systemRanking,verbosity,mPick)
%
% PROCESS
% plot the weight factors defined as wt.*pF for each mode in the order 
% they appear in the U matrix. It is possible that some of these vectors
% are dependent on other vectors because different modes can come from 
% different orthogonal sets of U all partitioning the vector space into a 
% discriminant subspace, but they are not identical, which effectively 
% expands the net span of the discriminant subspace and creates redundancy.
%
% OUTPUT (graphs to screen when verbosity = 0,1,3, and to file for 2,3)
% bar graph for the weight assigned to each mode.
%%                                             set global sploc parameters
global gvSPLOC                   % shares information across sploc toolset
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                             parse input
   if( nargin == 1 )
   verbosity = 0;                            % default => normal operation
   mPick = 1;
   elseif( nargin == 2 )
   verbosity = setVerbosity(verbosity); % returns allowed values {0,1,2,3}
   mPick = 1;
   else
       if( mPick < 1 )
       error('mPick must be > 0');
       end
       if( mPick > systemRanking.nXsystems )
       error(['mPick must not be greater than nXsystems = ', ...
             num2str(splocResults.nSols)]);
       end
   end
%%                                  extract information from datastructure
rankType = systemRanking.rankType;
baseFname = systemRanking.baseFname;
nModes = systemRanking.nModes;
systemXname1 = systemRanking.systemXname1{mPick};
pXFk = systemRanking.pXF(mPick,:);
qXNk = systemRanking.qXN(mPick,:);
pqXk = systemRanking.pqX(mPick,:);
% ---------------------------------
wmk = systemRanking.wm;
%%                                      set verbosity output level details
   if( verbosity > 1 )
   subFolder = 'classification';
   fName = [baseFname,'_',systemXname1,'_rankDecomp'];
   fName = getOutputFileName(subFolder,fName);
   end
   if( verbosity ~= 2 )
   figure;
   figure_number0 = get(gcf,'Number');
   end
%%                                                 plot rank decomposition
x = 1:nModes;
y_pXFk = 100*pXFk;
y_qXNk = 100*qXNk;
y_pqXk = 100*pqXk;
weight = 100*wmk;
minPFk = min(pqXk);
avePFk = mean(pqXk);
maxPFk = max(pqXk);
stdPFk = std(pqXk);
% ------------------------------------ prepare to plot mode decompositions
      if( verbosity == 2 )                    % prints graphs to file only
      rankDecomp = figure('visible','off');
      else
      rankDecomp = figure(figure_number0);        %  screen and maybe file
      end  
   clf
% --------------------------------------------------------- PROB(function)
   subplot(4,1,1);                               % apply title only at top
   bar(x,y_pqXk,'k'); 
   ylabel('function');
   ylim( [0,100] );
   xx = systemRanking.pFave(mPick); 
   str1 = [rankType,' ',num2str(mPick),': ',baseFname,'_', ... 
           systemXname1,'  pFunct= ',num2str(xx)];
   title(str1,'Interpreter','none');
% ------------------------------ PROB(to be similar to functional systems)
   subplot(4,1,2);
   bar(x,y_pXFk,'b'); 
   title('function similarity');
   ylabel('percent');
   ylim( [0,100] );
% ------------------------ PROB(to be dissimilar to nonfunctional systems)
   subplot(4,1,3);
   bar(x,y_qXNk,'r'); 
   title('nonfunction dissimilarity');
   ylabel('percent');
   ylim( [0,100] );
% ---------------------------------------- mode decomposition mode weights
   subplot(4,1,4);
   bar(x,weight,'m'); 
   xlabel('discriminant basis vector index');
   ylabel('percent');
   ylim( [0,100] );
   title(['mode weights for ',rankType]);  
      if( verbosity > 1 )
      saveas(rankDecomp,fName,gvSPLOC.gFileType);
      end
      if( verbosity == 3 )
      commandwindow
      figure(rankDecomp);
      pause
      close(rankDecomp);
      end
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['              reference name = ',baseFname]);
fprintf(fid,'%s \n',['   X-system being classified = ',systemXname1]);
fprintf(fid,'%s \n',['classifying method: rankType = ',rankType]);
fprintf(fid,'%s \n',['number of discriminant modes = ',num2str(nModes)]);
fprintf(fid,'%s \n',[' min basis vector likelihood = ',num2str(minPFk)]);
fprintf(fid,'%s \n',[' ave basis vector likelihood = ',num2str(avePFk)]);
fprintf(fid,'%s \n',[' max basis vector likelihood = ',num2str(maxPFk)]);
fprintf(fid,'%s \n',['   mode likelihood std. dev. = ',num2str(stdPFk)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end
