function plotCongruencySpectrum(splocResults,verbosity,dName)
% plot the selection power and concensus power summary bar graphs
%
% INPUT
% splocResults (a complicated data stucture)          
% verbosity:   {0,1,2,3} --> for plotting functions:
%              0,1 => print figures only to the screen.
%              2   => print figures to file only.
%              3   => print figures to file and screen. 
%
% OPERATIONAL CONTROL (using global variable switches)
% splocSplit = (-1,1) => (split graphs, no split graphs)
%                       split means to put indifference modes < 0
%                               and to put discriminant modes > 0
%                       No split means to put all modes > 0
%                       MCSPLOC requires split mode, SPLOC need not use it
% showGeometricMean = (0,1) => (no show, show) the geometric mean of the
%                              discriminant modes for visual reference
% ------------------------------------------------------------------------
%
% REMARK: Cind = congruency indicator (2,1,0,-1)
%                2 => discriminant and indifference congruences
%                1 => projections for discriminant congruences
%                0 => projections that are undetermined
%               -1 => projections for indifference congruences
%
% USAGE:                                             default verbosity = 0
%
%   plotCongruencySpectrum(splocResults)  
%   plotCongruencySpectrum(splocResults,verbosity)
%
% PROCESS
% Plot data using the standard split level format or the top format when
% the user specifies
%
% OUTPUT (graphs to screen when verbosity = 0,1,3, and to file for 2,3)
% two graphs, each with three panels stacked
% stacked bar graph 1: selection power, consensus power, quality factor
% stacked bar graph 2: efficacy, low & high selection within similarities
%%                                             set global sploc parameters
global gvSPLOC                   % shares information across sploc toolset
splocSplit = gvSPLOC.splocSplit;
showGeometricMean = gvSPLOC.showGeometricMean;
minScore1 = gvSPLOC.minScore1;       % score > minScore1 =>  "on" subspace
maxScore0 = gvSPLOC.maxScore0;       % score < maxScore0 => "off" subspace
minQF = gvSPLOC.qualityThreshold;      % minimum allowed quality threshold
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
% ---------------------------------------------------- dependent variables
lnMinScore1 = log(minScore1);
lnMaxScore0 = log(maxScore0);
lnMidScoreX = (lnMinScore1 + lnMaxScore0)/2;
%%                                                             parse input
 if( nargin == 1 )
 verbosity = 0;                              % default => normal operation
 elseif( nargin == 2 || 3)
 verbosity = setVerbosity(verbosity);   % returns allowed values {0,1,2,3}
 else
 disp('       ');
 disp('USAGE: ');
 disp('plotCongruencySpectrum(splocResults)');
 disp('plotCongruencySpectrum(splocResults,verbosity)'); 
 disp('                                    0,1 => plot to screen');
 disp('                                    2   => plot to file');
 disp('                                    3   => plot to file & screen');
 error('wrong usage/format');
 end
%%                                  extract information from datastructure
sType = splocResults.sType;                                % spectrum type
fname = splocResults.baseFname;
EEVd = splocResults.EEVd;
EEVi = splocResults.EEVi;
USVd = splocResults.USVd;
USVi = splocResults.USVi;
LSVd = splocResults.LSVd;
LSVi = splocResults.LSVi;
QEVd = splocResults.QEVd;
QEVi = splocResults.QEVi;
SEVd = splocResults.SEVd;
SEVi = splocResults.SEVi;
CEVd = splocResults.CEVd; 
CEVi = splocResults.CEVi;
Cind = splocResults.Cind;
vT = splocResults.vT;                                     % vote threshold
%%                                      set verbosity output level details
   if( verbosity > 1 )
   subFolder = dName;
   fNameC = [fname,'_congruencySpectrum'];
   fNameC = getOutputFileName(subFolder,fNameC);
   fName2 = [fname,'_propertiesSpectrum'];
   fName2 = getOutputFileName(subFolder,fName2);
   end
   if( verbosity ~= 2 )
   figure;
   figure_number0 = get(gcf,'Number');
   end
% ======================================================== old method hack
   if( strcmp(splocResults.sType,'SPLOC') && ...
       (splocSplit == 1) )
%%                                       map from new format to old format
EEV = zeros( size(EEVd) );
QEV = zeros( size(QEVd) );
USV = zeros( size(USVd) );
LSV = zeros( size(LSVd) );
SEV = zeros( size(SEVd) );
CEV = zeros( size(CEVd) );
CIP = Cind;
L3 = ( Cind > 0 );        % Note L3 = Ld for sploc  Dd = sum(Ld) = sum(L3)
D3 = sum(L3);
Ddi = -1;                                     % flags non-existence of Ddi
Ld = ( Cind == 1 );
CIP(Ld) = 1;
Lu = ( CIP == 0 );
Li = ( CIP == -1 );
Lupper = ( log(SEVd) > lnMidScoreX );
Llower = ~Lupper;
LuUpper = and(Lu,Lupper);
LuLower = and(Lu,Llower);
% ---------------------------
EEV(Ld) = EEVd(Ld);
EEV(LuUpper) = EEVd(LuUpper);
EEV(Li) = EEVi(Li);
EEV(LuLower) = EEVi(LuLower);
% ---------------------------
QEV(Ld) = QEVd(Ld);
QEV(LuUpper) = QEVd(LuUpper);
QEV(Li) = QEVi(Li);
QEV(LuLower) = QEVi(LuLower);
% ---------------------------
USV(Ld) = USVd(Ld);
USV(LuUpper) = USVd(LuUpper);
USV(Li) = USVi(Li);
USV(LuLower) = USVi(LuLower);
% ---------------------------
LSV(Ld) = LSVd(Ld);
LSV(LuUpper) = LSVd(LuUpper);
LSV(Li) = LSVi(Li);
LSV(LuLower) = LSVi(LuLower);
% ---------------------------
LSV(Ld) = LSVd(Ld);
LSV(LuUpper) = LSVd(LuUpper);
LSV(Li) = LSVi(Li);
LSV(LuLower) = LSVi(LuLower);
% ---------------------------
SEV(Ld) = SEVd(Ld);
SEV(LuUpper) = SEVd(LuUpper);
SEV(Li) = SEVi(Li);
SEV(LuLower) = SEVi(LuLower);
% ---------------------------
CEV(Ld) = CEVd(Ld);
CEV(LuUpper) = CEVd(LuUpper);
CEV(Li) = CEVi(Li);
CEV(LuLower) = CEVi(LuLower);
% ---------------------------
CEV(Ld) = CEVd(Ld);
CEV(LuUpper) = CEVd(LuUpper);
CEV(Li) = CEVi(Li);
CEV(LuLower) = CEVi(LuLower);
% ---------------------------
%%                      partition vector directions into congruency groups
% ------------------------------------------------------------ error check
Dd = sum(Ld);
Du = sum(Lu);
Di = sum(Li);
flagWrong = -1;                             % => to be wrong is impossible
   if( splocResults.Dd ~= Dd )
   flagWrong = 1;
   end
   if( splocResults.Du ~= Du )
   flagWrong = 1;
   end
   if( splocResults.Di ~= Di )
   flagWrong = 1;
   end
   if( flagWrong > 0 )
   disp( 'Cind is inconsistent with Dd, Du and Di');
   error('congruent subspace dimensions are inconsistent with Cind');
   end
dd = Dd + Du + Di;
% -------------------------------------------
   if( Dd > 0 )
   tempE = EEV(Ld);
   tempU = USV(Ld);
   tempL = LSV(Ld);
   tempQ = QEV(Ld);
   tempS = SEV(Ld); 
   tempC = CEV(Ld);
   [tempS,indx] = sort(tempS,'descend');
   dEfficacyFactor = tempE(indx);
   d_maxSimilarity = tempU(indx);
   d_minSimilarity = tempL(indx);
   dQualityFeature = tempQ(indx);
   dSelectionPower = tempS;
   dConsensusPower = tempC(indx);
   end
% -------------------------------------------
   if( Du > 0 )
   tempE = EEV(Lu);
   tempU = USV(Lu);
   tempL = LSV(Lu);
   tempQ = QEV(Lu);
   tempS = SEV(Lu);
   tempC = CEV(Lu);
   [tempS,indx] = sort(tempS,'descend');
   uEfficacyFactor = tempE(indx);
   u_maxSimilarity = tempU(indx);
   u_minSimilarity = tempL(indx);
   uQualityFeature = tempQ(indx);
   uSelectionPower = tempS;
   uConsensusPower = tempC(indx);
   end
% -------------------------------------------
   if( Di > 0 )
   tempE = EEV(Li);
   tempU = USV(Li);
   tempL = LSV(Li);
   tempQ = QEV(Li);
   tempS = SEV(Li);
   tempC = CEV(Li);
   [tempS,indx] = sort(tempS,'descend');
   iEfficacyFactor = tempE(indx);
   i_maxSimilarity = tempU(indx);
   i_minSimilarity = tempL(indx);
   iQualityFeature = tempQ(indx);
   iSelectionPower = tempS;
   iConsensusPower = tempC(indx);
   end
%%                                             parts of spectrum catenated
   Rank = 1:dd;
   dRank = 1:Dd;
   uRank = Dd+1:Dd+Du;
   iRank = Dd+Du+1:dd;
%%                  plot  selection - consensus - quality  spectrum bundle
      if( verbosity == 2 )                    % prints graphs to file only
      bCombo = figure('visible','off');
      else
      bCombo = figure(figure_number0);            %  screen and maybe file
      end  
   clf
%%                                                    plot selection power
   subplot(3,1,3);                                  
   hold on;
% ---------------------------------------------------------- place markers
   sI = 1:dd;
   SPbv = sqrt(minScore1*maxScore0);
   SPbvData = SPbv*ones(1,dd);
   xx1 = [0.6,dd+0.4,dd+0.4,0.6];                       % bin width is 0.8
   yy1 = [minScore1,minScore1,maxScore0,maxScore0]; 
   h5 = fill(xx1,yy1,[0.88,0.88,0.88],'LineStyle','none'); 
   flagDummyGot = -1;
   yMaxD = -1;
   yMaxU = -1;
   yMaxI = -1;
      if( Dd > 0 )
      yMaxD = ceil( max(dSelectionPower) );
      xData = SPbvData;
      xData(dRank) = dSelectionPower;
      h2 = bar(sI,xData,'y');
      h3 = bar(sI,xData,'b');
      h1 = bar(sI,xData,'r');
      h1(1).BaseValue = SPbv;
      flagDummyGot = 1;
      end
      if( Du > 0 )
      yMaxU = ceil( max(uSelectionPower) );
         if( flagDummyGot < 0 )
         xData = SPbvData;
         xData(uRank) = uSelectionPower;
         h1 = bar(sI,xData,'r');
         h3 = bar(sI,xData,'b');
         flagDummyGot = 1;
         end
      h2 = bar(uRank,uSelectionPower,'y');
      h2(1).BaseValue = SPbv;
      end
      if( Di > 0 )
      yMaxI = ceil( max(iSelectionPower) );
         if( flagDummyGot < 0 )
         xData = SPbvData;
         xData(iRank) = iSelectionPower;
         h1 = bar(sI,xData,'r');
         h2 = bar(sI,xData,'y');
         end
      xData = SPbvData;
      xData(iRank) = iSelectionPower;
      h3 = bar(sI,xData,'b'); 
      h3(1).BaseValue = SPbv;
      end
   xlabel('congruent basis vector index');
   ylabel('selection power');
   yMax = max( [3,yMaxD,yMaxU,yMaxI] );
   yMin = 1;
   ylim( [yMin,yMax] );
      if( (Dd > 0) && (showGeometricMean == 1) )
      aveLog1 = mean( log(dSelectionPower) );
      meanSelectionPower = exp(aveLog1);
      xt = [0.6,Dd+0.4];
      yt = [meanSelectionPower,meanSelectionPower];
      h4 = plot(xt,yt,'k','linewidth',2.0);
      legend([h4,h1,h2,h3,h5],{'geometric mean','discriminant', ...
      'undetermined','indifference','not-significant'},'Orientation','horizontal');
      else
      legend([h1,h2,h3,h5],{'discriminant','undetermined', ...
                            'indifference','not-significant'},'Orientation','horizontal');
      end
%%                                                    plot consensus power
   zData = zeros(1,dd);
   subplot(3,1,2);
   hold on;
   %xx1 = [0.6,d+0.4,d+0.4,0.6];             % seems like bin width is 0.8
   yy1 = [vT,vT,0,0]; 
   fill(xx1,yy1,[0.88,0.88,0.88],'LineStyle','none'); 
      if( Dd > 0 )
      rData = zData;
      rData(dRank) = dConsensusPower; 
      bar(sI,rData,'r');
      end
      if( Du > 0 )
      yData = zData;
      yData(uRank) = uConsensusPower; 
      bar(sI,yData,'y');
      end
      if( Di > 0 )
      bData = zData;
      bData(iRank) = iConsensusPower; 
      bar(sI,bData,'b');
      end
   ylabel('consensus power');
%%                                                   plot quality spectrum
   subplot(3,1,1);
   hold on;
   %xx1 = [0.6,d+0.4,d+0.4,0.6];             % seems like bin width is 0.8
   yy1 = [minQF,minQF,0,0]; 
   fill(xx1,yy1,[0.88,0.88,0.88],'LineStyle','none');
      if( Dd > 0 )
      rData = zData;
      rData(dRank) = dQualityFeature; 
      bar(sI,rData,'r');
      end
      if( Du > 0 )
      yData = zData;
      yData(uRank) = uQualityFeature; 
      bar(sI,yData,'y');
      end
      if( Di > 0 )
      bData = zData;
      bData(iRank) = iQualityFeature; 
      bar(sI,bData,'b');
      end
   ylabel('quality factor');
   totalEfficacy = splocResults.efficacy;
   totalEfficacy = round(totalEfficacy*10)/10;
   str1 = [fname,': ',sType,'-spectrum', ...
           '  net efficacy = ',num2str(totalEfficacy)];
   title(str1,'Interpreter','none');  
%%                                                       position subplots
   set(bCombo,'position',[400,250,520,600]); % [left bottom width height]
   sp1 = subplot(3,1,1);
   sp2 = subplot(3,1,2);
   sp3 = subplot(3,1,3);
   set(sp3,'position',[0.09 0.1 0.90 0.25]);
   set(sp2,'position',[0.09 0.4 0.90 0.25]);
   set(sp1,'position',[0.09 0.7 0.90 0.25]);
      if( verbosity > 1 )
      saveas(bCombo,fNameC,gvSPLOC.gFileType);
      end
% ========================================================================
%%         plot  efficacy - maxSimilarity - minSimilarity  spectrum bundle
      if( verbosity == 2 )                    % prints graphs to file only
      bCombo = figure('visible','off');
      else
      bCombo = figure(figure_number0 + 1);        %  screen and maybe file
      end  
   clf
%%                        plot maxSimilarity, (max - min) for discriminant
   xmin = 0.4;
   xmax = dd + 0.6;
   xBase = [xmin,xmax];
   yZero = [0,0];
% -------------------------------- determine which panel to use for legend
   Lupper = (SEV > SPbv);
   Llower = ~Lupper;
   LuUpper = and(Lu,Lupper);
   LuLower = and(Lu,Llower);
   indexBot = dd - min( [Rank(LuLower),Rank(Li)] );
   indexTop = max( [Rank(LuUpper),Rank(Ld)] );
      if( indexBot < indexTop )
      botLegend = 1;
      topLegend = 0;
      else
      botLegend = 0;
      topLegend = 1;
      end
   ytemp = 0.1*ones( size(SEV) );
   subplot(3,1,1);
   hold on;
% ------------------------------- plot dummy stuff for legend completeness
      if( topLegend == 1 )
      h1 = bar(1:dd,ytemp,'r');
      h2 = bar(1:dd,ytemp,'y');
      h3 = bar(1:dd,ytemp,'b');
      h4 = bar(1:dd,ytemp,'c');
      bar(1:dd,ytemp,'w','linewidth',1.5,'EdgeColor','w');
      end
% ---------------------------------------------------- plot the real stuff
      if( Dd > 0 )
      d_diffSimilarity = d_maxSimilarity - d_minSimilarity;
      rData = zData;
      cData = zData;
      rData(dRank) = d_maxSimilarity;
      cData(dRank) = d_diffSimilarity;
      bar(sI,rData,'r');
      bar(sI,cData,'c','linewidth',1.5,'EdgeColor','r');
      end
      if( Du > 0 )
      u_diffSimilarity = u_maxSimilarity - u_minSimilarity;
      yData = zData;
      cData = zData;
      yData( uRank(Lupper) ) = u_maxSimilarity(Lupper);
      cData( uRank(Lupper) ) = u_diffSimilarity(Lupper);
      bar(sI,yData,'y');
      bar(sI,cData,'c','LineWidth',1.5,'EdgeColor','y');
      end
   ylabel('within similarity');
   plot(xBase,yZero,'k');
      if( topLegend == 1)
      legend([h1,h2,h3,h4],{'max discriminant','max undetermined', ...
             'max indifference','(max - min)'},'location','northeast');
      end
   xlim([xmin,xmax]);
   str1 = [fname,': ',sType,'-spectrum', ...
           '  net efficacy = ',num2str(totalEfficacy)];
   title(str1,'Interpreter','none');
%%                        plot maxSimilarity, (max - min) for indifference
   subplot(3,1,2);
   hold on;
      if( botLegend == 1 )
      h1 = bar(1:dd,ytemp,'r');
      h2 = bar(1:dd,ytemp,'y');
      h3 = bar(1:dd,ytemp,'b');
      h4 = bar(1:dd,ytemp,'c');
      bar(1:dd,ytemp,'w','linewidth',1.5,'EdgeColor','w');
      end
      if( Di > 0 )
      i_diffSimilarity = i_maxSimilarity - i_minSimilarity;
      bData = zData;
      cData = zData;
      bData(iRank) = i_maxSimilarity;
      cData(iRank) = i_diffSimilarity;
      bar(sI,bData,'b');
      bar(sI,cData,'c','LineWidth',1.5,'EdgeColor','b');     
      end
      if( Du > 0 )
      u_diffSimilarity = u_maxSimilarity - u_minSimilarity;
      yData = zData;
      cData = zData;      
      yData(uRank) = u_maxSimilarity;
      cData(uRank) = u_diffSimilarity;
      bar(sI,yData,'y');
      bar(sI,cData,'c','LineWidth',1.5,'EdgeColor','y');
      end
   ylabel('within similarity');
   plot(xBase,yZero,'k');
   xlim([xmin,xmax]); 
      if( topLegend == 1)
      legend([h1,h2,h3,h4],{'max discriminant','max undetermined', ...
             'max indifference','(max - min)'},'location','northeast');
      end
%%                                                           plot efficacy
% check only
         Etemp = sum( EEVd( Cind == 1 ) ) ...
               + sum( EEVi( Cind == 2 ) ) ...
               + sum( EEVi( Cind == -1 ) );
            if( abs(Etemp - splocResults.efficacy) > 0.00001 )
            disp(['total Efficacy = ',num2str(Etemp), ...
                  ' by sum over EEVd and EEVi']);
            disp(['total Efficacy = ',num2str(splocResults.efficacy), ...
                  ' by prior summation']);
            error('error in sum rule detected!'); 
            end
   subplot(3,1,3);                                  
   hold on;
   flagDummyGot = -1;
      if( Dd > 0 )
      h3 = bar(dRank,dEfficacyFactor*0.1,'b');
      h2 = bar(dRank,dEfficacyFactor,'y');
      h1 = bar(dRank,dEfficacyFactor,'r');
      flagDummyGot = 1;
      end
      if( Du > 0 )
         if( flagDummyGot < 0 )
         h1 = bar(uRank,uEfficacyFactor*0.1,'r');
         h3 = bar(uRank,uEfficacyFactor*0.1,'b');
         flagDummyGot = 1;
         end
      h2 = bar(uRank,uEfficacyFactor,'y');
      end
      if( Di > 0 )
          if( flagDummyGot < 0 )
          h1 = bar(iRank,iEfficacyFactor*0.1,'r'); 
          h2 = bar(iRank,iEfficacyFactor*0.1,'y'); 
          end
      h3 = bar(iRank,iEfficacyFactor,'b');   
      end
   plot(xBase,yZero,'k');
   xlabel('congruent basis vector index');
   ylabel('efficacy');
   legend([h1,h2,h3],{'discriminant','undetermined', ...
                      'indifference'},'Location','southoutside','Orientation','horizontal');
   xlim([xmin,xmax]);
%%                                                       position subplots
   set(bCombo,'position',[920,250,520,600]); % [left bottom width height]
   sp1 = subplot(3,1,1);
   sp2 = subplot(3,1,2);
   sp3 = subplot(3,1,3);
   set(sp3,'position',[0.09 0.1 0.90 0.25]);
   set(sp2,'position',[0.09 0.4 0.90 0.25]);
   set(sp1,'position',[0.09 0.7 0.90 0.25]);
      if( verbosity > 1 )
      saveas(bCombo,fName2,gvSPLOC.gFileType);
      end     
% ======================================================== new method hack
   else                                    % allows for dual-purpose modes
% color scheme: r --> discriminant  b -> indifference   y --> undetermined
%                     gray = undetermined zone   black = thresholds/levels
% ------------------------------------------------------------------------
%     | R
%     | R R
%     | R R 
%     | R R R
%     | R R R 
%     | R R R R
%     | R R R R
%     |-R-R-R-R-------------Y-----Y------------- minScore1  Gray
%     | R R R R Y Y Y Y     Y Y   Y                         Gray
%     |=Y===========Y=Y=Y=Y=Y=Y===Y=Y=========== ScoreMid   Gray otherwise
%     | Y B   B     Y Y Y Y B B B B B B                     Gray
%     |---B---B-------Y---Y-B-B-B-B-B-B--------- maxScore0  Gray
%     |   B   B             B B B B B B
%     |       B               B B B B B
%     |___________________________B_B_B_________ 1
%
% The above notes and diagram defines the selection power indexing schema
% ------------------------------------------------------------------------
%  efficacy: put +Efficacy for d-subspace  and  -Efficacy for i-subspace
%   quality: put +Quality  for d-subspace  and  -Quality  for i-subspace
% consensus: put -vote     for d-subpace   and  -vote     for i-subspace
% ------------------------------------------------------------------------
% maxSimilaritity and minSimilarity
% For d-suspace:  make bar graph (red) for maxSimilarity
%                 as a sub-component (gray) fill bar with (max-min)
% For i-suspace:  make bar graph (blue) for maxSimilarity
%                 as a sub-component (gray) fill bar with (max-min)
% ------------------------------------------------------------------------
% ========================================================================
         Dd = splocResults.Dd;         % # of modes that discriminate only
         Ddi= splocResults.Ddi;             % # of modes with dual-purpose
         Du = splocResults.Du;          % # of modes that are undetermined
         Di = splocResults.Di;          % # of modes for indifference only
         D4 = Dd + Ddi + Du + Di; % 4 mode types, divided into 3-subspaces
%%                  plot  selection - consensus - quality  spectrum bundle
            if( verbosity == 2 )              % prints graphs to file only
            bCombo = figure('visible','off');
            else
            bCombo = figure(figure_number0);      %  screen and maybe file
            end  
         clf
         L3 = ( Cind > 0 );                          % captures Dd and Ddi
         L2 = ( Cind == 2 );         % Ddi: dual-purpose (d & i) subspaces
         Ld = ( Cind == 1 );                         % Dd: d subspace only
         Lu = ( Cind == 0 );                         % Du: u subspace only
         Li = ( Cind == -1);                         % Di: i subspace only
%        1st verification
%          disp( [L3; L2; Lu; Li; Cind] );
%          disp( dividerLine );
          
         D3 = sum(L3);
         sI = 1:D4;
         SPbv = exp(lnMidScoreX);
%%                                                    plot selection power
         selePlot = subplot(3,1,3);
         hold on;
         sPd = SEVd;   
         sPi = SEVi;
% ----------------------------------------------------------- build legend
         xmin = 0.4;
         xmax = D4 + 0.6;
         xx1 = [xmin,xmax,xmax,xmin];                % => bin width is 0.8
         ymax = 0.6*minScore1 + 0.4*maxScore0;
         ymin = 0.4*minScore1 + 0.6*maxScore0;
         yys = [ymax,ymax,ymin,ymin]; 
         h1 = fill(xx1,yys,'r','LineStyle','none');
         h2 = fill(xx1,yys,'y','LineStyle','none');
         h3 = fill(xx1,yys,'b','LineStyle','none');
         yys = [minScore1,minScore1,maxScore0,maxScore0]; 
         h4 = fill(xx1,yys,[0.88,0.88,0.88],'LineStyle','none'); 
            if( sum(L3) > 0 )
            xlim( [xmin,xmax] );
            aData = SPbv*ones( size(sI) );
            bData = aData;
            cData = aData;
            aData(L3) = sPd(L3);
            bData(L3) = sPi(L3);
            a1 = bar(sI,aData,'r');
            b1 = bar(sI,bData,'y');
            a1(1).BaseValue = SPbv;
            b1(1).BaseValue = SPbv;
               if( sum(L2) > 0 )
               cData(L2) = sPi(L2);
               c1 = bar(sI,cData,'b');
               c1(1).BaseValue = SPbv;

               bars.cData = cData;
               end
            end
            if( sum(Lu) > 0 )
            aData = SPbv*ones( size(sI) );
            bData = aData;
            aData(Lu) = sPd(Lu);
            bData(Lu) = sPi(Lu);
            a2 = bar(sI,aData,'y');
            b2 = bar(sI,bData,'y');
            a2(1).BaseValue = SPbv;
            b2(1).BaseValue = SPbv;
            end
            if( sum(Li) > 0 )
            aData = SPbv*ones( size(sI) );
            bData = aData;
            aData(Li) = sPd(Li);
            bData(Li) = sPi(Li);
            a3 = bar(sI,aData,'y');
            b3 = bar(sI,bData,'b');
            a3(1).BaseValue = SPbv;
            b3(1).BaseValue = SPbv;
            end
         ymax = max( ceil(sPd) );
         ymax = max(3,ymax);
         xlim( [xmin,xmax] );
         ylim( [1,ymax] );
         xlabel('selection index');
         ylabel('selection');
         
           if( (D3 > 0) && (showGeometricMean == 1) )
           aveLog1 = mean( log(sPd(L3) ) );
           yt1 = exp(aveLog1);
           xt = [0.6,D3+0.4];
           yt = [yt1,yt1];
           h0 = plot(xt,yt,'k','linewidth',2.0);
           legend([h0,h1,h2,h3,h4],{'geometric mean','discriminant', ...
           'undetermined','indifference','not-significant'}, ...
           'location','northeast','Orientation','vertical');
           else
           legend([h1,h2,h3,h4],{'discriminant','undetermined', ...
                           'indifference','not significant'}, ...
                           'location','northeast','Orientation','vertical');
           end
         hold off;
         % ----------------------------------------------------- consensus
         zData = zeros( size(sI) );
         subplot(3,1,2);
         hold on;
         yyv = [vT,vT,-vT,-vT]; 
         fill(xx1,yyv,[0.88,0.88,0.88],'LineStyle','none');
         cPd = CEVd;
         cPi = -CEVi;
            if( sum(L3) > 0 )
            rData = zData;
            yData = zData;
            bData = zData;
            rData(L3) = cPd(L3);
            yData(L3) = cPi(L3);
            bData(L2) = cPi(L2);
            bar(sI,rData,'r');
            bar(sI,yData,'y');
            bar(sI,bData,'b');
            end
            if( sum(Lu) > 0 )
            tData = zData;
            bData = zData;
            tData(Lu) = cPd(Lu);
            bData(Lu) = cPi(Lu);
            bar(sI,tData,'y');
            bar(sI,bData,'y');
            end
            if( sum(Li) > 0 )
            tData = zData;
            bData = zData;
            tData(Li) = cPd(Li);
            bData(Li) = cPi(Li);
            bar(sI,tData,'y');
            bar(sI,bData,'b');
            end
         ylabel('consensus');
         xlim( [xmin,xmax] );
         ylim( [-1,1] );
         yticks([-1 0 1]);
         yticklabels({'1','0','1'});
         hold off;
         % -------------------------------------------------- mode quality
         subplot(3,1,1);
         hold on;
         yyq = [minQF,minQF,-minQF,-minQF]; 
         fill(xx1,yyq,[0.88,0.88,0.88],'LineStyle','none');
         qPd = max(0,QEVd);
         qPi = -max(0,QEVi);
            if( sum(L3) > 0 )
            rData = zData;
            yData = zData;
            bData = zData;
            rData(L3) = qPd(L3);
            yData(L3) = qPi(L3);
            bData(L2) = qPi(L2);
            bar(sI,rData,'r');
            bar(sI,yData,'y');
            bar(sI,bData,'b');
            end
            if( sum(Lu) > 0 )
            tData = zData;
            bData = zData;
            tData(Lu) = qPd(Lu);
            bData(Lu) = qPi(Lu);
            bar(sI,tData,'y');
            bar(sI,bData,'y');
            end
            if( sum(Li) > 0 )
            tData = zData;
            bData = zData;
            tData(Li) = qPd(Li);
            bData(Li) = qPi(Li);
            bar(sI,tData,'y');
            bar(sI,bData,'b');
            end
         xlim( [xmin,xmax] );
         ylim( [-10,10] );
         yticks( [-10,-5,0,5,10] );
         yticklabels({'10','5','0','5','10'})
         ylabel('quality');
         totalEfficacy = splocResults.efficacy;
         totalEfficacy = round(totalEfficacy*10)/10;
         str1 = [fname,': ',sType,'-spectrum', ...
              '  net efficacy = ',num2str(totalEfficacy)];
         title(str1,'Interpreter','none'); 
%%                                                       position subplots
%                            [left bottom width height]
         set(bCombo,'position',[400,250,520,600]); 
         sp1 = subplot(3,1,1);
         sp2 = subplot(3,1,2);
         sp3 = subplot(3,1,3);
         set(sp3,'position',[0.09 0.1 0.90 0.25]);
         set(sp2,'position',[0.09 0.4 0.90 0.25]);
         set(sp1,'position',[0.09 0.7 0.90 0.25]);
            if( verbosity > 1 )
            saveas(bCombo,fNameC,gvSPLOC.gFileType);
            end
% ========================================================================
% {
%%         plot  efficacy - maxSimilarity - minSimilarity  spectrum bundle
            if( verbosity == 2 )              % prints graphs to file only
            bCombo = figure('visible','off');
            else
            bCombo = figure(figure_number0 + 1);  %  screen and maybe file
            end  
         clf
         xBase = [xmin,xmax];
         yZero = [0,0];
% ---------------------------------------------------------- mode efficacy
         % check only
         Etemp = sum( EEVd(L3) ) ...
               + sum( EEVi(L2) ) ...
               + sum( EEVi(Li) );
            if( abs(Etemp - splocResults.efficacy) > 0.00001 )
            disp(['total Efficacy = ',num2str(Etemp), ...
                  ' by sum over EEVd and EEVi']);
            disp(['total Efficacy = ',num2str(splocResults.efficacy), ...
                  ' by prior summation']);
            error('error in sum rule detected!'); 
            end
         subplot(3,1,3);
         hold on;
         ePd = max(0,EEVd);
         ePi = -max(0,EEVi);
            if( sum(L3) > 0 )
            rData = zData;
            yData = zData;
            bData = zData;
            rData(L3) = ePd(L3);
            yData(L3) = ePi(L3);
            bData(L2) = ePi(L2);
            bar(sI,rData,'r');
            bar(sI,yData,'y');
            bar(sI,bData,'b');
            end
            if( sum(Lu) > 0 )
            tData = zData;
            bData = zData;
            tData(Lu) = ePd(Lu);
            bData(Lu) = ePi(Lu);
            bar(sI,tData,'y');
            bar(sI,bData,'y');
            end
            if( sum(Li) > 0 )
            tData = zData;
            bData = zData;
            tData(Li) = ePd(Li);
            bData(Li) = ePi(Li);
            bar(sI,tData,'y');
            bar(sI,bData,'b');
            end
         ymax = max( ceil(0.001 + ePd) );
         ymin = min( floor(-0.001 + ePi) );
         ymax = ceil( 5*max(ymax, -ymin) )/5;
         ymax = min(100,ymax);
         ymin = -ymax;
            if( ymax < 5 )
            yticks( [-4,-2,0,2,4] );
            yticklabels({'4','2','0','2','4'})
            elseif( ymax < 10 )
            yticks( [-8,-4,0,4,8] );
            yticklabels({'8','4','0','4','8'}) 
            elseif( ymax < 20 )
            yticks( [-20, -10,  0, 10, 20] );
            yticklabels({'20','10','0','10','20'})
            elseif( ymax < 30 )
            yticks( [-25, -15, 0, 15, 25] );
            yticklabels({'25','15','0','15','25'})
            elseif( ymax < 40 )
            yticks( [-30, -20, -10, 0, 10, 20, 30] );
            yticklabels({'30','20','10','0','10','20','30'})
            elseif( ymax < 50 )
            yticks( [-40, -20, 0, 20, 40] );
            yticklabels({'40','20','0','20','40'})
            elseif( ymax < 60 )
            yticks( [-50, -25, 0, 25, 50] );
            yticklabels({'50','25','0','25','50'})
            elseif( ymax < 70 )
            yticks( [-60, -40, -20, 0, 20, 40, 60] );
            yticklabels({'60','40','20','0','20','40','60'})
            elseif( ymax < 80 )
            yticks( [-75, -50, -25, 0, 25, 50, 75] );
            yticklabels({'75','50','25','0','25','50','75'})
            elseif( ymax < 100 )
            yticks( [-90, -60, -30, 0, 30, 60, 90] );
            yticklabels({'90','60','30','0','30','60','90'})
            end
         plot(xBase,yZero,'k');
         ylim( [ymin,ymax] );
         xlim( [xmin,xmax] );
         ylabel('efficacy');
         hold off
% =========================================================== similarities
         Lupper = ( SEVd > SPbv );
         Llower = ~Lupper;
         LuUpper = and(Lu,Lupper);
         LuLower = and(Lu,Llower);
         B3 = or(Li,L2);                          % => captures Di and Ddi
% -------------------------------- determine which panel to use for legend
         indexBot = D4 - min( [sI(B3),sI(LuLower)] );
         indexTop = max( [sI(L3),sI(LuUpper)] );
            if( indexBot < indexTop )
            botLegend = 1;
            topLegend = 0;
            else
            botLegend = 0;
            topLegend = 1;
            end
% ------------------------ consider USVd & LSVd on L3=(Ld + Ldi) & LuUpper
         subplot(3,1,1);
         hold on;
         ytemp = 0.1*ones( size(USVd) );
            if( topLegend == 1 )
            ht1 = bar(sI,ytemp,'r');
            ht2 = bar(sI,ytemp,'y');
            ht3 = bar(sI,ytemp,'b');
            ht4 = bar(sI,ytemp,'c');
            bar(sI,ytemp,'w','linewidth',1.5,'EdgeColor','w');
            end
         dDif = USVd - LSVd;
         rData = zData;
         cData = zData;
         rData(L3) = USVd(L3);
         cData(L3) = dDif(L3);
         bar(sI,rData,'r');
         bar(sI,cData,'c','linewidth',1.5,'EdgeColor','r');
         yData = zData;
         cData = zData;
         yData(LuUpper) = USVd(LuUpper);
         cData(LuUpper) = dDif(LuUpper);
         bar(sI,yData,'y');
         bar(sI,cData,'c','linewidth',1.5,'EdgeColor','y');
         ylabel('within similarity');
         str1 = [fname,': ',sType,'-spectrum', ...
                 '  net efficacy = ',num2str(totalEfficacy)];
         title(str1,'Interpreter','none');
         plot(xBase,yZero,'k');
         xlim( [xmin,xmax] );
         ylim( [0,1] );
            if( topLegend == 1 )
            legend([ht1,ht2,ht3,ht4],{'max discriminant', ...
                'max undetermined','max indifference', ...
                '(max - min)'},'location','northeast');
            end
         hold off;
% ------------------------ consider USVi & LSVi on B3=(Li + Ldi) & LuLower
         subplot(3,1,2);
         hold on;
            if( botLegend == 1 )
            hb1 = bar(sI,ytemp,'r');
            hb2 = bar(sI,ytemp,'y');
            hb3 = bar(sI,ytemp,'b');
            hb4 = bar(sI,ytemp,'c');
            bar(sI,ytemp,'w','linewidth',1.5,'EdgeColor','w');
            end
         iDif = USVi - LSVi;
         yData = zData;
         cData = zData;
         yData(L3) = USVi(L3);
         cData(L3) = iDif(L3);
         bar(sI,yData,'y');
         bar(sI,cData,'c','linewidth',1.5,'EdgeColor','y');
         bData = zData;
         cData = zData;
         bData(B3) = USVi(B3);
         cData(B3) = iDif(B3);
         bar(sI,bData,'b');
         bar(sI,cData,'c','linewidth',1.5,'EdgeColor','b');
         yData = zData;
         cData = zData;
         yData(LuLower) = USVi(LuLower);
         cData(LuLower) = iDif(LuLower);
         bar(sI,yData,'y');
         bar(sI,cData,'c','linewidth',1.5,'EdgeColor','y');
         plot(xBase,yZero,'k');
         xlim( [xmin,xmax] );
         ylim( [0,1] );
         ylabel('within similarity');
            if( botLegend == 1 )
            legend([hb1,hb2,hb3,hb4],{'max discriminant', ...
                'max undetermined','max indifference', ...
                '(max - min)'},'location','northwest');
            end
         hold off;
%%                                                       position subplots
          set(bCombo,'position',[920, 250,   520,  600]); 
                              % [left bottom width height]
          sp1 = subplot(3,1,1);
          sp2 = subplot(3,1,2);
          sp3 = subplot(3,1,3);
          set(sp3,'position',[0.09 0.1 0.90 0.25]);
          set(sp2,'position',[0.09 0.4 0.90 0.25]);
          set(sp1,'position',[0.09 0.7 0.90 0.25]);
             if( verbosity > 1 )
             saveas(bCombo,fName2,gvSPLOC.gFileType);
             end     
   end
%%                                         record action in sploc log file
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['                reference name = ',fname]);
fprintf(fid,'%s \n',['                vote threshold = ',num2str(vT)]);
fprintf(fid,'%s \n',['                total efficacy = ', ...
                    num2str(totalEfficacy)]);
   if( D3 > 0 )
   aveQEVd = sum( QEVd(L3) )/D3;
   fprintf(fid,'%s \n',['mean quality for d-subspace = ', ...
                       num2str(aveQEVd)]);
   end
      if( Di > 0 )
      aveQEVi = sum( QEVi(Li) )/Di;
      fprintf(fid,'%s \n',['mean quality for i-subspace = ', ...
                          num2str(aveQEVi)]);
      end
fprintf(fid,'%s \n',['        discriminant modes: Dd = ',num2str(Dd)]);
  if( Ddi > -1 )
  fprintf(fid,'%s \n',['        dual-purpose modes: Ddi= ',num2str(Ddi)]);                                                        
  end
fprintf(fid,'%s \n',['        undetermined modes: du = ',num2str(Du)]);
fprintf(fid,'%s \n',['        indifference modes: di = ',num2str(Di)]);
   if( Ddi > -1 )
   fprintf(fid,'%s \n',[' # of variables = Dd+Ddi+Du+Di = ',num2str(D4)]);                      
   else
   fprintf(fid,'%s \n',[' # of variables = Dd + Du + Di = ',num2str(dd)]);
   end   
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end
