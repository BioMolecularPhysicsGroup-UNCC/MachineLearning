function [splocResults,mapU2S,mapS2U] = ...
         getBVSmultiClass(U,baseFname,traitL,cID,vT)
% Given a basis set of vectors determine the MCSPLOC congruency spectrum
% Different types of inputs and outputs are possible
%
% ------------------ three conditions -----------------       congruency
% score > minScore1  &  dVote > vT  &  dQuality > minQF  =>  discriminant
% score < maxScore0  &  iVote > vT  &  iQuality > minQF  =>  indifference
%                otherwise                               =>  undetermined
% 
% DEFINITIONS
% ------------------------------------------------------- trait definition
% trait => data structure  (required components)
%     n = # of statistical metrics
%    mu = cell array containing n vectors for a characteristic property.
%    cM = cell array containing information on pairwise correlations that
%         are symmetric (cM is any symmetric matrix deemed useful).
%    ---> nVariables and nDtotal = total # of data samples is information
%         needed to calculate the vote threshold using a heristic formula.
% -------------------------------------------------------- local variables
% INPUT:
% U = proposed set of basis vectors
% baseFname => all output file names spawn off this base file name
% ----------------------------------------------------------------- traitL
% nL = # of statistical metrics for unlabeled systems
% ML = cell array containing n0 mean vectors for unlabeled systems
% QL = cell array containing n0 covariance matrices for unlabeled systems
% nVL = # of variables in each unlabeled system
% nDL = total # of data samples = n0*samplesize0
% --------- 
% vT = optional vote threshold   (default is calculated based on sampling)
%
% USAGE:
% splocResults = getBVSmultiClass(U,fname,3,traitL,cID)
% splocResults = getBVSmultiClass(U,fname,3,traitL,cID,vT)
%
%         Notes: traitL defines classified systems => class labels = cID
%                  cID = list of numbers from 1 to nC = # of class labels
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% PROCESS
% Emmulate the process found in mcsploc() identically, except this is a
% one-time shot that quantifies the input vectors. No optimization is 
% performed. Furthermore, the process within   getBVStwoClasses()  is done 
% here identically, except now it is done for all pairs of classes.
%
% OUTPUT 
% splocResults.           => data structure
% splocResults.sType       = type of spectrum: MCSPLOC
% splocResults.pursuitType = 1,0,-1  => d, d&i, i  set DEFAULT VALUE = 0
% splocResults.pType       = packing format describing the vector space
% splocResults.dim         = # of components in a local vector
% splocResults.baseFname   = base file name
% splocResults.SBV         = selection basis vectors
% splocResults.vT          = voting threshold to establish consensus
% splocResults.efficacy    = quantifies the ability for SBV to cluster
% splocResults.Dd          = # of discriminant only modes
% splocResults.Ddi         = # of (discriminant & indifference) modes
% splocResults.Di          = # of indifference only modes
% splocResults.Du          = # of undetermined modes
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% splocResults.EEVd = efficacy eigenvalues for discriminant congruence
% splocResults.EEVi = efficacy eigenvalues for indifference congruence 
% splocResults.USVd = conditional upper similar discriminant-eigenvalues
% splocResults.USVi = conditional upper similar indifference-eigenvalues
% splocResults.LSVd = conditional lower similar discriminant-eigenvalues
% splocResults.LSVi = conditional lower similar indifference-eigenvalues
% splocResults.QEVd = conditional quality discriminant-eigenvalues
% splocResults.QEVi = conditional quality indifference-eigenvalues
% splocResults.SEVd = conditional selection discriminant-eigenvalues
% splocResults.SEVi = conditional selection indifference-eigenvalues 
% splocResults.CEVd = conditional consensus discriminant-eigenvalues 
% splocResults.CEVi = conditional consensus indifference-eigenvalues 
% splocResults.Cind = congruency indicator (2,1,0,-1)
%                      2 => discriminant and indifference congruences
%                      1 => projections for discriminant congruences
%                      0 => projections that are undetermined
%                     -1 => projections for indifference congruences
%            mapU2S = map original U indices to SPLOC mode indices
%                  => Sindex = mapU2S(Uindex)
%            mapS2U = map SPLOC mode indices to original U indices
%                  => Uindex = mapS2U(Sindex)
%%                                               check number of arguments
% ---------------------------------------- work with specific format cases
   switch nargin
     case 4
     useDefault_vT = true;   % default voting threshold will be calculated   
     case 5
     useDefault_vT = false;
   end  
%%                                             simple error check on input
   if( useDefault_vT == false )
      if( vT < 0.5 )
      error('A voting threshold < 50% is not a majority vote!');
      elseif( vT > 0.95 )
      error('95% is the maximum voting threshold used in sploc!');
      end
   voteThreshold = vT;
   end 
%%                                 internally renumber static class labels
newCL2oldCL = unique(cID);                       % Note: CL = class labels
nCL = length(newCL2oldCL);                             % # of class labels
   if( nCL < 2 )
   error('nCL < 2: must have at least 2 labeled classes to compare');
   end
%%                                                extract data from traitL
nL = traitL.n;
Mgot = traitL.mu;                                                  % => ML
Qgot = traitL.cM;                                                  % => QL
nVL= traitL.nVariables;                                   % # of variables
NsL= traitL.sampleSize;
nt = nL;
% --------------------------------------------- overload normal SPLOC case
Ns0 = NsL; 
Ns1 = NsL;
%%                                                  error check input data
   if( length(cID) ~= nL )
   error('array length of cID is not equal to # of systems in traitL');
   end
   if( nL < nCL )
   error('nL < nCL: must have at least 1 system per class to compare');
   end
nDOF = nVL;
% -------------------------------------------------- check properties of U
[mV,nModes] = size(U);
   if( mV ~= nDOF )
   error('number of DOF for U is not equal to number of variables');
   end
% nModes = # of modes based on user's request. Typically nModes = nDOF
nC = nCL;                         % use simpler name: => number of classes
nCpairs = nC*(nC - 1)/2;                  % number of distinct class-pairs
%%                                       assign systems to labeled classes
classSystemList = cell(1,nC);
nSystemsPerClass = zeros(1,nC);             % # of systems per class label
sampleList = 1:nt;
   for newCL=1:nC
   oldCL = newCL2oldCL(newCL);
   Lkeep = ( cID == oldCL );
   classSystemList{newCL} = sampleList(Lkeep);
   nSystemsPerClass(newCL) = sum(Lkeep);
   %disp( [newCL,nSystemsPerClass(newCL),classSystemList{newCL}] ); 
   end
%%                            assign probability weight to each class pair
pcwt = zeros(1,nCpairs);
cPair = 0;
   for cL0=1:nC    % start double-loop over all pairs of clusters
   n0 = nSystemsPerClass(cL0);
      for cL1=cL0+1:nC
      n1 = nSystemsPerClass(cL1);
      cPair = cPair + 1;
      pcwt(cPair) = n0*n1;
      end
   end
total_pcwt = sum(pcwt);
pcwt = pcwt/total_pcwt;
%disp( pcwt );
%%                                                calculate vote threshold
% REMARK: This is complicated because based on the heristic formula given
% for voteThresold in sploc(), which is a binary classifer, the number of 
% systems in state-0 (n0) and in state-1 (n1) are accounted for. Since 
% these numbers are changing across different pairwise comparisons, the
% voting threshold should change to accommodate variations in n0 and n1.
% Even if this was calculated on the fly, it would be hard to plot the 
% results. Therefore, the approach taken here as the default evaluation is 
% to calculate an average vote theshold over the static cluster labels. If
% a user inputs a vote threshold, this single number is also applied to
% the entire dataset. In any case, only one vote threshold is considered 
% across all pairs of classes.
   if( useDefault_vT )   % default uses heristic formula for voteThreshold
   % ------------------- calculate average vT over labeled pair of classes
   vTsum = 0;
     for cL0=1:nC
     n0 = nSystemsPerCL(cL0);
        for cL1=cL0+1:nC
        n1 = nSystemsPerCL(cL1);
        % ------------------------------------------- calculate average vT
        nS0 = Ns0*sqrt(n0);                   % nS0 = sqrt(n0)*sampleSize0
        nS1 = Ns1*sqrt(n1);                   % nS1 = sqrt(n1)*sampleSize1
        voteThreshold = getDefaultVoteThreshold(nS0,nS1,nDOF);
        vTsum = vTsum + voteThreshold;
        end
     end
   voteThreshold = vTsum/nCpairs;
   end
%%                                             initialize SPLOC parameters
global gvSPLOC                         % shared information across toolset
qType = gvSPLOC.qType;
% -------------------------------------------- scoring function parameters
minScore1 = gvSPLOC.minScore1;       % score > minScore1 =>  "on" subspace
maxScore0 = gvSPLOC.maxScore0;       % score < maxScore0 => "off" subspace
minQF = gvSPLOC.qualityThreshold;      % minimum allowed quality threshold
add0 = gvSPLOC.add0;                 % prevents variance ratios to diverge
% ------------------------------------------------ bias on cluster quality
fIm = gvSPLOC.fIm;
fDm = gvSPLOC.fDm;
%disp( [fIm,fDm] );
% ---------------------------------------------------- dependent variables
lnMinScore1 = log(minScore1);
lnMaxScore0 = log(maxScore0);
lnMidScoreX = (lnMinScore1 + lnMaxScore0)/2;
% ---------------------------------------------------- efficacy parameters
sVS = gvSPLOC.sVS;      % reduction in vecScore for i-modes w.r.t. d-modes
%%                                             Initialize weight functions
dWt = gvSPLOC.dWtS;               % discriminant probability weight factor
iWt = gvSPLOC.iWtS;               % indifference probability weight factor
%uWt = gvSPLOC.uWtS;              % undetermined probability weight factor
xMax = gvSPLOC.xMax;                     % maximum value of score required
dampLR = 1; 
%%                              initialize for cluster properties measures
rFeatureSoft = 0.2;        % set value of rFeature when using softGeometry   
  switch qType
     case 'hardLinear'
     % no special initialization is necessary
     case {'softGeometry','hybrid'}
% ---------------------- these parameters have been subjectively optimized
% ---------------------- these parameters have been subjectively optimized
     haloDiff = 0.1;       % + 0.2*4/sqrt(n0*n1);           % ~0.10 FIDDLE
     haloSame = 0.1*haloDiff;                               % ~0.01 FIDDLE
%    REMARK: haloSame larger than 0.02 gives more problems more frequently
%            Not sure why, but in theory it should be zero for infinitely
%            many samples. Thus too large is presumably not good. haloDiff
%            greater than 0.3 does not work well. 0.5 almost always does 
%            not work. As haloDiff gets larger, the results look more like
%            hardLinear. When haloDiff gets really small, and should be 
%            about 5 times bigger (or more) than haloSame, clustering for
%            small numbers of samples looses all meaning. Since the data
%            fits in a unit square, 0.1 and 0.01 give reasonable scales.
%            if these numbers are good enough to be universal this would
%            be great. I tried to modify these numbers based on n0 and n1,
%            but the simple idea of setting both as constants worked best.
% ------------------------------------------------------------------------
     log10QdThreshold = 0.25;
     log10QiThreshold = 0.25;
     minQ = 0.001;
     maxQ = 10;
% ------------------------------------------------------------------------
     log10QdMIN = 1.0e50;
     log10QiMIN = 1.0e50;
     log10QdMAX = -1.0e50;
     log10QiMAX = -1.0e50;
% ------------------------------------------------------------------------
     log10Qmax = log10(maxQ);
     log10Qmin = log10(minQ);
     % x' = A*x + B
     % 0 = A*log10Qref) + B
     % log10Qmin = -2A + B
     %                 => A = -c/(1 + log10Qref)    &    B = log10Qmin + A
     Ad = -log10Qmin/(1 + log10QdThreshold);
     Bd = log10Qmin + Ad;
     Ai = -log10Qmin/(1 + log10QiThreshold);
     Bi = log10Qmin + Ai;
% -------------------------------------------- initialize some work arrays
     xNN = zeros(n0);
     xFF = zeros(n1);
     xNF = zeros(n0,n1);
     yNN = zeros(n0);
     yFF = zeros(n1);
     yNF = zeros(n0,n1);
     rNN = zeros(n0);
     rFF = zeros(n1);
     rNF = zeros(n0,n1);
     sNN = zeros(n0);
     sFF = zeros(n1);
     sNF = zeros(n0,n1);     
     otherwise
     error('unknown qType');
  end
%%                                                     prepare projections
% -------------------------- precalculate projected averages and variances
gotAve = zeros(nt,nModes); 
gotVar = zeros(nt,nModes);
ave0 = zeros(nt,nModes);
var0 = zeros(nt,nModes);
ave1 = zeros(nt,nModes);
var1 = zeros(nt,nModes);
%%                                 make projections on all labeled systems
   for j=1:nModes                    % consider each j-th vector-direction
   vec = U(:,j);                                     % the selected vector
      for kL=1:nL                                  % all unlabeled systems
      gotAve(kL,j) = vec'*Mgot{kL};                   % projected  average
      gotVar(kL,j) = max(vec'*Qgot{kL}*vec,add0);     % projected variance
      end
   end
%%                                         initialize counters of interest
% REMARK1:      upperSimilarity = max(fIndeterminable00,fIndeterminable11)
% REMARK2:      lowerSimilarity = min(fIndeterminable00,fIndeterminable11)
d_efficacyFuncton = zeros(1,nModes);                       % mode efficacy
i_efficacyFuncton = zeros(1,nModes);                       % mode efficacy
d_upperSimilarity = zeros(1,nModes); % greatest similarity between 00 & 11
i_upperSimilarity = zeros(1,nModes); % greatest similarity between 00 & 11
d_lowerSimilarity = zeros(1,nModes); % lowest similarity between 00 and 11
i_lowerSimilarity = zeros(1,nModes); % lowest similarity between 00 and 11
d_qualityProperty = zeros(1,nModes);                        % mode quality
i_qualityProperty = zeros(1,nModes);                        % mode quality
d_lnScoreFunctn10 = zeros(1,nModes);       % => mean ln of selection power
i_lnScoreFunctn10 = zeros(1,nModes);       % => mean ln of selection power
d_VoteConsensus10 = zeros(1,nModes); %=> mean discriminant consensus power
i_VoteConsensus10 = zeros(1,nModes); %=> mean indifference consensus power
% -------------------------------------- initialize class pair data arrays
DorIindicator = zeros(1,nCpairs);               % => (1, -1) => X = (D, I)
efficacyVectX = zeros(1,nCpairs);                    % X-subspace efficacy
qualityLevelX = zeros(1,nCpairs);              % X-subspace quality factor
voteFractionX = zeros(1,nCpairs);               % X-subspace vote fraction
lnScoreFunctX = zeros(1,nCpairs);             % X-subspace lnScoreFunction
upperSimilarX = zeros(1,nCpairs); %=> X-subspace: max(similar11,similar00)
lowerSimilarX = zeros(1,nCpairs); %=> X-subspace: min(similar11,similar00)
%%                         calculate properties of interest per projection
               for j=1:nModes                  % j-th requested projection
%%                            apply normal SPLOC for all pairs of clusters
% ************************************************************* start hack
% REMARK: The SPLOC code will be reused in this section. The difference is
% that the scores for each comparion made between a pair of clusters are
% added together to form a grand total score. For a given pair of clusters
% the dummy labels are still called 1 and 0, since the underlying process
% continues to be a binary comparison. The nonfunctional and functional
% terms are also used in places, but recall SPLOC is symmetrical w.r.t.
% these labels. Only one comparison is performed in SPLOC because there is
% only one pair (0 to 1). Here, nC*(nC - 1)/2 pairwise comparisons are 
% performed and the results are tallied.
% ------------------------------------------------------------------------
         cPair = 0;
            for cL0=1:nC    % start double-loop over all pairs of clusters
            n0 = nSystemsPerClass(cL0);
            ave0(1:n0,j) = gotAve( classSystemList{cL0} , j);
            var0(1:n0,j) = gotVar( classSystemList{cL0} , j);
               for cL1=cL0+1:nC
               cPair = cPair + 1;
               n1 = nSystemsPerClass(cL1);
               ave1(1:n1,j) = gotAve( classSystemList{cL1} , j);
               var1(1:n1,j) = gotVar( classSystemList{cL1} , j);
               totPairs00 = n0*(n0-1)/2;
               totPairs11 = n1*(n1 - 1)/2; 
               totPairs = n0*n1;
% ----------------------------------------------------- for debugging only
%                strL = ['cL0= ',num2str(cL0),' cL1= ',num2str(cL1)];
%                disp( dividerLine(strL) );
%                disp( [totPairs00, totPairs11, totPairs] );
%                disp( [n0,n1] );
%                aa0 = ave0(1:n0,j)';
%                vv0 = var0(1:n0,j)';
%                aa1 = ave1(1:n1,j)';
%                vv1 = var1(1:n1,j)';
%                disp( num2str( [classSystemList{cL0}; aa0; vv0] ) );
%                disp('  ');
%                disp( num2str( [classSystemList{cL1}; aa1; vv1] ) );
%                disp('  ');
%                pause 
%%                                                 SPLOC each cluster pair
% --------------------------------- get <ln[scoringFunction]> for 10-pairs
               lnScoreFunctn10 = 0;    % => get mean ln of selection power
               dVoteFraction10 = 0;  % => get discriminant consensus power
               iVoteFraction10 = 0;  % => get indifference consensus power
               fIndeterminable = 0;        % => quantify the indeterminacy
              %------------------------------------ monitoring consistency
               NconsistentSBN0 = 0;
               NconsistentSBN1 = 0;
               NconsistentVAR0 = 0;
               NconsistentVAR1 = 0;
               Nindeterminable = 0;
                 for k1=1:n1                  % functional 1-state systems
                 av1 = ave1(k1,j);
                 vQ1 = var1(k1,j);
                   for k0=1:n0             % nonfunctional 0-state systems
                   av0 = ave0(k0,j); 
                   vQ0 = var0(k0,j);
% -------------------------------------------------------- std. dev. ratio
                   minV = min(vQ1,vQ0);
                   maxV = max(vQ1,vQ0);
                   r1 = sqrt(vQ0/vQ1);
                   r2 = 1/r1;
                   r = sqrt(maxV/minV);     % <= symmetrical w.r.t. labels
                   noise = sqrt(vQ1 + vQ0); 
                   signal = abs(av1 - av0);
% ------------------------------------------- scoring function calculation
                   snr = signal/noise;             % signal to noise ratio
                   sbn = max(snr - 1,0);             % signal beyond noise
% -------------------------------------------- monitor indeterminate cases
                     if( ( snr < 1 ) && ...
                         ( r < minScore1 ) )    % => no significant signal
                     fff = 1 - snr*(r - 1)/(minScore1 - 1);  % tiny signal
                     fIndeterminable = fIndeterminable + fff;
                     end
% -------------------------------- begin monitor statistical consistencies
                     if( sbn > 0.667 )         % count only robust signals
                        if( av1 > av0 )
                        NconsistentSBN0 = NconsistentSBN0 + 1;  
                        else
                        NconsistentSBN1 = NconsistentSBN1 + 1;
                        end
                     else    % => weak signal => no inconsistency detected
                     NconsistentSBN0 = NconsistentSBN0 + 1;
                     NconsistentSBN1 = NconsistentSBN1 + 1;
                     end
                     if( r1 > maxScore0 )                   % => vQ0 > vQ1
                     NconsistentVAR0 = NconsistentVAR0 + 1;
                     end
                     if( r2 > maxScore0 )                   % => vQ1 > vQ0
                     NconsistentVAR1 = NconsistentVAR1 + 1;
                     end
                     if( (sbn < 0.01) && (r < minScore1) )
                     Nindeterminable = Nindeterminable + 1;
                     end
% ---------------------------------- end monitor statistical consistencies
% --------------------------------------- calculate ln of scoring function
                   r = r - 1;    % must do for scoring function definition
                   r = r/sqrt(1 + (r/99)^2);       % apply limiting factor
                   lnScore1 = log(sqrt(sbn*sbn+r*r) + 1);  % add 1 back in
                   %                   ^^^^^^^---> harder for discriminant
                   % --> Note: sbn on its own is a viable scoring function
                   lnScore0 = log(sqrt(snr*snr+r*r) + 1);  % add 1 back in
                   %                   ^^^^^^^---> harder for indifference
                   % ---> Note: pick harder case among sbn and snr options
                     if( lnScore0 < lnMidScoreX )  % using harder criteria
                     x = lnScore0;         % => satisfies harder condition
                     elseif(lnScore1>lnMidScoreX)  % using harder criteria
                     x = lnScore1;         % => satisfies harder condition
                     else       % split case?: set lnScore as undetermined
                     x = lnMidScoreX; %set score at center of undetermined
                     end              
% --------------------------------------------------- assign voter weights
                     if( x > xMax )
                     d_wt = 1;
                     i_wt = 0;
                     else
                     ii = 1 + floor(x/0.001);
                     d_wt = dWt(ii);
                     i_wt = iWt(ii);
                     end                                     
% -------------------------------------------------------- tally all pairs
                   lnScoreFunctn10 = lnScoreFunctn10 + x;
                   dVoteFraction10 = dVoteFraction10 + d_wt;
                   iVoteFraction10 = iVoteFraction10 + i_wt;
                   end      
                 end                 
% ---------------------------- selection & consensus power and consistency
               lnScoreFunctn10 = lnScoreFunctn10/totPairs;  % => selection
               dVoteFraction10 = dVoteFraction10/totPairs;  % => consensus
               iVoteFraction10 = iVoteFraction10/totPairs;  % => consensus
% ------------------------- account for statistical consistency conditions
% Some heuristic measures related to fidelity in the consistency across
% all members of each class yielding similar answers. When consistency is
% lost penalties are incorporated into the measures. These factors get 
% built into the efficacy score downstream, and as such, when consistency
% is poor, more undetermined modes is likely to be uncovered. 
% ------------------------------------------------------ for checking only
%    disp(['mode: ',num2str(j),'  Nsbn0= ',num2str(NconsistentSBN0), ...
%          '  Nsbn1= ',num2str(NconsistentSBN1), ...
%          '  Nvar0= ',num2str(NconsistentVAR0), ...
%          '  Nvar1= ',num2str(NconsistentVAR1), ...
%          '  fff= ',num2str(1 - Nindeterminable/totPairs)]);
               Nconsistent = max( [NconsistentSBN0, NconsistentSBN1, ...
                                   NconsistentVAR0, NconsistentVAR1] );
               ff1 = 1 - Nindeterminable/totPairs;          % => downgrade
               ff2 = Nconsistent/totPairs;                  % => downgrade
               consensus = 0.5 + 2.5*(ff1*ff2)^2;            % range [0,1]
               fIndeterminable = fIndeterminable/totPairs; %=> consistency
% REMARK:      ^^^^^^^^^^^^^^^---> quantifies the scoring function but not
%              the clustering properties. See next section for clustering!
%%                                         calculate clustering properties
   stdv0 = sqrt( var0(:,j) );          % nonfunctional: standard deviation
   stdv1 = sqrt( var1(:,j) );             % functional: standard deviation
% ------------------------------------------------------------------ logic
   switch qType
     case 'hardLinear'            % O( log(n0)*n0 + log(n1)*n1 ) algorithm
     doHardLinear = true;
     doSoftGeometry = false;
     case 'softGeometry'
     doHardLinear = false;
     doSoftGeometry = true;
     case 'hybrid'
     doHardLinear = true;
     doSoftGeometry = true;
   end
% ---------------------------------------------- proceed with calculations
     if( doHardLinear )
% ---------------------------------------------------- nonfunctional cases
               aveMean0 = mean( ave0(:,j) );
               minMean0 =  min( ave0(:,j) );
               maxMean0 =  max( ave0(:,j) );
              %--------------------------------------
               tmpLnStDv0 = log(stdv0);        % log of standard deviation
               aveLnStDv0 = mean(tmpLnStDv0); 
               minLnStDv0 =  min(tmpLnStDv0);
               maxLnStDv0 =  max(tmpLnStDv0);
% ------------------------------------------------------- functional cases
               aveMean1 = mean( ave1(:,j) ); 
               minMean1 =  min( ave1(:,j) );
               maxMean1 =  max( ave1(:,j) );
              %--------------------------------------
               tmpLnStDv1 = log(stdv1);        % log of standard deviation
               aveLnStDv1 = mean(tmpLnStDv1); 
               minLnStDv1 =  min(tmpLnStDv1);
               maxLnStDv1 =  max(tmpLnStDv1);
% -------------------------------------------------------- modify extremes
                 if( lnScoreFunctn10 > lnMidScoreX )         % d-subspace
                 aveStd0 = mean(stdv0);
                 tadShift0 = aveStd0/sqrt(Ns0);   % =0 when Ns0 -> infinty
                 tadDeflate0 = 1 - 0.5/sqrt(Ns0); % =1 when Ns0 -> infinty
                 tadInflate0 = 1/tadDeflate0;
                 % --------------------------------
                 aveStd1 = mean(stdv1);
                 tadShift1 = aveStd1/sqrt(Ns1);   % =0 when Ns1 -> infinty
                 tadDeflate1 = 1 - 0.5/sqrt(Ns1); % =1 when Ns1 -> infinty
                 tadInflate1 = 1/tadDeflate1;
                 minMean0 = minMean0 - tadShift0;
                 maxMean0 = maxMean0 + tadShift0;
                 minMean1 = minMean1 - tadShift1;
                 maxMean1 = maxMean1 + tadShift1;
                 minLnStDv0 = minLnStDv0 + log(tadDeflate0);
                 maxLnStDv0 = maxLnStDv0 + log(tadInflate0);
                 minLnStDv1 = minLnStDv1 + log(tadDeflate1);
                 maxLnStDv1 = maxLnStDv1 + log(tadInflate1);
                 end
% ------------------------------------------- quantify clusters on 1D line
               rangeMean = max(maxMean0,maxMean1) ...  % 111       0000000
                         - min(minMean0,minMean1);     % <---- range ---->
               %                               range = span1 + gap + span0
               %
               %                                           <--- span0 --->
               spanMean0 = maxMean0 - minMean0;  % e.g. => 000000000000000
               spanMean1 = maxMean1 - minMean1;      % e.g. => 11111111111
               %                                               <- span1 ->
               %
               gap10 = minMean0 - maxMean1;    % e.g. => 11111 -gap- 00000
               gap01 = minMean1 - maxMean0;    % e.g. => 00000 -gap- 11111
               gapMean = max( [0,gap10,gap01] );
                 if( gapMean > 0 )                   %=> gap>0 & overlap=0
                 overlapMean = 0;
                 else                                %=> gap=0 & overlap>0
                 overlapList = [spanMean0,spanMean1, -max(gap10,gap01)];
                 overlapMean = min(overlapList);   
                 % Notes on overlap. example 1:     1110101100100110010000
                 %                                     <-- overlap -->
                 %                                  <------  range ------>
                 %
                 %                    <----- not overlap ------>
                 %                       <------ not overlap -------->
                 %         example 2: 00010101000100110010111011000000
                 %                       <------- overlap ------>
                 %                    <----------- range ------------>
                 %
                 %      =>  overlap = min([span0,span1,-max(gap10,gap01)])
                 end               
% ------------ see comments with pics above for range, span, gap & overlap
               rangeLnStDv = max(maxLnStDv0,maxLnStDv1) ...
                           - min(minLnStDv0,minLnStDv1);
               spanLnStDv0 = maxLnStDv0 - minLnStDv0;
               spanLnStDv1 = maxLnStDv1 - minLnStDv1;
               gap10 = minLnStDv0 - maxLnStDv1; 
               gap01 = minLnStDv1 - maxLnStDv0;
               gapLnStDv = max( [0,gap10,gap01] );
                 if( gapLnStDv > 0 )                 %=> gap>0 & overlap=0
                 overlapLnStDv = 0;
                 else                                %=> gap=0 & overlap>0
                 overlapList = [spanLnStDv0,spanLnStDv1, ...
                                -max(gap10,gap01)];
                 overlapLnStDv = min(overlapList); 
                 end                  
% ------------------------------------------- construct universal measures
               sig0 = exp(aveLnStDv0);   % => geometric mean for std. dev.
               sig1 = exp(aveLnStDv1);   % => geometric mean for std. dev.
               sigMeanMin = min(sig0,sig1);
               sigMeanMax = max(sig0,sig1);
               r = sigMeanMax/sigMeanMin;
               d = sigMeanMin*sqrt(1+r*r);  
               % --------------------------------------- 4 scaled measures
               rangeMean = max(1.0e-9,rangeMean/d);
               snr10Mean = abs( aveMean1 - aveMean0)/d;
               overlapMean = overlapMean/d;
               gapMean = gapMean/d;        %=> d sets 10-pair length scale
% ------------------------------------------------------------------------
               d = lnMinScore1;          % meaningful on an absolute scale
               %d = lnMaxScore0;         % meaningful on an absolute scale
               %d = lnMidScoreX;                    % split the difference
               % Seems all three alternatives work. The larger the value,
               % the more difficult it is to find discriminant solutions
               % but the clustering separation should be increased due to
               % being picky. But it is hard to notice differences.
               % --------------------------------------- 4 scaled measures
               rangeLnStDv = max(1.0e-9,rangeLnStDv/d);
               snr10LnStDv = abs( aveLnStDv1 - aveLnStDv0)/d;
               overlapLnStDv = overlapLnStDv/d;
               gapLnStDv = gapLnStDv/d;    %=> d set absolute length scale
% ----------------------------------------------- scale invariant measures
% REMARK: sig0 & sig1 set the 00-pair & 11-pair length scales respectively
               spanMean0 = max(sig0*0.000000001,spanMean0); 
               spanMean1 = max(sig1*0.000000001,spanMean1); 
               spanLnStDv0 =max(d*0.000000001,spanLnStDv0);
               spanLnStDv1 =max(d*0.000000001,spanLnStDv1); 
               % ---------------------------------------------------------
               rSpan = max(spanMean0,spanMean1)/min(spanMean0,spanMean1);
               rGap = 1 + gapMean/rangeMean;
               rOverlap = 1 + overlapMean/rangeMean;
               rMobius = (snr10Mean - 1)/(snr10Mean + 1); 
               cpMUg = gapMean*rGap*rSpan + rMobius;            % ^ Note 1
               cpMUo = overlapMean*rOverlap/rSpan - rMobius;    % ^ Note 2
               % ---------------------------------------------------------
               rSpan = max(spanLnStDv0,spanLnStDv1)/ ...
                       min(spanLnStDv0,spanLnStDv1);
               rGap = 1 + gapLnStDv/rangeLnStDv;
               rOverlap = 1 + overlapLnStDv/rangeLnStDv;
               rMobius = (snr10LnStDv - 1)/(snr10LnStDv + 1); 
               cpSDg = gapLnStDv*rGap*rSpan + rMobius;          % ^ Note 3
               cpSDo = overlapLnStDv*rOverlap/rSpan - rMobius;  % ^ Note 4
               rFeature = sqrt(gapMean^2 + gapLnStDv^2);
% ---------------------------------------------------------------- summary
% ---------------------------------------------- invoked cluster property?
%   MEAN                    STD. DEV.            d-subspace    i-subspace
% rangeMean                rangeLnStDv             yes           yes
% snr10Mean                snr10LnStDv             yes           yes
% gapMean                  gapLnStDv               yes           yes
% overlapMean              overlapLnStDv           yes           yes
% spanMean0                spanLnStDv0             yes           yes
% spanMean1                spanLnStDv1             yes           yes
% ---------------------------------------------
% vFrac                                            yes           yes
% fIndeterminable                                  yes           yes
% fIndeterminable00                                yes           yes
% fIndeterminable11                                yes           yes
% Notes:                                        cluster property functions
% 1) cpMUg = nonlinear function. When >> 1 =>     means separate very well
% 2) cpMUo = nonlinear function. When >> 1 =>    means virtually identical
%            Note: When cpMUg > 0, cpMUo = 0, and vice versa
% 3) cpSDg = nonlinear function. When >> 1 => std. dev. separate very well
% 4) cpSDo = nonlinear function. When >> 1 => std. dev. virtualy identical
%            Note: When cpSDg > 0, cpSDo = 0, and vice versa 
% ---------------------------------------------------------- calculate dcp
               dcp1 = cpMUg - cpMUo;
               dcp2 = cpSDg - cpSDo;
               dcp3 = dcp1 + dcp2;
               dcp = max( [dcp1,dcp2,dcp3] );   bb    % either one or both
% ---------------------------------------------------------- calculate icp
               icp1 = cpMUo - cpMUg;                                % MEAN
               icp2 = cpSDo - cpSDg;                           % std. dev.
               icp3 = icp1 + icp2;
               icp = min( [icp1,icp2,icp3] );    % require icp1 & icp2 > 0
% using an if-condition to split the 2 cases does not alter the calc-speed
% ------------------------------------------------------------------------
     else                             % no hardLinear results will be used
     dcp = -0.0000000001;                                  % use as a flag
     icp = -0.0000000001;                                  % use as a flag
     rFeature = 0;                                         % use as a flag
     end
     if( doSoftGeometry && (rFeature < rFeatureSoft) )
% ----------------------------------- apply new clustering quality measure
% Created by Dr. Jacobs in June 2021 in relation to analyzing EEG signals.
% The hard separation criteria is replaced by a soft criteria for cluster
% quality. This will allow some u-modes to be recognized as d-modes. It 
% will also be able to start with something not good, but with potential,
% to self-grow and become strong clustering. To understand the rationale
% for the geometry-based measure, see testClusterMetric_v4.m
% ------------------------------------------- normalize data onto ~[0,1]^2
               xxF = ave1(:,j)';
               xxN = ave0(:,j)';
               tAdd = [xxF,xxN];
               maxX = max(tAdd);
               minX = min(tAdd);
               xSpan = max(maxX - minX, 1.0e-10);
               xWidth = 1.2*xSpan;
               xxF = 0.1 + (xxF - minX)/xWidth;
               xxN = 0.1 + (xxN - minX)/xWidth;
% -----------------------------------------------
               yyF = stdv1';
               yyN = stdv0';
               tAdd = [yyF,yyN];
               maxY = max(tAdd);
               minY = min(tAdd);
               ySpan = max(maxY - minY, 1.0e-10);
               yWidth = 1.2*ySpan;
               yyF = 0.1 + (yyF - minY)/yWidth;
               yyN = 0.1 + (yyN - minY)/yWidth;
% ================================================ calculate quality level
                 for kNa=1:n0
                    for kNb=1:n0
                    axx = abs( xxN(kNa) - xxN(kNb) );
                    ayy = abs( yyN(kNa) - yyN(kNb) );
                    arr = sqrt( axx*axx + ayy*ayy );
                    xNN(kNa,kNb) = max(axx,haloSame);
                    yNN(kNa,kNb) = max(ayy,haloSame);
                    rNN(kNa,kNb) = max(arr,haloSame);
                    sNN(kNa,kNb) = max(arr,haloDiff);
                    end
                 xNN(kNa,kNa) = 1.0e50;
                 yNN(kNa,kNa) = 1.0e50;
                 rNN(kNa,kNa) = 1.0e50;
                 end
% ---------------------------------------------------------------
                 for kFa=1:n1
                    for kFb=1:n1
                    axx = abs( xxF(kFa) - xxF(kFb) );
                    ayy = abs( yyF(kFa) - yyF(kFb) );
                    arr = sqrt( axx*axx + ayy*ayy );
                    xFF(kFa,kFb) = max(axx,haloSame);
                    yFF(kFa,kFb) = max(ayy,haloSame);
                    rFF(kFa,kFb) = max(arr,haloSame);
                    sFF(kFa,kFb) = max(arr,haloDiff);
                    end
                 xFF(kFa,kFa) = 1.0e50;
                 yFF(kFa,kFa) = 1.0e50;
                 rFF(kFa,kFa) = 1.0e50;
                 end
% ---------------------------------------------------------------
                 for kFa=1:n1
                    for kNb=1:n0
                    axx = abs( xxF(kFa) - xxN(kNb) );
                    ayy = abs( yyF(kFa) - yyN(kNb) );
                    arr = sqrt( axx*axx + ayy*ayy );
                    xNF(kNb,kFa) = max( axx - haloDiff , 0);
                    yNF(kNb,kFa) = max( ayy - haloDiff , 0);
                    rNF(kNb,kFa) = max( arr - haloDiff , 0);
                    sNF(kNb,kFa) = max(arr,haloDiff);
                    end
                 end 
% ---------------------------------------------------------------
               dNNmin = min(xNN);
               dFFmin = min(xFF);
               dNFmin = min(xNF);
               dFNmin = min(xNF');
               log10qNx = 2*dFNmin./(dNNmin + dFNmin) - 1;
               log10qFx = 2*dNFmin./(dFFmin + dNFmin) - 1;
               log10Qx = sum([log10qNx,log10qFx])/(n1 + n0); 
% ---------------------------------------------------------------
               dNNmin = min(yNN);
               dFFmin = min(yFF);
               dNFmin = min(yNF);
               dFNmin = min(yNF');
               log10qNy = 2*dFNmin./(dNNmin + dFNmin) - 1;
               log10qFy = 2*dNFmin./(dFFmin + dNFmin) - 1;
               log10Qy = sum([log10qNy,log10qFy])/(n1 + n0);  
% ---------------------------------------------------------------
               dNNmin = min(rNN);
               dFFmin = min(rFF);
               dNFmin = min(rNF);
               dFNmin = min(rNF');
               log10qNr = 2*dFNmin./(dNNmin + dFNmin) - 1;
               log10qFr = 2*dNFmin./(dFFmin + dNFmin) - 1;
               log10Qr = sum([log10qNr,log10qFr])/(n1 + n0);       
% ------------------------------------------------------------------------
               log10Qd = max([log10Qx,log10Qy,log10Qr]);   % => for d-mode
               log10QdMIN = min(log10QdMIN,log10Qd);
               log10QdMAX = max(log10QdMAX,log10Qd);
               log10Qd = Ad*log10Qd + Bd;
               log10Qd = log10Qd*(1 + 10*log10Qd^2);
                  if( log10Qd > 0 )
                  log10Qd = log10Qd/sqrt(1 + (log10Qd/log10Qmax)^2 );
                  end
               Qd = 10^(log10Qd);                              % raw value
% ------------------------------------------------------------------------
               dNNmin = min(sNN);
               dFFmin = min(sFF);
               dNFmin = min(sNF);
               dFNmin = min(sNF');
               log10qFs = 4*min([dFFmin; dNFmin])./(dFFmin + dNFmin) - 1;
               log10qNs = 4*min([dNNmin; dFNmin])./(dNNmin + dFNmin) - 1;
               log10qF = mean(log10qFs);
               log10qN = mean(log10qNs);
               log10Qi = min(log10qF,log10qN);
               log10QiMIN = min(log10QiMIN,log10Qi);
               log10QiMAX = max(log10QiMAX,log10Qi);
               log10Qi = Ai*log10Qi + Bi;
               log10Qi = log10Qi*(1 + 10*log10Qi^2);
                  if( log10Qi > 0 )
                  log10Qi = log10Qi/sqrt(1 + (log10Qi/log10Qmax)^2 );
                  end
               Qi = 10^(log10Qi);                              % raw value
% ============================================ finalize clustering quality
               if( lnScoreFunctn10 > lnMidScoreX )  % => possibly a d-mode
                  if( (Qd - Qi) > 1 )                        % => a d-mode
                  Q = Qd;
                  else                                         % => u-mode
                  dQ = (Qd - Qi);
                     if( dQ < -1 )
                     Q = -Qi;                         % => maximum penalty
                     else
                     fd = 0.5*(dQ + 1);
                     fi = 1 - fd;
                     Q = fd*Qd - fi*Qi;          % include gradual penalty
                     end
                  end
               dcp2 = Q;  
                  if( (dcp < 0) && (dcp2 < 0) )
                  dcp = min(dcp2,dcp);
                  else
                  dcp = max(dcp2,dcp);
                  end
               else                                % => possibly an i-mode
                  if( (Qi - Qd) > 1 )                       % => an i-mode
                  Q = Qi;
                  else                                         % => u-mode
                  dQ = (Qi - Qd);
                     if( dQ < -1 )
                     Q = -Qd;                         % => maximum penalty
                     else
                     fi = 0.5*(dQ + 1);
                     fd = 1 - fi;
                     Q = fi*Qi - fd*Qd;          % include gradual penalty
                     end
                  end
               icp2 = Q; 
                 if( (icp < 0) && (icp2 < 0) )
                 icp = min(icp2,icp);
                 else
                 icp = max(icp2,icp);
                 end
               end
% ------------------------------------------------------------------------
     rFeature = rFeatureSoft;
     end
%%                           calculate modifiers for clustering properties
% ----------------------------------------------------- consistency checks
               fIndeterminable00 = 0;
                 for k0=1:n0
                 av0 = ave0(k0,j);
                 sd0 = stdv0(k0);
                  % ------------------------------------------------------
                   for k1=k0+1:n0                 % using k1 as a dummy k0
                   av1 = ave0(k1,j);
                   sd1 = stdv0(k1);
                   signal = abs(av1 - av0);
                   % ------------------------------------- std. dev. ratio
                   sigMin = min(sd0,sd1);
                   sigMax = max(sd0,sd1);
                   r = sigMax/sigMin;       % <= symmetrical w.r.t. labels
                   noise = sqrt(sigMin*sigMin + sigMax*sigMax);
                   snr = signal/noise;             % signal to noise ratio
% -------------------------------------------- monitor indeterminate cases
                     if( (snr < 1) && ...
                         (r < minScore1) )      % => no significant signal
                     fff = 1 - snr*(r - 1)/(minScore1 - 1);  % tiny signal
                     fIndeterminable00 = fIndeterminable00 + fff;
                     end                         
                   end
                 end
               % ---------------------------------------------------------
               fIndeterminable11 = 0;
                 for k1=1:n1
                 av1 = ave1(k1,j);
                 sd1 = stdv1(k1);
                   for k0=k1+1:n1                 % using k0 as a dummy k1
                   av0 = ave1(k0,j);
                   sd0 = stdv1(k0);
                   signal = abs(av1 - av0);
                   % ------------------------------------- std. dev. ratio
                   sigMin = min(sd0,sd1);
                   sigMax = max(sd0,sd1);
                   r = sigMax/sigMin;       % <= symmetrical w.r.t. labels
                   noise = sqrt(sigMin*sigMin + sigMax*sigMax);  
                   snr = signal/noise;             % signal to noise ratio
% -------------------------------------------- monitor indeterminate cases
                     if( ( snr < 1 ) && ...
                         ( r < minScore1 ) )    % => no significant signal
                     fff = 1 - snr*(r - 1)/(minScore1 - 1);  % tiny signal
                     fIndeterminable11 = fIndeterminable11 + fff;
                     end
                   end
                 end 
               fIndeterminable00 = fIndeterminable00/totPairs00;
               fIndeterminable11 = fIndeterminable11/totPairs11;
%%                                          congruent subspace bifurcation
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& start of bifurcation of subspace selection
               if( lnScoreFunctn10 > lnMidScoreX ) % discriminant subspace
% +++++++++++++++++++++++++++++++++++++++++++++++ calculate dQualityFactor
               vFraction = dVoteFraction10;
                 %if( pursuitType ~= -1 )     % => requires dQualityFactor
                 %disp( [j,gapMean,gapLnStDv] );
                 %rFeature = sqrt(gapMean^2 + gapLnStDv^2);
                 vFrac = vFraction - voteThreshold;
                    if( vFrac < 0 )                % => poor quality level
                    vFrac = 0;
                    else                           % => good quality level
                    vFrac = vFrac/(1 - voteThreshold);  % normalize +range
                    end
                 fImax = max(fIndeterminable00,fIndeterminable11);
                %^^^^^-------------------------- only one needs to be good
                    if( (fImax > 0.5) && (dcp > 0) )   % qualify for bonus
                    fImin = min(fIndeterminable00,fIndeterminable11);
                    %^^^^-----> the bigger fImin is the better as a bonus!
                    tmp = 0.1*(0.1 + 1.9*vFrac)* ...
                          max(0,0.5 - fIndeterminable); 
                    dAdd = tmp/max(0.01,1 - fImax) ...
                         + tmp/max(0.01,1 - fImin) ...
                         + 2*(fImax - 0.5) ...
                         + max(0, (4*fImax-3)/(3 - 2*fImax) ) ...  % [0,1]
                         + max(0, (4*fImin-3)/(3 - 2*fImin) ) ...  % [0,1]
                         + max(0, (1 - 5*fIndeterminable)/ ...       
                           (1 + 3*fIndeterminable) ) ...    % range: [0,1]
                         + vFrac*consensus ...
                         + vFrac;
                    else  % bonus means there is only gain without penalty
                    dAdd = 0;                          % => nothing to add
                    end
% ------------------------------------------------ positive/negative cases
                   if( dcp > 0 )                % => forward learning bias
                   required1 = max(-0.99, 8*(fImax - 0.5)^3 );
                   required2 = max(-0.99, 8*(0.5 - fIndeterminable)^3 );
                   required = min(required1,required2);
                   boost = required + consensus + consensus*vFraction;
                   %disp([dcp,boost,dAdd,rFeature,fImax])
                   dcp = dcp*vFraction*boost/(1 + boost);
                   dcp = dcp*fImax*( (1 - (1-fImax)^8 )^6 );
                   dcp = dcp*dAdd/(1 + dAdd);
                   dcp = dcp*rFeature/rFeatureSoft;
                     if( dcp < 1 )
                     dcp = dcp^6;
                     end
                   %disp(dcp);
                   %pause
                   dQualityFactor = dampLR*dcp/sqrt(1 + (dcp/10)^2);
                   else      % => reverse bias w.r.t. minimum forward bias
                   dQualityFactor = dcp;
                   end 
                 % -------------------------------------------------------
                 % REMARK: Although vFraction has been included in the
                 % calculation of xx dQuality tends to reach its maximum
                 % of 10 relatively easily. By including the factor
                 % that is a function of "v", and the vFraction again, 
                 % dQualityFactor has some variation, where a value of
                 % 10 indicating v = 1 and vFraction = 1. 
                 % Note that "v" is a nonlinear term both related to 
                 % clustering. The two classes need to be separated in 
                 % at least one feature, and at least of the two classes
                 % must have good pooling. Both should be good, leading
                 % to the multiplication as an "and" operation. Also, a
                 % maximum of 1 is placed on this factor since the 
                 % non-linear function maps v on [0,1] to [0,1]. 
                 %--------------------------------------------------------
                 %else            % => looking for indifference modes only
                 %dQualityFactor = -0.1;    % <= default anti-discriminant
                 %end
               span = (lnScoreFunctn10 - lnMidScoreX);
               sWeight = sqrt(span) + span*(1 + span);
               %         ^^^^^^^^^^-------> rapidly bifurcates from span=0
                  if( dQualityFactor > 0 )
                      if( dQualityFactor > minQF )
                      vecScore = dQualityFactor*sWeight;
                      else                          % tends to nucleate by
                      vecScore = minQF*sWeight;  % encouraging bifurcation
                      end      % ^^^^^------------------> boosts the score
                  else
                  vecScore = dQualityFactor*sWeight;
                  end
               vecScore = vecScore*fDm;
               DorIindicator(cPair) = 1;                        % => X = D
               efficacyVectX(cPair) = vecScore;      % D-subspace efficacy
               qualityLevelX(cPair) = dQualityFactor;
               voteFractionX(cPair) = vFraction;
               lnScoreFunctX(cPair) = lnScoreFunctn10; 
               upperSimilarD = max(fIndeterminable00,fIndeterminable11);
               lowerSimilarD = min(fIndeterminable00,fIndeterminable11);
               upperSimilarX(cPair) = upperSimilarD;
               lowerSimilarX(cPair) = lowerSimilarD;
               else % -----------------------------> indifference subspace
% +++++++++++++++++++++++++++++++++++++++++++++++ calculate iQualityFactor
               vFraction = iVoteFraction10;
                 %if( pursuitType ~= 1 )      % => requires iQualityFactor
                 vFrac = vFraction - voteThreshold;
                    if( vFrac < 0 )                % => poor quality level
                    vFrac = 0;
                    else                           % => good quality level
                    vFrac = vFrac/(1 - voteThreshold);  % normalize +range
                    end
                 fImin = min( [fIndeterminable00, fIndeterminable, ...
                               fIndeterminable11] );
                 %^^^^^----------------> all three conditions are required
                    if( (fImin > 0.5) && (icp > 0) )   % qualify for bonus
                    iAdd = (0.1 + 1.9*vFraction)* ...
                            max(0,fIndeterminable - 0.5)/ ...
                            max(0.01,1 - fImin) ... 
                         + 2*(fImin - 0.5) ...
                         + vFrac*consensus ...
                         + vFrac;
                    else  % bonus means there is only gain without penalty
                    iAdd = 0;                          % => nothing to add
                    end
% ------------------------------------------------ positive/negative cases
                   if( icp > 0 )                % => forward learning bias
                   boost = 1 + consensus + vFraction ...
                             + consensus*vFraction ...
                             + max(-0.45, 8*(fImin - 0.5)^3 ) ...
                             + max(-0.45, 8*(fIndeterminable - 0.5)^3);
                   boost = sqrt(boost);
                   %disp([icp,boost,iAdd])
                   icp = icp*(min(2,0.4*iAdd) + min(4,boost*0.4) + 0.5);
                   %disp(icp);
                   %pause
                   iQualityFactor = dampLR*min(icp,10);
                   else      % => reverse bias w.r.t. minimum forward bias
                   iQualityFactor = icp;
                   end 
                % --------------------------------------------------------
                 %else            % => looking for discriminant modes only
                 %iQualityFactor = -0.1;    % <= default anti-indifference
                 %end
               span1 = (lnMidScoreX - lnScoreFunctn10); 
               span2 = span1*span1;
               span = span1 + 6.4027*span2 + 28.0462*span2*span2;
              %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^--> Note 1
              % Note 1: This allows the throw of the indifference subspace
              % to be about the same as the throw of discriminant subspace
              % lnMidScoreX = ( log(2) + log(1.3) )/2 = 0.4778
              % maximum selection power = 100 for discriminant situations.
              % log(100)= 4.6052 = maximum selection power for d-subspace
              % maximum throw for d-subspace = 4.6052 - 0.4778 = 4.1274
              % maximum throw for i-subspace = 0.4778 - log(1) = 0.4778
              % set maximum selection power for i-subspace = log(30)
              % Not as strong as d-modes, but has an appreciable impact.
              % After non-linear transformation i-subspace throw = 3.4012
              % Hence: span is an "effective span" as a function of span1
              % The function is approximately = to span1 for small span1.
              % ----------------------------------------------------------
               sWeight = sqrt(span) + span*(1 + span);
              %          ^^^^^^^^^^-------> rapidly bifurcates from span=0
                  if( iQualityFactor > 0 )
                      if( iQualityFactor > minQF )
                      vecScore = iQualityFactor*sWeight;
                      else                          % tends to nucleate by
                      vecScore = minQF*sWeight;  % encouraging bifurcation
                      %          ^^^^^------------------> boosts the score
                      end
                  else
                  vecScore = iQualityFactor*sWeight;
                  end
               vecScore = sVS*(vecScore*fIm);     % => bias toward d-modes
              %      This ^^^-> factor is incorporated for these reasons:
              % Units on efficacy are arbitrary. Any positive constant can
              % multiply efficacy as a scale. When comparing d-modes to 
              % d-modes this scale factor is irrelevant. When comparing 
              % i-modes to i-modes this scale factor is irrelevant. When
              % comparing d-modes to i-modes, the d-modes will win in more
              % head-2-head vector pair spins since d-modes have sVSx more
              % weight relative to i-modes. With the objective for using
              % SPLOC is to discriminate two systems, a 10-fold bias is
              % used. This asymmetry may reduce # of multiple solutions.
              % ----------------------------------------------------------
               DorIindicator(cPair) = -1;                       % => X = I
               efficacyVectX(cPair) = vecScore;      % I-subspace efficacy
               qualityLevelX(cPair) = iQualityFactor;
               voteFractionX(cPair) = vFraction;
               lnScoreFunctX(cPair) = lnScoreFunctn10; 
               upperSimilarI = max(fIndeterminable00,fIndeterminable11); 
               lowerSimilarI = min(fIndeterminable00,fIndeterminable11); 
               upperSimilarX(cPair) = upperSimilarI;
               lowerSimilarX(cPair) = lowerSimilarI;
               end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& end of bifurcation of subspace selection
               end
            end
% -------------------------------- combine quantities from all class pairs
         d_efficacyVectD = 0;
         d_qualityLevelD = 0;
         d_voteFractionD = 0;
         d_lnScoreFunctD = 0;
         d_upperSimilarD = 0;
         d_lowerSimilarD = 0;
         d_nConditionalD = 0;
         d_pcwtD = 0;
        % ------------------------ get weakest link
         u_efficacyVectD = 10000;
         u_qualityLevelD = 10000;
         u_voteFractionD = 10000;
         u_lnScoreFunctD = 10000;
         u_upperSimilarD = 10000;
         u_lowerSimilarD = 10000;
         u_nConditionalD = 0;
        % ------------------------ get weakest link
         u_efficacyVectI = 10000;
         u_qualityLevelI = 10000;
         u_voteFractionI = 10000;
         u_lnScoreFunctI = 0;
         u_upperSimilarI = 10000;
         u_lowerSimilarI = 10000;
         u_nConditionalI = 0;
        %------------------
         i_efficacyVectI = 0;
         i_qualityLevelI = 0;
         i_voteFractionI = 0;
         i_lnScoreFunctI = 0;
         i_upperSimilarI = 0;
         i_lowerSimilarI = 0;
         i_nConditionalI = 0;
         i_pcwtI = 0;
            for cPair=1:nCpairs
                if( DorIindicator(cPair) > 0 )          % => favors d-mode
                efficacyVectD = efficacyVectX(cPair);
                lnScoreFunctD = lnScoreFunctX(cPair);
                voteFractionD = voteFractionX(cPair);
                qualityLevelD = qualityLevelX(cPair);
                upperSimilarD = upperSimilarX(cPair);
                lowerSimilarD = lowerSimilarX(cPair);
                p = pcwt(cPair);
                   if( (lnScoreFunctD > lnMinScore1) && ...
                       (voteFractionD > voteThreshold) && ...
                       (qualityLevelD > minQF) ) %=> discriminant subspace
                   d_pcwtD         = d_pcwtD         + p;
                   d_efficacyVectD = d_efficacyVectD + p*efficacyVectD;
                   d_qualityLevelD = d_qualityLevelD + p*qualityLevelD;
                   d_lnScoreFunctD = d_lnScoreFunctD + p*lnScoreFunctD;
                   d_voteFractionD = d_voteFractionD + p*voteFractionD;
                   d_upperSimilarD = d_upperSimilarD + p*upperSimilarD;
                   d_lowerSimilarD = d_lowerSimilarD + p*lowerSimilarD;
                   d_nConditionalD = d_nConditionalD + 1;
                   else                         % => undetermined subspace
                   u_efficacyVectD = min(u_efficacyVectD,efficacyVectD);
                   u_qualityLevelD = min(u_qualityLevelD,qualityLevelD);
                   u_lnScoreFunctD = min(u_lnScoreFunctD,lnScoreFunctD);
                   u_voteFractionD = min(u_voteFractionD,voteFractionD);
                   u_upperSimilarD = min(u_upperSimilarD,upperSimilarD);
                   u_lowerSimilarD = min(u_lowerSimilarD,lowerSimilarD);
                   u_nConditionalD = u_nConditionalD + 1;
                   end  
                else  % DorIindicator < 0                 => favors i-mode
                efficacyVectI = efficacyVectX(cPair);
                lnScoreFunctI = lnScoreFunctX(cPair);
                voteFractionI = voteFractionX(cPair);
                qualityLevelI = qualityLevelX(cPair);
                upperSimilarI = upperSimilarX(cPair);
                lowerSimilarI = lowerSimilarX(cPair);
                p = pcwt(cPair);
                   if( (lnScoreFunctI < lnMaxScore0) && ...
                       (voteFractionI > voteThreshold) && ...
                       (qualityLevelI > minQF) ) %=> indifference subspace
                   i_pcwtI         = i_pcwtI         + p;
                   i_efficacyVectI = i_efficacyVectI + p*efficacyVectI;
                   i_qualityLevelI = i_qualityLevelI + p*qualityLevelI;
                   i_lnScoreFunctI = i_lnScoreFunctI + p*lnScoreFunctI;
                   i_voteFractionI = i_voteFractionI + p*voteFractionI;
                   i_upperSimilarI = i_upperSimilarI + p*upperSimilarI;
                   i_lowerSimilarI = i_lowerSimilarI + p*lowerSimilarI;
                   i_nConditionalI = i_nConditionalI + 1;
                   else                         % => undetermined subspace
                   u_efficacyVectI = min(u_efficacyVectI,efficacyVectI);
                   u_qualityLevelI = min(u_qualityLevelI,qualityLevelI);
                   u_lnScoreFunctI = max(u_lnScoreFunctI,lnScoreFunctI);
                   u_voteFractionI = min(u_voteFractionI,voteFractionI);
                   u_upperSimilarI = min(u_upperSimilarI,upperSimilarI);
                   u_lowerSimilarI = min(u_lowerSimilarI,lowerSimilarI);                   
                   u_nConditionalI = u_nConditionalI + 1;
                   end
                end
            end
% ----------------------------------------- stratify all possible outcomes
            if( d_nConditionalD > 0 )
            d_efficacyVectD = d_efficacyVectD/d_pcwtD;
            d_qualityLevelD = d_qualityLevelD/d_pcwtD;
            d_lnScoreFunctD = d_lnScoreFunctD/d_pcwtD;
            d_voteFractionD = d_voteFractionD/d_pcwtD;
            d_upperSimilarD = d_upperSimilarD/d_pcwtD;
            d_lowerSimilarD = d_lowerSimilarD/d_pcwtD;
            %dLd = true;                              % for debugging only
            else                       % some quantities must reset from 0
            d_lnScoreFunctD = lnMidScoreX;
            %dLd = false;                             % for debugging only
            end
% ----------------------------------------
            if( u_nConditionalD > 0 )
            uLd = true;
            else                       % some quantities must reset from 0
            u_efficacyVectD = 0;
            u_qualityLevelD = 0;
            u_lnScoreFunctD = lnMidScoreX;
            u_voteFractionD = 0;
            u_upperSimilarD = 0;
            u_lowerSimilarD = 0;
            uLd = false;
            end
% ----------------------------------------
            if( i_nConditionalI > 0 )
            i_efficacyVectI = i_efficacyVectI/i_pcwtI;
            i_qualityLevelI = i_qualityLevelI/i_pcwtI;
            i_lnScoreFunctI = i_lnScoreFunctI/i_pcwtI;
            i_voteFractionI = i_voteFractionI/i_pcwtI;
            i_upperSimilarI = i_upperSimilarI/i_pcwtI;
            i_lowerSimilarI = i_lowerSimilarI/i_pcwtI;
            %iLi = true;                              % for debugging only
            else                       % some quantities must reset from 0
            i_lnScoreFunctI = lnMidScoreX; 
            %iLi = false;                             % for debugging only
            end
% ----------------------------------------
            if( u_nConditionalI > 0 )
            uLi = true;
            else                       % some quantities must reset from 0
            u_efficacyVectI = 0;
            u_qualityLevelI = 0;
            u_lnScoreFunctI = lnMidScoreX; 
            u_voteFractionI = 0;
            u_upperSimilarI = 0;
            u_lowerSimilarI = 0;
            uLi = false;
            end
% -------------------------------------- partition results into categories
% REMARK: Note that except for within the undetermined subspace, if there
%         is implication that the variables represent the result of the Ld
%         or Li or Ldi subspaces, these results can be null. In particular
%         the results need not be good. However, this is the proper way to
%         partition all condition averages for all possible outcomes.
          uL = or(uLd,uLi);
            if( uL )
            uBoth = and(uLd,uLi);
              if( uBoth )     % => Lu subspace
              efficacyVectD = u_efficacyVectD;
              lnScoreFunctD = u_lnScoreFunctD;
              voteFractionD = u_voteFractionD;
              qualityLevelD = u_qualityLevelD;
              upperSimilarD = u_upperSimilarD;
              lowerSimilarD = u_lowerSimilarD;
             % -------------------------------
              efficacyVectI = u_efficacyVectI;
              lnScoreFunctI = u_lnScoreFunctI;
              voteFractionI = u_voteFractionI;
              qualityLevelI = u_qualityLevelI;
              upperSimilarI = u_upperSimilarI;
              lowerSimilarI = u_lowerSimilarI;
              elseif( uLd )   % => Li subspace
              efficacyVectD = u_efficacyVectD;
              lnScoreFunctD = u_lnScoreFunctD;
              voteFractionD = u_voteFractionD;
              qualityLevelD = u_qualityLevelD;
              upperSimilarD = u_upperSimilarD;
              lowerSimilarD = u_lowerSimilarD;
              iPenalty = (1 - u_nConditionalD/nCpairs)^2;
             % -------------------------------
              efficacyVectI = i_efficacyVectI*i_pcwtI*iPenalty;
              lnScoreFunctI = i_lnScoreFunctI;
              voteFractionI = i_voteFractionI;
              qualityLevelI = i_qualityLevelI;
              upperSimilarI = i_upperSimilarI;
              lowerSimilarI = i_lowerSimilarI;
              else % uLi = true => Ld subspace
              efficacyVectI = u_efficacyVectI;
              lnScoreFunctI = u_lnScoreFunctI;
              voteFractionI = u_voteFractionI;
              qualityLevelI = u_qualityLevelI;
              upperSimilarI = u_upperSimilarI;
              lowerSimilarI = u_lowerSimilarI;
              dPenalty = (1 - u_nConditionalI/nCpairs)^2;
              % -------------------------------
              efficacyVectD = d_efficacyVectD*d_pcwtD*dPenalty;
              lnScoreFunctD = d_lnScoreFunctD;
              voteFractionD = d_voteFractionD;
              qualityLevelD = d_qualityLevelD;
              upperSimilarD = d_upperSimilarD;
              lowerSimilarD = d_lowerSimilarD;
              end
            else             % => Ldi subspace
              efficacyVectD = d_efficacyVectD;
              lnScoreFunctD = d_lnScoreFunctD;
              voteFractionD = d_voteFractionD;
              qualityLevelD = d_qualityLevelD;
              upperSimilarD = d_upperSimilarD;
              lowerSimilarD = d_lowerSimilarD;
             % -------------------------------
              efficacyVectI = i_efficacyVectI;
              lnScoreFunctI = i_lnScoreFunctI;
              voteFractionI = i_voteFractionI;
              qualityLevelI = i_qualityLevelI;
              upperSimilarI = i_upperSimilarI;
              lowerSimilarI = i_lowerSimilarI;
            end  
% ----------------------------------------------------- for debugging only
%       disp( [dLd,uLd,uLi,iLi] ); 
%       disp( [lnScoreFunctD,voteFractionD,qualityLevelD,efficacyVectD] );
%       disp( [lnScoreFunctI,voteFractionI,qualityLevelI,efficacyVectI] );
%       pause     
% ----------------------------------------------------- transfer to arrays
         d_lnScoreFunctn10(j) = lnScoreFunctD;
         d_VoteConsensus10(j) = voteFractionD;
         d_qualityProperty(j) = qualityLevelD;
         d_efficacyFuncton(j) = efficacyVectD;
         d_lowerSimilarity(j) = lowerSimilarD;
         d_upperSimilarity(j) = upperSimilarD;
         %------------------------------------      
         i_lnScoreFunctn10(j) = lnScoreFunctI;
         i_VoteConsensus10(j) = voteFractionI;
         i_efficacyFuncton(j) = efficacyVectI;
         i_qualityProperty(j) = qualityLevelI;
         i_lowerSimilarity(j) = lowerSimilarI;
         i_upperSimilarity(j) = upperSimilarI;    
% *************************************************************** end hack        
               end                                % end of loop over modes            
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
% -------------------- sort by properties using logic to define partitions 
% Step 1: determine discriminant subspace:  3 conditions must be satisfied
Sd = ( d_lnScoreFunctn10 > lnMinScore1 );
Vd = ( d_VoteConsensus10 > voteThreshold );    % meets consensus threshold
Qd = ( d_qualityProperty > minQF );    % => exceeds minimum quality factor
Vd = and(Vd,Qd);                                    % a quality consensus!
Ld0 = and(Vd,Sd);                           % selects discriminant vectors
% Step 2: determine indifference subspace:  3 conditions must be satisfied
Si = ( i_lnScoreFunctn10 < lnMaxScore0 );
Vi = ( i_VoteConsensus10 > voteThreshold );    % meets consensus threshold
Qi = ( i_qualityProperty > minQF );    % => exceeds minimum quality factor
Vi = and(Vi,Qi);                                    % a quality consensus!
Li0 = and(Vi,Si);                           % selects discriminant vectors
% Step 3: define orthogonal sets of congruent subspaces
Ld = and(Ld0,~Li0);
Ldi = and(Ld0,Li0);
Lu = and(~Ld0,~Li0);
Li = and(~Ld0,Li0);
% Step 4: determine upper and lower levels within undetermined subspace
Lupper = ( d_lnScoreFunctn10 > lnMidScoreX );
Llower = ~Lupper;
LuUpper = and(Lu,Lupper);
LuLower = and(Lu,Llower);
Dd = sum(Ld);
Ddi= sum(Ldi);
Dd0= sum(Ld0);
DuU = sum(LuUpper);
DuL = sum(LuLower);
Du = DuU + DuL;
Di = sum(Li);
% step 5: rearrange congruent sections and sort them appropriately
EEVd = zeros(1,nModes);
EEVi = zeros(1,nModes);
QEVd = zeros(1,nModes);
QEVi = zeros(1,nModes);
USVd = zeros(1,nModes);
USVi = zeros(1,nModes);
LSVd = zeros(1,nModes);
LSVi = zeros(1,nModes);
SEVd = zeros(1,nModes);
SEVi = zeros(1,nModes);
CEVd = zeros(1,nModes);
CEVi = zeros(1,nModes);
SBV = zeros( size(U) );
Cind = zeros(1,nModes);
mapS2U = zeros(1,nModes);                     % => Uindex = mapS2U(Sindex)
indicator = zeros(1,nModes);
indicator(Ldi) = 2;
indicator(Ld) = 1;
indicator(Lu) = 0;
indicator(Li) = -1;
oldIndex = 1:nVL;
dPast = 0;
   if( Dd0 > 0 )
   d1st = dPast + 1;
   dEnd = dPast + Dd0;
   modeList = d1st:dEnd;
   mapIndex = oldIndex(Ld0);
   [~,newIndx] = sort( d_lnScoreFunctn10(Ld0), 'descend');
   indx = mapIndex(newIndx);
   EEVd(modeList) = d_efficacyFuncton(indx);
   QEVd(modeList) = d_qualityProperty(indx);
   USVd(modeList) = d_upperSimilarity(indx);
   LSVd(modeList) = d_lowerSimilarity(indx);
   SEVd(modeList) = exp( d_lnScoreFunctn10(indx) );
   CEVd(modeList) = d_VoteConsensus10(indx);
%-------------------------------------------  
   EEVi(modeList) = i_efficacyFuncton(indx);
   QEVi(modeList) = i_qualityProperty(indx);
   USVi(modeList) = i_upperSimilarity(indx);
   LSVi(modeList) = i_lowerSimilarity(indx);
   SEVi(modeList) = exp( i_lnScoreFunctn10(indx) );
   CEVi(modeList) = i_VoteConsensus10(indx);
   SBV(:,modeList) = U(:,indx);
   Cind(modeList) = indicator(indx);
%       if( sum( Cind(modeList) > 1 ) ~= Ddi )    % for debugging
%       disp( [Ld0; Li0; Ldi, indicator] );
%       disp( dividerLine );
%       disp( [indx; indicator(indx)] );
%       error('invalid sum-rule detected!');
%       end
   mapS2U(modeList) = indx;
   dPast = dEnd;
   end
% ------------------------------------------------------------------------
   if( DuU > 0 )
   d1st = dPast + 1;
   dEnd = dPast + DuU;
   modeList = d1st:dEnd;
   mapIndex = oldIndex(LuUpper);
   [~,newIndx] = sort( d_lnScoreFunctn10(LuUpper), 'descend');
   indx = mapIndex(newIndx);
   EEVd(modeList) = d_efficacyFuncton(indx);
   QEVd(modeList) = d_qualityProperty(indx);
   USVd(modeList) = d_upperSimilarity(indx);
   LSVd(modeList) = d_lowerSimilarity(indx);
   SEVd(modeList) = exp( d_lnScoreFunctn10(indx) );
   CEVd(modeList) = d_VoteConsensus10(indx);
%-------------------------------------------  
   EEVi(modeList) = i_efficacyFuncton(indx);
   QEVi(modeList) = i_qualityProperty(indx);
   USVi(modeList) = i_upperSimilarity(indx);
   LSVi(modeList) = i_lowerSimilarity(indx);
   SEVi(modeList) = exp( i_lnScoreFunctn10(indx) );
   CEVi(modeList) = i_VoteConsensus10(indx);
   SBV(:,modeList) = U(:,indx);
   Cind(modeList) = indicator(indx);
   mapS2U(modeList) = indx;
   dPast = dEnd;
   end
% ------------------------------------------------------------------------
   if( DuL > 0 )
   d1st = dPast + 1;
   dEnd = dPast + DuL;
   modeList = d1st:dEnd;
   mapIndex = oldIndex(LuLower);
   [~,newIndx] = sort( i_lnScoreFunctn10(LuLower), 'descend');
   indx = mapIndex(newIndx);
   EEVd(modeList) = d_efficacyFuncton(indx);
   QEVd(modeList) = d_qualityProperty(indx);
   USVd(modeList) = d_upperSimilarity(indx);
   LSVd(modeList) = d_lowerSimilarity(indx);
   SEVd(modeList) = exp( d_lnScoreFunctn10(indx) );
   CEVd(modeList) = d_VoteConsensus10(indx);
%-------------------------------------------  
   EEVi(modeList) = i_efficacyFuncton(indx);
   QEVi(modeList) = i_qualityProperty(indx);
   USVi(modeList) = i_upperSimilarity(indx);
   LSVi(modeList) = i_lowerSimilarity(indx);
   SEVi(modeList) = exp( i_lnScoreFunctn10(indx) );
   CEVi(modeList) = i_VoteConsensus10(indx);
   SBV(:,modeList) = U(:,indx);
   Cind(modeList) = indicator(indx);
   mapS2U(modeList) = indx;
   dPast = dEnd;
   end
% ------------------------------------------------------------------------
   if( Di > 0 )
   d1st = dPast + 1;
   dEnd = dPast + Di;
   modeList = d1st:dEnd;
   mapIndex = oldIndex(Li);
   [~,newIndx] = sort( i_lnScoreFunctn10(Li), 'descend');
   indx = mapIndex(newIndx);
   EEVd(modeList) = d_efficacyFuncton(indx);
   QEVd(modeList) = d_qualityProperty(indx);
   USVd(modeList) = d_upperSimilarity(indx);
   LSVd(modeList) = d_lowerSimilarity(indx);
   SEVd(modeList) = exp( d_lnScoreFunctn10(indx) );
   CEVd(modeList) = d_VoteConsensus10(indx);
%-------------------------------------------  
   EEVi(modeList) = i_efficacyFuncton(indx);
   QEVi(modeList) = i_qualityProperty(indx);
   USVi(modeList) = i_upperSimilarity(indx);
   LSVi(modeList) = i_lowerSimilarity(indx);
   SEVi(modeList) = exp( i_lnScoreFunctn10(indx) );
   CEVi(modeList) = i_VoteConsensus10(indx);
   SBV(:,modeList) = U(:,indx);
   Cind(modeList) = indicator(indx);
   mapS2U(modeList) = indx;
   end
% ----------------------------------------------------- for debugging only
% % % D4 = Dd + Ddi + Du + Di;
% % %    if( D4 ~= nVL )
% % %    error('2nd round: D4 is unequal to nVL'); 
% % %    end
% % %    if( Dd0 ~= Dd + Ddi )
% % %    error('dd is unequal to Dd + Ddi');
% % %    end

totalEfficacy = sum( d_efficacyFuncton(Ld) ) ...
              + sum( d_efficacyFuncton(Ldi) ) ...
              + sum( i_efficacyFuncton(Ldi) ) ...
              + sum( i_efficacyFuncton(Li) );
%%                                      package output into data structure
gvSPLOC.sVS = gvSPLOC.std_sVS;   % reset gvSPLOC.sVS to the standard value
splocResults = struct;
splocResults.sType = 'MCSPLOC';
splocResults.pursuitType = 0;             % This function always assumes 0
splocResults.pType = traitL.pType; 
splocResults.dim = traitL.dim; 
splocResults.baseFname = baseFname;
splocResults.SBV = SBV;                           % can have nModes < nDOF
splocResults.vT = voteThreshold;
splocResults.efficacy = totalEfficacy;
splocResults.Dd = Dd;        
splocResults.Ddi=Ddi;                
splocResults.Du = Du;
splocResults.Di = Di;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
splocResults.EEVd = EEVd;
splocResults.EEVi = EEVi;
splocResults.USVd = USVd;
splocResults.USVi = USVi;
splocResults.LSVd = LSVd;
splocResults.LSVi = LSVi;
splocResults.QEVd = QEVd;
splocResults.QEVi = QEVi;
splocResults.SEVd = SEVd;
splocResults.SEVi = SEVi;
splocResults.CEVd = CEVd; 
splocResults.CEVi = CEVi;
splocResults.Cind = Cind;
%%                                                  get mapU2S from mapS2U
mapU2S = zeros( size(mapS2U) );
   for Sindex=1:nModes
   Uindex = mapS2U(Sindex);
   mapU2S(Uindex) = Sindex;
   end
%%                                    summarize sploc results in sploc log
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
fid = fopen(splocLogFile,'a');
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('calculation summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['                base file name = ',baseFname]);
fprintf(fid,'%s \n',['                vote threshold = ', ...
                                          num2str(voteThreshold)]);
fprintf(fid,'%s \n',['                total efficacy = ', ...
                                          num2str(totalEfficacy)]);
D3 = Dd + Ddi;
   if( D3 > 0 )
   L3 = ( Cind > 0 );
   aveQEVd = sum( QEVd(L3) )/D3;
   fprintf(fid,'%s \n',['   mean quality for d-subspace = ', ...
                       num2str(aveQEVd)]);
   end
Di0 = sum(Li0);
   if( Di0 > 0 )
   aveQEVi = sum( QEVi(Li0) )/Di0;
   fprintf(fid,'%s \n',['   mean quality for i-subspace = ', ...
                       num2str(aveQEVi)]);
   end
fprintf(fid,'%s \n',['        discriminant modes: Dd = ',num2str(Dd)]);
fprintf(fid,'%s \n',['        dual-purpose modes: Ddi= ',num2str(Ddi)]);
fprintf(fid,'%s \n',['        undetermined modes: Du = ',num2str(Du)]);
fprintf(fid,'%s \n',['        indifference modes: Di = ',num2str(Di)]);
fprintf(fid,'%s \n',[' # of variables = Dd+Ddi+Du+Di = ',num2str(nVL)]);
fprintf(fid,'%s \n',['      # of labeled systems: nL = ',num2str(nL)]);
m0 = traitL.sampleSize; 
fprintf(fid,'%s \n',['sample size per labeled system = ',num2str(m0)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end
