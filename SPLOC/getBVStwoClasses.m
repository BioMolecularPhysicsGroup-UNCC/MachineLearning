function [splocResults,mapU2S,mapS2U] = ...
         getBVStwoClasses(U,baseFname,trait1,trait0,vT)
% Given a basis set of vectors U, determine its congruency spectrum
% Using SPLOC and binary classification
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
% -------------------------------------- define local language for sploc()
% binary states: 1 => "on" => Function    and    0 => "off" => Nonfunction
% M => mean column vector
% Q => symmetric matrix, such as a covariance matrix. 
% -------------------------------------------------------- local variables
% INPUT:
% U = proposed set of basis vectors
% baseFname => from which all output file names spawn
% ------------------------------------------------------------ from trait1
% Functional systems
% n1 = # of statistical metrics for Functional systems
% M1 = cell array containing n1 mean vectors for Functional systems
% Q1 = cell array containing n1 covariance matrices for Functional systems
% nV1 = # of variables for all Functional systems
% nD1 = total # of data samples = n1*samplesize1
% ------------------------------------------------------------ from trait0
% Nonfunctional systems
% n0 = # of statistical metrics for Nonfunctional systems
% M0 = cell array containing n0 mean vectors 
% Q0 = cell array containing n0 covariance matrices 
% nV0 = # of variables for all nonfunctional systems
% nD0 = total # of data samples = n0*samplesize0  
% --------- 
% vT = optional vote threshold   (default is calculated based on sampling)
%
% USAGE:
% ------------------------------------------------------------------------
%   splocResults = getBVS2(U,fname,trait1,trait0)
%   splocResults = getBVS2(U,fname,trait1,trait0,vT)
%           Notes: trait1 defines functional systems only    => class = 1
%                  trait0 defines nonfunctional systems only => class = 0
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   sploc ----> supervised projective learning for orthogonal congruences
% 
%   splocResults => quantifies the effectiveness of a set of basis vectors
%                   to characterize differences and similarities between a
%                   specified set of functional systems with respect to a
%                   specified set of nonfunctional systems
%
% PROCESS
% Emmulate the process found in sploc() identically, except the difference
% is that this is a one-time shot that quantifies the input vectors. No
% optimization is performed. 
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
%%                                     extract data from trait1 and trait0
n1 = trait1.n;
M1 = trait1.mu;
Q1 = trait1.cM;
nV1= trait1.nVariables;                                   % # of variables
Ns1= trait1.sampleSize;
nD1= trait1.nDtotal;      % nD1 = n1*sampleSize1 = total # of data samples
% ---------------------
n0 = trait0.n;
M0 = trait0.mu;
Q0 = trait0.cM;
nV0= trait0.nVariables;                                   % # of variables
Ns0= trait0.sampleSize;
nD0= trait0.nDtotal;      % nD0 = n0*sampleSize0 = total # of data samples
totPairs00 = n0*(n0-1)/2; % # of distinct nonfunction to nonfunction-pairs
totPairs11 = n1*(n1 - 1)/2;     % # of distinct function to function-pairs
totPairs = n0*n1;                   % # of [function to nonfunction]-pairs
%%                                                  error check input data
   if( strcmp(trait1.pType,trait0.pType) == 0 )
   error('packing type of trait1 is not the same as trait0');
   end
     if( trait1.dim ~= trait0.dim )
     error('local vector dimension of trait1 is not the same as trait0');
     end 
   if( n1 < 1 )
   error('n1 < 1: must have at least 1 representation of the 1-state');
   end
     if( n0 < 1 )
     error('n0 < 1: must have at least 1 representation of the 0-state');
     end
   if( nV1 ~= nV0 )
   error('input matrices are not the same size');
   end
nDOF = nV1;
% -------------------------------------------------- check properties of U
[mV,nModes] = size(U);
   if( mV ~= nDOF )
   error('number of DOF for U is not equal to number of variables');
   end
% nModes = # of modes based on user's request. Typically nModes = nDOF
%%                                                calculate vote threshold
   if( useDefault_vT )   % default uses heristic formula for voteThreshold
   % nD1 = n1*samplesize      and   nD0 = n0*sampleSize
   nS0 = nD0/sqrt(n0);          % nS0 = sqrt(n0)*sampleSize = nD0/sqrt(n0)
   nS1 = nD1/sqrt(n1);          % nS1 = sqrt(n1)*sampleSize = nD1/sqrt(n1)
   voteThreshold = getDefaultVoteThreshold(nS0,nS1,nDOF);
   end
%%                                             initialize SPLOC parameters
global gvSPLOC                         % shared information across toolset
% -------------------------------------------- scoring function parameters
minScore1 = gvSPLOC.minScore1;       % score > minScore1 =>  "on" subspace
maxScore0 = gvSPLOC.maxScore0;       % score < maxScore0 => "off" subspace
minQF = gvSPLOC.qualityThreshold;      % minimum allowed quality threshold
add0 = gvSPLOC.add0;                 % prevents variance ratios to diverge
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
%%                                                     prepare projections
% -------------------------- precalculate projected averages and variances
ave0 = zeros(n0,nModes);
var0 = zeros(n0,nModes);
ave1 = zeros(n1,nModes);
var1 = zeros(n1,nModes);
   for j=1:nModes                    % consider each j-th vector-direction
   vec = U(:,j);                                     % the selected vector
      for k0=1:n0                      % all nonfunctional 0-state systems
      ave0(k0,j) = vec'*M0{k0};                       % projected  average
      var0(k0,j) = max(vec'*Q0{k0}*vec,add0);         % projected variance
      end
      for k1=1:n1                         % all functional 1-state systems
      ave1(k1,j) = vec'*M1{k1};                       % projected  average
      var1(k1,j) = max(vec'*Q1{k1}*vec,add0);         % projected variance
      end
   end
%%                                         initialize counters of interest
% REMARK1:      upperSimilarity = max(fIndeterminable00,fIndeterminable11)
% REMARK2:      lowerSimilarity = min(fIndeterminable00,fIndeterminable11)
efficacyQuality = zeros(1,nModes);                         % mode efficacy
upperSimilarity = zeros(1,nModes); % greatest similarity between 00 and 11
lowerSimilarity = zeros(1,nModes);   % lowest similarity between 00 and 11
qualityProperty = zeros(1,nModes);                          % mode quality
lnScoreFunctn10 = zeros(1,nModes);         % => mean ln of selection power
dVoteFraction10 = zeros(1,nModes);  % => mean discriminant consensus power
iVoteFraction10 = zeros(1,nModes);  % => mean indifference consensus power
%%                         calculate properties of interest per projection
               for j=1:nModes                  % j-th requested projection
% --------------------------------- get <ln[scoringFunction]> for 10-pairs
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
                   lnScoreFunctn10(j) = lnScoreFunctn10(j) + x;
                   dVoteFraction10(j) = dVoteFraction10(j) + d_wt;
                   iVoteFraction10(j) = iVoteFraction10(j) + i_wt;
                   end      
                 end
% ---------------------------- selection & consensus power and consistency
         lnScoreFunctn10(j) = lnScoreFunctn10(j)/totPairs;  % => selection
         dVoteFraction10(j) = dVoteFraction10(j)/totPairs;  % => consensus
         iVoteFraction10(j) = iVoteFraction10(j)/totPairs;  % => consensus
         lnScoreFunction = lnScoreFunctn10(j);
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
% ---------------------------------------------------- nonfunctional cases
               stdv0 = sqrt( var0(:,j) );             % standard deviation
               aveMean0 = mean( ave0(:,j) );
               minMean0 =  min( ave0(:,j) );
               maxMean0 =  max( ave0(:,j) );
              %--------------------------------------
               tmpLnStDv0 = log(stdv0);        % log of standard deviation
               aveLnStDv0 = mean(tmpLnStDv0); 
               minLnStDv0 =  min(tmpLnStDv0);
               maxLnStDv0 =  max(tmpLnStDv0);
% ------------------------------------------------------- functional cases
               stdv1 = sqrt( var1(:,j) );             % standard deviation
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
%%                                          congruent subspace bifurcation
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& start of bifurcation of subspace selection
               if( lnScoreFunction > lnMidScoreX ) % discriminant subspace
% +++++++++++++++++++++++++++++++++++++++++++++++ calculate dQualityFactor
               vFraction = dVoteFraction10(j);
                %if( pursuitType ~= -1 )      % => requires dQualityFactor
                %disp( [j,gapMean,gapLnStDv] );
                 rFeature = sqrt(gapMean^2 + gapLnStDv^2);
                 vFrac = vFraction - voteThreshold;
                    if( vFrac < 0 )                % => poor quality level
                    vFrac = 0;
                    else                           % => good quality level
                    vFrac = vFrac/(1 - voteThreshold);  % normalize +range
                    end
                 fImax = max(fIndeterminable00,fIndeterminable11);
                %^^^^^-------------------------- only one needs to be good
                 required1 = max(-0.99, 8*(fImax - 0.5)^3 );
                 required2 = max(-0.99, 8*(0.5 - fIndeterminable)^3 );
                 required = min(required1,required2); %critical conditions
                 boost = 1 + required + consensus + consensus*vFraction;
                 dcp1 = cpMUg - cpMUo;
                 dcp2 = cpSDg - cpSDo;
                 dcp3 = dcp1 + dcp2;
                 dcp = max( [dcp1,dcp2,dcp3] );       % either one or both
                    if( dcp > 0 )
                    dcp = boost*dcp;
                    end
                 dAdd = 0;                      % => assume nothing to add     
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
                    end
                 dcp = dcp + dAdd;                % focus on the best case
                 dcp = dampLR*vFraction*dcp; 
                   if( dcp > 0 )                % => forward learning bias
                   xx = 0.1 + ...                      % sets minimum bias
                      + dcp/sqrt(1 + (dcp/9.9)^2);     % sets maximum bias
                   else      % => reverse bias w.r.t. minimum forward bias
                   xx = 0.1 + ...                % <= minimum forward bias
                      + dcp/sqrt(1 + (dcp/0.2)^2);  % maximum reverse bias
                   end
                 v = min(1,fImax*rFeature);
                 dQualityFactor = xx*( (1 - (1-v)^8 )^6 )*(0.5 + 0.5*v);
                 dQualityFactor = dQualityFactor*vFraction;
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
                %else             % => looking for indifference modes only
                %dQualityFactor = -0.1;     % <= default anti-discriminant
                %end
               qualityFactor = dQualityFactor;
               span = (lnScoreFunction - lnMidScoreX); 
               sWeight = sqrt(span) + span*(1 + span);
               %         ^^^^^^^^^^-------> rapidly bifurcates from span=0
                  if( qualityFactor > 0 )
                      if( qualityFactor > minQF )
                      vecScore = qualityFactor*sWeight;
                      else                          % tends to nucleate by
                      vecScore = minQF*sWeight;  % encouraging bifurcation
                      end      % ^^^^^------------------> boosts the score
                  else
                  vecScore = qualityFactor*sWeight;
                  end
               else % -----------------------------> indifference subspace
% +++++++++++++++++++++++++++++++++++++++++++++++ calculate iQualityFactor
               vFraction = iVoteFraction10(j);
                %if( pursuitType ~= 1 )       % => requires iQualityFactor
                 vFrac = vFraction - voteThreshold;
                    if( vFrac < 0 )                % => poor quality level
                    vFrac = 0;
                    else                           % => good quality level
                    vFrac = vFrac/(1 - voteThreshold);  % normalize +range
                    end
                 fImin = min( [fIndeterminable00, fIndeterminable, ...
                               fIndeterminable11] );
                 %^^^^^----------------> all three conditions are required
                 boost = 1 + consensus + vFraction ...
                       + consensus*vFraction ...
                       + max(-0.45, 8*(fImin - 0.5)^3 ) ...
                       + max(-0.45, 8*(fIndeterminable - 0.5)^3);
                 icp1 = cpMUo - cpMUg;                              % MEAN
                 icp2 = cpSDo - cpSDg;                         % std. dev.
                 icp3 = icp1 + icp2;
                 icp = min( [icp1,icp2,icp3] );  % require icp1 & icp2 > 0
                    if( icp > 0 )
                    icp = boost*icp;
                    end
                 iAdd = 0;                      % => assume nothing to add
                    if( (fImin > 0.5) && (icp > 0) )   % qualify for bonus
                    iAdd = (0.1 + 1.9*vFraction)* ...
                            max(0,fIndeterminable - 0.5)/ ...
                            max(0.01,1 - fImin) ... 
                         + 2*(fImin - 0.5) ...
                         + vFrac*consensus ...
                         + vFrac;                         
                    end
                 icp = icp + iAdd;               % focus on the worse case
                 icp = dampLR*vFraction*icp;
                   if( icp > 0 )                % => forward learning bias
                   xx = 0.1 + ...                      % sets minimum bias
                      + icp/sqrt(1 + (icp/9.9)^2);     % sets maximum bias
                   else      % => reverse bias w.r.t. minimum forward bias
                   xx = 0.1 + ...                % <= minimum forward bias
                      + icp/sqrt(1 + (icp/0.2)^2);  % maximum reverse bias
                   end
                 iQualityFactor = vFraction*xx;                 
                 % The non-linear term related to ensuring good clustering
                 % is not necessary for indifference modes because to be
                 % good everything has to be the same. But, multiplying by
                 % vFraction has the same effect as when done for d-modes.
                % --------------------------------------------------------
                %else             % => looking for discriminant modes only
                %iQualityFactor = -0.1;     % <= default anti-indifference
                %end                
               qualityFactor = iQualityFactor;
               span1 = (lnMidScoreX - lnScoreFunction); 
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
               %         ^^^^^^^^^^-------> rapidly bifurcates from span=0
                  if( qualityFactor > 0 )
                      if( qualityFactor > minQF )
                      vecScore = qualityFactor*sWeight;
                      else                          % tends to nucleate by
                      vecScore = minQF*sWeight;  % encouraging bifurcation
                      %          ^^^^^------------------> boosts the score
                      end
                  else
                  vecScore = qualityFactor*sWeight;
                  end
               vecScore = sVS*vecScore;           % => bias toward d-modes
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
               end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& end of bifurcation of subspace selection
               efficacyQuality(j) = vecScore;
               upperS = max(fIndeterminable00,fIndeterminable11); 
               lowerS = min(fIndeterminable00,fIndeterminable11);
               upperSimilarity(j) = upperS; 
               lowerSimilarity(j) = lowerS;
               qualityProperty(j) = qualityFactor;
               end                                % end of loop over modes
% ---------------------------------------------------------- sort by score
[lnScoreFunction,indx] = sort(lnScoreFunctn10,'descend');
efficacyScore = efficacyQuality(indx);
maxSimilarity = upperSimilarity(indx);
minSimilarity = lowerSimilarity(indx);
qualityFactor = qualityProperty(indx);
dVoteFraction = dVoteFraction10(indx);
iVoteFraction = iVoteFraction10(indx);
mapS2Uindices = indx;                 %=> Uindex = mapS2Uindices(1:nModes)
U = U(:,indx);                                 % to have consistent output
% --------------------------------------------------- assign voterFraction
voterFraction = dVoteFraction;        % assume all modes are in d-subspace
Si = (lnScoreFunction < lnMidScoreX);
voterFraction(Si) = iVoteFraction(Si);          % replace values as needed
Sd = ~Si;
% ------------------------------------------ using logic define partitions
Vf = ( voterFraction > voteThreshold );        % meets consensus threshold
qT = ( qualityFactor > minQF );        % => exceeds minimum quality factor
Vf = and(Vf,qT);                                    % a quality consensus!
Ld = and(Vf,Sd);                            % selects discriminant vectors
dd = sum(Ld);                                     % discriminant-dimension
   if( dd > 0 )
   dEfficacyFactor = efficacyScore(Ld);
   d_maxSimilarity = maxSimilarity(Ld);
   d_minSimilarity = minSimilarity(Ld);
   dQualityFeature = qualityFactor(Ld);
   dSelectionPower = exp( lnScoreFunction(Ld) ); 
   dConsensusPower = voterFraction(Ld);
   d_mapS2Uindices = mapS2Uindices(Ld);
   end
Li = and(Vf,Si);                            % selects indifference vectors
di = sum(Li); 
   if( di > 0 )
   iEfficacyFactor = efficacyScore(Li);
   i_maxSimilarity = maxSimilarity(Li);
   i_minSimilarity = minSimilarity(Li);
   iQualityFeature = qualityFactor(Li);
   iSelectionPower = exp( lnScoreFunction(Li) ); 
   iConsensusPower = voterFraction(Li);
   i_mapS2Uindices = mapS2Uindices(Li);
   end
Lu = and(~Ld,~Li);                 % not discriminant and not indifference 
du = sum(Lu);
   if( du > 0 )
   uEfficacyFactor = efficacyScore(Lu);
   uQualityFeature = qualityFactor(Lu);
   u_maxSimilarity = maxSimilarity(Lu);
   u_minSimilarity = minSimilarity(Lu);
   uSelectionPower = exp( lnScoreFunction(Lu) ); 
   uConsensusPower = voterFraction(Lu);
   u_mapS2Uindices = mapS2Uindices(Lu);
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                            prepare output information for basis vectors
CIP = zeros(1,nModes);
jF1 = 1;
jF2 = dd;
   if( jF2 >= jF1 )
   CIP(jF1:jF2) = 1;
   end
jU1 = jF2 + 1;
jU2 = jF2 + du;
   if( jU2 >= jU1 )
   CIP(jU1:jU2) = 0;
   end
jI1 = jU2 + 1;
jI2 = jU2 + di;
   if( jI2 >= jI1 )
   CIP(jI1:jI2) = -1;
   end
% disp( num2str([dd,du,di]) );
% disp( num2str(jI2) );
% disp( num2str(nDOF) );
   if( jI2 ~= nModes )               % this condition should be impossible
   error('Number of selection vectors is not equal to nModes');
   end
%%                                      combine all sections of selections
SBV = [];
   if( dd > 0 )
   EEV = dEfficacyFactor;
   QEV = dQualityFeature;
   USV = d_maxSimilarity;
   LSV = d_minSimilarity;
   SEV = dSelectionPower;
   CEV = dConsensusPower;
   SBV = U(:,Ld);
   mapS2U = d_mapS2Uindices;
      if( du > 0 )
      EEV = [EEV,uEfficacyFactor];
      QEV = [QEV,uQualityFeature];
      USV = [USV,u_maxSimilarity];
      LSV = [LSV,u_minSimilarity];
      SEV = [SEV,uSelectionPower];
      CEV = [CEV,uConsensusPower];
      SBV = [SBV,U(:,Lu)];
      mapS2U = [mapS2U,u_mapS2Uindices];
      end
      if( di > 0 )
      EEV = [EEV,iEfficacyFactor];
      QEV = [QEV,iQualityFeature];
      USV = [USV,i_maxSimilarity];
      LSV = [LSV,i_minSimilarity];
      SEV = [SEV,iSelectionPower];
      CEV = [CEV,iConsensusPower];
      SBV = [SBV,U(:,Li)];
      mapS2U = [mapS2U,i_mapS2Uindices];
      end
   else                                             % ------------> dd = 0
      if( du > 0 )
      EEV = uEfficacyFactor;
      QEV = uQualityFeature;
      USV = u_maxSimilarity;
      LSV = u_minSimilarity;
      SEV = uSelectionPower;
      CEV = uConsensusPower;
      SBV = U(:,Lu);
      mapS2U = u_mapS2Uindices;
         if( di > 0 )
         EEV = [EEV,iEfficacyFactor];
         QEV = [QEV,iQualityFeature];
         USV = [USV,i_maxSimilarity];
         LSV = [LSV,i_minSimilarity];
         SEV = [SEV,iSelectionPower];
         CEV = [CEV,iConsensusPower];
         SBV = [SBV,U(:,Li)];
         mapS2U = [mapS2U,i_mapS2Uindices];
         end
      else                                     % ------------> dd = du = 0
      EEV = iEfficacyFactor;
      QEV = iQualityFeature;
      USV = i_maxSimilarity;
      LSV = i_minSimilarity;
      SEV = iSelectionPower;
      CEV = iConsensusPower;
      SBV = U(:,Li);
      mapS2U = i_mapS2Uindices;
      end
   end
%%                                          map old format into new format
Ld = ( CIP == 1 );  
Lu = ( CIP == 0 );
Li = ( CIP == -1 );
refSEV = exp(lnMidScoreX);
Lupper = ( SEV > refSEV );
Llower = ~Lupper;
LuUpper = and(Lu,Lupper);
LuLower = and(Lu,Llower);
% -----------------------
EEVd = zeros(1,nModes);
EEVi = zeros(1,nModes);
QEVd = zeros(1,nModes);
QEVi = zeros(1,nModes);
USVd = zeros(1,nModes);
USVi = zeros(1,nModes);
LSVd = zeros(1,nModes);
LSVi = zeros(1,nModes);
SEVd = refSEV*ones(1,nModes);
SEVi = refSEV*ones(1,nModes);
CEVd = zeros(1,nModes);
CEVi = zeros(1,nModes);
% ---------------------------
EEVd(Ld) = EEV(Ld);
EEVd(LuUpper) = EEV(LuUpper);
EEVi(Li) = EEV(Li);
EEVi(LuLower) = EEV(LuLower);
% ---------------------------
QEVd(Ld) = QEV(Ld);
QEVd(LuUpper) = QEV(LuUpper);
QEVi(Li) = QEV(Li);
QEVi(LuLower) = QEV(LuLower);
% ---------------------------
USVd(Ld) = USV(Ld);
USVd(LuUpper) = USV(LuUpper);
USVi(Li) = USV(Li);
USVi(LuLower) = USV(LuLower);
% ---------------------------
LSVd(Ld) = LSV(Ld);
LSVd(LuUpper) = LSV(LuUpper);
LSVi(Li) = LSV(Li);
LSVi(LuLower) = LSV(LuLower);
% ---------------------------
LSVd(Ld) = LSV(Ld);
LSVd(LuUpper) = LSV(LuUpper);
LSVi(Li) = LSV(Li);
LSVi(LuLower) = LSV(LuLower);
% ---------------------------
SEVd(Ld) = SEV(Ld);
SEVd(LuUpper) = SEV(LuUpper);
SEVi(Li) = SEV(Li);
SEVi(LuLower) = SEV(LuLower);
% ---------------------------
CEVd(Ld) = CEV(Ld);
CEVd(LuUpper) = CEV(LuUpper);
CEVi(Li) = CEV(Li);
CEVi(LuLower) = CEV(LuLower);
% ---------------------------
CEVd(Ld) = CEV(Ld);
CEVd(LuUpper) = CEV(LuUpper);
CEVi(Li) = CEV(Li);
CEVi(LuLower) = CEV(LuLower);
% ---------------------------
Cind = CIP;
Lefficacy = ~Lu;
totalEfficacy = sum( EEV(Lefficacy) );
%%                                      package output into data structure
gvSPLOC.sVS = gvSPLOC.std_sVS;   % reset gvSPLOC.sVS to the standard value
splocResults = struct;
splocResults.sType = 'SPLOC';
splocResults.pursuitType = 0;             % This function always assumes 0
splocResults.pType = trait1.pType; 
splocResults.dim = trait1.dim; 
splocResults.baseFname = baseFname;
splocResults.SBV = SBV;                           % can have nModes < nDOF
splocResults.vT = voteThreshold;
splocResults.efficacy = totalEfficacy;
splocResults.Dd = dd;
splocResults.Ddi = 0;
splocResults.Du = du;
splocResults.Di = di;
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

% % % % % % % % % % % % % % %%                                                  get mapS2U from mapU2S
% % % % % % % % % % % % % % mapS2U = zeros( size(mapU2S) );
% % % % % % % % % % % % % %    for Uindex=1:nModes
% % % % % % % % % % % % % %    Sindex = mapU2S(Uindex);
% % % % % % % % % % % % % %    mapS2U(Sindex) = Uindex;
% % % % % % % % % % % % % %    end

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
   if( dd > 0 )
   aveQEVd = sum( QEVd(Ld) )/dd;
   fprintf(fid,'%s \n',['   mean quality for d-subspace = ', ...
                       num2str(aveQEVd)]);
   end
   if( di > 0 )
   aveQEVi = sum( QEVi(Li) )/di;
   fprintf(fid,'%s \n',['   mean quality for i-subspace = ', ...
                       num2str(aveQEVi)]);
   end                                                                 
fprintf(fid,'%s \n',['        discriminant modes: Dd = ',num2str(dd)]);
fprintf(fid,'%s \n',['        undetermined modes: Du = ',num2str(du)]);
fprintf(fid,'%s \n',['        indifference modes: Di = ',num2str(di)]);
fprintf(fid,'%s \n',[' # of variables = Dd + Du + Di = ',num2str(nV1)]);
% ---------------------- information on orgin of nonfunctional system data
fprintf(fid,'%s \n',['  nonfunctional system dataset = ', ...
                        trait0.dataRefName]);
fprintf(fid,'%s \n',['  # of nonfunction systems: n0 = ',num2str(n0)]);
m0 = trait0.sampleSize; 
fprintf(fid,'%s \n',['      sample size per system-0 = ',num2str(m0)]);
fprintf(fid,'%s \n',['total sample size for system-0 = ',num2str(nD0)]);
% ------------------------- information on orgin of functional system data
fprintf(fid,'%s \n',['     functional system dataset = ', ...
                        trait1.dataRefName]);
fprintf(fid,'%s \n',['     # of function systems: n1 =  ',num2str(n1)]);
m1 = trait1.sampleSize; 
fprintf(fid,'%s \n',['      sample size per system-1 = ',num2str(m1)]);
fprintf(fid,'%s \n',['total sample size for system-1 = ',num2str(nD1)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end
