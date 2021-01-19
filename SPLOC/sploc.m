function splocResults = sploc(o,U0,baseFname,trait1,trait0,verbosity,vT) 
% Supervised projective Learning with Orthogonal Completeness
% Written by Donald J. Jacobs at UNC Charlotte (initial code Oct 24, 2018)
% SPLOC is a binary classifer. It requires two sets of labeled data.
% 
% DEFINITIONS
% ------------------------------------------------------- trait definition
% trait => data structure  (required components for sploc)
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
% o = 1 or 0 or -1 => type of SPLOC => d, or (d & i) or i
% U0 = initial complete basis set of vectors or the simple number 0.
% baseFname = all output file names spawn off this base file name.
% --------- from trait1
% n1 = # of statistical metrics for Functional systems
% M1 = cell array containing n1 mean vectors for Functional systems
% Q1 = cell array containing n1 covariance matrices for Functional systems
% nV1 = # of variables in each Functional system
% Ns1 = trait1.sampleSize;
% nD1 = total # of data samples = n1*samplesize1 for Functional systems
% --------- from trait0
% n0 = # of statistical metrics for Nonfunctional systems
% M0 = cell array containing n0 mean vectors for Nonfunctional systems
% Q0 = cell array containing n0 cov-matrices for NonFunctional systems
% nV0 = # of variables in each Nonfunctional system
% Ns0 = trait0.sampleSize;
% nD0 = total # of data samples = n0*samplesize0 for Nonfunctional systems
% --------- 
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
% vT = optional vote threshold   (default is calculated based on sampling)
%
% USAGE:
%      splocResults = sploc(o,0,baseFname,trait1,trait0)
%      splocResults = sploc(o,0,baseFname,trait1,trait0,verbosity)
%      splocResults = sploc(o,0,baseFname,trait1,trait0,verbosity,vT)
%                           ^ ^---> no initial guess for basis set
%                           |-----> SPLOC pursuit objective
%
%      splocResults = sploc(o,1,baseFname,trait1,trait0)
%      splocResults = sploc(o,1,baseFname,trait1,trait0,verbosity)
%      splocResults = sploc(o,1,baseFname,trait1,trait0,verbosity,vT)
%                           ^ ^---> identity matrix: coordinate basis set
%                           |-----> SPLOC pursuit objective
%
%      splocResults = sploc(o,U0,baseFname,trait1,trait0)
%      splocResults = sploc(o,U0,baseFname,trait1,trait0,verbosity)
%      splocResults = sploc(o,U0,baseFname,trait1,trait0,verbosity,vT)
%                           ^ ^^--> user defined initial basis set
%                           |-----> SPLOC pursuit objective
%
%                                  AVAILABLE MODES
%       o = 2 => maximize discriminant modes & MINIMIZE indifference modes
%                RANU is modified to give inversion of i-modes
%                sets reference value for sVS = 0.1 = set_sVS (standard)
%       o = 1 => focus on maximizing discriminant modes
%                sets the reference value for sVS = 0.01 = set_sVS
%                use adaptive sVS: 0.1 < sVS/set_sVS < 10  depending on Dd
%NORMAL o = 0 => maximize discriminant & indifference modes simultaneously
%                sets reference value for sVS = 0.1 = set_sVS (standard)
%                use adaptive sVS: 0.1 < sVS/set_sVS < 10  depending on Dd
%       o= -1 => focus on maximizing indifference modes
%                sets the reference value for sVS = 1.0 = set_sVS
%                use adaptive sVS: 0.1 < sVS/set_sVS < 10  depending on Di
%       o= -2 => maximize indifference modes & MINIMIZE discriminant modes
%                RANU is modified to give inversion of d-modes
%                sets reference value for sVS = 0.1 = set_sVS (standard)
% 
% RARE^ use o = 1 or -1 when a bias toward d-modes or i-modes is desired.
%               When the bias for d-modes increases, the bias for i-modes
%               is unaltered. Thus, i-modes only get removed because a 
%               d-mode can be created. Vice versa, when the bias for an
%               i-mode increases, the bias for d-modes is unaltered. Thus
%               i-modes only get removed because a i-mode can be created.
%
% RARE^ use o = 2 when a bias for d-modes and AGAINST i-modes is desired.
%               Generating i-modes decresses efficacy.
%
% RARE^ use o = -2 when a bias for i-modes and AGAINST d-modes is desired.
%               Generating i-modes decresses efficacy.
%
% ^=> This would be used if the NORMAL case does not work. Which way is 
%     best to rescue the normal case is up to the user. 
%
% PROCESS
% Apply projective learning to discriminant two categoric classifications.
% Random vector-directions are generated efficiently,  
%    For each vector-direction considered individually calculate a score:
%    score = scoring_function[vector-direction]
%    augmented score attempts to disperse vector-directions optimally.
%    selection power and consensus power is optimized using a generalized
%    Jacobi method to optimize the scoring function. 
%
%    Classification of vector directions is done via three conditions
%
% ------------------ three conditions -----------------       congruency
% score > minScore1  &  dVote > vT  &  dQuality > minQF  =>  discriminant
% score < maxScore0  &  iVote > vT  &  iQuality > minQF  =>  indifference
%                otherwise                               =>  undetermined
%                            
% cip = congruency indicator partitions
%       1 => congruent discriminant subspace     
%       0 => congruent undetermined subspace
%      -1 => congruent indifference subspace
%
% OUTPUT splocResults <-- data structure
% ========================================================================
% splocResults.sType       = type of spectrum: SPLOC
% splocResults.nSols       = number of solutions in this data structure
% splocResults.rankScore   = ranks multiple solutions based on consistency
% splocResults.pursuitType = 2,1,0,-1,-2 
% splocResults.pType       = packing format describing the vector space
% splocResults.dim         = # of components in a local vector
% splocResults.baseFname   = base file name
% splocResults.SBV         = selection basis vectors
% splocResults.vT          = voting threshold to establish consensus
% splocResults.efficacy    = quantifies the ability for SBV to cluster
% splocResults.Dd          = # of discriminant only modes
% splocResults.Du          = # of undetermined modes
% splocResults.Di          = # of indifference only modes
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% splocResults.EEV = efficacy eigenvalues
% splocResults.USV = upper similar eigenvalues => max(similar11,similar00)
% splocResults.LSV = lower similar eigenvalues => min(similar11,similar00)
% splocResults.QEV = quality eigenvalues
% splocResults.SEV = selection eigenvalues       ** used for sort ordering
% splocResults.CEV = consensus eigenvalues
% splocResults.CIP = congruency indicator partition (1,0,-1)
%                    1 => projections for discriminant congruences
%                    0 => projections that are undetermined
%                   -1 => projections for indifference congruences
% ========================================================================

%%                                          General history if development
% Modified substantially for first released version, March 28, 2019.
% Optimization Method was initially a random search over basis rotations.
% Modified version uses 2D subspace rotation method for optimization of
% selection and consensus spectrums, which is a generalized Jacobi method.
% Furthermore, the scoring function was simplified considerably, allowing
% for separability between modes. Although many changes were made to 
% improve performance, the paradigm remains identical. Finished 04/14/19.
% Modified in Oct 13-20, 2019 to include initial basis set, U0, to start 
% from. In addition, the scoring function was substantially changed to a 
% Rectified Adaptive Nonlinear Unit (RANU) for learning to obtain efficacy
% as an objective. A pursuit type was added to modify the shape of the
% RANU to generate d, d&i or i modes. The initial guess, if U0 is not 
% supplied, is selected based on PCA for different scenarios. More work 
% was done on the sploc() function to fix a bug Chris Avery spotted. This
% bug [see patience] was removed, but was also thought to be responsible
% for occasional unexplained excessive slowness. Searching through the
% code further revealed a few more bugs that were fixed while fine tunning 
% was applied in the optimization algorithms. Much of the code was then
% streamlined, and a [gremlin] was expelled. This was another place that
% sploc() could get excessively slow. Most important changes are the 
% introduction of mode-type probabilities, initial probabilities are
% randomized, sequential resets are made in a simpler way, including 
% the stopping of anemic performance, and the gremlin was found related
% to skipping rotations. By adding a new concept of quorum the gremlin was
% expelled. a LOT of MINOR & MAJOR changes were made [~800 lines] keeping
% the basic ideas the same, but adding several NEW features that help the
% importance sampling. The front end of sploc() was modified slighly too.
% The most important added concept was instead of trying to get all the 
% modes to converge accurately at the same time, the program focuses down
% on the side of spectrum MOST important to converge first, and the rest
% of the modes follow like a domino effect. This concept also includes 
% an adaptive multiplicaive factor on the relative weighting between the
% d-modes and i-modes. This was always a source of confusion, but now the
% program internally adaptively optimizes this factor. At the end, sploc
% has much faster convergence, runs faster, while achieves much greater
% accuracy. Moreover, sensitively to initial conditions has been greatly 
% reduced. This latest version [called 2020b release] was finished on 
% April 28, 2020. 
%%                                               parse input into function
useDefault_vT = true;        % default voting threshold will be calculated
   switch nargin
     case 5
     verbosity = 1;                     % => default value (0 is reserved)
     case 6
     verbosity = setVerbosity(verbosity);  % restricts values to {0,1,2,3}
     case 7
     verbosity = setVerbosity(verbosity);  % restricts values to {0,1,2,3}
        if( vT < 0.5 )
        error('A voting threshold < 50% is not a majority vote!');
        elseif( vT > 0.95 )
        error('95% is the maximum voting threshold used in sploc!');
        end
     useDefault_vT = false;
     voteThreshold = vT;
     otherwise
       if( nargin < 5 )
       error('must specify: {o, U0, baseFileName, trait1, trait0}');
       end
   end
%%                                       determine SPLOC pursuit objective
   if( o < 0 )
   pursuitType = -1;
   o = -2;                                         % define true o setting
   elseif( o > 1 )
   pursuitType = 1;
   o = 2;                                          % define true o setting
   else
   pursuitType = 0;
   end
%%                                     extract data from trait1 and trait0
n1 = trait1.n;
M1 = trait1.mu;
Q1 = trait1.cM;
nV1= trait1.nVariables;                                   % # of variables
nD1= trait1.nDtotal;      % nD1 = n1*sampleSize1 = total # of data samples
Ns1 = trait1.sampleSize;                               % number of samples
% ---------------------
n0 = trait0.n;
M0 = trait0.mu;
Q0 = trait0.cM;
nV0= trait0.nVariables;                                   % # of variables
Ns0 = trait0.sampleSize;                               % number of samples
nD0= trait0.nDtotal;      % nD0 = n0*sampleSize0 = total # of data samples
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
% -------------------------------
   if( nV1 ~= nV0 )
   error('input matrices are not the same size');
   end
nDOF = nV1;
%%                                                  check properties of U0
flagUseU0 = -1;      % assume user provides no initial guess for basis set
[nU0,mU0] = size(U0);
   if( (nU0 == 1) && (mU0 == 1) )
      if( U0 == 1 )
      U0 = eye(nDOF);       % => user wishes to use identity matrix for U0
      flagUseU0 = 1; %=> user provides a valid initial guess for basis set
      msgInitial_U = 'initial basis set is the identity matrix';
      elseif( U0 ~= 0 )
      error('U0=0 => default,  U0=1 => identity,  U0=orthonormal basis');
      end
   else   
      if( (nU0 ~= nDOF) || (mU0 ~= nDOF) )
      error('initial guess for basis set has wrong dimensions');
      end
   temp = U0*U0';       
   t1 = ones(1,nDOF);
   temp1 = diag(t1);                         % this is the identity matrix
   errorTolerance = mean( mean( abs(temp - temp1) ) );
      if( errorTolerance > 0.000002 ) 
      disp(['error level = ',num2str(errorTolerance)]);
      error(['given U basis set is not orthonormal ', ... 
             'within required precision']); 
      else
      msgInitial_U = 'initial basis set given as input';
      end
   flagUseU0 = 1;   % => user provides a valid initial guess for basis set
   end
%%                                                calculate vote threshold
   if( useDefault_vT )   % default uses heristic formula for voteThreshold
   % nD1 = n1*samplesize  and   nD0 = n0*sampleSize
   nS0 = nD0/sqrt(n0);          % nS0 = sqrt(n0)*sampleSize = nD0/sqrt(n0)
   nS1 = nD1/sqrt(n1);          % nS1 = sqrt(n1)*sampleSize = nD1/sqrt(n1)
   voteThreshold = getDefaultVoteThreshold(nS0,nS1,nDOF);
   end
vT = voteThreshold;
%%                                                      set figure_number0
   if( verbosity == 3 )
   figure;
   figure_number0 = get(gcf,'Number');
   clf
   end
%%                                             initialize SPLOC parameters
t0 = cputime;
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
middleScore = exp(lnMidScoreX);
gap_lnScore = lnMinScore1 - lnMidScoreX;
% --------------------------------------------- set relative score weights
   switch o
       case -2
       set_sVS = gvSPLOC.std_sVS;
       case -1
       set_sVS = 1.00;     %0.25    % i-modes are 2.5 times more effective
       case 0                           
       set_sVS = 0.10;                        % should be equal to std_sVS
       case 1
       set_sVS = 0.01;     %0.04    % d-modes are 2.5 times more effective
       otherwise   % 
       set_sVS = gvSPLOC.std_sVS;
   end
   if( abs( gvSPLOC.std_sVS - 0.1 ) > 0.0000001 )
   error('playing around with global variables is dangerous!');
   end
%%                                             Initialize weight functions
dWt = gvSPLOC.dWtS;               % discriminant probability weight factor
iWt = gvSPLOC.iWtS;               % indifference probability weight factor
%uWt = gvSPLOC.uWtS;              % undetermined probability weight factor
xMax = gvSPLOC.xMax;                     % maximum value of score required
%%                         prepare for projection pursuit over the dataset
totPairs00 = n0*(n0-1)/2; % # of distinct nonfunction to nonfunction-pairs
totPairs11 = n1*(n1 - 1)/2;     % # of distinct function to function-pairs
totPairs = n0*n1;                   % # of [function to nonfunction]-pairs
% ------------------------------- convergence and random search parameters
% applied to randomly selecting distinct pairs of orthonormal vectors that
% are optimally rotated in a 2D subspace to maximize an objective scoring 
% function.  Funnel Diffusion is not applied in this optimization process.
nvp = nDOF*(nDOF-1)/2;                            % number of vector pairs
nRounds = 15 + max(0,37 - nDOF);     % smaller systems require more times;
minFrac = 1/sqrt(2*nDOF);
targetLesson = 5*max(25,nDOF);
maxRotations = 12*nvp + ceil( nvp*28/sqrt(nDOF/2) );
minRotations0 = max(2*targetLesson, ceil(minFrac*nvp) );
maxLessonBatch = 5*minRotations0/targetLesson;     % maximum lesson length
minConvergenceRatio = 0.005;          % => a 0.5% change in consensus vote
%                     ^^^^^---> percent error = 100*0.005/(consensus vote)
writeCheckPoint = max( ceil(0.1*nDOF),25);
% =========================================== growth factor considerations
growFactor = 2^0.2;    % => 0.5*[2^(1/5)]^5 = 1 => 5 skips per 1 hit @ 0.5
endGameRatio = 2*minConvergenceRatio;
nTimes = 1 + ceil(0.20*sqrt(nDOF+25)*nRounds);           % a heristic only
%disp(nTimes);
%%                                             setup output file base name
logFileName = [baseFname,'_sploc.log'];    % => for a specialized log file
subFolder = 'training';
logFileName = getOutputFileName(subFolder,logFileName);
SPLOCresultsFileName = [baseFname,'_splocResults.dlm'];
splocFileName = getOutputFileName(subFolder,SPLOCresultsFileName);
%%                                        write to training sploc log file
fid = fopen(logFileName,'w');
msgLine = ['----------------------------------------', ...
           '----------------------------------------'];
msg ='Supervised Projective Learning with Orthogonal Completeness (SPLOC)';
fprintf(fid,'%s \n',msg);
msg ='sploc()';
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',msgLine);
msg = '0-state: (mean, std) = (ave0,sigma0)';
fprintf(fid,'%s \n',msg);
msg = '1-state: (mean, std) = (ave1,sigma1)';
fprintf(fid,'%s \n',msg);
msg = 'signal = |ave1-ave0|';
fprintf(fid,'%s \n',msg);
msg = 'noise = sqrt[ sigma0^2 + sigma1^2 ]';
fprintf(fid,'%s \n',msg);
msg = 'snr= signal to noise ratio = signal/noise';
fprintf(fid,'%s \n',msg);
msg = 'sbn = signal beyond noise ratio = max(0,snr-1)';
fprintf(fid,'%s \n',msg);
msg = 'r = max( sigma0/sigma1, sigma1/sigma0 ) - 1';
fprintf(fid,'%s \n',msg);
msg = ['----------------------------------------', ...
       '----------- scoring function definition:'];
fprintf(fid,'%s \n',msg); 
msg = [num2str(gvSPLOC.maxScore0),' < score < ', ... 
       num2str(gvSPLOC.minScore1),' => undetermined congruency'];
fprintf(fid,'%s \n',msg);
msg = ['lnMaxScore0 = maximum score for indifference = ln(', ... 
      num2str(gvSPLOC.maxScore0),')'];
fprintf(fid,'%s \n',msg);
msg = ['lnMinScore1 = minimum score for discriminant = ln(', ...
      num2str(gvSPLOC.minScore1),')'];
fprintf(fid,'%s \n',msg);
msg = 'lnMidScoreX = (lnMaxScore0 + lnMinScore1)/2';
fprintf(fid,'%s \n',msg);
msg = 'lnScore0 = ln(sqrt(r*r + snr*snr) + 1)';
fprintf(fid,'%s \n',msg);
msg = 'lnScore1 = ln(sqrt(r*r + sbr*sbr) + 1)';
fprintf(fid,'%s \n',msg);
msg = '   if( lnScore0 < lnMidScoreX )';
fprintf(fid,'%s \n',msg);
msg = '   lnScore = lnScore0';
fprintf(fid,'%s \n',msg);
msg = '   elseif( lnScore1 > lnMidScoreX )';
fprintf(fid,'%s \n',msg);
msg = '   lnScore = lnScore1';
fprintf(fid,'%s \n',msg);
msg = '   else';
fprintf(fid,'%s \n',msg);
msg = '   lnScore = lnMidScoreX';
fprintf(fid,'%s \n',msg);
msg = '   end';
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',msgLine);
%%                                         initialize 2D rotation matrices
a_upper = pi/4;
a_lower = -a_upper;
Np2 = 2^14;              
da = (a_upper - a_lower)/(1 + Np2);
aR = round(Np2/2) + 520;             % must add 512 but an extra 8 is fine
aL = -aR ;
a = aL:aR;
a = da*a;
[~,iaZero] = min( abs(a) );
Na = length(a);
dampFactor = 0.5 + 0.598659*abs(a);  % range: 0.5 to 1 for 0 < |a| < a_max
% plot(a,dampFactor)
% hold on
% plot(a,dampFactor*growFactor^5);
% error('stop here for now');
%%                                                    initialize zoomLevel
% ======================================================== search schedule
%                                            REMARKS:  dA = 0.0055 degrees
iaSkip = 1024;                     % => skip resolution = dA*1024 = 5.6320
aR1 = floor(iaSkip/2) + iaSkip*7;          % => span1 = dA*2*aR1 = 84.4800
aR2 = iaSkip + iaSkip*6;                     % => interlace between span 1
shiftL = [aR1,aR2,128, 64,32,16, 8, 4,2, 1]; %levels (1&2)=> resol = 2.816
skipL = [1024,1024,256,128,64,32,16,8,4, 2];  % resolution = dA*1 = 0.0055
% zLevel:  1    2   3   4   5  6  7 8 9 10 # of evals = (8 + 7 + 8)*2 = 46
zoomLevelEnd = 10;                        % level 1 evals alone = 8*2 = 16
% ========================================================================
% --------------------------------------- build/store 2D rotation matrices
fs = sin(a);                                             % sine   function
fc = cos(a);                                             % cosine function
TRcell = cell(1,Na);     % Transformation from Right side: rotation matrix
TLcell = cell(1,Na);     % Transformation from Left  side: rotation matrix
% ---------------------------------- construct rotation matrices in memory
   for ia=1:Na
   TR = [ fc(ia), fs(ia); -fs(ia), fc(ia) ];          % => rotation matrix
   TRcell{ia} = TR; 
   TLcell{ia} = TR';                    % TL = transpose(TR) = inverse(TR)
   %      NOTE: ^^^-----> takes ~2x more time to take transpose vs. recall
   %disp( TRcell{ia} );
   %disp( TLcell{ia} );
   %pause
   end
%%                   initialize variables for random vector-pair selection
nVectorPairs = nDOF - 1;    % for fixed j1: nDOF - 1 distinct vector pairs
j2Pair = 1:nDOF;   % => each vector is covered once: order does not matter
nDOF1 = nDOF + 1;
% -------------------------------------------- define waiting time weights
vIndx = (1:nDOF)/nDOF1;
vIndx = vIndx';
vPROB = zeros(nDOF,2);
vPROB(:,1) = sqrt(vIndx);
vPROB(:,2) = 1 - sqrt(vIndx); 
clear vIndx;
% ---------------------------------------- define j1 selection probabilies
   switch o
     case -2
     prob4j12Bd = 0;  % probability for j1 vector to be selected as d-mode
     case -1
     prob4j12Bd = 0.1;
     case 0
     prob4j12Bd = 0.5;
     case 1
     prob4j12Bd = 0.9;
     otherwise
     prob4j12Bd = 1;
   end
% ----------------------------------------- define cummulative probabilies
   switch o
       case -2
       cPROBii = 0.60;
       cPROBid = 0.80;  
       cPROBdi = 1.00;  
       %cPROBdd = 1.00;
       overRideRank = nDOF;
       case -1
       cPROBii = 0.60;
       cPROBid = 0.80;  
       cPROBdi = 1.00;  
       %cPROBdd = 1.00;
       overRideRank = nDOF;
       case 0                             % ii-10%  id-35%  di-45%  dd-10%
       cPROBii = 0.10;                  % runs SLOW generates lot of skips
       cPROBid = 0.45;  % runs with moderate speed via moderate # of skips
       cPROBdi = 0.90;  % runs with moderate speed via moderate # of skips
       %cPROBdd = 1.00;                    % runs fast generates few skips
       overRideRank = 1;
       case 1
       cPROBii = 0.00;
       cPROBid = 0.20;  
       cPROBdi = 0.40;  
       %cPROBdd = 1.00;
       overRideRank = 1;
       otherwise   % 
       cPROBii = 0.00;
       cPROBid = 0.20;  
       cPROBdi = 0.40;  
       %cPROBdd = 1.00;
       overRideRank = 1;
   end
overRideZone = 0.05;              % defines when extra scrutiny is applied
                          % whenever a vote is close to the vote threshold
% ------------------------------ set reference parameters for adaptive sVS
% adaptive_sVS = 10^pow    |pow| = (1 - r)^p4t  at   ro = sqrt(nDOF)/nDOF1
%                             t0 = 1 - ro                Let  t0^p4t = 0.1  
%                            p4t = log(0.1)/log(t0) 
ttt = 1 - sqrt(nDOF)/nDOF1;
p4t = log(0.1)/log(ttt);
%disp(['adaptive exponent = ',num2str(p4t)]); 
%%               initialize parameters for control in optimization process
pFmaxStartValue = 0.8;                             % 0.7 to 0.8 works well
sWeightReference = sqrt(gap_lnScore) + gap_lnScore*(1 + gap_lnScore);
bestVector2DscoreREF = minQF*sWeightReference;
pREF = 2 + bestVector2DscoreREF;
pREF2 = 2*pREF;
iWaitTime10percent = 1 + ceil(0.1*nDOF);    %=> looks at ~10% of the modes
upperDrotateFraction = 1.0;
lowerDrotateFraction = 0.3;
ref_randomMatrixAmplitude = 0.05;           % reference rotation amplitude
maxRMA = 5*ref_randomMatrixAmplitude;       % need not go larger than this
minRMA = ref_randomMatrixAmplitude/5;      % need not go smaller than this
Drotate = 0;
sumDrotate = 0;
Nrotate = 1.0e-20;
quorum = round( sqrt(nDOF-1) );
%%                                            create initial basis vectors
   if( flagUseU0 < 0 )
% -------------------------- ave statistical matrix for functional systems
   a0 = zeros(nDOF,1);
      for k0=1:n0
      a0 = a0 + M0{k0};
      end
   a0 = a0/n0;
   a1 = zeros(nDOF,1);
      for k1=1:n1
      a1 = a1 + M1{k1};
      end
   a1 = a1/n1;
   aa = (a0 + a1)/2;
% ------------------------------
   R1 = zeros(nDOF); 
   S1 = zeros(nDOF);
      for k1=1:n1
      R1 = R1 + Q1{k1};
      S1 = S1 + (M1{k1} - aa)*(M1{k1} - aa)';
      end
   R1 = R1/n1;
   S1 = S1/n1;
% ------------------------------
   R0 = zeros(nDOF); 
   S0 = zeros(nDOF);
      for k0=1:n0
      R0 = R0 + Q0{k0}; 
      S0 = S0 + (M0{k0} - aa)*(M0{k0} - aa)';
      end
   R0 = R0/n0;  
   S0 = S0/n0;
% ------------------------------------------------------------------------
   R2 = R0 + (a0 - aa)*(a0 - aa)' ...
      + R1 + (a1 - aa)*(a1 - aa)';
   T0 = R0 + S0;
   T1 = R1 + S1;
   T2 = T0 + T1;
   R3 = R1 + (a1 - a0)*(a1 - a0)';
   R4 = R0 + (a0 - a1)*(a0 - a1)';
% ---------------------------------- try 8 complete orthonormal basis sets
   [U1,D1] = eig(R0);                  % pooled over nonfunctional systems
   [U2,D2] = eig(R1);                     % pooled over functional systems
   [U3,~] = eig(R2);   % averaged over pooled functional and nonfunctional
   [U4,D4] = eig(T0);                              % more combinations ...
   [U5,D5] = eig(T1); 
   [U6,~] = eig(T2);
   [U7,~] = eig(R3); 
   [U8,~] = eig(R4);
% ------------------------------------ check against numerical instability
   U0 = eye(nDOF);
     if( ~isreal(U1) )
     U1 = U0;
     else
     temp = U1*U1';       
     t1 = ones(1,nDOF);
     temp1 = diag(t1);                       % this is the identity matrix
     errorTolerance = mean( mean( abs(temp - temp1) ) );
       if( errorTolerance > 0.000002 )
       U1 = U0;
       end
     end
   % ------------------
     if( ~isreal(U2) )
     U2 = U0;
     else
     temp = U2*U2';       
     t1 = ones(1,nDOF);
     temp1 = diag(t1);                       % this is the identity matrix
     errorTolerance = mean( mean( abs(temp - temp1) ) );
       if( errorTolerance > 0.000002 )
       U2 = U0;
       end
     end
   % ------------------
     if( ~isreal(U3) )
     U3 = U0;
     else
     temp = U3*U3';       
     t1 = ones(1,nDOF);
     temp1 = diag(t1);                       % this is the identity matrix
     errorTolerance = mean( mean( abs(temp - temp1) ) );
       if( errorTolerance > 0.000002 )
       U3 = U0;
       end
     end
   % ------------------
     if( ~isreal(U4) )
     U4 = U0;
     else
     temp = U4*U4';       
     t1 = ones(1,nDOF);
     temp1 = diag(t1);                       % this is the identity matrix
     errorTolerance = mean( mean( abs(temp - temp1) ) );
       if( errorTolerance > 0.000002 )
       U4 = U0;
       end
     end
   % ------------------
     if( ~isreal(U5) )
     U5 = U0;
     else
     temp = U5*U5';       
     t1 = ones(1,nDOF);
     temp1 = diag(t1);                       % this is the identity matrix
     errorTolerance = mean( mean( abs(temp - temp1) ) );
       if( errorTolerance > 0.000002 )
       U5 = U0;
       end
     end
   % ------------------
     if( ~isreal(U6) )
     U6 = U0;
     else
     temp = U6*U6';       
     t1 = ones(1,nDOF);
     temp1 = diag(t1);                       % this is the identity matrix
     errorTolerance = mean( mean( abs(temp - temp1) ) );
       if( errorTolerance > 0.000002 )
       U6 = U0;
       end
     end
   % ------------------
     if( ~isreal(U7) )
     U7 = U0;
     else
     temp = U7*U7';       
     t1 = ones(1,nDOF);
     temp1 = diag(t1);                       % this is the identity matrix
     errorTolerance = mean( mean( abs(temp - temp1) ) );
       if( errorTolerance > 0.000002 )
       U7 = U0;
       end
     end
   % ------------------
     if( ~isreal(U8) )
     U8 = U0;
     else
     temp = U8*U8';       
     t1 = ones(1,nDOF);
     temp1 = diag(t1);                       % this is the identity matrix
     errorTolerance = mean( mean( abs(temp - temp1) ) );
       if( errorTolerance > 0.000002 )
       U8 = U0;
       end
     end
% ++++++++++++++++++++++++++++++++++++++++++++++++ do not use: Problematic
% ----------------------------------- try 2 common spatial pattern methods
% % % % % % % % % %    diagLambda = diag(D1);
% % % % % % % % % %    lowLambda = max( diagLambda )/10000000;
% % % % % % % % % %    D1 = diag( max(lowLambda,diagLambda) );
% % % % % % % % % %    diagLambda = diag(D2);
% % % % % % % % % %    lowLambda = max( diagLambda )/10000000;
% % % % % % % % % %    D2 = diag( max(lowLambda,diagLambda) );
% % % % % % % % % %    diagLambda = diag(D4);
% % % % % % % % % %    lowLambda = max( diagLambda )/10000000;
% % % % % % % % % %    D4 = diag( max(lowLambda,diagLambda) );
% % % % % % % % % %    diagLambda = diag(D5);
% % % % % % % % % %    lowLambda = max( diagLambda )/10000000;
% % % % % % % % % %    D5 = diag( max(lowLambda,diagLambda) );   
% % % % % % % % % %    R0 = U1*D1*U1';
% % % % % % % % % %    R1 = U2*D2*U2';
% % % % % % % % % %    T0 = U4*D4*U4';
% % % % % % % % % %    T1 = U5*D5*U5';
% % % % % % % % % %    CR = R0 + R1;
% % % % % % % % % %    CT = T0 + T1;
% % % % % % % % % %    [WR,sigR] = eig(CR);
% % % % % % % % % %    sigSqrtInv = diag( 1./sqrt( diag(sigR) ) );
% % % % % % % % % %    PR = sigSqrtInv*WR';
% % % % % % % % % %    [WT,sigT] = eig(CT);
% % % % % % % % % %    sigSqrtInv = diag( 1./sqrt( diag(sigT) ) );
% % % % % % % % % %    PT = sigSqrtInv*WT';
% % % % % % % % % %    RR0 = PR*R0*PR';
% % % % % % % % % %    %RR1 = PR*R1*PR';   
% % % % % % % % % %    TT0 = PT*T0*PT';
% % % % % % % % % %    %TT1 = PT*T1*PT';
% % % % % % % % % %    [U9,~] = eig(RR0);
% % % % % % % % % %    [U0,~] = eig(TT0);
% ---------------------------------------
%    isreal(U9)
%    isreal(U0)
% %    U9'*RR0*U9
% %    pause
% %    U9'*RR1*U9
% %    pause
% %    U0'*TT0*U0
% %    pause
% %    U0'*TT1*U0
% %    pause
% % % %    U0'*U0 
% % % %    pause
% % % %    transpose( U0'*U0 ) 
% % % %    pause
% % % %    U0*U0'
% % %    pause
% % %    U9'*U9 
% % %    pause
% % %    U9*U9'
% % %    pause
%%                             compare different initialization conditions
%vT = voteThreshold; 
gvSPLOC.sVS = set_sVS;                             % must reset every time
splocResults1 = getBasisVecSpectrum(U1,baseFname,trait1,trait0,vT); % try1
                                  % ^^-> nonfunctional systm PCA basis set
gvSPLOC.sVS = set_sVS;                             % must reset every time
splocResults2 = getBasisVecSpectrum(U2,baseFname,trait1,trait0,vT); % try2
                                  % ^^---> functional system PCA basis set
gvSPLOC.sVS = set_sVS;                             % must reset every time
splocResults3 = getBasisVecSpectrum(U3,baseFname,trait1,trait0,vT); % try3
                                  % ^^------> pooled systems PCA basis set
gvSPLOC.sVS = set_sVS;                             % must reset every time
splocResults4 = getBasisVecSpectrum(U4,baseFname,trait1,trait0,vT);
                                  % ^^------------------------------> try4
gvSPLOC.sVS = set_sVS;                             % must reset every time
splocResults5 = getBasisVecSpectrum(U5,baseFname,trait1,trait0,vT);
                                  % ^^------------------------------> try5
gvSPLOC.sVS = set_sVS;                             % must reset every time
splocResults6 = getBasisVecSpectrum(U6,baseFname,trait1,trait0,vT);
                                  % ^^------------------------------> try6
gvSPLOC.sVS = set_sVS;                             % must reset every time
splocResults7 = getBasisVecSpectrum(U7,baseFname,trait1,trait0,vT);
                                  % ^^------------------------------> try7
gvSPLOC.sVS = set_sVS;                             % must reset every time
splocResults8 = getBasisVecSpectrum(U8,baseFname,trait1,trait0,vT);
                                  % ^^------------------------------> try8
gvSPLOC.sVS = set_sVS;                             % must reset every time
% ----------------------------------------------------------- tricky cases
% % % REMARK: spatial pattern methods are often unstable
% %    try
% %    splocResults9 = getBasisVecSpectrum(U9,baseFname,trait1,trait0,vT);
% %    passed09 = true;                  % ^^-----------------------> try9
% %    catch
% %    passed09 = false;
% %    end
% % gvSPLOC.sVS = set_sVS;                         % must reset every time
% %    try
% %    splocResults0 = getBasisVecSpectrum(U0,baseFname,trait1,trait0,vT);
% %    passed10 = true;                  % ^^-----------------------> try0
% %    catch
% %    passed10 = false;      
% %    end                              
% ------------------------------------------------------------------------
   if( pursuitType == 1 )
   eee1 = sum(splocResults1.EEVd);
   eee2 = sum(splocResults2.EEVd);
   eee3 = sum(splocResults3.EEVd);
   eee4 = sum(splocResults4.EEVd);
   eee5 = sum(splocResults5.EEVd);
   eee6 = sum(splocResults6.EEVd);
   eee7 = sum(splocResults7.EEVd);
   eee8 = sum(splocResults8.EEVd);
%       if( passed09 )
%       eee9 = sum(splocResults9.EEVd);
%       else
%       eee9 = -1;
%       end
%       if( passed10 )
%       eee0 = sum(splocResults0.EEVd);
%       else
%       eee0 = -1;
%       end
   elseif( pursuitType == -1 )
   eee1 = sum(splocResults1.EEVi);
   eee2 = sum(splocResults2.EEVi);
   eee3 = sum(splocResults3.EEVi);
   eee4 = sum(splocResults4.EEVi);
   eee5 = sum(splocResults5.EEVi);
   eee6 = sum(splocResults6.EEVi);
   eee7 = sum(splocResults7.EEVi);
   eee8 = sum(splocResults8.EEVi);
%       if( passed09 )
%       eee9 = sum(splocResults9.EEVi);
%       else
%       eee9 = -1;
%       end
%       if( passed10 )
%       eee0 = sum(splocResults0.EEVi);
%       else
%       eee0 = -1;
%       end
   else
   eee1 = splocResults1.efficacy;
   eee2 = splocResults2.efficacy;
   eee3 = splocResults3.efficacy;
   eee4 = splocResults4.efficacy;
   eee5 = splocResults5.efficacy;
   eee6 = splocResults6.efficacy;
   eee7 = splocResults7.efficacy;
   eee8 = splocResults8.efficacy;
%       if( passed09 )
%       eee9 = splocResults9.efficacy;
%       else
%       eee9 = -1;
%       end
%       if( passed10 )
%       eee0 = splocResults0.efficacy;
%       else
%       eee0 = -1;
%       end
   end 
% ------------------------------------------------------------------------
%    if( passed09 )
%    Dd09 = splocResults9.Dd;
%    Di09 = splocResults9.Di;
%    else
%    Dd09 = 0;
%    Di09 = 0;
%    end
%    if( passed10 )
%    Dd10 = splocResults0.Dd;
%    Di10 = splocResults0.Di;
%    else
%    Dd10 = 0;
%    Di10 = 0;
%    end  
% % %    disp( pursuitType );
% % %    disp( [eee1,eee2,eee3,eee4,eee5,eee6,eee7,eee8,eee9,eee0] );
% % %    disp( [splocResults1.Dd, splocResults2.Dd, splocResults3.Dd, ...
% % %           splocResults4.Dd, splocResults5.Dd, splocResults6.Dd, ...
% % %           splocResults7.Dd, splocResults8.Dd, Dd09,  Dd10] );
% % %    disp( [splocResults1.Di, splocResults2.Di, splocResults3.Di, ...
% % %           splocResults4.Di, splocResults5.Di, splocResults6.Di, ...
% % %           splocResults7.Di, splocResults8.Di, Di09,  Di10] );


%    disp( pursuitType );
%    disp( [eee1,eee2,eee3,eee4,eee5,eee6,eee7,eee8] );
%    disp( [splocResults1.Dd, splocResults2.Dd, splocResults3.Dd, ...
%           splocResults4.Dd, splocResults5.Dd, splocResults6.Dd, ...
%           splocResults7.Dd, splocResults8.Dd] );
%    disp( [splocResults1.Di, splocResults2.Di, splocResults3.Di, ...
%           splocResults4.Di, splocResults5.Di, splocResults6.Di, ...
%           splocResults7.Di, splocResults8.Di] );
      
% ------------------------------------------------------------------------
%[~,jsort] = sort( [eee1,eee2,eee3,eee4,eee5,eee6,eee7,eee8,eee9,eee0]);
[~,jsort] = sort( [eee1,eee2,eee3,eee4,eee5,eee6,eee7,eee8]);
breakDegeneracy = zeros(1,10);
breakDegeneracy( jsort(1) )  = 0.05;
breakDegeneracy( jsort(2) )  = 0.10;
breakDegeneracy( jsort(3) )  = 0.15;
breakDegeneracy( jsort(4) )  = 0.20;
breakDegeneracy( jsort(5) )  = 0.25;
breakDegeneracy( jsort(6) )  = 0.30;
breakDegeneracy( jsort(7) )  = 0.35;
breakDegeneracy( jsort(8) )  = 0.40;
% % % breakDegeneracy( jsort(9) )  = 0.45;
% % % breakDegeneracy( jsort(10) ) = 0.50;
   if( pursuitType == -1 )                        % <= focus on i-subspace
   eee1 = splocResults1.Di + breakDegeneracy(1);
   eee2 = splocResults2.Di + breakDegeneracy(2);
   eee3 = splocResults3.Di + breakDegeneracy(3);
   eee4 = splocResults4.Di + breakDegeneracy(4);
   eee5 = splocResults5.Di + breakDegeneracy(5);
   eee6 = splocResults6.Di + breakDegeneracy(6);
   eee7 = splocResults7.Di + breakDegeneracy(7);
   eee8 = splocResults8.Di + breakDegeneracy(8);
% % %    eee9 = Di09             + breakDegeneracy(9);
% % %    eee0 = Di10             + breakDegeneracy(10);
   else                                           % <= focus on d-subspace
   eee1 = splocResults1.Dd + breakDegeneracy(1);
   eee2 = splocResults2.Dd + breakDegeneracy(2);
   eee3 = splocResults3.Dd + breakDegeneracy(3);
   eee4 = splocResults4.Dd + breakDegeneracy(4);
   eee5 = splocResults5.Dd + breakDegeneracy(5);
   eee6 = splocResults6.Dd + breakDegeneracy(6);
   eee7 = splocResults7.Dd + breakDegeneracy(7);
   eee8 = splocResults8.Dd + breakDegeneracy(8);
% % %    eee9 = Dd09             + breakDegeneracy(9);
% % %    eee0 = Dd10             + breakDegeneracy(10);
   end
%[~,jsort] = max( [eee1,eee2,eee3,eee4,eee5,eee6,eee7,eee8,eee9,eee0] );
[~,jsort] = max( [eee1,eee2,eee3,eee4,eee5,eee6,eee7,eee8] );
jSelect = jsort(1);
      switch jSelect
          case 1     % nonfunctional PCA basis set   try 1
          U = U1;
          msgInitial_U = 'initial basis set given by: PCA nonfunctional';
          %disp('case U = U1');
          case 2     % functional PCA basis set   try 2
          U = U2;
          msgInitial_U = 'initial basis set given by: PCA functional';
          %disp('case U = U2');
          case 3     % pooled average PCA basis set  try 3
          U = U3;
          msgInitial_U = 'initial basis set given by: PCA pooled';
          %disp('case U = U3');
          case 4     % generalized eigenvalue try 4
          U = U4;
          msgInitial_U = 'initial basis set given by: try4';
          %disp('case U = U4');
          case 5     % generalized eigenvalue try 5
          U = U5;
          msgInitial_U = 'initial basis set given by: try5';
          %disp('case U = U5');
          case 6     % generalized eigenvalue try 6
          U = U6;
          msgInitial_U = 'initial basis set given by: try6';
          %disp('case U = U6');
          case 7     % generalized eigenvalue try 7
          U = U7;
          msgInitial_U = 'initial basis set given by: try7';
          %disp('case U = U7');
          case 8     % generalized eigenvalue try 8
          U = U8;
          msgInitial_U = 'initial basis set given by: try8';
          %disp('case U = U8');
% % %           case 9     % generalized eigenvalue try 9
% % %           U = U9;
% % %           msgInitial_U = 'initial basis set given by: try9';
% % %           %disp('case U = U9');
% % %           case 10     % generalized eigenvalue try 0
% % %           U = U0;
% % %           msgInitial_U = 'initial basis set given by: try0';
% % %           %disp('case U = U0');
      end  
% ---------------------------- an un-need check if by-pass was implemented
% DO NOT REMOVE:  Orthogonality is required. 
   temp = U*U';       
   t1 = ones(1,nDOF);
   temp1 = diag(t1);                         % this is the identity matrix
   errorTolerance = mean( mean( abs(temp - temp1) ) );
      if( errorTolerance > 0.000002 )
      disp(['error level = ',num2str(errorTolerance)]);
      error('U basis set is not orthonormal within required precision');
      end
   else
   U = U0;                       % initial guess provided by user: Thanks!
   end  
% disp(msgInitial_U);
% % error('stop here for now');
%%                                             maximize objective function
% ------------------- initialize arrays for projected averages & variances
ave0 = zeros(n0,2);   % only 2D projections are needed, not nDOF dimension
var0 = zeros(n0,2);
ave1 = zeros(n1,2);
var1 = zeros(n1,2);
% --------------------- initialize cell arrays for trait projections in 2D
m0 = cell(1,n0);
q0 = cell(1,n0);
m1 = cell(1,n1);
q1 = cell(1,n1);
%%                                         initialize counters of interest
% -------------------------------------- tracking rotation characteristics
debugEfficacyVecJ = zeros(Na,2);       % efficacy score for vector J = 1,2
debugQualitiesPwJ = zeros(Na,2);       % quality factor for vector J = 1,2
debugConsensusPwJ = zeros(Na,2);      % consensus power for vector J = 1,2
debugSelectionPwJ = zeros(Na,2);      % selection power for vector J = 1,2
debugSubspaceType = zeros(Na,2);                % 1,2,3 => d,i,u subspaces
count_j1 = zeros(1,nDOF);   % # of times j1 is 1st vector in a vector pair
% --------------------------------------------- one time initialize arrays
LoverRide = false(1,nDOF);
efficacyArray = ones(1,nDOF);                  % controls overall efficacy
voterFraction = zeros(1,nDOF);   % each basis vector gets a consensus vote
got2DEfficacy = ones(1,2);                  % initialize 2D Efficacy score
aDevLnScore2D = zeros(1,2);         % = abs(LnScoreFunction - lnMidScoreX)
voter2DFrcton = zeros(1,2);
modeSelectionTime = zeros(1,nDOF);
modeIncrementTime = ones(1,nDOF);
passoverCount = zeros(nDOF);   % controls converged vector pair scheduling
% ------------------------------------------- projection pursuit iteration
runningAveScoreIncrement = 1.0;
convergenceRatio = 1;      %=> requires some effort to get below threshold
efficacySwapPerc = 1;
waitTime = 1.0e-10;                      % a small number should work well
pFmax = pFmaxStartValue;
nVectorPairsLess1 = nVectorPairs - 1;
DebugON = (verbosity == 3);
initialScore = -1.0; 
countRotations = 0;
consensusVote = 0;
dampLR = 0;
countSucc = 0;                        % count number of non-zero rotations
countSkip = 0;               % count number of times rotations are skipped
waitStep = 1;
CayleyCPU = 0;          
xCayley = 0;              % count relative to reference rotation amplitude
nextCayley = nDOF;
checkPointWrite = 0;
vT = voteThreshold;
refEfficacy = 0;
nEfficacyCount = 0;
flag_converged = -1;
skipRate = 0;
nFailed = 0;
nResets_final = nTimes;
% histSucc = zeros(1,501);                                         % check
% --------------------------------- reset nTimes is an iterative procedure
for nResets=1:nTimes
gvSPLOC.sVS = set_sVS;                             % must reset every time
gotBVS = getBasisVecSpectrum(U,'internal',trait1,trait0,vT);
newEfficacy = gotBVS.efficacy;
% ------------------------------------------------------ check if finished
   if( nResets > 1 )                             % must do at least 1 run!
   tempNumerator = abs(newEfficacy - refEfficacy);
   tempDenominator = abs(refEfficacy) + 0.00001;
   temp = tempNumerator/tempDenominator;
% ------------------------------------------------ instantaneous condition 
      if( temp < 0.01 )
      flag_converged = 1;
      nResets_final = nResets - 1;
      break;
      end
% ------------------------------------------------ count down: 12%, 8%, 4%
      if( temp < (0.15 - 0.05*nEfficacyCount) )  
      nEfficacyCount = nEfficacyCount + 1;
      end
% ----------------------------------------------------------- exceed count
      if( nEfficacyCount > 2 )
      flag_converged = 1;
      nResets_final = nResets - 1;
      break;
      end
% -------------------------------------------------- prevent infinite loop
      if( countRotations > maxRotations )                % emergency break
      break;
      end
   end
% --------------------------------------------------------------- continue
refEfficacy = newEfficacy;
oldEfficacy = refEfficacy;
% ----------------------------------------
   if( pursuitType == 1 )
   efficacyScore = sum(gotBVS.EEVd);
   elseif( pursuitType == -1 )
   efficacyScore = sum(gotBVS.EEVi);
   else
   efficacyScore = sum(gotBVS.EEVd + gotBVS.EEVi);
   end
% --------------------------------- control variables for Cayley rotations
maxDrotateFraction = upperDrotateFraction;
randomMatrixAmplitude = ref_randomMatrixAmplitude;
% ---------------------------------- control variables for vector rotation
% REMARK: The efficacyArray is a conglomerate evaluation that informs how
%         effective each projection is in terms of selection, consensus 
%         and clustering qualities. 
efficacyArray(:) = initialScore;
voterFraction(:) = 0;
got2DEfficacy(:) = initialScore;
aDevLnScore2D(:) = 0;
voter2DFrcton(:) = 0;
oldCountRotations = 0;
oldCountSkip = 0;
% ------------------------------------------ initialize relevant variables
consensusVote = 0;                             % => initially 0% consensus
old_vScore = 0;                              % used to monitor convergence
penaltyFactor = pFmax*ones(nDOF);     % schedules low quality vector pairs
% REMARK:       ^^^^^^^^^^^^^^^^--> replace various pF values with mean pF
rotateProbability = 0.125 + 0.875*rand(nDOF).^(0.6*log(nDOF) - 0.4); %nice
rotateProbability = min( rotateProbability , rotateProbability' );
skipFraction = 0;
skipRate = 0;
nFailed = 0;
modeSelectionTime(:) = 0;
modeIncrementTime(:) = 0;
passoverCount(:) = 0;
showWaitTimePlot = 0;                 % counter to show wait time bar plot
mVPc = 0;                                   % # of vector pairs considered
minRotations = minRotations0 + countRotations;          % w.r.t last count
             % ^^^^^^^^^^^^^-------> check at least this many vector pairs
nLessons = 0;                  % # of lessons used to adjust learning rate
tryAgain = true;
   while( tryAgain )            % => take indefinite # of lessons to learn
   ppp = nLessons/targetLesson;
   dampLR = min(1,0.1 + 0.9*ppp);        % cap => no accelerated learning!
% REMARK: Upon each restart, the learning rate is slower. This has the
% effect of slowing developing memory without jumping to conclusions.  
  %disp( [nLessons,dampLR] );
%%                                      set control flags for optimization
      if( pursuitType < -0.1 )
      rrr = gotBVS.Di/nDOF1;
      mFpow = (1 - rrr)^p4t;
      adaptive_sVS = set_sVS*10^mFpow;
      else
      rrr = gotBVS.Dd/nDOF1;
      mFpow = (1 - rrr)^p4t;
      adaptive_sVS = set_sVS/10^mFpow;
      end
   sVS = adaptive_sVS;
%    disp( num2str([gotBVS.Dd,gotBVS.Du,gotBVS.Di,adaptive_sVS]) );
%    pause(0.3)
   gvSPLOC.sVS = adaptive_sVS;                     % must reset every time
   [gotBVS,mapU2S,mapS2U] = ...
                       getBasisVecSpectrum(U,'internal',trait1,trait0,vT);
   %^^^^^^^---------> current characteristics of the basis vector spectrum
% {
% ------------------------------------- check criteria for Cayley rotation
   addCayleyCPU = cputime; %=> includes all prep steps even if no rotation
% REMARK 1: There are two reasons for appying Cayley rotations. 
%     1) The convergence ratio is large, allowing for time to destablize
%        the undetermined subspace so as to explore possibilities at the
%        early stage of projection pursuit. 
%     2) A certain number of Jacobi rotations must be made before another
%        Cayley rotation is made, otherwise the system does not have time
%        to relax from the random perturbations. 
% REMARK 2: The average deviation angle made by a Cayley rotation should 
%        increase up to a maximum value as long as large gains in efficacy
%        are being made. Of course, the average deviation should not be 
%        extreme, although extreme is okay if large enough gains are being
%        made. But it was noted that being too aggressive with high degree
%        of rotations can confuse the Jacobi rotations. Confusion means 
%        that the convergence criteria stalls, and the process will run
%        "forever", analogous to waiting for water to turn into ice at
%        30 degrees Celsuis at 1 ATM. Although POSSIBLE for H2O to take
%        the form of ice at 30 degrees Celsuis at 1 ATM presure, it is 
%        NEVER going to happen thanks to entropy and massive fluctuations
%        that occur at relatively high temperature. Likewise here, it is 
%        important not to apply too large of random fluctuations. This is
%        accomplished by repeating moderate to small size rotations in 
%        order to create a diffusive process. The Jacobi-rotations select
%        the most fit vector directions, and thus sets up a bias. Even for
%        small degree rotations, the diffusive process works very well, 
%        allowing vector space directions to be explored very efficiently.
%        Changing the amplitude of the Cayley rotation matrix changes the
%        effective time scale of this diffusive process. Relative to a 
%        reference rotation amplitude the time scale of the diffusive 
%        process can be reduced or increased by a factor of 1/sqrt(5) or
%        sqrt(5) respectively, yielding an effective range of 5 to set the
%        rate of the Cayley randomization method of u-basis vectors.
%        However, upon convergence the typical gain in efficacy for a spin
%        decreases. As convergence is approached, randomization of u-modes
%        will potentially substantially slow down convergence, or worse, 
%        prevent convergence. This process is analogous to simulated
%        annealing, where it is always the case that as the ground state
%        is aproached, the thermal fluctuations must go to zero. That is,
%        T must change from high to 0 value, which controls fluctuations
%        from high to low (zero). This processess is controlled here in 
%        a fully automated data driven adaptive approach.
% ------------------------------------------------------------------------
     if( convergenceRatio > endGameRatio )
        if( countRotations > nextCayley )
       % -------------------------------------------------- efficacy level
       % REMARK: Argued that efficacy = 0 for u-modes, but note precursor
       %         efficacy is being tracked here, which is not zero. This
       %         information is being used in an informative way. 
        pRotate1 = max(0, 1 - sqrt(0.0000001 + gotBVS.EEVd/0.1) );   
        pRotate2 = max(0, 1 - sqrt(0.0000001 + gotBVS.EEVi/0.1) );   
        pRotate = pRotate1.*pRotate2;
       % ---------------------------------------- statistical significance
        pRotate1 = max(0, 1 - (gotBVS.CEVd/0.5).^2 );
        pRotate2 = max(0, 1 - (gotBVS.CEVi/0.5).^2 );
        pRotate = pRotate.*pRotate1.*pRotate2;
       % ------------------------------------------------- signal to noise
           if( pursuitType == 1 )
           tempSEVd = gotBVS.SEVd - middleScore;
           pRotate1 = max(0, 1 - tempSEVd/0.25);
           elseif( pursuitType == -1 )
           tempSEVi = middleScore - gotBVS.SEVi;
           pRotate1 = max(0, 1 - tempSEVi/0.25);
           else
           tempSEVd = gotBVS.SEVd - middleScore;
           tempSEVi = middleScore - gotBVS.SEVi;
           absDiffSEV = abs( tempSEVd - tempSEVi );
           pRotate1 = 2*max(0, 1 - ( absDiffSEV/0.25).^2 );
           end       %^-------------> greater than 1 will not be a problem
        pRotate = maxDrotateFraction*pRotate.*pRotate1;
       % ----------------------------- determine which u-vectors to rotate
        randPick = rand(1,nDOF); 
        Su = ( gotBVS.Cind == 0 );              % => must be in u-subspace
        Su = and( Su, (randPick < pRotate) );  % random rotation selection
        Drotate = sum(Su);     %=> # of vectors selected for a random walk
% ----------------------------------------------------- for debugging only
% %         tempSEVd = gotBVS.SEVd - middleScore;
% %         tempSEVi = middleScore - gotBVS.SEVi;
% %         absDiffSEV = abs( tempSEVd - tempSEVi );
% %         disp( [gotBVS.EEVd; gotBVS.CEVd; tempSEVd; ...
% %                gotBVS.EEVi; gotBVS.CEVi; tempSEVi; ...
% %                absDiffSEV;      pRotate; randPick; Su] );
% %         disp( Drotate );
% %         pause
           if( Drotate > 2 )   % 3 is minimum # to make effort worth doing
           doCayleyRotation = true;
           nextCayley = nextCayley + nDOF;     % set relaxation delay time
           sumDrotate = sumDrotate + Drotate;
           Nrotate = Nrotate + 1;
           else
           Drotate = 0;      % not needed, but used for optional reporting
           doCayleyRotation = false;    % skip diffusion process this time
           end
        else
        doCayleyRotation = false;
        end 
     else
     doCayleyRotation = false;
     end
%%                          randomize vectors in the undetermined subspace
% --------------------------------------------- undue what was done before
     Lrandom = ( modeIncrementTime < 1.0e-20 );
     modeSelectionTime(Lrandom) = waitTime;
     modeIncrementTime(Lrandom) = 2*waitTime;
%%                          randomize vectors in the undetermined subspace
        if( doCayleyRotation )
% REMARK: all u-basis vectors have a degenerate efficacy value of 0. This
%         means we can spin them randomly. By considering only pairs of
%         u-basis vectors any random spin is as good as any other. This 
%         will not change the efficacy. It could change other metrics such
%         as quality, consensus and selection power, but efficacy is what
%         we are trying to maximimize, and undetermined basis vectors have
%         0-efficacy since they cannot be used in any useful way! As such,
%         this means that the directions of the u-basis vectors can be 
%         randomized using a random walk approach. In time, successive
%         random rotations will allow the u-space to be explored more
%         efficiently than the alternatve of spinning all possible pairs
%         of u-basis vectors -- one pair at a time. As such, this global
%         randomization by a single random Cayley rotation 
%         SEE: https://en.wikipedia.org/wiki/Cayley_transform     enhances
%         the probability of spinning a d-basis vector with a u-basis
%         vector or an i-basis vector with a u-basis vector, or a d-basis
%         vector with an i-basis vector. Since the u-basis vectors are in
%         a sea of randomness, this will allow the d&i-subspace to extract
%         information from the u-subspace abyss. The spinning of u-modes 
%         will not mix the subspaces d with u or i with u since all types
%         of subspaces are orthogonal. Nevertheless, due to spinning of a 
%         u-mode with either a d- or i-mode, the u-subspace changes within
%         the total vector space. In some sense, the entire u-subspace is
%         diffusing through the total vector space. Greater effectiveness
%         to extract information promotes faster convergence and greater
%         accuracy. This is the heristic for this section of the code. But
%         there is a problem. If most basis vectors have efficacy of zero,
%         excessive random spinning will prevent a signal to be found. The
%         indiscriminate spinning of vectors is like T->infinity thermal 
%         bath. Therefore, the number of random spins and the amplitude of
%         a typical random spin angle must be limited. This idea is like
%         having a signal + noise. If noise >> signal, there is no hope to
%         capture a signal since the vectors are going haywire! If the
%         signal > noise, where the noise does not do much in 1 iteration,
%         but when the noise is applied over long times it is seen that
%         the vector-directions of u-modes are "diffusing" as a "blob" in
%         a high dimensional space. This diffusion is slow, but it allows
%         for the possibility to extract the maximum amount of information
%         out of the u-subspace, and empty its information content as it
%         is sucked dry (as much as possible). Furthermore, observations
%         show that having too much randomness prevents  NUCLEATION  to
%         occur, meaning a d- or i-mode appears for the first time. Thus,
%         the number of vectors to randomize over is also contrained by 
%         enforcing a maximum cap so there are vectors that can nucleate.
% REMARK: In practice the typical deviation angle of a random rotation
%         should be small (say ~8 degrees). This will not allow the plane
%         of a previous vector pair to be much different than what it was
%         before the rotation between a d-mode with a u-mode or a i-mode
%         with a u-mode. When this vector pair spins again, a slightly
%         different answer will be obtained, but the answer will be a 
%         small perturbation from the original value. Over many such spins
%         the vector directions over u-modes will diffuse, and the d- and
%         i-modes will wiggle in such a way as to maximize the efficacy.
% ------------------------------------------------------------------------
%         tempCheck1 = gotBVS.efficacy;                 % debug check only
        UmodeIndices = mapS2U(Su);
        Uundetermined = U(:,UmodeIndices);
% ------------------------------------------------ begin a Cayley rotation
        xRotation = randomMatrixAmplitude/ref_randomMatrixAmplitude; 
        xCayley = xCayley + sqrt(xRotation)*Drotate; %<= diffusive process
                          % ^^^^^^^^^^^^^^^--> weight by angular dev. size
        D0 = min(Drotate,10);
% ---------------------------------------- mathematics for Cayley rotation
        indexREF = randperm(Drotate);
        jLower = 1;
        jUpper = D0;
        nRmax = Drotate - D0;
           for nR=0:D0:nRmax
           index = indexREF(jLower:jUpper);
           V0 = Uundetermined(:,index);
% ---------------------------------------------- construct rotation matrix
           Id = diag( ones(1,D0) );
           S0 = randomMatrixAmplitude*rand(D0);      % 0.05 => small angle
%          S0 = 0*S0;                                   % debug check only
           S = S0 - S0';                       % => S is now a skew matrix
           R = (Id - S)/(Id + S);    % => random rotation matrix in D0 DIM
           Uundetermined(:,index) = V0*R;
           jLower = jLower + D0;
           jUpper = jUpper + D0;
           end
        U(:,UmodeIndices) = Uundetermined;
% -------------------------------- modify waitTime for this set of u-modes
        modeSelectionTime(UmodeIndices) = 1.4*waitTime;
        modeIncrementTime(UmodeIndices) = 1.0e-25;  % needs to be a tiny #
        gvSPLOC.sVS = adaptive_sVS;                % must reset every time
        [gotBVS,mapU2S,mapS2U] = ...
                       getBasisVecSpectrum(U,'internal',trait1,trait0,vT);
        %^^^^^--------> must update in case a u-mode leaves the u-subspace
        % ----------------------------------------------------------------
%         tempCheck2 = gotBVS.efficacy;                 % debug check only
%            if( abs(tempCheck2 - tempCheck1) > 1.0e-9 )
%            error('IMPOSSIBLE error occurred!');
%            end
        % ----------------------------------------------------------------
        Su = ( gotBVS.Cind == 0 ); 
           if( sum(~Su) > 0 )
           notUmodeIndices = mapS2U(~Su);
           modeSelectionTime(notUmodeIndices) = waitTime;
           modeIncrementTime(notUmodeIndices) = 2*waitTime;
           end
% ------------------------------------- adaptively set rotation parameters
% REMARK: pp0 defines a reference point, where 50% success rate sets pp0=1
%         which is found to be a reasonably good value for pp0 if it were
%         a constant. When the success rate dips below 50% randomization
%         in the degenerate undetermined subspace is applied more strongly
%         via more randomization, and less randomization is applied as
%         the success rate of Jacobi spins increase. This process attempts
%         to explore options better, but it also makes the solutions that
%         are obtained more random. 
        pp0 = (2*countSucc/countRotations)^2;
        pp = pp0 - min(pp0 + 10,scoreChange);
        adaptiveFactor = (0.985)^pp;
        maxDrotateFraction = adaptiveFactor*maxDrotateFraction;
        maxDrotateFraction = max(lowerDrotateFraction,maxDrotateFraction);
        maxDrotateFraction = min(upperDrotateFraction,maxDrotateFraction);
        randomMatrixAmplitude = adaptiveFactor*randomMatrixAmplitude;
        randomMatrixAmplitude = max(minRMA,randomMatrixAmplitude);
        randomMatrixAmplitude = min(maxRMA,randomMatrixAmplitude);
        %disp([pp,scoreChange,randomMatrixAmplitude,maxDrotateFraction]);
        %pause(0.5)
        end 
     addCayleyCPU = cputime - addCayleyCPU;
     CayleyCPU = CayleyCPU + addCayleyCPU;    % includes prep steps, which
     %                                  may include not doing any rotation
%}  
%%                                        select importance weight factors
   selectPRBcase = rand;                        % 12
      if( selectPRBcase < cPROBii )             % ii pairs are most likely
      v1 = 2;
      v2 = 2;
      elseif( selectPRBcase < cPROBid )         % id pairs are most likely
      v1 = 2;
      v2 = 1;
      elseif( selectPRBcase < cPROBdi )         % dd pairs are most likely
      v1 = 1;
      v2 = 2;
      else               % cPROBdd = 1          % di pairs are most likely
      v1 = 1;
      v2 = 1;
      end  
   %disp( [mapU2S; mapS2U] );
%   disp( mapU2S );
    mapU2S(LoverRide) = overRideRank;
%   disp([LoverRide ; mapU2S]);
%   pause 
%%                                     work with one vector-pair at a time
      if( mVPc == 0 )                            % mVPc = 0 => a new start
      max2DscoreDiff = -1;  % => flags no change in score since last start
      end
   j1 = j2Pair(nDOF);   % j1 is initially set to j2 that cannot be reached
   count_j1(j1) = count_j1(j1) + 1;      % for non-essential tracking only
   prb1 = vPROB(mapU2S(j1),v1);
% -------------------------------------------------------- secure a quorum
% REMARK: Found the gremlin: This section that secures a quorum is needed
% to prevent anemic performance simply because nothing was being done.
   old_rotateProbability = rotateProbability(:,j1);
   new_rotateProbability = old_rotateProbability;
   prb = prb1*vPROB(mapU2S(1:nDOF),v2);
   nPresent = 0;
      while( nPresent < quorum )
      new_rotateProbability = growFactor*new_rotateProbability;
      Lselect = ( new_rotateProbability > prb );
      nPresent = sum(Lselect);
% ------------------------------------------------------------------ check
%       disp( [nPresent,quorum] );
%          if( nPresent < 1 )
%          beep;
%          pause
%          end
      end
   rotateProbability(:,j1) = growFactor*old_rotateProbability;
   rotateProbability(Lselect,j1) = new_rotateProbability(Lselect);
   new_rotateProbability = rotateProbability(:,j1);
   rotateProbability(j1,:) = new_rotateProbability'; 
% --------------------------------------------------- select a 2D subspace
   scoreChange = 0;           % accumulate total change in score per sweep
   countSpin = 1.0e-25;
     for iVP=1:nVectorPairs        % j1 cannot reach last j2 value in list
     j2 = j2Pair(iVP);        % note that it is impossible to have j1 = j2
     %===================================================================|
     %prb2 = vPROB(mapU2S(j2),v2);  % NO NEED TO DO THIS: precalculated! |
     %prbTEST = prb1*prb2;          % NO NEED TO DO THIS: look at quorum |
     %       if( abs( prbTEST - prb(j2) ) > 0.00001 )                  % |
     %       beep                                                      % |
     %       disp( [prbTest, prb(j2)] );                               % |
     %       end                                                       % |
     % This box explains the relationship between the vectorized and the |
     % nonvectorized versions. The reason for implementing the quorum is |
     % to prevent senseless skips that eat up cpu time ineffectively. But|
     % simply shifting all probabilities up with as many growFactor      |
     % multiplications destroys the correlations. So, only the ones that |
     % pop above the thresholds are accelerated, and the vast majority of|
     % vector pairs remain the same. All this conditioning is done above |
     % in the quorum section, and can be done with vectorized code using |
     % the MATLAB way.                       ==> Nice!                   |
     %===================================================================|
      if( rotateProbability(j2,j1) < prb(j2) )  % skip vector pair (j1,j2)
      %prob2rotate = rotateProbability(j1,j2);      precalulated in quorum
      %prob2rotate = growFactor*prob2rotate;        precalulated in quorum
      %rotateProbability(j1,j2) = prob2rotate;      precalulated in quorum
      %rotateProbability(j2,j1) = prob2rotate;      precalulated in quorum
      modeSelectionTime(j2) = modeSelectionTime(j2) + waitTime;
      countSkip = countSkip + 1;
      else                                   % => spin vector pair (j1,j2)
      nLessons = nLessons + 1;                 % each rotation is a lesson
      countRotations = countRotations + 1;
      ShowSpin = (nDOF*floor(countRotations/nDOF) == countRotations);
      plotRotationInfo = and(DebugON,ShowSpin);   
%           if( mapU2S(j1) == 1 )
%           disp( mapU2S );    
%           plotRotationInfo = true;
%           disp(['sVS = ',num2str(sVS)]);
%           pause
%           end
      mVPc = mVPc + 1;
      modeList = [j1,j2];
      PR = U(:,modeList);      % define projection operator in 2D subspace
      PL = PR';      % PR is a nDOF x 2 matrix  &  PL is a 2 x nDOF matrix
      %                              REMARK: PL*PR = 2 x 2 identity matrix
      % ---------------------------- apply projection on M0, Q0, M1 and Q1
         for k0=1:n0
         m0{k0} = PL*M0{k0};
         q0{k0} = PL*Q0{k0}*PR;
         end
      % ------------------------
         for k1=1:n1
         m1{k1} = PL*M1{k1};
         q1{k1} = PL*Q1{k1}*PR;
         end
%%                                                optimize rotation matrix
% ------------------------------------------ initialize tracking registers
      oldScore2D = efficacyArray(j1) + efficacyArray(j2);  % value to beat
      bestScore2D = oldScore2D;     % => current best value is where it is
      bestVector2Dscore1 = efficacyArray(j1);
      bestVector2Dscore2 = efficacyArray(j2);
      bestAbsDevLnScore1 = 0;                   % set to 0 in lieu of span
      bestAbsDevLnScore2 = 0;                   % set to 0 in lieu of span
      bestVoterFraction1 = voterFraction(j1);
      bestVoterFraction2 = voterFraction(j2);
      bestInertIncrement = 1;
      iaBest = iaZero;                      % => always start at 0 degrees
%%                                             apply numerical Jacobi spin
         for zoomLevel=1:zoomLevelEnd                    % trial rotations
            if( zoomLevel > 2 )
            iaMIN = iaBest - shiftL(zoomLevel);
            iaMAX = iaBest + shiftL(zoomLevel);
            skip = skipL(zoomLevel);   
            else
            iaMIN = iaZero - shiftL(zoomLevel);
            iaMAX = iaZero + shiftL(zoomLevel);
            skip = skipL(zoomLevel);
            end
%          disp(iaMIN)
%          disp(iaMAX)
%          disp(skip)
%          pause
% ---------------------------------- apply hierarchical search in rotation
            for ia=iaMIN:skip:iaMAX       % specified range and resolution
            TR = TRcell{ia};             % 2D rotation matrix (right side)
            TL = TLcell{ia};    % TL = transpose(TR) = inv(TR) (left side)
            % Note:              similarity transform: A_new = TL*A_old*TR
% -------------------------- reference to machine learning jargon: COMMENT
% REMARK: Within the selected 2D subspace, there are two projections that
% are made. Using the "mean" vector, and "covariance" matrix, the exact
% sample mean & variance is rapidly calculated along projected directions.
% Note that  ave0(:,1) = 1st moment FEATURE (as mean) for projection 1
%            ave0(:,2) = 1st moment FEATURE (as mean) for projection 2
%            var0(:,1) = 2nd moment FEATURE (as variance) for projection 1
%            var0(:,2) = 2nd moment FEATURE (as variance) for projection 2
%               ^-^-------> for each of the nonfunctional systems
% Similarly: ave1(:,1), ave1(:,2), var1(:,1), var1(:,2) are 2 FEATURES per
% projection (e.g. 1 and 2) for functional systems. 
% In short, SPLOC is analyzing TWO FEATURES per DOF per system.
% Total number of features = 2xnDOF.
% The feature space is calculated using brute force (no kernel trick).
% The scoring function uses the Euclidean Pythagoras formula as selection
% power. It does not matter which feature of the two (mean or variance) is
% separating the data, as long as one is, or even by their combination via
% an Euclidean distance where standard deviation and mean define two axes
% in FEATURE space. 
% ------------------------------- calculate 2 features for two projections
               for k0=1:n0                           % nonfunctional cases
               vec2D = TL*m0{k0};
               var2D = TL*q0{k0}*TR;
               ave0(k0,1) = vec2D(1);                 % projected  average
               ave0(k0,2) = vec2D(2);                 % projected  average
               var0(k0,1) = max(var2D(1,1),add0);     % projected variance
               var0(k0,2) = max(var2D(2,2),add0);     % projected variance
               end
            % ------------------------------------------------------------
               for k1=1:n1                              % functional cases
               vec2D = TL*m1{k1};
               var2D = TL*q1{k1}*TR;
               ave1(k1,1) = vec2D(1);                 % projected  average
               ave1(k1,2) = vec2D(2);                 % projected  average
               var1(k1,1) = max(var2D(1,1),add0);     % projected variance
               var1(k1,2) = max(var2D(2,2),add0);     % projected variance
               end
%%                                       calculate efficacy per projection
               for j=1:2             % j is dummy vec-index in 2D subspace
               inertAddition = 1;                 % counts converged cases
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
               if( lnScoreFunctn10 > lnMidScoreX ) % discriminant subspace
% +++++++++++++++++++++++++++++++++++++++++++++++ calculate dQualityFactor
               vFraction = dVoteFraction10;
                 if( pursuitType ~= -1 )      % => requires dQualityFactor
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
                 % ----------------------------------------- for debugging
                 %{
              if( nDOF*floor(countRotations/nDOF) == countRotations )  
                 dMM0 = [dcp,vFrac,vFraction,consensus,dAdd,dampLR];
                 dMM1 = [fIndeterminable, fIndeterminable11, ...
                         fIndeterminable00,fImax];
                 dMM2 = [snr10Mean,rangeMean,gapMean,overlapMean, ...
                         spanMean0,spanMean1];
                 dMM3 = [snr10LnStDv,rangeLnStDv,gapLnStDv, ...
                         overlapLnStDv,spanLnStDv0,spanLnStDv1];
                 dMM4 = [cpMUg,cpMUo,cpSDg,cpSDo];
                 disp(['vector pair indices = (', ...
                      num2str(j1),',',num2str(j2),')']);
                 disp(['lnScoreFunction = ',num2str(lnScoreFunctn10)]);
                 disp(['              discriminant quality factor = ', ...
                       num2str(dQualityFactor)] );
                 disp(['dcp,vFrac,vFraction,consensus,dAdd,dampLR = ', ...
                       num2str(dMM0)] );
                 disp(['         fI10, fI11, fI00, max(fI00,fI11) = ', ...
                       num2str(dMM1)] );
                 disp(['  Mean: snr10,range,gap,ovlap,span0,span1 = ', ...
                       num2str(dMM2)] );  
                 disp(['LnStDv: snr10,range,gap,ovlap,span0,span1 = ', ...
                       num2str(dMM3)] );
                 disp(['               cpMUg, cpMUo, cpSDg, cpSDo = ', ...
                       num2str(dMM4)] );
                 disp( dividerLine );
              pause(0.25)
              end   
                 %}
                 else             % => looking for indifference modes only
                 dQualityFactor = -0.1;     % <= default anti-discriminant
                 end
               qualityFactor = dQualityFactor;
               span = (lnScoreFunctn10 - lnMidScoreX); 
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
               %inertAddition = inertAddition + 0;  NO addition for d-mode
               else % -----------------------------> indifference subspace
% +++++++++++++++++++++++++++++++++++++++++++++++ calculate iQualityFactor
               vFraction = iVoteFraction10;
                 if( pursuitType ~= 1 )       % => requires iQualityFactor
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
                 % ----------------------------------------- for debugging
                 %{ 
             %if( vFrac > 0.0 )
             %if( (icp > 0.0) && (nLessons > nDOF) )
             %if( (fIndeterminable11 > 0.5) || (fIndeterminable00 > 0.5) )
              if( nDOF*floor(countRotations/nDOF) == countRotations )
                 dMM0 = [icp,vFraction,consensus,iAdd,dampLR];
                 dMM1 = [fIndeterminable, fIndeterminable11, ...
                         fIndeterminable00];
                 dMM2 = [snr10Mean,rangeMean,gapMean,overlapMean, ...
                         spanMean0,spanMean1];
                 dMM3 = [snr10LnStDv,rangeLnStDv,gapLnStDv, ...
                         overlapLnStDv,spanLnStDv0,spanLnStDv1];
                 dMM4 = [cpMUg,cpMUo,cpSDg,cpSDo];
                 disp(['vector pair indices = (', ...
                       num2str(j1),',',num2str(j2),')']);
                 disp(['lnScoreFunction = ',num2str(lnScoreFunctn10)]);
                 disp(['              indifference quality factor = ', ...
                       num2str(iQualityFactor)] );
                 disp(['  icp, vFraction, consensus, iAdd, dampLR = ', ...
                       num2str(dMM0)] );
                 disp(['                         fI10, fI11, fI00 = ', ...
                       num2str(dMM1)] );
                 disp(['  Mean: snr10,range,gap,ovlap,span0,span1 = ', ...
                       num2str(dMM2)] );  
                 disp(['LnStDv: snr10,range,gap,ovlap,span0,span1 = ', ...
                       num2str(dMM3)] );
                 disp(['               cpMUg, cpMUo, cpSDg, cpSDo = ', ...
                       num2str(dMM4)] );
                 disp( dividerLine );
               pause(0.25)
               end
                 %}
                 else             % => looking for discriminant modes only
                 iQualityFactor = -0.1;     % <= default anti-indifference
                 end                
               qualityFactor = iQualityFactor;
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
%                disp(['proof in the pudding: sVS = ',num2str(sVS), ...
%                      '  set_sVS = ',num2str(set_sVS)]);
               vecScore = sVS*vecScore;           % => bias toward d-modes
              %      This ^^^--> factor is incorporated for these reasons:
              % Units on efficacy are arbitrary. Any positive constant can
              % multiply efficacy as a scale. When comparing d-modes to 
              % d-modes this scale factor is irrelevant. When comparing 
              % i-modes to i-modes this scale factor is irrelevant. When
              % comparing d-modes to i-modes, the d-modes will win in more
              % head-2-head vector pair spins since d-modes have sVSx more
              % weight relative to i-modes. With the objective for using
              % SPLOC is to discriminant two systems, a 10-fold bias is
              % used. This asymmetry may reduce # of multiple solutions.
              % ----------------------------------------------------------
               inertAddition = inertAddition + 1;    % <= add 1 for i-mode
               end
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& end of bifurcation of subspace selection
               aDevLnScore2D(j) = span;                          % control
               voter2DFrcton(j) = vFraction;                     % control
               got2DEfficacy(j) = vecScore;                      % control
                 if(plotRotationInfo) % => record rotation characteristics
% --------------------------------------------------- 
                 debugEfficacyVecJ(ia,j) = vecScore;
                 debugQualitiesPwJ(ia,j) = qualityFactor;
                 debugConsensusPwJ(ia,j) = vFraction;
                 debugSelectionPwJ(ia,j) = exp(lnScoreFunctn10); 
                   if( (lnScoreFunctn10 > lnMinScore1) && ...
                       (dVoteFraction10 > voteThreshold) && ...
                       (qualityFactor > minQF) ) 
                   debugSubspaceType(ia,j) = 1;  %=> discriminant subspace
                   elseif( (lnScoreFunctn10 < lnMaxScore0) && ... 
                           (iVoteFraction10 > voteThreshold) && ...
                           (qualityFactor > minQF) )
                   debugSubspaceType(ia,j) = 2;  %=> indifference subspace
                   else
                   debugSubspaceType(ia,j) = 3;  %=> undetermined subspace
                   end
                 end
               end                   % end loop for vector pair evaluation
            score2D = got2DEfficacy(1) + got2DEfficacy(2);
               if( score2D > bestScore2D )
               bestScore2D = score2D;
               bestVector2Dscore1 = got2DEfficacy(1);
               bestVector2Dscore2 = got2DEfficacy(2);
               bestAbsDevLnScore1 = aDevLnScore2D(1);
               bestAbsDevLnScore2 = aDevLnScore2D(2);
               bestVoterFraction1 = voter2DFrcton(1);
               bestVoterFraction2 = voter2DFrcton(2);
               bestInertIncrement = inertAddition;
               iaBest = ia;
               end  
            end                          % end of angle subset exploration
         end                                                 % end of zoom
         %pause
% {
% ------------------------------------------ plot rotation characteristics
         if( plotRotationInfo )
         %disp(['best angle = ',num2str( a(iaBest)*180/pi )]);
         %disp(['best score = ',num2str(bestScore2D)]);
         %disp('   ');
         La = (debugSubspaceType(:,1) > 0 );                   % all cases
         xa = a(La)*180/pi;
         % ------------------------------------------------------- index 1
         Ld1 = (debugSubspaceType(:,1) == 1 );
         Li1 = (debugSubspaceType(:,1) == 2 );
         Lu1 = (debugSubspaceType(:,1) == 3 );
         xd1 = a(Ld1)*180/pi;
         xu1 = a(Lu1)*180/pi;
         xi1 = a(Li1)*180/pi;
         % ------------------------------------------------------- index 2
         Ld2 = (debugSubspaceType(:,2) == 1 );
         Li2 = (debugSubspaceType(:,2) == 2 );
         Lu2 = (debugSubspaceType(:,2) == 3 );
         xd2 = a(Ld2)*180/pi;
         xu2 = a(Lu2)*180/pi;
         xi2 = a(Li2)*180/pi;
         % ---------------------------------------------- set up multiplot
         hfig = figure(figure_number0 + 3);
         clf; 
         % --------------------------------------------- 2D efficacy score
         ya_ScoreN = debugEfficacyVecJ(La,1) ...
                   + debugEfficacyVecJ(La,2);
         % ------------------------------------- 2D efficacy decomposition
         yd_Score1 = debugEfficacyVecJ(Ld1,1);
         yd_Score2 = debugEfficacyVecJ(Ld2,2);
         % --------------------------------
         yi_Score1 = debugEfficacyVecJ(Li1,1);
         yi_Score2 = debugEfficacyVecJ(Li2,2);
         % --------------------------------
         yu_Score1 = debugEfficacyVecJ(Lu1,1);
         yu_Score2 = debugEfficacyVecJ(Lu2,2);
         % -------------------------------------
         ymax = 5*ceil( 0.0000000001 + max( max(ya_ScoreN) )/5 );
         subplot(3,3,1);
         hold off;
         plot(xa,ya_ScoreN,'g','linewidth',1.7);
         ylim( [-1,ymax] );
         xlim( [-45,45] );
         ylabel('2D efficacy');
         title11 = ['index pair: (',num2str(j1),',',num2str(j2),')'];
         title(title11);    
         % ------------------------------------------------------ index j1
         subplot(3,3,4);
         hold on;
         plot(xd1,yd_Score1,'^r');
         plot(xu1,yu_Score1,'xk');
         plot(xi1,yi_Score1,'vb');
         ylim( [-1,ymax] );
         xlim( [-45,45] );
         ylabel(['efficacy: ',num2str(j1)]);
         hold off;
         % ------------------------------------------------------ index j2
         subplot(3,3,7);
         hold on;
         plot(xd2,yd_Score2,'^r');
         plot(xu2,yu_Score2,'xk');
         plot(xi2,yi_Score2,'vb');
         ylim( [-1,ymax] );
         xlim( [-45,45] );
         xlabel('angle sweep (degrees)');
         ylabel(['efficacy: ',num2str(j2)]);
         hold off;
         % ------------------------------------------------------- quality
         ydQuality1 = debugQualitiesPwJ(Ld1,1);
         ydQuality2 = debugQualitiesPwJ(Ld2,2);
         % ------------------------------------
         yiQuality1 = debugQualitiesPwJ(Li1,1);
         yiQuality2 = debugQualitiesPwJ(Li2,2);
         % ------------------------------------
         yuQuality1 = debugQualitiesPwJ(Lu1,1);
         yuQuality2 = debugQualitiesPwJ(Lu2,2);
         ymax = max( max( debugQualitiesPwJ(La,:) ) );
         ymax = max( ceil(ymax ),3);
         % ------------------------------------------------------ index j1
         subplot(3,3,2);
         hold on;
         xRectangle = [-44.9,44.9,44.9,-44.9];
         yRectangle = [minQF,minQF,0,0]; 
         fill(xRectangle,yRectangle,[0.88,0.88,0.88],'LineStyle','none');
         % --------------------------------------------------------------
         plot(xd1,ydQuality1,'^r');
         plot(xu1,yuQuality1,'xk');
         plot(xi1,yiQuality1,'vb');
         ylim( [-1,ymax] );
         xlim( [-45,45] );
         ylabel('quality');
         title(['index: ',num2str(j1)]);
         hold off;
         % ------------------------------------------------------ index j2
         subplot(3,3,3);
         hold on;
         xRectangle = [-44.9,44.9,44.9,-44.9];
         yRectangle = [minQF,minQF,0,0]; 
         fill(xRectangle,yRectangle,[0.88,0.88,0.88],'LineStyle','none');
         % --------------------------------------------------------------
         plot(xd2,ydQuality2,'^r');
         plot(xu2,yuQuality2,'xk');
         plot(xi2,yiQuality2,'vb');
         title(['index: ',num2str(j2)]);
         ylim( [-1,ymax] );
         xlim( [-45,45] );
         hold off;
         % ----------------------------------------------------- consensus
         yd_Vote1 = debugConsensusPwJ(Ld1,1);
         yd_Vote2 = debugConsensusPwJ(Ld2,2);
         % ----------------------------------
         yi_Vote1 = debugConsensusPwJ(Li1,1);
         yi_Vote2 = debugConsensusPwJ(Li2,2);
         % ----------------------------------
         yu_Vote1 = debugConsensusPwJ(Lu1,1);
         yu_Vote2 = debugConsensusPwJ(Lu2,2);
         % ------------------------------------------------------ index j1
         subplot(3,3,5);
         hold on;
         xRectangle = [-44.9,44.9,44.9,-44.9];
         yRectangle = [vT,vT,0,0]; 
         fill(xRectangle,yRectangle,[0.88,0.88,0.88],'LineStyle','none');
         % --------------------------------------------------------------
         plot(xd1,yd_Vote1,'^r');
         plot(xu1,yu_Vote1,'xk');
         plot(xi1,yi_Vote1,'vb');
         ylim( [0,1] );
         xlim( [-45,45] );
         ylabel('consensus');
         hold off;
         % ------------------------------------------------------ index j2
         subplot(3,3,6);
         hold on;
         xRectangle = [-44.9,44.9,44.9,-44.9];
         yRectangle = [vT,vT,0,0]; 
         fill(xRectangle,yRectangle,[0.88,0.88,0.88],'LineStyle','none');
         % --------------------------------------------------------------
         plot(xd2,yd_Vote2,'^r');
         plot(xu2,yu_Vote2,'xk');
         plot(xi2,yi_Vote2,'vb');
         ylim( [0,1] );
         xlim( [-45,45] );
         hold off;
         % ----------------------------------------------------- selection
         yd_Spow1 = debugSelectionPwJ(Ld1,1);
         yd_Spow2 = debugSelectionPwJ(Ld2,2);
         % --------------------------------
         yi_Spow1 = debugSelectionPwJ(Li1,1);
         yi_Spow2 = debugSelectionPwJ(Li2,2);
         % --------------------------------
         yu_Spow1 = debugSelectionPwJ(Lu1,1);
         yu_Spow2 = debugSelectionPwJ(Lu2,2);
         ymax = max( max( debugSelectionPwJ(La,:) ) );
         ymax = ceil(ymax);
         ymax = max(ymax,3);
         % ------------------------------------------------------ index j1
         subplot(3,3,8);
         hold on;
         xRectangle = [-44.9,44.9,44.9,-44.9];
         yRectangle = [minScore1,minScore1,maxScore0,maxScore0]; 
         fill(xRectangle,yRectangle,[0.88,0.88,0.88],'LineStyle','none');
         % ---------------------------------------------------------------
         plot(xd1,yd_Spow1,'^r');
         plot(xu1,yu_Spow1,'xk');
         plot(xi1,yi_Spow1,'vb');
         xlabel('angle sweep (degrees)');
         ylabel('selection');
         ylim( [0,ymax] );
         xlim( [-45,45] );
         hold off;
         % ------------------------------------------------------ index j2
         subplot(3,3,9);
         hold on;
         xRectangle = [-44.9,44.9,44.9,-44.9];
         yRectangle = [minScore1,minScore1,maxScore0,maxScore0]; 
         fill(xRectangle,yRectangle,[0.88,0.88,0.88],'LineStyle','none');
         % ---------------------------------------------------------------
         plot(xd2,yd_Spow2,'^r');
         plot(xu2,yu_Spow2,'xk');
         plot(xi2,yi_Spow2,'vb');
         xlabel('angle sweep (degrees)');
         ylim( [0,ymax] );
         xlim( [-45,45] );
         hold off;
% --------------------------- [left bottom width height] position subplots
         set(hfig,'position',  [ 50,  500,  800,  400]);
         sp11 = subplot(3,3,1);
         sp12 = subplot(3,3,2);
         sp13 = subplot(3,3,3);
         %---------------------
         sp21 = subplot(3,3,4);
         sp22 = subplot(3,3,5);
         sp23 = subplot(3,3,6);
         %---------------------
         sp31 = subplot(3,3,7);
         sp32 = subplot(3,3,8);
         sp33 = subplot(3,3,9);
         %--------------------- [left bottom width height] ------- 3rd row
         set(sp31, 'position',  [0.06  0.1   0.25  0.25]);
         set(sp32, 'position',  [0.36  0.1   0.29  0.25]);
         set(sp33, 'position',  [0.68  0.1   0.29  0.25]);
         % ------------------------------------------------------- 2nd row
         set(sp21, 'position',  [0.06  0.4   0.25  0.25]);
         set(sp22, 'position',  [0.36  0.4   0.29  0.25]);
         set(sp23, 'position',  [0.68  0.4   0.29  0.25]);
         % ------------------------------------------------------- 1st row
         set(sp11, 'position',  [0.06  0.7   0.25  0.25]);
         set(sp12, 'position',  [0.36  0.7   0.29  0.25]);
         set(sp13, 'position',  [0.68  0.7   0.29  0.25]);
         % -----------------------------------------------
         %pause(0.5);
% ---------------------------------------------------------- for debugging
%          eScore = sum(efficacyArray);
%          disp(['best 2D efficacy score = ',num2str(bestScore2D)]);
%          disp(['      convergenceRatio = ',num2str(convergenceRatio)]);
%          disp(['                    j1 = ',num2str(j1)]);
%          disp(['     # of vector pairs = ',num2str(mVPc)]);
%          disp(['        max2DscoreDiff = ',num2str(max2DscoreDiff)]);
%          disp(['internal efficacyScore = ',num2str(eScore)]); 
%          disp(['        consensus vote = ',num2str(consensusVote)]); 
%          disp('   ');
%          disp('----------------------------------------------------- ');
% ==================================== information on vector pair sampling
         hfig = figure(figure_number0 + 4);
         clf; 
% -------------------------- plot vector pair rotation importance sampling
         flatProbRotate = sort( rotateProbability(:) ); %sort looks better
         %flatProbRotate = rotateProbability(:);              % DO NOT USE
         set(hfig,'position',[850,   0,  600,  200]); 
                          % [left bottom width height]
         sp21 = subplot(2,1,1);
         bar(flatProbRotate);
         tmp = round(10*gotBVS.efficacy)/10;
         tpt = round(1000*countSucc/countRotations)/10;
         title(['efficacy= ',num2str(tmp), ...
                '    PP-type= ',num2str(pursuitType), ...
                '    resets= ',num2str(nResets), ...
                '    spins= ',num2str(countRotations), ...
                '    %succ= ',num2str(tpt), ...
                '    %skip= ',num2str( ceil(skipRate) )]);
         xlabel('vector pairs (sort ordered)');        % sort looks better
         %xlabel('vector pairs (fixed order)');               % DO NOT USE
         ylabel('spin probability');
         ylim([0,1]);
         sp22 = subplot(2,1,2);
         histogram( penaltyFactor(:),'FaceColor','g', ...  %[0.15,0.1,0.9]
                   'Normalization','probability' );
         xlabel('penalty factor');
         ylabel('probability')
         %--------------------- [left bottom width height] ---------------
         set(sp21, 'position',  [0.08  0.19   0.89  0.72]);
         set(sp22, 'position',  [0.15  0.52   0.3  0.35]);
         %pause
         debugSubspaceType = 0*debugSubspaceType;      % required to reset
         end
% ----------------------------------------------------------- end debuging
%}
%%                               calculate vector pair selection reduction
% The conceptul idea of the algorithm that implements IMPORTANCE SAMPLING.
%
% This section deals with importance sampling to prioritize which pairs of
% vectors should be considered to spin in their 2D plane. Spinning a pair
% of vectors recurrently at too fast of a rate is not efficient because if
% the vectors lie in the same (or approximately the same) plane, the same
% answer will be obtained upon subsequent spinning. If each of the vectors
% that define a pair have time to partner with other vectors, then when
% they come back together for a spin, they have higher chance of being in
% a different plane than before. In that case, there is a chance that the
% overall efficacy can improve. Even if the efficacy does not improve, it 
% should be checked if the plane represents a new orientation than was not
% checked previously. In principle, one could store the normal vector to a
% plane, and then recalulate a new normal vector, and see the dot-product
% of the current normal to the previous normal.  This would be a nice way
% to make the determination whether to consider a particular vector pair 
% over again or not. However, there are technical difficulties to do this.
% Unfortunately, for nDOF, nDOF - 2 vectors will be orthogonal to this
% plane. Unlike in 3D where the normal to a plane is uniquely defined, at
% least to a +/- sign --> in higher dimensions the normal is not defined. 
% Alternatively the use RMSIP comparing before and after vecotr could be 
% done instead. However, this strategy has baggage because each vector in 
% must be stored using nDOF components, and taking dot products between 
% these nDOF dimensional vectors could be computationally costly. Instead,
% the idea is to place all vector pairs in a queue. A current vector pair
% just checked need not go back to the end of the line. Instead, it is
% assigned a priority. Less priority is given to vectors that have low
% quality properties, or to vectors just checked, since opportunity must
% be given to each of the vectors to EXPLORE spinning with other vectors.
%
% The IMPORTANCE SAMPLING assigns probability for spinning a vector-pair.
% The probability decreases if a vector pair is well converged, or has low
% quality. The probability for spinning a vector-pair increases when the
% two vectors forming the pair are of high quality, not checked recently, 
% and, are volatile thereby not converged. By monitoring all vector pairs,
% it is possible to let probability assignments via importance sampling to
% speed up calculations by not wasting resources on duds, while increasing
% accuracy for the same # of spins because each rotation is relevant for
% potentially improving the efficacy score.
% ------------------------------------------------------------------------
%     REMARK: Consider the function used to calculate Efficacy:
%
%             sWeight = sqrt(span) + span*(1 + span);
%
%             At the threshold just after nucleation, span = gap_lnScore
%             sWeightReference = sWeight(gap_lnScore)
%             bestVector2DscoreREF = minQF*sWeightReference
%             pREF = 1 + 1 + bestVector2DscoreREF
%             cases: A   B   C
%             For index #1:
%          A) When bestVoterFraction1 = voteThreshold/1.25  => cutoff of 0
%          B) When bestAbsDevLnScore1 = gap_lnScore         => cutoff of 0
%          C) When bestVector2Dscore1 = bestVector2DscoreREF =>cutoff of 0
%             Likewise the same applies to index #2. 
%
%             When p1 < 0 => likelihood j1 vector is INTERESTING or at 
%                            least some effort should be made to see if it
%                            becomes more interesting with more iterations
%             When p1 > 0 => surely the current state of the j1-vector is
%                            NOT INTERESTING, and as such less effort will
%                            will be expended on it because the interest
%                            is to find interesting things and not waste
%                            valuable resources on uninteresting things.
%       DELAY JUDGEMENTS --> making judgements about a vector early in the
%                            process is too risky because as vectors pair
%                            up and spin, it is possible to extract some 
%                            interesting information out of a vector as it
%                            changes direction. Also, transfer of some
%                            interesting information can come out of one 
%                            vector and into the other vector such that 
%                            the rich get richer and the poor get poorer.
%                            When the SAME pair is re-encountered after 
%                            some delay, the vector directions are LIKELY
%                            to be different. Therefore a spin of the same
%                            two vectors will be confined to a different 
%                            plane in a high dimensional space, possibly 
%                            orthogonal to the plane from before.
%               STRATEGY --> Try to delay spinning a low-quality vector so
%                            it can make a difference when paired up to a
%                            different vector that has time to change its
%                            direction. When both vectors are of little to
%                            no interest, the delay should be longer. Thus
%                            the probability for a vector pair to spin
%                            must be based on the quality of both vectors.
%               RECOVERY --> When a vector pair is re-encountered, but has
%                            become interesting, reset the penaltyFactor
%                            to one so that it can continue to contribute
%                            to greater efficacy, or interestingness! 
% ------------------------------------------------------------------------
      %vpsr = 1;                                 % => with no optimization
      % {
      p1 = pREF - bestVector2Dscore1 - bestAbsDevLnScore1/gap_lnScore ...
         - (1.25*bestVoterFraction1/voteThreshold)^2;
      p2 = pREF - bestVector2Dscore2 - bestAbsDevLnScore2/gap_lnScore ...
         - (1.25*bestVoterFraction2/voteThreshold)^2;
      p3 = p1 + p2;
      maxPx = max( [p1,p2,p3] );
      %disp(maxPx); 
         if( maxPx > 0 )
         % pPower = maxPx;               less agressive original algorithm
         pPower = min(pREF2,sqrt(maxPx) + maxPx);
         pF = penaltyFactor(j1,j2);
         vpsr = pF^pPower;       % vpsr => vector pair selection reduction
         % --------------------------------------------------- not so fast
         % EXCEPTION: slow down cases when consensus vote is very close to
         %            but not quite at the voteThreshold. These cases must
         %            remain active, but it is expected that there won't 
         %            be many. Hence at this point, there is a chance to
         %            switch to the pPower < 0 case.
         dV1 = bestVoterFraction1 - voteThreshold;
         Lexception1 = and( (dV1 < 0) , (dV1 > -0.075) );
         dV2 = bestVoterFraction2 - voteThreshold;
         Lexception2 = and( (dV2 < 0) , (dV2 > -0.075) );
         sumException = sum(Lexception1 + Lexception2);
           if( sumException > 0 )               % => inspect the exception
           aaa = passoverCount(j1,j2) + bestInertIncrement;
           aaa = max(10,aaa);
           vpsrNEW = 0.9^aaa;
              if( vpsrNEW > vpsr )                      % exception caught
              vpsr = vpsrNEW;
              %disp( [0,pF,vpsr,aaa] );
              % pF stays the same
              passoverCount(j1,j2) = passoverCount(j1,j2) ...
                                   + bestInertIncrement;
              passoverCount(j2,j1) = passoverCount(j1,j2);
            % REMARK: The exception reduces activity to allow for other
            % vector pairs to spin before returning to this pair. A small
            % reduction of activity at first retains moderate frequency of
            % rechecks. As this process continues the reduction increases,
            % moving to faster termination due to no progress.
              else                          % => FALSE ALARM: no exception
              %disp( [pPower,pF,vpsr] );
              pFshift = min(0.05,0.05*pPower/pREF2);     % 0.05 works well
              pF = pF - pFshift*pF*pF;
              penaltyFactor(j1,j2) = pF;
              penaltyFactor(j2,j1) = pF;
              passoverCount(j1,j2) = 0; % <= detected changes => not inert
              passoverCount(j2,j1) = 0; % <= detected changes => not inert
              end
           else                       % => continue as usual: no exception
           %disp( [pPower,pF,vpsr] );
           pFshift = min(0.05,0.05*pPower/pREF2); 
           pF = pF - pFshift*pF*pF;
           penaltyFactor(j1,j2) = pF;
           penaltyFactor(j2,j1) = pF;
           passoverCount(j1,j2) = 0;    % <= detected changes => not inert
           passoverCount(j2,j1) = 0;    % <= detected changes => not inert
           end
         else       % as passoverCount increases => ~converged vector pair
         vpsr = 0.9^passoverCount(j1,j2);       
         %disp( [0,vpsr,passoverCount(j1,j2)] );
         penaltyFactor(j1,j2) = 1;
         penaltyFactor(j2,j1) = 1;
         passoverCount(j1,j2) = passoverCount(j1,j2) + bestInertIncrement;
         passoverCount(j2,j1) = passoverCount(j1,j2);
         end 
      %}
%%                                                record results: update U
         if( iaBest ~= iaZero )     % => rotate vectors within 2D subspace
% -------------------------------- make rotation in high dimensional space
         TR = TRcell{iaBest};
         newU = PR*TR;                              % generate new vectors
         U(:,modeList) = newU;
         %U'*U
         %U
% ------------------------------------------------- record rotated vectors
         efficacyArray(j1) = bestVector2Dscore1;
         efficacyArray(j2) = bestVector2Dscore2;
         voterFraction(j1) = bestVoterFraction1;
         voterFraction(j2) = bestVoterFraction2;
         diff2Dscore = bestScore2D - oldScore2D;
         %disp( diff2Dscore );
% --------------------------------------------------- track maximum change
            if( diff2Dscore > max2DscoreDiff )
            max2DscoreDiff = diff2Dscore;
            end
            if( diff2Dscore < 0.01 )             % => ineffective rotation
            spinProbability = min(1, 0.10*vpsr*rotateProbability(j1,j2) );
            rotateProbability(j1,j2) = spinProbability;
            rotateProbability(j2,j1) = spinProbability; 
            else
            ccc = min(1,100*diff2Dscore);
            countSucc = countSucc + ccc;                      % successful
            spinProbability = min(1, 0.5*vpsr*dampFactor(iaBest)* ...
                                    rotateProbability(j1,j2) );
            rotateProbability(j1,j2) = spinProbability;
            rotateProbability(j2,j1) = spinProbability;
            end
         scoreChange = scoreChange + diff2Dscore;
         % --------------------------------------------- for checking only
%          kSucc = max(min(round( 1000*diff2Dscore ),500),0);
%          kSucc = kSucc + 1;
%          histSucc(kSucc) = histSucc(kSucc) + 1;   
         else                     % |a| = 0 => appy maximum damping factor
%          histSucc(1) = histSucc(1) + 1;                          % check
         spinProbability = min(1, 0.25*vpsr*rotateProbability(j1,j2) );
         rotateProbability(j1,j2) = spinProbability;
         rotateProbability(j2,j1) = spinProbability; 
         end
      countSpin = countSpin + 1;
      end                          % end of rotating specified vector pair
     end                                  % finished with vector pair list
% --------------------- determine most important mode for a complete sweep
%         Note that a complete sweep means spin a reference mode against 
%         all other (nDOF - 1) modes. The question is, which mode should
%         be selected as the reference mode for an entire sweep?
% REMARK: The modes that will be selected most frequently are those that
%         are on the ends of the spectrum. Also, the modes with shorter
%         waiting times will likely be spun more frequently than the modes
%         with longer waiting times. The greater # of successful spins in
%         a single complete sweep will shorten the waiting time of a mode.
%         The greater # of failed spins will increase the waiting time of
%         a mode. This creates a bias for higher probability to spin the
%         modes that will most likely increase the efficacy of the basis 
%         vectors both more quickly and more efficiently. A mode should 
%         not be spun too frequently. Therefore the algorithm requires at
%         least one sit-out. If a mode is successfull and gets selected
%         frequently, it will become less successfull since other vectors
%         need time to spin too. As such, the algorithm self-adapts in
%         that no mode will spin excessively too much or too little based
%         on observation and feedback. This process creates a heristic 
%         for an automatated data driven sampling method.
% --------------------------------------------- update delay time schedule
   waitStep_j1 = scoreChange/countSpin;     % mean time step for j1 vector
   waitStep = 0.9*waitStep + 0.1*waitStep_j1; % defines relevant time unit
% REMARK: To understand the waiting time algorithm, set waitStep = 1.
%         Then, realize that as convergence is approached, the change in
%         score --> 0, so this means setting waitStep = 1 is NOT good.
%         Rather, one must consider relative gains of one vector compared
%         to the average gains of all vectors in the system.
   waitTime = waitStep*nVectorPairs;
   modeSelectionTime(j1) = waitTime;
   waitStep_j1 = max(waitStep_j1,waitTime/nDOF); % prevents complete stall
   modeIncrementTime(j1) = waitStep_j1;        % this mode's effectiveness
   sortedWaitTimes = sort(modeSelectionTime);
   shiftWaitTime = sortedWaitTimes(iWaitTime10percent);
   shiftWaitTime = max(0,shiftWaitTime);      % only shorten, not lengthen
   shiftWaitTime = shiftWaitTime + waitTime/nDOF;   % => must sink below 0
                               % ^^^^^^^^^^^^^^^--> patch preventing stall
   modeSelectionTime = modeSelectionTime - shiftWaitTime;
% ------------------------------------------- select next j1 vector to use
   j1Skip = j2Pair(nDOF);                       % this was the previous j1
   tfactor = countSucc/countRotations;
   probRS = 0.8*(1 - tfactor*tfactor); 
      if( rand < probRS )      % randomize selection order 10% of the time
      lineUp = randperm(nDOF);  % attempts to break potential cyclic traps
      else
      lineUp = mapS2U;
      end
      if( rand < prob4j12Bd )    % order in preference of d-mode selection
      gotNext_j1 = -1;                                              % flag 
        for maxDepth=1:nVectorPairs
          for iDepth = 1:maxDepth
          j1_next = lineUp(iDepth);
            if( j1_next ~= j1Skip )
            %disp( num2str([iDepth,j1_next, ...
            %               modeSelectionTime(j1_next), ...
            %               modeIncrementTime(j1_next)]) );
              if( modeSelectionTime(j1_next) < 0 )
              modeSelectionTime(j1_next) = waitTime;
              gotNext_j1 = 1;
              break;
              else
              modeSelectionTime(j1_next) = modeSelectionTime(j1_next) ...
                                         - modeIncrementTime(j1_next);
              end
            end              
          end
          if( gotNext_j1 > 0 )
          break
          end
        end
      else                       % order in preference of i-mode selection
      gotNext_j1 = -1;                                              % flag
        for maxDepth=nVectorPairs:-1:1
          for iDepth = nVectorPairs:-1:maxDepth
          j1_next = lineUp(iDepth);
            if( j1_next ~= j1Skip )
            %disp( num2str([iDepth,j1_next, ...
            %               modeSelectionTime(j1_next), ...
            %               modeIncrementTime(j1_next)]) );
              if( modeSelectionTime(j1_next) < 0 )
              modeSelectionTime(j1_next) = waitTime;
              gotNext_j1 = 1;
              break;
              else
              modeSelectionTime(j1_next) = modeSelectionTime(j1_next) ...
                                         - modeIncrementTime(j1_next);
              end
            end              
          end
          if( gotNext_j1 > 0 )
          break
          end
        end
      end
   Lget_iVP_next = ( j2Pair == j1_next );
   iVP_next = find(Lget_iVP_next); 
%    disp('------------------------------------------------');
%    disp(['iVP_next = ',num2str(iVP_next)]);
%    disp( num2str([j1_next,j2Pair(iVP_next),j2Pair(nDOF)]) );
%    disp('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^');
%    disp('   ');
%    pause
    j2Pair(iVP_next) = j2Pair(nDOF);         % swap out previous j1 vector
   %^^^^^^--------------> not an ordered indexed variable, must use find()
    j2Pair(nDOF) = j1_next;                       % next j1 vector in line
% ----------------------------------------------------- for debugging only
% {
   if( verbosity == 3 )
   showWaitTimePlot = showWaitTimePlot + 1;
      if( mod(showWaitTimePlot,10) == 0 )
      nfig = figure(figure_number0+5);
      clf;
      sortedModeSelectionTime = modeSelectionTime(mapS2U);
      sortedModeIncrementTime = modeIncrementTime(mapS2U);
      subplot(2,1,1);                                  
      ymax = 1.5*waitTime;
      ymin = 1.25*min(sortedModeSelectionTime);
      ymin = max(-0.1,ymin);
         if( doCayleyRotation )
         hold on;
         tempArray = sortedModeSelectionTime;
         Lrandom = ( sortedModeIncrementTime < 1.0e-20 );
         tempArray(Lrandom) = ymax;
         tempArray(~Lrandom) = 0;
         bar(tempArray,'y');
         tempArray = 0*tempArray;
         tempArray(~Lrandom) = sortedModeSelectionTime(~Lrandom);
         tempArray(Lrandom) = 0;
         bar(tempArray,'k');
         else
         bar(sortedModeSelectionTime,'k');
         end
      xlabel('mode index');
      ylabel('delay time');
      cr = ceil(10000*convergenceRatio)/10000;         % convergence ratio
      title(['convR= ', num2str(cr), ...
             '  nCayley= ',num2str( round(xCayley) ), ...
             '  waitTime = ',num2str(waitTime)],'FontSize',12);
      ylim( [ymin,ymax] );
      hold off;
% --------------------------------------------------------- second subplot
      subplot(2,1,2);  
      tempArray = count_j1(mapS2U);
      bar(tempArray,'b');
      %xlabel('mode index')
      ylabel('j1-selection count');
      set(nfig,'position',  [ 50,  0,   350,  340]);
                         % [left bottom width height] 
      sp1 = subplot(2,1,1);
      sp2 = subplot(2,1,2);      
      set(sp1,'position',[0.14 0.11 0.81 0.37]);
      set(sp2,'position',[0.14 0.63 0.81 0.35]);
%       pause(0.5)
%       disp('   ');
%       disp('   ');
%       disp( mapS2U )
%       beep;
%       pause
      end
   end
%}
% ------------------------------------------- select next j1 vector to use
% REMARK: This is the old code before importance sampling via delay time
%         was added by tracking the effectiveness of rotating each vector.
%         In simplest terms, the concept is to swap out one reference
%         vector for another vector. The optimized algorithm builds a
%         queue, with the more successful vectors called more frequently.
%  j1_next = j2Pair(iVP_next);       % gives next vector (no optimization)
%  j2Pair(iVP_next) = j2Pair(nDOF);                  % next vector in line
%  j2Pair(nDOF) = j1_next;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                           get net score
     if( mVPc > nVectorPairsLess1 )  % => evaluation period finishes round
     mVPc = 0;                             % reset counter for fresh start
%%                                       get updated basis vector spectrum
    checkPointWrite = checkPointWrite + 1;
    gvSPLOC.sVS = set_sVS;                         % must reset every time
    gotBVS = getBasisVecSpectrum(U,'checkpoint',trait1,trait0,vT); 
    LoverRideH = ( (voterFraction - voteThreshold) < 0.000001 );
    LoverRideL = ( (voteThreshold - voterFraction) < overRideZone );
    LoverRide = and(LoverRideH,LoverRideL);
    LscoreGood = ( abs(gotBVS.SEVd - gotBVS.SEVi) > 0.3 );
    LqualityGood = ( (gotBVS.QEVd + gotBVS.QEVi) > minQF );
    LoverRide = and(LoverRide,LscoreGood);
    LoverRide = and(LoverRide,LqualityGood);
%%                                                        monitor progress
    nOverRides = sum(LoverRide);
% -------------------------------------------------- various sanity checks
%        if( nOverRides > 0 )
%        %%dScoreArray = abs( gotBVS.SEVd - gotBVS.SEVi );
%        %%tmpArray = [gotBVS.QEVd; gotBVS.QEVi; dScoreArray; LoverRide];
%        LscoreGood = ( abs(gotBVS.SEVd - gotBVS.SEVi) > 0.3 );
%        LqualityGood = ( (gotBVS.QEVd + gotBVS.QEVi) > minQF );
%        tmpArray = [LqualityGood; LscoreGood; LoverRide];
%        disp( tmpArray );
%        pause(0.3)
%        end
    newEfficacy = gotBVS.efficacy;  
    efficacySwapPerc0 = 100*min(1, scoreChange/(newEfficacy + 1.0e-20) );
    efficacySwapPerc = 0.1*efficacySwapPerc0 + 0.9*efficacySwapPerc;
       if( pursuitType == 1 )
       new_efficacyScore = sum(gotBVS.EEVd);
       elseif( pursuitType == -1 )
       new_efficacyScore = sum(gotBVS.EEVi);
       else
       new_efficacyScore = sum(gotBVS.EEVd + gotBVS.EEVi);
       end  
    consensusVec = (voterFraction - voteThreshold)/(1 - voteThreshold);
    Lvote = (consensusVec > 0);
    nVotes = sum(Lvote);
       if( nVotes > 0 )
       consensusVec(~Lvote) = 0;
       tmpYEA = sum(consensusVec); 
       consensusVote = voteThreshold ...
                     + (tmpYEA/nVotes)*(1 - voteThreshold);
       else
       consensusVote = 0;
       end
    tmp_vF = (voterFraction/voteThreshold).^2;
    vScore = sum(tmp_vF);                    % weight scores above vT more
    scoreIncrement = min(abs(vScore - old_vScore),1);          % <= actual
    old_vScore = vScore;
    runningAveScoreIncrement = 0.1*scoreIncrement ...
                             + 0.9*runningAveScoreIncrement;
    convergenceRatio = runningAveScoreIncrement/vScore;
    convergenceRatio = min(1,convergenceRatio);
% --------------------------------- calculate running average of skip rate
    dCountRotations = countRotations - oldCountRotations;
       if( dCountRotations > 100 )
       dCountSkip = countSkip - oldCountSkip;
       incrementalSkipFraction = dCountSkip/dCountRotations;
       skipFraction = 0.1*incrementalSkipFraction + 0.9*skipFraction;
       oldCountSkip = countSkip;
       oldCountRotations = countRotations;
       end
    skipRate = ceil(1000*skipFraction)/10;           % calculate skip rate
%%                                                      record information
      if( checkPointWrite > writeCheckPoint )      % => intermediate write
      %writeSPLOCresults(splocFileName,gotBVS);
      % this was commented out to reduce risk that MATLAB slows down when
      % there are a lot of successive writing with large files. This may
      % not improve speed, but this feature isn't being used anyway.
      checkPointWrite = 0;
% --------------------------------------- write to screen some information
         if( verbosity > 1 )
         dd = gotBVS.Dd;
         du = gotBVS.Du;
         di = gotBVS.Di;
         d3 = dd + du + di;
         disp('   ');
         disp(['              pursuit type = ',num2str(pursuitType)]);
         disp(msgInitial_U);
         disp(['            vote threshold = ',num2str(voteThreshold)]);
         disp([' discriminant subspace DIM = ',num2str(dd)]);
         disp([' undetermined subspace DIM = ',num2str(du)]);
         disp([' indifference subspace DIM = ',num2str(di)]);
         disp(['    vector space dimension = ',num2str(d3)]);
         disp(['distinct # of vector pairs = ',num2str(nvp)]); 
         disp([' total # of pair rotations = ',num2str(countRotations)]);
         aveSpinN = countRotations/nvp;
         disp([' mean spin # per Vec. pair = ',num2str(aveSpinN)]);
         disp(['     # of successful spins = ',num2str(countSucc)]);
         successRate = 100*countSucc/countRotations;
         disp(['      percent success rate = ',num2str(successRate)]);
         disp(['    # of skipped rotations = ',num2str(countSkip)]);
         disp(['runningAverage % skip rate = ',num2str(skipRate)]);
         disp('----------------------------------------- control');
         disp(['               # of resets = ',num2str(nResets)]);
         disp(['      damped learning rate = ',num2str(dampLR)]);
         disp(['   random matrix amplitude = ', ...
               num2str(randomMatrixAmplitude)]);
         disp([' maximum rotation fraction = ', ...
               num2str(maxDrotateFraction)]);
         aveDrotate = sumDrotate/Nrotate;
         disp([' <# of u-vector rotations> = ',num2str(aveDrotate)]);
         disp(['            scoreIncrement = ',num2str(scoreIncrement)]);
         disp(['  runningAveScoreIncrement = ', ...
               num2str(runningAveScoreIncrement)]);
         disp('-------------------------------------------------');
         disp(['    # of failed gains > 1% = ',num2str(nFailed)]);
         disp(['            total efficacy = ',num2str(gotBVS.efficacy)]);
         disp(['            consensus vote = ',num2str(consensusVote)]);
         disp(['effective diffusive spin # = ',num2str(round(xCayley) )]);
         disp(['%efficacy change per sweep = ', ...
              num2str(efficacySwapPerc)]);
         disp(['         convergence ratio = ', ...
               num2str(convergenceRatio)]);
         disp('   ');
         %pause(2.5) 
         end
         if( verbosity > 0 )
         dd = gotBVS.Dd;
         du = gotBVS.Du;
         di = gotBVS.Di;
         d3 = dd + du + di;
         msg = '    ';
         fprintf(fid,'%s \n',msg);
         msg = ['              pursuit type = ',num2str(pursuitType)];
         fprintf(fid,'%s \n',msg);
         fprintf(fid,'%s \n',msgInitial_U);
         msg = ['            vote threshold = ',num2str(voteThreshold)];
         fprintf(fid,'%s \n',msg);
         msg = [' discriminant subspace DIM = ',num2str(dd)];
         fprintf(fid,'%s \n',msg);
         msg = [' undetermined subspace DIM = ',num2str(du)];
         fprintf(fid,'%s \n',msg);
         msg = [' indifference subspace DIM = ',num2str(di)];
         fprintf(fid,'%s \n',msg);
         msg = ['    vector space dimension = ',num2str(d3)];
         fprintf(fid,'%s \n',msg);
         msg = [' distinct # of vector pairs = ',num2str(nvp)];
         fprintf(fid,'%s \n',msg);
         msg = [' total # of pair rotations = ',num2str(countRotations)];
         fprintf(fid,'%s \n',msg);
         aveSpinN = countRotations/nvp;
         msg = [' mean spin # per Vec. pair = ',num2str(aveSpinN)];
         fprintf(fid,'%s \n',msg);    
         msg = ['     # of successful spins = ',num2str(countSucc)];
         fprintf(fid,'%s \n',msg);
         successRate = 100*countSucc/countRotations;
         msg = ['      percent success rate = ',num2str(successRate)];
         fprintf(fid,'%s \n',msg);
         msg = ['    # of skipped rotations = ',num2str(countSkip)];
         fprintf(fid,'%s \n',msg);
         msg = ['runningAverage % skip rate = ',num2str(skipRate)];
         fprintf(fid,'%s \n',msg);
         msg = ['      damped learning rate = ',num2str(dampLR)];
         fprintf(fid,'%s \n',msg);
         msg = ['    # of failed gains > 1% = ',num2str(nFailed)];
         fprintf(fid,'%s \n',msg);
         msg = ['            total efficacy = ',num2str(gotBVS.efficacy)];
         fprintf(fid,'%s \n',msg);
         msg = ['            consensus vote = ',num2str(consensusVote)];
         fprintf(fid,'%s \n',msg);
         msg = ['               # of resets = ',num2str(nResets)];
         fprintf(fid,'%s \n',msg);
         msg = ['effective diffusive spin # = ',num2str(round(xCayley))];
         fprintf(fid,'%s \n',msg);
         msg = ['%efficacy change per sweep = ', ...
                num2str(efficacySwapPerc)];
         fprintf(fid,'%s \n',msg);
         msg = ['         convergence ratio = ', ...
                num2str(convergenceRatio)];
         fprintf(fid,'%s \n',msg);
         dt = cputime - t0;
         msg = ['    cumulative CPU seconds used = ',num2str(dt)];
         fprintf(fid,'%s \n',msg);                     
         end 
% ---------------------------------------------- plot intermediate results
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
% efficacy  -> stacked bar graph based on the selection power indexing
%           1st put total efficacy as red (discriminant on top)
%           Put contribution of i-subspace as blue or u-subpace as yellow
%  if no discriminant subspace contribution
%           1st put total efficacy as blue (indifference on top)
%           Put contribution of u-subspace as yellow 
%  if no d-subspace or i-subspace
%           make total efficacy bar as yellow.
%                              **This example corresponds to above example
%                                     B
%             R                   B   B
%         R   R                   B B B
%         R R R                 B B B B
%       R R R R               B B B B B
%       R R R B             B B B B B B
%       R R R B       Y     B B B B B B
%       R B R B     Y Y   Y Y B B B B B
%       Y B R B     Y Y Y Y Y Y B Y B B
%       Y B R B Y Y Y Y Y Y Y Y B Y Y B
% ------------------------------------------------------------------------
%   quality: put +Quality for d-subspace  and  -Quality for i-subspace
% consensus: put -vote    for d-subpace   and  -vote    for i-subspace
% ------------------------------------------------------------------------
% maxSimilaritity and minSimilarity
% For d-suspace:  make bar graph (red) for maxSimilarity
%                 as a sub-component (gray) fill bar with (max-min)
% For i-suspace:  make bar graph (blue) for maxSimilarity
%                 as a sub-component (gray) fill bar with (max-min)
% ------------------------------------------------------------------------
         if( verbosity == 3 )
         Cind = gotBVS.Cind;                       % got gotBVS from above
         L3 = ( Cind > 0 );                          % captures Dd and Ddi
         L2 = ( Cind == 2 );         % Ddi: dual-purpose (d & i) subspaces
        %Ld = ( Cind == 1 );                         % Dd: d subspace only
         Lu = ( Cind == 0 );                         % Du: u subspace only
         Li = ( Cind == -1);                         % Di: i subspace only
         sI = 1:d3;
         SPbv = exp(lnMidScoreX);
         SPbvData = SPbv*ones(1,d3);
         hfig = figure(figure_number0+1);
         clf;
         % ----------------------------------------------------- selection
         subplot(4,1,4);
         hold on;
         sPd = gotBVS.SEVd;   
         sPi = gotBVS.SEVi;
% ----------------------------------------------------------- build legend
         xx1 = [0.6,d3+0.4,d3+0.4,0.6];                 % bin width is 0.8
         ymax = 0.6*minScore1 + 0.4*maxScore0;
         ymin = 0.4*minScore1 + 0.6*maxScore0;
         yys = [ymax,ymax,ymin,ymin]; 
         h1 = fill(xx1,yys,'r','LineStyle','none');
         h2 = fill(xx1,yys,'y','LineStyle','none');
         h3 = fill(xx1,yys,'b','LineStyle','none');
         yys = [minScore1,minScore1,maxScore0,maxScore0]; 
         h4 = fill(xx1,yys,[0.88,0.88,0.88],'LineStyle','none'); 
            if( sum(L3) > 0 )
            rData = SPbvData;
            rData(L3) = sPd(L3);
            a1 = bar(sI,rData,'r');
            a1(1).BaseValue = SPbv;
            end
            if( sum(Lu) > 0 )
            tData = SPbvData;
            bData = SPbvData;
            tData(Lu) = sPd(Lu);
            bData(Lu) = sPi(Lu);
            a2 = bar(sI,tData,'y');
            b2 = bar(sI,bData,'y'); 
            a2(1).BaseValue = SPbv;
            b2(1).BaseValue = SPbv;    
            end
            if( sum(Li) > 0 )
            bData = SPbvData;
            bData(Li) = sPi(Li);
            b3 = bar(sI,bData,'b');
            b3(1).BaseValue = SPbv;
            end
         ymax = max( ceil(sPd) );
         ymax = max(3,ymax);
         ylim( [1,ymax] );
         xlabel('selection index');
         ylabel('selection');
         legend([h1,h2,h3,h4],{'discriminant','undetermined', ...
                           'indifference','not significant'}, ...
                           'location','northeast');
         hold off;
         % ----------------------------------------------------- consensus
         zData = zeros( size(sI) );
         subplot(4,1,3);
         hold on;
         yyv = [vT,vT,-vT,-vT]; 
         fill(xx1,yyv,[0.88,0.88,0.88],'LineStyle','none');
         cPd = gotBVS.CEVd;
         cPi = -gotBVS.CEVi;
            if( sum(L3) > 0 )
            rData = zData;
            rData(L3) = cPd(L3);
            bar(sI,rData,'r');
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
            bData = zData;
            bData(Li) = cPi(Li);
            bar(sI,bData,'b');
            end
         ylabel('consensus');
         ylim( [-1,1] );
         yticks([-1 0 1]);
         yticklabels({'1','0','1'});
         hold off;
         % -------------------------------------------------- mode quality
         subplot(4,1,2);
         hold on;
         yyq = [minQF,minQF,-minQF,-minQF]; 
         fill(xx1,yyq,[0.88,0.88,0.88],'LineStyle','none');
         qPd = max(0,gotBVS.QEVd);
         qPi = -max(0,gotBVS.QEVi);
            if( sum(L3) > 0 )
            rData = zData;
            rData(L3) = qPd(L3);
            bar(sI,rData,'r');
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
            bData = zData;
            bData(Li) = qPi(Li);
            bar(sI,bData,'b');
            end
         ylim( [-10,10] );
         yticks( [-10,-5,0,5,10] );
         yticklabels({'10','5','0','5','10'})
         ylabel('quality');
         % ------------------------------------------------- mode efficacy
         subplot(4,1,1);
         hold on;
         ePd = max(0,gotBVS.EEVd);
         ePi = -max(0,gotBVS.EEVi);
            if( sum(L3) > 0 )
            rData = zData;
            rData(L3) = ePd(L3);
            bar(sI,rData,'r');
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
            bData = zData;
            bData(Li) = ePi(Li);
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
         ylim( [ymin,ymax] );
         ylabel('efficacy');
         hold off
         % --------------------------------------------- position subplots
         set(hfig,'position',[850,  370,  600, 530]);
                          % [left bottom width height]
         sp1 = subplot(4,1,1);
         sp2 = subplot(4,1,2);
         sp3 = subplot(4,1,3);
         sp4 = subplot(4,1,4);
         set(sp4,'position',[0.07 0.07 0.89 0.20]);
         set(sp3,'position',[0.07 0.31 0.89 0.20]);
         set(sp2,'position',[0.07 0.55 0.89 0.20]);
         set(sp1,'position',[0.07 0.79 0.89 0.20]);
% =========================================================== similarities
       % {
          if( gotBVS.Du < nDOF )
          hfig = figure(figure_number0+2);
          clf;
          % --------------------------------------------------------------
          B3 = or(Li,L2);                         % => captures Di and Ddi
          xdE3 = gotBVS.EEVd(L3);      % => d subspace or (d & i)-subspace
          xiE3 = gotBVS.EEVi(B3);      % => i subspace or (d & i)-subspace
          xmin = max(0, floor( min( [xdE3,xiE3] ) ) );
          xmax =  ceil( max( [xdE3,xiE3] ) );
             if( xmax - xmin < 1 )
             xmax = xmin + 1;
             end
             if( xmin < 5 )
             xmin = 0;
             end
          % --------------------------------------------------------------
          subplot(2,2,1);
          hold on;
          ydC3 = gotBVS.CEVd(L3);      % => d subspace or (d & i)-subspace
          yiC3 = gotBVS.CEVi(B3);      % => d subspace or (d & i)-subspace
          scatter(xiE3,yiC3,30,'b','filled'); 
          scatter(xdE3,ydC3,50,'r');         
          ylabel('vFraction')
          xlim( [xmin,xmax] );
          hold off;
          % --------------------------------------------------------------
          subplot(2,2,2);
          hold on;
          ydU3 = gotBVS.USVd(L3);      % => d subspace or (d & i)-subspace
          yiU3 = gotBVS.USVi(B3);      % => d subspace or (d & i)-subspace
          scatter(xiE3,yiU3,30,'b','filled'); 
          scatter(xdE3,ydU3,50,'r');  
          ylabel('maxSimilar')
          xlim( [xmin,xmax] );
          ylim( [0,1.1] );
          hold off;
          % --------------------------------------------------------------
          subplot(2,2,3);
          hold on;
          ydS3 = gotBVS.SEVd(L3);      % => d subspace or (d & i)-subspace
          yiS3 = gotBVS.SEVi(B3);      % => d subspace or (d & i)-subspace
          scatter(xiE3,yiS3,30,'b','filled'); 
          scatter(xdE3,ydS3,50,'r');
          xlabel('efficacy');
          ylabel('selection')
          xlim( [xmin,xmax] );
          hold off;
          % --------------------------------------------------------------
          subplot(2,2,4);
          hold on;
          ydL3 = gotBVS.LSVd(L3);      % => d subspace or (d & i)-subspace
          yiL3 = gotBVS.LSVi(B3);      % => d subspace or (d & i)-subspace
          scatter(xiE3,yiL3,30,'b','filled'); 
          scatter(xdE3,ydL3,50,'r'); 
          xlabel('efficacy');
          ylabel('minSimilar');
          xlim( [xmin,xmax] );
          ylim( [0,1.1] );
          hold off;      
          % ---------------- [left bottom width height] ------------------
          set(hfig,'position',[ 400,   0,  450,  340]);
          sp1 = subplot(2,2,1);
          sp2 = subplot(2,2,2);
          sp3 = subplot(2,2,3);
          sp4 = subplot(2,2,4);
          set(sp4,'position',[0.58 0.10 0.40 0.40]);
          set(sp3,'position',[0.10 0.10 0.40 0.40]);
          set(sp2,'position',[0.58 0.56 0.40 0.40]);
          set(sp1,'position',[0.10 0.56 0.40 0.40]);
          pause(0.5)
          end
         %}   
         end
      end                                         % end intermediate write
%%                                                   check for improvement
% ------------------------------------------------------------ condition 1
        if( new_efficacyScore > 1.01*efficacyScore )       
        fail1 = 0;                                          % => succeeded
        efficacyScore = new_efficacyScore;
        else
        fail1 = 1;                                            % => failure
        end
% ------------------------------------------------------------ condition 2
        tmp = abs(newEfficacy - oldEfficacy)/(0.00001 + newEfficacy);
           if( tmp > 0.01 ) 
           fail2 = 0;                                       % => succeeded
           oldEfficacy = newEfficacy;
           else
           fail2 = 1;                                         % => failure
           end
     failLevel = fail1 + fail2;
        switch failLevel
          case 0                            % => both conditions satisfied
          nFailed = 0;                      % => reset successive failures
          case 1                  % => 1 out of 2 conditions are satisfied
          accelerate = ppp*sqrt( max(1,skipFraction) );
          nFailed = nFailed + 0.1*accelerate;
          otherwise  
          nFailed = nFailed + 1;
        end
     end
%%                                               check for stop conditions
       if( countRotations > minRotations )   % necessary condition to stop
          if( (nFailed - 5*nOverRides) > 10 )        
          tryAgain = false;
          end  
          if( convergenceRatio < minConvergenceRatio )
          tryAgain = false;
          end 
          if( ppp > maxLessonBatch )
          tryAgain = false;
          end
       else
       nFailed = 0;
       end
   %disp( [tryAgain, convergenceRatio, minConvergenceRatio, ...
   %       nFailed, dampLR] );
   %disp( num2str([nLessons,dampLR,convergenceRatio, ...
   %               minRotations0,nFailed]) );
   %pause(0.1)
   end                                                     % end of lesson
% ----------------------------------- be more agressive on penalty factors
Lnot1 = ( penaltyFactor < 0.999 );
totalElements = 1.0e-20 + sum( Lnot1(:) );
totalPenalty = 0 + sum( penaltyFactor( Lnot1(:) ) ); 
pFmax = totalPenalty/totalElements;                           % mean PFmax
%disp(pFmax)
%pause
% % figure(101);                                                   % check
% % clf;
% % tttt = histSucc/1000;
% % bar(0:0.001:0.5,tttt);
% % pause(0.75)
end                                                 % BIG LOOP over resets
%%                                      warning for imcomplete convergence
   if( flag_converged < 0 )                             % => not converged
   disp('   ');
   disp(dividerLine);
   disp('WARNING: Not Converged');  
   disp(['        convergenceRatio = ',num2str(convergenceRatio)]);
   end
%%                                                     iterations finished
% sVS = std_sVS is what getBasisVecSpectrum() uses unless forced otherwise
%                                  |----------------------> best basis set
splocResults = getBasisVecSpectrum(U,baseFname,trait1,trait0,vT);
splocResults.pursuitType = pursuitType;                 % retain objective
%writeSPLOCresults(splocFileName,splocResults);
% this was commented out to reduce risk that MATLAB slows down when
      % there are a lot of successive writing with large files. This may
      % not improve speed, but this feature isn't being used anyway.
      % A user can use this command outside of sploc.m when results from
      % sploc are to be recorded. 
dd = splocResults.Dd;
du = splocResults.Du;
di = splocResults.Di;
%%                                             write information on screen
   if( verbosity > 0 )
   disp('   ');
   disp(['              pursuit type = ',num2str(pursuitType)]);
   disp(msgInitial_U);
   disp(['            vote threshold = ',num2str(voteThreshold)]);
   disp([' discriminant subspace DIM = ',num2str(dd)]);
   disp([' undetermined subspace DIM = ',num2str(du)]);
   disp([' indifference subspace DIM = ',num2str(di)]);
   disp(['    vector space dimension = ',num2str(dd+du+di)]);
   disp(['distinct # of vector pairs = ',num2str(nvp)]); 
   disp([' total # of pair rotations = ',num2str(countRotations)]);
   aveSpinN = countRotations/nvp;
   disp([' mean spin # per Vec. pair = ',num2str(aveSpinN)]);
   disp(['     # of successful spins = ',num2str(countSucc)]);
   successRate = 100*countSucc/countRotations;
   disp(['      percent success rate = ',num2str(successRate)]);   
   disp(['    # of skipped rotations = ',num2str(countSkip)]);
   disp(['runningAverage % skip rate = ',num2str(skipRate)]);
   disp(['      damped learning rate = ',num2str(dampLR)]);
   disp(['    # of failed gains > 1% = ',num2str(nFailed)]);
   disp(['            total efficacy = ',num2str(splocResults.efficacy)]);
   disp(['            consensus vote = ',num2str(consensusVote)]);
   disp(['               # of resets = ',num2str(nResets_final)]);
   disp(['effective diffusive spin # = ',num2str( round(xCayley) )]);
   disp(['%efficacy change per sweep = ',num2str(efficacySwapPerc)]);
   disp(['         convergence ratio = ',num2str(convergenceRatio)]);
   dt = cputime - t0;
   disp(['   Cayley CPU seconds used = ',num2str(CayleyCPU)]);
   disp(['    total CPU seconds used = ',num2str(dt)]);
   disp('   ');
   end
%%                                         write information into log file
msg = '   ';
fprintf(fid,'%s \n',msg);
msg = ['----------------------------------------', ...
       '-------------------------- final status:']; 
fprintf(fid,'%s \n',msg);
msg = ['              pursuit type = ',num2str(pursuitType)];
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',msgInitial_U);
msg = ['            vote threshold = ',num2str(voteThreshold)];
fprintf(fid,'%s \n',msg);
msg = [' discriminant subspace DIM = ',num2str(dd)];
fprintf(fid,'%s \n',msg);
msg = [' undetermined subspace DIM = ',num2str(du)];
fprintf(fid,'%s \n',msg);
msg = [' indifference subspace DIM = ',num2str(di)];
fprintf(fid,'%s \n',msg);
msg = ['    vector space dimension = ',num2str(dd+du+di)];
fprintf(fid,'%s \n',msg);
msg = [' distinct # of vector pairs = ',num2str(nvp)];
fprintf(fid,'%s \n',msg);
msg = [' total # of pair rotations = ',num2str(countRotations)];
fprintf(fid,'%s \n',msg);
aveSpinN = countRotations/nvp;
msg = [' mean spin # per Vec. pair = ',num2str(aveSpinN)];
fprintf(fid,'%s \n',msg);    
msg = ['     # of successful spins = ',num2str(countSucc)];
fprintf(fid,'%s \n',msg);
successRate = 100*countSucc/countRotations;
msg = ['      percent success rate = ',num2str(successRate)];
fprintf(fid,'%s \n',msg);
msg = ['    # of skipped rotations = ',num2str(countSkip)];
fprintf(fid,'%s \n',msg);
msg = ['runningAverage % skip rate = ',num2str(skipRate)];
fprintf(fid,'%s \n',msg);
msg = ['      damped learning rate = ',num2str(dampLR)];
fprintf(fid,'%s \n',msg);
msg = ['    # of failed gains > 1% = ',num2str(nFailed)];
fprintf(fid,'%s \n',msg);
msg = ['            total efficacy = ',num2str(splocResults.efficacy)];
fprintf(fid,'%s \n',msg);
msg = ['            consensus vote = ',num2str(consensusVote)];
fprintf(fid,'%s \n',msg);
msg = ['               # of resets = ',num2str(nResets_final)];
fprintf(fid,'%s \n',msg);
msg = ['effective diffusive spin # = ',num2str( round(xCayley) )];
fprintf(fid,'%s \n',msg);
msg = ['%efficacy change per sweep = ',num2str(efficacySwapPerc)];
fprintf(fid,'%s \n',msg);
msg = ['         convergence ratio = ',num2str(convergenceRatio)];
fprintf(fid,'%s \n',msg);
dt = cputime - t0;
msg = ['    total CPU seconds used = ',num2str(dt)];
fprintf(fid,'%s \n',msg);
fclose(fid);
%%                                    summarize sploc results in sploc log
  if( verbosity > 0 )    % do not write sploc results within splocEnsemble
  splocLogFile = gvSPLOC.splocLogFile;
  fid = fopen(splocLogFile,'a');
  fprintf(fid,'%s \n','  ');
  fprintf(fid,'%s \n',[mfilename,'()']);
  msg = dividerLine('calculation summary');
  fprintf(fid,'%s \n',msg);
  fprintf(fid,'%s \n',['               sploc file name = ',baseFname]);
  fprintf(fid,'%s \n',['                  pursuit type = ', ...
                                            num2str(pursuitType)]);
  fprintf(fid,'%s \n',msgInitial_U);
  fprintf(fid,'%s \n',['                vote threshold = ', ...
                                            num2str(voteThreshold)]);
  fprintf(fid,'%s \n',['        discriminant modes: dd = ',num2str(dd)]);
  fprintf(fid,'%s \n',['        undetermined modes: du = ',num2str(du)]);
  fprintf(fid,'%s \n',['        indifference modes: di = ',num2str(di)]);
  fprintf(fid,'%s \n',[' # of variables = dd+du+di = d = ',num2str(nV1)]);
% ---------------------- information on orgin of nonfunctional system data
  fprintf(fid,'%s \n',['  nonfunctional system dataset = ', ...
                       trait0.dataRefName]);
  fprintf(fid,'%s \n',['# of nonfunctional systems: n0 = ',num2str(n0)]); 
  fprintf(fid,'%s \n',['      sample size per system-0 = ',num2str(Ns0)]);
  fprintf(fid,'%s \n',['total sample size for system-0 = ',num2str(nD0)]);
% ------------------------- information on orgin of functional system data
  fprintf(fid,'%s \n',['     functional system dataset = ', ...
                       trait1.dataRefName]);
  fprintf(fid,'%s \n',['   # of functional systems: n1 =  ',num2str(n1)]);
  fprintf(fid,'%s \n',['      sample size per system-1 = ',num2str(Ns1)]);
  fprintf(fid,'%s \n',['total sample size for system-1 = ',num2str(nD1)]);
  fprintf(fid,'%s \n',dividerLine);
  fclose(fid);
  end
end
