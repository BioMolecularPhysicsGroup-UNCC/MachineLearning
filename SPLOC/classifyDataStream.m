function [ranking,T] = classifyDataStream(fname,dataU,dataF,dataN, ...
                                          U,verbosity)
% calculates PROB[unclassified system is functional] based on overlaps of
% probability density functions together with activation functions to  
% calculate the probability that a system viewed along a particular basis
% vector direction has functional properties. This calculated probability
% is context dependent, with an elaborate procedure that is implemented in
% this function. The probabilities calculated per basis vector depend on 
% the separation between known functional and nonfunctional systems and on
% the unknown system. 
%
% A brief summary of the math is given that is appied to a set of PDFs for
% systems that function {F} and systems that do not function {N}. Each PDF
% is for a particlar unit-vector direction. The basis vector describing a 
% direction is called a mode. Overlaps of PDFs are then translated into
% probabilities that {F} and {N} systems are similar or dissimilar to one
% another as well as a system with unknown function {U} is similar or 
% dissimilar to either an {F}-system or a {N}-system. Note that a system
% that is unclassified is called unknown (U): Based on all this ... 
% MATH SUMMARY:
%
%  pUfunct = PROB[unclassified system is functional]        (Final result)
% ----------------------
%      pUk = prob. U is        functional w.r.t. mode k based on all {F}.
%      qUk = prob. U is NOT nonfunctional w.r.t. mode k based on all {N}.
%            Treat pUk and qUk as independent probabilities since the 
%            systems {F} and {N} are unrelated. 
% ----------------------
%   pUFk_F = conditional probabilty that U mimicks F w.r.t. mode k
%            given that F is functional. 
%      wFk = probability that the system labeled as F is distinct from all
%            nonfunctional systems in mode k.
%     pUFk = pUFk_F*wFk    => probability U functions like F w.r.t. mode k
% 1 - pUFk = probability U is nonfunctional in mode k.
%  1 - pUk = product over all F of (1 - pUFk)
%            This product represents an AND operation for all F to be
%            distinct from U for non-function to take place. Note that
%            if U mimics only 1 F very well, nonfunction is null. 
%*-->  pUk = 1 - product over all F of (1 - pUFk) 
% ----------------------
%   qUNk_N = conditional probabilty that U mimicks N w.r.t. mode k
%            given that N is nonfunctional. 
%      wNk = probability that the system labeled as N is distinct from all
%            functional systems in mode k. 
%     qUNk = qUNk_N*wNk => probabilty U does not function like N in mode k
% 1 - qUNk = probability U is functional in mode k. 
%*-->  qUk = product over all N of (1 - qUNk)  
%            This product represents an AND operation for all N to be
%            distinct from U for function to take place. Note that
%            if U mimics only 1 N very well, function is null. 
% ----------------------
%      wUk = pUk*qUk = probability U is functional in mode k.
%            This product represents an AND operation saying that U should
%            have characteristics that are similar or the same as those 
%            characteristics shared by systems that function while at the 
%            same time be as dissimilar from those characteristics found
%            in systems having nonfunction.  
%  method A:                                  % works, but is conservative
%  pUFprod = product over all k of wUk 
%  pUfunct = (pUFprod)^(1/Nmodes)  => geometric mean of wUk, since this is
%            an AND operation. Reporting this geometrical mean removes the
%            relative number of discriminant modes as a factor.
%            For numerical stability: let wUk --> max(wUk,1.0e-125)
%  method B: 
%  pUfunct = root mean square average of wUk      works, but has high risk
%  method C: pUfunct = hybridMethodFunction(pUfunctA,pUfunctB)
%            pwA = log10(pUfunctA);           % typically too conservative
%            pwB = log10(pUfunctB);                  % typically too risky
%            pwD = min(10*pUfunctB + (pwA - pwB), 0); %powers are negative
%            pwD = pwD/(1 + 10*pUfunctB);   % variable scale for reduction
%            pUfunct = pUfunctB*10^pwD;         % pwD => reduction in risk
% ------------------------------------------------------------------------
%
% DEFINITIONS
% dataMatrixInfo. <-- data structure
% dataRefName    = reference name for A matrix data with similar traits
% dataMatrixName = cell array for file names that store the A matrix data 
%     nVariables = number of variables defining size of vector space
%          pType = packing type (format) applied to the vector space
%            dim = # of components in local vector (e.g. x,y,z for 1 atom)
%              n = number of data matrices
%              A = data matrix in the rowVar format
%       nSamples = number of samples for each data matrix 
%
% -------------------------------------------------------- local variables
% INPUT
% fname = file name to specify origin of data and/or identify the output.
% dataU = dataMatrixInfo for unclassified systems.
% dataF = dataMatrixInfo for functional systems.
% dataN = dataMatrixInfo for nonfunctional systems.
% U = basis vectors that are to be used for the classification process.
%     Note that U can be catenated from an ensemble of viable sets of U.
%     In the case that U represents multiple U sets, the vectors contained
%     within U (ie the modes) will not satisfy all-pairwise orthogonality.
%     This means some redundant information will be present in the product
%     probability. In principle the probabilities employed in the product
%     probability calculation should pertain to independent events, which 
%     is the case for orthogonal modes. With redundancy present, the final
%     likelihoods will decrease as more redundant information is added.
%     This reduction is uniformly applied to all predictions. Therefore,
%     rank ordering may not be affected. Wether classification will become
%     better or worse using redundant information is unknown. 
%
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
%
% PROCESS
%++  Per basis vector (or per mode):
% 1) calculate Fpdf using known pooled data matrix for F-systems.
% 2) calculate Npdf using known pooled data matrix for N-systems.
% 3) calculate Updf using test data matrix for U-systems. 
% 4) calculate overlap integrals for all NF, UF and UN pairs per mode.
%   --> NF and UF pairs use discriminant-sigmodal since interested in how 
%       different F is to either N or U. F should be different to N.
%   --> UN pairs use indifference-sigmodal since interested in how similar
%       U is to N. Only N-states with significant similarity are relevant.
% 5) Store all overlap integrals in arrays. 
% 6) perturb stored overlap values with added small random Gaussian noise. 
% 7) calculate {wFk, pUFk_F, pUFk, pUk} and {wNk, qUNk_N, qUNk, qUk} 
% 8) combine pUk and qUk via product to arrive at wUk
% 9) take geometric mean of wUk to arrive a final probability estimate
% 10) for error analysis perform 626 calculations with randomized values
% 11) determine mean value and standard deviation
% 12) report mean and std. dev. for classification prediction in summary
%
%
% OUTPUT systemRanking. <-- data structure
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
%%                                                    set sploc parameters
global gvSPLOC
% ---------------------------- define discriminant weighting from overlaps
dWtOVL = gvSPLOC.dWtOVL;          % discriminant probability weight factor
%iWtOVL = gvSPLOC.iWtOVL;         % indifference probability weight factor
dxOVL = gvSPLOC.dxOVL;
% --------------------------------------------------- for recording action
rankType = 'dataStream';
splocLogFile = gvSPLOC.splocLogFile;              % to record in sploc log
%%                                                         check verbosity
   if( nargin == 5 )
   verbosity = 1;                                    % => normal operation
   elseif( nargin < 5 )
   error('minimum specification is: fname, dataU, dataF, dataN, U');
   else
   verbosity = setVerbosity(verbosity);
   end
%%                                                 unpackage data matrices
Nu = dataU.n;                                       % # of unknown systems
Au = dataU.A;
aMname = dataU.dataMatrixName;           % cell array of data matrix names
% -------------------------------------
N1 = dataF.n;                                    % # of functional systems
A1 = dataF.A;
% ---------------
N0 = dataN.n;                                 % # of nonfunctional systems
A0 = dataN.A;
%%                                                      simple error check
[nDOF,Nsu] = size(Au{1}); % nDOF= # of degrees of freedom = # of variables
[nV,nModes] = size(U);         % nModes = # of basis vectors to be checked
   if( nDOF ~= nV )
   disp(['Au:   nDOF = ',num2str(nDOF)]);
   disp([' U:   nDOF = ',num2str(nV)]);
   disp([' U: nmodes = ',num2str(nModes)]);
   error('number of variables differ between Au and basis set dimension');
   end
% ------------------------------------------------------- get sample sizes 
[~,Ns1] = size(A1{1});
[~,Ns0] = size(A0{1});
%%                                      set verbosity output level details
 if( verbosity > 0 )
 subFolder = 'classification';
 fName = [fname,'_pdf.log']; 
 %               ^^^^^^^^^^^--------> user controls uniqueness in basename
 iLogFileName = getOutputFileName(subFolder,fName);
%  disp(iLogFileName);                                     % for debugging
   if( isfile(iLogFileName) )
   fid = fopen(iLogFileName,'a');     % => do not need header: append data
   else                                % => virgin file: needs header info
   fid = fopen(iLogFileName,'w');          % => will overwrite old results
   msg ='layout for how probabilities are calculated for classification';
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n','   ');
   msg ='pUfunct = PROB[unclassified system will function]  Final result';
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n',dividerLine);
   msg ='    pUk = prob. U is functional w.r.t. mode k based on all {F}';
   fprintf(fid,'%s \n',msg);
   msg = ['    qUk = prob. U is NOT nonfunctional w.r.t. mode ', ... 
          'k based on all {N}.'];
   fprintf(fid,'%s \n',msg);
   msg ='          Treat pUk and qUk as independent probabilities since ';
   fprintf(fid,'%s \n',msg);
   msg ='          the systems {F} and {N} are unrelated.';
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n',dividerLine);
   msg =' pUFk_F = conditional probabilty that U mimicks F w.r.t. mode k';
   fprintf(fid,'%s \n',msg);
   msg ='          given that F is functional.';
   fprintf(fid,'%s \n',msg);
   msg ='    wFk = probability that the system labeled as F is distinct ';
   fprintf(fid,'%s \n',msg);
   msg ='          from all nonfunctional systems in mode k.';
   fprintf(fid,'%s \n',msg);
   msg ='   pUFk = pUFk_F*wFk  => Prob. U functions like F w.r.t. mode k';
   fprintf(fid,'%s \n',msg);
   msg =' 1-pUFk = probability U is nonfunctional in mode k.';
   fprintf(fid,'%s \n',msg);
   msg ='1 - pUk = product over all F of (1 - pUFk)';
   fprintf(fid,'%s \n',msg);
   msg ='          This product represents an AND operation for all F to';
   fprintf(fid,'%s \n',msg);
   msg ='          be distinct from U for non-function to take place.';
   fprintf(fid,'%s \n',msg);
   msg ='          If U mimics only 1 F very well, nonfunction is null.';
   fprintf(fid,'%s \n',msg);
   msg ='*-> pUk = 1 - product over all F of (1 - pUFk)';
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n',dividerLine);
   msg =' qUNk_N = conditional probabilty that U mimicks N w.r.t. mode k';
   fprintf(fid,'%s \n',msg);
   msg ='          given that N is nonfunctional.';
   fprintf(fid,'%s \n',msg);
   msg ='    wNk = probability that the system labeled as N is distinct';
   fprintf(fid,'%s \n',msg);
   msg ='          from all functional systems in mode k.';
   fprintf(fid,'%s \n',msg);
   msg ='   qUNk = qUNk_N*wNk => Prob. U does not function like N in mode k';
   fprintf(fid,'%s \n',msg);
   msg =' 1-qUNk = probability U is functional in mode k. ';
   fprintf(fid,'%s \n',msg);
   msg ='*-> qUk = product over all N of (1 - qUNk)';
   fprintf(fid,'%s \n',msg);
   msg ='          This product represents an AND operation for all N';
   fprintf(fid,'%s \n',msg);
   msg ='          to be distinct from U for function to take place.';
   fprintf(fid,'%s \n',msg);
   msg ='          If U mimics only 1 N very well, function is null.';
   fprintf(fid,'%s \n',msg);
   msg ='    wUk = pUk*qUk = probability U is functional in mode k.';
   fprintf(fid,'%s \n',msg);
   msg ='          This product represents an AND operation saying that U';
   fprintf(fid,'%s \n',msg);
   msg ='          has characteristics that are similar or the same as';
   fprintf(fid,'%s \n',msg);
   msg ='          the characteristics shared by systems that function ';
   fprintf(fid,'%s \n',msg);
   msg ='          while at the same time be as dissimilar from those ';
   fprintf(fid,'%s \n',msg);
   msg ='          characteristics found in systems having nonfunction.';
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n',dividerLine);
   msg ='pUFprod = product over all k of wUk ';
   fprintf(fid,'%s \n',msg);
   msg ='pUfunct = (pUFprod)^(1/Nmodes)  => geometric mean of wUk, since';
   fprintf(fid,'%s \n',msg);
   msg ='          this is an AND operation. Reporting this geometrical';
   fprintf(fid,'%s \n',msg);
   msg ='          mean removes the relative number of discriminant ';
   fprintf(fid,'%s \n',msg);
   msg ='          modes as a factor that depends on system size.';
   fprintf(fid,'%s \n',msg);
   msg ='          For numerical stability: let wUk --> max(wUk,1e-125)';
   fprintf(fid,'%s \n',msg);
   fprintf(fid,'%s \n',dividerLine);
   msg ='          modify likelihood for improvisational exploration';
   fprintf(fid,'%s \n',msg);
   msg ='          Let A = 1 - sqrt(pUfunct) ';
   fprintf(fid,'%s \n',msg);
   msg ='          Let B = 0.0001*sqrt( <pUk^2> ) ';
   fprintf(fid,'%s \n',msg);
   msg ='          pUfunct --> pUfunct +  A*B';
   fprintf(fid,'%s \n',msg);  
   fprintf(fid,'%s \n',dividerLine);
   end
 end
   if( verbosity > 1 )
   figure;
   figure_number0 = get(gcf,'Number');
   end
%%                                             pre-calculate PDFs per mode
% --------------------------------------------- calculate Fpdf, Npdf, Updf
PDFu = cell(Nu,nModes);
PRJu = cell(Nu,nModes);
PDF1 = cell(N1,nModes);
PRJ1 = cell(N1,nModes);
PDF0 = cell(N0,nModes);
PRJ0 = cell(N0,nModes);
% ------------------------------------------- project into unknown systems
   for ju=1:Nu 
   [modePROJ,modePDF] = generateModePDF(Au{ju},U);
      for j=1:nModes
      PRJu{ju,j} = modePROJ{j};
      PDFu{ju,j} = modePDF{j};
      end
   end
% ---------------------------------------- project into functional systems
   for j1=1:N1 
   [modePROJ,modePDF] = generateModePDF(A1{j1},U);
      for j=1:nModes
      PRJ1{j1,j} = modePROJ{j};
      PDF1{j1,j} = modePDF{j};
      end
   end
% ------------------------------------- project into nonfunctional systems
   for j0=1:N0 
   [modePROJ,modePDF] = generateModePDF(A0{j0},U);
      for j=1:nModes
      PRJ0{j0,j} = modePROJ{j};
      PDF0{j0,j} = modePDF{j};
      end
   end
%%           classification probabilities based on reference training data
% REMARK 1: Classification predictions can only be as good as the training
% data is. The training data can be used to rank the effectiveness of the
% basis vectors that are being used for discrimination. The basis vectors
% can come from anywhere, so the first step is to define a probability
% for a given mode to be effective in contrasting a F-system from being
% a N-system. This is done mathematically as follows: 
% For complete details of all mathematics, see above, MATH SUMMARY:
% ---------------------------------- calculate and store overlap integrals
overlapFNk = zeros(N1,N0,nModes);
overlapUNk = zeros(Nu,N0,nModes);
overlapUFk = zeros(Nu,N1,nModes);
   for k=1:nModes                    % consider each j-th vector-direction
% --------------------------------------------------------------- FN pairs
      for j1=1:N1                % consider each functional 1-state system
      x1 = PRJ1{j1,k};
      pdf1 = PDF1{j1,k};
         for j0=1:N0          % consider each nonfunctional 0-state system
         x0 = PRJ0{j0,k};            % Note that j0 labels 1-state systems 
         pdf0 = PDF0{j0,k};                
% ---------------------------------------------------- calculate overlapFN
         overlapFNk(j1,j0,k) = getOverlap(x0,pdf0,x1,pdf1); 
         end
      end
% ------------------------------------------------------------- consider U
      for ju=1:Nu                     % consider each unknown-state system
      xu = PRJu{ju,k};
      pdfu = PDFu{ju,k};
% --------------------------------------------------------------- UN pairs
         for j0=1:N0          % consider each nonfunctional 0-state system
         x0 = PRJ0{j0,k};            % Note that j0 labels 0-state systems
         pdf0 = PDF0{j0,k};                
% ---------------------------------------------------- calculate overlapUN
         overlapUNk(ju,j0,k) = getOverlap(x0,pdf0,xu,pdfu); 
         end
% --------------------------------------------------------------- UF pairs
         for j1=1:N1             % consider each functional 1-state system
         x1 = PRJ1{j1,k};            % Note that j1 labels 1-state systems
         pdf1 = PDF1{j1,k};                
% ---------------------------------------------------- calculate overlapUF
         overlapUFk(ju,j1,k) = getOverlap(x1,pdf1,xu,pdfu); 
         end
      end                                              
   end                                         % <---- end loop over modes
% -------------------------------------------------------- debugging check
% disp( ['N1 = ',num2str(N1),'  N0 = ',num2str(N0)] );
% disp( overlapFNk );
% disp('   ');
% disp( ['Nu = ',num2str(Nu),'  N0 = ',num2str(N0)] );
% disp( overlapUNk );
% disp('   ');
% disp( ['Nu = ',num2str(Nu),'  N1 = ',num2str(N1)] );
% disp( overlapUFk );
% error('stop here for now'); 
% ------------------------------------------------ calculate probabilities
pUk = zeros(Nu,nModes);
pUFk_F = zeros(Nu,N1,nModes);
wFk = zeros(N1,nModes);
pUFk = zeros(Nu,N1,nModes);
% ----------------------------
qUk = zeros(Nu,nModes);
qUNk_N = zeros(Nu,N0,nModes);
wNk = zeros(N0,nModes);
qUNk = zeros(Nu,N0,nModes);
% ----------------------------
wUk = zeros(Nu,nModes);
% ----------------------- mean probability for mode classification quality
ave_wFk = zeros(nModes,1);
ave_wNk = zeros(nModes,1);
% ------------------------------- moments for probability that U functions
pUfunct1 = zeros(Nu,1);
pUfunct2 = zeros(Nu,1);
sigFN = 0.75/sqrt( sqrt(Ns1*Ns0) );      % using 0.75 is somewhat generous
sigUF = 0.75/sqrt( sqrt(Nsu*Ns1) );
sigUN = 0.75/sqrt( sqrt(Nsu*Ns0) );
nTrials = 626; 
for nt=1:nTrials
   for k=1:nModes              % consider each mode: k-th vector-direction
% ---------------------------------- step 1: calculate prior probabilities
%      wFk = probability that the system labeled as F is distinct from all
%            nonfunctional systems in mode k.
%      wNk = probability that the system labeled as N is distinct from all
%            functional systems in mode k. 
% --------------------------------------------------------------- FN pairs
      for j1=1:N1                % consider each functional 1-state system
      wFk(j1,k) = 1;              % compare this 1-state with all 0-states
         for j0=1:N0          % consider each nonfunctional 0-state system
         xx = overlapFNk(j1,j0,k) + sigFN*randn;
         overlapFN = min(1, max(0,xx) );
         jFN = max(round(overlapFN/dxOVL),1);
        %wFk(j1,k) = wFk(j1,k)*( 1 - iWtOVL(jFN) ); % not sensitive enough
         wFk(j1,k) = wFk(j1,k)*dWtOVL(jFN);
%          disp(['trial= ',num2str(nt),'  k=',num2str(k), ...
%                '  j1=',num2str(j1),'  j0=',num2str(j0), ...
%                '  overlapFN=',num2str(overlapFN), ...
%                '  wFk=',num2str( wFk(j1,k) )]);
         end
%       disp('  ');                                        % for debugging
      end
   ave_wFk(k) = ave_wFk(k) + mean( wFk(:,k) );
      for j0=1:N0             % consider each nonfunctional 0-state system
      wNk(j0,k) = 1;              % compare this 0-state with all 1-states
         for j1=1:N1             % consider each functional 1-state system
         xx = overlapFNk(j1,j0,k) + sigFN*randn;
         overlapFN = min(1, max(0,xx) );
         jFN = max(round(overlapFN/dxOVL),1);
        %wNk(j0,k) = wNk(j0,k)*( 1 - iWtOVL(jFN) ); % not sensitive enough
         wNk(j0,k) = wNk(j0,k)*dWtOVL(jFN);
%          disp(['trial= ',num2str(nt),'  k=',num2str(k), ...
%                '  j0=',num2str(j0),'  j1=',num2str(j1), ...
%                '  overlapFN=',num2str(overlapFN), ...
%                '  wNk=',num2str( wNk(j0,k) )]);
         end
%       disp('  ');                                        % for debugging
      end
   ave_wNk(k) = ave_wNk(k) + mean( wNk(:,k) );
%    disp(['trial= ',num2str(nt),'  k=',num2str(k), ...    % for debugging
%          '  <wFk>=',num2str( mean( wFk(:,k) ) ), ... 
%          '  <wNk>=',num2str( mean( wNk(:,k) ) )]);
% ---------------------------- step 2: calculate conditional probabilities
%   pUFk_F = conditional probabilty that U mimicks F w.r.t. mode k
%            given that F is functional.
%   qUNk_N = conditional probabilty that U mimicks N w.r.t. mode k
%            given that N is nonfunctional. 
% ------------------------------------------------------------- consider U
      for ju=1:Nu                     % consider each unknown-state system
% --------------------------------------------------------------- UF pairs
         for j1=1:N1             % consider each functional 1-state system
% ---------------------------------------------------- calculate overlapUF
         xx = overlapUFk(ju,j1,k) + sigUF*randn;
         overlapUF = min(1, max(0,xx) );   
         jUF = max(round(overlapUF/dxOVL),1);
         pUFk_F(ju,j1,k) = 1 - dWtOVL(jUF);
%          disp(['trial= ',num2str(nt),'  k=',num2str(k), ...   % to debug
%                '  ju=',num2str(ju),'  j1=',num2str(j1), ...
%                '  overlapUF=',num2str(overlapUF), ...
%                '  pUFk_F=',num2str( pUFk_F(ju,j1,k) )]);
         end
%       disp('  ')                                         % for debugging
% --------------------------------------------------------------- UN pairs
         for j0=1:N0          % consider each nonfunctional 0-state system
% ---------------------------------------------------- calculate overlapUN
         xx = overlapUNk(ju,j0,k) + sigUN*randn;
         overlapUN = min(1, max(0,xx) );         
         jUN = max(round(overlapUN/dxOVL),1);
         qUNk_N(ju,j0,k) = 1 - dWtOVL(jUN);
%          disp(['trial= ',num2str(nt),'  k=',num2str(k), ...   % to debug
%                '  ju=',num2str(ju),'  j0=',num2str(j0), ...
%                '  overlapUN=',num2str(overlapUN), ...
%                '  qUNk_N=',num2str( qUNk_N(ju,j0,k) )]);
         end
%       disp( dividerLine('next ju') );                    % for debugging
      end
% ----------------------------- step 3: calculate similarity probabilities
%     pUFk = pUFk_F*wFk    => probability U functions like F w.r.t. mode k
%     qUNk = qUNk_N*wNk => probabilty U does not function like N in mode k
      for ju=1:Nu                     % consider each unknown-state system
         for j1=1:N1             % consider each functional 1-state system
         pUFk(ju,j1,k) = pUFk_F(ju,j1,k)*wFk(j1,k);
%          disp(['trial= ',num2str(nt),'  k=',num2str(k), ...   % to debug
%                '  ju=',num2str(ju),'  j1=',num2str(j1), ...
%                '  pUFk=',num2str( pUFk(ju,j1,k) )]);
         end
         for j0=1:N0          % consider each nonfunctional 0-state system
         qUNk(ju,j0,k) = qUNk_N(ju,j0,k)*wNk(j0,k);
%          disp(['trial= ',num2str(nt),'  k=',num2str(k), ...   % to debug
%                '  ju=',num2str(ju),'  j0=',num2str(j0), ...
%                '  qUNk=',num2str( qUNk(ju,j0,k) )]);
         end
%       disp( dividerLine('next ju') );                    % for debugging
      end
% --------------------- step 4a: calculate probability U is NOT functional
% 1 - pUFk = probability U is nonfunctional in mode k.
%  1 - pUk = product over all F of (1 - pUFk)
%            This product represents an AND operation for all F to be
%            distinct from U for non-function to take place. Note that
%            if U mimics only 1 F very well, nonfunction is null. 
%*-->  pUk = 1 - product over all F of (1 - pUFk)
      for ju=1:Nu                     % consider each unknown-state system
      temp = 1;
         for j1=1:N1             % consider each functional 1-state system
         temp = temp*(1 - pUFk(ju,j1,k) ); % AND => U does not look like F
         end
      pUk(ju,k) = 1 - temp;  % note: if temp=1 => U cannot be F => pUk = 0
      end
% ------------------ step 4b: calculate probability U is NOT nonfunctional
% 1 - qUNk = probability U is functional in mode k. 
%*-->  qUk = product over all N of (1 - qUNk)  
%            This product represents an AND operation for all N to be
%            distinct from U for function to take place. Note that
%            if U mimics only 1 N very well, function is null. 
      for ju=1:Nu                     % consider each unknown-state system
      temp = 1;
         for j0=1:N0          % consider each nonfunctional 0-state system
         temp = temp*( 1 - qUNk(ju,j0,k) );
         end
      qUk(ju,k) = temp; % note: if temp=1 => U is surely NOT nonfunctional
      end
% ------------ step 5: calculate probability U is functional w.r.t. mode k
%      wUk = pUk*qUk = probability U is functional in mode k.
%            This product represents an AND operation saying that U should
%            have characteristics that are similar or the same as those 
%            characteristics shared by systems that function while at the 
%            same time be as dissimilar from those characteristics found
%            in systems having nonfunction.

% % % %        wUk(:,k) = pUk(:,k).*qUk(:,k);                    % orginal
% % % %        wUk(:,k) = (0.99*pUk(:,k) + 0.01).*qUk(:,k);   % saturation
% ----------------------------------------------- another saturation model
%      Let pUk = tpF   and   qUk = tqN
%      tpF =          probability unknown system is functional
%      tqF = 1 - tpF = probability unknown system is not functional
%      tqN =          probability unknown system is not nonfunctional
%      tpN = 1 - tqN = probability unknown system is nonfunctional
%      (tpN*tpF + tqN) assigned likelihood of unknown system is functional
%            as a function of the probability that the unknown system is
%            is not nonfunctional. If qN --> 1, the system might function 
%            disregarding what is known about functional systems.
%            As qN --> 0 the unknown system is likely nonfunctional but
%            the knowledge of what is known about functional systems is
%            included.  This modifies pUk, which would normally be tpF.
%      (1 - tpN*tqF) = (tpN*tpF + tqN) identically. 
% ------------------------------------------------------------------------
       tqF = 1 - pUk(:,k);
       tpN = 1 - qUk(:,k);
       wUk(:,k) = (1 - tqF.*tpN).*qUk(:,k);    
       
%       for ju=1:Nu                                        % for debugging
%       disp(['trial= ',num2str(nt),'  k=',num2str(k), ...  
%             '  ju=',num2str(ju), ...
%             '  pUk=',num2str( pUk(ju,k) ), ... 
%             '  qUk=',num2str( qUk(ju,k) ), ...
%             '  wUk=',num2str( wUk(ju,k) )]);
%       end
%                                       for debugging on any section above
%    disp(dividerLine);
%    pause
%    disp('         ');
   end 
% ---------- step 6: calculate probability U functions as a geometric mean
%  pUFprod = product over all k of wUk 
%  pUfunct = (pUFprod)^(1/Nmodes)  => geometric mean of wUk, since this is
%            an AND operation. Reporting this geometrical mean removes the
%            relative number of discriminant modes as a factor.
%            For numerical stability: let wUk --> max(wUk,1.0e-125)
   for ju=1:Nu
% --------------------------------------------------------------- method A
   sum_ln_pUFprod = 0;
      for k=1:nModes
      ln_wUk = log( max( wUk(ju,k), 1.0e-125) );
      sum_ln_pUFprod = sum_ln_pUFprod + ln_wUk;
      end
   pUfunctA = exp( sum_ln_pUFprod/nModes );           % => exp( <ln_wUk> )
   
% --------------------------------------------------------------- method B
   pUfunctB = sqrt( mean( wUk(ju,:).^2 ) );       
   
% --------------------------------------------------------------- method C
%  pUfunct = hybridMethodFunction(pUfunctA,pUfunctB)
   pwA = log10(pUfunctA);                     % typically too conservative
   pwB = log10(pUfunctB);                            % typically too risky
   pwD = min(10*pUfunctB + (pwA - pwB), 0); %NOTE: all powers are negative
   pwD = pwD/(1 + 10*pUfunctB);      % variable scale factor for reduction
   pUfunct = pUfunctB*10^pwD;                   % pwD => reduction in risk
   %disp( [log10(pUfunctA),pUfunctB,pUfunct] );
   
   pUfunct1(ju) = pUfunct1(ju) + pUfunct;
   pUfunct2(ju) = pUfunct2(ju) + pUfunct*pUfunct;
   end
end
ave_wFk = ave_wFk/nTrials;
ave_wNk = ave_wNk/nTrials;
% ------------------------- note: using geometrical averages <wFk> = <wNk>
% REMARK: here a normal average was used for ave_wFk and ave_wNk
% since these would be equal had we used geometrical averages, and because
% both estimates are using different noise factors .... it makes sense to
% report only one number for mode quality, as the average of the averages.
wm = 0.5*(ave_wFk + ave_wNk);         % switch from k-th mode to m-th mode
pUfunct1 = pUfunct1/nTrials;
pUfunct2 = pUfunct2/nTrials;
sdPrbUfunct = sqrt( max(1.0e-125, pUfunct2 - pUfunct1.^2 ) );
% ------------------------------------------------------ debugging check 1
% disp('   ');
% disp( dividerLine('print: [wm,ave_wFm,ave_wNm]') );      % for debugging
% disp( [wm,ave_wFk,ave_wNk] );
% ------------------------------------------------------ debugging check 2
% disp('   ');
% disp( dividerLine('print: [pUfunct,sdPrbUfunct]') );
% format long
% disp( aMname );
% disp('  p1    sig_p1    -log10(p1) ');
% disp( [pUfunct1,sdPrbUfunct, -log10(pUfunct1)] );
% if( nModes > -1000 )
% error('stop here for now');
% end
%%                                         visualize the above calculation
padd = '%03i';
ndigits = min(1 + floor( log10(0.0001 + nModes) ),9);
padd(3) = num2str(ndigits);
if( verbosity > 1 )
  for ju=1:Nu
    for k=1:nModes
% ------------------------------------------------ determine range of plot
    xu = PRJu{ju,k}; 
    xmin = xu(1);  
    xmax = xu(end);
       for j1=1:N1
       minx1 = PRJ1{j1,k}(1);
       maxx1 = PRJ1{j1,k}(end);
       xmin = min(xmin,minx1);
       xmax = max(xmax,maxx1);
       end
% -------------------------------
       for j0=1:N0
       minx0 = PRJ0{j0,k}(1);
       maxx0 = PRJ0{j0,k}(end);
       xmin = min(xmin,minx0);
       xmax = max(xmax,maxx0);
       end
    pdfu = PDFu{ju,k};
       if( verbosity == 2 )
       H = figure('visible','off');
       else                                             % => verbosity = 3
       H = figure(figure_number0);
       end
    clf
    hold on;
% ------------------------------------------------------------------------
       for j0=1:N0
       x0 = PRJ0{j0,k};
       pdf0 = PDF0{j0,k};
       h0 = plot(x0,pdf0,'r','linewidth',1.8);
       end
% -------------------------------
       for j1=1:N1
       x1 = PRJ1{j1,k};
       pdf1 = PDF1{j1,k};
       h1 = plot(x1,pdf1,'b','linewidth',1.8);
       end
% -------------------------------
    hu = plot(xu,pdfu,'k','linewidth',2.5);
% --------------------------------------------------- set figure character
    xlabel('mode projection');
    ylabel('probability density');
    legend([h1,h0,hu],{'functional','nonfunctional','unclassified'});
    %legend([h1,h0],{'functional','nonfunctional'});
    temp0 = wm(k);                              % mean weight of k-th mode
    temp1 = pUk(ju,k);
    temp2 = qUk(ju,k);
    temp3 = wUk(ju,k); 
    title([aMname{ju},'d-mode=',num2str(k), ...
                      ' wm=',num2str(temp0,3), ...
                      ' p=',num2str(temp1,3), ...
                      ' q=',num2str(temp2,3), ...
                      ' pq=',num2str(temp3,3)], ...
                      'Interpreter','none','FontSize',12);
    str = ['_',aMname{ju},'_',num2str(k,padd)];
    fName = [fname,str];                          % overwrites old results
    fnamePlot = getOutputFileName(subFolder,fName);
    %disp(fName)
    saveas(H,fnamePlot,gvSPLOC.gFileType);    % FIX ME (remove comment)
       if( verbosity == 3 )
       commandwindow
       figure(H)
       pause
       end
    close(H)    
    end
%     if( verbosity == 3 )
%     commandwindow
%     sound(cos(1:7000).*sin(1:7000) + cos(1:7000));
%     pause(2);
%     commandwindow
%     disp('------------------------');
%     sound( sin(1:3000) )
%     ttt = input('enter stop to stop: ','s');   
%       if( strcmp(ttt,'stop') == 1 )
%       break;
%       end
%     else
%     break;
%     end
  end
end
%error('stop here for now'); 
%%                        package output into systemRanking data structure
ranking = struct;
% ------------------------------------------------------- single variables
ranking.rankType = rankType;
ranking.baseFname = fname;
ranking.nXsystems = Nu;
ranking.dataRefNameX = dataU.dataRefName;
ranking.dataRefNameF = dataF.dataRefName;
ranking.dataRefNameN = dataN.dataRefName;
ranking.nModes = nModes;
% ------------------------------------ cell arrays with 1 element per cell
ranking.systemXname1 = aMname;                               % data matrix
ranking.systemXname2 = aMname;                            % give same name
%%                                                  rank order the ranking
[~,indx] = sort(pUfunct1,'descend');
% --------------------------------------
pFave = pUfunct1(indx);
pFstd = sdPrbUfunct(indx);
% ------------------------
ranking.pFave = pFave;
ranking.pFstd = pFstd;
% ------------------------
ranking.wm = wm;
% ------------------------
ranking.pXF = pUk(indx,:);      % probability X is functional based on {F}
ranking.qXN = qUk(indx,:);     % prob X is not non-functional based on {N}
ranking.pqX = wUk(indx,:);                   % probability X is functional
% ---------------------------------------------------------- single matrix
ranking.U = U;                % a collection of vectors possibly dependent
%%                                      change name based on rank ordering
   for j=1:Nu
   j0 = indx(j);
   ranking.systemXname1{j} = aMname{j0}; 
   ranking.systemXname2{j} = aMname{j0};                   % use same name
   end
%%                                                            create table
aMname = ranking.systemXname1;       % => this puts aMname in sorted order
log10inv_pFave = -log10(pFave);
j = 1:Nu;
T = table(j',aMname',pFave,pFstd,log10inv_pFave, ...
         'VariableNames', ...
         {'rank','aMatrixName','pFave','pFstd','log10inv_pFave'});
% ------------------------------------------------------- sort locally too
   if( verbosity > 0 )
   fclose(fid);
   end
% format short;
% disp(T);
% error('stop here for now');
%%                                         record action in sploc log file
nM = nModes;                                       % nM is a temp variable
maxP = max(pFave);
aveP = mean(pFave);
minP = min(pFave);
%disp( splocLogFile );                                     % for debugging
fid = fopen(splocLogFile,'a');                    % append new information
fprintf(fid,'%s \n','  ');
fprintf(fid,'%s \n',[mfilename,'()']);
msg = dividerLine('basic summary');
fprintf(fid,'%s \n',msg);
fprintf(fid,'%s \n',['  classifying method: rankType = ',rankType]);
fprintf(fid,'%s \n',['                reference name = ',fname]);
fprintf(fid,'%s \n',['number of X-systems classified = ',num2str(Nu)]);
fprintf(fid,'%s \n',['   number of example 1-systems = ',num2str(N1)]);
fprintf(fid,'%s \n',['   number of example 0-systems = ',num2str(N0)]);
fprintf(fid,'%s \n',['  number of discriminant modes = ',num2str(nM)]);
fprintf(fid,'%s \n',['           minimum probability = ',num2str(minP)]);
fprintf(fid,'%s \n',['           average probability = ',num2str(aveP)]);
fprintf(fid,'%s \n',['           maximum probability = ',num2str(maxP)]);
fprintf(fid,'%s \n',dividerLine);
fclose(fid);
end
