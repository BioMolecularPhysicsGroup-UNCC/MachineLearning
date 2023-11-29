function [splocResults,varargout] = ...
         getBasisVecSpectrum(U,fname,trait1,traitORcID,vT)
% Given a basis set of vectors determine its congruency spectrum
% Different types of inputs and outputs are possible
%
% ------------------ three conditions -----------------       congruency
% score > minScore1  &  dVote > vT  &  dQuality > minQF  =>  discriminant
% score < maxScore0  &  iVote > vT  &  iQuality > minQF  =>  indifference
%                otherwise                               =>  undetermined
%
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
% U = proposed set of basis vectors as a [nDOF,nModes] size matrix
% fname = baseFname => all output file names spawn off this base file name
%
% ------------------------------------------------------------ from trait1
% n1 = # of statistical metrics for Functional or labeled systems
% M1 = cell array containing n1 mean vectors for Functional systems
% Q1 = cell array containing n1 covariance matrices for Functional systems
% nV1 = # of variables in each Functional system
% nD1 = total # of data samples = n1*samplesize1 for Functional systems
% 
% ------------------------------------------------------------ from trait0
% n0 = # of statistical metrics for Nonfunctional or unlabeled systems
% M0 = cell array containing n0 mean vectors for Nonfunctional systems
% Q0 = cell array containing n0 covariance matrices for Functional systems
% nV0 = # of variables in each Nonfunctional system
% nD0 = total # of data samples = n0*samplesize0 for Nonfunctional systems
%
% cID = class indentifiers or labels for multiclass case
% 
% vT = optional vote threshold   (default is calculated based on sampling)
%
% USAGE:
%
%   splocResults = getBasisVecSpectrum(U,fname,trait1,trait0)
%   splocResults = getBasisVecSpectrum(U,fname,trait1,trait0,vT)
%           Notes: trait1 defines functional systems only    => class = 1
%                  trait0 defines nonfunctional systems only => class = 0
% 
% mcsplocResults = getBasisVecSpectrum(U,fname,traitL,cID)
% mcsplocResults = getBasisVecSpectrum(U,fname,traitL,cID,vT)
%           Notes: traitL defines classified systems => class labels = cID
%                  cID = list of numbers from 1 to nC = # of class labels
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - Optional output
%  [splocResults, mapU2S] = getBasisVecSpectrum( ... )
%  [splocResults, mapU2S, mapS2U] = getBasisVecSpectrum( ... )
%
%  REMARKS:    mapU2S =>    Sindex = mapU2S(Uindex)
%              mapS2U =>    Uindex = mapS2U(Sindex)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Meaning of function names:
%   sploc ----> supervised projective learning for orthogonal congruences
% mcsploc ----> multiclass-SPLOC
% 
% List of distinct data structures that are OUTPUT for this program
%
%   splocResults => generalized data structure for sploc or mcsploc
%
% PROCESS
% Emmulate the process found in sploc() or mcsploc(), except this is a one
% time shot that quantifies the input vectors supplied. No optimization is
% performed. 
%
%
% OUTPUT (datastructure)    
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
%
% OPTIONAL output
%
%            mapU2S = map original U indices to SPLOC mode indices
%                  => Sindex = mapU2S(Uindex)
%            mapS2U = map SPLOC mode indices to original U indices
%                  => Uindex = mapS2U(Sindex)
iFormat = -1;                      % => incorrect format for this function
%%                                               parse input into function
%  iFormat = 1  => getBasisVecSpectrum(U,fname,trait1,trait0)            1
%  iFormat = 2  => getBasisVecSpectrum(U,fname,trait1,trait0,vT)         2
% ------------------------------------------------------------------------
%  iFormat = 3  => getBasisVecSpectrum(U,fname,traitL,cID)               3
%  iFormat = 4  => getBasisVecSpectrum(U,fname,traitL,cID,vT)            4
%%                            first only check for consistent input format
   switch nargin
       case 4
          if( isnumeric(U) && ischar(fname) && ...
              ~isnumeric(trait1) && isnumeric(traitORcID) )
          iFormat = 3;
          elseif( isnumeric(U) && ischar(fname) && ...
              ~isnumeric(trait1) && ~isnumeric(traitORcID) )
          iFormat = 1;
          end
       case 5
          if( isnumeric(U) && ischar(fname) && ...
              ~isnumeric(trait1) && isnumeric(traitORcID) && ...
               isnumeric(vT) )
          iFormat = 4;
          elseif( isnumeric(U) && ischar(fname) && ...
              ~isnumeric(trait1) && ~isnumeric(traitORcID) && ...
               isnumeric(vT) )
          iFormat = 2;
          end       
   end
  if( iFormat < 0 )
  disp('  ');
  disp('USAGE options:')
  disp('--------------------------------------------------------------');
  disp('  sploc: getBasisVecSpectrum(U,fname,trait1,trait0)');          %1
  disp('         getBasisVecSpectrum(U,fname,trait1,trait0,vT)');       %2
  disp('--------------------------------------------------------------');
  disp('mcsploc: getBasisVecSpectrum(U,fname,traitL,cID)');             %3
  disp('         getBasisVecSpectrum(U,fname,traitL,cID,vT)');          %4
  disp('--------------------------------------------------------------');
  error('wrong input format');
  end
% -------------------- distribute to proper function based on format cases
  switch iFormat
    case 1
    trait0 = traitORcID;
    [splocResults,mapU2S,mapS2U] = ...
    getBVStwoClasses(U,fname,trait1,trait0);
    case 2
    trait0 = traitORcID;
    [splocResults,mapU2S,mapS2U] = ...
    getBVStwoClasses(U,fname,trait1,trait0,vT);
    case 3
    cID = traitORcID;
    L = ( cID < 1 );
       if( sum(L) > 0 )
       error('class labels must be positive integers');
       end
      % ------
       try
       [splocResults,mapU2S,mapS2U] = ...
       getBVSmultiClass(U,fname,trait1,cID);
       catch me
       disp('   ');
       disp('multiclass SPLOC is not ready to be released to public');
       disp('contact Dr. Jacobs by Email for more information');
       disp(' Email: djacobs1@uncc.edu')
       error('Implementation is still in development');
       end  
    case 4
    cID = traitORcID;
    L = ( cID < 1 );
       if( sum(L) > 0 )
       error('class labels must be positive integers');
       end 
      % -----
       try
       [splocResults,mapU2S,mapS2U] = ...
       getBVSmultiClass(U,fname,trait1,cID,vT);
       catch me
       disp('   ');
       disp('multiclass SPLOC is not ready to be released to public');
       disp('contact Dr. Jacobs by Email for more information');
       disp(' Email: djacobs1@uncc.edu')
       error('Implementation is still in development');
       end 
  end 
%%                                       deal with variable output request
   switch nargout
       case 2
       varargout{1} = mapU2S;
       case 3
       varargout{1} = mapU2S;
       varargout{2} = mapS2U;
   end
end
