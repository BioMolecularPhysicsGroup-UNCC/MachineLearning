function newBVS = randomizeSBV(BVS,subspace,A)
% Given a Basis Vector Set (BVS): randomize a subspace to yeild a new BVS
%
% INPUT
% BVS      = nV x nV matrix defining a complete basis vector set
% subspace = integer array defining a set of m vectors to randomly rotate
% A        = amplitude used in Cayley transformation formula.
%            Generally 0.001 < A < 1 is a reasonable range. 
%            A = 0.001 => mean angle shifts ~ 0.05 degrees
%            A = 0.005 => mean angle shifts ~ 1 to 2 degrees
%            A = 0.01  => mean angle shifts ~ 2 to 3 degrees
%            A = 0.05  => mean angle shifts ~ 5 to 30 degrees   nV < 20
%                                             10 to 20 degrees  nV > 20
%            A = 0.2   => mean angle shifts ~ 10 to 70 degrees  nV < 20
%                                             30 to 65 degrees  nV > 20
%            A = 1     => mean angle shifts ~ 20 to 160 degrees nV < 20
%                                             60 to 120 degrees nV > 20
%                                             75 to 105 degrees nV > 100
%                                             85 to 95  degrees nV > 1000
% In general as nV increases the distribution gets sharper and lower and 
% upper limits shrink about the average more. However, there is overall
% consistency in the range of angle deviations independent of nV.
%
% PROCESS 
% Apply a series of Cayley rotations to yield a randomized rotation within
% a specified subspace hold all other vectors outside the subspace fixed.
% The degree of randomization is controled by aa.
%
% OUTPUT
% newBVS = nV x nV matrix defining final basis vector set
%%                                                               set input
[nV,mV] = size(BVS);
   if( nV ~= mV )
   error('BVS is not a square matrix');
   end
% ----------------------------------------------
flagSubspace = 1;                                     % subspace specified
   if( nargin < 2 )               
   flagSubspace = -1;                              % no subspace specified
   end
% ----------------------------------------------
   if( nargin == 1 )
      if( nV > 1 )
      U = BVS;
      else
      newBVS = BVS;
      return
      end
   else
   m = length(subspace);
      if( m > nV )
      error('subspace list exeeds dimension of space');
      elseif( m > 1 )
      smax = max(subspace);
         if( smax > nV )
         error('vector # in subspace list exeeds maximum dimension');
         end
      smin = min(subspace);
         if( smin < 1 )
         error('vector # in subspace is less than 1');
         end
      m2 = length( unique(subspace) );
         if( m ~= m2 )
         error('not all vector # in subspace is are unique');
         end
      U = BVS(:,subspace);
      else
      newBVS = BVS;
      return
      end
   end
% -------------------------------------------------------- set amplitude A
   if( nargin ~= 3 )
   A = 0.04;                 % <= default: angle range is between 5 and 15
   else
   A = max(min(1,A),0.001);
   end
%%                                      apply a series of Cayley Rotations
D0 = min(m,10); 
Id = diag( ones(1,D0) );
   for jj=1:4
   indexREF = randperm(m);
   jLower = 1;
   jUpper = D0;
   nRmax = m - D0;
      for nR=0:D0:nRmax
      index = indexREF(jLower:jUpper);
      U0 = U(:,index);
% ---------------------------------------------- construct rotation matrix
      S0 = A*rand(D0);
      S = S0 - S0';                            % => S is now a skew matrix
      R = (Id - S)/(Id + S);         % => random rotation matrix in D0 DIM
      newU0 = U0*R;
      U(:,index) = newU0;
      jLower = jLower + D0;
      jUpper = jUpper + D0;
      end
   end
%%                                                          package output
   if( flagSubspace < 0 )         % subspace not specified, use full space
   newBVS = U;
   else                                               % subspace specified
   newBVS = BVS;
   newBVS(:,subspace) = U;
   end
end

