function U = getDiscriminantSBV(splocResults)
% extracts entire discriminant selection basis vectors from splocResults
% 
% INPUT splocResults.   <-- data structure 
% baseFname = base file name useful for spawning more file names
%       SBV = selection basis vectors
%       SEV = selection eigenvalues        reporting: (average of logs)
%       CEV = consensus eigenvalues
%       CIP = congruency indicators (1,0,-1)
%        vT = voting threshold used to establish consensus
% rankScore = score that ranks multiple solutions based on consistency
%             -1 => nonexistent score: no ensemble to rank against
% 
% OUTPUT
% U = selection basis vectors with CIP = 1 => discriminant. The size of
%     the U matrix is nV x nMode, where nV = dimension of the vector space
%     and nMode is the cumulative number of modes collected that span the
%     discriminant subspace in addition to horizontal concatenation across
%     the ensemble of sploc results.
%   = [] when the discriminant subspace is empty. 
%%                                          apply horizontal concatenation
I = 1;                               % indicator for discriminant subspace
% -------------------------------- remaining code is subspace type neutral
CIP = splocResults.Cind;
Ldi = (CIP == 2);                    % identifies dual-purpose (d&i)-modes
CIP(Ldi) = I;              % includes d-modes and dual-purpose (d&i)-modes
L = (CIP == I);
nModes = sum(L);
   if( nModes == 0 )
   U = [];
   else
   U = splocResults.SBV(:,L);
   end
end