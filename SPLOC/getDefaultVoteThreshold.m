function voteThreshold = getDefaultVoteThreshold(nS0,nS1,nDOF)
% applies a heristic formula to calculate the vote threshold
%
% INPUT
% nS0 =         precalculated as: nD0/sqrt(n0) = sampleSize*sqrt(n0)
% nS1 =         precalculated as: nD1/sqrt(n1) = sampleSize*sqrt(n1)
% REMARKS:      n0 = number of systems of type 0
%               n1 = number of systems of type 1
%               nD0 = sampleSize*n0 => system 0 sampling statistics
%               nD1 = sampleSize*n1 => system 1 sampling statistics
%               Assumption is that sampleSize is the same for all systems
%               being compared. This is important to do to remove biases.
%
% OUTPUT
% voteThresold = vote threshold used as default if not user-specified
%
% ----------------------------------------------------- start calculations
% nD1 = n1*samplesize  and   nD0 = n0*sampleSize
% nS0 = nD0/sqrt(n0);           % nS0 = sqrt(n0)*sampleSize = nD0/sqrt(n0)
% nS1 = nD1/sqrt(n1);           % nS1 = sqrt(n1)*sampleSize = nD1/sqrt(n1)
nSmin = min(nS1,nS0);
nSmax = max(nS1,nS0);
Neff = 0.5*( nSmin + nSmin* ( 2*nSmax/(nSmin + nSmax) )^2 );
pw = max(0,0.5 - 0.49999*sqrt(nDOF/Neff) );
voteThreshold = min(0.95,0.5 + 0.7*(nDOF/Neff)^pw);
end

