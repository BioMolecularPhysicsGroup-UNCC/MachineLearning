function rmsip = getRMSIP(U,V)
% calculation of RMSIP between basis sets
% use for general purposes or as a simple tool for the SPLOC toolbox
% 
% INPUT 
% U is a nDOF x nU basis set matrix  where nU = # of modes in the U matrix
% V is a nDOF x nV basis set matrix  where nV = # of modes in the V matrix
%
% PROCESS 
% calculate all pairs of inner products and average their squares. 
% take the square root of the mean of the squares. 
%
% OUTPUT
% rmsip
%%                                                      do the calculation
% [nDOF,nU] = size(U);
% [~,nV] = size(V);
allInnerProductsSq = (U'*V).^2;                   % gives a nU x nV matrix
rmsip = sqrt( mean( sum(allInnerProductsSq,2) ) );
end
