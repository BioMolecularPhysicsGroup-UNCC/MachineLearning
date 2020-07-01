function rmsipSequence = getRMSIPsequence(U,V)
% yeilds sequence of root mean squared inner products as DIM(U) is scanned 
m = length(U);
rmsipSequence = zeros(1,m);
   for k=1:m
   subset = 1:k;
   rmsipSequence(k) = getRMSIP(U(:,subset),V);
   end
end 
