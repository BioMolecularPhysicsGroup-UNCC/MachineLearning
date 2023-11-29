function reducedC = reduceCmatrix(C,pType,dim)
% reduces a matrix based on vectorized data to one without vectorization
% Utility function
% -------------------------------------------------------- local variables 
% INPUT
% C = input matrix
% pType = packing type (format) that is applied to the vector space
% dim = # of components in a local vector (e.g. x,y,z for 1 atom)  
%
% PROCESS
% reduces a (dimxN) by (dimxN) matrix down to a NxN matrix by grouping 
% corresponding vectorized components. For instance, the packing in the 
% input Cmatrix can be of the form xyz-xyz-xyz or xxx-yyy-zzz. At the end,
% the reduced matrix is really of the form: xx + yy + zz per element. This 
% function generalizes the number of components in a vector from 3 to dim.
%
% OUTPUT  
% reducedC = the reduced matrix. 
%%                                                      simple error check
[n,m] = size(C);
   if( n ~= m )
   error('C is not square: method only works for square matrices');
   end
%%                                                          calculate RMSF
   switch pType
      case 'notVectored'
      reducedC = C;                         % easy case: requires no work!
      case 'xyz-xyz-xyz'
      n = round(m/dim); 
      reducedC = zeros(n);
      reduceMap = zeros(1,m);
      ii1 = 1:dim:(m-1);
      i1 = 1:n;
         for k=1:dim
         ii = ii1 + (k-1);
         reduceMap(ii) = i1;
         end
% -------------------------------
         for ii=1:m
         i = reduceMap(ii);
         ki = ii - (i-1)*dim;
            for jj=1:m
            j = reduceMap(jj);
            kj = jj - (j-1)*dim;
               if( ki == kj )
               reducedC(i,j) = reducedC(i,j) + C(ii,jj);
%                disp(['index ',num2str([ii,jj]),' --> [', ...
%                               num2str(i),': ',num2str(ki),'  ', ...
%                               num2str(j),': ',num2str(kj),']']);
               end
            end
         end
      case 'xxx-yyy-zzz'
      n = round(m/dim);
      reducedC = zeros(n);
      reduceMap = zeros(1,m); 
      i1 = 1:n;
         for k=1:dim
         ii = i1 + (k - 1)*n;
         reduceMap(ii) = i1;
         end
% -------------------------------      
         for i=1:n
             for j=1:n
                for k=1:dim
                ii = i + (k-1)*n;
                jj = j + (k-1)*n;
                reducedC(i,j) = reducedC(i,j) + C(ii,jj);
%                 disp(['index ',num2str([ii,jj]),' --> [', ...
%                                num2str(i),': ',num2str(k),'  ', ...
%                                num2str(j),': ',num2str(k),']']);
                end
             end
         end
      otherwise
      error('pType is unknown');
   end
end