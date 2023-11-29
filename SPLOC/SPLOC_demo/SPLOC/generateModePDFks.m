function [modePROJ,modePDF] = generateModePDFks(A,basisVectors)
% Project data matrix onto basis vectors to construct PDF per mode.
%
% INPUT
% A = data matrix
% basisVectors = list of column vectors to generate projections
%
% PROCESS
% Using each vector in basis set calculate PDF as a projection.
% 
% OUTPUT
% modePROJ{:} = projection of data matrix to define independent variable
%  modePDF{:} = probablity density function for the projected data
%         {:} => per mode
% ------------------------------------------------------------------------
[nDOF,nVec] = size(basisVectors);
[iDOF,~] = size(A);
   if( iDOF ~= nDOF )
   error('inconsistent number of DOF with basis vectors and data matrix');
   end
modePROJ = cell(1,nVec);
modePDF =  cell(1,nVec);
   for k=1:nVec
   V = basisVectors(:,k)';
   yU = V*A;
   max_y = max(yU);
   min_y = min(yU);
% -----------------------------------------------------------------------
   dy = 0.05*(max_y - min_y);
   max_y = max_y + dy/2;
   min_y = min_y - dy/2;
   dy = (max_y - min_y)/200;
   pts = min_y:dy:max_y;   
   [modePDF{k},modePROJ{k}] = ksdensity(yU,pts);   
   end
end
