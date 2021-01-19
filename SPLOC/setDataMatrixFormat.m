function matrixFormat = setDataMatrixFormat(pType,dim)
% defines the way a data matrix is formated in terms of vector components
%
% INPUT
% ---- mode 1:
% pType = string => packing type = {xyz-xyz-xyz, xxx-yyy-zzz, notVectored}
% dim   = # of components in local vector property.
% ---- mode 2:
% mFormat = data stucture => mFormat.pType  and  mFormat.dim 
%                            mFormat.pType = packing type
%                            mFormat.dim = # of components of local vector
%
% USAGE: matrixFormat = setDataMatrixFormat()
%        matrixFormat = setDataMatrixFormat('xyz-xyz-xyz',3)
%        matrixFormat = setDataMatrixFormat(mFormat)
%
% PROCESS 
% queries user about the characteristics of the vectorized data. Checks
% if input is admissible. 
%
% OUTPUT: data structure
% matrixFormat.pType = {xyz-xyz-xyz, xxx-yyy-zzz, notVectored};
% matrixFormat.dim   = # of components in local vector property
%%                                                transfer input to output
  matrixFormat = struct;
  if( nargin < 1 )                                  % user query is needed
% ------------------------------------------------------------------------
  disp('  ');
  disp('------------------------ allowed options for data matrix format');
  disp(' 1. not vectorized');
  disp(' 2. vectorized with confined packing format. (e.g. xyz-xyz-xyz)');
  disp(' 3. vectorized with extended packing format. (e.g. xxx-yyy-zzz)');
  disp('  ');
  pckType = input('    Enter option (1 through 3) for packing format: ');
     switch pckType
     case 1
     pType = 'notVectored';
     case 2
     pType = 'xyz-xyz-xyz';
     case 3
     pType = 'xxx-yyy-zzz';
     otherwise
     error('Option selected is outside allowed range');
     end
  disp('  ');
    if( pckType ~= 1 )
    dim = input('    Enter number of components in vector: ');
       if( dim < 2 )
       error(['dim = ',num2str(dim),' => Not vectorized data']);
       elseif( dim > 9999 )
       error(['dim = ',num2str(dim),' => currently too large to handle']);
       end
    else
    dim = 1;
    end
  matrixFormat.pType = pType;
  matrixFormat.dim = dim;
  return;
  elseif( nargin == 1 )
      if( isstruct(pType) )
      dim = pType.dim;
      temp = pType.pType;
      pType = temp;
      elseif( ischar(pType) )
         if( strcmp(pType,'notVectored') )
         dim = 1;
         else
         error('expecting:  pType,dim');
         end
      else
      error('expecting:  pType,dim');
      end
  end
%%                                                     % apply error check 
    switch pType
      case 'xyz-xyz-xyz'
      matrixFormat.pType = pType;
      case 'xxx-yyy-zzz'
      matrixFormat.pType = pType;
      case 'notVectored'
      matrixFormat.pType = pType;
      otherwise
    error('unknown format');
    end
    if( strcmp(pType,'notVectored') == 1 )
    matrixFormat.dim = 1;
    else
       if( dim < 2 )
       error(['dim = ',num2str(dim),' => Not vectorized data']);
       elseif( dim > 9999 )
       error(['dim = ',num2str(dim),' => currently too large to handle']);
       end
    matrixFormat.dim = dim;
    end
end