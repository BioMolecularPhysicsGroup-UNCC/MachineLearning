%qType = 'hardLinear';
qType = 'softGeometry';
% ------------------------------------------------------------------------
for i=1:20
    
    
   switch qType
     case 'hardLinear'
     doHardLinear = true;
     doSoftGeometry = false;
     case 'softGeometry'
        if( rand < 0.2 )
        doHardLinear = true;
        else
        doHardLinear = false;
        end
     doSoftGeometry = true;
     otherwise
     error('unknown clustering property measure');
   end
% ------------------------------------------------------------------------
   if( doHardLinear )
   disp('do 1 ');
   end
% ------------------------------------------------------------------------
   if( doSoftGeometry )
      if( doHardLinear )
      disp('do 2 after 1');
      else
      disp('do 2 ');
      end
   end
   
   
end

% % a = 1;
% % b = true;
% % switch a
% %     case 1 
% %     
% %     case 2
% %        if( b )
% %        disp('do 1 when b is true');
% %        end
% % end
% % % ------------------------------------------------------------------------
% % switch a
% %     case 2
% %     disp('do 2 ');
% % end