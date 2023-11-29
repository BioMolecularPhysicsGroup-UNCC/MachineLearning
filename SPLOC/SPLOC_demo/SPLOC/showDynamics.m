function Xshow = showDynamics(X,U)
% From X trajectory show simulated dynamics with an option to filter on U
% 
% INPUT:
% X = conformations tracked by a simulation of interest
% U = basis vectors for a subspace of interest
% 
% PROCESS:
% apply a uniform sampling of conformations from X in such a way that it 
% looks like a periodic system.
%
% OUTPUT
% Xshow = time dynamics for the moleculr motion, perhaps filtered by U
%
%%                                          select subset of conformations
[~,nt] = size(X);
   if( nt < 40 )
   dn = 1;
   elseif( nt < 200 )
   dn = 2 + round(nt/100);
   elseif( nt < 2000 )
   dn = 10 + round(nt/100);
   elseif( nt < 20000 )
   dn = 25 + round(nt/200);
   elseif( nt < 40000 )
   dn = 40 + round(nt/250);
   elseif( nt < 80000 )
   dn = 100 + round(nt/400);
   elseif( nt < 400000 )
   dn = 250 + round(nt/500);
   else
   dn = 1000 + round(nt/1000);
   end
t1 = 1:dn:nt;
t2 = (nt + round(dn/2) ):dn:(2*nt);
t = [t1,t2];
index = nt:-1:1;
Xrev = X(:,index);
Xtemp = [X,Xrev];
clear X Xrev
Xshow = Xtemp(:,t);
clear Xtemp
%%                                             filter out unwanted motions
   if( nargin == 2 )
   aveX = mean(Xshow,2);
   A = Xshow - aveX;
   MSD = mean(A.*A);
   [~,indx] = min(MSD);
   Xref = Xshow(:,indx);           % use this frame as reference structure
   dv = Xshow - Xref;     % displacement vector w.r.t. to reference struct
   P = U*U';
   dv = P*dv;
   Xshow = Xref + dv;
   end
end

