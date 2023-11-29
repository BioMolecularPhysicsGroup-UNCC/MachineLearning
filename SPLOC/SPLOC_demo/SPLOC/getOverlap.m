function overlap = getOverlap(x1,pdf1,x2,pdf2)
% calculates overlap between two probability density functions
% overlap = [2 - integral |pdf2 - pdf1| *dx ]/2
xmin1 = min(x1);
xmin2 = min(x2);
xmax1 = max(x1);
xmax2 = max(x2);
xmin = min(xmin1,xmin2);
xmax = max(xmax1,xmax2);
dx = (xmax - xmin)/200;
x = xmin:dx:xmax;
p1 = zeros( size(x) );
p2 = zeros( size(x) );
% ------------------------------------------------- create new x1 and pdf1
   if( xmin < xmin1 )
   L = (x < xmin1 );
   xL = x(L);
   pL = p1(L);
   pdf1 = [pL,pdf1];
   x1 = [xL,x1];
   end
% --------------------------
   if( xmax > xmax1 )
   L = (x > xmax1 );
   xR = x(L);
   pR = p1(L);
   pdf1 = [pdf1,pR];
   x1 = [x1,xR];
   end
% ------------------------------------------------- create new x2 and pdf2
   if( xmin < xmin2 )
   L = (x < xmin2 );
   xL = x(L);
   pL = p2(L);
   pdf2 = [pL,pdf2];
   x2 = [xL,x2];
   end
% --------------------------
   if( xmax > xmax2 )
   L = (x > xmax2 );
   xR = x(L);
   pR = p2(L);
   pdf2 = [pdf2,pR];
   x2 = [x2,xR];
   end
% ------------------------------------------------------------ interpolate
f1 = interp1(x1,pdf1,x);
f2 = interp1(x2,pdf2,x);
% --------------------------------------------- calculate overlap integral
temp = dx*sum( abs(f2 - f1) );
overlap = 0.5*(2 - temp);
   if( overlap > 0.999 )
   overlap = 1;
   end
% -----------------------
   if( overlap < 0.001 )
   overlap = 0.0;
   end
end
