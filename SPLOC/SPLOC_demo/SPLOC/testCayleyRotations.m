        A = input(' Enter A: ');
Du = input('       Enter Du: ');
% Conclusion: Use second method for all Du since the first method is 
%             faster for Du < ~20 and approximately equal until Du = ~40
%             After this, the second method is hands down faster.
%             but, even when the first method is the fastest relative to
%             to the second for Du ~10 .... it is ~6 times faster. But, at
%             this point both methods are sufficiently fast that a factor
%             of 6 will not matter, plus the second method has better 
%             distribution of angles, which is similar for all Du. The
%             first method has variable distributions for different Du.
% ------------------------------------------------------------------------
Ntries = ceil(50000/Du);
%Ntries = input('Enter Ntries: '); 
% ------------------------------------------------- for small Du: Du < 100
flagMethod1 = -1;
%{
if( Du < 200 )                               % quit before going too large
flagMethod1 = 1;
V = diag( ones(1,Du) );
t1 = cputime;
% ---------------------- error check
   if( Du < 2 )
   error('Du must be 2 or greater');
   end

%A = A/sqrt(Du);           no good
% A = 0.1*A/sqrt( sqrt(Du) );
A = max(0.2,1.2 - Du/100 );
   for i=1:Ntries
   Id = diag( ones(1,Du) );
   S0 = A*rand(Du);
   S = S0 - S0';                          % => S is now a skew matrix
   R = (Id - S)/(Id + S);       % => random rotation matrix in Du DIM
   newV = V*R;                           % 100x7 = (100x7) x (7x7)
% %    M = newV'*V;                          % 7x7   = (7x100) x (100x7)
% %    anglesM = diag( acosd(M) );
   
% %    disp( anglesM )
% %    cosAngles = diag(M);
% %    angles = acosd( cosAngles );
% %    disp( sort(angles') );
% %    
% %    M2 = V'*newV;
% %    cosAngles = diag(M2);
% %    angles = acosd( cosAngles );
% %    disp( sort(angles') );

%    figure(1);
%    clf;
%    histogram( anglesM(:) );
%    title(['test # = ',num2str(i),'  A = ',num2str(A)]);
%    pause(0.5)
   
   V = newV;
   end
t1 = cputime - t1;
disp(['small Du: Ntries = ',num2str(Ntries), ...
      '  CPUseconds = ',num2str(t1)]);
end
%}
% ------------------------------------------------- for large Du: Du > 100
%{
t2 = cputime;
D0 = min(Du,10);
V = diag( ones(1,Du) );
itest = Ntries*Du;
   if( itest > 2000000 )
   error('Ntries*Du is too large to record data');
   else
   anglesM = zeros(Ntries,Du);
   cosAnglesM = zeros(1,Du);
   end
%A = 0.2;                                           % good for 5 rotations
A = 0.05;                    % for tickle --> small amplitude random noise
   for i=1:Ntries
   Vo = V;
   
   %for jj=1:5                                               % 5 rotations
   %for jj=1:1                                                % 1 rotation
   for jj=1:25
   %for jj=1:100                                           % 100 rotations
   %for jj=1:400                                            % 400 rotations
   indexREF = randperm(Du);
   jLower = 1;
   jUpper = D0;
   nRmax = Du - D0;
      for nR=0:D0:nRmax
      index = indexREF(jLower:jUpper);
      V0 = V(:,index);
% ---------------------------------------------- construct rotation matrix
      Id = diag( ones(1,D0) );
      S0 = A*rand(D0);
      S = S0 - S0';                            % => S is now a skew matrix
      R = (Id - S)/(Id + S);         % => random rotation matrix in D0 DIM
      newV0 = V0*R;
      V(:,index) = newV0;
      jLower = jLower + D0;
      jUpper = jUpper + D0;
      end
   end
%  % REMARK: M = V'*Vo;                        % 7x7   = (7x100) x (100x7)
%            From M, we need to check angles along the diagonal
% ----------------------------------------------------------- for checking
   
      for k=1:Du
      cosAnglesM(k) = sum( V(:,k).*Vo(:,k) );
      end
   anglesM(i,:) = acosd(0.999999999*cosAnglesM );
   end
t2 = cputime - t2;
disp(['large Du: Ntries = ',num2str(Ntries), ...
      '  CPUseconds = ',num2str(t2)]);
  if( flagMethod1 > 0 )
  disp(['ratio of large-Du/small-Du methods: ratio = ',num2str(t2/t1)]);
  end
figure(1);
   clf;
   histogram( anglesM(:) );
   title(['large Du test: Ntries = ',num2str(Ntries)]);
   pause(0.5)
%}
%%                                           test function  randomizeSBV.m
% ------------------------------------------------- for large Du: Du > 100
% {
t2 = cputime;
Vo = diag( ones(1,Du) );
itest = Ntries*Du;
   if( itest > 2000000 )
   error('Ntries*Du is too large to record data');
   else
   anglesM = zeros(Ntries,Du);
   cosAnglesM = zeros(1,Du);
   end
% ---------------------------------
%m = ceil( sqrt(Du) );
m = ceil(0.2*Du);
   for i=1:Ntries
   %V = randomizeSBV(Vo,subspace,A);
   %V = randomizeSBV(Vo,subspace);
   subspace = randperm(Du,m);
   V = randomizeSBV(Vo,subspace);
% ----------------------------------------------------------- for checking
      for k=1:Du
      cosAnglesM(k) = sum( V(:,k).*Vo(:,k) );
      end
   anglesM(i,:) = acosd(0.999999999*cosAnglesM );
   end
t2 = cputime - t2;
disp(['large Du: Ntries = ',num2str(Ntries), ...
      '  CPUseconds = ',num2str(t2)]);
  if( flagMethod1 > 0 )
  disp(['ratio of large-Du/small-Du methods: ratio = ',num2str(t2/t1)]);
  end
figure(1);
   clf;
   histogram( anglesM(:) );
   title(['large Du test: Ntries = ',num2str(Ntries)]);
   pause(0.5)
%}



