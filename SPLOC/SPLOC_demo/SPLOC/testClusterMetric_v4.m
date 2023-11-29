%%                                                                set mode
plotScatterPlot = true;
%%                                     define clustering metric parameters
haloSame = 0.01;
haloDiff0 = 0.5;
log10QdThreshold = 0.25;
log10QiThreshold = 0.25;
minQ = 0.001;
maxQ = 10;
%%                                             specify sampling properties
nBoost = 6;                                % # of replicas for each system
nR = 200;                                         % number of realizations
Ns = 20000;                                            % size of data bank
%%                                                   specify scale factors
scaleFaveX = 1.0;
scaleFstdX = 1.0;
% ------------------
scaleFaveY = 1.0;
scaleFstdY = 1.0;
% ------------------------
scaleNaveX = 1.0;
scaleNstdX = 1.0;
% ------------------
scaleNaveY = 1.0;
scaleNstdY = 1.0;


% scaleFaveX = 1.0;
% scaleFstdX = 0.01;
% % ------------------
% scaleFaveY = 10.0;
% scaleFstdY = 0.01;
% % ------------------------
% scaleNaveX = 1.0;
% scaleNstdX = 0.01;
% % ------------------
% scaleNaveY = 10.0;
% scaleNstdY = 0.01;


% scaleFaveX = 5.0;
% scaleFstdX = 0.04;
% % ------------------
% scaleFaveY = 5.0;
% scaleFstdY = 0.04;
% % ------------------------
% scaleNaveX = 1.0;
% scaleNstdX = 0.01;
% % ------------------
% scaleNaveY = 1.0;
% scaleNstdY = 0.01;
%%                                                       specify test case
% ----------------------------------------------------------------- case 1
% {
% --------------------------------------------------- 2 functional members
copulaF = 'Gaussian'; 
pho_F = [0.00, 0.00];
% --------------------------------------------------------- x-axis on MFPS
aveXF = [2.50, 2.50]; 
stdXF = [0.25, 0.25];
% --------------------------------------------------------- y-axis on MFPS
aveYF = [0.2, 0.2];
stdYF = [0.10, 0.10]; 

% ------------------------------------------------ 4 nonfunctional members
copulaN = 'Gaussian'; 
pho_N = [0.00, 0.00, 0.00, 0.00];
% --------------------------------------------------------- x-axis on MFPS
aveXN = [2.00, 2.00, 2.00, 2.00]; 
stdXN = [0.20, 0.20, 0.20, 0.20];
% --------------------------------------------------------- y-axis on MFPS
aveYN = [0.50, 0.50, 0.50, 0.50];
stdYN = [0.10, 0.10, 0.10, 0.10];
%}
% ----------------------------------------------------------------- case 2
%{
% --------------------------------------------------- 2 functional members
copulaF = 'Gaussian'; 
pho_F = [0.00, -0.00];
% --------------------------------------------------------- x-axis on MFPS
aveXF = [2.00, 2.50]; 
stdXF = [0.25, 0.20];
% --------------------------------------------------------- y-axis on MFPS
aveYF = [0.30, 0.40];
stdYF = [0.05, 0.10]; 

% ------------------------------------------------ 4 nonfunctional members
copulaN = 'Gaussian'; 
pho_N = [0.00, 0.00, 0.00, 0.00];
% --------------------------------------------------------- x-axis on MFPS
aveXN = [1.00, 2.10, 3.00, 5.00]; 
stdXN = [0.40, 0.50, 0.30, 0.45];
% --------------------------------------------------------- y-axis on MFPS
aveYN = [0.70, 0.80, 0.85, 0.75];
stdYN = [0.15, 0.20, 0.10, 0.10];
%}
% ----------------------------------------------------------------- case 3
%{
nF = 6;
nN = 6;
% ------------------------------------------------------------------------
dAlphaF = (2*pi/3)/(nF + 1);
dAlphaN = (2*pi/3)/(nN + 1);
alphaF = -(pi/3) - pi/2 + dAlphaF*(1:nF);
alphaN = (pi/3) - pi/6 + dAlphaN*(1:nN);

% --------------------------------------------------- 2 functional members
copulaF = 'Gaussian'; 
pho_F = zeros(1,nF);
% --------------------------------------------------------- x-axis on MFPS
aveXF = cos(alphaF); 
stdXF = 0.5*ones(1,nF);
% --------------------------------------------------------- y-axis on MFPS
aveYF = sin(alphaF);
stdYF = 0.5*ones(1,nF);

% ------------------------------------------------ 4 nonfunctional members
copulaN = 'Gaussian'; 
pho_N = zeros(1,nN);
% --------------------------------------------------------- x-axis on MFPS
aveXN = 5*cos(alphaN);
stdXN = 2*ones(1,nN);
% --------------------------------------------------------- y-axis on MFPS
aveYN = 5*sin(alphaN);
stdYN = 2*ones(1,nN);
%}
%%                                                          saniety checks
nF = length(aveXF);
nN = length(aveXN);
   if( length(pho_F) ~= nF )
   error('length(pho_F) is not equal to length(aveXF)');
   end
   if( length(stdXF) ~= nF )
   error('length(stdXF) is not equal to length(aveXF)');
   end
   if( length(aveYF) ~= nF )
   error('length(aveYF) is not equal to length(aveXF)');
   end
   if( length(stdYF) ~= nF )
   error('length(stdYF) is not equal to length(aveXF)');
   end
   if( length(pho_N) ~= nN )
   error('length(pho_N) is not equal to length(aveXN)');
   end
   if( length(stdXN) ~= nN )
   error('length(stdXN) is not equal to length(aveXN)');
   end
   if( length(aveYN) ~= nN )
   error('length(aveYN) is not equal to length(aveXN)');
   end
   if( length(stdYN) ~= nN )
   error('length(stdYN) is not equal to length(aveXN)');
   end
% ---------------------------------------------------------------- phase 2
   if( nBoost < 1 )
   error('nBoost must be at least 1');
   end
   if( nF < 1 )
   error('nF must be at least 1');
   end
   if( nN < 1 )
   error('nN must be at least 1');
   end
   if( nF*nBoost < 2 )
   error('nF*nBoost must be at least 2');
   end
   if( nN*nBoost < 2 )
   error('nN*nBoost must be at least 2');
   end
   if( nF*nBoost > 600 )
   error('nF*nBoost should not exceed 600');
   end
   if( nN*nBoost > 600 )
   error('nN*nBoost should not exceed 600');
   end
   if( nR < 1 )
   error('nR must be at least 1');
   end
   minScaleFactor = min([scaleFaveX,scaleFstdX, ...
                         scaleFaveY,scaleFstdY, ...
                         scaleNaveX,scaleNstdX, ...
                         scaleNaveY,scaleNstdY]);
   if( minScaleFactor < 1.0e-20 )
   error('all scale factors must be greater than 1.0E-20');
   end
   maxScaleFactor = max([scaleFaveX,scaleFstdX, ...
                         scaleFaveY,scaleFstdY, ...
                         scaleNaveX,scaleNstdX, ...
                         scaleNaveY,scaleNstdY]);
   if( maxScaleFactor > 1.0e20 )
   error('all scale factors must be less than 1.0E20');
   end
   if( Ns < 5*nR*nBoost )
   error('Ns must be at least 5*nR*nBoost');
   end
%%                                                    generate random data
xF = zeros(Ns,nF);
yF = zeros(Ns,nF);
xN = zeros(Ns,nN);
yN = zeros(Ns,nN);
str1 = ['  nF= ',num2str(nF),'  nN= ',num2str(nN)];
% ------------------------------------------------ scale numbers uniformly
aveXF = scaleFaveX*aveXF;
stdXF = scaleFstdX*stdXF;
aveYF = scaleFaveY*aveYF;
stdYF = scaleFstdY*stdYF;
% -----------------------
aveXN = scaleNaveX*aveXN;
stdXN = scaleNstdX*stdXN;
aveYN = scaleNaveY*aveYN;
stdYN = scaleNstdY*stdYN;
% -------------------------------------------- shift origin of coordinates
% ymin = min(aveYF - 5*stdYF);
% aveYF = aveYF - ymin;
% ymin = min(aveYN - 5*stdYN);
% aveYN = aveYN - ymin;
%error('stop here for now');
% ------------------------------------------------------------- functional
   for jF=1:nF
% --------------------------------------------------------------- copula F
   pho = pho_F(jF);
   u = copularnd(copulaF,[ [1.0 pho]; [pho 1.0] ],Ns);
% ------------------------------------------------------------ marginal xF
   ave = aveXF(jF);
   sig = stdXF(jF);
   xF(:,jF) = norminv(u(:,1),ave,sig);
% ------------------------------------------------------------ marginal yF
   ave = aveYF(jF);
   sig = stdYF(jF);
   yF(:,jF) = norminv(u(:,2),ave,sig);
   end
% ---------------------------------------------------------- nonfunctional
   for jN=1:nN
% --------------------------------------------------------------- copula N
   pho = pho_N(jN);
   u = copularnd(copulaN,[ [1.0 pho]; [pho 1.0] ],Ns);
% ------------------------------------------------------------ marginal xF
   ave = aveXN(jN);
   sig = stdXN(jN);
   xN(:,jN) = norminv(u(:,1),ave,sig);                 % representing mean
% ------------------------------------------------------------ marginal yF
   ave = aveYN(jN);
   sig = stdYN(jN);
   yN(:,jN) = norminv(u(:,2),ave,sig);                  % representing STD
   end 
% --------------------------------- force yF and yN to be > 0 as variances
yF = abs(yF);
yN = abs(yN);
% --------------------------------------------------------- initialization
nFig = 0;
sd = 2;
si = 1.3;
lnMidScoreX = ( log(sd) + log(si) )/2;
lnSgap = log(sd) - lnMidScoreX; 
scoreDivide = exp(lnMidScoreX);
%%                                                       plot data scatter
nFig = nFig + 1;
figure(nFig);
clf;
hold on;
   for jF=1:nF
   hF = scatter(xF(:,jF),yF(:,jF),6,'b','filled');
   end
   for jN=1:nN
   hN = scatter(xN(:,jN),yN(:,jN),6,'r','filled');
   end
legend([hF,hN],{'F','N'},'location','best');
xlabel('mean');
ylabel('STD');
title(['MFSP: Ns= ',num2str(Ns),str1]);
%error('stop here for now');
%%                                            apply unit-box normalization
maxXF = max(xF(:));
maxXN = max(xN(:));
maxX = max(maxXF,maxXN);
% -----------------------
minXF = min(xF(:));
minXN = min(xN(:));
minX = min(minXF,minXN);
% -----------------------
xSpan = max(maxX - minX, 1.0e-10);
xWidth = 1.2*xSpan;
xxF = 0.1 + (xF - minX)/xWidth;
xxN = 0.1 + (xN - minX)/xWidth;
% ------------------------------------------------------------------------
maxYF = max(yF(:));
maxYN = max(yN(:));
maxY = max(maxYF,maxYN);
% -----------------------
minYF = min(yF(:));
minYN = min(yN(:));
minY = min(minYF,minYN);
% -----------------------
ySpan = max(maxY - minY, 1.0e-10);
yWidth = 1.2*ySpan;
yyF = 0.1 + (yF - minY)/yWidth;
yyN = 0.1 + (yN - minY)/yWidth;
%error('stop here for now');
%%                                            plot normalized data scatter
nFig = nFig + 1;
figure(nFig);
clf;
hold on;
   for jF=1:nF
   hF = scatter(xxF(:,jF),yyF(:,jF),6,'b','filled');
   end
   for jN=1:nN
   hN = scatter(xxN(:,jN),yyN(:,jN),6,'r','filled');
   end
legend([hF,hN],{'F','N'},'location','best');
xlabel('mean');
ylabel('STD');
title(['normalized MFSP: Ns= ',num2str(Ns),str1]);
xlim( [0,1] );
ylim( [0,1] );
% ------------------------------------------------------------------------
% REMARK: xxF and yyF are dummy variables that refect the scaling of the 
% xF and yF data. yyF and xxF will be recalculated on the fly as it will
% be done in the actual program. 
%error('stop here for now');
%%                                 quantitative analysis of quality factor
haloDiff = 5*haloSame + haloDiff0/(sqrt(nF*nN)*nBoost);
log10Qmax = log10(maxQ);
log10Qmin = log10(minQ);
% x' = A*x + B
% 0 = A*log10Qref) + B
% log10Qmin = -2A + B
%                      => A = -c/(1 + log10Qref)    &    B = log10Qmin + A
Ad = -log10Qmin/(1 + log10QdThreshold);
Bd = log10Qmin + Ad;
Ai = -log10Qmin/(1 + log10QiThreshold);
Bi = log10Qmin + Ai;
% ------------------------------------------------------------------------
nF0 = nF;
nN0 = nN;
nF = nBoost*nF0;
nN = nBoost*nN0;
% -------------------------------------- use for on the fly dummy variable
xxF = zeros(1,nF); 
yyF = zeros(1,nF);
xxN = zeros(1,nN);
yyN = zeros(1,nN);
% --------------------------------------
score = zeros(1,nR);
rawQd = zeros(1,nR);
rawQi = zeros(1,nR);
dbF = nF0 - 1;
dbN = nN0 - 1;
log10QdMIN = 1.0e50;
log10QiMIN = 1.0e50;
log10QdMAX = -1.0e50;
log10QiMAX = -1.0e50;
   for m=1:nR    
% --------------------------------------------------------- create dataset
   indxF = randi(Ns,1,nBoost);
   indxN = randi(Ns,1,nBoost);
   bEndF = 0;
   bEndN = 0;
      for b=1:nBoost
      bStartF = bEndF + 1;
      bEndF = bStartF + dbF;
      xxF(bStartF:bEndF) = xF(indxF(b),:);
      yyF(bStartF:bEndF) = yF(indxF(b),:);
     %------------------------------------
      bStartN = bEndN + 1;
      bEndN = bStartN + dbN;
      xxN(bStartN:bEndN) = xN(indxN(b),:);
      yyN(bStartN:bEndN) = yN(indxN(b),:);
      end
%error('stop here for now');
   %{
   disp(['------------------------------------------------------ m= ', ...
         num2str(m)]);
   disp( [xxF; yyF]);
   disp( [xxN; yyN]);
   figure(50)
   clf
   hold on;
   scatter(xxF,yyF,25,'b','filled');
   scatter(xxN,yyN,25,'r','filled');
   xlabel('mean');
   ylabel('STD');
   legend('F','N');
   pause
   disp('  ');
   %}
% +++++++++++++++++++++++++++++++++++++++++++++++ calculate score function
      lnScoreFunctn10 = 0;
      n1 = nF;
      n0 = nN;
      totPairs = n1*n0;
                 for k1=1:n1                  % functional 1-state systems
                 av1 = xxF(k1);
                 vQ1 = yyF(k1)^2;
                   for k0=1:n0             % nonfunctional 0-state systems
                   av0 = xxN(k0);
                   vQ0 = yyN(k0)^2;
% -------------------------------------------------------- std. dev. ratio
                   minV = min(vQ1,vQ0);
                   maxV = max(vQ1,vQ0);
                   r = sqrt(maxV/minV);     % <= symmetrical w.r.t. labels
                   noise = sqrt(vQ1 + vQ0);
                   signal = abs(av1 - av0);
% ------------------------------------------- scoring function calculation
                   snr = signal/noise;             % signal to noise ratio
                   sbn = max(snr - 1,0);             % signal beyond noise
% --------------------------------------- calculate ln of scoring function
                   r = r - 1;    % must do for scoring function definition
                   r = r/sqrt(1 + (r/99)^2);       % apply limiting factor
                   lnScore1 = log(sqrt(sbn*sbn+r*r) + 1);  % add 1 back in
                   %                   ^^^^^^^---> harder for discriminant
                   % --> Note: sbn on its own is a viable scoring function
                   lnScore0 = log(sqrt(snr*snr+r*r) + 1);  % add 1 back in
                   %                   ^^^^^^^---> harder for indifference
                   % ---> Note: pick harder case among sbn and snr options
                     if( lnScore0 < lnMidScoreX )  % using harder criteria
                     x = lnScore0;         % => satisfies harder condition
                     elseif(lnScore1>lnMidScoreX)  % using harder criteria
                     x = lnScore1;         % => satisfies harder condition
                     else       % split case?: set lnScore as undetermined
                     x = lnMidScoreX; %set score at center of undetermined
                     end
% -------------------------------------------------------- tally all pairs
                   lnScoreFunctn10 = lnScoreFunctn10 + x;
                   end
                 end
      lnScoreFunctn10 = lnScoreFunctn10/totPairs;           % => selection
      s = exp(lnScoreFunctn10);
      score(m) = s;                          % raw score (can be modified)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------- normalize data to make scale invariant
      maxX = max([xxF,xxN]);
      minX = min([xxF,xxN]);
      xSpan = max(maxX - minX, 1.0e-10);
      xWidth = 1.2*xSpan;
      xxF = 0.1 + (xxF - minX)/xWidth;
      xxN = 0.1 + (xxN - minX)/xWidth;
% --------------------------------------
      maxY = max([yyF,yyN]);
      minY = min([yyF,yyN]);
      ySpan = max(maxY - minY, 1.0e-10);
      yWidth = 1.2*ySpan;
      yyF = 0.1 + (yyF - minY)/yWidth;
      yyN = 0.1 + (yyN - minY)/yWidth;
% ================================================ calculate quality level
      xNN = zeros(nN);
      xFF = zeros(nF);
      xNF = zeros(nN,nF);
      yNN = zeros(nN);
      yFF = zeros(nF);
      yNF = zeros(nN,nF);
      rNN = zeros(nN);
      rFF = zeros(nF);
      rNF = zeros(nN,nF);
      sNN = zeros(nN);
      sFF = zeros(nF);
      sNF = zeros(nN,nF);
         for kNa=1:nN
            for kNb=1:nN
            axx = abs( xxN(kNa) - xxN(kNb) );
            ayy = abs( yyN(kNa) - yyN(kNb) );
            arr = sqrt( axx*axx + ayy*ayy );
            xNN(kNa,kNb) = max(axx,haloSame);
            yNN(kNa,kNb) = max(ayy,haloSame);
            rNN(kNa,kNb) = max(arr,haloSame);
            sNN(kNa,kNb) = max(arr,haloDiff);
            end
         xNN(kNa,kNa) = 1.0e50;
         yNN(kNa,kNa) = 1.0e50;
         rNN(kNa,kNa) = 1.0e50;
         end
% ------------------------------------------------------
         for kFa=1:nF
            for kFb=1:nF
            axx = abs( xxF(kFa) - xxF(kFb) );
            ayy = abs( yyF(kFa) - yyF(kFb) );
            arr = sqrt( axx*axx + ayy*ayy );
            xFF(kFa,kFb) = max(axx,haloSame);
            yFF(kFa,kFb) = max(ayy,haloSame);
            rFF(kFa,kFb) = max(arr,haloSame);
            sFF(kFa,kFb) = max(arr,haloDiff);
            end
         xFF(kFa,kFa) = 1.0e50;
         yFF(kFa,kFa) = 1.0e50;
         rFF(kFa,kFa) = 1.0e50;
         end
% ------------------------------------------------------ 
         for kFa=1:nF
            for kNb=1:nN
            axx = abs( xxF(kFa) - xxN(kNb) );
            ayy = abs( yyF(kFa) - yyN(kNb) );
            arr = sqrt( axx*axx + ayy*ayy );
            xNF(kNb,kFa) = max( axx - haloDiff , 0);
            yNF(kNb,kFa) = max( ayy - haloDiff , 0);
            rNF(kNb,kFa) = max( arr - haloDiff , 0);
            sNF(kNb,kFa) = max(arr,haloDiff);
            end
         end 
% ------------------------------------------------------
      dNNmin = min(xNN);
      dFFmin = min(xFF);
      dNFmin = min(xNF);
      dFNmin = min(xNF');
      log10qNx = 2*dFNmin./(dNNmin + dFNmin) - 1;
      log10qFx = 2*dNFmin./(dFFmin + dNFmin) - 1;
%       log10qF = mean(log10qFx);
%       log10qN = mean(log10qNx);
%       %disp([log10qF,log10qN]);
%       log10Qx = min(log10qF,log10qN);
      log10Qx = sum([log10qNx,log10qFx])/(nF + nN); 
% ------------------------------------------------------
      dNNmin = min(yNN);
      dFFmin = min(yFF);
      dNFmin = min(yNF);
      dFNmin = min(yNF');
      log10qNy = 2*dFNmin./(dNNmin + dFNmin) - 1;
      log10qFy = 2*dNFmin./(dFFmin + dNFmin) - 1;
%       log10qF = mean(log10qFy);
%       log10qN = mean(log10qNy);
%       %disp([log10qF,log10qN]);
%       log10Qy = min(log10qF,log10qN);
      log10Qy = sum([log10qNy,log10qFy])/(nF + nN);  
% ------------------------------------------------------
      dNNmin = min(rNN);
      dFFmin = min(rFF);
      dNFmin = min(rNF);
      dFNmin = min(rNF');
      log10qNr = 2*dFNmin./(dNNmin + dFNmin) - 1;
      log10qFr = 2*dNFmin./(dFFmin + dNFmin) - 1;
%       log10qF = mean(log10qFr);
%       log10qN = mean(log10qNr);
%       %disp([log10qF,log10qN]);
%       log10Qr = min(log10qF,log10qN);
      log10Qr = sum([log10qNr,log10qFr])/(nF + nN);       
% ------------------------------------------------------------------------
      log10Qd = max([log10Qx,log10Qy,log10Qr]);     % => discriminant mode
      log10QdMIN = min(log10QdMIN,log10Qd);
      log10QdMAX = max(log10QdMAX,log10Qd);
      %disp(log10Qd)
      log10Qd = Ad*log10Qd + Bd;
      log10Qd = log10Qd*(1 + 10*log10Qd^2);
        if( log10Qd > 0 )
        log10Qd = log10Qd/sqrt(1 + (log10Qd/log10Qmax)^2 );            
        end
      Qd = 10^(log10Qd);                                       % raw value
% ------------------------------------------------------------------------
      dNNmin = min(sNN);
      dFFmin = min(sFF);
      dNFmin = min(sNF);
      dFNmin = min(sNF');
      log10qFs = 4*min( [dFFmin; dNFmin] )./(dFFmin + dNFmin) - 1;
      log10qNs = 4*min( [dNNmin; dFNmin] )./(dNNmin + dFNmin) - 1;
      log10qF = mean(log10qFs);
      log10qN = mean(log10qNs);
      log10Qi = min(log10qF,log10qN);
      log10QiMIN = min(log10QiMIN,log10Qi);
      log10QiMAX = max(log10QiMAX,log10Qi);
      log10Qi = Ai*log10Qi + Bi;
      log10Qi = log10Qi*(1 + 10*log10Qi^2);
        if( log10Qi > 0 )
        log10Qi = log10Qi/sqrt(1 + (log10Qi/log10Qmax)^2 );
        end
      Qi = 10^(log10Qi);                                       % raw value
% ---------------------------------------------------- visualization check
         if( plotScatterPlot && (m < 16) )
            if( s < scoreDivide )                  % => possibly an i-mode
               if( (Qi - Qd) > 1 )                          % => an i-mode
               modeType = 'i-mode';
               Q = Qi;
               else                                            % => u-mode
               modeType = 'u-mode';
               dQ = (Qi - Qd);
                  if( dQ < -1 )
                  Q = -Qd;                            % => maximum penalty
                  else
                  fi = 0.5*(dQ + 1);
                  fd = 1 - fi;
                  Q = fi*Qi - fd*Qd;             % include gradual penalty
                  end
               end
            else                                    % => possibly a d-mode
               if( (Qd - Qi) > 1 )                           % => a d-mode
               modeType = 'd-mode';
               Q = Qd;
               else                                            % => u-mode
               modeType = 'u-mode';
               dQ = (Qd - Qi);
                  if( dQ < -1 )
                  Q = -Qi;                            % => maximum penalty
                  else
                  fd = 0.5*(dQ + 1);
                  fi = 1 - fd;
                  Q = fd*Qd - fi*Qi;             % include gradual penalty
                  end
               end
            end
         figure(50);
         clf;
         hold on;
         scatter(xxF,yyF,'b','filled');
         scatter(xxN,yyN,'r','filled');
         legend('F','N','location','best');
         xlabel('mean');
         ylabel('STD');
         title(['normalized MFSP: (nF,nN)= (',num2str(nF),',', ...
                 num2str(nN),')   ',modeType,'   Qd= ', ...
                 num2str(Qd),'   Qi= ',num2str(Qi), ...
                 '  score= ',num2str(s),'   Q= ',num2str(Q)]);
         pause
         end
% -------------------------------------------------------- data collection
      rawQd(m) = Qd;
      rawQi(m) = Qi;
   end
nFig = nFig + 1;
%%                                                           data analysis
disp(['log10QdMIN = ',num2str(log10QdMIN)]);
disp(['log10QiMIN = ',num2str(log10QiMIN)]);
disp(['log10QdMAX = ',num2str(log10QdMAX)]);
disp(['log10QiMAX = ',num2str(log10QiMAX)]);
figure(nFig);
clf;
hold on;
scatter( score(:), rawQd(:), 30,'r','filled');
scatter( score(:), rawQi(:), 30,'b','filled');
ydMax = max( rawQd(:) );
yiMax = max( rawQi(:) );
yMax = max( ydMax, yiMax );
x = 2*ones(1,11);
y = 1.05*yMax*(0:0.1:1);
plot(x,y,'c','linewidth',1.8);
x = 1.3*ones(1,11);
plot(x,y,'m','linewidth',1.8);
x = scoreDivide*ones(1,11);
plot(x,y,'k--','linewidth',2.2);
ylim( [0,1.05*yMax] );
xlabel('selection score');
ylabel('raw cluster quality');
title(['nR= ',num2str(nR),'   nBoost= ',num2str(nBoost),str1]);
legend('Qd','Qi','location','best');
% ---------------------------------------------------------- final Q vaues
nFig = nFig + 1;
%%
figure(nFig);
clf;
hold on;
finalQ = zeros(1,nR);
uMode = false(1,nR);
   for m=1:nR   
% % % %    s = score(m);
% % % %    lnS = log(s);
% % % %    sFactor = abs(lnS - lnMidScoreX)/lnSgap;
% % % %    sFactor = min(1,sFactor);
   Qi = rawQi(m);
   Qd = rawQd(m);
      if( score(m) < scoreDivide )                 % => possibly an i-mode
         if( (Qi - Qd) > 1 )                                % => an i-mode
         %modeType = 'i-mode';
         Q = Qi;
         else                                                  % => u-mode
         %modeType = 'u-mode';
         uMode(m) = true;
         dQ = (Qi - Qd);
            if( dQ < -1 )
            Q = -Qd;                                  % => maximum penalty
            else
            fi = 0.5*(dQ + 1);
            fd = 1 - fi;
            Q = fi*Qi - fd*Qd;                   % include gradual penalty
            end
         end
      else                                          % => possibly a d-mode
         if( (Qd - Qi) > 1 )                                 % => a d-mode
         %modeType = 'd-mode';
         Q = Qd;
         else                                                  % => u-mode
         %modeType = 'u-mode';
         uMode(m) = true;
         dQ = (Qd - Qi);
            if( dQ < -1 )
            Q = -Qi;                                  % => maximum penalty
            else
            fd = 0.5*(dQ + 1);
            fi = 1 - fd;
            Q = fd*Qd - fi*Qi;                   % include gradual penalty
            end
         end
      end
   finalQ(m) = Q;
   end
uModeNot = ~uMode;
scatter( score(uModeNot), finalQ(uModeNot), 30,'k','filled');
scatter( score(uMode), finalQ(uMode), 30,'g','filled');
yMax = max(1.05*max(finalQ),10.001);
yMin = 1.05*min(min(finalQ),0);
x = 2*ones(1,11);
y = yMin + (yMax - yMin)*(0:0.1:1);
plot(x,y,'c','linewidth',1.8);
x = 1.3*ones(1,11);
plot(x,y,'m','linewidth',1.8);
x = scoreDivide*ones(1,11);
plot(x,y,'k--','linewidth',2.2);
ylim( [yMin,yMax]);
xlabel('selection score');
ylabel('cluster quality');
title(['nR= ',num2str(nR),'   nBoost= ',num2str(nBoost),str1]);
legend('Nay u-mode','Yay u-mode','location','best');

% boneyard
% ------------------------------------------------------------------------
% % % % %%
% % % % nFig = nFig + 1;
% % % % figure(nFig);
% % % % clf;
% % % % hold on;
% % % %   for m=1:Ms
% % % %     for jF=1:MF
% % % %     nF = nFarray(jF);
% % % %       for jN=1:MN
% % % %       nN = nNarray(jN);
% % % %       s = score(jN,jF,m);
% % % %          if( s < scoreDivide )
% % % %          scatter(nF,nN,25,'g');
% % % %          end
% % % %       end
% % % %     end
% % % %   end
% % % % xlabel('nF');
% % % % ylabel('nN');

% % aveQd = zeros(MN,MF);
% % stdQd = zeros(MN,MF);
% % aveQi = zeros(MN,MF);
% % stdQi = zeros(MN,MF);
  
%{  
aveQd = aveQd/Ms;
av2Quality = stdQd/Ms;
stdQd = sqrt( 1.0e-12 + av2Quality - aveQd.^2 );
   for jF=1:MF
   %errorbar(nNarray,aveQuality(jF,:),stdQuality(jF,:),'b-');
   plot(nNarray,aveQd(jF,:),'b-');
   end
xlabel('number of nonfunctional samples');
ylabel('quality level');
title(['MFSP: quality contours: nF= ', ...
       num2str(nFarray(1)),' to ',num2str(nFarray(end))]);
%error('stop here for now');

nFig = nFig + 1;
figure(nFig);
clf;
hold on;
% --------------------------------------------- slice through different nF
%}


%{
% ------------------------------------------- final stretch to calculate Q
         if( Qi > Qd )
         else
         end
        
      Q = max( min(sFactor*10^(1.05*log10Q),10),0.1);
%}


%{
% --------------------------------------------- calculate scoring function
         if( s > scoreDivide )
         xs = (s - scoreDivide)/(2 - scoreDivide);
         modeType = 'd-mode';
         dMode = 1;
         else
         xs = (scoreDivide - s)/(scoreDivide - 1.3);
         modeType = 'i-mode';
         dMode = -1;
         end
      sFactor = ( min(1,xs) )^2;
 %}
      
     