function [msipMap] = subspaceComparison(cSBV1,cCIP1, ...
                                        cSBV2,cCIP2,verbosity,fname)
% calculates MSIP comparisons & plots results on a 2D grid and histograms 
% 
% INPUT (two modes)
% ----------------------------------------------------------------- mode 1
% subspaceComparison(cSBV1,cCIP1,verbosity,fname)
% subspaceComparison(cSBV1,cCIP1,verbosity)
% subspaceComparison(cSBV1,cCIP1)
%
% cSBV1 = cell array of selection basis vectors
% cCIP1 = cell array of congruency indicator partitions
%         1 => congruent discriminant subspace
%         0 => congruent undetermined subspace
%        -1 => congruent indifference subspace
% verbosity: 0 => output msipMap only
%            1 => output msipMap and print figures to the screen.
%            2 => output msipMap and print figures to file only.
%            3 => output msipMap and print figures to file and screen.
% fname = base file name to give to output files
%
% ----------------------------------------------------------------- mode 2
% subspaceComparison(cSBV1,cCIP1,cSBV2,cCIP2,verbosity,fname)
% subspaceComparison(cSBV1,cCIP1,cSBV2,cCIP2,verbosity)
% subspaceComparison(cSBV1,cCIP1,cSBV2,cCIP2)
%
% Input variables are the same as mode 1 plus:
% cSBV2 = cell array of selection basis vectors for second block
% cCIP2 = cell array of congruency indicator partitions for second block
%
% PROCESS
% Each selection basis vector set has three congruent types of subspaces:
% d => discriminant, u => undetermined, i => indifference
% MSIP pairwise calculations compare similar and dissimilar subspaces.
% similar    --> pairwise subspace comparisons are dd, uu and ii.
% dissimilar --> pairwise subspace comparisons are du, di and ui. 
%
% Based on mode, the entire msipMap is calculated, which gives the MSIP
% values of each system stratisfied across the d,u,i congruent subspaces. 
%
% Plots a heat map of all pairwise MSIP comparisons and histograms of the
% MSIP pairs. There are two modes of operation that can be considered. A
% comparison over one block of data defined by cSBV1,cCIP1 or two blocks
% cSBV1,cCIP1 and cSBV2,cCIP2. For one block of data, one histrogram is
% plotted for self-comparisons and if two blocks are specified, a self-
% comparison is made per block and a cross-comparison between the two 
% blocks of data is made. The size of msipMap depends on the mode. 
%
% OUTPUT  (no plot is made when verbosity = 0)
% msipMap = the summary of all calculations that are performed
%
% PLOTS (graphs to screen when verbosity = 1,3, and to file for 2,3)
% mode 1: 
% histogram of MSIP values for block 1 data (intra-comparison).
% a heat map of all intra-pairwise comparisons
%
% mode 2:
% histogram of MSIP values per block (two intra-comparisons).
% histogram of MSIP values for cross comparisons.
% one heat map that gives all intra- and inter-pairwise comparisons
%
% Note: The format for a histogram shows MSIP values for expected similar
% subspaces that are clustered as (dd+uu+ii) and for expected orthogonal
% subspaces that are clustered as (du+di+ui). Information about the origin 
% of MSIP values from either similar or dissimilar subspaces is lost. The
% heat map showing all individual MSIP pairwise comparisions will allow a
% user to extract this information.
%
% heat maps are arranged as:
%          data1-d  data1-u  data1-i  data2-d  data2-u  data2-i
% data1-d
% data1-u
% data1-i
% data2-d
% data2-u
% data2-i
%%                                             set global sploc parameters
global gvSPLOC                   % shares information across sploc toolset
%%                                check input and determine operation mode
   if( nargin == 2 )
   verbosity = 0;                                          % default value
   fname = 'unnamed';                                      % default value
   flagMode = 1;
   elseif( nargin == 3 )
      if( isnumeric(cSBV2) )                  % => this value is verbosity
      verbosity = cSBV2;
      fname = 'unnamed';                                   % default value
      flagMode = 1;
      else
      flagMode = -1;                               % => input format error
      end
   elseif( nargin == 4 )          % =====> interface between modes 1 and 2
      if( isnumeric(cSBV2) )                  % => this value is verbosity
      verbosity = cSBV2;
      fname = cCIP2;                              % => this value is fname
      flagMode = 1;
      else                     
         if( iscell(cSBV2) )                                   % => mode 2
         verbosity = 0;                                    % default value
         fname = 'unnamed';                                % default value
         flagMode = 2; 
         else
         flagMode = -1;                            % => input format error
         end
      end
   elseif( nargin == 5 )                               % => must be mode 2
     if( iscell(cSBV2) )             % => this is a requirement for mode 2
     fname = 'unnamed';                                    % default value
     flagMode = 2;
     else
     flagMode = -1;                                % => input format error
     end
   elseif( nargin == 6 )  
   flagMode = 2;
   else
   flagMode = -1;                                 % cannot have nargin < 2
   end
% -------------------------------------- terminate if an error is detected
  if( flagMode < 0 )
  disp('   ');
  disp('allowed usage:');
  disp('subspaceComparison(cSBV1,cCIP1)');
  disp('subspaceComparison(cSBV1,cCIP1,verbosity)');
  disp('subspaceComparison(cSBV1,cCIP1,verbosity,fname)');
  disp('subspaceComparison(cSBV1,cCIP1,cSBV2,cCIP2)');
  disp('subspaceComparison(cSBV1,cCIP1,cSBV2,cCIP2,verbosity)');
  disp('subspaceComparison(cSBV1,cCIP1,cSBV2,cCIP2,verbosity,fname)');
  error('Input format error: Incorrect usage detected');
  end
verbosity = setVerbosity(verbosity);
% disp(['     file name = ',fname])                  % for debugging input
% disp(['     verbosity = ',num2str(verbosity)])     % for debugging input
% disp(['operation mode = ',num2str(flagMode)])      % for debugging input
% disp('-----------------------------------------'); % for debugging input
% disp('   ');                                       % for debugging input
%%                                         partition subspaces for block 1
Z = 0*cSBV1{1}(:,1);      % DIM of vector space is same for blocks 1 and 2
% Z= zero vector ^------> any vector from any block is just as good to use
% ---------------------------------------------------------------- block 1
N1 = length(cSBV1);
Ld1 = cell(1,N1);
Lu1 = cell(1,N1);
Li1 = cell(1,N1);
dd1 = zeros(1,N1);
du1 = zeros(1,N1);
di1 = zeros(1,N1);
   for k=1:N1
   Ld1{k} = ( cCIP1{k} == 1 );
   Lu1{k} = ( cCIP1{k} == 0 );
   Li1{k} = ( cCIP1{k} == -1);
   dd1(k) = sum( Ld1{k} );
   du1(k) = sum( Lu1{k} );
   di1(k) = sum( Li1{k} );
   end
%%                           mean squared inner products between subspaces
msip_dd1 = zeros(N1);
msip_du1 = zeros(N1);
msip_di1 = zeros(N1);
msip_uu1 = zeros(N1);
msip_ui1 = zeros(N1);
msip_ii1 = zeros(N1);
   for j=1:N1
% ----------------------------
      if( dd1(j) > 0 )
      Ud = cSBV1{j}(:,Ld1{j});
      else
      Ud = Z;
      end
% ----------------------------
      if( du1(j) > 0 )
      Uu = cSBV1{j}(:,Lu1{j});
      else
      Uu = Z;
      end
% ----------------------------
      if( di1(j) > 0 )
      Ui = cSBV1{j}(:,Li1{j});
      else
      Ui = Z;
      end
% ----------------------------------------------- construct the V matrices
      for k=j:N1                  % k=j => must calculate the diagonal too
% -------------------------------
         if( dd1(k) > 0 )
         Vd = cSBV1{k}(:,Ld1{k});
         else
         Vd = Z;
         end
% -------------------------------
         if( du1(k) > 0 )
         Vu = cSBV1{k}(:,Lu1{k});
         else
         Vu = Z;
         end
% -------------------------------
         if( di1(k) > 0 )
         Vi = cSBV1{k}(:,Li1{k});
         else
         Vi = Z;
         end
% ---------------------------------------------- calculate matrix elements
      msip_dd1(k,j) = mean( sum( (Ud'*Vd).^2 ,2) );
      msip_du1(k,j) = mean( sum( (Ud'*Vu).^2 ,2) );
      msip_di1(k,j) = mean( sum( (Ud'*Vi).^2 ,2) );
      msip_uu1(k,j) = mean( sum( (Uu'*Vu).^2 ,2) );
      msip_ui1(k,j) = mean( sum( (Uu'*Vi).^2 ,2) );
      msip_ii1(k,j) = mean( sum( (Ui'*Vi).^2 ,2) );
      % -------------------------------------------- get transposed values
      msip_dd1(j,k) = msip_dd1(k,j);
      msip_du1(j,k) = msip_du1(k,j);
      msip_di1(j,k) = msip_di1(k,j);
      msip_uu1(j,k) = msip_uu1(k,j);
      msip_ui1(j,k) = msip_ui1(k,j);
      msip_ii1(j,k) = msip_ii1(k,j);
      end
   end
a1 = [msip_dd1,msip_du1,msip_di1];
b1 = [msip_du1,msip_uu1,msip_ui1];
c1 = [msip_di1,msip_ui1,msip_ii1];
%%                                            process data based on mode 1
  if( flagMode == 1 )                             % work only with block 1
  msipMap = [a1; b1; c1];
     if( verbosity > 0 )
        if( 2*floor(verbosity/2) ~= verbosity )               % for screen
        showFig = 'yes';
        figure;
        figure_number = get(gcf,'Number');
        h_msipS11 = figure(figure_number);                       % => Same
        else                               % => do not show plot on screen
        showFig = 'no';
        h_msipS11 = figure('visible','off');
        end
% -------------------------------------------------------------
     Ldd1 = ( triu(ones( size(msip_dd1) ),1) == 1 );
     ydd1 = msip_dd1(Ldd1);
     Luu1 = ( triu(ones( size(msip_uu1) ),1) == 1 );
     yuu1 = msip_uu1(Luu1);
     Lii1 = ( triu(ones( size(msip_ii1) ),1) == 1 );
     yii1 = msip_ii1(Lii1);
% -------------------------------------------------------------
     ySame1 = [ydd1', yuu1', yii1'];
     edg = 0:0.05:1;
     histogram(ySame1,edg);
     xlim([0,1]);
     xlabel('MSIP of dd, uu, ii subspaces');
     ylabel('distribution of counts');
     str1 = [num2str(length(ySame1)),' # of squared inner product pairs'];
     title(str1,'Interpreter','none');
% ------------------------------------------------------------------------
        if( 2*floor(verbosity/2) ~= verbosity )               % for screen
        figure_number = figure_number + 1;
        h_msipD11 = figure(figure_number);                  % => Different
        else                               % => do not show plot on screen
        h_msipD11 = figure('visible','off');
        end
     Ldu1 = ( triu(ones( size(msip_du1) ),1) == 1 );
     ydu1 = msip_du1(Ldu1);
     Ldi1 = ( triu(ones( size(msip_di1) ),1) == 1 );
     ydi1 = msip_di1(Ldi1);
     Lui1 = ( triu(ones( size(msip_ui1) ),1) == 1 );
     yui1 = msip_ui1(Lui1);
     yDiff1 = [ydu1', ydi1', yui1'];
     histogram(yDiff1,edg);
     xlim([0,1]);
     xlabel('MSIP of du, di, ui subspaces');
     ylabel('distribution of counts');
     str2 = [num2str(length(yDiff1)),' # of squared inner product pairs'];
     title(str2,'Interpreter','none');
% ------------------------------------------------------------------------
     tLstr = [fname,'  Intra-comparisons across subsaces 1x1'];
        if( verbosity > 1 )                                     % for file
        subFolder = 'basisComparison';
        fNameS11 = [fname,'_MSIP11Shist'];      % => 1-1 intra-comparisons
        fNameS11 = getOutputFileName(subFolder,fNameS11);   
        fNameD11 = [fname,'_MSIP11Dhist'];      % => 1-1 intra-comparisons
        fNameD11 = getOutputFileName(subFolder,fNameD11);
        fNameMap = [fname,'_MSIPmatrix'];
        fNameMap = getOutputFileName(subFolder,fNameMap);
        disp(fNameS11)
        disp(fNameD11)
        disp(fNameMap)
        colorMatrixTool(msipMap,0,'clean','no', ...
                   'xLabel','(vector indices): (d+u+i)', ...
                   'yLabel','(vector indices): (d+u+i)', ...
                   'scale','value','cutlower',0,'cutupper',1, ...
                   'fShow',showFig,'fName',fNameMap, ...
                   'fType',gvSPLOC.gFileType, ...
                   'tLabel',tLstr,'cType','0b','bLabel','MSIP');
        saveas(h_msipS11,fNameS11,gvSPLOC.gFileType);
        saveas(h_msipD11,fNameD11,gvSPLOC.gFileType);
        else                                                % not for file
        colorMatrixTool(msipMap,0,'clean','no', ...
                   'xLabel','(vector indices): (d+u+i)', ...
                   'yLabel','(vector indices): (d+u+i)', ...
                   'scale','value','cutlower',0,'cutupper',1, ...
                   'fShow',showFig, ...
                   'tLabel',tLstr,'cType','0b','bLabel','MSIP');
        end
     end
  return
  end
%%                                         partition subspaces for block 2
N2 = length(cSBV2);
Ld2 = cell(1,N2);
Lu2 = cell(1,N2);
Li2 = cell(1,N2);
dd2 = zeros(1,N2);
du2 = zeros(1,N2);
di2 = zeros(1,N2);
   for k=1:N2
   Ld2{k} = ( cCIP2{k} == 1 );
   Lu2{k} = ( cCIP2{k} == 0 );
   Li2{k} = ( cCIP2{k} == -1);
   dd2(k) = sum( Ld2{k} );
   du2(k) = sum( Lu2{k} );
   di2(k) = sum( Li2{k} );
   end
%%                           mean squared inner products between subspaces
msip_dd2 = zeros(N2);
msip_du2 = zeros(N2);
msip_di2 = zeros(N2);
msip_uu2 = zeros(N2);
msip_ui2 = zeros(N2);
msip_ii2 = zeros(N2);
   for j=1:N2
% ----------------------------
      if( dd2(j) > 0 )
      Ud = cSBV2{j}(:,Ld2{j});
      else
      Ud = Z;
      end
% ----------------------------
      if( du2(j) > 0 )
      Uu = cSBV2{j}(:,Lu2{j});
      else
      Uu = Z;
      end
% ----------------------------
      if( di2(j) > 0 )
      Ui = cSBV2{j}(:,Li2{j});
      else
      Ui = Z;
      end
% ----------------------------------------------- construct the V matrices
      for k=j:N2                  % k=j => must calculate the diagonal too
% -------------------------------
         if( dd2(k) > 0 )
         Vd = cSBV2{k}(:,Ld2{k});
         else
         Vd = Z;
         end
% -------------------------------
         if( du2(k) > 0 )
         Vu = cSBV2{k}(:,Lu2{k});
         else
         Vu = Z;
         end
% -------------------------------
         if( di2(k) > 0 )
         Vi = cSBV2{k}(:,Li2{k});
         else
         Vi = Z;
         end
% ---------------------------------------------- calculate matrix elements
      msip_dd2(k,j) = mean( sum( (Ud'*Vd).^2 ,2) );
      msip_du2(k,j) = mean( sum( (Ud'*Vu).^2 ,2) );
      msip_di2(k,j) = mean( sum( (Ud'*Vi).^2 ,2) );
      msip_uu2(k,j) = mean( sum( (Uu'*Vu).^2 ,2) );
      msip_ui2(k,j) = mean( sum( (Uu'*Vi).^2 ,2) );
      msip_ii2(k,j) = mean( sum( (Ui'*Vi).^2 ,2) );
      % -------------------------------------------- get transposed values
      msip_dd2(j,k) = msip_dd2(k,j);
      msip_du2(j,k) = msip_du2(k,j);
      msip_di2(j,k) = msip_di2(k,j);
      msip_uu2(j,k) = msip_uu2(k,j);
      msip_ui2(j,k) = msip_ui2(k,j);
      msip_ii2(j,k) = msip_ii2(k,j);
      end
   end
%%                 mean squared inner products between different subspaces
msip_dd12 = zeros(N1,N2);
msip_du12 = zeros(N1,N2);
msip_di12 = zeros(N1,N2);
msip_uu12 = zeros(N1,N2);
msip_ui12 = zeros(N1,N2);
msip_ii12 = zeros(N1,N2);
   for j=1:N1
% ----------------------------
      if( dd1(j) > 0 )
      Ud = cSBV1{j}(:,Ld1{j});
      else
      Ud = Z;
      end
% ----------------------------
      if( du1(j) > 0 )
      Uu = cSBV1{j}(:,Lu1{j});
      else
      Uu = Z;
      end
% ----------------------------
      if( di1(j) > 0 )
      Ui = cSBV1{j}(:,Li1{j});
      else
      Ui = Z;
      end
% ----------------------------------------------- construct the V matrices
      for k=1:N2
% -------------------------------
         if( dd2(k) > 0 )
         Vd = cSBV2{k}(:,Ld2{k});
         else
         Vd = Z;
         end
% -------------------------------
         if( du2(k) > 0 )
         Vu = cSBV2{k}(:,Lu2{k});
         else
         Vu = Z;
         end
% -------------------------------
         if( di2(k) > 0 )
         Vi = cSBV2{k}(:,Li2{k});
         else
         Vi = Z;
         end
% ---------------------------------------------- calculate matrix elements
      msip_dd12(j,k) = mean( sum( (Ud'*Vd).^2 ,2) );
      msip_du12(j,k) = mean( sum( (Ud'*Vu).^2 ,2) );
      msip_di12(j,k) = mean( sum( (Ud'*Vi).^2 ,2) );
      msip_uu12(j,k) = mean( sum( (Uu'*Vu).^2 ,2) );
      msip_ui12(j,k) = mean( sum( (Uu'*Vi).^2 ,2) );
      msip_ii12(j,k) = mean( sum( (Ui'*Vi).^2 ,2) );
      end
   end
% % ----------------------------------------- assemble the parts together
% % REMARK: Packing format: (d1 + u1 + i1)-(d2 + u2 + i2)
% % a1,b1,c1 are of size:                                          N1 x N1
% a2 = [msip_dd2,msip_du2,msip_di2];                             % N2 x N2
% b2 = [msip_du2,msip_uu2,msip_ui2];
% c2 = [msip_di2,msip_ui2,msip_ii2];
% a12 = [msip_dd12,msip_du12,msip_di12];                         % N1 x N2
% b12 = [msip_du12,msip_uu12,msip_ui12];
% c12 = [msip_di12,msip_ui12,msip_ii12];
% M11 = [a1;  b1;  c1];
% M12 = [a12; b12; c12];
% M21 = M12';
% M22 = [a2;  b2;  c2];
% m1 = [M11, M12];
% m2 = [M21, M22];
% msipMap = [m1; m2];
% -------------------------------------------- assemble the parts together
% packing format: (d1+d2)-(u1+u2)-(i1+i2) with block size= (N1+N2)x(N1+N2)
msip_dd21 = msip_dd12';
block_dd = [ [msip_dd1, msip_dd12]; [msip_dd21, msip_dd2] ];
%
msip_du21 = msip_du12';
block_du = [ [msip_du1, msip_du12]; [msip_du21, msip_du2] ];  
%
msip_di21 = msip_di12';
block_di = [ [msip_di1, msip_di12]; [msip_di21, msip_di2] ];
% --------------------------------------------------------------
msip_uu21 = msip_uu12';
block_uu = [ [msip_uu1, msip_uu12]; [msip_uu21, msip_uu2] ];  
msip_ui21 = msip_ui12';
block_ui = [ [msip_ui1, msip_ui12]; [msip_ui21, msip_ui2] ];
% --------------------------------------------------------------
msip_ii21 = msip_ii12';
block_ii = [ [msip_ii1, msip_ii12]; [msip_ii21, msip_ii2] ]; 
% size(block_dd)
% size(block_du)
% size(block_di)
% size(block_uu)
% size(block_ui)
% size(block_ii)
% ========================================================================
block_ud = block_du';
block_id = block_di';
block_iu = block_ui';
msipMap = [ [block_dd, block_du, block_di]; ...
            [block_ud, block_uu, block_ui]; ...
            [block_id, block_iu, block_ii] ];
     if( verbosity > 0 )
% ================================================================ same 11
        if( 2*floor(verbosity/2) ~= verbosity )               % for screen
        showFig = 'yes';
        figure;
        figure_number = get(gcf,'Number');
        h_msipS11 = figure(figure_number);                       % => Same
        else                               % => do not show plot on screen
        showFig = 'no';
        h_msipS11 = figure('visible','off');
        end
% -------------------------------------------------------------
     Ldd1 = ( triu(ones( size(msip_dd1) ),1) == 1 );
     ydd1 = msip_dd1(Ldd1);
     Luu1 = ( triu(ones( size(msip_uu1) ),1) == 1 );
     yuu1 = msip_uu1(Luu1);
     Lii1 = ( triu(ones( size(msip_ii1) ),1) == 1 );
     yii1 = msip_ii1(Lii1);
% -------------------------------------------------------------
     ySame1 = [ydd1', yuu1', yii1'];
     edg = 0:0.05:1;
     histogram(ySame1,edg);
     xlim([0,1]);
     xlabel('MSIP of dd, uu, ii subspaces 11');
     ylabel('distribution of counts');
     str3 = [num2str(length(ySame1)),' # of squared inner product pairs'];
     title(str3,'Interpreter','none');
% ================================================================ diff 11
        if( 2*floor(verbosity/2) ~= verbosity )               % for screen
        figure_number = figure_number + 1;
        h_msipD11 = figure(figure_number);                  % => Different
        else                               % => do not show plot on screen
        h_msipD11 = figure('visible','off');
        end
     Ldu1 = ( triu(ones( size(msip_du1) ),1) == 1 );
     ydu1 = msip_du1(Ldu1);
     Ldi1 = ( triu(ones( size(msip_di1) ),1) == 1 );
     ydi1 = msip_di1(Ldi1);
     Lui1 = ( triu(ones( size(msip_ui1) ),1) == 1 );
     yui1 = msip_ui1(Lui1);
     yDiff1 = [ydu1', ydi1', yui1'];
     mmm = length(yDiff1);
     histogram(yDiff1,edg);
     xlim([0,1]);
     xlabel('MSIP of du, di, ui subspaces 11');
     ylabel('distribution of counts');
     str4 = [num2str(mmm),' # of squared inner product pairs'];
     title(str4,'Interpreter','none');
% ================================================================ same 22
        if( 2*floor(verbosity/2) ~= verbosity )               % for screen
        figure_number = figure_number + 1;
        h_msipS22 = figure(figure_number);                       % => Same
        else                               % => do not show plot on screen
        h_msipS22 = figure('visible','off');
        end
% -------------------------------------------------------------
     Ldd2 = ( triu(ones( size(msip_dd2) ),1) == 1 );
     ydd2 = msip_dd2(Ldd2);
     Luu2 = ( triu(ones( size(msip_uu2) ),1) == 1 );
     yuu2 = msip_uu2(Luu2);
     Lii2 = ( triu(ones( size(msip_ii2) ),1) == 1 );
     yii2 = msip_ii2(Lii2);
% -------------------------------------------------------------
     ySame2 = [ydd2', yuu2', yii2'];
     mmm = length(ySame2);
     histogram(ySame2,edg);
     xlim([0,1]);
     xlabel('MSIP of dd, uu, ii subspaces 22');
     ylabel('distribution of counts');
     str5 = [num2str(mmm),' # of squared inner product pairs'];
     title(str5,'Interpreter','none');
% ================================================================ diff 22
        if( 2*floor(verbosity/2) ~= verbosity )               % for screen
        figure_number = figure_number + 1;
        h_msipD22 = figure(figure_number);                  % => Different
        else                               % => do not show plot on screen
        h_msipD22 = figure('visible','off');
        end
     Ldu2 = ( triu(ones( size(msip_du2) ),1) == 1 );
     ydu2 = msip_du2(Ldu2);
     Ldi2 = ( triu(ones( size(msip_di2) ),1) == 1 );
     ydi2 = msip_di2(Ldi2);
     Lui2 = ( triu(ones( size(msip_ui2) ),1) == 1 );
     yui2 = msip_ui2(Lui2);
     yDiff2 = [ydu2', ydi2', yui2'];
     mmm = length(yDiff2); 
     histogram(yDiff2,edg);
     xlim([0,1]);
     xlabel('MSIP of du, di, ui subspaces 22');
     ylabel('distribution of counts');
     str6 = [num2str(mmm),' # of squared inner product pairs'];
     title(str6,'Interpreter','none');
% ================================================================ same 12
        if( 2*floor(verbosity/2) ~= verbosity )               % for screen
        figure_number = figure_number + 1;
        h_msipS12 = figure(figure_number);                       % => Same
        else                               % => do not show plot on screen
        h_msipS12 = figure('visible','off');
        end
% -------------------------------------------------------------
     ydd12 = msip_dd12(1:end);
     yuu12 = msip_uu12(1:end);
     yii12 = msip_ii12(1:end);
% -------------------------------------------------------------
     ySame12 = [ydd12, yuu12, yii12];
     mmm = length(ySame12);
     histogram(ySame12,edg);
     xlim([0,1]);
     xlabel('MSIP of dd, uu, ii subspaces 12');
     ylabel('distribution of counts');
     str7 = [num2str(mmm),' # of squared inner product pairs'];
     title(str7,'Interpreter','none');
% ================================================================ diff 12
        if( 2*floor(verbosity/2) ~= verbosity )               % for screen
        figure_number = figure_number + 1;
        h_msipD12 = figure(figure_number);                  % => Different
        else                               % => do not show plot on screen
        h_msipD12 = figure('visible','off');
        end
     ydu12 = msip_du12(1:end);
     ydi12 = msip_di12(1:end);
     yui12 = msip_ui12(1:end);
     yDiff12 = [ydu12, ydi12, yui12];
     mmm = length(yDiff12);
     histogram(yDiff12,edg);
     xlim([0,1]);
     xlabel('MSIP of du, di, ui subspaces 12');
     ylabel('distribution of counts');
     str8 = [num2str(mmm),' # of squared inner product pairs'];
     title(str8,'Interpreter','none');
% =============================================================== MSIP map
     tLstr = [fname,'  comparisons across subsaces (1+2)x(1+2)'];
        if( verbosity > 1 )                                     % for file
        subFolder = 'basisComparison';
        fNameS11 = [fname,'_MSIP11Shist'];             % => 1-1 comparison
        fNameS11 = getOutputFileName(subFolder,fNameS11);   
        fNameD11 = [fname,'_MSIP11Dhist'];             % => 1-1 comparison
        fNameD11 = getOutputFileName(subFolder,fNameD11);
        % ----------------
        fNameS22 = [fname,'_MSIP22Shist'];             % => 2-2 comparison
        fNameS22 = getOutputFileName(subFolder,fNameS22);   
        fNameD22 = [fname,'_MSIP22Dhist'];             % => 2-2 comparison
        fNameD22 = getOutputFileName(subFolder,fNameD22);
        % ----------------
        fNameS12 = [fname,'_MSIP12Shist'];             % => 1-2 comparison
        fNameS12 = getOutputFileName(subFolder,fNameS12);   
        fNameD12 = [fname,'_MSIP12Dhist'];             % => 1-2 comparison
        fNameD12 = getOutputFileName(subFolder,fNameD12);
        % ----------------
        fNameMap = [fname,'_MSIPmatrix'];
        fNameMap = getOutputFileName(subFolder,fNameMap);
        colorMatrixTool(msipMap,0,'fType',gvSPLOC.gFileType, ...
                'xLabel','(vector indices): (d1+d2)-(u1+u2)-(i1+i2)', ...
                'yLabel','(vector indices): (d1+d2)-(u1+u2)-(i1+i2)', ...
                'scale','value','cutlower',0,'cutupper',1, ...
                'clean','no','fShow',showFig,'fName',fNameMap, ... 
                'tLabel',tLstr,'cType','0b','bLabel','MSIP');
        saveas(h_msipS11,fNameS11,gvSPLOC.gFileType);
        saveas(h_msipD11,fNameD11,gvSPLOC.gFileType);
        saveas(h_msipS22,fNameS22,gvSPLOC.gFileType);
        saveas(h_msipD22,fNameD22,gvSPLOC.gFileType);
        saveas(h_msipS12,fNameS12,gvSPLOC.gFileType);
        saveas(h_msipD12,fNameD12,gvSPLOC.gFileType);
        else                                                % not for file
        colorMatrixTool(msipMap,0,'clean','no','fShow',showFig, ...
                'xLabel','(vector indices): (d1+d2)-(u1+u2)-(i1+i2)', ...
                'yLabel','(vector indices): (d1+d2)-(u1+u2)-(i1+i2)', ...
                'scale','value','cutlower',0,'cutupper',1, ...
                'tLabel',tLstr,'cType','0b','bLabel','MSIP');
        end
     end
end
