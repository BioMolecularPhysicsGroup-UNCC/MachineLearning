function colorMatrixTool(M,varargin)
% Color 2D matrix with varying thresholding rules and coloring schemes
% Written by Donald J. Jacobs at UNC Charlotte (initial code Nov 19, 2018)
% This code is based on the earlier version: ColorMatrixElementsLINEAR()
% general purpose code, used in SPLOC toolbox
%
% INPUT
% M = 2D matrix to plot as a heat map or it is a cell array of matrices, 
%     each of which will be plotted as a heat map using the exact same 
%     thresholding and coloring rules. This is the only required input. If
%     no other arguements are specified, then all default settings will be
%     applied to the processing.
%
% verbosity: {0,1,2,3} --> for plotting functions:
%             0   => quiet operation
%             1,2 => write basic information about plot to screen.
%             3   => same as (1,2) plus specify the input received.
%
% commands      purpose/options
% 'scale'       {percentile, std, value}
% 'cutlower'    sets lower threshold based on specified scale
% 'cutupper'    sets upper threshold based on specified scale
% 'morezero'    sets differential above and below zero that is set to zero
% 'yType'       {normal,reverse} to change direction of y-axis.
% 'tType'       transformation type: {'linear','rank','log10','nonlinear'}
% 'fType'       file type: {'eps','fig','png','jpg','pdf','bmp','tif'}    
% 'zSymm'       {yes,no} to force mirror symmetry about zero on color bar
% 'bLabel'      label for color bar: {no => none, tType => transformation} 
% 'dShow'       {NO,yes} to plot probability density for all input data
% 'fShow'       {YES,no} to show or not to show figures on the screen
% 'clean'       {NO,yes} to close all figures first or just append figures
% 'cType'       color spectrum options that enforce boundary conditions:
%               supports mirror symmetry: {'rwb','bwr'} 
%               defaults mirror symmetry: {'r0b','b0r'} 
%               0 at bottom  of spectrum: {'wb','0b','wr','0r'}
%               0 at the top of spectrum: {'rw','r0','bw','b0'}
% ------------------------------------------------ optional as cell arrays
% 'tLabel'      title of heat maps
% 'xLabel'      label for the x-axis
% 'yLabel'      label for the y-axis
% 'xLimit'      [xmin,xmax] to define range of INTEGER counting for x-axis
% 'yLimit'      [ymin,ymax] to define range of INTEGER counting for y-axis
% 'fName'       file name for heat map output
%
%
% PROCESS
% 1. Transform the numerical value of each pixel of the input data to a 
%    color based on thresholding and coloring rules.
% 2. If more than one matrix is given in a cell array, process all pixels
%    as a single group without assignment to a specific matrix.
% 3. from generic data, determine parameters of a transformation function.
% 4. color[pixel] = function[(numerical value of pixel)|parameters]
%
% OUTPUT
% Write to screen and/or files for each input matrix given.
%%                                                              read input
p = inputParser;
% ------------------------------------------------- define function checks
msg = 'Input the absolute value: Do not consider the sign.';
validVerbosity = @(x) isnumeric(x) && isscalar(x) && (x > -1) && (x<4);
validPercent = @(x) assert( isnumeric(x) && isscalar(x) ...
                    && (x > -1.0e-11) , msg);
expectedAnswers = {'yes','no','YES','NO'};
% -------------------------------------------------------------- verbosity
addOptional(p,'verbosity',0,validVerbosity)   % default prints out nothing
% ------------------------------------------------------------------------
expectedScales = {'percentile','std','value'};
addParameter(p,'scale','percentile', ...
             @(x) any(validatestring(x,expectedScales)));
% ------------------------------------------------------------- thresholds
addParameter(p,'cutlower',-1,validPercent);     % -1 is used as a flag
addParameter(p,'cutupper',-1,validPercent);     % -1 is used as a flag
addParameter(p,'morezero', 0 ,validPercent);        %   0 => no cutoff
% ------------------------------------------------------------------------
expectedTransform = {'linear','rank','log10','nonlinear'};
addParameter(p,'tType','linear', ...
             @(x) any(validatestring(x,expectedTransform)));
% ------------------------------------------------------------------------
% yType = normal or reverse. Y-axis points up or down respectively.
expectedYtype = {'normal','reverse'};
addParameter(p,'yType','normal', ...
             @(x) any(validatestring(x,expectedYtype)));
% ------------------------------------------------------------------------
expectedFtype = {'eps','fig','png','jpg','pdf','bmp','tif'};
addParameter(p,'fType','fig', ...
             @(x) any(validatestring(x,expectedFtype)));
% ------------------------------------------------------------------------
% zSymm = yes or no. Forces color scheme to be symmetric 
addParameter(p,'zSymm','0',@(x) any(validatestring(x,expectedAnswers)));
% ------------------------------------------------------------------------
% bLabel = tType or no. Showing the tType label on the color bar
addParameter(p,'bLabel','no');
% ------------------------------------------------------------------------
% dShow = yes or no. Showing the statistical distribution of all the data.
addParameter(p,'dShow','no',@(x) any(validatestring(x,expectedAnswers)));
% ------------------------------------------------------------------------
% fShow = yes or no. Showing figure on the screen.
addParameter(p,'fShow','yes',@(x) any(validatestring(x,expectedAnswers)));
% ------------------------------------------------------------------------
% clean = yes or no. Closing all figures first or append figures.
addParameter(p,'clean','no',@(x) any(validatestring(x,expectedAnswers)));
% ------------------------------------------------------------------------
% cType = Color scheme: 'rwb', 'bwr', 'wb', 'bw', 'wr', 'rw' (see lbmap.m)
%         Force true 0: 'r0b', 'b0r', '0b', 'b0', '0r', 'r0'  Note 0 => w. 
expectedColour = {'rwb','bwr','wb','bw','wr','rw', ...
                  'r0b','b0r','0b','b0','0r','r0'};          % 0 => weight
addParameter(p,'cType','bwr', ...
             @(x) any(validatestring(x,expectedColour)));
% ------------------------------------------------------------------------
addParameter(p,'tLabel','    ');            % initialize as no title given
% ------------------------------------------------------------------------
addParameter(p,'xLabel','    ');          % initialize as no x-label given
% ------------------------------------------------------------------------
addParameter(p,'yLabel','    ');          % initialize as no y-label given
% ------------------------------------------------- treated as cell arrays
addParameter(p,'xLimit',0);                % initialize as no xLimit given
% ------------------------------------------------------------------------
addParameter(p,'yLimit',0);                % initialize as no yLimit given
% ------------------------------------------------------------------------
addParameter(p,'fName','noWrite');  % assign default name if no name given
% ------------------------------------------------------------------------
%%                                                  parse input parameters
parse(p,varargin{:});
%%                                                  define local variables
verbosity = round(p.Results.verbosity);
scale = p.Results.scale;
cutupper = p.Results.cutupper;
cutlower = p.Results.cutlower;
morezero = p.Results.morezero;
yType = p.Results.yType;
tType = p.Results.tType;
fType = p.Results.fType;
zSymm = p.Results.zSymm;
dShow = p.Results.dShow;
fShow = p.Results.fShow;
clean = p.Results.clean;
bLabel = p.Results.bLabel;
cType = p.Results.cType;
xLimit = p.Results.xLimit;
yLimit = p.Results.yLimit;
tLabel = p.Results.tLabel;
xLabel = p.Results.xLabel;
yLabel = p.Results.yLabel;
fName = p.Results.fName;
% --------------------------------------------- translate yes/no to on/off
   if( strcmp(fShow,'yes') )
   fShow2 = 'on';
   else
   fShow2 = 'off';
   end
% --------------------------------------------------------- reset defaults
   if( strcmp(scale,'percentile') )
      if( abs(cutupper + 1) < 0.1 )
      cutupper = 100;
      end
% ---------------------------------
      if( abs(cutlower + 1) < 0.1 )
      cutlower = 100;
      end
   elseif( strcmp(scale,'std') ) % =======================================
      if( abs(cutupper + 1) < 0.1 )
      cutupper = 3;
      end
% ---------------------------------
      if( abs(cutlower + 1) < 0.1 )
      cutlower = 3;
      end
   else   % by value =====================================================
      if( abs(cutupper + 1) < 0.1 )
      error('cutupper is by value: must specify value!');
      end
% ---------------------------------------------------
      if( abs(cutlower + 1) < 0.1 )
      error('cutlower is by value: must specify value!');
      end      
   end
% ---------------------------------------------------
   if( strcmp(scale,'percentile') )
   cutupper = min(cutupper,100);
   cutlower = min(cutlower,100);
   morezero = min(morezero,100);
   end
%%                                                            error checks
   if( morezero > min(cutlower,cutupper) )
   disp('   ');
   disp('---------------------------------------');
   disp(['                cutupper = ',num2str(cutupper)]);
   disp(['                cutlower = ',num2str(cutlower)]);
   disp(['                morezero = ',num2str(morezero)]);
   disp('   ');
   error('morezero is too large');
   end
% ------------------------------------- check zSymm against color spectrum
   if( strcmp(zSymm,'yes') )
       switch cType
       case {'r0b','b0r','rwb','bwr'}       % 0 or w must be in the middle
       % do nothing -> all is consistent --> zSymm will be reconciled next
       otherwise
       error('zSymm = yes --> requires cType = {r0b, rwb, b0r, bwr}');
       end
   end
%%                                      set zSymm if undetermined thus far
   if( contains(cType,'0') )
% cType = {'r0b','b0r','0b','b0','0r','r0'}
   L = ( cType == '0' );
   cType(L) = 'w'; 
       if( contains(zSymm,'0') )
       zSymm = 'yes';                % assume YES, could change back to NO
       %else zSymm is specified as 'yes' or 'no' already.
       end
   else
       if( contains(zSymm,'0') )
       zSymm = 'no';                                   % => must remain no
       end      
   end   
% ------------------------------------- determine if zero is forced or not
   if( strcmp(zSymm,'yes') )       % also determine if zSymm changes to NO
   forceZero = true;  
%  cType = 'rwb','bwr','wb','bw','wr','rw'
      switch cType
      case {'rwb','bwr'}
      symFactor = 0;                      % zero is centered with symmetry
      case {'wb' , 'wr'}
      symFactor = -1;                      % zero is at bottom of spectrum
      zSymm = 'no';
      case {'bw' , 'rw'}
      symFactor = 1;                          % zero is at top of spectrum
      zSymm = 'no';
      otherwise
      error('unknown color spectrum');
      end
   else
   forceZero = false;  
   end   
%%                                                  set flag for color bar
   if( strcmp(bLabel,'no') )
   addColarBarString = false;
   else
   addColarBarString = true;
   end
%%                                convert to standard format of cell array
   if( iscell(M) )           % => multiple figures at once need processing
   nM = length(M);
% --------------------------------------------------------- follow through
      if( iscell(tLabel) == false )           % => duplicate default value
      str1 = tLabel;
      tLabel = cell(1,nM);
         for k=1:nM
         tLabel{k} = str1; 
         end
      else                            % make sure we have a row cell array
      temp = tLabel;
      imax = length(tLabel);
      tLabel = cell(1,imax);
         for i=1:imax
         tLabel{i} = temp{i};
         end

      end
% ------------------------------------------------ duplicate default value 
      if( iscell(xLabel) == false )     
      str1 = xLabel;
      xLabel = cell(1,nM);
         for k=1:nM
         xLabel{k} = str1; 
         end
      else                            % make sure we have a row cell array
      temp = xLabel;
      imax = length(xLabel);
      xLabel = cell(1,imax);
         for i=1:imax
         xLabel{i} = temp{i};
         end
      end
% ------------------------------------------------ duplicate default value 
      if( iscell(yLabel) == false )     
      str1 = yLabel;
      yLabel = cell(1,nM);
         for k=1:nM
         yLabel{k} = str1; 
         end
      else                            % make sure we have a row cell array
      temp = yLabel;
      imax = length(yLabel);
      yLabel = cell(1,imax);
         for i=1:imax
         yLabel{i} = temp{i};
         end

      end
% ------------------------------------------------ duplicate default value 
      if( iscell(xLimit) == false )  
      temp = xLimit;
      xLimit = cell(1,nM);
         for k=1:nM
         xLimit{k} = temp; 
         end
      else                            % make sure we have a row cell array
      temp = xLimit;
      imax = length(xLimit);
      xLimit = cell(1,imax);
         for i=1:imax
         xLimit{i} = temp{i};
         end
      end
% ------------------------------------------------ duplicate default value 
      if( iscell(yLimit) == false )  
      temp = yLimit;
      yLimit = cell(1,nM);
         for k=1:nM
         yLimit{k} = temp; 
         end
      else                            % make sure we have a row cell array
      temp = yLimit;
      imax = length(yLimit);
      yLimit = cell(1,imax);
         for i=1:imax
         yLimit{i} = temp{i};
         end
      end
% ------------------------------------------------ duplicate default value 
      if( iscell(fName) == false )  
      baseName = fName;
      fName = cell(1,nM);
         for k=1:nM
         fName{k} = [baseName,'_',num2str(k)]; 
         end
      else                            % make sure we have a row cell array
      temp = fName;
      imax = length(fName);
      fName = cell(1,imax);
         for i=1:imax
         fName{i} = temp{i};
         end
      end
   clear temp str1
   else               % => convert M & associated variables to cell arrays
   temp = M;
   M = cell(1,1);
   M{1} = temp;
% --------------------------------------------------------- follow through 
   temp = tLabel;
   tLabel = cell(1,1);
   tLabel{1} = temp; 
% ---------------------------------------------------------
   temp = xLabel;
   xLabel = cell(1,1);
   xLabel{1} = temp; 
% ---------------------------------------------------------
   temp = yLabel;
   yLabel = cell(1,1);
   yLabel{1} = temp; 
% ---------------------------------------------------------
   temp = xLimit;
   xLimit = cell(1,1);
   xLimit{1} = temp; 
% ---------------------------------------------------------
   temp = yLimit;
   yLimit = cell(1,1);
   yLimit{1} = temp;
% ---------------------------------------------------------
   temp = fName;
   fName = cell(1,1);
   fName{1} = temp;
   clear temp
   nM = 1;
   end
%%                                         summarize parameters by request
   if( verbosity == 3 )
   disp(p.Results);
   end
   if( verbosity > 0 )
   disp('   ');
   disp('----------------------------------------');
   disp(['                verbosity = ',num2str(verbosity)]);
   disp(['                    scale = ',scale]);
   disp(['                 cutupper = ',num2str(cutupper)]);
   disp(['                 cutlower = ',num2str(cutlower)]);
   disp(['                 morezero = ',num2str(morezero)]);
   disp(['    y-axis direction type = ',yType]);
   disp(['           transformation = ',tType]);
   disp(['                file type = ',fType]);
   disp(['color bar mirror symmetry = ',zSymm]);
   disp(['           show bar label = ',bLabel]);   
   disp(['         show data as pdf = ',dShow]);
   disp(['   show figures on screen = ',fShow]);
   disp(['  clean figures on screen = ',clean]);
   disp(['             color scheme = ',cType]);
% ------------------------------------------- create string for x,y limits
   xL = cell(1,nM);
   yL = cell(1,nM);
      for k=1:nM
        if( isscalar(xLimit{k}) )
        xL{k} = 'indices';
        else
        xL{k} = ['[',num2str(xLimit{k}(1)),',',num2str(xLimit{k}(2)),']'];
        end
        if( isscalar(yLimit{k}) )
        yL{k} = 'indices';
        else
        yL{k} = ['[',num2str(yLimit{k}(1)),',',num2str(yLimit{k}(2)),']'];
        end   
      end
% ------------------------------------------- 
      if( nM == 1 )
      disp(['                 x limits = ',xL{1}]);
      disp(['                 y limits = ',yL{1}]);
      disp(['              title label = ',tLabel{1}]);
      disp(['                  x label = ',xLabel{1}]);
      disp(['                  y label = ',yLabel{1}]);
      disp(['                file name = ',fName{1}]);
      else
      disp('----------------------------------------');   
      cellArray = [fName',tLabel',xLabel',yLabel',xL',yL'];
      T = cell2table(cellArray, ...
          'VariableNames',{'fName','tLabel','xLabel','yLabel', ...
                           'xLimits','yLimits'});
      disp(T);
      end
   disp('----------------------------------------');
   end
%%                               consider cleaning out all current figures
   if( strcmp(clean,'yes') )
   close all
   end
%%                                                 combine all data into q
q = [];
   for k=1:nM
   q = horzcat(M{k}(1:end),q);
   end
q = sort(q);
nq = length(q);
%%                                   write relevant information by request
   if( verbosity > 0 )
   disp(['        number of figures = ',num2str(nM)]);
   disp(['         number of pixels = ',num2str(nq)]);
   disp('----------------------------------------');
   end
%%                                                               process q
qave = mean(q);
qstd = std(q);
% ----------------------- need absolute values to work with data near zero
abs_q = sort(abs(q));
% --------------------------------------------------- define the zero band
setZero = 1.0e-12;                           % this is considered a true 0
   if( strcmp(scale,'percentile') )
   nzmax = round(morezero*nq/100);         % maximum # for negative values
      if( nzmax > 0 )
      max_abs_qzero = abs_q(nzmax);
      else
      max_abs_qzero = 0;
      end
   elseif( strcmp(scale,'std') )
   max_abs_qzero = morezero*qstd;
   else                                                   % => value based
   max_abs_qzero = abs(morezero);
   end
max_abs_qzero = max(max_abs_qzero,setZero);                         
% ---------------------------------------------- split data into 3 regions
Lp = ( q > max_abs_qzero );
qp = q(Lp);
Ln = ( q < -max_abs_qzero ); 
qn = q(Ln);
Lz = and(~Lp,~Ln);
qz = q(Lz);
% --------------------------------- number of positive and negative pixels
npp = length(qp);
nnp = length(qn);
nzp = length(qz);
% ----------------------------------------------- determine special values
   if( nzp > 0 )
   qzAve = mean(qz);
   qzStd = std(qz);
   end
   if( strcmp(scale,'percentile') )
      if( npp > 0 )
      npmax = round(cutupper*npp/100);     % maximum # for positive values 
      qmax = qp(npmax);                  % already sorted from low to high
      qpAve = mean(qp);
      qpStd = std(qp);
      else
      qmax = 0;
      end
% -----------------------------------
      if( nnp > 0 )
      nnmax = round(cutlower*nnp/100);     % maximum # for negative values
      sneg = sort(abs(qn));            % sneg changes the sign and resorts
      qmin = -sneg(nnmax);
      qnAve = -mean(sneg);
      qnStd = std(sneg);
      else
      qmin = 0;
      end
   elseif( strcmp(scale,'std') )
      if( npp > 0 )
      qpAve = mean(qp);
      qpStd = std(qp);
      qmax = cutupper*qstd;
      else
      qmax = 0;
      end
% --------------------------------
      if( nnp > 0 )
      sneg = sort(abs(qn));            % sneg changes the sign and resorts
      qnAve = -mean(sneg);
      qnStd = std(sneg);
      qmin = -cutlower*qstd;    
      else
      qmin = 0;
      end  
   else                                                         % by value
      if( npp > 0 )
      qpAve = mean(qp);
      qpStd = std(qp);
      qmax = cutupper;                   % already sorted from low to high
      else
      qmax = 0;
      end
% --------------------------------
      if( nnp > 0 )
      sneg = sort(abs(qn));            % sneg changes the sign and resorts
      qnAve = -mean(sneg);
      qnStd = std(sneg);
      qmin = -cutlower;                         % value has to be negative
      else
      qmin = 0;
      end
   end
%%                                   write relevant information by request
   if( verbosity > 0 )
   disp(['         mean pixel value = ',num2str(qave)]);
   disp(['       std of pixel value = ',num2str(qstd)]);
   disp('   ');
   disp(['     # of (+)pixel values = ',num2str(npp)]);
      if( npp > 0 )
      disp(['     mean (+)pixel values = ',num2str(qpAve)]);
      disp(['      std (+)pixel values = ',num2str(qpStd)]);
      end
   disp([' cutupper (+)pixel  value = ',num2str(qmax)]);
   disp('   ');
   disp(['    # of (~0)pixel values = ',num2str(nzp)]);
      if( nzp > 0 )
      disp(['    mean (~0)pixel values = ',num2str(qzAve)]);
      disp(['     std (~0)pixel values = ',num2str(qzStd)]);
      end
   disp(['morezero (~0)pixel  value = ',num2str(max_abs_qzero)]);
   disp('   ');
   disp(['     # of (-)pixel values = ',num2str(nnp)]);
      if( nnp > 0 )
      disp(['     mean (-)pixel values = ',num2str(qnAve)]);
      disp(['      std (-)pixel values = ',num2str(qnStd)]);
      end  
   disp([' cutlower (-)pixel  value = ',num2str(qmin)]);
   disp('----------------------------------------');
      if( forceZero == true )
      disp('            force 0 level = TRUE');
      else
      disp('            force 0 level = FALSE');
      end
   disp('   ');
   end
%%                                  plot pdf for q when dShow is requested
   if( strcmp(dShow,'yes') )
   figure;
   figure_number = get(gcf,'Number');
   nv = nq/5;
      if( nv > 4 )
      q2 = max(q);
      q1 = min(q);
      dq = 0.04*(q2 - q1);
      q1 = q1 - dq;
      q2 = q2 + dq;
      dq = (q2 - q1)/(nq/5);
      x = q1:dq:q2;
      [pdf,x] = ksdensity(q,x);
      figure(figure_number);
      clf;
      plot(x,pdf,'b','linewidth',1.5);
      xlabel('pixel value');
      ylabel('probability density');
      title('calculation is now paused');
      beep
      pause
      end
   end
%%                                      transform M -> Cplot -> plot image
figure;
figure_number = get(gcf,'Number');
   for k=1:nM    
   hf = figure(figure_number);
   set(hf, 'Visible',fShow2);
   clf
   Cplot = M{k};     % elements of Cplot will be transformed based on rule
% ===================================================== start thresholding
% WARNING: thresholding takes precedence before performing transformations
   L0 = and( (Cplot < max_abs_qzero) , ...          % define band of zeros
             (Cplot > -max_abs_qzero) );
   Cplot(L0) = 0;                         % apply thresholding around zero
   Lp = ( Cplot > qmax );       % define positive limit not to be exceeded
   Cplot(Lp) = qmax;                            % apply upper thresholding
   Ln = ( Cplot < qmin );       % define negative limit not to be exceeded
   Cplot(Ln) = qmin;                            % apply lower thresholding
% ======================================================= end thresholding
% 'tType'       transformation type: {'linear','rank','log10','nonlinear'}
      if( forceZero )
% --------------------------------------------------- apply transformation
        switch tType
        case 'linear'
        cmax = qmax;
        cmin = qmin;
        cbStr = 'linear';
        case 'rank'
        xu = unique( abs(q) );
        xu = horzcat( sort(-xu), xu );
        ju = length(xu);
        yu = (0:ju-1)/ju;
        yu = 2*(yu - 0.5);
        yu = 99.9*yu;
        Cplot = interp1(xu,yu,Cplot,'linear','extrap');
        cmin = -100; 
        cmax = 100;
        cbStr = 'rank';
        case 'log10'
        Lp = ( Cplot > 0 );
        Cplot(Lp) = log10(Cplot(Lp) + 1);
        Ln = ( Cplot < 0 );
        Cplot(Ln) = -log10( abs(Cplot(Ln)) + 1);
        cmax = log10(qmax + 1);
        cmin = -log10( abs(qmin) + 1);
        cbStr = 'piecewise: -log10( 1+|m| ) _0_ log10( 1+m )';
        case 'nonlinear'
        qmax = max(abs(qmin),qmax);
        fixedPoint = 0.7616*qmax; 
        Cplot = qmax*tanh(Cplot/fixedPoint);    % 0.7616 is a ~fixed point
        cmax = qmax;
        cmin = qmin;
        cbStr = 'nonlinear: Tanh()';        
        otherwise
        error('unknown transformation is being requested')
        end
% ------------------------------------------------------ force zero at ...
        if( symFactor < 0 )                      % ... bottom of sprectrum
           if( cmin < 0 )
           disp(['    title = ',tLabel{k}]);
           disp(['file name = ',fName{k}]);
           disp(['cmax = ',num2str(cmax)]);
           disp(['cmin = ',num2str(cmin)]);
           error('negative values detected: use rwb or bwr spectrum');
           end
        cmin = 0;
        elseif( symFactor > 0 )                      % ... top of spectrum
           if( cmax > 0 )
           disp(['    title = ',tLabel{k}]);
           disp(['file name = ',fName{k}]);
           disp(['cmin = ',num2str(cmin)]);
           disp(['cmax = ',num2str(cmax)]);
           error('positive values detected: use rwb or bwr spectrum');
           end
        cmax = 0;
        else                                   % ... in middle of spectrum
        temp = max( cmax, abs(cmin) );
        cmax = temp;
        cmin = -temp;
        end
      else                                            % let the zero float
% --------------------------------------------------- apply transformation
        switch tType
        case 'linear'
        cmax = qmax;
        cmin = qmin;
        cbStr = 'linear';
        case 'rank'
        xu = unique(q);
        ju = length(xu);
        yu = 1:ju;
        yu = 99.9*(yu - 0.5)/ju;
        Cplot = interp1(xu,yu,Cplot,'linear','extrap');
        cmin = 0; 
        cmax = 100;
        cbStr = 'rank';
        case 'log10'
        Lp = ( Cplot > 0 );
        Cplot(Lp) = log10(Cplot(Lp) + 1);
        Ln = ( Cplot < 0 );
        Cplot(Ln) = -log10( abs(Cplot(Ln)) + 1);
        cmax = log10(qmax + 1);
        cmin = -log10( abs(qmin) + 1);
        cbStr = 'piecewise: -log10( 1+|m| ) _0_ log10( 1+m )';
        case 'nonlinear'
        qmax = max(abs(qmin),qmax);
        fixedPoint = 0.7616*qmax; 
        Cplot = qmax*tanh(Cplot/fixedPoint);    % 0.7616 is a ~fixed point
        cmax = qmax;
        cmin = qmin;
        cbStr = 'nonlinear: Tanh()'; 
        otherwise
        error('unknown transformation is being requested')
        end   
      end
% --------------------------------------------- reset cbStr when necessary
      if( strcmp(bLabel,'tType') )
      bLabel = cbStr;                      % sets bar label to method name
      end     
% ---------------------------------------------------------- make the plot     
      if( isscalar( xLimit{k}) )            % => no specified x-axis range
         if( isscalar( yLimit{k}) )         % => no specified y-axis range
         [jy,jx] = size(Cplot);                 % rows => y   columns => x
         xx = [1,jx];
         yy = [1,jy];
         else                                  % => y-axis range specified
         yy = yLimit{k};
         [~,jx] = size(Cplot);                  % rows => y   columns => x
         xx = [1,jx];
         end
      else                                     % => x-axis range specified
      xx = xLimit{k};
         if( isscalar( yLimit{k}) )         % => no specified y-axis range
         [jy,~] = size(Cplot);                  % rows => y   columns => x
         yy = [1,jy];
         else                                  % => y-axis range specified
         yy = yLimit{k};
         end
      end 
   imagesc('XData',xx,'YData',yy,'CData',Cplot); 
   %imagesc(Cplot);                            
   nc = 256;
      if( npp == 0 )
      nc = nc/2;
      end
      if( nnp == 0 )
      nc = nc/2;
      end
   caxis([cmin,cmax]);
   set(gca,'YDir',yType);
   %colormap( lbmap(nc,cType) )           % cType defines the color scheme
% --------------------------------------------------------- generalization
   % example linear transformation on basic color scale
   lighten = 0.0;                                                 % FIX ME
   cRule = lighten*[1,1,1] + (1 - lighten)*lbmap(nc,cType);       % FIX ME
   % cRule can be used here to make any type of transformation on color
   % also should modify lbmap ..... to make relevant.               FIX ME
   colormap(cRule);
      if( addColarBarString )
      c = colorbar;
      c.Label.String = bLabel;
      else
      colorbar;
      end
%    figure(figure_number)
%    pause
% ------------------------------------------------------------------------
   title(tLabel{k},'Interpreter','none');
   xlabel(xLabel{k},'Interpreter','none');
   ylabel(yLabel{k},'Interpreter','none');
%    figure(figure_number)
%    pause
   xx0 = xx;
   yy0 = yy;
   xx(1) = min(xx0) - 0.5;
   xx(2) = max(xx0) + 0.5;
   yy(1) = min(yy0) - 0.5;
   yy(2) = max(yy0) + 0.5;
   xlim(xx);
   ylim(yy);
%    figure(figure_number)
%    pause
      if( contains(fName,'noWrite') == 0 )
          if( strcmp(fType,'eps') )
          saveas(gca,fName{k},'psc2')
          else
          saveas(gca,fName{k},fType)
          end
      end
   figure_number = figure_number + 1;
   end 
%%
function map = lbmap(n,scheme)
%LBMAP Returns specified Light-Bertlein colormap.
%
%   LBMAP(N,SCHEME) returns an Nx3 colormap. SCHEME can be one of the
%   following strings:
%
%    'Blue'       Single-hue progression to purlish-blue (default)
%    'BlueGray'   Diverging progression from blue to gray
%    'BrownBlue'  Orange-white-purple diverging scheme
%    'RedBlue'    Modified spectral scheme
%================================================= modifed by Dr. Jacobs
%    'BWR'        Blue-White-Red diverging scheme
%    'RWB'        Red-White-Blue diverging scheme
%    'WB'         White-Blue diverging scheme
%    'BW'         Blue-White diverging scheme
%    'RW'         Red-White diverging scheme
%    'WR'         White-Red diverging scheme
%
%   If N is not specified, the size of the colormap is determined by the
%   current figure. If no figure exists, MATLAB creates one.
%
%Example 1: 7-color single-hue blue (default)
%   load penny
%   imagesc(P)
%   colormap(lbmap(7))
%   colorbar
%
%Example 2: 11-color modified spectrum
%   load penny
%   imagesc(P)
%   colormap(lbmap(11,'RedBlue'))
%   colorbar
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, COLORMAP, RGBPLOT.

% Reference:
% A. Light & P.J. Bartlein, "The End of the Rainbow? Color Schemes for
% Improved Data Graphics," Eos,Vol. 85, No. 40, 5 October 2004.
% http://geography.uoregon.edu/datagraphics/EOS/Light&Bartlein_EOS2004.pdf

% Copyright 2007-2010 The MathWorks, Inc.

%defensive programming
narginchk(0,2)
nargoutchk(0,1)

%defaults
if nargin<2
  scheme = 'Blue';
end
if nargin<1
  n = size(get(gcf,'colormap'),1);
end

%valid schemes
switch lower(scheme)
  case 'blue'
    baseMap = BlueMap;
  case 'bluegray'
    baseMap = BlueGrayMap;
  case 'brownblue'
    baseMap = BrownBlueMap;
  case 'redblue'
    baseMap = RedBlueMap;
  case 'bwr'
    baseMap = BlueWhiteRedMap;
  case 'rwb'
    baseMap = RedWhiteBlueMap;
  case 'wr'
    baseMap = WhiteRedMap;
  case 'rw'
    baseMap = RedWhiteMap;
  case 'wb'
    baseMap = WhiteBlueMap;
  case 'bw'
    baseMap = BlueWhiteMap;
  otherwise
    error(['Invalid scheme ' scheme])
end
idx1 = linspace(0,1,size(baseMap,1));
idx2 = linspace(0,1,n);
map = interp1(idx1,baseMap,idx2);
% ------------------------------------------------------------------ RGB

function baseMap = BlueMap
baseMap = [243 246 248;
           224 232 240;
           171 209 236;
           115 180 224;
            35 157 213;
             0 142 205;
             0 122 192]/255;

function baseMap = BlueGrayMap
%DivergingBlueGray
baseMap = [  0 170 227;
            53 196 238;
           133 212 234;
           190 230 242;
           217 224 230;
           146 161 170;
           109 122 129;
            65  79  81]/255;

function baseMap = BrownBlueMap
baseMap = [144 100  44;
           187 120  54;
           225 146  65;
           248 184 139;
           244 218 200;
           241 244 245;
           207 226 240;
           160 190 225;
           109 153 206;
            70  99 174;
            24  79 162]/255;

function baseMap = RedBlueMap
baseMap = [175  53  71;
           216  82  88;
           239 133 122;
           245 177 139;
           249 216 168;
           242 238 197;
           216 236 241;
           154 217 238;
            68 199 239;
             0 170 226;
             0 116 188]/255;
         
function baseMap = BlueWhiteRedMap
baseMap = [  0   0  70;
             0   0 130;
             0   0 180;
             0   0 220;
             3   7 240;
            11  21 250;
            41  51 255;
            61 101 255;
           102 140 255;
           153 200 255;
           204 229 255;
           225 245 255;
           245 250 255;
           255 255 255;
           255 250 245;
           255 245 225;
           255 229 204;
           255 200 153;
           255 140 102;
           255 101  61;
           255  51  41;
           250  21  11;
           240   7   3;
           220   0   0;
           180   0   0;
           130   0   0]/255;
       
function baseMap = RedWhiteBlueMap
baseMap = [130   0   0;
           180   0   0
           220   0   0;
           240   7   3;
           250  21  11;
           255  51  41;
           255 101  61;
           255 140 102;
           255 200 153;
           255 229 204;
           255 245 225;
           255 250 245;
           255 255 255;
           245 250 255;
           225 245 255;
           204 229 255;
           153 200 255;
           102 140 255;
            61 101 255;
            41  51 255;
            11  21 250;
             3   7 240;
             0   0 220;
             0   0 180;
             0   0 130;   
             0   0  70]/255;
%------------------------------------
function baseMap = WhiteRedMap
baseMap = [255 255 255;
           255 250 245;
           255 245 225;
           255 229 204;
           255 200 153;
           255 140 102;
           255 101  61;
           255  51  41;
           250  21  11;
           240   7   3;
           220   0   0;
           180   0   0;
           130   0   0]/255;
       
function baseMap = RedWhiteMap
baseMap = [130   0   0;
           180   0   0
           220   0   0;
           240   7   3;
           250  21  11;
           255  51  41;
           255 101  61;
           255 140 102;
           255 200 153;
           255 229 204;
           255 245 225;
           255 250 245;
           255 255 255;]/255;
             
function baseMap = WhiteBlueMap
baseMap = [255 255 255;
           245 250 255;
           225 245 255;
           204 229 255;
           153 200 255;
           102 140 255;
            61 101 255;
            41  51 255;
            11  21 250;
             3   7 240;
             0   0 220;
             0   0 180;
             0   0 130;   
             0   0  70]/255;
       
function baseMap = BlueWhiteMap
baseMap = [  0   0  70;
             0   0 130;
             0   0 180;
             0   0 220;             
             3   7 240;
            11  21 250;
            41  51 255;
            61 101 255;
           102 140 255;
           153 200 255;
           204 229 255;
           225 245 255;   
           245 250 255;
           255 255 255]/255;
