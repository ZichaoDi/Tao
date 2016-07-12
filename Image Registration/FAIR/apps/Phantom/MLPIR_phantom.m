
function [wOpt,his] = MLPIR_phantom(MLdata,minLevel,maxLevel,varargin)

% setup default parameter for parametric pre-registration
PIRopt      = @GaussNewton; % optimizer
PIRLS       = @Armijo;      % linesearch scheme
w0          = trafo('w0');  % starting guess
%  w0=[0.2,0,0]';
wStop       = w0;           % global stopping for PIR
wRef        = w0;           % regularization: (w-wRef)'*M*(w-wRef)
M           = [];           %
beta        = 0;            % regularization: H -> H + beta*I
maxIter     = 100;           % maximum number of iterations for PIR
solver      = '';           % solver for PIR
getGrid     = @getCellCenteredGrid;
% setup additional default parameter
pause       = 'off';            % flag for pauses
plots       = 1;            % flag for plots
plotIter    = 0;            % flag for output of iteration history each level
plotMLiter  = 0;            % flag for output of summarized iteration history
dimstr      = @(m) sprintf('m = [%s]',sprintf(' %d',m));

for k=1:2:length(varargin),   % overwrites default parameter
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end;
% build parameter cell for PIR
PIRpara = setOptPara(...
    'optimizer',  PIRopt,...
    'lineSearch', PIRLS,...
    'maxIter',    maxIter,...
    'Plots',      @FAIRplots2,...
    'yStop',      wStop,...
    'solver',     solver,...
    'tolG',       1e-2,...
    varargin{:});

reportStatus
FAIRmessage(sprintf('%s, minLevel=%d:maxLevel=%d',...
    'MultiLevel Parametric Image Registration',minLevel,maxLevel));
wOpt  = w0;
his   = [];
tic;

% -- for loop over all levels ---------------------------------------------
for level=minLevel:maxLevel;
    
    FAIRmessage(sprintf('%s: level %d from %d to %d, %s',...
        mfilename,level,minLevel,maxLevel,dimstr(MLdata{level}.m)));
    
    
    % get data for current level, compute interpolation coefficients
    m     = MLdata{level}.m;
    mt    = MLdata{level}.mt;
    omega=MLdata{level}.omega;
    omegat=MLdata{level}.omegat;
    T = inter('coefficients',MLdata{level}.T,[],omegat,'regularizer','moments');%,'theta',1e-1
    R = inter('coefficients',MLdata{level}.R,[],omega,'regularizer','moments');%,'theta',1e-1
    
    % update transformation
    trafo('set','omega',omega,'m',m,'omegat',omegat,'mt',mt);
    
    % initialize plots
    FAIRplots2('reset','mode','PIR-multi level','fig',level,'plots',plots);
    FAIRplots2('init',struct('Tc',T,'Rc',R,'omega',omega,'omegat',omegat,'m',m,'mt',mt));
    
    % ----- call PIR ------------------------------------
    xc   = getGrid(omega,m);
    Rc   = inter(R,omega,center(xc,m));
    xt=getCellCenteredGrid(omegat,mt);
    fctn = @(wc)  PIRobjFctnTbRs(T,Rc,omega,omegat,m,mt,beta,M,wRef,xc,xt,wc);
    if level == minLevel,
        fctn([]);   % report status
    else
        w0 = [1,2,2]'.*wOpt;  % update starting guess
    end;
    %===============================================================
%     w=linspace(-10,5,100);
%     Jc=zeros(size(w));
%     for i=1:length(w)
%         Jc(i)=feval(fctn,[w(i),w(i)]');
%     end
%     drawnow;
%     figure(12);plot(w,Jc,'r.-');hold on;
%     return;
    %===============================================================
    [wOpt,hisPIR] = optim(PIRpara{:},'objFctn',fctn,'yc',w0,'yStop',wStop);
    % ----- end PIR --------------------------------------
    
    if plotIter,
        plotIterationHistory(hisPIR,'J',[1,2,5],'fig',20+level);
    end;
    
    
    % update iteration history
    if level == minLevel,
        his.str = hisPIR.str;
        his.his = hisPIR.his;
    else
        his.his = [his.his;hisPIR.his];
    end;
    doPause(pause)
    
end;%for level
% -- for loop over all levels ---------------------------------------------
his.time = toc;

if plotMLiter,
  plotMLIterationHistory(his,'fig',30);
end;
if isempty(wStop), wStop = w0; end
% his.reduction = fctn(wOpt)/fctn(wStop);
% J = find(his.his(:,1)==-1);
% his.iter(minLevel:maxLevel) = his.his([J(2:end)-1;size(his.his,1)],1)';

FAIRmessage([mfilename,' : done !']);

%==============================================================================

function doPause(p)
if strcmp(p,'on'),
    pause;
elseif p>0,
    pause(p);
end;

