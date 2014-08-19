% ==================================================================================
% (c) Lars Ruthotto 2011/02/08, see FAIR.2 and FAIRcopyright.m.
% http://www.mic.uni-luebeck.de/people/lars-ruthotto.html
%
% function [G,dG] = geometry(yc,m,flag,varargin)
%
% area, volume and det computation for tetrahedral partition of yc
%
% Input:
% ------
%    yc     	nodal grid
%    m     		number of discretization points
%    flag		can be {'A','V','Jac'}
%    varargin   optional parameter (see below)
%
% Output:
% ------
%    G   		value
%    dG  		derivative 
%   if ~matrixFree, dG is sparse matrix; else, dG is struct endif		
%
% See also
%  @article{2011-BMR,
%	Author = {Burger M., Modersitzki J., Ruthotto L. },
%	Publisher = {University of Muenster},
%	Title = {A hyperelastic regularization energy for image registration},
%	Year = {2011}
%  }
% 
% see also hyperelatic and BigTutorialRegularizer
% ==================================================================================
function [G,dG] = geometry(yc,m,flag,varargin)
if nargin==0,
    testThisMfile;
    return;
end
matrixFree = regularizer('get','matrixFree');
doDerivative = (nargout>1);
dG = [];

for k=1:2:length(varargin) % overwrites defaults
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

dim = length(m);
if matrixFree
  flag = [flag '-mf-' num2str(dim)];
else
  flag = [flag '-mb-' num2str(dim)];
end

n  = prod(m);
nn = prod(m+1);
if exist('omega','var') && not(isempty(omega)),
    hd = prod((omega(2:2:end)-omega(1:2:end)) ./m);
end
  
switch flag
    case 'V-mb-2'
        S  = speye(prod(m+1),prod(m+1));
        PA = S(1:end-(m(1)+1),:); PA((m(1)+1)*(1:m(2)),:) = [];
        PB = S(1:end-(m(1)+1),:); PB((m(1)+1)*(0:m(2)-1)+1,:) = [];
        PC = S(m(1)+2:end,:);     PC((m(1)+1)*(1:m(2)),:) = [];
        PD = S(m(1)+2:end,:);     PD((m(1)+1)*(0:m(2)-1)+1,:) = [];
        
        PM = (1/4) * (PA + PB + PC + PD);
        
        [volABM dvolABM] = arTriangle2D(PA,PB,PM,yc,doDerivative);
        [volBDM dvolBDM] = arTriangle2D(PB,PD,PM,yc,doDerivative);
        [volDCM dvolDCM] = arTriangle2D(PD,PC,PM,yc,doDerivative);
        [volCAM dvolCAM] = arTriangle2D(PC,PA,PM,yc,doDerivative);
        
        G = [volABM ; volBDM; volDCM; volCAM];
        if doDerivative,
            dG =[dvolABM ; dvolBDM; dvolDCM; dvolCAM];
        end
    case 'V-mb-3'
        % build edge projectors
        PA = Pedge3D('A',m);
        PB = Pedge3D('B',m);
        PC = Pedge3D('C',m);
        PD = Pedge3D('D',m);
        PE = Pedge3D('E',m);
        PF = Pedge3D('F',m);
        PG = Pedge3D('G',m);
        PH = Pedge3D('H',m);
        PABDC = Pedge3D('ABDC',m);
        PBDHF = Pedge3D('BDHF',m);
        PABFE = Pedge3D('ABFE',m);
        PCAEG = Pedge3D('CAEG',m);
        PDCGH = Pedge3D('DCGH',m);
        PEFHG = Pedge3D('EFHG',m);
        PM = Pedge3D('M',m);
        
        % compute volumes of 24 tetras per voxel
        [V1, dV1] = volTetra3D(yc,PA,PB,PABDC,PM,doDerivative);
        [V2, dV2] = volTetra3D(yc,PB,PD,PABDC,PM,doDerivative);
        [V3, dV3] = volTetra3D(yc,PD,PC,PABDC,PM,doDerivative);
        [V4, dV4] = volTetra3D(yc,PC,PA,PABDC,PM,doDerivative);
        
        [V5, dV5] = volTetra3D(yc,PB,PD,PBDHF,PM,doDerivative);
        [V6, dV6] = volTetra3D(yc,PD,PH,PBDHF,PM,doDerivative);
        [V7, dV7] = volTetra3D(yc,PH,PF,PBDHF,PM,doDerivative);
        [V8, dV8] = volTetra3D(yc,PF,PB,PBDHF,PM,doDerivative);
        
        [V9, dV9] = volTetra3D(yc,PA,PB,PABFE,PM,doDerivative);
        [V10, dV10] = volTetra3D(yc,PB,PF,PABFE,PM,doDerivative);
        [V11, dV11] = volTetra3D(yc,PF,PE,PABFE,PM,doDerivative);
        [V12, dV12] = volTetra3D(yc,PE,PA,PABFE,PM,doDerivative);
        
        [V13, dV13] = volTetra3D(yc,PC,PA,PCAEG,PM,doDerivative);
        [V14, dV14] = volTetra3D(yc,PA,PE,PCAEG,PM,doDerivative);
        [V15, dV15] = volTetra3D(yc,PE,PG,PCAEG,PM,doDerivative);
        [V16, dV16] = volTetra3D(yc,PG,PC,PCAEG,PM,doDerivative);
        
        [V17, dV17] = volTetra3D(yc,PD,PC,PDCGH,PM,doDerivative);
        [V18, dV18] = volTetra3D(yc,PC,PG,PDCGH,PM,doDerivative);
        [V19, dV19] = volTetra3D(yc,PG,PH,PDCGH,PM,doDerivative);
        [V20, dV20] = volTetra3D(yc,PH,PD,PDCGH,PM,doDerivative);
        
        [V21, dV21] = volTetra3D(yc,PE,PF,PEFHG,PM,doDerivative);
        [V22, dV22] = volTetra3D(yc,PF,PH,PEFHG,PM,doDerivative);
        [V23 dV23] = volTetra3D(yc,PH,PG,PEFHG,PM,doDerivative);
        [V24, dV24] = volTetra3D(yc,PG,PE,PEFHG,PM,doDerivative);
        
        G = [-V1;-V2;-V3;-V4;V5;V6;V7;V8;V9;V10;V11;V12;V13;V14;V15;V16;V17;V18;V19;V20;V21;V22;V23;V24];
        if doDerivative,
            dG = [-dV1;-dV2;-dV3;-dV4;dV5;dV6;dV7;dV8;dV9;dV10;dV11;dV12;dV13;dV14;...
                dV15;dV16;dV17;dV18;dV19;dV20;dV21;dV22;dV23;dV24];
        end
    case {'V-mf-2','V-mf-3'}
        G = geometrymexC(yc(:),m,'V');
    case {'VRange-mb-2','VRange-mb-3'}
        vol = geometry(yc,m,'V');
        G = zeros(2,1);
        G(1) = min(vol(:));
        G(2) = max(vol(:));
    case {'VRange-mf-2','VRange-mf-3'}
        G = geometrymexC(yc(:),m,'VRange');
     case 'Jac-mb-2'
        S  = speye(prod(m+1),prod(m+1));
        PA = S(1:end-(m(1)+1),:); PA((m(1)+1)*(1:m(2)),:) = [];
        PB = S(1:end-(m(1)+1),:); PB((m(1)+1)*(0:m(2)-1)+1,:) = [];
        PC = S(m(1)+2:end,:);     PC((m(1)+1)*(1:m(2)),:) = [];
        PD = S(m(1)+2:end,:);     PD((m(1)+1)*(0:m(2)-1)+1,:) = [];
        
        PM = (1/4) * (PA + PB + PC + PD);
        
        [volABM dvolABM] = arTriangle2D(PA,PB,PM,yc,doDerivative);
        [volBDM dvolBDM] = arTriangle2D(PB,PD,PM,yc,doDerivative);
        [volDCM dvolDCM] = arTriangle2D(PD,PC,PM,yc,doDerivative);
        [volCAM dvolCAM] = arTriangle2D(PC,PA,PM,yc,doDerivative);
        
        G = (volABM + volBDM + volDCM + volCAM) /hd;
        if doDerivative,
            dG =(dvolABM + dvolBDM + dvolDCM + dvolCAM) / hd;
        end
    case {'Jac-mf-2','Jac-mf-3'}
        G = geometrymexC(yc(:),m,'Jac')/hd;
        if doDerivative,
            dG.dJac = @(y,m,x) geometrymexC(y(:),m,'dJacx',x(:))/hd;
            dG.dJacadj = @(y,m,x) geometrymexC(y(:),m,'dJacadjx',x(:))/hd;
        end
    case {'dJacx-mf-2','dJacx-mf-3'}
        G = geometrymexC(yc(:),m,'dJacx',x(:));
    case {'dJacadjx-mf-2','dJacadjx-mf-3'}
        G = geometrymexC(yc(:),m,'dJacadjx',x(:));
       
    case 'Jac-mb-3'
        % build edge projectors
        PA = Pedge3D('A',m);
        PB = Pedge3D('B',m);
        PC = Pedge3D('C',m);
        PD = Pedge3D('D',m);
        PE = Pedge3D('E',m);
        PF = Pedge3D('F',m);
        PG = Pedge3D('G',m);
        PH = Pedge3D('H',m);
        PABDC = Pedge3D('ABDC',m);
        PBDHF = Pedge3D('BDHF',m);
        PABFE = Pedge3D('ABFE',m);
        PCAEG = Pedge3D('CAEG',m);
        PDCGH = Pedge3D('DCGH',m);
        PEFHG = Pedge3D('EFHG',m);
        PM = Pedge3D('M',m);
        
        % compute volumes of 24 tetras per voxel
        [V1, dV1] = volTetra3D(yc,PA,PB,PABDC,PM,doDerivative);
        [V2, dV2] = volTetra3D(yc,PB,PD,PABDC,PM,doDerivative);
        [V3, dV3] = volTetra3D(yc,PD,PC,PABDC,PM,doDerivative);
        [V4, dV4] = volTetra3D(yc,PC,PA,PABDC,PM,doDerivative);
        
        [V5, dV5] = volTetra3D(yc,PB,PD,PBDHF,PM,doDerivative);
        [V6, dV6] = volTetra3D(yc,PD,PH,PBDHF,PM,doDerivative);
        [V7, dV7] = volTetra3D(yc,PH,PF,PBDHF,PM,doDerivative);
        [V8, dV8] = volTetra3D(yc,PF,PB,PBDHF,PM,doDerivative);
        
        [V9, dV9] = volTetra3D(yc,PA,PB,PABFE,PM,doDerivative);
        [V10, dV10] = volTetra3D(yc,PB,PF,PABFE,PM,doDerivative);
        [V11, dV11] = volTetra3D(yc,PF,PE,PABFE,PM,doDerivative);
        [V12, dV12] = volTetra3D(yc,PE,PA,PABFE,PM,doDerivative);
        
        [V13, dV13] = volTetra3D(yc,PC,PA,PCAEG,PM,doDerivative);
        [V14, dV14] = volTetra3D(yc,PA,PE,PCAEG,PM,doDerivative);
        [V15, dV15] = volTetra3D(yc,PE,PG,PCAEG,PM,doDerivative);
        [V16, dV16] = volTetra3D(yc,PG,PC,PCAEG,PM,doDerivative);
        
        [V17, dV17] = volTetra3D(yc,PD,PC,PDCGH,PM,doDerivative);
        [V18, dV18] = volTetra3D(yc,PC,PG,PDCGH,PM,doDerivative);
        [V19, dV19] = volTetra3D(yc,PG,PH,PDCGH,PM,doDerivative);
        [V20, dV20] = volTetra3D(yc,PH,PD,PDCGH,PM,doDerivative);
        
        [V21, dV21] = volTetra3D(yc,PE,PF,PEFHG,PM,doDerivative);
        [V22, dV22] = volTetra3D(yc,PF,PH,PEFHG,PM,doDerivative);
        [V23 dV23] = volTetra3D(yc,PH,PG,PEFHG,PM,doDerivative);
        [V24, dV24] = volTetra3D(yc,PG,PE,PEFHG,PM,doDerivative);
        
        G = (-V1-V2-V3-V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14+V15+V16+V17+V18+V19+V20+V21+V22+V23+V24) /hd;
        if doDerivative,
            dG = (-dV1-dV2-dV3-dV4+dV5+dV6+dV7+dV8+dV9+dV10+dV11+dV12+dV13+dV14+...
                dV15+dV16+dV17+dV18+dV19+dV20+dV21+dV22+dV23+dV24) /hd;
        end
     case 'A-mb-3'
        PA = Pedge3D('A',m);
        PB = Pedge3D('B',m);
        PC = Pedge3D('C',m);
        PD = Pedge3D('D',m);
        PE = Pedge3D('E',m);
        PF = Pedge3D('F',m);
        PG = Pedge3D('G',m);
        PH = Pedge3D('H',m);
        PABDC = Pedge3D('ABDC',m);
        PBDHF = Pedge3D('BDHF',m);
        PABFE = Pedge3D('ABFE',m);
        PCAEG = Pedge3D('CAEG',m);
        PDCGH = Pedge3D('DCGH',m);
        PEFHG = Pedge3D('EFHG',m);
        
        [A1,  dA1] = arTriangle3D(yc,PA,PB,PABDC,doDerivative);
        [A2,  dA2] = arTriangle3D(yc,PB,PD,PABDC,doDerivative);
        [A3,  dA3] = arTriangle3D(yc,PD,PC,PABDC,doDerivative);
        [A4,  dA4] = arTriangle3D(yc,PC,PA,PABDC,doDerivative);
        
        [A5,  dA5] = arTriangle3D(yc,PB,PD,PBDHF,doDerivative);
        [A6,  dA6] = arTriangle3D(yc,PD,PH,PBDHF,doDerivative);
        [A7,  dA7] = arTriangle3D(yc,PH,PF,PBDHF,doDerivative);
        [A8,  dA8] = arTriangle3D(yc,PF,PB,PBDHF,doDerivative);
        
        [A9,  dA9]  = arTriangle3D(yc,PA,PB,PABFE,doDerivative);
        [A10, dA10] = arTriangle3D(yc,PB,PF,PABFE,doDerivative);
        [A11, dA11] = arTriangle3D(yc,PF,PE,PABFE,doDerivative);
        [A12, dA12] = arTriangle3D(yc,PE,PA,PABFE,doDerivative);
        
        [A13, dA13] = arTriangle3D(yc,PC,PA,PCAEG,doDerivative);
        [A14, dA14] = arTriangle3D(yc,PA,PE,PCAEG,doDerivative);
        [A15, dA15] = arTriangle3D(yc,PE,PG,PCAEG,doDerivative);
        [A16, dA16] = arTriangle3D(yc,PG,PC,PCAEG,doDerivative);
        
        [A17, dA17] = arTriangle3D(yc,PD,PC,PDCGH,doDerivative);
        [A18, dA18] = arTriangle3D(yc,PC,PG,PDCGH,doDerivative);
        [A19, dA19] = arTriangle3D(yc,PG,PH,PDCGH,doDerivative);
        [A20, dA20] = arTriangle3D(yc,PH,PD,PDCGH,doDerivative);
        
        [A21, dA21] = arTriangle3D(yc,PE,PF,PEFHG,doDerivative);
        [A22, dA22] = arTriangle3D(yc,PF,PH,PEFHG,doDerivative);
        [A23, dA23] = arTriangle3D(yc,PH,PG,PEFHG,doDerivative);
        [A24, dA24] = arTriangle3D(yc,PG,PE,PEFHG,doDerivative);
        
        G = [A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12;A13;A14;A15;A16;A17;A18;A19;A20;A21;A22;A23;A24];
        if doDerivative,
            dG = [dA1;dA2;dA3;dA4;dA5;dA6;dA7;dA8;dA9;dA10;dA11;dA12;dA13;dA14;...
                  dA15;dA16;dA17;dA18;dA19;dA20;dA21;dA22;dA23;dA24];
        end
    case 'A-mf-3'
        G = geometrymexC(yc(:),m,'A');
   otherwise
        error('Error in file geometry.m : flag '' %s '' nyi!',flag);
end

%
%---------- matrix based helper --------------------
%

function [vol dVol] = volTetra3D(yc,P1,P2,PF,PM,doDerivative)
% computes volume of tetrahedra with edges P1 and P2, face-staggered point
% PF and cell-center PM
%
% volume is given by the determinant of the matrix
%
%     |       |       |       |
% M = | P1-PM | P2-PM | PF-PM |
%     |       |       |       |
%       col1  col2  col3
%
dVol = [];
yc = reshape(yc,[],3);
P1 = P1-PM;
P2 = P2-PM;
P3 = PF-PM;
col1 = P1*yc;
col2 = P2*yc;
col3 = P3*yc;

vol  = (1/6)*(col1(:,1) .* col2(:,2) .* col3(:,3) ...
    + col2(:,1) .* col3(:,2) .* col1(:,3) ...
    + col3(:,1) .* col1(:,2) .* col2(:,3) ...
    - col1(:,3) .* col2(:,2) .* col3(:,1) ...
    - col2(:,3) .* col3(:,2) .* col1(:,1) ...
    - col3(:,3) .* col1(:,2) .* col2(:,1));
if ~doDerivative,
    return;
end

dVol = (1/6)*[...
    sdiag(col2(:,2) .* col3(:,3) - col2(:,3) .* col3(:,2)) * P1 ...
    + sdiag(col3(:,2) .* col1(:,3) - col3(:,3) .* col1(:,2)) * P2 ...
    + sdiag(col1(:,2) .* col2(:,3) - col1(:,3) .* col2(:,2)) * P3,...
    sdiag(col1(:,1) .* col3(:,3) - col1(:,3) .* col3(:,1)) * P2...
    + sdiag(col2(:,1) .* col1(:,3) - col2(:,3) .* col1(:,1)) * P3...
    + sdiag(col3(:,1) .* col2(:,3) - col3(:,3) .* col2(:,1)) * P1,...
    + sdiag(col1(:,1) .* col2(:,2) - col1(:,2) .* col2(:,1)) * P3...
    + sdiag(col2(:,1) .* col3(:,2) - col2(:,2) .* col3(:,1)) * P1...
    + sdiag(col3(:,1) .* col1(:,2) - col3(:,2) .* col1(:,1)) * P2...
    ];

function [V,dV] =  arTriangle2D(P1,P2,PM,yc,doDerivative)
% computes area of triangle with P1 and P2 and cell-center PM
%
% volume is given by the determinant of the matrix
% 
%     |       |       |
% M = | P1-PM | P2-PM | 
%     |       |       |
%       col1    col2 

dV= [];
yc = reshape(yc,[],2);
V = ((P1*yc(:,1)-PM*yc(:,1)).*(P2*yc(:,2)-PM*yc(:,2)) ...
    - (P2*yc(:,1)-PM*yc(:,1)).*(P1*yc(:,2)-PM*yc(:,2)))/2;
if doDerivative,
    dV =[
        sdiag(P2*yc(:,2)-PM*yc(:,2))*(P1-PM) ...
        - sdiag(P1*yc(:,2)-PM*yc(:,2))*(P2-PM),...
        sdiag(P1*yc(:,1)-PM*yc(:,1))*(P2-PM) ...
        - sdiag(P2*yc(:,1)-PM*yc(:,1))*(P1-PM)
        ]/2;
end
function [ar dAr] = arTriangle3D(yc,P1,P2,PF,doDerivative)
% computes (squared) volume of triangle with edges P1 and P2 and the
% face-staggered point PF
%
% area is given by the euclidean length of the cross product
%
% ( B-A ) x ( C-A )
% =: v      =: w
%
dAr = [];
yc = reshape(yc,[],3);
Pv = (P1-PF);
Pw = (P2-PF);
v = Pv*yc;
w = Pw*yc;

vec = [(v(:,2) .* w(:,3) - v(:,3) .* w(:,2)), ...
    (v(:,1) .* w(:,3) - v(:,3) .* w(:,1)), ...
    (v(:,1) .* w(:,2) - v(:,2) .* w(:,1))];
ar = ( vec(:,1) .* vec(:,1) + vec(:,2) .* vec(:,2) + vec(:,3) .* vec(:,3));


if ~doDerivative,
    return;
end

dAr = 2* [
    sdiag( (v(:,1) .* w(:,3) - v(:,3) .* w(:,1)).* w(:,3) ...
    +(v(:,1) .* w(:,2) - v(:,2) .* w(:,1)).* w(:,2)) * Pv ...
    - sdiag((v(:,1) .* w(:,3) - v(:,3) .* w(:,1)).* v(:,3) ...
    +(v(:,1) .* w(:,2) - v(:,2) .* w(:,1)).* v(:,2))* Pw,...
    sdiag( (v(:,2) .* w(:,3) - v(:,3) .* w(:,2)).* w(:,3) ...
    -(v(:,1) .* w(:,2) - v(:,2) .* w(:,1)).* w(:,1)) * Pv ...
    - sdiag((v(:,2) .* w(:,3) - v(:,3) .* w(:,2)).* v(:,3)...
    -((v(:,1) .* w(:,2) - v(:,2) .* w(:,1)).* v(:,1)))* Pw,...
    - sdiag((v(:,2) .* w(:,3) - v(:,3) .* w(:,2)).* w(:,2) ...
    +(v(:,1) .* w(:,3) - v(:,3) .* w(:,1)).* w(:,1) )* Pv...
    + sdiag( (v(:,2) .* w(:,3) - v(:,3) .* w(:,2)).* v(:,2)...
    +(v(:,1) .* w(:,3) - v(:,3) .* w(:,1)).* v(:,1)) * Pw ...
    ];



function P = Pedge3D(edge,m)
% function P = Pedge3D(edge,m)
%
% projector for edge points of 3D nodal data
%
% INPUT
% edge - String specifying the corner
% m
%
persistent S;
persistent PA PB PC PD PE PF PG PH;
if isempty(S) || size(S,1)~=prod(m+1),
    S  = speye(prod(m+1),prod(m+1));
end
ox = ones(m(1)+1,1);
oy = ones(m(2)+1,1);
oz = ones(m(3)+1,1);

switch edge
    case 'A'
        if isempty(PA) || size(PA,2)~=prod(m+1),
            vecx = ox; vecx(end) = 0;
            vecy = oy; vecy(end) = 0;
            vecz = oz; vecz(end) = 0;
            PA = S(kron(kron(vecz,vecy),vecx)==1,:);
        end
        P = PA;
    case 'B'
        if isempty(PB) || size(PB,2)~=prod(m+1),
            vecx = ox; vecx(1) = 0;
            vecy = oy; vecy(end) = 0;
            vecz = oz; vecz(end) = 0;
            PB = S(kron(kron(vecz,vecy),vecx)==1,:);
        end
        P = PB;
    case 'C'
        if isempty(PC) || size(PC,2)~=prod(m+1),
            vecx = ox; vecx(end) = 0;
            vecy = oy; vecy(1) = 0;
            vecz = oz; vecz(end) = 0;
            PC = S(kron(kron(vecz,vecy),vecx)==1,:);
        end
        P = PC;
    case 'D'
        if isempty(PD) || size(PD,2)~=prod(m+1),
            vecx = ox; vecx(1) = 0;
            vecy = oy; vecy(1) = 0;
            vecz = oz; vecz(end) = 0;
            PD = S(kron(kron(vecz,vecy),vecx)==1,:);
        end
        P = PD;
    case 'E'
        if isempty(PE) || size(PE,2)~=prod(m+1),
            vecx = ox; vecx(end) = 0;
            vecy = oy; vecy(end) = 0;
            vecz = oz; vecz(1) = 0;
            PE = S(kron(kron(vecz,vecy),vecx)==1,:);
        end
        P = PE;
    case 'F'
        if isempty(PF) || size(PF,2)~=prod(m+1),
            vecx = ox; vecx(1) = 0;
            vecy = oy; vecy(end) = 0;
            vecz = oz; vecz(1) = 0;
            PF = S(kron(kron(vecz,vecy),vecx)==1,:);
        end
        P = PF;
    case 'G'
        if isempty(PG) || size(PG,2)~=prod(m+1),
            vecx = ox; vecx(end) = 0;
            vecy = oy; vecy(1) = 0;
            vecz = oz; vecz(1) = 0;
            PG = S(kron(kron(vecz,vecy),vecx)==1,:);
        end
        P = PG;
    case 'H'
        if isempty(PH) || size(PH,2)~=prod(m+1),
            vecx = ox; vecx(1) = 0;
            vecy = oy; vecy(1) = 0;
            vecz = oz; vecz(1) = 0;
            PH = S(kron(kron(vecz,vecy),vecx)==1,:);
        end
        P = PH;
    case 'ABDC'
        P = (Pedge3D('A',m) + Pedge3D('B',m) + Pedge3D('D',m) + Pedge3D('C',m)) /4;
    case 'BDHF'
        P = (Pedge3D('B',m) + Pedge3D('D',m) + Pedge3D('H',m) + Pedge3D('F',m)) /4;
    case 'ABFE'
        P = (Pedge3D('A',m) + Pedge3D('B',m) + Pedge3D('F',m) + Pedge3D('E',m)) /4;
    case 'CAEG'
        P = (Pedge3D('C',m) + Pedge3D('A',m) + Pedge3D('E',m) + Pedge3D('G',m)) /4;
    case 'DCGH'
        P = (Pedge3D('D',m) + Pedge3D('C',m) + Pedge3D('G',m) + Pedge3D('H',m)) /4;
    case 'EFHG'
        P = (Pedge3D('E',m) + Pedge3D('F',m) + Pedge3D('H',m) + Pedge3D('G',m)) /4;
    case 'M'
        P = (Pedge3D('A',m) + Pedge3D('B',m) + Pedge3D('D',m) + Pedge3D('C',m) + ...
            Pedge3D('E',m) + Pedge3D('F',m) + Pedge3D('H',m) + Pedge3D('G',m)) /8;
    otherwise
        error('Invalid Edge');
        
end
% shortcut for sparse diagonal matrices
function a = sdiag(a)
a = spdiags(reshape(a,[],1),0,length(a),length(a));

function testThisMfile
clc;
help(mfilename)

fprintf('Testing 2D \n');
omega = [0 2 1 3];
m = [4 9];
tol = 1e-14;
xc = getNodalGrid(omega,m);
mbFctn = @(y,m,flag,x) geometry(y,m,flag,'x',x,'matrixFree',0);
mfFctn = @(y,m,flag,x) geometry(y,m,flag,'x',x,'matrixFree',1);
fprintf('Test V \t\t ');
[mbV, mbdV] = mbFctn(xc,m,'V',[]);
[mfV, mfdV] = mfFctn(xc,m,'V',[]);
RE = norm(mbV-mfV)/norm(mbV);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
fprintf('Test VRange \t ');
[mbVRange] = mbFctn(xc,m,'VRange',[]);
[mfVRange] = mfFctn(xc,m,'VRange',[]);
RE = norm(mbVRange-mfVRange)/norm(mbVRange);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
x = randn(2*prod(m+1),1);
fprintf('Test Jac \t ');
[mbJac, mbdJac] = geometry(xc,m,'Jac','x',x,'matrixFree',0,'omega',omega);
[mfJac, mfdJac] = geometry(xc,m,'Jac','x',x,'matrixFree',1,'omega',omega);
RE =norm(mbJac-mfJac)/norm(mbJac);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
fprintf('Test dJac \t ');
x = 1e2*randn(2*prod(m+1),1);
RE = norm(mbdJac*x-mfdJac.dJac(xc,m,x)) / norm(mbdJac*x);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
fprintf('Test dJacadj \t ');
x = 1e2*randn(prod(m),1);
RE = norm(mbdJac'*x-mfdJac.dJacadj(xc,m,x)) /norm(mbdJac'*x);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);

fprintf('\nTesting 3D \n');
omega = [0 2 1 3 4 6];
m = [2 3 5];
xc = getNodalGrid(omega,m);
xc = xc + 1e-2 * randn(size(xc));
mbFctn = @(y,m,flag,x) geometry(y,m,flag,'x',x,'matrixFree',0);
mfFctn = @(y,m,flag,x) geometry(y,m,flag,'x',x,'matrixFree',1);
fprintf('Test V \t\t ');
[mbV, mbdV] = mbFctn(xc,m,'V',[]);
[mfV, mfdV] = mfFctn(xc,m,'V',[]);
RE = norm(mbV-mfV)/norm(mbV);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
fprintf('Test VRange \t ');
[mbVRange] = mbFctn(xc,m,'VRange',[]);
[mfVRange] = mfFctn(xc,m,'VRange',[]);
RE = norm(mbVRange-mfVRange)/norm(mbVRange);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
fprintf('Test Jac \t ');
[mbJac, mbdJac] = geometry(xc,m,'Jac','x',x,'matrixFree',0,'omega',omega);
[mfJac, mfdJac] = geometry(xc,m,'Jac','x',x,'matrixFree',1,'omega',omega);
RE = norm(mbJac-mfJac)/norm(mbJac);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
fprintf('Test dJac \t ');
x = randn(3*prod(m+1),1);
RE = norm(mbdJac*x-mfdJac.dJac(xc,m,x)) / norm(mbdJac*x);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
fprintf('Test dJacadj \t ');
x = randn(prod(m),1);
RE = norm(mbdJac'*x-mfdJac.dJacadj(xc,m,x)) /norm(mbdJac'*x);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);
fprintf('Test A \t\t ');
[mbA] = geometry(xc,m,'A','matrixFree',0);
[mfA] = geometry(xc,m,'A','matrixFree',1);
RE = norm(mbA-mfA) /norm(mbA);
fprintf(' error: %1.2e \t passed ? %d \n',RE,RE<tol);


%{ 
	=======================================================================================
	FAIR: Flexible Algorithms for Image Registration, Version 2011
	Copyright (c): Jan Modersitzki
	Maria-Goeppert-Str. 1a, D-23562 Luebeck, Germany
	Email: jan.modersitzki@mic.uni-luebeck.de
	URL:   http://www.mic.uni-luebeck.de/people/jan-modersitzki.html
	=======================================================================================
	No part of this code may be reproduced, stored in a retrieval system,
	translated, transcribed, transmitted, or distributed in any form
	or by any means, means, manual, electric, electronic, electro-magnetic,
	mechanical, chemical, optical, photocopying, recording, or otherwise,
	without the prior explicit written permission of the authors or their
	designated proxies. In no event shall the above copyright notice be
	removed or altered in any way.

	This code is provided "as is", without any warranty of any kind, either
	expressed or implied, including but not limited to, any implied warranty
	of merchantibility or fitness for any purpose. In no event will any party
	who distributed the code be liable for damages or for any claim(s) by
	any other party, including but not limited to, any lost profits, lost
	monies, lost data or data rendered inaccurate, losses sustained by
	third parties, or any other special, incidental or consequential damages
	arrising out of the use or inability to use the program, even if the
	possibility of such damages has been advised against. The entire risk
	as to the quality, the performace, and the fitness of the program for any
	particular purpose lies with the party using the code.
	=======================================================================================
	Any use of this code constitutes acceptance of the terms of the above statements
	=======================================================================================
%}