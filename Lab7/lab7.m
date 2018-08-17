function varargout = lab7(varargin)
% lab7 MATLAB code for lab7.fig
%      lab7, by itself, creates a new lab7 or raises the existing
%      singleton*.
%
%      H = lab7 returns the handle to a new lab7 or the handle to
%      the existing singleton*.
%
%      lab7('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in lab7.M with the given input arguments.
%
%      lab7('Property','Value',...) creates a new lab7 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lab7_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lab7_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lab7

% Last Modified by GUIDE v2.5 01-May-2018 15:44:38

% Begin initialization code - DO NOT EDIT


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lab7_OpeningFcn, ...
                   'gui_OutputFcn',  @lab7_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



function [m,n,newE] = grValidation(E);
% The validation of array E - auxiliary function for GrTheory Toolbox.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if ~isnumeric(E),
  error('The array E must be numeric!') 
end
if ~isreal(E),
  error('The array E must be real!') 
end
se=size(E); % size of array E
if length(se)~=2,
  error('The array E must be 2D!') 
end
if (se(2)<2),
  error('The array E must have 2 or 3 columns!'), 
end
if ~all(all(E(:,1:2)>0)),
  error('1st and 2nd columns of the array E must be positive!')
end
if ~all(all((E(:,1:2)==round(E(:,1:2))))),
  error('1st and 2nd columns of the array E must be integer!')
end
m=se(1);
if se(2)<3, % not set the weight
  E(:,3)=1; % all weights =1
end
newE=E(:,1:3);
n=max(max(newE(:,1:2))); % number of vertexes
return

function [pTS,fmin]=grTravSale(C)
% Function [pTS,fmin]=grTravSale(C) solve the nonsymmetrical
% traveling salesman problem.
% Input parameter: 
%   C(n,n) - matrix of distances between cities, 
%     maybe, nonsymmetrical;
%     n - number of cities.
% Output parameters: 
%   pTS(n) - the order of cities;
%   fmin - length of way.
% Uses the reduction to integer LP-problem:
% Look: Miller C.E., Tucker A. W., Zemlin R. A. 
% Integer Programming Formulation of Traveling Salesman Problems. 
% J.ACM, 1960, Vol.7, p. 326-329.
% Needed other products: MIQP.M.
% This software may be free downloaded from site:
% http://control.ee.ethz.ch/~hybrid/miqp/
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
if ~isnumeric(C),
  error('The array C must be numeric!') 
end
if ~isreal(C),
  error('The array C must be real!') 
end
s=size(C); % size of array C
if length(s)~=2,
  error('The array C must be 2D!') 
end
if s(1)~=s(2),
  error('Matrix C must be square!')
end
if s(1)<3,
  error('Must be not less than 3 cities!')
end

% ============ Size of problem ====================
n=s(1); % number of vertexes
m=n*(n-1); % number of arrows

% ============ Parameters of integer LP problem ========
Aeq=[]; % for the matrix of the boundary equations
for k1=1:n,
  z1=zeros(n);
  z1(k1,:)=1;
  z2=[z1;eye(n)];
  Aeq=[Aeq z2([1:2*n-1],setdiff([1:n],k1))];
end
Aeq=[Aeq zeros(2*n-1,n-1)];
A=[]; % for the matrix of the boundary inequations
for k1=2:n,
  z1=[];
  for k2=1:n,
    z2=eye(n)*(n-1)*(k2==k1);
    z1=[z1 z2(setdiff([2:n],k1),setdiff([1:n],k2))];
  end
  z2=-eye(n);
  z2(:,k1)=z2(:,k1)+1;
  A=[A;[z1 z2(setdiff([2:n],k1),2:n)]];
end
beq=ones(2*n-1,1); % the right parts of the boundary equations
b=ones((n-1)*(n-2),1)*(n-2); % the right parts of the boundary inequations
C1=C'+diag(ones(1,n)*NaN);
C2=C1(:);
c=[C2(~isnan(C2));zeros(n-1,1)]; % the factors for objective function
vlb=[zeros(m,1);-inf*ones(n-1,1)]; % the lower bounds
vub=[ones(m,1);inf*ones(n-1,1)]; % the upper bounds
H=zeros(n^2-1); % Hessian

% ============= We solve the MIQP problem ==========
[xmin,fmin]=MIQP(H,c,A,b,Aeq,beq,[1:m],vlb,vub);

% ============= We return the results ==============
eik=round(xmin(1:m)); % the arrows of the way
e1=[zeros(1,n-1);reshape(eik,n,n-1)];
e2=[e1(:);0]; % we add zero to a diagonal
e3=(reshape(e2,n,n))'; % the matrix of the way
pTS=[1 find(e3(1,:))]; % we build the way
while pTS(end)>1, % we add the city to the way
  pTS=[pTS find(e3(pTS(end),:))];
end
return

function Etc=grTranClos(E)
% Function Etc=grTranClos(E) built 
% the transitive closure for the digraph E.
% Input parameter: 
%   E(m,2) - the arrows of digraph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of arrows.
% Output parameter:
%   Etc(mtc,2) - the arrows of digraph with trasitive closure for E;
%     mtc - number of arrows.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % data validation
% ============ Atc=A+A^2+A^3+A^4+... ==========================
A1=zeros(n);
A1((E(:,2)-1)*n+E(:,1))=1; % adjacency matrix A
A0=A1;
An=(A1*A1)>0; % A^2
A=(A1+An)>0; % A+A^2
while ~all(all(A==A1)),
  A1=A;
  An=((double(An))*A0)>0;
  A=(A+An)>0;
end
[r,c]=find(A);
Etc=[r c]; % all arrows of the transitive closure
return

function [dMWP,ssp]=grShortVerPath(E,Wv)
% Function dMVP=grShortVerPath(E,Wv) for digraph with weighted vertexes 
% solve the problem about the path with minimal weight of verticies.
% Input parameters: 
%   E(m,2) - the arrows of digraph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of arrows.
%   Wv(n,1) - the weights of verticies; if this parameter omitted, 
%             then all vertices have weight 1.
%     n - number of verticies.
% Output parameter:
%   dMWP(k) - the vector with number of verticies included
%     to path with minimal weight.
% [dMWP,ssp]=grMaxVerPath(E,Wv) return also
%   ssp - sum of vertexes weights on this path.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru
% Acknowledgements to Dr. Albert Niel (Austria)
% for testing of this algorithm.

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation
if nargin<2,
  Wv = ones(n,1);
end
if n~=length(Wv),
  error('Length of Wv not equal to n!')
end
%============ Duplication of vertices ==================
E1 = [[E(:,1)+n, E(:,2), zeros(m,1)]; [(1:n)', (n+1:2*n)', Wv(:)]];
%============ Finding of the base and contrabase =======
BG = grBase(E1); % base of digraph
CBG = grCoBase(E1); % contrabase of digraph
S1 = unique(BG(:))'; % all vertices-sources
T1 = unique(CBG(:))'; % all vertices-tails
ns = length(S1); % number of vertices in base
nt = length(T1); % number of vertices in contrabase
%============ Add source and tail vertexes =============
s = 2*n+1; % source vertex
t = 2*n+2; % tail vertex
E2=[[ones(1,ns)*s;S1;zeros(1,ns)]';E1;[T1;ones(1,nt)*t;zeros(1,nt)]'];
%============ Shortest Path Problem ====================
[dSP,sp] = grShortPath(E2,s,t);
dMWP = sp(2:2:end-1);
ssp = sum(Wv(dMWP));
return

function [dSP,sp]=grShortPath(E,s,t)
% Function dSP=grShortPath(E) solve the problem about
% the shortest path between any vertexes of digraph.
% Input parameter: 
%   E(m,2) or (m,3) - the arrows of digraph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of arrow;
%     m - number of arrows.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   dSP(n,n) - the matrix of shortest path.
%   Each element dSP(i,j) is the shortest path 
%   from vertex i to vertex j (may be inf,
%   if vertex j is not accessible from vertex i).
% Uses the algorithm of R.W.Floyd, S.A.Warshall.
% [dSP,sp]=grShortPath(E,s,t) - find also
% the shortest paths from vertex s (source) to vertex t (tail).
% In this case output parameter sp is vector with numbers 
% of vertexes, included to shortest path from s to t.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru
% Acknowledgements to Prof. Gerard Biau (France)
% for testing of this algorithm.

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

% ================ Initial values ===============
dSP=ones(n)*inf; % initial distances
dSP((E(:,2)-1)*n+E(:,1))=E(:,3);
%dSP0=dSP;
% ========= The main cycle of Floyd-Warshall algorithm =========
for j=1:n,
  i=setdiff((1:n),j);
  dSP(i,i)=min(dSP(i,i),repmat(dSP(i,j),1,n-1)+repmat(dSP(j,i),n-1,1));
end
sp=[];
if (nargin<3)|(isempty(s))|(isempty(t)),
  return
end
s=s(1);
t=t(1);
if (~(s==round(s)))|(~(t==round(t)))|(s<1)|(s>n)|(t<1)|(t>n),
  error(['s and t must be integer from 1 to ' num2str(n)])
end
if isinf(dSP(s,t)), % t is not accessible from s
  return
end
dSP1=dSP;
dSP1(1:n+1:n^2)=0; % modified dSP
l=ones(m,1); % label for each arrow
sp=t; % final vertex
while ~(sp(1)==s),
  nv=find((E(:,2)==sp(1))&l); % all labeled arrows to sp(1)
  vnv=abs((dSP1(s,sp(1))-dSP1(s,E(nv,1)))'-E(nv,3))<eps*1e3; % valided arrows
  l(nv(~vnv))=0; % labels of not valided arrows
  if all(~vnv), % invalided arrows
    l(find((E(:,1)==sp(1))&(E(:,2)==sp(2))))=0; 
    sp=sp(2:end); % one step back
  else
    nv=nv(vnv); % rested valided arrows
    sp=[E(nv(1),1) sp]; % add one vertex to shortest path
  end
end
return

function h=grPlot(V,E,kind,vkind,ekind,sa)
% Function h=grPlot(V,E,kind,vkind,ekind,sa) 
% draw the plot of the graph (digraph).
% Input parameters: 
%   V(n,2) or (n,3) - the coordinates of vertexes
%     (1st column - x, 2nd - y) and, maybe, 3rd - the weights;
%     n - number of vertexes.
%     If V(n,2), we write labels: numbers of vertexes,
%     if V(n,3), we write labels: the weights of vertexes.
%     If V=[], use regular n-angle.
%   E(m,2) or (m,3) - the edges of graph (arrows of digraph)
%     and their weight; 1st and 2nd elements of each row 
%     is numbers of vertexes;
%     3rd elements of each row is weight of arrow;
%     m - number of arrows.
%     If E(m,2), we write labels: numbers of edges (arrows);
%     if E(m,3), we write labels: weights of edges (arrows).
%     For disconnected graph use E=[] or h=PlotGraph(V).
%   kind - the kind of graph.
%   kind = 'g' (to draw the graph) or 'd' (to draw digraph);
%   (optional, 'g' default).
%   vkind - kind of labels for vertexes (optional).
%   ekind - kind of labels for edges or arrows (optional);
%   sa - size of arrows (optional, default 1).
%   For vkind and ekind use the format of function FPRINTF,
%   for example, '%8.3f', '%14.10f' etc. Default value is '%d'.
%   Use '' (empty string) for don't draw labels.
% Output parameter:
%   h - handle of figure (optional).
% See also GPLOT.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru
% Acknowledgements to Mr.Howard (howardz@cc.gatech.edu)
% for testing of this algorithm.

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
if (nargin==1) & isempty(V),
  error('V is empty and E is not determined!')
end
if (nargin==2) & isempty(V) & isempty(E),
  error('V and E are empty!')
end
if ~isempty(V),
  if ~isnumeric(V),
    error('The array V must be numeric!') 
  end
  sv=size(V); % size of array V
  if length(sv)~=2,
    error('The array V must be 2D!') 
  end
  if (sv(2)<2),
    error('The array V must have 2 or 3 columns!'), 
  end
  if nargin==1, % disconnected graph
    E=[];
  end
end
if ~isempty(E), % for connected graph
  [m,n,newE]=grValidation(E);
  we=min(3,size(E,2)); % 3 for weighed edges
  E=newE;
  if isempty(V), % regular n-angle
    V=[cos(2*pi*[1:n]'/n) sin(2*pi*[1:n]'/n)];
    sv=size(V); % size of array V
  end
  if n>sv(1),
    error('Several vertexes is not determined!');
  end
else
  m=0;
end

% ============= Other arguments ==================
n=sv(1); % real number of vertexes
wv=min(3,sv(2)); % 3 for weighted vertexes
if nargin<3, % only 2 input parameters
  kind1='g';
else
  if isempty(kind),
    kind='g';
    kind1='g';
  end
  if ~ischar(kind),
    error('The argument kind must be a string!')
  else
    kind1=lower(kind(1));
  end
end
if nargin<4,
  vkind1='%d';
else
  if ~ischar(vkind),
    error('The argument vkind must be a string!')
  else
    vkind1=lower(vkind);
  end
end
if nargin<5,
  ekind1='%d';
else
  if ~ischar(ekind),
    error('The argument ekind must be a string!')
  else
    ekind1=lower(ekind);
  end
end
if nargin<6,
  sa=1;
end
md=inf; % the minimal distance between vertexes
for k1=1:n-1,
  for k2=k1+1:n,
    md=min(md,sum((V(k1,:)-V(k2,:)).^2)^0.5);
  end
end
if md<eps, % identical vertexes
  error('The array V have identical rows!')
else
  V(:,1:2)=V(:,1:2)/md; % normalization
end
r=0.1; % for multiple edges
tr=linspace(pi/4,3*pi/4);
xr=0.5-cos(tr)/2^0.5;
yr=(sin(tr)-2^0.5/2)/(1-2^0.5/2);
t=linspace(-pi/2,3*pi/2); % for loops
xc=0.1*cos(t);
yc=0.1*sin(t);

% we sort the edges
if ~isempty(E),
  E=[zeros(m,1),[1:m]',E]; % 1st column for change, 2nd column is edge number
  need2=find(E(:,4)<E(:,3)); % for replace v1<->v2
  tmp=E(need2,3);
  E(need2,3)=E(need2,4);
  E(need2,4)=tmp;
  E(need2,1)=1; % 1, if v1<->v2
  [e1,ie1]=sort(E(:,3)); % sort by 1st vertex
  E1=E(ie1,:);
  for k2=E1(1,3):E1(end,3),
    num2=find(E1(:,3)==k2);
    if ~isempty(num2), % sort by 2nd vertex
      E3=E1(num2,:);
      [e3,ie3]=sort(E3(:,4));
      E4=E3(ie3,:);
      E1(num2,:)=E4;
    end
  end
  ip=find(E1(:,3)==E1(:,4)); % we find loops
  Ep=E1(ip,:); % loops
  E2=E1(setdiff([1:m],ip),:); % edges without loops
end

% we paint the graph
hh=figure;
hold on
plot(V(:,1),V(:,2),'k.','MarkerSize',20)
axis equal
h1=get(gca); % handle of current figure
if ~isempty(vkind1), % labels of vertexes
  for k=1:n,
    if wv==3,
      s=sprintf(vkind1,V(k,3));
    else
      s=sprintf(vkind1,k);
    end
    text(V(k,1)+0.05,V(k,2)-0.07,s);
  end
end

% edges (arrows)
if ~isempty(E),
  k=0;
  m2=size(E2,1); % number of edges without loops
  while k<m2,
    k=k+1; % current edge
    MyE=V(E2(k,3:4),1:2); % numbers of vertexes 1, 2
    k1=1; % we find the multiple edges
    if k<m2,
      while all(E2(k,3:4)==E2(k+k1,3:4)),
        k1=k1+1;
        if k+k1>m2,
          break;
        end
      end
    end
    ry=r*[1:k1];
    ry=ry-mean(ry); % radius
    l=norm(MyE(1,:)-MyE(2,:)); % lenght of line
    dx=MyE(2,1)-MyE(1,1);
    dy=MyE(2,2)-MyE(1,2);
    alpha=atan2(dy,dx); % angle of rotation
    cosa=cos(alpha);
    sina=sin(alpha);
    MyX=xr*l;
    for k2=1:k1, % we draw the edges (arrows)
      MyY=yr*ry(k2);
      MyXg=MyX*cosa-MyY*sina+MyE(1,1);
      MyYg=MyX*sina+MyY*cosa+MyE(1,2);
     if E2(k+k2-1,5)<0
         plot(MyXg,MyYg,'k-','Color','blue');
     else
      plot(MyXg,MyYg,'k-','Color','red');
     end
      if kind1=='d', % digraph with arrows
        if E2(k+k2-1,1)==1,
          [xa,ya]=CreateArrow(MyXg(1:2),MyYg(1:2),sa);
          fill(xa,ya,'k');
        else
          [xa,ya]=CreateArrow(MyXg(end:-1:end-1),MyYg(end:-1:end-1),sa);
          fill(xa,ya,'k');
        end
      end
      if ~isempty(ekind1), % labels of edges (arrows)
        if we==3,
          s=sprintf(ekind1,E2(k+k2-1,5));
        else
          s=sprintf(ekind1,E2(k+k2-1,2));
        end
        text(MyXg(length(MyXg)/2),MyYg(length(MyYg)/2),s);
      end
    end
    k=k+k1-1;
  end
  % we draw the loops
  k=0;
  ml=size(Ep,1); % number of loops
  while k<ml,
    k=k+1; % current loop
    MyV=V(Ep(k,3),1:2); % vertexes
    k1=1; % we find the multiple loops
    if k<ml,
      while all(Ep(k,3:4)==Ep(k+k1,3:4)),
        k1=k1+1;
        if k+k1>ml,
          break;
        end
      end
    end
    ry=[1:k1]+1; % radius
    for k2=1:k1, % we draw the loop
      MyX=xc*ry(k2)+MyV(1);
      MyY=(yc+r)*ry(k2)+MyV(2);
      plot(MyX,MyY,'k-');
      if kind1=='d',
        [xa,ya]=CreateArrow(MyX([1 10]),MyY([1 10]),sa);
        fill(xa,ya,'k');
      end
      if ~isempty(ekind1), % labels of edges (arrows)
        if we==3,
          s=sprintf(ekind1,Ep(k+k2-1,5));
        else
          s=sprintf(ekind1,Ep(k+k2-1,2));
        end
        text(MyX(length(MyX)/2),MyY(length(MyY)/2),s);
      end
    end
    k=k+k1-1;
  end
end
hold off
axis off
if nargout==1,
  h=hh;
end
return

function [xa,ya]=CreateArrow(x,y,sa)
% create arrow with length 0.1*sa with tip x(1), y(1) 
% and direction from x(2), y(2)

xa1=sa*[0 0.1 0.08 0.1 0]';
ya1=sa*[0 0.03 0 -0.03 0]';
dx=diff(x);
dy=diff(y);
alpha=atan2(dy,dx); % angle of rotation
cosa=cos(alpha);
sina=sin(alpha);
xa=xa1*cosa-ya1*sina+x(1);
ya=xa1*sina+ya1*cosa+y(1);
return

function [CrP,Ts,Td]=grPERT(E)
% Function [CrP,Ts,Td]=grPERT(E) solve
% the project evaluation research task.
% Input parameter: 
%   E(m,2) or (m,3) - the arrows of digraph 
%     and their weight (working time);
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of arrow;
%     m - number of arrows.
%     If we set the array E(m,2), then all weights is 1.
% The digraph E must be acyclic.
% Output parameters:
%   CrP - the critical path (vector-row with numbers of vertexes);
%   Ts(1,n) - the start times for each vertex (event);
%   max(Ts) is length of critical path;
%   Td(m,1) - the delay times for each arrow (work).
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

[Dec,Ord]=grDecOrd(E); % the decomposition and partial ordering
n=max(max(E(:,1:2))); % number of vertexes (events)
m=size(E,1); % number of arrows (works)
if size(Dec,2)<n,
  error('The digraph is cyclic!')
end
s=find(Dec(:,1)); % the first event
t=find(Dec(:,end)); % the last event
Em=[E(:,1:2) -E(:,3)]; % weight with minus
[dSP,CrP]=grShortPath(Em,s,t); % the shortest path
Ts=-dSP(s,:); % the start time for each vertex (event)
Ts(find(isinf(Ts)))=0; % change inf to 0
Td=(Ts(E(:,2))-Ts(E(:,1)))'-E(:,3); % the time of delay for each work
return


function nMC=grMinVerCover(E,d)
% Function nMC=grMinVerCover(E,d) solve the minimal vertex cover problem.
% Input parameters: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
%   d(n) (optional) - the weights of vertexes,
%     n - number of vertexes.
%     If we have only 1st parameter E, then all d=1.
% Output parameter:
%   nMC - the list of the numbers of vertexes included 
%     in the minimal (weighted) vertex cover.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation
if nargin<2, % we may only 1st parameter
  d=ones(n,1); % all weights =1
else
  d=d(:); % reshape to vector-column
  if length(d)<n, % the poor length
    error('The length of the vector d is poor!')
  else
    n=length(d); % Number of Vertexes
  end
end

% ============= Parameters of integer LP problem ==========
A=zeros(n,m); % for incidence matrix
A(E(:,1:2)+repmat(([1:m]'-1)*n,1,2))=1; % we fill the incidence matrix
options=optimset('bintprog'); % the default options
options.Display='off'; % we change the output

% ============= We solve the MILP problem ==========
xmin=bintprog(d,-A',-ones(m,1),[],[],[],options);
nMC=find(round(xmin)); % the answer - numbers of vertexes
return

function nMST=grMinSpanTree(E)
% Function nMST=grMinSpanTree(E) solve 
% the minimal spanning tree problem for a connected graph.
% Input parameter: 
%   E(m,2) or (m,3) - the edges of graph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of edge;
%     m - number of edges.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   nMST(n-1,1) - the list of the numbers of edges included 
%     in the minimal (weighted) spanning tree in the including order.
% Uses the greedy algorithm.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

% ============= The data preparation ==================
En=[(1:m)',E]; % we add the numbers
En(:,2:3)=sort(En(:,2:3)')'; % edges on increase order
ln=find(En(:,2)==En(:,3)); % the loops numbers
En=En(setdiff([1:size(En,1)]',ln),:); % we delete the loops
[w,iw]=sort(En(:,4)); % sort by weight
Ens=En(iw,:); % sorted edges

% === We build the minimal spanning tree by the greedy algorithm ===
Emst=Ens(1,:); % 1st edge include to minimal spanning tree
Ens=Ens(2:end,:); % rested edges
while (size(Emst,1)<n-1)&(~isempty(Ens)),
  Emst=[Emst;Ens(1,:)]; % we add next edge to spanning tree
  Ens=Ens(2:end,:); % rested edges
  if any((Emst(end,2)==Emst(1:end-1,2))&...
         (Emst(end,3)==Emst(1:end-1,3))) | ...
     IsCycle(Emst(:,2:3)), % the multiple edge or cycle
    Emst=Emst(1:end-1,:); % we delete the last added edge
  end
end
nMST=Emst(:,1); % numbers of edges
return

function ic=IsCycle(E); % true, if graph E have cycle
n=max(max(E)); % number of vertexes
A=zeros(n);
A((E(:,1)-1)*n+E(:,2))=1;
A=A+A'; % the connectivity matrix
p=sum(A); % the vertexes power
ic=false;
while any(p<=1), % we delete all tails
  nc=find(p>1); % rested vertexes
  if isempty(nc),
    return
  end
  A=A(nc,nc); % new connectivity matrix
  p=sum(A); % new powers
end
ic=true;
return

function nMC=grMinEdgeCover(E)
% Function nMC=grMinEdgeCover(E) solve the minimal edge cover problem.
% Input parameter: 
%   E(m,2) or (m,3) - the edges of graph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of edge;
%     m - number of edges.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   nMC - the list of the numbers of edges included 
%     in the minimal (weighted) edge cover.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

% ============= Parameters of integer LP problem ==========
A=zeros(n,m); % for incidence matrix
A(E(:,1:2)+repmat(([1:m]'-1)*n,1,2))=1; % we fill the incidence matrix
options=optimset('bintprog'); % the default options
options.Display='off'; % we change the output

% ============= We solve the integer LP problem ==========
xmin=bintprog(E(:,3),-A,-ones(n,1),[],[],[],options);
nMC=find(round(xmin)); % the answer - numbers of edges
return

function [nMCS,mf]=grMinCutSet(E,s,t)
% Function [nMCS,mf]=grMinCutSet(E,s,t) find the first 
% minimal cut-sets of the network.
% Input parameters: 
%   E(m,2) or (m,3) - the arrows of digraph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of arrow;
%     m - number of arrows.
%     If we set the array E(m,2), then all weights is 1.
%   s - input (source) of the network (number of vertex);
%   t - output (sink) of the network (number of vertex).
% Output parameters: 
%   nMCS(ncs) - the list of the numbers of arrows included 
%     in first minimal cut-set; ncs - number of of arrows
%     of the minimal cut-sets.
%   mf - the total flow through each minimal cut-set.
% Uses the reduction to maximal flow problem.
% Required the Optimization Toolbox
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

[v,mf]=grMaxFlows(E,s,t); % the maximal flow
[m,n,E] = grValidation(E); % E data validation
ecs=find((abs(v)<1e-8)|(abs(E(:,3)-v)<1e-8)); % v=0 or v=E(:,3);
E1=E(setdiff([1:m]',ecs),1:2);
E1=[E1;fliplr(E1)]; % all arrows
d=grDecOrd(E1); % strongly connected components
nd1=n-size(d,1); % number of terminal vertexes
d(size(d,1)+1:size(d,1)+nd1,size(d,2)+1:size(d,2)+nd1)=eye(nd1);
nz=rem(find(d'),size(d,2)); % number of zone for each vertex
nz(find(nz==0))=size(d,2);
Ecsn=sort(nz(E(ecs,1:2))')'; % numbers of connected zones
nMCS=ecs(find((Ecsn(:,1)==nz(s))|(Ecsn(:,2)==nz(s)))); % first cut-set
return

function nMS=grMinAbsVerSet(E,d)
% Function nMS=grMinAbsVerSet(E,d) solve the minimal absorbant set problem
%   for the graph vertexes.
% Input parameters: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
%   d(n) (optional) - the weights of vertexes,
%     n - number of vertexes.
%     If we have only 1st parameter E, then all d=1.
% Output parameter:
%   nMS - the list of the numbers of vertexes included 
%     in the minimal (weighted) absorbant set of vertexes.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation
if nargin<2, % we may only 1st parameter
  d=ones(n,1); % all weights =1
else
  d=d(:); % reshape to vector-column
  if length(d)<n, % the poor length
    error('The length of the vector d is poor!')
  else
    n=length(d); % Number of Vertexes
  end
end

% ============= Parameters of integer LP problem ==========
A=eye(n);
A((E(:,2)-1)*n+E(:,1))=1; % adjacency matrix + main diagonal
A=double(A+A'>0); % symmetrical
options=optimset('bintprog'); % the default options
options.Display='off'; % we change the output

% ============= We solve the MILP problem ==========
xmin=bintprog(d,-A,-ones(n,1),[],[],[],options);
nMS=find(round(xmin)); % the answer - numbers of vertexes
return

function nMS=grMinAbsEdgeSet(E)
% Function nMS=grMinAbsEdgeSet(E) solve the minimal absorbant set problem
%   for the graph edges.
% Input parameter: 
%   E(m,2) or (m,3) - the edges of graph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of edge;
%     m - number of edges.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   nMS - the list of the numbers of edges included 
%     in the minimal (weighted) absorbant set of edges.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

% ============= Parameters of integer LP problem ==========
B=zeros(m);
for k=1:m,
  nn=find((E(:,1)==E(k,1))|(E(:,1)==E(k,2))| ...
          (E(:,2)==E(k,1))|(E(:,2)==E(k,2)));
  B(nn,k)=1;
  B(k,nn)=1;
end % adjacency matrix for edges + main diagonal
options=optimset('bintprog'); % the default options
options.Display='off'; % we change the output

% ============= We solve the MILP problem ==========
xmin=bintprog(E(:,3),-B,-ones(m,1),[],[],[],options);
nMS=find(round(xmin)); % the answer - numbers of vertexes
return

function nMS=grMaxStabSet(E,d)
% Function nMS=grMaxStabSet(E,d) solve the maximal stable set problem.
% Input parameters: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
%   d(n) (optional) - the weights of vertexes,
%     n - number of vertexes.
%     If we have only 1st parameter E, then all d=1.
% Output parameter:
%   nMS - the list of the numbers of vertexes included 
%     in the maximal (weighted) stable set.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation
if nargin<2, % we may only 1st parameter
  d=ones(n,1); % all weights =1
else
  d=d(:); % reshape to vector-column
  if length(d)<n, % the poor length
    error('The length of the vector d is poor!')
  else
    n=length(d); % Number of Vertexes
  end
end

% ============= Parameters of integer LP problem ==========
A=zeros(n,m); % for incidence matrix
A(E(:,1:2)+repmat(([1:m]'-1)*n,1,2))=1; % we fill the incidence matrix
options=optimset('bintprog'); % the default options
options.Display='off'; % we change the output

% ============= We solve the MILP problem ==========
xmin=bintprog(-d,A',ones(m,1),[],[],[],options);
nMS=find(round(xmin)); % the answer - numbers of vertexes
return

function nMM=grMaxMatch(E)
% Function nMM=grMaxMath(E) solve the maximal matching problem.
% Input parameter: 
%   E(m,2) or (m,3) - the edges of graph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of edge;
%     m - number of edges.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   nMM - the list of the numbers of edges included 
%     in the maximal (weighted) matching.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

% ============= Parameters of integer LP problem ==========
A=zeros(n,m); % for incidence matrix
A(E(:,1:2)+repmat(([1:m]'-1)*n,1,2))=1; % we fill the incidence matrix
options=optimset('bintprog'); % the default options
options.Display='off'; % we change the output

% ============= We solve the integer LP problem ==========
xmin=bintprog(-E(:,3),A,ones(n,1),[],[],[],options);
nMM=find(round(xmin)); % the answer - numbers of edges
return

function [v,mf]=grMaxFlows(E,s,t)
% Function [v,mf]=grMaxFlows(E,s,t) solve the problem 
% about the maximal flow in the network.
% Input parameters: 
%   E(m,2) or (m,3) - the arrows of digraph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of arrow;
%     m - number of arrows.
%     If we set the array E(m,2), then all weights is 1.
%   s - input (source) of the network (number of vertex);
%   t - output (sink) of the network (number of vertex).
% Output parameters: 
%   v(m,1) - vector-column of flows in the arrows;,
%   mf - the total maximal flow in the network.
% Uses the reduction to LP-problem.
% Required the Optimization Toolbox
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<3,
  error('Required three input data: E, s, t.')
end
[m,n,E] = grValidation(E); % E data validation
if ~ismember(s,E(:,1)),
  error(['From vertex number s=' num2str(s) ' does not go out any arrow!'])
end
if ~ismember(t,E(:,2)),
  error(['To vertex number t=' num2str(t) ' does not go in any arrow!'])
end

% ============= Parameters of LP problem ==========
A=zeros(n,m); % for incidence matrix
A(E(:,1)+([1:m]'-1)*n)=1; % we fill the incidence matrix
A(E(:,2)+([1:m]'-1)*n)=-1;
Aeq=A(setdiff([1:n],[s t]),:); % the balance of flows
options=optimset('linprog'); % the default options
options.Display='off'; % we change the output
options.TolX=1e-12; % we change the output
% ============= We solve the LP problem ==========
v=linprog(-A(s,:)',[],[],Aeq,zeros(size(Aeq,1),1),...
  zeros(m,1),E(:,3),[],options); % we solve the LP-problem
mf=A(s,:)*v; % the total flow
return

function nMS=grMaxComSu(E,d)
% Function nMS=grMaxComSu(E,d) solve 
% the maximal complete subgraph (clique) problem.
% Input parameters: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
%   d(n) (optional) - the weights of vertexes,
%     n - number of vertexes.
%     If we have only 1st parameter E, then all d=1.
% Output parameter:
%   nMS - the list of the numbers of vertexes included 
%     in the maximal (weighted) complete sugraph.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation
E=sort(E(:,1:2)')'; % each row in ascending order
E=unique(E,'rows'); % we delete multiple edges
E=E(setdiff([1:size(E,1)]',find((E(:,1)==E(:,2)))),:); % we delete loops
m=size(E,1);
n=max(max(E));
k1=repmat([1:n],n,1);
k2=k1';
K=[k1(:) k2(:)];
K1=reshape(K(2:end,:),n+1,2*(n-1));
K=unique(sort(reshape(K1(1:end-1,:),n*(n-1),2),2),'rows'); % clique
E1=setdiff(K,E,'rows'); % rested
if nargin<2, % we have only 1st parameter
  d=ones(n,1); % all weights =1
else
  d=d(:); % reshape to vector-column
  if length(d)<n, % the poor length
    error('The length of the vector d is poor!')
  end
end
nMS=grMaxStabSet(E1,d); % the maximal stable set problem for graph G(V,~E)
return

function [IsIsomorph,Permut]=grIsomorph(E1,E2,n)
% Function [IsIsomorph,Permut]=grIsomorph(E1,E2,n) 
% compare two simple graphs (without loops and multiple edges) 
% about their isomorphism.
% Input parameters: 
%   E1(m1,2) - the edges of 1st graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m1 - number of edges;
%   E2(m2,2) - the edges of 2nd graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m2 - number of edges;
%   n - number of vertexes (optional); if this parameter 
%     are omitted, by default n1=max(E1(:)); n2=max(E2(:));
%     if n1~=n2, the graphs are nonisomorphic.
% Output parameters:
%   IsIsomorph = true for isomorphic graphs;
%   IsIsomorph = false otherwise;
%   For nonisomorphic graphs Permut = [];
%   For isomorphic graphs Permut(1:n,1) is permutation 
%   from vertexes numbers G1 to vertexes numbers of G2.
% Needed other products: ZerOne Toolbox.
% This software may be free downloaded from site:
% http://www.apmath.spbu.ru/grafomann/download.html
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<2,
  error('There are no input data!')
end
[m1,n1,E1] = grValidation(E1); % E1 data validation
[m2,n2,E2] = grValidation(E2); % E2 data validation
IsIsomorph=false;
Permut=[];
if nargin>=3,
  n1=n;
  n2=n;
end
if (m1~=m2)|(n1~=n2),
  return;
end
A1=zeros(n1);
A1((E1(:,2)-1)*n1+E1(:,1))=1;
A1=A1+A1'; % adjacency matrix A1
A2=zeros(n2);
A2((E2(:,2)-1)*n2+E2(:,1))=1;
A2=A2+A2'; % adjacency matrix A2
[ef,P]=psimilar(A1,A2);
IsIsomorph=(ef>0);
if IsIsomorph,
  [Permut,j]=find(P');
end

function [eu,cEu]=grIsEulerian(E)
% Function eu=grIsEulerian(E) returns 1 for Eulerian graph, 
% 0.5 for semi-Eulerian and 0 otherwise.
% Input parameter: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
% Output parameter:
%   eu = 1 for Eulerian graph;
%   eu = 0.5 for semi-Eulerian graph;
%   eu = 0 otherwise.
% The graph is Eulerian, if it's connected and powers of 
% all vertexes is even.
% The graph is semi-Eulerian, if it's connected and only 
% two vertexes have odd powers, and other powers of vertexes is even.
% [eu,cEu]=grIsEulerian(E) return also 
%   cEu(m,1) - the vector-column with numbers of edges
%     included to Eulerian cycle (if eu=1) or Eulerian path (if eu=0.5)
%     and cEu=[] otherwise. The Fleury algorithm is used.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if nargin<1,
  error('There are no input data!')
end
eu = 0; % the graph is not Eulerian
cEu=[];
ncV=grComp(E); % the number of components
if max(ncV)>1, % 2 or more components
  return % the graph is not Eulerian
end
n=max(max(E(:,1:2))); % number of vertexes
E1=find([diff(sort(E(:)));1]);
p=[E1(1);diff(E1)]; % powers of vertexes
rp=rem(p,2); % remainder after division to 2
srp=sum(rp); % summa of remainders
switch srp
  case 0, % Eulerian graph
    eu=1;
  case 2, % semi-Eulerian graph
    eu=0.5;
  otherwise, % not Eulerian graph
    return
end
%=========== we find the Eulerian cycle or way =============
if nargout>1,
  if srp==0, % Eulerian graph
    v1=1; % first vertex of Eulerian cycle
  else % semi-Eulerian graph
    v1=find(rp);
    v1=v1(1); % first vertex of Eulerian way
  end
  vc=v1; % the current vertex
  m=size(E,1); % number of edges
  E1=[E(:,1:2), [1:m]']; % all edges with numbers
  while ~isempty(E1), % the Fleury algorithm
    evc=find((E1(:,1)==vc)|(E1(:,2)==vc)); % edges connected with vc
    levc=length(evc); % number of edges connected with vertex vc
    if levc==1, % only one way
      cEu=[cEu;E1(evc,3)]; % we add new edge to Eulerian cycle (way)
      vcold=vc;
      vc=sum(E1(evc,1:2))-vc; % new current vertex
      E1=E1(setdiff([1:size(E1,1)],evc),:); % we delete isolated vertex
      E2=E1(:,1:2);
      E2gv=E2>vcold;
      E2(E2gv)=E2(E2gv)-1;
      E1(:,1:2)=E2;
      if vc>vcold,
        vc=vc-1;
      end
      if v1>vcold,
        v1=v1-1;
      end
    else % several ways from vertex vc
      for k=1:levc,
        E2=E1(setdiff([1:size(E1,1)],evc(k)),:);
        ncv=grComp(E2); % number of components
        nco=max(ncv);
        if (max(ncv)==1), % this edge is not bridge
          cEu=[cEu;E1(evc(k),3)]; % we add new edge to Eulerian cycle (way)
          vc=sum(E1(evc(k),1:2))-vc; % new current vertex
          E1=E2;
          break;
        end
      end
    end
  end
end
return

function [Ec,Rad,Diam,Cv,Pv]=grEccentricity(E)
% Function Ec=grEccentricity(E) find the (weighted) 
% eccentricity of all vertexes of graph.
% Input parameter: 
%   E(m,2) or (m,3) - the edges of graph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of arrow;
%     m - number of arrows.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   Ec(1,n) - the (weighted) eccentricity of all vertexes.
% [Ec,Rad,Diam]=grEccentricity(E) find also the radius Rad and 
%   diameter Diam of the graph (the scalars).
% [Ec,Rad,Diam,Cv,Pv]=grEccentricity(E) find also 
%   the center vertexes Cv and the periphery vertexes Pv 
%   of the graph (the vector-rows with numbers of vertexes).
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if nargin<1,
  error('There are no input data!')
end
dSP=grDistances(E); % the matrix of distances
Ec=max(dSP); % the eccentricity of all vertexes
if nargout>=2,
  Rad=min(Ec); % the radius
  if nargout>=3,
    Diam=max(Ec); % the diameter
    if nargout>=4,
      Cv=find(Ec==Rad); % the center vertexes of the graph
      if nargout>=5,
        Pv=find(Ec==Diam); % the periphery vertexes of the graph
      end
    end
  end
end
return

function [dSP,sp]=grDistances(E,s,t)
% Function dSP=grDistances(E) find the distances
% between any vertexes of graph.
% Input parameter: 
%   E(m,2) or (m,3) - the edges of graph and their weight;
%     1st and 2nd elements of each row is numbers of vertexes;
%     3rd elements of each row is weight of arrow;
%     m - number of arrows.
%     If we set the array E(m,2), then all weights is 1.
% Output parameter:
%   dSP(n,n) - the symmetric matrix of destances between 
%     all vertexes (may be dSP(i,j)=inf for disconnected graph).
% [dSP,sp]=grDistances(E,s,t) - find also
% the shortest path between vertexes s and t.
% In this case output parameter sp is vector with numbers 
% of vertexes, included to shortest path between s and t.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

Ev=[E;E(:,[2 1 3])]; % all arrows and vice versa
sp=[];
if (nargin<3)|(isempty(s))|(isempty(t)),
  dSP=grShortPath(Ev); % the shortest path
else
  s=s(1);
  t=t(1);
  if s==t, % the trivial way
    dSP=grShortPath(Ev); % the shortest path
    sp=s;
  else
    [dSP,sp]=grShortPath(Ev,s,t);
  end
end
dSP=dSP-diag(diag(dSP)); % we delete the main diagonal
return

function [Dec,Ord]=grDecOrd(E)
% Function Dec=grDecOrd(E) solve 
% the problem about decomposition of the digraph
% to the sections with mutually accessed vertexes
% (strongly connected components).
% Input parameter: 
%   E(m,2) - the arrows of digraph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of arrows.
% Output parameter:
%   Dec(n,ns) - the Boolean array with numbers of vertexes.
%     n - number of vertexes;
%     ns - number of sections with mutually accessed vertexes.
%     In each column of the array Dec True value have
%     numbers of vertexes of this section.
% Other syntax: [Dec,Ord]=grDecOrd(E) also ordered all sections
%   in partial ordering. Second output parameter:
%   Ord(ns,ns) - the Boolean matrix of partial ordering
%     of sections. This matrix have right-up triangle structure.
%     If Ord(i,j)=True, then we have access 
%     from the section i (vertexes of i-st column of Dec)
%     to the section j (vertexes of j-st column of Dec).
%     In this syntax all columns of Dec(n,ns) 
%     is in partial ordering.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru
% Acknowledgements to Dr.Eduardo D. Sontag (sontag@math.rutgers.edu)
% for testing of this algorithm.

if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % data validation
% ============ Decomposition ==========================
A=eye(n);
A((E(:,2)-1)*n+E(:,1))=1; % adjacency matrix with main diagonal
A2=(A*A)>0;
while ~all(all(A2==A)),
  A=double(A2);
  A2=(A*A)>0;
end
[T,ir,jc]=unique(A.*A','rows');
Dec=T';

% ============== Partial ordering =========================
if nargout>1, 
  Ord=A(ir,ir);
  ns=size(Ord,1); % number of sections
  for it=1:ns*(ns-1)/2, % the iterations for partial ordering
    Mlow=tril(Ord,-1); % the left down trialgle
    [is,js,Mw]=find(Mlow); % we find not ordering elements
    if isempty(is), % all ordered
      break; % exit from loop for
    end
    num=[1:ns];
    num(is(1))=js(1); % we change two numbers
    num(js(1))=is(1);
    Ord=Ord(num,num);
    Dec=Dec(:,num);
  end
end  
return

function Cycles=grCycleBasis(E)
% Function Cycles=grCycleBasis(E) find 
% all independent cycles for a connected simple graph
% without loops and multiple edges
% (fundamental set of circuits).
% For loops and multiple edges you need add new vertexes.
% Input parameter: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
% Output parameter:
%   Cycles(m,m-n+1) - the Boolean array with numbers of edges.
%     n - number of vertexes;
%     m-n+1 - number of independent cycles.
%     In each column of the array Cycles True value have
%     numbers of edges of this cycle.
% Uses the addition of one edge to the spanning tree and deleting of tails.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

nMST=grMinSpanTree(E); % data validation and minimal spanning tree
E=sort(E(:,1:2)')'; % only numbers of vertexes
m=size(E,1); % number of edges
n=max(max(E)); % number of vertexes
Erest=E(setdiff([1:m],nMST),:); % rested edges
nr=m-n+1; % number of rested edges
Cycles=zeros(m,nr); % array for cycles
for k1=1:nr, % we add one independent cycle
  Ecurr=[E(nMST,:);Erest(k1,:)]; % spanning tree + one edge
  A=zeros(n);
  A((Ecurr(:,1)-1)*n+Ecurr(:,2))=1;
  A=A+A'; % the connectivity matrix
  p=sum(A); % the vertexes power
  nv=[1:n]; % numbers of vertexes
  while any(p==1), % we delete all tails
    nc=find(p>1); % rested vertexes
    A=A(nc,nc); % new connectivity matrix
    nv=nv(nc); % rested numbers of vertexes
    p=sum(A); % new powers
  end
  [i1,j1]=find(A);
  incedg=nv(unique(sort([i1 j1]')','rows')); % included edges
  Cycles(:,k1)=ismember(E,incedg,'rows'); % current column
end
return

function ncV=grComp(E,n)
% Function ncV=grComp(E,n) find all components of the graph.
% Input parameter: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
%   n - number of vertexes ( optional, by default n=max(max(E)) ).
%     This input parameter is needed, if last vertexes is isolated.
% Output parameter:
%   ncV(n,1) - the the vector-column with the number of component 
%     for each vertex;
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

if nargin<1,
  error('There are no input data!')
end
[m,n1,E1] = grValidation(E); % data validation
E2=[E1(:,1:2);E1(:,[2 1])]; % all arrows and vice versa
Dec=grDecOrd(E2); % the components
ncV=sum(Dec*diag([1:size(Dec,2)]),2); % the numbers of components
if (nargin>1)&(n>n1), % last isolated vertexes
  ncV=[ncV;[1:n-n1]'+max(ncV)];
end
return

function nCol=grColVerOld(E)
% function nCol=grColVer(E) solve the color graph problem
% for vertexes of the graph.
% Input parameter: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
% Output parameter:
%   nCol(n,1) - the list of the colors of vertexes.
% Uses the sequential deleting of the maximal stable sets.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation
E=sort(E(:,1:2)')'; % each row in ascending order
E=unique(E,'rows'); % we delete multiple edges
E=E(setdiff([1:size(E,1)]',find((E(:,1)==E(:,2)))),:); % we delete loops
nCol=zeros(n,1); % initial value
% ============= Main cycle with MaxStabSet deleting ====
while any(nCol==0),
  nv=find(nCol==0); % uncolored vertexes
  E1=E(find(ismember(E(:,1),nv)&ismember(E(:,2),nv)),:); % it's edges
  if isempty(E1),
    nCol(find(nCol==0))=max(nCol)+1; % the last color
    break;
  end
  nvs=unique(E1(:)); % all vertexes
  for kk=1:length(nvs),
    E1(find(E1==nvs(kk)))=kk;
  end
  nMS=grMaxStabSet(E1); % the maximal stable set
  nCol(nvs(nMS))=max(nCol)+1; % the next color
end
return

function nCol=grColVer(E)
% function nCol=grColVer(E) solve the color graph problem
% for vertexes of the graph.
% Input parameter: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
% Output parameter:
%   nCol(n,1) - the list of the colors of vertexes.
% Uses the reduction to integer LP-problem.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation

E=sort(E(:,1:2)')'; % each row in ascending order
E=unique(E,'rows'); % we delete multiple edges
E=E(setdiff([1:size(E,1)]',find((E(:,1)==E(:,2)))),:); % we delete loops

% ============= Parameters of integer LP problem ==========
A1=zeros(2*m,n*n);
for kk=1:m,
  A1(2*kk-1,(E(kk,1)-1)*n+1:E(kk,1)*n)=1;
  A1(2*kk-1,(E(kk,2)-1)*n+1:E(kk,2)*n)=-1;
  A1(2*kk,(E(kk,1)-1)*n+1:E(kk,1)*n)=-1;
  A1(2*kk,(E(kk,2)-1)*n+1:E(kk,2)*n)=1;
end
A=[zeros(2*m,n),A1,...
  reshape([-n*reshape(eye(m),1,m*m);n*reshape(eye(m),1,m*m)],2*m,m);...
  -ones(n),reshape(repmat(reshape(eye(n),1,n*n),n,1),n*n,n)',zeros(n,m)];
b=[reshape([-ones(1,m);(n-1)*ones(1,m)],2*m,1);zeros(n,1)];
c=[ones(n,1);zeros(n*n+m,1)];
options=optimset('bintprog'); % the default options
options.Display='off'; % we change the output

% ============= We solve the integer LP problem ==========
xmin=round(bintprog(c,A,b,[],[],[],options));
nCol=(sum(reshape(xmin(n+1:n*(n+1)),n,n))+1)'
return

function mCol=grColEdge(E)
% function mCol=grColEdge(E) solve the color graph problem
% for edges of the graph.
% Input parameter: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
% Output parameter:
%   mCol(m,1) - the list of the colors of edges.
% Uses the sequential deleting of the maximal matching sets.
% Required the Optimization Toolbox v.3.0.1 or over.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

% ============= Input data validation ==================
if nargin<1,
  error('There are no input data!')
end
[m,n,E] = grValidation(E); % E data validation
E=[E(:,1:2),[1:m]']; % numbers of vertexes and numbers of edges
mCol=zeros(m,1); % initial value
% ============= Main cycle with MaxMatch deleting ====
while any(mCol==0),
  ne=find(mCol==0); % uncolored edges
  E1=E(ne,:); % it's edges
  nMM=grMaxMatch(E1(:,1:2)); % the maximal matching
  mCol(E1(nMM,3))=max(mCol)+1; % the next colorend
end
return


function CoCycles=grCoCycleBasis(E)
% Function CoCycles=grCoCycleBasis(E) find 
% the cocycle basis for a connected simple graph
% without loops and multiple edges
% (fundamental set of cut-sets).
% For loops and multiple edges you need add new vertexes.
% Input parameter: 
%   E(m,2) - the edges of graph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of edges.
% Output parameter:
%   CoCycles(m,n-1) - the Boolean array with numbers of edges.
%     n - number of vertexes;
%     n-1 - number of independent cocycles.
%     In each column of the array CoCycles True value have
%     numbers of edges of this cocycle.
% Uses the deletion one edge from spanning tree.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

nMST=grMinSpanTree(E); % data validation and minimal spanning tree
E=E(:,1:2); % only numbers of vertexes
Emst=E(nMST,:); % edges of minimal spanning tree
m=size(E,1); % number of edges
n=max(max(E)); % number of vertexes
CoCycles=zeros(m,n-1); % array for cocycles
for k1=1:n-1, % we add one independent cocycle
  ncV=grComp(Emst(setdiff([1:n-1],k1),:),n); % two components
  n1=find(ncV==1); % the vertexes of 1st component
  n2=find(ncV==2); % the vertexes of 2nd component
  CoCycles(find((ismember(E(:,1),n1)&ismember(E(:,2),n2))|...
    (ismember(E(:,1),n2)&ismember(E(:,2),n1))),k1)=1;
end
return


function BG=grBase(E)
% Function BG=grBase(E) find all bases of digraph. 
% Input parameter: 
%   E(m,2) - the arrows of digraph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of arrows.
% Output parameter:
%   BG(nb,nv) - the array with numbers of vertexes.
%     nb - number of bases;
%     nv - number of vertexes in each base.
%     In each row of the array BG is numbers 
%     of vertexes of this base.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

[d,po]=grDecOrd(E); % decomposition
ngr=1; % numbers of groups
ua=po(1,:);
while ~all(ua),
  ung=find(~ua);
  ngr=[ngr ung(1)];
  ua=ua | po(ung(1),:);
end
levels=sum(d(:,ngr)); % full-factorial designs
ssize = prod(levels);
ncycles = ssize;
cols = length(levels);
design = zeros(ssize,cols);
for k = 1:cols,
  settings = (1:levels(k));
  ncycles = ncycles./levels(k);
  nreps = ssize./(ncycles*levels(k));
  settings = settings(ones(1,nreps),:);
  settings = settings(:);
  settings = settings(:,ones(1,ncycles));
  design(:,k) = settings(:);
end
for k=1:size(design,2), % we change
  ff=find(d(:,ngr(k))); % number of factor
  BG(:,k)=ff(design(:,k)); % to number of vertex
end
return

function CBG=grCoBase(E)
% Function CBG=grCoBase(E) find all contrabases of digraph. 
% Input parameter: 
%   E(m,2) - the arrows of digraph;
%     1st and 2nd elements of each row is numbers of vertexes;
%     m - number of arrows.
% Output parameter:
%   CBG(ncb,nv) - the array with numbers of vertexes.
%     ncb - number of contrabasis;
%     nv - number of vertexes in each contrabasis.
%     In each row of the array BG is numbers 
%     of vertexes of this contrabasis.
% Author: Sergii Iglin
% e-mail: siglin@yandex.ru
% personal page: http://iglin.exponenta.ru

E=E(:,[2 1]);
CBG=grBase(E);
return


function fc=findcycles(G)

numNodes = size(G, 1); 
ofc=zeros(60,20);

j=1;
for n = 1: 25
    nn = int16(n);
    G
   [D,P] = graphtraverse(G, n);
   for d = D
       if G(d,n)
           graphpred2path(P,d);
           s=size(graphpred2path(P,d));
           fc=graphpred2path(P,d);
           for i=1:s(2)
               ofc(j,i)=fc(i);
           end
           j=j+1;
%            fc(j,15)=graphpred2path(P,d);
%            j=j+1;
           
       end
   end
   G(n,:)=0; 
   fc=ofc;
end



% --- Executes just before lab7 is made visible.
function lab7_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lab7 (see VARARGIN)

% Choose default command line output for lab7
handles.output = hObject;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
filename = 'Map.xlsx';
sheet = 1;
xlRange = 'B2:L12';

% A=xlsread(filename,sheet,xlRange);
A = load('values.txt');
set(handles.uitable3,'Data',A);

%xlRange = 'A15:A25';
%XX=xlsread(filename,sheet,xlRange);
%set(handles.listbox5,'UserData',XX);

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = lab7_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes when entered data in editable cell(s) in uitable3.
function uitable3_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable3 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
%data=get(handles.uitable3,'Data');
%A_n=data;
%     A_n=cell2mat(data);




% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data=get(handles.uitable3,'Data');
A_n=data;

%create the matrix of the connections without zeros

global E sA;
sA=size(A_n);
E=zeros(1,3);
V=zeros(1,2);

iv=1;
for i=1:sA(1)
    for j=1:sA(2)
        if A_n(i,j) ~=0
            E(iv,1)=i;
            E(iv,2)=j;
            E(iv,3)=A_n(i,j);
         iv=iv+1;
        end
    end
end
 R=1;
%R=500;
for i=1:sA(1)
    V(i,1)=R*cos(i*2*pi/sA(1));
    V(i,2)=R*sin(i*2*pi/sA(1));
end 

   
kind='d';
sa=1.1;
vkind='%d';
ekind='%3.3f';
h=grPlot(V,E,kind,vkind,ekind,sa);
                



% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

    


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global E sA;

data=get(handles.uitable3,'Data');
    A_n=data;
    sA=size(A_n);
    E=zeros(1,3);
    
    iv=1;
    for i=1:sA(1)
        for j=1:sA(2)
            if A_n(i,j) ~=0
                E(iv,1)=i;
                E(iv,2)=j;
                E(iv,3)=A_n(i,j);
             iv=iv+1;
            end
        end
    end
    
 sE=size(E);
       
v1=E(1:sE(1),1)';
v2=E(1:sE(1),2)';

G = sparse( v1,v2,true,sE(1),sE(1));

fc = findcycles(G) %array 60x20 of all cycles
 
 i=1;
 c=[];
 m=1; %indicator of even cycle
 evenc=zeros(40,10); %all even cycles
 u=0;
 ec=0; %counter of even cycles
 while fc(i,1) ~=0
i=i+1;
     ind = find(fc(i,:), 1, 'last');
     c=[c,fc(i,1:ind)] %one cycle

      sC=size(c);
     for i2=1:sC(2)
         for k=1:sE(1)
             if i2==sC(2)
                 if E(k,1)==c(i2) && E(k,2)==c(1)
                 m=m*E(k,3); % + or - relationship
                 end
             else
                if E(k,1)==c(i2) && E(k,2)==c(i2+1)
                 m=m*E(k,3); % + or - relationship
                end
             end
         end
     end
     
     if m>0 && ~isempty(c)
         m;
          u=u+1;
          ec=ec+1;
         for i3=1:sC(2)
             evenc(u,i3)=c(i3);
         end
   end
     
     m=1;      
     c=[];
     f=false;
 end
  ecS=num2str(ec);
 odc=i-1-ec;
 odcS=num2str(odc);
 set(handles.edit1, 'String',ecS);
 set(handles.edit2, 'String',odcS);

 indEV=find(evenc(:,1),1,'last');

 s0=' ->';
 s='';
 
 if evenc(1,1)~=0
 for y=1:indEV
     indJ=find(evenc(y,:),1,'last');
     for q=1:indJ
     s1=strcat(num2str(evenc(y,q)),s0);
     s=strcat(s,s1);
     end
     if y==1
         old_str='';
     else
    old_str = get(handles.lb3, 'String' );
     end
     if isempty(evenc)
         new_str='0';
          set(handles.lb3, 'String', new_str);
     else
	new_str = strvcat(old_str,s); %#ok<VCAT>
    set(handles.lb3, 'String', new_str);
     end
     s='';
 end
 else
       set(handles.lb3, 'String', '0');
 end


d=eig(A_n);
da=abs(d);
mxl=0;
sD=size(da);
 for y1=1:sD(1)
     s1=strcat(num2str(da(y1)));
 if da(y1)>mxl
     mxl=da(y1);
 end
if y1==1
    old_str='';
else
    old_str = get(handles.listbox4, 'String' );
end
    new_str = strvcat(old_str,s1); %#ok<VCAT>
    set(handles.listbox4, 'String', new_str);
 end 
 mxlS=num2str(mxl);
set(handles.edit4, 'String', mxlS);

iOK = true;

if mxl<=1
    set(handles.checkbox2, 'Value', 1);
    if mxl<1
        set(handles.checkbox3, 'Value', 1);
        set(handles.checkbox5, 'Value', 1);
    else
        set(handles.checkbox3, 'Value', 0);
        set(handles.checkbox5, 'Value', 0);
        iOK = false;
    end
else
    set(handles.checkbox2, 'Value', 0);
    set(handles.checkbox3, 'Value', 0);
    set(handles.checkbox5, 'Value', 0);
    iOK = false;
end
if ec>1
    set(handles.checkbox1, 'Value', 0);
    iOK = false;
else
     set(handles.checkbox1, 'Value', 1);
end
% if iOK 
%     set(handles.OKtext, 'Visible', 'on');
% else
%     set(handles.OKtext, 'Visible', 'off');
% end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lb3.
function lb3_Callback(hObject, eventdata, handles)
% hObject    handle to lb3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb3


% --- Executes during object creation, after setting all properties.
function lb3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.~~~~~~~~~~~~~~~~~~~~~~~~~~
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% nm=14;
% nms=num2str(nm);
%%%dados=get(handles.uitable3, 'Data');
% c=get(handles.uitable3,'ColumnName');
% r=get(handles.uitable3, 'RowName');
% 
% 
% c(size(c)+1,1)=nms;
% r(size(r)+1,1)=nms;
% size(dados,1)
%%%linha=dados (size(dados,1),:);  %last line
%%%Usuario_1=cat(1,dados,linha);
%%%colha=Usuario_1 (:,size(dados,2)); %last column
%%%Usuario_2= cat(2,Usuario_1,colha);
% colnha= Usuario_1 (:,size(Usuario_1,2));
% Usuario_2=cat(1,Usuario_1,colnha);
% nm=nm+1;
%%%set(handles.uitable3,'data', Usuario_2, 'ColumnEditable',true);


% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5


% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
