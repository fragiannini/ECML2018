% This script produces the solution of Example 2 in the paper. We make use 
% of the MATLAB function 'quadprog' that is specific for quadratic 
% programming problems and provides as additional output a vector of 
% Lagrange multipliers satisfying the KKT-conditions for the given solution

%% INITIALIZATION

% Supervisions and a domain 'dom' for the evaluation are loaded. The matrix 'L'
% contains the points whose labels (for every predicate) are collected in
% the matrix 'target'
load('data2.mat')  


U=[]; % This set contains the unsupervised sample

l=size(L,2);
u=size(U,2);
sdom=size(dom,2);

S=[L,U]; % This set contains the overall sample for the predicates
s=size(S,2);

print=1; % If this parameter is equal to 1 then the optimal solutions are plotted
J=3; % Number of learnable predicates
v=3; % Number of considered constrains
eta=0.00001; % To avoid numerical issues

% The following parameters define a common kernel for the predicates and we
% use the script calckernel to get the Gram matrix on the sets 'dom' and
% 'S'
degree=2;

options.Kernel='poly';
options.KernelParam=degree;

K=calckernel(options,S');

Kdom=calckernel(options,dom',S');

%% PARAMETERS FOR QUADPROG

lb=[];
ub=[];

H=zeros(J*s,J*s);
c=0;
for i=1:J
    H(1+c:s+c,1+c:s+c)=K;
    c=c+s;
end

f=zeros(J*s,1);

A=zeros();
b=zeros();

%% POINTWISE CONSTRAINTS
% These constraints enforce the available supervisions

c=0;
for i=1:J
    A(1+c:c+l,1+(i-1)*s:i*s)=-2*repmat(target(i,:)',1,s).*K(1:l,:);
    b(1+c:c+l,1)=-1-target(i,:)';
    c=c+l;
end

%% CONSISTENCY CONSTRAINTS
% These constraints enforce the predicate to be limited between 0 and 1

  for i=1:J
      A(1+c:s+c,1+(i-1)*s:i*s)=K;
      A(1+s+c:2*s+c,1+(i-1)*s:i*s)=-K;
      b(1+c:c+s,1)=1;
      b(1+s+c:2*s+c,1)=0;
      c=c+2*s;
  end

%% LOGICAL CONSTRAINTS
% These are the linear constraints corresponding to the logical formulas
% p_1(x)-> p_2(x), p_2(x)-> p_3(x), p_1(x)-> p_3(x)

    A(c+1:c+s,:)=[K,-K,zeros(s,s)]; %1
    b(c+1:c+s,1)=zeros(s,1);
    c=c+s;
    A(c+1:c+s,:)=[zeros(s,s),K,-K]; %2
    b(c+1:c+s,1)=zeros(s,1);
    c=c+s;
    A(c+1:c+s,:)=[K,zeros(s,s),-K]; %1
    b(c+1:c+s,1)=zeros(s,1);
    c=c+s;
    
%% OPTIMIZATION AND EVALUATION


[x,fval,exitflag,output,lambda] = quadprog(H+eye(size(H,1))*eta,f,A,b,[],[],lb,ub);


% Evaluation of the predicates on the domain
valDom=zeros(J,sdom);

c=0;
for i=1:J
    valDom(i,:)=x(1+c:s+c)'*Kdom;
    c=c+s;
end

%% PRINT 

if print==1

figure

grand=100;
ax1 = subplot(3,1,1);
scatter(dom(1,:),dom(2,:),grand,valDom(1,:),'s','filled')
colormap(ax1,parula)
colorbar
xlabel('p_1','FontSize',12)
caxis([0,1])

hold on
scatter(S(1,1:s),S(2,1:s),grand,'black','s','filled')


ax2 = subplot(3,1,2);
scatter(dom(1,:),dom(2,:),grand,valDom(2,:),'s','filled')
colormap(ax2,parula)
colorbar
xlabel('p_2','FontSize',12)
caxis([0,1])

hold on
scatter(S(1,1:s),S(2,1:s),grand,'black','s','filled')


ax3 = subplot(3,1,3);
scatter(dom(1,:),dom(2,:),grand,valDom(3,:),'s','filled')
colormap(ax3,parula)
colorbar
xlabel('p_3','FontSize',12)
caxis([0,1])

hold on
scatter(S(1,1:s),S(2,1:s),grand,'black','s','filled')

axis([0, 1, 0, 1])
end


    
    