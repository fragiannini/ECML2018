%supp sia tot sup, cioè quando so una sup per un pred la so per tutti gli
%altri

%% PARAMETERS


%examples

global dom valDom S

comp=4;
%L=[];
L=L(:,2:comp);
%U=U(:,1:4);
U=[];
%U=U(:,1);   %version with only one point and one supervision
% U=[U,[1;0]];
target=target(:,2:comp);

l=size(L,2);
u=size(U,2);

%target=round(rand(3,l));
%target=ones(3,l); %sono da metterci le vere sup per ogni predicato

%U=round(rand(2,l)*10)*1/10;
S=[L,U]; %SUPPONENDO L e U DISGIUNTI
%S=unique([L,U]',t'rows','stable')';
s=size(S,2);

%dom

n=0.01;
[domainx,domainy]=meshgrid(0:n:1,0:n:1);
sizedom=size(domainx,1)*size(domainy,2);
dom=zeros(2,sizedom);
for i=1:sizedom
dom(:,i)=[domainx(i);domainy(i)];
end


stampa=0;
t=3;
v=3;
eta=0.00001;


degree=2;
sigma=1;

options.Kernel='poly';
options.KernelParam=degree;

% options.Kernel='rbf';
% options.KernelParam=sigma;


K=calckernel(options,S');

Kdom=calckernel(options,dom',S');

%% INITIALIZZATIONS


lb=[];
ub=[];

H=zeros(t*s,t*s);
c=0;
for i=1:t
    H(1+c:s+c,1+c:s+c)=K;
    c=c+s;
end

f=zeros(t*s,1);

A=zeros();
b=zeros();

%A=zeros(t*l+2*t*s+v*s,t*s);
%b=zeros(t*l+2*t*s+v*s,1);

%% P-CONSTRAINTS

c=0;
for i=1:t
    A(1+c:c+l,1+(i-1)*s:i*s)=-2*repmat(target(i,:)',1,s).*K(1:l,:);
    b(1+c:c+l,1)=-1-target(i,:)';
    c=c+l;
end

%% C-CONSTRAINTS

if 1
  for i=1:t
      A(1+c:s+c,1+(i-1)*s:i*s)=K;
      A(1+s+c:2*s+c,1+(i-1)*s:i*s)=-K;
      b(1+c:c+s,1)=1;
      b(1+s+c:2*s+c,1)=0;
      c=c+2*s;
  end
end

%% L-CONSTRAINTS

if 1
    A(c+1:c+s,:)=[K,-K,zeros(s,s)]; %1
    b(c+1:c+s,1)=zeros(s,1);
    c=c+s;
    A(c+1:c+s,:)=[zeros(s,s),K,-K]; %2
    b(c+1:c+s,1)=zeros(s,1);
    c=c+s;
    A(c+1:c+s,:)=[K,zeros(s,s),-K]; %1
    b(c+1:c+s,1)=zeros(s,1);
    c=c+s;
end
    
%% OPTIM + VAL



[x,fval,exitflag,output,lambda] = quadprog(H+eye(size(H,1))*eta,f,A,b,[],[],lb,ub);


valDom=zeros(t,sizedom);

%valutazione sui punti di DOM
c=0;
for i=1:t
    valDom(i,:)=x(1+c:s+c)'*Kdom;
    c=c+s;
end

%% PRINT 

if stampa==1

figure

%legend({'NoLog','LogConv','LogNonConv'},'FontSize',12,'FontWeight','bold')

grand=100;
ax1 = subplot(3,1,1);
scatter(dom(1,:),dom(2,:),grand,valDom(1,:),'s','filled')
colormap(ax1,parula)
colorbar
xlabel('p_1','FontSize',12)
caxis([0,1])

%title('Best Objective Functions \it','FontSize',16,'FontWeight','bold')


% hold on
% scatter(S(1,1),S(2,1),grand,'black','s','filled')
hold on
scatter(S(1,1:size(S,2)),S(2,1:size(S,2)),grand,'black','s','filled')

ax2 = subplot(3,1,2);
scatter(dom(1,:),dom(2,:),grand,valDom(2,:),'s','filled')
colormap(ax2,parula)
colorbar
xlabel('p_2','FontSize',12)
caxis([0,1])

% hold on
% scatter(S(1,1),S(2,1),grand,'black','s','filled')
hold on
scatter(S(1,1:size(S,2)),S(2,1:size(S,2)),grand,'black','s','filled')

ax3 = subplot(3,1,3);
scatter(dom(1,:),dom(2,:),grand,valDom(3,:),'s','filled')
colormap(ax3,parula)
colorbar
xlabel('p_3','FontSize',12)
caxis([0,1])

% hold on
% scatter(S(1,1:1),S(2,1:1),grand,'black','s','filled')
hold on
scatter(S(1,1:size(S,2)),S(2,1:size(S,2)),grand,'black','s','filled')


axis([0, 1, 0, 1])
end


    
    