close all, clear all, clc

%% Fecha de Modificación 10 de Mayo de 2014

%% Dimensiones del Mecanismo

r1=5;
r2=2;
r3=5;
r4=4.5;
theta1=pi;
theta2=0:pi/180:2*pi;
dtheta1=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Valores de Entrada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

dr1=0.02;
dr2=0.01;
dr3=0.02;
dr4=0.015;

%%%%%%%%%%%%%%%%%%%%%

c=-0.0017; %Variation in theta2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Valores de Entrada 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

rd1=r1+dr1;
rd2=r2+dr2;
rd3=r3+dr3;
rd4=r4+dr4;
theta2d=theta2+c;

for j=1:1
   
    
 dtheta2=c;
 dt2(j,1:361)=c;


for i=1:361
  
     
 %% Cálculo de los valores nominales de theta3 y theta4
 %% Sin variaciones

 theta_2(i)=theta2(i);
  [theta_3 theta_4]=angulos4barras(r2,r3,r4,r1,theta_2(i));
  theta3(i)=theta_3;
  theta4(i)=theta_4;
   
  %% Con variaciones
  theta_2d(i)=theta2d(i);
  [theta_3d theta_4d]=angulos4barras(rd2,rd3,rd4,rd1,theta_2d(i));
  theta3D(i)=theta_3d;
  theta4D(i)=theta_4d;
    
    
%% Método de Leishman R.C. & Chase K.W.(2010) 
%% Matriz A    
      
A=[cos(theta1) cos(theta2(i)) cos(theta3(i)) -cos(theta4(i)) -r1*sin(theta1) -r2*sin(theta2(i)); sin(theta1) sin(theta2(i)) sin(theta3(i)) -sin(theta4(i)) r1*cos(theta1) r2*cos(theta2(i))];
B=[-r3*sin(theta3(i)) r4*sin(theta4(i));r3*cos(theta3(i)) -r4*cos(theta4(i))];
dX=[dr1;dr2;dr3;dr4;dtheta1;dtheta2];


%% Results of variations dteta3 y dteta4
dUl(:,i)=((-inv(B)*A))*dX;
H=((-inv(B))*A);

Lamaxl(i)=max(abs(eig(B)));
Laminl(i)=min(abs(eig(B)));



%% Inicializando la variable dU2l
dU2l(1:2,i)=0;

%% Suma Cuaddratica
for k=1:4 
dU2l(1,i)=dU2l(1,i)+(H(1,k)*dX(k,1)).^2;
dU2l(2,i)=dU2l(2,i)+(H(2,k)*dX(k,1)).^2;
end

%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Metodo Propuesto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
AA=[-sin(theta1) -sin(theta2(i)) -sin(theta3(i)) sin(theta4(i)) -r1*cos(theta1)-sin(theta1) -r2*cos(theta2(i))-sin(theta2(i));cos(theta1) cos(theta2(i)) cos(theta3(i)) -cos(theta4(i)) -r1*sin(theta1)+cos(theta1) -r2*sin(theta2(i))+cos(theta2(i))];
BB=[-sin(theta3(i))-r3*cos(theta3(i)) sin(theta4(i))+r4*cos(theta4(i));-r3*sin(theta3(i))+cos(theta3(i)) -cos(theta4(i))+r4*sin(theta4(i))];
dXX=[dr1;dr2;dr3;dr4;dtheta1;dtheta2];

%% Variaciones
dU(:,i)=(-inv(BB))*(AA*dXX);
dU2(1:2,i)=0;
H=(inv(BB)*AA);
Lamax(i)=max(abs(eig(BB)));
Lamin(i)=min(abs(eig(BB)));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:4 
dU2(1,i)=dU2(1,i)+(H(1,k)*dXX(k,1)).^2;
dU2(2,i)=dU2(2,i)+(H(2,k)*dXX(k,1)).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Minimization Problem Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EQUALITIES

  
Anl=[-sin(theta1) -sin(theta2(i)) -sin(theta3(i)) sin(theta4(i)) -r1*cos(theta1)-sin(theta1);...
    cos(theta1) cos(theta2(i)) cos(theta3(i)) -cos(theta4(i)) -r1*sin(theta1)+cos(theta1)];
Bnl=[-r2*cos(theta2(i))-sin(theta2(i)) -sin(theta3(i))-r3*cos(theta3(i)) sin(theta4(i))+r4*cos(theta4(i)) ;...
    -r2*sin(theta2(i))+cos(theta2(i)) -r3*sin(theta3(i))+cos(theta3(i)) -cos(theta4(i))+r4*sin(theta4(i)) ];
%Anl=[cos(theta1) cos(theta2(i)) cos(theta3(i)) -cos(theta4(i)) -r1*sin(theta1) ; sin(theta1) sin(theta2(i)) sin(theta3(i)) -sin(theta4(i)) r1*cos(theta1)];
%Bnl=[-r2*sin(theta2(i)) -r3*sin(theta3(i)) r4*sin(theta4(i)); r2*cos(theta2(i)) r3*cos(theta3(i)) -r4*cos(theta4(i))];
dXnl=[dr1;dr2;dr3;dr4;dtheta1];

%% Equalities Variables Bnl=bb

bb=-Anl*dXnl
Bnl(:,4)=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Minimization Problem Part 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Function
AAn=[cos(theta1) cos(theta2(i)) cos(theta3(i)) -cos(theta4(i)) -r1*sin(theta1) ;sin(theta1) sin(theta2(i)) sin(theta3(i)) -sin(theta4(i)) r1*cos(theta1)];
BBn=[-r2*sin(theta2(i)) -r3*sin(theta3(i)) r4*sin(theta4(i));r2*cos(theta2(i)) r3*cos(theta3(i)) -r4*cos(theta4(i))];

dXXn=[dr1;dr2;dr3;dr4;dtheta1];


%% bT*X
ccn(1:4)=[BBn(1,:)+BBn(2,:) 1];
% Four variable X4=Vn
Vn=AAn(1,:)*dXXn+AAn(2,:)*dXXn;
%% Completing Equalities Matrix 
bb(3,1)=Vn;
Bnl(3,4)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% optimisation Solution
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ct= ccn';
a=Bnl;
blc=bb;
buc=bb;
blx =[-abs(dtheta2);-2*abs(dtheta2);-2*abs(dtheta2);Vn];
bux = [abs(dtheta2);2*abs(dtheta2);2*abs(dtheta2);Vn];

%% Maximization and Minimization
[res] = msklpopt(-ct,a,blc,buc,blx,bux);
sol=res.sol;
XX(:,i)=sol.itr.xx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

dt3(j,:)=dU(1,:);
dt4(j,:)=dU(2,:);
dU2=sqrt(dU2);
dt3l(j,:)=dUl(1,:);
dt4l(j,:)=dUl(2,:);
dU2l=sqrt(dU2l);

end
%% print Results
Results_12_Abril_2015
