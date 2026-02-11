clear all
a=30* 1e-3; %单位 m
b=3* 1e-3;  %单位 m
rhoab=[1:0.25:16]';
a=b.*rhoab;

%sa 气体的占比
sa=0.7;
%上部分水分占比
sw1=(1-sa)/2;
%下部分水分占比
sw2=(1-sa)/2;

dt=0.1;
t=[0:dt:9]';
%t=5;
q=1.028;
%冻结温度
Ti=-0;
betaT=0.09-0.045*Ti/22;
si=betaT.*(1-exp(-q*t));
ET=66*(1-0.012*Ti);
%冻结压力
Pi=ET.*si+0.001;
plot(t,Pi);
%切向摩擦力
oi=72.86-0.2476*(273.15+Ti);%冰的单轴压缩强度 MPa
or=60;%岩石的单轴压缩强度 MPa
JRC=0.1;
Rs=sqrt(0.0001^2+0.003^2)/0.003;

taus=(0.52+0.089*sqrt(abs(Ti) ))*Rs;
tauP=Pi.*tan( 0.04*exp(1)^ (4.689*( 1-(273.15+Ti)/273.15 ) ) + JRC .*    log10( (oi^2+or^2)./ ( (oi+or).*Pi )  ) )  ;
taui=taus+tauP;
n=size(Pi,1);

Pa=1.01;
A1=Pi(5,1);
A2=taui(n,1);

lmin=sqrt( (A1-Pa)^2.*a.^2.*b^2 ./( a.^2*A2^2 + (A1-Pa)^2.*b^2      ))            ;

% a(37,1)=30mm
j=37;



L=a(j,1)-lmin(j,1);
%冻结前 上部分水分的体积
L1=L*sw1;
%冻结前 气体的体积长度
L2=L*sa;
%冻结前 下部分水分的体积
L3=L*sw2;

%冻结前 气体的体积长度
X1=L1+lmin(j,1); 
X2=L1+lmin(j,1)+L2; 
f = @(x) 2*b/a(j,1)*sqrt(a(j,1)^2-x.^2);  % 函数句柄

Va1 = integral(f, X1, X2);  % 对 f(x) 从 a 到 b 积分，冻结前 气体的体积长度

%上部分水分的 冰压力、切向力作用长度 
L4=L1.*(1+si);
%下部分水分的 冰压力、切向力作用长度 
L6=L3.*(1+si);
%气体压力z作用长度、
L5=a(j,1)-(L4+lmin(j,1))-L6;

%上部分水分的 冰压力、切向力作用坐标终点
X3=(L4+lmin(j,1)); %冰压力、切向力作用长度 
%气体压力作用坐标终点
X4=a(j,1)-L6;
%下部分水分的 冰压力、切向力作用坐标终点
X5=a(j,1); 


%冻结后 气体的体积
for i=1:1:n
Va2(i,1)= integral(f, X3(i,1),X4(i,1) );
end

%气体参数
T0=273.15;  %0℃
Pa0=101325; %1 atm标准大气压强 单位Pa
R=8.314;
Ma=28.96;
ma=Pa0*Va1*Ma/(R*T0);

T1=273.15+Ti;  %0℃
ref3=ma*R*T1/Ma/1e5;
%气体方程求解气体压力
Pa1=ref3./Va2;


%平衡方程求解气体压力
du=2*b/a(j,1).*sqrt(a(j,1)^2-X3.^2);
dp=2*b/a(j,1).*sqrt(a(j,1)^2-(X3+L5).^2);

%Pa2=Pi-a(j,1)/b.*taui.*X2./sqrt(a(j,1)^2-X2.^2);
% taui=(Pi.*de-Pa1.*de)./X2/2;
%平衡方程8
P1=Pi.*du;
P2=3*taui.*L4;
P3=Pa1.*du;
Z1=P1-P2-P3;
%平衡方程9
P4=Pi.*dp;
P5=5*taui.*L6;
P6=Pa1.*dp;
Z2=P4-P5-P6;

%Pa3=Pi(n,1)-a(j,1)/b*taui(n,1)*(L2+lmin(j,1))/sqrt( a(j,1)^2-(L2+lmin(j,1))^2 );

plot(t,Z1);
hold on 
plot(t,Z2);
