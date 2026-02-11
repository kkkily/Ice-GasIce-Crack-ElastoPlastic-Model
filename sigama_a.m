clear all
% 裂隙长度
a=[10:0.25:40]'.*1e-3; %单位mm
b=3e-3;
sa=0.5;
dt=0.05;
Ti=0;
sw1=0.25;%上部分水分占比
sa= 0.5; %sa 气体的占比
sw2=0.25;%下部分水分占比


%saa=[0.1:0.005:0.7]';
n=size(a,1);
for i=1:1:n
[L4(:,i),L5(:,i),L6(:,i),Pi(:,i),taui(:,i),Pa(:,i),t]=refs(a(i,1),b,dt,Ti,sw1,sa,sw2);

[KI(:,i),KII(:,i),Ki(:,i),Ka(:,i),Kt(:,i),s4]=KK(a(i,1),n,L4(:,i),L5(:,i),L6(:,i),Pi(:,i),taui(:,i),Pa(:,i));
end


oys=2e8; %岩石屈服强度
% rp随a的变化，取t=6h；
KIa=real(1*KI(n,:)');
KIIa=real(1*KII(n,:)');
for i=1:1:n
    Rpa(:,i)=rp(KIa(i,1),KIIa(i,1),oys);
end
sigam=[0:1:360]'.*(pi/180);

Rpa=real(Rpa);

% rp随t的变化，取a=40mm；
KIt=KI(:,n);
KIIt=KII(:,n);
for i=1:1:n
    Rpt(:,i)=rp(KIt(i,1),KIIt(i,1),oys);
end

R1=Rpa(3,:)';



%周向应力o1、起裂角sita_sa、起裂应力sigama_sa

rc=0.2;
sigam0=[0:1:360]';
MN=size(KIa,1);
oo=[0:1:360]'.*(pi/180);
MM=size(oo,1);

for j=1:1:MN
    KI=KIa(j,1);
    KII=KIIa(j,1);
   for i=1:1:MM   
    o1(i,1)=KI/( 2*sqrt(2*pi*rc) )*cos(oo(i,1)/2)*(1+ cos(oo(i,1)) )-3*KII/( 2*sqrt(2*pi*rc) )*cos(oo(i,1)/2)*sin(oo(i,1));
   end
   
    [~,B]=max(abs(o1));
    sita_sa1(j,1)=sigam0(B,1);
    sigama_sa1(j,1)=o1(B,1);
    clear B
end
o1=real(o1);
sigama_sa1=abs(sigama_sa1);

figure
plot(sigam0,o1)
figure
plot(a,sita_sa1)
% 周向应力云图

dN=300;
dn=300;

jj=121;
rpa=Rpa(:,jj);%j=1,a=10mm
sigam=[0:1:360]'.*(pi/180);

z1=LNZ(KIa,KIIa,rpa,dN,dn,sigam,oys,jj);




function [X3,L5,L6,Pi,taui,Pa,t]=refs(a,b,dt,Ti,sw1,sa,sw2)

t=[0:dt:6]';
q=1.028;
%冻结温度
%Ti=-0;
betaT=0.09-0.045*Ti/22;
si=betaT.*(1-exp(-q*t));
ET=66*(1-0.012*Ti);

%冻结压力
Pi=ET.*si+0.001;
%切向摩擦力
oi=72.86-0.2476*(273.15+Ti);%冰的单轴压缩强度 MPa
or=50;%岩石的单轴压缩强度 MPa
JRC=0.1;
Rs=sqrt(0.0001^2+0.003^2)/0.003;
taus=(0.52+0.089*sqrt(abs(Ti) ))*Rs;

tauP=Pi.*tan( 0.04*exp(1)^ (4.689*( 1-(273.15+Ti)/273.15 ) ) + JRC .*    log10( (oi^2+or^2)./ ( (oi+or).*Pi )  ) )  ;
taui=taus+tauP;
% 切向力作用最小段长度
n=size(Pi,1);
Pa=1.01;
A1=Pi(5,1);
A2=taui(n,1);
lmin=sqrt( (A1-Pa)^2*a^2*b^2 /( a^2*A2^2 + (A1-Pa)^2*b^2      ))       ;
%sa 气体的占比
%sa=0.2;




L=a-lmin;
%冻结前 上部分水分的体积
L1=L*sw1;
%冻结前 气体的体积长度
L2=L*sa;
%冻结前 下部分水分的体积
L3=L*sw2;
%冻结前 气体的体积长度
X1=L1+lmin; 
X2=L1+lmin+L2; 
f = @(x) 2*b/a*sqrt(a^2-x.^2);  % 函数句柄
Va1 = integral(f, X1, X2);  % 对 f(x) 从 a 到 b 积分，冻结前 气体的体积长度

%上部分水分的 冰压力、切向力作用长度 
L4=L1.*(1+si);
%下部分水分的 冰压力、切向力作用长度 
L6=L3.*(1+si);
%气体压力z作用长度、
L5=a-(L4+lmin)-L6;

%上部分水分的 冰压力、切向力作用坐标终点
X3=(L4+lmin); %冰压力、切向力作用长度 
%气体压力作用坐标终点
X4=a-L6;
%下部分水分的 冰压力、切向力作用坐标终点
X5=a; 


for i=1:1:n
Va2(i,1)= integral(f, X3(i,1), X4(i,1));
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
Pa=ref3./Va2;


end

function [KI,KII,Ki,Ka,Kt,s4]=KK(a,n,L4,L5,L6,Pi,taui,Pa)

% 冻结压力作用段长度
% Li=0.1; %单位mm
% % 气体压力作用段长度
% La=0.1; %单位mm
% % 切向力作用段长度
% Lt=0.1; %单位mm

%起点坐标
s1=0;
%上部分冰压力、摩擦力作用终点坐标
s2=L4./a;
%气体压力作用终点坐标
s3=(L4+L5)./a;
%下部分冰压力、摩擦力作用终点坐标
s4=(L4+L5+L6)./a;

%权函数参数
b1=2;
b2=0.980246;
b3=1.105334;
b4=-0.314547;
b5=-0.103356;

for j=1:1:n
    %上部分冰压力、摩擦力作用无量纲应力因子
f1=sqrt(2)/pi*( (1*b1+ 1/3*b2 + 1/5*b3  +1/7*b4 +1/9*b5 )-...
    ( 1*b1*(1-s2(j,1))^(1/2) + 1/3*b2*(1-s2(j,1))^(3/2) + 1/5*b3*(1-s2(j,1))^(5/2)  +1/7*b4*(1-s2(j,1))^(7/2) +1/9*b5*(1-s2(j,1))^(9/2)        )   );

%气体压力作用无量纲应力因子
f2=sqrt(2)/pi*( ( 1*b1*(1-s2(j,1))^(1/2) + 1/3*b2*(1-s2(j,1))^(3/2) + 1/5*b3*(1-s2(j,1))^(5/2)  +1/7*b4*(1-s2(j,1))^(7/2) +1/9*b5*(1-s2(j,1))^(9/2) )-...
    ( 1*b1*(1-s3(j,1))^(1/2) + 1/3*b2*(1-s3(j,1))^(3/2) + 1/5*b3*(1-s3(j,1))^(5/2)  +1/7*b4*(1-s3(j,1))^(7/2) +1/9*b5*(1-s3(j,1))^(9/2)        )   );
 %下部分冰压力、摩擦力作用无量纲应力因子
f3=sqrt(2)/pi*( ( 1*b1*(1-s3(j,1))^(1/2) + 1/3*b2*(1-s3(j,1))^(3/2) + 1/5*b3*(1-s3(j,1))^(5/2)  +1/7*b4*(1-s3(j,1))^(7/2) +1/9*b5*(1-s3(j,1))^(9/2)        ) -...
    ( 1*b1*(1-s4(j,1))^(1/2) + 1/3*b2*(1-s4(j,1))^(3/2) + 1/5*b3*(1-s4(j,1))^(5/2)  +1/7*b4*(1-s4(j,1))^(7/2) +1/9*b5*(1-s4(j,1))^(9/2)        )   );
end

%冻结压力应力因子
Ki=Pi.*f1.*sqrt(pi*a)+Pi.*f3.*sqrt(pi*a);
%气体压力应力因子
Ka=Pa.*f2.*sqrt(pi*a);

%切向压力应力因子
Kt=taui.*f1.*sqrt(pi*a)+taui.*f3.*sqrt(pi*a);

KI=Ki+Ka;
KII=Kt;

end
function [Rp]=rp(KI,KII,oys)
sigam0=[0:1:360]'.*(pi/180);
NN=size(sigam0,1); 
for i=1:1:NN
    sigam=sigam0(i,1);
A1=6*KI^2/(2*pi*oys^2) *cos(sigam/2)*cos(sigam/2) *sin(sigam/2)*sin(sigam/2)*...
    sin(3*sigam/2)*sin(3*sigam/2);
A21=cos(sigam/2)*sin(sigam/2)*( 1+cos(sigam/2)*cos(3*sigam/2)  ) ;
A22=cos(sigam/2)*cos(sigam/2)*( 1-sin(sigam/2)*sin(3*sigam/2)  ) ;

A2=3*KI*KII/(2*pi*oys^2)*sin(sigam/2)*sin(3*sigam/2)*...
    ( A21 + A22          );

A31=sin(sigam/2)*sin(sigam/2)*( 1+cos(sigam/2)*cos(3*sigam/2) )^2;
A32=cos(sigam/2)*cos(sigam/2)*( 1-sin(sigam/2)*sin(3*sigam/2)  )^2;

A3=3*KII^2/(2*pi*oys^2)*(A31+A32 );
A4=KI^2/(2*pi*oys^2) *cos(sigam/2)*cos(sigam/2)-KI*KII/(pi*oys^2)*cos(sigam/2)*sin(sigam/2)+...
    KI*KII/(pi*oys^2)*sin(sigam/2)*sin(sigam/2);
Rp(i,1)=A1+A2+A3+A4;
end

end

function [Z1]=LNZ(K1a,K2a,rpa,dN,dn,sigam,oys,jj)

KI=K1a;  
KII=K2a;
% RR1 ?右端点 r距离
[NN,~]=size(rpa);
RR1=zeros(NN,dN);
RR1(:,1)=rpa;
R1=zeros(NN,dn);
R1(:,1)=rpa;

% dn=300;

% oys=2e8;
dnn=dn/2;
a=oys/(149*2);
b=oys/2-a;
x=[1:1:dnn];
oy(:,1:dnn)=a*x+b;
oy(:,dnn+1:dn)=oys;
clear oys
for i=1:1:NN
A=RR1(i,1);
rend=30*A;
Rend=0;

dr=(rend-A)/(dN-1);
RR1(i,:)=(RR1(i,1):dr:rend)  ;

dR=(A-Rend)/(dn-1);
R1(i,:)=(R1(i,1):-dR:Rend)  ;
oys(i,:)=oy;

end

 for jkL=1:1:dN
for i=1:1:NN   
    o1=or(K1a(jj,1),K2a(jj,1),RR1(:,jkL)); 
 end
o_R(:,jkL)=o1;
 end
 
 
 
 for i=1:1:NN
 [X(i,:),Y(i,:),Z(i,:)]=pol2cart(sigam(i,1),RR1(i,:),o_R(i,:)) ;
 [X1(i,:),Y1(i,:),Z1(i,:)]=pol2cart(sigam(i,1),R1(i,:),oys(i,:)) ;
  end
 
 for i=2:1:NN
 AXX(1:dN,:)=X(1,:)';
 AXX((i-1)*dN+1:i*dN,:)=X(i,:)'; 
 BXX(1:dN,:)=X1(1,:)';
 BXX((i-1)*dN+1:i*dN,:)=X1(i,:)'; 

 
  AYY(1:dN,:)=Y(1,:)';
  AYY((i-1)*dN+1:i*dN,:)=Y(i,:)'; 
  BYY(1:dN,:)=Y1(1,:)';
  BYY((i-1)*dN+1:i*dN,:)=Y1(i,:)'; 
  
  
  
  AZZ(1:dN,:)=Z(1,:)';
  AZZ((i-1)*dN+1:i*dN,:)=Z(i,:)'; 
  BZZ(1:dN,:)=Z1(1,:)';
  BZZ((i-1)*dN+1:i*dN,:)=Z1(1,:)'; 
 
 end
 XX=[AXX;BXX];
 YY=[AYY;BYY];
 ZZ=[AZZ;BZZ];

Z1=[XX YY ZZ];
end

function o_sigam=or(KI,KII,rp)
sigam0=[0:1:360]'.*(pi/180);
NN=size(sigam0,1); 
[m,~]=size(rp);
for i=1:1:NN
  o=sigam0(i,1);
if m>1
    rc=rp(i,1);
elseif m==1
    rc=rp;
end
o_sigam(i,1)=KI/( 2*sqrt(2*pi*rc) )*cos(o/2)*(1+ cos(o) )-3*KII/( 2*sqrt(2*pi*rc) )*cos(o/2)*sin(o);
end
end



