clc;
clear;
%==========================================================================
tt=0.2;                              %!Total time
%1=fixed
%2=Free
BCv=2;
%=====================================================[Initial Information]
ge=9.81;%9.806;           %!gravitational acceleration(m/s2)
k=2.1e9;            %!Fluid's Bulk modulus,[pa]
vp=0.3; %!0.5!0.29  %!Possion ratio,
r=(0.797/2);%(0.0506/2);%(0.797/2);         %!inner Radius of pipe,[m]
e=2.11D11 %1.43e9;             %!Young modulus of pipe wall material,[pa]
ee=0.008 %0.0063;           %!pipe wall thickness,[m]
r_f=1000;           %!fluid mass density,[kg/m3]
r_t=7870.0D0 %5000.;          %!Tube mass density,[kg/m3]
pres=0;%1.e5;          %!reservoir pressure
hres=0.0;%pres/(r_f*ge);
v0=1.002 %0.5;           %!initial velosity
%!========================================================================= 
ks=1/((1/k)+((1-vp^2)*((2*r)/(e*ee))));                          %!effictive bulk modulus,[pa]
cf=(ks/r_f)^0.5;                                                 %!Classical Fluid wave speed,[m/s]
ct=(e/r_t)^0.5;                                                  %!Classical Tube wave speed,[m/s]
g=(((1+((2*vp^2)*(r_f/r_t)*(r/ee)))*cf^2)+ct^2)^0.5;             %!constant,[m/s]
cfh=0.5*sqrt(2.)*(g^2-(g^4-(4*cf^2)*ct^2)^0.5)^0.5            %!FSI fluid wave speed,[m/s]
cth=0.5*sqrt(2.)*(g^2+(g^4-(4*cf^2)*ct^2)^0.5)^0.5            %!FSI Tube wave speed,[m/s]
%cfh=1024.7;
%cth=5280.5;
%rc=cth/cfh
ddd3=cth/cfh
%!=========================================================================
a=1/(r_f*cfh);
b=(2*vp)*(((cfh/ct)^2)/(1-(cfh/ct)^2));
c=((2*vp)/(r_t*cfh))*(((cfh/ct)^2)/(1-(cfh/ct)^2));
%!=========================================================================
a2=1/(r_t*cth);
b2=((-vp*r*r_f)/(e*r_t))*(((cf/cth)^2)/(1-(cf/cth)^2));
c2=((vp*r)/(e*r_t*cth))*(((cf/cth)^2)/(1-(cf/cth)^2));
%!=========================================================================
as=(3.14/4)*(2*r)^2;
at=(3.14/4)*(2*(r-ee))^2;
az=as-at;
q0=as*v0;
%!=========================================================================
l=20;                                %!length
ne=20;                               %!number of element
%tt=0.1;                              %!Total time
nx=ne+1;                             %!number of node
dxz=l/ne;                            
dt=dxz/cfh
nt=floor(tt/dt)+1;                   %!number of time steps
%!=========================================================================
t=zeros(4,4);

tr=zeros(4,1);

vi(:,:)=zeros(nx,nt+1);
pi(:,:)=zeros(nx,nt+1);
ui(:,:)=zeros(nx,nt+1);
zi(:,:)=zeros(nx,nt+1);

%!=========================================================================
t(1,1)=1/as;
t(1,2)=a*r_f*ge;
t(1,3)=b;
t(1,4)=-c;
%!----------
t(2,1)=1/as;
t(2,2)=-a*r_f*ge;
t(2,3)=b;
t(2,4)=c;
%!----------
t(3,1)=-b2/as;
t(3,2)=-c2*r_f*ge;
t(3,3)=1;
t(3,4)=-a2;
%!----------
t(4,1)=-b2/as;
t(4,2)=c2*(r_f*ge);
t(4,3)=1;
t(4,4)=a2;
%!=========================================================================[initial boundry]
%hres=0;
vi(:,:)=0;
pi(:,:)=0;
pi(1,:)=hres;
ui(:,:)=0;
zi(:,:)=0;

vi(:,1)=q0;
vi(nx,:)=0;


pi(:,1)=hres;
pi(nx,1)=0;%hres;

ui(:,1)=0;
zi(:,1)=hres*(as/az);

%!=========================================================================[initial boundry]
  for tx=1:(nt-1)
%!=========================================================================[boundry Condition]
%!=========================================================================[reservoir]

[vir,zir,o1,o2,ox,vi,zi]=fixedR(r_f,ge,hres,as,tx,vi,pi,zi,ui,dxz,a,b,c,a2,b2,c2,cth,dt);

%!=========================================================================[valve]

[vi,pi,ui,zi]=boundryCondition(BCv,o1,o2,r_f,ge,as,tx,vi,pi,zi,ui,ox,a,b,c,a2,b2,c2,dxz,cth,dt,az,hres,nx);

%!=========================================================================[boundry Condition] 
for i=2:nx-1
    
if i<=5            %!......................................................FIRST NODEs:
  zt=dt-(((i-1)*dxz)/cth)
  zz=dt-zt
%%!-----------------------------------------[opproximation of V,P,U,Z]   
  gx(1)=((zz*vi(1,tx))+(zt*vi(1,tx+1)))/dt;
  gx(2)=((zz*pi(1,tx))+(zt*pi(1,tx+1)))/dt;
  gx(3)=((zz*ui(1,tx))+(zt*ui(1,tx+1)))/dt;
  gx(4)=((zz*zi(1,tx))+(zt*zi(1,tx+1)))/dt;
%%!-----------------------------------------[opproximation of V,P,U,Z]
 o1=abs((5*dxz)-(cth*dt))
 o2=dxz-o1

%%!-----------------------------------------[opproximation of V,P,U,Z]   
  ox(1)=((o1*vi(i+4,tx))+(o2*vi(i+5,tx)))/dxz;
  ox(2)=((o1*pi(i+4,tx))+(o2*pi(i+5,tx)))/dxz;
  ox(3)=((o1*ui(i+4,tx))+(o2*ui(i+5,tx)))/dxz;
  ox(4)=((o1*zi(i+4,tx))+(o2*zi(i+5,tx)))/dxz;
%%!-----------------------------------------[opproximation of V,P,U,Z]
 
tr(3)=((-b2/as)*gx(1))-((c2*r_f*ge)*gx(2))+(gx(3))-(a2*gx(4));
tr(4)=((-b2/as)*ox(1))+((c2*r_f*ge)*ox(2))+(ox(3))+(a2*ox(3));

elseif i>=nx-5      %!.....................................................END NODES:
  zt=dt-(((nx-(i-1))*dxz)/cth);
  zz=dt-zt;
%%!-----------------------------------------[opproximation of V,P,U,Z]   
  gx(1)=((zz*vi(nx,tx))+(zt*vi(nx,tx+1)))/dt;
  gx(2)=((zz*pi(nx,tx))+(zt*pi(nx,tx+1)))/dt;
  gx(3)=((zz*ui(nx,tx))+(zt*ui(nx,tx+1)))/dt;
  gx(4)=((zz*zi(nx,tx))+(zt*zi(nx,tx+1)))/dt;
%%!-----------------------------------------[opproximation of V,P,U,Z]
 o1=abs((5*dxz)-(cth*dt));
 o2=dxz-o1;
%%!-----------------------------------------[opproximation of V,P,U,Z]   
  ox(1)=((o1*vi(i-4,tx))+(o2*vi(i-5,tx)))/dxz;
  ox(2)=((o1*pi(i-4,tx))+(o2*pi(i-5,tx)))/dxz;
  ox(3)=((o1*ui(i-4,tx))+(o2*ui(i-5,tx)))/dxz;
  ox(4)=((o1*zi(i-4,tx))+(o2*zi(i-5,tx)))/dxz;
%%!-----------------------------------------[opproximation of V,P,U,Z]
 
tr(3)=((-b2/as)*ox(1))-((c2*r_f*ge)*ox(2))+(ox(3))-(a2*ox(4));
tr(4)=((-b2/as)*gx(1))+((c2*r_f*ge)*gx(2))+(gx(3))+(a2*gx(4));

else                     %!................................................MIDDEL NODES

o1=abs((5*dxz)-(cth*dt));
o2=dxz-o1;

%%!-----------------------------------------[opproximation of V,P,U,Z]   
  gx(1)=((o1*vi(i+4,tx))+(o2*vi(i+5,tx)))/dxz;
  gx(2)=((o1*pi(i+4,tx))+(o2*pi(i+5,tx)))/dxz;
  gx(3)=((o1*ui(i+4,tx))+(o2*ui(i+5,tx)))/dxz;
  gx(4)=((o1*zi(i+4,tx))+(o2*zi(i+5,tx)))/dxz;
%%!-----------------------------------------[opproximation of V,P,U,Z]
  %o1=(5*dxz)-(cth*dt);
  %o2=dxz-o1;
%%!-----------------------------------------[opproximation of V,P,U,Z]   
  ox(1)=((o1*vi(i-4,tx))+(o2*vi(i-5,tx)))/dxz;
  ox(2)=((o1*pi(i-4,tx))+(o2*pi(i-5,tx)))/dxz;
  ox(3)=((o1*ui(i-4,tx))+(o2*ui(i-5,tx)))/dxz;
  ox(4)=((o1*zi(i-4,tx))+(o2*zi(i-5,tx)))/dxz;
%%!-----------------------------------------[opproximation of V,P,U,Z]

tr(3)=((-b2/as)*ox(1))-((c2*r_f*ge)*ox(2))+(ox(3))-(a2*ox(4));
tr(4)=((-b2/as)*gx(1))+((c2*r_f*ge)*gx(2))+(gx(3))+(a2*gx(4));

end 

tr(1)=((1/as)*vi(i-1,tx))+((a*r_f*ge)*pi(i-1,tx))+((b)*ui(i-1,tx))-((c)*zi(i-1,tx));
tr(2)=((1/as)*vi(i+1,tx))-((a*r_f*ge)*pi(i+1,tx))+((b)*ui(i+1,tx))+((c)*zi(i+1,tx));
%%!--------------------------------------------------[guess]
ng=4;
[xg]=guesswaterhammer(t,tr,ng);
%%!--------------------------------------------------[guess]

vi(i,tx+1)=xg(1);
pi(i,tx+1)=xg(2);
ui(i,tx+1)=xg(3);
zi(i,tx+1)=xg(4);

end

  end
 time=0:dt:dt*(nt);
[vi,pi,ui,zi]=print(vi,pi,zi,ui,as,time);
[vi,pi,ui,zi]=plotali(vi,pi,zi,ui,BCv,ne,time);
  