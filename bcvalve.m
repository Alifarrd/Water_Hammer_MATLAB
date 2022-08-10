function[m,ansf,vi,pi,ui,zi]=bcvalve(o1,o2,r_f,ge,as,ii,az,tx,vi,pi,zi,ui,ox,a,b,c,a2,b2,c2,dxz,hres)




m=zeros(4,4);
n=zeros(4,1);
ansf=zeros(4,1);


%=======================
n(1)=((1/as)*vi(ii-1,tx))+((a*r_f*ge)*pi(ii-1,tx))+(b*ui(ii-1,tx))-(c*zi(ii-1,tx));
%!-----------------------------------------[opproximation of V,P,U,Z]   
   ox(1)=((o1*vi(ii-4,tx))+(o2*vi(ii-5,tx)))/dxz;
   ox(2)=((o1*pi(ii-4,tx))+(o2*pi(ii-5,tx)))/dxz;
   ox(3)=((o1*ui(ii-4,tx))+(o2*ui(ii-5,tx)))/dxz;
   ox(4)=((o1*zi(ii-4,tx))+(o2*zi(ii-5,tx)))/dxz;
%!-----------------------------------------[opproximation of V,P,U,Z]
n(2)=((-b2/as)*ox(1))-((c2*r_f*ge)*ox(2))+(ox(3))-(a2*ox(4));
%=============================================================================================================
%============================

m(1,1)=1/as;
m(1,2)=a*r_f*ge;
m(1,3)=b;
m(1,4)=-c;

m(2,1)=-b2/as;
m(2,2)=-c2*r_f*ge;
m(2,3)=1;
m(2,4)=-a2;

m(3,1)=1;
m(3,2)=0;
m(3,3)=-as;
m(3,4)=0;

m(4,1)=0;
m(4,2)=-(r_f*ge)*(as/az);
m(4,3)=0;
m(4,4)=1;

n(3)=0;
n(4)=0;

%=======================
ansf=m\n;
%==========================================================================

vi(ii,tx+1)=ansf(1);
pi(ii,tx+1)=ansf(2);
ui(ii,tx+1)=ansf(3);
zi(ii,tx+1)=ansf(4);

