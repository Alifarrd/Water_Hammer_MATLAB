function[pend,vend,vi,pi,ui,zi]=freevalve(o1,o2,r_f,ge,as,ii,tx,vi,pi,zi,ui,ox,a,b,c,a2,b2,c2,dxz,cth,dt,az,hres)
pi(ii,tx+1)=(((1/as)*vi(ii-1,tx))+((a*r_f*ge)*pi(ii-1,tx))+(b*ui(ii-1,tx))-(c*zi(ii-1,tx))+(c*r_f*ge*(as/az)*hres))/((r_f*ge)*(a-(c*(as/az))));

zi(ii,tx+1)=((r_f*ge)*(pi(ii,tx+1)-hres)*(as/az));
%=========================================================================[3]
y1=((1/as)*vi(ii-1,tx))+((a*r_f*ge)*pi(ii-1,tx))+(b*ui(ii-1,tx))-(c*zi(ii-1,tx))-(c*r_f*ge*(as/az)*hres);
%gg1=((1/as)*vi(ii-1,tx))+((a*r_f*ge)*pi(ii-1,tx))+(b*ui(ii-1,tx))-(c*zi(ii-1,tx))

o1=abs((5*dxz)-(cth*dt));
o2=dxz-o1;

%!-----------------------------------------[opproximation of V,P,U,Z]   
   ox(1)=((o1*vi(ii-4,tx))+(o2*vi(ii-5,tx)))/dxz;
   ox(2)=((o1*pi(ii-4,tx))+(o2*pi(ii-5,tx)))/dxz;
   ox(3)=((o1*ui(ii-4,tx))+(o2*ui(ii-5,tx)))/dxz;
   ox(4)=((o1*zi(ii-4,tx))+(o2*zi(ii-5,tx)))/dxz;
%!-----------------------------------------[opproximation of V,P,U,Z]
 y2=((-b2/as)*ox(1))-((c2*r_f*ge)*ox(2))+(ox(3))-(a2*ox(4))-(a2*r_f*ge*(as/az)*hres);
%gg2=((-b2/as)*ox(1))-((c2*r_f*ge)*ox(2))+(ox(3))-(a2*ox(4))

 l1=(1+b)/as;
 l2=(r_f*ge)*(a-(c*(as/az)));
 l3=(1-b2)/as;
 l4=(r_f*ge)*(-c2-(a2*(as/az)));
% 
 ff=(l1*l4)-(l2*l3);
% 
vend=((l4*y1)+(-l2*y2))/ff;
pend=((-l3*y1)+(l1*y2))/ff;


vi(ii,tx+1)=vend;
pi(ii,tx+1)=pend;

ui(ii,tx+1)=vi(ii,tx+1)/as;
zi(ii,tx+1)=((r_f*ge)*(pi(ii,tx+1)-hres)*(as/az));

