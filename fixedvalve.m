function[pend,zend,pi,zi]=fixedvalve(o1,o2,r_f,ge,as,ii,tx,vi,pi,zi,ui,ox,a,b,c,a2,b2,c2,dxz,cth,dt,hres,az)

vdot=0;
udot=0;

y1=((1/as)*vi(ii-1,tx))+((a*r_f*ge)*pi(ii-1,tx))+(b*ui(ii-1,tx))-(c*zi(ii-1,tx))-((1/as)*vdot)-(b*udot);

  o1=(5*dxz)-(cth*dt);
  o2=dxz-o1;
%!-----------------------------------------[opproximation of V,P,U,Z]   
  ox(1)=((o1*vi(ii-4,tx))+(o2*vi(ii-5,tx)))/dxz;
  ox(2)=((o1*pi(ii-4,tx))+(o2*pi(ii-5,tx)))/dxz;
  ox(3)=((o1*ui(ii-4,tx))+(o2*ui(ii-5,tx)))/dxz;
  ox(4)=((o1*zi(ii-4,tx))+(o2*zi(ii-5,tx)))/dxz;
%!-----------------------------------------[opproximation of V,P,U,Z]
y2=((-b2/as)*ox(1))-((c2*r_f*ge)*ox(2))+(ox(3))-(a2*ox(4))-((-b2/as)*vdot)-(udot);

l1=a*r_f*ge;
l2=-c;
l3=-c2*r_f*ge;
l4=-a2;

ff=(l1*l4)-(l2*l3);

pend=((l4*y1)+(-l2*y2))/ff;
zend=((-l3*y1)+(l1*y2))/ff;
zend=((r_f*ge)*(pend-hres)*(as/az));

 pi(ii,tx+1)=pend;
 zi(ii,tx+1)=zend;

