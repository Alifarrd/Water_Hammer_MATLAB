function [vir,zir,o1,o2,ox,vi,zi]=fixedR(r_f,ge,hres,as,tx,vi,pi,zi,ui,dxz,a,b,c,a2,b2,c2,cth,dt)

x1=((1/as)*vi(2,tx))-((a*r_f*ge)*pi(2,tx))+(a*r_f*ge*hres)+(b*ui(2,tx))+(c*zi(2,tx));
  o1=(5*dxz)-(cth*dt);
  o2=dxz-o1;
% !-----------------------------------------[opproximation of V,P,U,Z]   
  ox(1)=((o1*vi(5,tx))+(o2*vi(6,tx)))/dxz;
  ox(2)=((o1*pi(5,tx))+(o2*pi(6,tx)))/dxz;
  ox(3)=((o1*ui(5,tx))+(o2*ui(6,tx)))/dxz;
  ox(4)=((o1*zi(5,tx))+(o2*zi(6,tx)))/dxz;
%!-----------------------------------------[opproximation of V,P,U,Z]
x2=((-b2/as)*ox(1))-((c2*r_f*ge)*ox(2))-(c2*r_f*ge*hres)+(ox(3))+(a2*ox(4));


l1=1/as;
l2=c;
l3=-b2/as;
l4=a2;

ff=(l1*l4)-(l2*l3);

vir=((l4*x1)+(-l2*x2))/ff;
zir=((-l3*x1)+(l1*x2))/ff;

vi(1,tx+1)=vir;
zi(1,tx+1)=zir;
return


