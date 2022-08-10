function[vi,pi,ui,zi]=boundryCondition(BCv,o1,o2,r_f,ge,as,tx,vi,pi,zi,ui,ox,a,b,c,a2,b2,c2,dxz,cth,dt,az,hres,nx)
ii=nx;
%==========================================================================[1]fixedvalve
if BCv==1 
[pend,zend,pi,zi]=fixedvalve(o1,o2,r_f,ge,as,ii,tx,vi,pi,zi,ui,ox,a,b,c,a2,b2,c2,dxz,cth,dt,hres,az);

%==========================================================================[2]freevalve
elseif BCv==2
[pend,vend,vi,pi,ui,zi]=freevalve(o1,o2,r_f,ge,as,ii,tx,vi,pi,zi,ui,ox,a,b,c,a2,b2,c2,dxz,cth,dt,az,hres);

%==========================================================================[free valve dastgah]
    else 
[m,ansf,vi,pi,ui,zi]=bcvalve(o1,o2,r_f,ge,as,ii,az,tx,vi,pi,zi,ui,ox,a,b,c,a2,b2,c2,dxz,hres);

%==========================================================================
end
    
return




