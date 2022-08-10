function[xg]=guesswaterhammer(t,tr,ng);

for kx=1:4
for    px=1:4
    tnew(kx,px)=t(kx,px);
end
end

for kx=1:4
    tnew(kx,5)=tr(kx);
end


for kg=1:ng
ag(kg,:)=tnew(kg,:);
end 

for ig=1:ng
for jg1=ig+1:ng
landag=-ag(jg1,ig)/ag(ig,ig);
ag(jg1,:)=(landag*ag(ig,:))+ag(jg1,:);
end 
end 
 xg(ng)=ag(ng,ng+1)/ag(ng,ng);
 for ig3=ng-1:-1:1
 sg=0.;
 for jg3=ig3+1:ng
 sg=sg+(ag(ig3,jg3)*xg(jg3));
 end 
 xg(ig3)=(ag(ig3,ng+1)-sg)/ag(ig3,ig3);
 end
end



