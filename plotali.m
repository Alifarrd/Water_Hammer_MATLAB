function[vi,pi,ui,zi]=plotali(vi,pi,zi,ui,BCv,ne,time)
hold on;
nx=floor(ne);
if(BCv==1)
plot(time,pi(nx,:),'b')% 'DisplayName', strcat('Head at the valve') , 'YDataSource', 'h(nx,1:nt)'); figure(gcf)    
legend('Fixed')
else
plot(time,pi(nx,:),'r')% 'DisplayName', strcat('Head at the valve') , 'YDataSource', 'h(nx,1:nt)'); figure(gcf)    
legend('Free')
end