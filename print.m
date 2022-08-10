function[vi,pi,ui,zi]=print(vi,pi,zi,ui,as,time)
% % xlswrite(filename,A,sheet,xlRange) 
 xlswrite('C:\Users\Ali Fard\Desktop\ans-WH_FSI4.xlsx',vi/as,'v');
 xlswrite('C:\Users\Ali Fard\Desktop\ans-WH_FSI4.xlsx',pi,'h');
 xlswrite('C:\Users\Ali Fard\Desktop\ans-WH_FSI4.xlsx',ui,'udot');
 xlswrite('C:\Users\Ali Fard\Desktop\ans-WH_FSI4.xlsx',zi,'zigma');
xlswrite('C:\Users\Ali Fard\Desktop\ans-WH_FSI4.xlsx',time,'time');
