

if exist('data_summary.txt')
  save_data=load('data_summary.txt');
  Nsave=size(save_data,1);
else
  save_data=zeros(1,5); 
  Nsave=0;
end

for Vnum=0.5:0.5:5
  for Tnum=10:5:95
    for Snum=5:5:90
      Sfrac=Snum/100;
      Vfrac=Vnum/100;
      Tfrac=Tnum/100;
      dat=[Vnum,Tnum,Snum]
      [aver,seed]=oxy_tissue_make_aver(Vfrac,Tfrac,Sfrac);
    
      if (aver>0)
        Nsave=Nsave+1;
        save_data(Nsave,1:5)=[Vnum,Tnum,Snum,aver,seed];
        save('data_summary.txt','save_data','-ascii')
      end
      close all
    end
  end
end
%save('data_summary.txt','save_data','-ascii')

