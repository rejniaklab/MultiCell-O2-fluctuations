function script_blue_influx

oxy_init=load('Exp_data/blue_exp.txt');
minLimit=oxy_init(1,2)-3.5;   % lower bound for initial oxy saturation
maxLimit=oxy_init(1,2)+3.5;   % upper bound for initial oxy saturation


tissueIn=['tissue'];    % dir with tissue data  
tissueOut=['tiss_fluc_blue'];% dir to save fluctuations


file_to_save=['data_blue_influx.txt'];


if exist(file_to_save)
  save_data=load(file_to_save);
  Nsave=size(save_data,1);
else
  save_data=zeros(1,7); %[V,T,S,oxy,seed,r2,norm2]
  Nsave=0;
end


data_files=load([tissueIn,'/data_summary.txt']);
ind=find((data_files(:,4)>=minLimit)&(data_files(:,4)<=maxLimit));
Nfiles=size(ind,1); 


for ii=1:Nfiles
  Vfrac=data_files(ind(ii),1)/100; 
  Tfrac=data_files(ind(ii),2)/100; 
  Sfrac=data_files(ind(ii),3)/100; 
  
  [r2,norm2]=oxy_tissue_fluc_blue_script(Vfrac,Tfrac,Sfrac,tissueIn,...
             tissueOut,1);

   if (r2>0)||(norm2<100)
     Nsave=Nsave+1;
     save_data(Nsave,1:5)=data_files(ind(ii),1:5);
     save_data(Nsave,6:7)=[r2,norm2];   
     save([file_to_save],'save_data','-ascii')
   end
   close all
end

save([file_to_save],'save_data','-ascii')

end