function script_influx

  clear all
  close all
  
 
  global pathOutTiss textExp Vfrac Tfrac Sfrac
  
  pathOut='tissue_fluctuations';  % create a new general directory
  if ~exist(pathOut, 'dir') mkdir(pathOut); end


  for num_exp=1:4      % select the tissue for each ROI
    if     num_exp==1  Vfrac=4;   Tfrac=15; Sfrac=75; textExp='black';
    elseif num_exp==2  Vfrac=1.5; Tfrac=15; Sfrac=60; textExp='red';
    elseif num_exp==3  Vfrac=1.5; Tfrac=40; Sfrac=45; textExp='blue';
    elseif num_exp==4  Vfrac=0.5; Tfrac=20; Sfrac=75; textExp='magenta';
    end
    pathOutTiss=[pathOut,'/data_V',num2str(Vfrac),'_T',num2str(Tfrac),...
                 '_S',num2str(Sfrac)]   % create new directory
    if ~exist(pathOutTiss, 'dir') mkdir(pathOutTiss); end        
             
    
    
    % determine the optimal schedule for influx for each segment
    Main_fluct_optim_influx(Vfrac,Tfrac,Sfrac,textExp,pathOutTiss)
     
    
    
    % run the determined optimal schedule for all segments 
    [r2,norm2]=oxy_tissue_fluct_influx(Vfrac,Tfrac,Sfrac,'tissue',...
               pathOutTiss,num_exp);
    
           
    % save final data including [Vfrac,Tfrac,Sfrac,stable oxy,seed,r2,norm2 
    % read seed and stable oxy form the saved data
      pathIn=['tissue/data_V',num2str(Vfrac),'_T',num2str(Tfrac),...
              '_S',num2str(Sfrac)]; 
      param=load([pathIn,'/parameters.txt']);
      oxy=load([pathIn,'/Data/oxy.txt']);  
      
    % final data to save  
      save_data(1,1:7)=[Vfrac,Tfrac,Sfrac,mean(mean(oxy)),param(1),r2,norm2];       
      save([pathOutTiss,'/final_',textExp,'_influx.txt'],'save_data','-ascii')

    
      close all
    
  end % for
end % function

