function [r2,norm_res]=oxy_tissue_fluc_magenta_script(Vfrac,Tfrac,Sfrac,tiss_pathIn,...
tiss_pathOut,case_fluct)
% Vfrac -- fraction of tissue that is vasculature, i.e. 0.05 means 5%
% Tfrac -- fraction of tissue that is tumor, i.e. 0.2 means 20%
% Sfrac -- fraction of tissue that is stromal cells, i.e. 0.1 means 10%

% tiss_path -- name of the directory to save data
% case_fluct -- 1: influcx only; 2: uptake only; 3: both
r2=0; norm_res=100;



dt=0.05;        % time step [sec]
Niter=28*60/dt; % 30 minutes
Ntime=60/dt;    %  1 minute interval 


toSave=1;         % switch =1 will save data ==0 no data save
toSavelast=1;
Nsave=floor(60/dt);  % frequency of saving the data
Ndraw=floor(60/dt);  % frequency of drawing the data


pathDir=[tiss_pathOut,'/data_V',num2str(Vfrac*100),'_T',num2str(Tfrac*100),...
         '_S',num2str(Sfrac*100)];
pathData=[pathDir,'/DataFluc',num2str(case_fluct)]; % name of the directiory for saving data, 2 for influx, 3 for uptake
pathFigs=[pathDir,'/FigsFluc',num2str(case_fluct)]; % name of the directiory for saving figs    



%~~~~~~~~~~~~~
exp_oxy=load('Exp_data/magenta_exp.txt');

if case_fluct==1
  influx=load(['Opt_data/magenta/rateInflux_opt_influx.txt']);
  uptake=load(['Opt_data/magenta/rateUptake_opt_influx.txt']);
elseif case_fluct==2
  influx=load(['Opt_data/magenta/rateInflux_opt_uptake.txt']);
  uptake=load(['Opt_data/magenta/rateUptake_opt_uptake.txt']);
else
  influx=load(['Opt_data/magenta/rateInflux_opt_both.txt']);
  uptake=load(['Opt_data/magenta/rateUptake_opt_both.txt']);
end
%~~~~~~~~~~~~


% make output directory
if (toSave==1)  mkdir(pathDir); mkdir(pathData); mkdir(pathFigs); end
if (toSavelast==1)  mkdir(pathDir); mkdir(pathData); mkdir(pathFigs); end


% tissue not to dense
no_go=0;        
if (Vfrac+Tfrac+Sfrac>1) % if more than 100% do not simulate
  disp('more than 100% cellualrity & vascularity'); no_go=1;
end


% parameters
stiff=50;        % how hard the cells are pushing away from each other
nu=  250/10;     % viscosity of the medium
force_thr=5;     % threshold for stopping the forces
error_thr=1e-6;  % threshold for stopping the oxygen stability

diamVess   =40;        % vessel diameter [microns]
diamCell   =15;        % tumor cell diameter [microns]
diamStroma = 7.5;      % stromal cell diameter

% domain 
xmin=-500; xmax=-xmin; ymin=xmin; ymax=-xmin;  % 1000x1000 micron^2=1mm^2
hg=5;           % grid width [microns]


Doxy   =100;    % diffusion coefficient for oxygen [micron^2/s]
oxy_max= 60;    % level of oxygen in the vessel [mmHg=sigma-gram/micron^3]
oxy_init= .25; % initial oxygen level

oxy_upTum0  =9.54/(hg*hg);   % max rate of oxy uptake by tumor cells V
oxy_upStrom=9.54/(hg*hg);   % max rate of oxy uptake by stromal cells V
km         =1.134;          % Michealis-Menten coefficient, V*O/(km+O) 


cfl_cond=Doxy*dt/(hg*hg);   % stability condition
if (cfl_cond>0.25)
  disp('stability condition is not satisfied'); pause;
end


% to remember the average oxygen and maximal force
avg_oxy=zeros(1,1);
max_for=zeros(1,1);
Nav=0;   % index for the above vectors



if no_go==0    % if the tissue is not too dense 
  
  % domain 
  [xx,yy]=meshgrid(xmin:hg:xmax,ymin:hg:ymax); % grid for oxygen 
  xxg=xx(1,:);  yyg=yy(:,1);  Ngx=length(xxg);  Ngy=length(yyg);

  % create the figure
  figure('position',[750,500,1500,500])
  pause(0.1)

  
calc_oxy=zeros(9,1);

%read intital data from outside files
xyVess=load([tiss_pathIn,'/data_V',num2str(Vfrac*100),'_T',...
      num2str(Tfrac*100),'_S',num2str(Sfrac*100),'/Data/xyVess.txt']);
xyTum=load([tiss_pathIn,'/data_V',num2str(Vfrac*100),'_T',...
      num2str(Tfrac*100),'_S',num2str(Sfrac*100),'/Data/xyTum.txt']);
xyStrom=load([tiss_pathIn,'/data_V',num2str(Vfrac*100),'_T',...
      num2str(Tfrac*100),'_S',num2str(Sfrac*100),'/Data/xyStrom.txt']);
oxy=load([tiss_pathIn,'/data_V',num2str(Vfrac*100),'_T',...
      num2str(Tfrac*100),'_S',num2str(Sfrac*100),'/Data/oxy.txt']);
Nvess=size(xyVess,1);
Ntum=size(xyTum,1);
Nstroma=size(xyStrom,1);
vess_oxy=oxy_max*ones(Nvess,1);


% -------------- main code ---------------- %
  
  iter=-1;                            % initial iteration  
   
  % conditions for cell relocation
  force_ok=0;                                       
  
% ----- main loop ----- %  
  
 tic


 
iter=4*Ntime-1;
 % while (force_ok==0)||(err_oxy>error_thr) %force & error too large  
  while (iter<Niter)
    iter=iter+1;   % increase iteration
    
    
    if (iter<=4*Ntime) %300 incerements because the paper recorded fluc every 3 minutes
      if iter==4*Ntime
        calc_oxy(1,1)=mean(mean(oxy));
      end
      vess_oxy =influx(1)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(1)*oxy_upTum0;        
    elseif (iter<=7*Ntime)
      if iter==7*Ntime
        calc_oxy(2,1)=mean(mean(oxy));
      end
      vess_oxy =influx(2)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(2)*oxy_upTum0;        
    elseif (iter<=10*Ntime)
      if iter==10*Ntime
        calc_oxy(3,1)=mean(mean(oxy));
      end
      vess_oxy =influx(3)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(3)*oxy_upTum0;        
    elseif (iter<=13*Ntime)
      if iter==13*Ntime
        calc_oxy(4,1)=mean(mean(oxy));
      end
      vess_oxy =influx(4)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(4)*oxy_upTum0;       
    elseif (iter<=16*Ntime)
      if iter==16*Ntime
        calc_oxy(5,1)=mean(mean(oxy));
      end
      vess_oxy =influx(5)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(5)*oxy_upTum0;       
    elseif (iter<=19*Ntime)
      if iter==19*Ntime
        calc_oxy(6,1)=mean(mean(oxy));
      end
      vess_oxy =influx(6)*oxy_max*ones(Nvess,1);
      oxy_upTum =uptake(6)*oxy_upTum0;        
    elseif (iter<=22*Ntime)
      if iter==22*Ntime
        calc_oxy(7,1)=mean(mean(oxy));
      end
      vess_oxy =influx(7)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(7)*oxy_upTum0;        
    elseif (iter<=25*Ntime)
      if iter==25*Ntime
        calc_oxy(8,1)=mean(mean(oxy));
      end
      vess_oxy =influx(8)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(8)*oxy_upTum0;       
    elseif (iter<=28*Ntime)
      if iter==28*Ntime
        calc_oxy(9,1)=mean(mean(oxy));
      end
      vess_oxy =influx(9)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(9)*oxy_upTum0;        
    else
      vess_oxy =influx(10)*oxy_max*ones(Nvess,1); 
      oxy_upTum =uptake(10)*oxy_upTum0;       
    end
    
           
      
    % define influx of oxy from veins
    for ii=1:Nvess
      Nx=1+floor((xyVess(ii,1)-xmin)/hg);  % oxy grid closest to the vessel
      Ny=1+floor((xyVess(ii,2)-ymin)/hg);  
      for ix=-5:5                  % inspect near-by grid points
        for iy=-5:5                % to find this inside the vein
          if (Nx+ix>0)&&(Nx+ix<=Ngx)&&(Ny+iy>0)&&(Ny+iy<=Ngy)  
            ixy=sqrt((xyVess(ii,1)-(xmin+(Nx+ix)*hg))^2+...
                     (xyVess(ii,2)-(ymin+(Ny+iy)*hg))^2);
            if (ixy<0.5*diamVess)   % is inside the vein? 
              oxy(2+Nx+ix,2+Ny+iy)=vess_oxy(ii);
            end                    % influx rate
          end
        end
      end
    end
     
    % define diffusion of oxygen
    oxy(1:Ngx+2,1)      =oxy(1:Ngx+2,2);       % boundary conditions
    oxy(1:Ngx+2,Ngy+2)  =oxy(1:Ngx+2,Ngy+1);
    oxy(      1,1:Ngx+2)=oxy(      2,1:Ngx+2);
    oxy(  Ngx+2,1:Ngx+2)=oxy(  Ngx+1,1:Ngx+2);
    for ii=2:Ngx+1                             % diffusion loop 
      for jj=2:Ngy+1
        oxy(ii,jj)=oxy(ii,jj)+Doxy*(dt/hg^2)*(oxy(ii-1,jj)+oxy(ii+1,jj)+...
                   oxy(ii,jj-1)+oxy(ii,jj+1)-4*oxy(ii,jj));
      end
    end
    
    % define oxygen uptake by tumor cells
    for ii=1:Ntum
      Nx=1+floor((xyTum(ii,1)-xmin)/hg);  % grid closest to the tumor cell
      Ny=1+floor((xyTum(ii,2)-ymin)/hg);  
      for ix=-3:3                         % inspect near-by grid points
        for iy=-3:3                       % to find this inside the tumor
          if (Nx+ix>0)&&(Nx+ix<=Ngx)&&(Ny+iy>0)&&(Ny+iy<=Ngy)  
            ixy=sqrt((xyTum(ii,1)-(xmin+(Nx+ix)*hg))^2+...
                     (xyTum(ii,2)-(ymin+(Ny+iy)*hg))^2);
            if (ixy<0.5*diamCell)     
              oxy(2+Nx+ix,2+Ny+iy)=max(0,oxy(2+Nx+ix,2+Ny+iy)*...
                              (1-dt*oxy_upTum/(km+oxy(2+Nx+ix,2+Ny+iy))));
            end                 % M-M method V*o/(km+o)
          end
        end   
      end
    end
    
    %oxy uptake by stromal cells
    for ii=1:Nstroma
      Nx=1+floor((xyStrom(ii,1)-xmin)/hg);  % grid closest to the cell
      Ny=1+floor((xyStrom(ii,2)-ymin)/hg);  
      for ix=-3:3                          % inspect near-by grid points
        for iy=-3:3                        % to find this inside the cell
          if (Nx+ix>0)&&(Nx+ix<=Ngx)&&(Ny+iy>0)&&(Ny+iy<=Ngy)  
            ixy=sqrt((xyStrom(ii,1)-(xmin+(Nx+ix)*hg))^2+...
                     (xyStrom(ii,2)-(ymin+(Ny+iy)*hg))^2);
            if (ixy<0.5*diamStroma)     
              oxy(2+Nx+ix,2+Ny+iy)=max(0,oxy(2+Nx+ix,2+Ny+iy)*...
                            (1-dt*oxy_upStrom/(km+oxy(2+Nx+ix,2+Ny+iy))));
            end             % uptake rate M-M method V*o/(km+o)
          end
        end   
      end   
    end
    
     
    
  % draw the results every Ndraw and at the end of simulation 
  if (mod(iter,Ndraw)==0)
      
    % add average oxygen and maxmal force to the history vector  
    Nav=Nav+1;  
    avg_oxy(Nav)=mean(mean(oxy));
    
    
    % this is a color map like in teh experiments
    oxy_map=[ 0,255,255;  0,255,255;  0,215,255;  0,174,255;  0,134,255;...
           0, 94,255;  0, 40,255;  0, 13,255;  5, 11,226; 12, 27,183;...
          47, 29, 93;102,  0,  0;168,  0,  0;233,  0,  0;255, 26,  0;...
         255, 64,  0;255,102,  0;255,153,  0;255,178,  0;255,217,  0;...
         255,255,  0;255,255,225]/255;

     
     
    % clear the figure 
    clf
           
    subplot(1,2,1)  %4 windows for cells and oxygen
      axis([xmin,xmax,ymin,ymax])
      axis equal
      hold on
  
      imagesc(xxg,yyg,oxy(2:Ngx+1,2:Ngy+1)') 
      colorbar
      colormap(oxy_map)
      caxis([0,oxy_max])
      
      viscircles(xyVess,.5*diamVess*ones(Nvess,1),    'color',[0.5,0,0]);
      viscircles(xyStrom,.5*diamStroma*ones(Nstroma,1),'color',[1,1,1]*0.25);
      viscircles(xyTum, .5*diamCell*ones(Ntum,1),     'color', 'k');
   
      axis([xmin,xmax,ymin,ymax])
      axis square
      
      xlabel('[microns]')
      ylabel('[microns]')
      title(['[Nv,Nt,Ns]=[',num2str(Nvess),',',num2str(Ntum),',',...
        num2str(Nstroma),'];  Tum frac=',num2str(Tfrac*100),...
        '%; Str frac=',num2str(Sfrac*100),'%; Vasc frac=',...
        num2str(Vfrac*100),'%'])       
      
            
    subplot(1,2,2)  % 2 windows for average oxygen history
    
       plot(exp_oxy(1,1):Ndraw/Ntime:exp_oxy(1,1)+(Nav-1)*Ndraw/Ntime,...
            avg_oxy(1:Nav), 'bo') 
       hold on
       
       plot(exp_oxy(:,1),exp_oxy(:,2), 'k')
       plot(exp_oxy(:,1),exp_oxy(:,2), 'k*', 'markersize', 12)
      
       title(['Average O2: ',num2str(mean(mean(oxy))),';   iter=',...
         num2str(iter),'; case: ',num2str(case_fluct)]) 
       
       axis([exp_oxy(1,1)-1,exp_oxy(end,1)+1,0,60])
       xticks([4 7 10 13 16 19 22 25 28])
       xlabel('time [minutes]')
       ylabel('oxygen [mmHg]')
       
       grid on 
       axis square  


    pause(0.01) 
    
    end  % end drawing  
    
    
    if (iter==Niter)
      %calculate norm and r2
      norm_res=norm(exp_oxy(:,2)-calc_oxy)/length(calc_oxy);
      r2= rsquare(calc_oxy,exp_oxy(:,2));
        
      title(['Average O2: ',num2str(mean(mean(oxy))),';    R2= ',...
          num2str(r2),';    norm= ',num2str(norm_res),' case: ',...
          num2str(case_fluct)])
    end
  
    if ((toSave==1)||(toSavelast==1))&&(mod(iter,Nsave)==0)
      % save data
      
      print('-djpeg',[pathFigs,'/fig_V',num2str(round(Vfrac*100)),...
         '_T',num2str(round(Tfrac*100)),'_S',num2str(round(Sfrac*100)),...
         '_',num2str(iter)])
      
      save([pathData,'/oxy_',num2str(iter),'.txt'],'oxy','-ascii')
      save([pathData,'/xyVess_',num2str(iter),'.txt'],'xyVess','-ascii')
      save([pathData,'/xyTum_',num2str(iter),'.txt'],'xyTum','-ascii')
      save([pathData,'/xyStrom_',num2str(iter),'.txt'],'xyStrom','-ascii')
      avg_oxy_tosave=avg_oxy';
      save([pathData,'/avg_oxy_',num2str(iter),'.txt'],'avg_oxy_tosave','-ascii')
    end
    
  end  % end of while loop
  
 
toc

if (toSavelast==1)   
  print('-djpeg','-r300',[pathFigs,'/fig_V',num2str(round(Vfrac*100)),'_T',...
        num2str(round(Tfrac*100)),'_S',num2str(round(Sfrac*100))])
   
  save([pathData,'/oxy.txt'],'oxy','-ascii')
  save([pathData,'/xyVess.txt'],'xyVess','-ascii')
  save([pathData,'/xyTum.txt'],'xyTum','-ascii')  
  save([pathData,'/xyStrom.txt'],'xyStrom','-ascii')
  avg_oxy_tosave=avg_oxy';
    
  save([pathData,'/avg_oxy.txt'],'avg_oxy_tosave','-ascii')
    
  parameters=[xmin,xmax,ymin,ymax,eps,hg,dt,Nav,Nvess,Ntum,...
     Nstroma,Ngx,Ngy,stiff,nu,force_thr,error_thr,Vfrac,Tfrac,Sfrac,...
     diamVess,diamCell,diamStroma,Doxy,oxy_max,oxy_upTum,oxy_upStrom,...
     case_fluct,r2,norm_res]';
  save([pathDir,'/parametersFluct',num2str(case_fluct),'.txt'],'parameters','-ascii')
  save([pathDir,'/influx_',num2str(case_fluct),'.txt'],'influx','-ascii')
  save([pathDir,'/uptake_',num2str(case_fluct),'.txt'],'uptake','-ascii')
  
end
  
  
end % if no_go==1

end %main function oxy_tissue_fluc

%----------------------------------------------------------
function r2 = rsquare(x,y)
  
  n=length(y);
  sx=sum(x);  sy=sum(y);  % elements of corr coef calculations
  sxy=sum(x.*y);  sxx=sum(x.*x);  syy=sum(y.*y);

  r=(sxy-sx*sy/n) / sqrt((sxx-(sx*sx)/n)*(syy-(sy*sy)/n) );
  r2=r*r;
end

%----------------------------------------------------------
%----------------------------------------------------------




