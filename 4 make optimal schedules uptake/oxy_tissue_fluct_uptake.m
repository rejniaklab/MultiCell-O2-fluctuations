function [r2,norm_res]=oxy_tissue_fluct_uptake(Vfrac,Tfrac,Sfrac,...
tiss_pathIn,tiss_pathOut,exp_num)
% Vfrac -- fraction of tissue that is vasculature, i.e. 0.05 means 5%
% Tfrac -- fraction of tissue that is tumor, i.e. 0.2 means 20%
% Sfrac -- fraction of tissue that is stromal cells, i.e. 0.1 means 10%
% tiss_pathIn  -- name of the directory with simualted data
% tiss_pathOut -- name of the directory to save data
% exp_num -- determines expeirmental data to read



  % default data
  r2=0; norm_res=100;

  % read experimental data
  if     exp_num==1 exp_oxy=load('Exp_data/black_exp.txt');
  elseif exp_num==2 exp_oxy=load('Exp_data/red_exp.txt');  
  elseif exp_num==3 exp_oxy=load('Exp_data/blue_exp.txt');
  elseif exp_num==4 exp_oxy=load('Exp_data/magenta_exp.txt');
  end

  % read optimal influx and uptake rates 
  influx=load([tiss_pathOut,'/rateInflux_opt_uptake.txt']);
  uptake=load([tiss_pathOut,'/rateUptake_opt_uptake.txt']);


  % parameters
  dt=0.05;        % time step [sec]
  Niter=28*60/dt; % 30 minutes
  Ntime=60/dt;    %  1 minute interval 


  toSave=1;         % switch =1 will save data ==0 no data save
  toSavelast=1;
  Nsave=floor(60/dt);  % frequency of saving the data
  Ndraw=floor(60/dt);  % frequency of drawing the data


  pathData=[tiss_pathOut,'/DataFinalUptake']; % name of the directiory for saving data, 2 for influx, 3 for uptake
  pathFigs=[tiss_pathOut,'/FigsFinalUptake']; % name of the directiory for saving figs    


  % make output directory
  if (toSave==1)||(toSavelast==1)   
    if ~exist(pathData,'dir') mkdir(pathData); end
    if ~exist(pathFigs,'dir') mkdir(pathFigs); end
  end


  % parameters
  stiff=50;        % how hard the cells are pushing away from each other
  nu=  250/10;     % viscosity of the medium

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


  % to remember the average oxygen and maximal force
  avg_oxy=zeros(1,1);
  Nav=0;   % index for the above vectors


  
  % domain 
  [xx,yy]=meshgrid(xmin:hg:xmax,ymin:hg:ymax); % grid for oxygen 
  xxg=xx(1,:);  yyg=yy(:,1);  Ngx=length(xxg);  Ngy=length(yyg);

  % create the figure
  figure('position',[750,500,1500,500])
  pause(0.1)

  
  calc_oxy=zeros(9,1);

  % read intital data from outside files
  tiss_pathIn=[tiss_pathIn,'/data_V',num2str(Vfrac),'_T',num2str(Tfrac),...
                 '_S',num2str(Sfrac),'/Data/'];
  xyVess =load([tiss_pathIn,'xyVess.txt']);
  xyTum  =load([tiss_pathIn,'xyTum.txt']);
  xyStrom=load([tiss_pathIn,'xyStrom.txt']);
  oxy    =load([tiss_pathOut,'/Data_uptake/oxy_1.txt']); %%KAR 11/25
  Nvess=size(xyVess,1);
  Ntum=size(xyTum,1);
  Nstroma=size(xyStrom,1);
  vess_oxy=oxy_max*ones(Nvess,1);


  
% -------------- main code ---------------- %
  
  iter=-1;                            % initial iteration  
   
  
% ----- main loop ----- %  
 tic
 
 iter=4*Ntime-1;
 while (iter<Niter)
    iter=iter+1;   % increase iteration
       
    if (iter<=4*Ntime) %300 incerements because the paper recorded fluc every 3 minutes
      if iter==4*Ntime
        calc_oxy(1,1)=mean(mean(oxy));
      end
      vess_oxy =influx(1)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(1)*oxy_upTum0;        
    elseif (iter<=7*Ntime)
      if iter==7*Ntime
        calc_oxy(2,1)=mean(mean(oxy));
      end
      vess_oxy =influx(2)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(2)*oxy_upTum0;        
    elseif (iter<=10*Ntime)
      if iter==10*Ntime
        calc_oxy(3,1)=mean(mean(oxy));
      end
      vess_oxy =influx(3)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(3)*oxy_upTum0;        
    elseif (iter<=13*Ntime)
      if iter==13*Ntime
        calc_oxy(4,1)=mean(mean(oxy));
      end
      vess_oxy =influx(4)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(4)*oxy_upTum0;       
    elseif (iter<=16*Ntime)
      if iter==16*Ntime
        calc_oxy(5,1)=mean(mean(oxy));
      end
      vess_oxy =influx(5)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(5)*oxy_upTum0;       
    elseif (iter<=19*Ntime)
      if iter==19*Ntime
        calc_oxy(6,1)=mean(mean(oxy));
      end
      vess_oxy =influx(6)*oxy_max*ones(Nvess,1);
      oxy_upTum=uptake(6)*oxy_upTum0;        
    elseif (iter<=22*Ntime)
      if iter==22*Ntime
        calc_oxy(7,1)=mean(mean(oxy));
      end
      vess_oxy =influx(7)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(7)*oxy_upTum0;        
    elseif (iter<=25*Ntime)
      if iter==25*Ntime
        calc_oxy(8,1)=mean(mean(oxy));
      end
      vess_oxy =influx(8)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(8)*oxy_upTum0;       
    elseif (iter<=28*Ntime)
      if iter==28*Ntime
        calc_oxy(9,1)=mean(mean(oxy));
      end
      vess_oxy =influx(9)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(9)*oxy_upTum0;        
    else
      vess_oxy =influx(9)*oxy_max*ones(Nvess,1); 
      oxy_upTum=uptake(9)*oxy_upTum0;       
    end
    
           
      
    % define influx of oxy from vessels
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
        num2str(Nstroma),'];  Tum frac=',num2str(Tfrac),'%; Str frac=',...
        num2str(Sfrac),'%; Vasc frac=',num2str(Vfrac),'%'])       
      
            
    subplot(1,2,2)  % 2 windows for average oxygen history
    
       plot(exp_oxy(1,1):Ndraw/Ntime:exp_oxy(1,1)+(Nav-1)*Ndraw/Ntime,...
            avg_oxy(1:Nav), 'bo') 
       hold on
       
       plot(exp_oxy(:,1),exp_oxy(:,2), 'k')
       plot(exp_oxy(:,1),exp_oxy(:,2), 'k*', 'markersize', 12)
      
       title(['Average O2: ',num2str(mean(mean(oxy))),';   iter=',...
         num2str(iter),'; case: uptake']) 
       
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
          num2str(r2),';    norm= ',num2str(norm_res),' case: uptake'])
    end
  
    if ((toSave==1)||(toSavelast==1))&&(mod(iter,Nsave)==0)
      % save data
      
      print('-djpeg',[pathFigs,'/fig_V',num2str(Vfrac),'_T',...
        num2str(Tfrac),'_S',num2str(Sfrac),'_',num2str(iter),'.jpg'])
      
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
  print('-djpeg','-r300',[pathFigs,'/fig_V',num2str(Vfrac),'_T',...
        num2str(Tfrac),'_S',num2str(Sfrac),'.jpg'])
   
  save([pathData,'/oxy.txt'],'oxy','-ascii')
  save([pathData,'/xyVess.txt'],'xyVess','-ascii')
  save([pathData,'/xyTum.txt'],'xyTum','-ascii')  
  save([pathData,'/xyStrom.txt'],'xyStrom','-ascii')
  avg_oxy_tosave=avg_oxy';
    
  save([pathData,'/avg_oxy.txt'],'avg_oxy_tosave','-ascii')
    
  parameters=[xmin,xmax,ymin,ymax,eps,hg,dt,Nav,Nvess,Ntum,Nstroma,...
     Ngx,Ngy,stiff,nu,Vfrac,Tfrac,Sfrac,diamVess,diamCell,diamStroma,...
     Doxy,oxy_max,oxy_upTum,oxy_upStrom,r2,norm_res]';
  save([tiss_pathOut,'/parametersFluct_uptake.txt'],'parameters','-ascii')
  
end
  
  

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




