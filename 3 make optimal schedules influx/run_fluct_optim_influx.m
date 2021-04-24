function oxy_out=run_fluct_optim_influx(rate,time,Vfrac,Tfrac,Sfrac,exp_oxy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the mathematical oxygenation model for determining the average   %
% oxygen level 'oxy_out' for the time segment [(tiem-1)*3,time*3] with % 
% rate 'rate' for tissue characterized by [Vfrac,Tfrac,Sfra], and for  %
% exp_oxy corresponding to ROIs:                                       %
% exp_num=  #1 black; #2 red; #3 blue and #4 magenta ROIs              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  dt=0.05;                    % time step [sec]
  Niter=(1+time*3)*60/dt;     % 3 minutes final time
  Ntime=60/dt;                % 1 minute interval 
  iter=(1+(time-1)*3)*60/dt;
 
  Ndraw=floor(60/dt);  % frequency of drawing the data

  tiss_pathIn=['tissue/data_V',num2str(Vfrac),'_T',...
               num2str(Tfrac),'_S',num2str(Sfrac)];           
  tiss_pathOut=['tissue_fluctuations/data_V',num2str(Vfrac),...
       '_T',num2str(Tfrac),'_S',num2str(Sfrac),'/Data_influx'];  
  
  
  
% parameters
  stiff=50;        % how hard the cells are pushing away from each other
  nu=  250/10;     % viscosity of the medium

  diamVess   =40;        % vessel diameter [microns]
  diamCell   =15;        % tumor cell diameter [microns]
  diamStroma = 7.5;      % stromal cell diameter

% domain 
  xmin=-500; xmax=-xmin; ymin=xmin; ymax=-xmin;  % 1000x1000 micron^2=1mm^2
  hg=5;           % grid width [microns]

% oxygen 
  Doxy   =100;    % diffusion coefficient for oxygen [micron^2/s]
  oxy_max= 60;    % level of oxygen in the vessel [mmHg=sigma-gram/micron^3]
  oxy_init= .25; % initial oxygen level

  oxy_upTum0 =9.54/(hg*hg);   % max rate of oxy uptake by tumor cells V
  oxy_upStrom=9.54/(hg*hg);   % max rate of oxy uptake by stromal cells V
  km         =1.134;          % Michealis-Menten coefficient, V*O/(km+O) 

  cfl_cond=Doxy*dt/(hg*hg);   % stability condition
  if (cfl_cond>0.25)
    disp('stability condition is not satisfied'); pause;
  end


% to remember the average oxygen and maximal force
  avg_oxy=zeros(1,1);
  Nav=0;   % index for the above vectors



% domain 
  [xx,yy]=meshgrid(xmin:hg:xmax,ymin:hg:ymax); % grid for oxygen 
  xxg=xx(1,:);  yyg=yy(:,1);  Ngx=length(xxg);  Ngy=length(yyg);

  

%read intital data from outside files
  xyVess =load([tiss_pathIn,'/Data/xyVess.txt']);
  xyTum  =load([tiss_pathIn,'/Data/xyTum.txt']);
  xyStrom=load([tiss_pathIn,'/Data/xyStrom.txt']);
  
  if (time==1)
    oxy=load([tiss_pathIn,'/Data/oxy.txt']);
    save([tiss_pathOut,'/oxy_0.txt'],'oxy','-ascii')
  end
  oxy=load([tiss_pathOut,'/oxy_',num2str(time-1),'.txt']);

  Nvess=size(xyVess,1); Ntum=size(xyTum,1); Nstroma=size(xyStrom,1);
  vess_oxy=oxy_max*ones(Nvess,1);

  oxy_out=mean(mean(oxy));


  % create the figure
  fg=figure('position',[750,500,1000,400]);

  
  
  
% -------------- main code ---------------- % 
% ----- main loop ----- %  
  


  while (iter<Niter)
    iter=iter+1;   % increase iteration
    
    % return final oxygenation level
    if iter==Niter
      oxy_out=mean(mean(oxy));
    end
    % define the rate of oxygen influx 'rate' and oxygen uptake '1'  
    vess_oxy =rate*oxy_max*ones(Nvess,1); 
    oxy_upTum =1*oxy_upTum0;        
    
    
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
    clf(fg);
           
    subplot(1,2,1)  %4 windows for cells and oxygen
      axis([xmin,xmax,ymin,ymax])
      axis equal
      hold on
  
      imagesc(xxg,yyg,oxy(2:Ngx+1,2:Ngy+1)') 
      colorbar
      colormap(oxy_map)
      caxis([0,oxy_max])
      
      axis([xmin,xmax,ymin,ymax])
      axis square
      
      xlabel('[microns]')
      ylabel('[microns]')
      title(['[Nv,Nt,Ns]=[',num2str(Nvess),',',num2str(Ntum),',',...
         num2str(Nstroma),'];  Tum frac=',num2str(Tfrac),'%; Str frac=',...
         num2str(Sfrac),'%; Vasc frac=',num2str(Vfrac),'%'])       
      
            
    subplot(1,2,2)  % 2 windows for average oxygen history
    
       t0=(1+(time-1)*3);
       tn=(1+time*3);  
       plot(t0+Ndraw/Ntime:Ndraw/Ntime:t0+(Ndraw/Ntime)*Nav,avg_oxy(1:Nav),'bo')      
       hold on
       
       plot(exp_oxy(:,1),exp_oxy(:,2), 'k')
       plot(exp_oxy(:,1),exp_oxy(:,2), 'k*', 'markersize', 12)
      
       title(['Average O2: ',num2str(mean(mean(oxy))),';   iter=',num2str(iter)]) 
       
       axis([exp_oxy(1,1)-4,exp_oxy(end,1)+1,0,60])
       xticks([4 7 10 13 16 19 22 25 28])
       xlabel('time [minutes]')
       ylabel('oxygen [mmHg]')
       
       grid on 
       axis square  
       pause(0.1) 
    
    end  % end drawing  
         
  end  % end of while loop
  
  
save([tiss_pathOut,'/oxy_',num2str(time),'.txt'],'oxy','-ascii')
oxy_out=mean(mean(oxy));
close(fg);
end %main function oxy_tissue_fluc

%----------------------------------------------------------
%----------------------------------------------------------




