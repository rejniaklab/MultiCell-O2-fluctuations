function Main_fluct_optim_uptake(Vfrac,Tfrac,Sfrac,textExp,pathOutTiss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimizing uptake schedule to match experimental oxygen fluctuation %
% data for textExp:                                                   %
%   exp_num=  #1 black; #2 red; #3 blue and #4 magenta ROIs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% global variables  
global exp_oxy time iter_opt res_seq
global pathOutTiss pathOutSch Vfrac Tfrac Sfrac textExp 


  % read experimental data of oxygen level
  exp_oxy=load(['Exp_data/',textExp,'_exp.txt']); 

  
  % create a directory for saving the data 
  pathOutSch=[pathOutTiss,'/Data_uptake']  
  if ~exist(pathOutSch, 'dir') mkdir(pathOutSch); end  

  
  % initial guess and lower & upper limits
  x0 =[1];    % initial guess for uptake
  lb =[0];    % lower bound for uptake
  ub =[50];   % upper bound for uptake
  
  
  %  objective function handle (defined below)
  obj=@objective_function;


  % options for stopping optimization process
  opts=optimoptions('patternsearch','MeshTolerance',0.001,...
       'StepTolerance',0.001,'FunctionTolerance',0.001);

   
  % to remember optimal rates for all 9 fluctuation segments  
  rate_influx_tosave=ones(9,1);
  rate_uptake_tosave=zeros(9,1);
  oxy_level_tosave=zeros(9,1);
  
  
  % optimize 9 rates for all fluctuation segments  
  for time=1:9
    iter_opt=0;
    res_seq=[];
    
    tic
    % for objective function 'obj', initial guess 'x0', and lower and 
    % upper bounds 'lb', 'ub', and xsearch options 'opts', this routine 
    % will minimize for optimal difference between simulated and 
    % expeirmental values of oxygen 'oxy_diff', and return the optimal 
    % influx rate 'rate_opt'
    [rate_opt,oxy_diff] = patternsearch(obj,x0,[],[],[],[],lb,ub,opts);
    toc
    
    
    
    % dave the optimal data for this segment as a vector
    rate_uptake_tosave(time,1)=rate_opt; 
      oxy_level_tosave(time,1)=exp_oxy(time,2)-oxy_diff; 
   
    % save the whole data vector for the current results 
    save([pathOutTiss,'/rateInflux_opt_uptake.txt'],'rate_influx_tosave','-ascii')
    save([pathOutTiss,'/rateUptake_opt_uptake.txt'],'rate_uptake_tosave','-ascii')
    save([pathOutTiss,'/oxyLevel_opt_uptake.txt'],'oxy_level_tosave','-ascii')
   
    close all
  end

  
  % save the whole data vector for the final results
  save([pathOutTiss,'/rateInflux_opt_uptake.txt'],'rate_influx_tosave','-ascii')
  save([pathOutTiss,'/rateUptake_opt_uptake.txt'],'rate_uptake_tosave','-ascii')
  save([pathOutTiss,'/oxyLevel_opt_uptake.txt'],'oxy_level_tosave','-ascii')
  
  
end   % Main_fluct_optim_uptake

% ------------------------------------------------------------
% ------------------------------------------------------------

function oxy_diff = objective_function(rate_in)

% global variables
global iter_opt time exp_oxy Vfrac Tfrac Sfrac res_seq
global pathOutSch 


  % count iterations for drawing the current result of optimization
  iter_opt=iter_opt+1;  
  disp(['time: ',num2str(time),' iteration: ',num2str(iter_opt),...
          ' uptake: ',num2str(rate_in)])
  
  % run our oxygenation model to return average level of oxygen 'oxy_out'
  oxy_out = run_fluct_optim_uptake(rate_in,time,Vfrac,Tfrac,Sfrac,exp_oxy);
  oxy_diff = abs((exp_oxy(time,2)-oxy_out));
  % for optimization purpose, we need to return the difference between 
  % simulated and experimental oxygen level -- this value will be mninimized 
  
  
  res_seq=[res_seq;rate_in,oxy_diff];
  save([pathOutSch,'/res_seq_',num2str(time),'.txt'],'res_seq','-ascii')
  
  
  % drawing and saving results
    figure(10)    
    subplot(2,1,1)
      stem(iter_opt,rate_in,'filled','r')
      axis([0,iter_opt+1,0,50])
      hold on  
      grid on
      xlabel('iteration')
      ylabel('uptake rate')
     
      title(['Uptake rate=',num2str(rate_in),';  optim=',num2str(oxy_out),...
           '  exp=',num2str(exp_oxy(time,2)),';  time segment: ',...
           num2str(time)],'fontsize',12)
       

    subplot(2,1,2)
      stem(iter_opt,oxy_out-exp_oxy(time,2),'filled','b')
      axis([0,iter_opt+1,-5,5])
      hold on  
      grid on
      xlabel('iteration')
      ylabel('fitting data')
    
      
   figure(10)
   print('-djpeg',[pathOutSch,'/fig_uptake_',num2str(time),'.jpg'])
   
end  % objective_function

% ------------------------------------------------------------
% ------------------------------------------------------------
