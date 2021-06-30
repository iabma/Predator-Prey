% This is a basic predator-prey script that is intended to show
% how to organize a test code.
% The predator and prey strategies are very basic.

function predator_prey

   close all

    Initial_fuel_r = 500000; % Max stored energy for predator
    Initial_fuel_y = 50000;  % Max stored energy for prey   
   
   force_table_prey = rand(51,2)-0.5;
   force_table_predator = rand(51,2)-0.5;

   options = odeset('Events',@event,'RelTol',0.005);
   
   initial_w = [500,0,0,0,0,0,0,0,Initial_fuel_r,Initial_fuel_y]; % Initial position/velocity/energy  
   
   continue_running = true;
   time_vals = [];
   sol_vals = [];
   start_time = 0;
   
   while (continue_running)
       
       % ODE113 will continue running until either 250s is up; a catch occurs; or predator or prey hit the ground
       if (250-start_time>2)
         tspan = start_time:1:250; % This spaces the time intervals by 1s for a smooth animation
       else
         tspan = [start_time,250];
       end
%      ODE45 behaves strangely for this problem so we use the more accurate ode113       
       [time_stage,sol_stage,time_event,sol_event,index_event] = ode113(@(t,w) eom(t,w,force_table_predator,force_table_prey), ...
           tspan,initial_w,options);
       % Store the solution produced by ode113 for plotting
       time_vals = [time_vals;time_stage];
       sol_vals = [sol_vals;sol_stage];
       % Check what happened at the event that terminated ode113 and restart ode113 if predator or prey landed safely
       if (time_stage(end)<250)
          [continue_running,initial_w,start_time] = handle_event(time_event,sol_event,index_event);
       else
           disp('Prey escaped!')
           break;
       end
   end
   
   animate_projectiles(time_vals,sol_vals);
   
   % You might find it helpful to add some code below this line
   % to plot graphs showing what happened during the contest.
   % A few ideas are:
   %(1) Plot the distance between predator & prey as function of time
   %(2) Plot the altitude of predator & prey as function of time
   %(3) Plot speed of predator & prey as functions of time
   %(4) Plot energy of predator & prey as functions of time



end
function dwdt = eom(t,w,force_table_predator,force_table_prey)

    % Extract the position and velocity variables from the vector w
    % Note that this assumes the variables are stored in a particular order in w.
 
    pr = w(1:2); % Predator position (2D column vetor)
    vr = w(5:6); % Predator velocity 
    Er = w(9); % Energy remaining for predator
    py = w(3:4); % Prey position 
    vy = w(7:8); % Prey velocity
    Ey = w(10); % Energy remaining for prey

    %      Constants given in the project description
    g = 9.81;
    mr = 100; % Mass of predator, in kg
    my = 10.; % Mass of prey, in kg
    Frmax = 1.3*mr*g; % Max force on predator, in Newtons
    Fymax = 1.4*my*g; % Max force on prey, in Newtons
    c = 0.2; % Viscous drag coeft, in N s/m
    Eburnrate_r = 0.1;
    Eburnrate_y = 0.2;
    Frrand_magnitude = 0.4*mr*g; % Magnitude of random force on predator
    Fyrand_magnitude = 0.4*my*g; % Magnitude of random force on prey

    % Compute all the forces on the predator
    amiapredator = true;
    Fr = compute_f_mygroupname(t,Frmax,Fymax,amiapredator,pr,vr,Er,py,vy,Ey);
    Frmag = sqrt(dot(Fr,Fr)); % Prevent prey from cheating....
    if (Frmag>Frmax)
        Fr=Fr*Frmax/Frmag;
    end
    if (Er<=0)  % Out of fuel!
        Fr = [0;0];
    end

    Frrand = Frrand_magnitude*compute_random_force(t,force_table_predator); % Random force on predator
    Frvisc = -norm(vr)*vr*c;   % Drag force on predator
    Frgrav = -mr*g*[0;1];      % Gravity force on predator
    Frtotal = Fr+Frrand+Frvisc+Frgrav;  % Total force on predator

    %       If predator is on ground and stationary, and resultant vertical force < 0, set force on predator to zero
    if (pr(2)<=0 && vr(2)<=0 && Frtotal(2)<0)
        Frtotal = [0;0];
    end

    dErdt = -Eburnrate_r*norm(Fr)^(3/2);

    % Write similar code below to call your compute_f_groupname function to
    % compute the force on the prey, determine the random forces on the prey,
    % and determine the viscous forces on the prey
    
    


    amiapredator = false;
    Fy = compute_f_mygroupname(t,Frmax,Fymax,amiapredator,pr,vr,Er,py,vy,Ey);
    Fymag = sqrt(dot(Fy,Fy));
    if (Fymag>Fymax)
        Fy=Fy*Fymax/Fymag;
    end
    if (Ey<=0)  % Out of fuel!
        Fy = [0;0];
    end
    
    Fyrand = Fyrand_magnitude*compute_random_force(t,force_table_prey); % Random force on prey
    Fyvisc = -norm(vy)*vy*c; 
    Fygrav = -my*g*[0;1];      % Gravity force on predator
    Fytotal = Fy+Fyrand+Fyvisc+Fygrav;  % Total force on predator

    %       If predator is on ground and stationary, and resultant vertical force < 0, set force on predator to zero
    if (py(2)<=0 && vy(2)<=0 && Fytotal(2)<0)
        Fytotal = [0;0];
    end

    dEydt = -Eburnrate_y*norm(Fy)^(3/2);
    
 
    


    dwdt = [vr;vy;Frtotal/mr;Fytotal/my;dErdt;dEydt];

    %      This displays a message every time 10% of the computation
    %      is completed, so you can see if your code hangs up.
    %      If the code does hang up, it is because your predator or prey
    %      algorithm is making the applied forces fluctuate too rapidly.
    %      To fix this you will need to modify your strategy.
    %
    persistent percentcompleted % A persistent variable retains its value after the function has completed
    if (isempty(percentcompleted))
        percentcompleted = false(1,11);
    end
    ii = floor(t/25)+1;
    if mod(floor(t/2.5),10)==0
        if (~percentcompleted(ii))
            fprintf('%d %s \n',floor(t/2.5),'% complete')
            percentcompleted(ii) = true;
        end
    end


end
    function [ev,s,dir] = event(t,w)
        pr = w(1:2); py = w(3:4); % Positions of predator and prey
%       ev is a vector - ev(1) = 0 for a catch; ev(2)=0 for predator landing; ev(3) = 0 for prey landing        
        ev = [norm(pr-py)-1.;pr(2);py(2)];
%       All three events will stop the ODE solver
        s = [1;1;1];
%       We want all three to cross zero from above.
        dir = [-1;-1;-1];
    end
    
function [continue_running,initial_w,start_time] = handle_event(event_time,event_sol,event_index)

    
    predator_crash_limit = 15; % Predator max landing speed to survive
    prey_crash_limit = 8; % Prey max landing speed to survive
    Max_fuel_r = 500000; % Max stored energy for predator
    Max_fuel_y = 50000;  % Max stored energy for prey
    
    if (isempty(event_index))
        continue_running = false;
        start_time = event_time;
        initial_w = event_sol;
        return
    end
    continue_running = true;
    initial_w = event_sol;
    start_time = event_time;
    vr = event_sol(5:6); vy = event_sol(7:8);


    if (event_index==1) % Catch 
        continue_running = false;
        disp('Prey was caught!')
    elseif (event_index==2) % Predator landed or crashed
        vrmag = norm(vr);
        if (vrmag>predator_crash_limit ) % Crash
            continue_running = false;
            disp('Predator crashed!')
        else
            initial_w = [initial_w(1),0.0,initial_w(3:4),0.0,0.0,initial_w(7:8),Max_fuel_r,initial_w(10)]; % Set predator velocity to zero; refuel
            disp('Predator landed & refueled!')
        end
    elseif (event_index==3)  % Prey landed or crashed
        vymag = norm(vy);
        if (vymag>prey_crash_limit ) % Crash
            continue_running = false;
            disp('Prey crashed!')
        else % Safe landing
            initial_w = [initial_w(1:2),initial_w(3),0.0,initial_w(5:6),0.0,0.0,initial_w(9),Max_fuel_y]; % Set prey velocity to zero; refuel
            disp('Prey landed & refueled!')
        end
    end
    
    
end
    
%% CHANGE THE NAME OF THE FUNCTION TO A UNIQUE GROUP NAME BEFORE SUBMITTING    
function F = compute_f_mygroupname(t,Frmax,Fymax,amiapredator,pr,vr,Er,py,vy,Ey)


% PLEASE FILL OUT THE INFORMATION BELOW WHEN YOU SUBMIT YOUR CODE
% Test time and place: Enter the time and room for your test here 
% Group members: list the names of your group members here


%   t: Time
%   Frmax: Max force that can act on the predator
%   Fymax: Max force that can act on th eprey
%   amiapredator: Logical variable - if amiapredator is true,
%   the function must compute forces acting on a predator.
%   If false, code must compute forces acting on a prey.
%   pr - 2D column vector with current position of predator eg pr = [x_r;y_r]
%   vr - 2D column vector with current velocity of predator eg vr= [vx_r;vy_r]
%   Er - energy remaining for predator
%   py - 2D column vector with current position of prey py = [x_prey;y_prey]
%   vy - 2D column vector with current velocity of prey py = [vx_prey;vy_prey]
%   Ey - energy remaining for prey
%   F - 2D column vector specifying the force to be applied to the object
%   that you wish to control F = [Fx;Fy]
%   The direction of the force is arbitrary, but if the
%   magnitude you specify exceeds the maximum allowable
%   value its magnitude will be reduced to this value
%   (without changing direction)

    g = 9.81;
    mr = 100; % Mass of predator, in kg
    my = 10.; % Mass of prey, in kg
    predator_crash_limit = 15; % Predator max landing speed to survive
    prey_crash_limit = 8; % Prey max landing speed to survive
    Max_fuel_r = 500000; % Max stored energy for predator
    Max_fuel_y = 50000;  % Max stored energy for prey

  %
  if (amiapredator)
    % Code to compute the force to be applied to the predator
   
     
    dt=8;
    if (norm(py-pr)<15)
        dt=2;
    end
%     if (py(2)==0&&py(1)==0)
%         F=[-1;0];
%     end
    F=py+dt*vy-(pr+dt*vr);%
    F=Frmax*F/norm(F);
   
  
   
        
      
 
  else
%     Code to compute the force to be applied to the prey
    if (norm(py-pr)<50)&&(norm(py-pr)>20)
       
        
        F=Fymax*[-vy(1);-vy(2)];
        
    else
    F=[sin(t);2+cos(t)];
    F=Fymax*(F/norm(F));
    end

  
   end
  
end
%%
function F = compute_random_force(t,force_table)
% Computes value of fluctuating random force at time t, where 0<t<250.
% The variable force_table is a 251x2 matrix of pseudo-random
% numbers between -0.5 and 0.5, computed using
% force_table = rand(251,2)-0.5;
% The force is in Newtons – if you use another system of units you
% must convert.
F = [interp1([0:5:250],force_table(:,1),t);interp1([0:5:250],force_table(:,2),t)];
end

function animate_projectiles(t,sols)

Max_fuel_r = 500000; % Max stored energy for predator
Max_fuel_y = 50000;  % Max stored energy for prey

scrsize = get(0,'ScreenSize');
nn = min(scrsize(3:4));
x0 = scrsize(3)/2-nn/3;
y0 = scrsize(4)/3;
figure1 = figure('Position',[x0,y0,2*nn/3,nn/3]);

xmax = max(max(sols(:,3)),max(sols(:,1)));
xmin = min(min(sols(:,3)),min(sols(:,1)));
ymax = max(max(sols(:,4)),max(sols(:,2)));
ymin = min(min(sols(:,4)),min(sols(:,2)));

dx = 0.1*(xmax-xmin)+0.5;
dy = 0.1*(ymax-ymin)+0.5;

for i = 1:length(t)
    clf
    axes1 = axes('Parent',figure1,'Position',[0.08 0.06914 0.44 0.8]);

    plot(axes1,sols(1:i,3),sols(1:i,4),'LineWidth',2,'LineStyle',...
    ':','Color',[0 0 1]);
    ylim(axes1,[ymin-dy ymax+dy]);
    xlim(axes1,[xmin-dx xmax+dx]);
    hold on
    plot(axes1,sols(1:i,1),sols(1:i,2),'LineWidth',2,'LineStyle',':',...
    'Color',[1 0 0]);
    plot(axes1,sols(i,1),sols(i,2),'ro','MarkerSize',11,'MarkerFaceColor','r');
    plot(axes1,sols(i,3),sols(i,4),'ro','MarkerSize',5,'MarkerFaceColor','g');
    if (ymin-dy<0) 
        fill([xmin-dx;xmax+dx;xmax+dx;xmin-dx;xmin-dx],[0;0;-dy;-dy;0],'g');
    end
    axes2 = axes('Parent',figure1,'Position',[0.62 0.06914 0.3 0.8]);
    bar(axes2,[100*sols(i,9)/Max_fuel_r,0],'FaceColor','r');
    hold on
    bar(axes2,[0,100*sols(i,10)/Max_fuel_y],'FaceColor','b');
    ylim(axes2, [0 120 ]);
    ylabel('%')
    title({'Energy'},'FontSize',14);
    pause(0.1);
end
end



