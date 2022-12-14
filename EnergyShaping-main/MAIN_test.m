% This code is used to develop a energy shaping controller for
% quadruped locomotion

% Reference code used from preprint available at: https://arxiv.org/abs/2012.10002
% Author: Yanran Ding
% Modified by: Sriram SKS Narayanan

% Notes on control
% force and moment are epressed in inertial frame

% Notes on ref traj:
% for walking in a straight line, x(t) is obtained using kinematics
% robot accelerates with constant acc p.acc_d until it reaches the desired
% final velocity p.vel_d
% for trotting gait, u_ref_i is taken as [f_x, f_y, f_z]_i
% f_z = mg/(no of feet in contact), 
% f_x and f_y are set to zero. 
% A contact sequence is hard coded to get desired gait pattern setting u_ref_i of
% swing feet to zero

% Notes on energy shaping:
% optimization is done using gurobi
% objective is norm(u-uref)
% Adding Selection Matrix to Optimization to snap forces of swing feet to
% zero makes the problem infeasible.
% Adding tight bounds on f_z makes the problem infeasible
% the controller overshoots the desired velocity setpoint initllay
% controller is stable over long periods of time

% Updates:
% Working animaiton with Energy Shaping control
% Gains are tuned to the best possible value manually
% 1norm - non zero alphas, 2norm or fro - zero alphas

% To-do:
% fix contact selection problem
% change to circular traj with MPC controller (without ES)

% To check
% add same Q and R weights from MPC

%% initialization
clear all;close all;clc
addpath fcns fcns_MPC

%% --- parameters ---
% ---- gait ----
% 0-trot; 1-bound; 2-pacing 3-gallop; 4-trot run; 5-crawl
gait = 0;
p = get_params(gait);
p.playSpeed = 1;
p.flag_movie = 1;       % 1 - make movie
use_qpSWIFT = 0;        % 0 - quadprog, 1 - qpSWIFT (external)

dt_sim = p.simTimeStep;
SimTimeDuration = 2;  % [sec]
MAX_ITER = floor(SimTimeDuration/p.simTimeStep);

% desired trajectory
p.acc_d = 0.5;
p.vel_d = [0.5;0];
p.yaw_d = 0;

% --- initial condition ---
% Xt = [pc dpc vR wb pf]': [30,1]
if gait == 1
    [p,Xt,Ut] = fcn_bound_ref_traj(p);
else
    [Xt,Ut] = fcn_gen_XdUd(0,[],[1;1;1;1],p);
end

% --- logging ---
tstart = 0;
tend = dt_sim;

[tout,Xout,Uout,Xdout,Udout,Uext,FSMout] = deal([]);

%% --- simulation ----
h_waitbar = waitbar(0,'Calculating...');
tic

r=[];
Coeff=[]; B1=[];
U_ES = []; U_ES_out = [];
Ut_ref = []; Ud_ref=[];
Selection = [];

for ii = 1:MAX_ITER
    % --- time vector ---
    t_ = dt_sim * (ii-1) + p.Tmpc * (0:p.predHorizon-1);
    
    % --- FSM ---
    if gait == 1
        [FSM,Xd,Ud,Xt] = fcn_FSM_bound(t_,Xt,p);
    else
        [FSM,Xd,Ud,Xt] = fcn_FSM(t_,Xt,p);
    end
    %% --- MPC ---- (used only for comparison)
    % form QP
    [H,g,Aineq,bineq,Aeq,beq] = fcn_get_QP_form_eta(Xt,Ut,Xd,Ud,p);

    if ~use_qpSWIFT
        % solve QP using quadprog
        [zval] = quadprog(H,g,Aineq,bineq,Aeq,beq,[],[]);
    else
        % interface with the QP solver qpSWIFT
        [zval,basic_info] = qpSWIFT(sparse(H),g,sparse(Aeq),beq,sparse(Aineq),bineq);
    end

    %get the foot forces value for first time step and repeat
    Ut = Ut + zval(1:12);
    
    %% --- Energy shaping --- (used for simulating dynamics and plot animations)
    % parameterize Ut as r*alpha + beta
    % get desired values
    i_hor=1; %gives desired values for current timestep
    Ut_ref = [Ut_ref, Ut]; % collect Ut_ref from MPC output to plot
    Ud_ref = [Ud_ref, Ud(:,i_hor)]; % collect desried Ud_ref from FSM to plot
    
    xd = Xd(1:3,i_hor);
    vd = Xd(4:6,i_hor);
    Rd = reshape(Xd(7:15,i_hor),[3,3]);
    wd = Xd(16:18,i_hor);
    X_des = [xd;vd;Xd(7:15,i_hor);wd];

    % get current values
    w = Xt(16:18); % get angular velocity in body frame wb
    x = Xt(1:3); % get position x
    v = Xt(4:6); % get velocity v
    R = reshape(Xt(7:15),[3,3]); % Rotation matrix R
    X_cur = [x;v;Xt(7:15);w];

    % foot positions r_pos in R^3x1 and r in R^3x4
    pf34 = reshape(Xt(19:30),[3,4]);
    pc = reshape(Xt(1:3),[3,1]);
    r34 = pf34 - repmat(pc,[1,4]);
    r_pos = reshape(r34,12,1);
    r =[r;r34];
    r1_hat = hatMap(r34(:,1));
    r2_hat = hatMap(r34(:,2));
    r3_hat = hatMap(r34(:,3));
    r4_hat = hatMap(r34(:,4));
    r21 = veeMap(r2_hat-r1_hat); %vector(r2-r1)
    r31 = veeMap(r3_hat-r1_hat); %vector(r3-r1)
    r41 = veeMap(r4_hat-r1_hat); %vector(r4-r1)

    % setup gains
    Kr = [10 0 0;0 10 0;0 0 10]; % check ki + kj > 0
    Kw = [1 0 0;0 1 0;0 0 1]; % positive definite
    Kp = [1e2 0 0;0 1e2 0;0 0 1e2];
    Kd = [1e2 0 0;0 1e2 0;0 0 1e2];

    % define errors
    R_error_right = Rd'*R;
    R_error_left = R*Rd';
    eR_right = 1/2 * veeMap(R_error_right - R_error_right');
    eR_left = 1/2 * veeMap(R_error_left - R_error_left');
    eW_right = w - R_error_right' * wd;
    eW_left = w - wd;
    ex = x - xd;
    ev = v - vd;

    % integral error
    KI_pos = 1; KI_ang = 1;
    p_pos = 1; p_ang = 1;
    eI_pos = ev + p_pos * ex;
    eI_ang = eW_right + p_ang * eR_right;

    % set max values of fi_z (workaround)
    %fi_z_lb = -[40;40;40;40];
    %fi_z_ub = [40;40;40;40];

    % this constraint makes solution infeasible
    % But, this is the right way to add bounds on fz
    Fzd = Ud([3 6 9 12],i_hor);
    fi_z_lb = -10 * max(Fzd);
    fi_z_ub = 10 * max(Fzd);

    % set up selection matrix (makes problem infeasible)
    % S = 1 if feet not in contant
    S = Ud(:,i_hor);
    S(S~=0)=1;  

    % YALMIP for Energy Shpaing optimiation   
    fi = sdpvar(3,4);
    
    % get only the stance feet
    % this leads to infeasible solution
    %idx = any(reshape(S,3,4));
    %stance_feet = ones(3,4).*idx;
    %fi = fi.*stance_feet; %select only the feet in stance 
    fi = reshape(fi,12,1);

    fi_x = fi([1,4,7,10]);
    fi_y = fi([2,5,8,11]);
    fi_z = fi([3,6,9,12]);

    f1 = fi(1:3);
    f2 = fi(4:6);
    f3 = fi(7:9);
    f4 = fi(10:12);
    f_net = f1 + f2 + f3 + f4;

    M_net = R'*(r1_hat*f1 + r2_hat*f2 + r3_hat*f3 + r4_hat*f4);
    
    % warm start
    assign(fi,Ud(:,i_hor));

    Objective = norm(Ud(:,i_hor)-fi, 2);
    
    Constraints = [
        [Rd'*Kr*eR_left + Kw*eW_left + KI_ang*eI_ang + M_net == hatMap(wd)*p.J*w]:'ES constaint 1 left',

%         [Kr*eR_right + Kw*eW_right + M_net ...        
%         == hatMap(w)*p.J*(R_error_right'*wd)]:'ES constaint 1 right',
       
        [f_net + Kp*ex + Kd*ev + KI_pos*eI_pos == p.mass*([p.acc_d;0;0]) + p.mass*[0;0;p.g]]:'ES constraint 2'

        % Constraints to make f_z positive for all 4 legs
        %[fi_z >= 0]:'Enforce positivity of fi_z',

        % Apply selection matrix for f_i (to force gait pattern)
        % leads to infeasible solution
        %[~S([3,6,9,12]).*fi_z == 0],

        % Friction constraints (using linear friction pyramid)
        [-p.mu*fi_z <= fi_x <= p.mu*fi_z]:'Friction Constraint on fi_x 1',
        [-p.mu*fi_z <= -fi_x <= p.mu*fi_z]:'Friction Constraint on fi_x 2',
        [-p.mu*fi_z <= fi_y <= p.mu*fi_z]:'Friction Constraint on fi_y 1',
        [-p.mu*fi_z <= -fi_y <= p.mu*fi_z]:'Friction Constraint on fi_y 2',

        % Lower and upper bounds on f_z
        [fi_z_lb <= fi_z <= fi_z_ub]: 'Bounds on fi_z',

        % bound z position to be > 0 always
        %[x(3) >= 0]
        ];   
        
    opt = sdpsettings('solver','mosek','verbose',2,'cachesolvers',1);
    sol = optimize(Constraints,Objective,opt);

    % collect all the coeff
    u_ES = value(fi);

    % collect all the U_ES
    U_ES = [U_ES,u_ES];

    % collect selection 
    Selection = [Selection, S];

    %% --- external disturbance ---
    [u_ext,p_ext] = fcn_get_disturbance(tstart,p);
    p.p_ext = p_ext;        % position of external force
    u_ext = 0*u_ext;        % set externam force to zero for now
    
    %% --- simulate without any external disturbances ---
    % simulate nonlinear SRB dynamics with ENERGY SHAPING control u_ES
    [t,X] = ode45(@(t,X)dynamics_SRB(t,X,u_ES,Xd,u_ext,p),[tstart,tend],Xt);
    
    
    %% --- update ---
    Xt = X(end,:)';
    tstart = tend;
    tend = tstart + dt_sim;
    
    % --- log ---  
    lent = length(t(2:end));
    tout = [tout;t(2:end)];
    Xout = [Xout;X(2:end,:)];
    Uout = [Uout;repmat(Ut',[lent,1])];
    U_ES_out = [U_ES_out;repmat(u_ES',[lent,1])]; 
    Xdout = [Xdout;repmat(Xd(:,1)',[lent,1])];
    Udout = [Udout;repmat(Ud(:,1)',[lent,1])];
    Uext = [Uext;repmat(u_ext',[lent,1])];
    FSMout = [FSMout;repmat(FSM',[lent,1])];
    
    waitbar(ii/MAX_ITER,h_waitbar,'Calculating...');
end
close(h_waitbar)
fprintf('Calculation Complete!\n')
toc

%% Animation
% generate animation with ENERGY SHAPING Control with u_ES_out
[t,EA,EAd] = fig_animate(tout,Xout,U_ES_out,Xdout,Udout,Uext,p);

%% Compare U_ES and Ut_ref
%% Z direction foot forces
figure()
subplot(2,2,1)
plot(Ud_ref(3,:),'--'); hold on;
plot(U_ES(3,:)); hold on;
plot(Ut_ref(3,:))
title('Z direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','energy shaping','MPC')
hold off

subplot(2,2,2)
plot(Ud_ref(6,:),'--'); hold on;
plot(U_ES(6,:)); hold on;
plot(Ut_ref(6,:))
title('Z direction foot force of front right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','energy shaping','MPC')
hold off

subplot(2,2,3)
plot(Ud_ref(9,:),'--'); hold on;
plot(U_ES(9,:)); hold on;
plot(Ut_ref(9,:))
title('Z direction foot force of rear left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','energy shaping','MPC')
hold off

subplot(2,2,4)
plot(Ud_ref(12,:),'--'); hold on;
plot(U_ES(12,:)); hold on;
plot(Ut_ref(12,:))
title('Z direction foot force of rear right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','energy shaping','MPC')
hold off

%% X direction foot forces
figure()
subplot(2,2,1)
plot(U_ES(1,:)); hold on;
plot(Ut_ref(1,:))
title('X direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,2)
plot(U_ES(4,:)); hold on;
plot(Ut_ref(4,:))
title('X direction foot force of front right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,3)
plot(U_ES(7,:)); hold on;
plot(Ut_ref(7,:))
title('X direction foot force of rear left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,4)
plot(U_ES(10,:)); hold on;
plot(Ut_ref(10,:))
title('X direction foot force of rear right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

%% Y direction foot forces
figure()
subplot(2,2,1)
plot(U_ES(2,:)); hold on;
plot(Ut_ref(2,:))
title('Y direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,2)
plot(U_ES(5,:)); hold on;
plot(Ut_ref(5,:))
title('Y direction foot force of front right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,3)
plot(U_ES(8,:)); hold on;
plot(Ut_ref(8,:))
title('Y direction foot force of rear left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,4)
plot(U_ES(11,:)); hold on;
plot(Ut_ref(11,:))
title('Y direction foot force of rear right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

%% Robot path
figure()
plot3(Xout(:,1), Xout(:,2), Xout(:,3)); grid on
%plot(Xdout(:,1), Xdout(:,2)); grid on
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
title('Robot path in 3D')
