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
% optimization is done using mosek
% objective is norm(u-uref)
% Adding tight bounds on f_z makes the problem infeasible
% the controller overshoots the desired velocity setpoint initllay
% controller is stable over long periods of time

% Updates:
% Working animaiton with Energy Shaping control
% Gains are are added to the optimization variable

% To-do:
% change to circular traj with MPC controller (without ES)

% To check
% add same Q and R weights from RF-MPC

%% initialization
clear all;close all;clc
addpath fcns fcns_MPC

%% --- parameters ---
% ---- gait ----
% 0-trot; 1-bound; 2-pacing 3-gallop; 4-trot run; 5-crawl
gait = 0;
p = get_params(gait);
p.playSpeed = 3;
p.flag_movie = 1;       % 1 - make movie
use_qpSWIFT = 0;        % 0 - quadprog, 1 - qpSWIFT (external)

dt_sim = p.simTimeStep;
SimTimeDuration = 1;  % [sec]
MAX_ITER = floor(SimTimeDuration/p.simTimeStep);

% desired trajectory
p.acc_d = 0.5;
p.vel_d = [0.25; 0];
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
KP = []; KD = []; 
KW = []; KR = [];
Slack = []; STORAGE = 10;
ER = []; EW = []; Change = [];
V = [0;0;0];
eI_pos = 0; eI_ang = 0;

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
    V = [V,v]; 
    a = (V(:,end) - V(:,end-1)) / dt_sim;
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
    p_pos = 1; p_ang = 1;
    eI_pos = eI_pos + ev + p_pos * ex;
    eI_ang = eI_ang + eW_right + p_ang * eR_right;

    % set max values of fi_z
    Fzd = Ud([3 6 9 12],i_hor);
    fi_z_lb = -2 * Fzd;
    fi_z_ub = 2 * Fzd;

    % set up selection matrix (makes problem infeasible)
    % S = 0 if feet not in contant
    S = Ud(:,i_hor);
    S(S~=0)=1;  

    %% YALMIP for Energy Shpaing optimiation   
    fi = sdpvar(3,4);

    %Kp_nom = [1e2 0 0;0 1e2 0;0 0 1e2];
    %Kr_nom = [10 10 10]'; % check ki + kj > 0
    Kp = sdpvar(3,3,'full'); % fix 
    Kd = sdpvar(3,3,'full');
    Kr = sdpvar(3,1); % fix
    Kw = sdpvar(3,3,'full');
    Ki_pos = 1e1*[1;1;1];
    Ki_ang = 1e1*[1;1;1];
    slack = sdpvar(7,1);
   
    KP = [KP, Kp];
    KR = [KR, Kr];

    % change Kp Kr only when contact changes
    %change = 0;
    %if(ii>1)
    %    change = 1e6;
    %    if(any(Ud_ref(:,end-1)~=Ud_ref(:,end)))
    %        change = 0; % if ref does not change, dont penalize
    %    end
    %end
    %Change = [Change, change];

    % get only the stance feet
    idx = any(reshape(S,3,4));
    stance_feet = ones(3,4).*idx;
    fi = fi.*stance_feet; %select only the feet in stance 
    fi = reshape(fi,12,1);

    % get net force in each direction
    fi_x = fi([1,4,7,10]);
    fi_y = fi([2,5,8,11]);
    fi_z = fi([3,6,9,12]);
    
    % get force of each foot
    f1 = fi(1:3);
    f2 = fi(4:6);
    f3 = fi(7:9);
    f4 = fi(10:12);
    
    % calculate overall net force and moment based on geometric control
    f_net = f1 + f2 + f3 + f4;
    %f_net = dot(f_net,R*[0;0;1]); %body frame
    %f_net = R'*f_net;

    M_net_right = -diag(Kr)*eR_right - Kw*eW_right - diag(Ki_ang)*eI_ang + hatMap(w)*p.J*(R_error_right'*wd);
    %M_net_left = -Rd'*diag(Kr)*eR_left - Kw*eW_left + hatMap(wd)*p.J*w;

    % M from geometric control of quadrotor
%     M_net_right = - diag(Kr)*eR_right ...
%                   - Kw*eW_right ...
%                   - diag(Ki_ang)*eI_ang ...
%                   + hatMap(R_error_right'*wd)*p.J*(R_error_right'*wd);
   
    % calculate storage function as a function of Kp and Kr
    %S1 = 0.5*norm(eW_right,2)^2;
    %S2 = 0.5*trace(diag(Kr)*(eye(3)-eR_right)); % function of Kr
    %S3 = 0.5*p.mass*norm(v-vd)^2;
    %S4 = 0.5*x'*Kp*x;% function of Kp
    %Storage = S1 + S2 + S3 + S4;
    %STORAGE = [STORAGE, Storage];
    
    Qx = 1e6*eye(3);
    Qv = 1e6*eye(3);
    Qr = 1e6*eye(9);
    Qw = 1e6*eye(3);
    Q = blkdiag(Qx, Qv, Qr, Qw);
    Ru = diag([1e-6;1e-6;1e-2]); % Ru for one f
    Ru = blkdiag(Ru, Ru, Ru, Ru);

    % warm start
    %assign(fi,Ud(:,i_hor));
    %assign(Kp,Kp_nom);
    %assign(Kr,Kr_nom);   

    Objective = (X_cur-X_des)'*Q*(X_cur-X_des) ...
        + (Ud(:,i_hor)-fi)'*(Ru)*(Ud(:,i_hor)-fi) ...
        + 1*(norm(Kp,1)+norm(Kd,1))...
        + 1*(norm(Kr,1)+norm(Kw,1))...
        + 1e3*norm(slack,1); ...
        %+ change*(norm(Kp(end)-Kp(end-1),1)) + change*(norm(Kr(end)-Kr(end-1),1));

    Constraints = [
        % ES force constraint (in inertial frame)
        f_net + Kp*ex + Kd*ev + diag(Ki_pos)*eI_pos ...
        == -p.mass*([p.acc_d;0;0]) + p.mass*[0;0;p.g],

        % ES moment constraint (in body frame) 
        R'*(r1_hat*f1 + r2_hat*f2 + r3_hat*f3 + r4_hat*f4) == M_net_right,

        % Storage function constraint
        % STORAGE(end) <= value(STORAGE(end-1)), 
  
        % Friction constraints (using linear friction pyramid)
        -p.mu*fi_z <= fi_x <= p.mu*fi_z,
        -p.mu*fi_z <= -fi_x <= p.mu*fi_z,
        -p.mu*fi_z <= fi_y <= p.mu*fi_z,
        -p.mu*fi_z <= -fi_y <= p.mu*fi_z,

        % Lower and upper bounds on f_z
        fi_z_lb <= fi_z <= fi_z_ub,

        % Positive Defininte constraint on gains 
        % to do - fix strict inequality using slack variables
        Kp >= slack(1),
        Kd >= slack(2),
        Kr >= slack(3),
        Kr(1) + Kr(2) >= slack(4),
        Kr(1) + Kr(3) >= slack(5),
        Kr(2) + Kr(3) >= slack(6),
        Kw >= slack(7),
        slack >= 1e-6,
        ];
        
    opt = sdpsettings('solver','mosek','verbose',2,'cachesolvers',1);
    sol = optimize(Constraints,Objective,opt);
    
    u_ES = value(fi);
    
    % collect all the U_ES
    U_ES = [U_ES,u_ES];
    
    % collect all gains 
    KD = [KD, value(Kd)];
    KW = [KW, value(Kw)];

    % collect errors
    ER = [ER, eR_left];
    EW = [EW, eW_left];

    % collect all slack
    Slack = [Slack, value(slack)];

    % collect selection matrix
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
    Xdout = [Xdout;repmat(Xd(:,i_hor)',[lent,1])];
    Udout = [Udout;repmat(Ud(:,i_hor)',[lent,1])];
    Uext = [Uext;repmat(u_ext',[lent,1])];
    FSMout = [FSMout;repmat(FSM',[lent,1])];
    
    waitbar(ii/MAX_ITER,h_waitbar,'Calculating...');
end
close(h_waitbar)
fprintf('Calculation Complete!\n')
toc

%% offline control
% U_ES_offline = U_ES(:,250:end);
% for i = 1:length(U_ES_offline)
%     [t,X] = ode45(@(t,X)dynamics_SRB(t,X,U_ES_offline(:,i),Xd,u_ext,p),[tstart,tend],Xt);
%  %% --- update ---
%     Xt = X(end,:)';
%     tstart = tend;
%     tend = tstart + dt_sim;
%     
%     % --- log ---  
%     lent = length(t(2:end));
%     tout = [tout;t(2:end)];
%     Xout = [Xout;X(2:end,:)];
%     Uout = [Uout;repmat(Ut',[lent,1])];
%     U_ES_out = [U_ES_out;repmat(U_ES_offline(:,i)',[lent,1])]; 
%     Xdout = [Xdout;repmat(Xd(:,1)',[lent,1])];
%     Udout = [Udout;repmat(Ud(:,1)',[lent,1])];
%     Uext = [Uext;repmat(u_ext',[lent,1])];
%     FSMout = [FSMout;repmat(FSM',[lent,1])];
% end

% past constraint which are not useful
%         % ES plane constraint 1 (in body frame) 
%         r21' * (M_net_right - R'*r1_hat*f_net) == ...
%         r21' * R'*((r2_hat-r1_hat)*f2 + (r3_hat-r1_hat)*f3 + (r4_hat-r1_hat)*f4),
%         
%         % ES plane constraint 2 (in body frame) 
%         r31' * (M_net_right - R'*r1_hat*f_net) == ...
%         r31' * R'*((r2_hat-r1_hat)*f2 + (r3_hat-r1_hat)*f3 + (r4_hat-r1_hat)*f4),
%         
%         % ES plane constraint 3 (in body frame) 
%         r41' * (M_net_right - R'*r1_hat*f_net) == ...
%         r41' * R'*((r2_hat-r1_hat)*f2 + (r3_hat-r1_hat)*f3 + (r4_hat-r1_hat)*f4),
%% Animation
% generate animation with ENERGY SHAPING Control with u_ES_out
[t,EA,EAd] = fig_animate(tout,Xout,U_ES_out,Xdout,Udout,Uext,p);

%% Compare U_ES and Ut_ref
%% Z direction foot forces
figure()
subplot(2,2,1)
%plot(Ud_ref(3,:),'--'); hold on;
plot(U_ES(3,:)); hold on;
%plot(Ut_ref(3,:))
title('Z direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','energy shaping')
hold off

subplot(2,2,2)
%plot(Ud_ref(6,:),'--'); hold on;
plot(U_ES(6,:)); hold on;
%plot(Ut_ref(6,:))
title('Z direction foot force of front right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','energy shaping')
hold off

subplot(2,2,3)
%plot(Ud_ref(9,:),'--'); hold on;
plot(U_ES(9,:)); hold on;
%plot(Ut_ref(9,:))
title('Z direction foot force of rear left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','energy shaping')
hold off

subplot(2,2,4)
%plot(Ud_ref(12,:),'--'); hold on;
plot(U_ES(12,:)); hold on;
%plot(Ut_ref(12,:))
title('Z direction foot force of rear right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','energy shaping')
hold off

%% X direction foot forces
figure()
subplot(2,2,1)
plot(U_ES(1,:)); hold on;
%plot(Ut_ref(1,:))
title('X direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,2)
plot(U_ES(4,:)); hold on;
%plot(Ut_ref(4,:))
title('X direction foot force of front right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,3)
plot(U_ES(7,:)); hold on;
%plot(Ut_ref(7,:))
title('X direction foot force of rear left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,4)
plot(U_ES(10,:)); hold on;
%plot(Ut_ref(10,:))
title('X direction foot force of rear right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

%% Y direction foot forces
figure()
subplot(2,2,1)
plot(U_ES(2,:)); hold on;
%plot(Ut_ref(2,:))
title('Y direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,2)
plot(U_ES(5,:)); hold on;
%plot(Ut_ref(5,:))
title('Y direction foot force of front right foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,3)
plot(U_ES(8,:)); hold on;
%plot(Ut_ref(8,:))
title('Y direction foot force of rear left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('energy shaping','MPC')
hold off

subplot(2,2,4)
plot(U_ES(11,:)); hold on;
%plot(Ut_ref(11,:))
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
