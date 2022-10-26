% This code is used to develop a energy shaping controller for
% quadruped locomotion

% Reference code used from preprint available at: https://arxiv.org/abs/2012.10002
% video available at: https://www.youtube.com/watch?v=iMacEwQisoQ&t=101s

% Author: Yanran Ding
% Modified by: Sriram SKS Narayanan

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
SimTimeDuration = 5;  % [sec]
MAX_ITER = floor(SimTimeDuration/p.simTimeStep);

% desired trajectory
p.acc_d = 1;
p.vel_d = [0.5;0];
p.yaw_d = 0;

%% Model Predictive Control
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

% --- simulation ----
h_waitbar = waitbar(0,'Calculating...');
tic
Ut_ref = []; Ud_ref = [];

for ii = 1:MAX_ITER
    % --- time vector ---
    t_ = dt_sim * (ii-1) + p.Tmpc * (0:p.predHorizon-1);
    
    % --- FSM ---
    if gait == 1
        [FSM,Xd,Ud,Xt] = fcn_FSM_bound(t_,Xt,p);
    else
        [FSM,Xd,Ud,Xt] = fcn_FSM(t_,Xt,p);
    end
    %% --- MPC ----
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
    i_hor = 1;
    Ut_ref = [Ut_ref, Ut];
    Ud_ref = [Ud_ref, Ud(:,i_hor)]; 
    %% --- external disturbance ---
    [u_ext,p_ext] = fcn_get_disturbance(tstart,p);
    p.p_ext = p_ext;        % position of external force
    u_ext = u_ext;
    
    %% --- simulate without any external disturbances ---
    [t,X] = ode45(@(t,X)dynamics_SRB(t,X,Ut,Xd,u_ext,p),[tstart,tend],Xt);
    
    
    %% --- update ---
    Xt = X(end,:)';
    tstart = tend;
    tend = tstart + dt_sim;
    
    %% --- log ---  
    lent = length(t(2:end));
    tout = [tout;t(2:end)];
    Xout = [Xout;X(2:end,:)];
    Uout = [Uout;repmat(Ut',[lent,1])];
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
[t,EA,EAd] = fig_animate(tout,Xout,Uout,Xdout,Udout,Uext,p);

%% Robot path
figure()
plot3(Xout(:,1), Xout(:,2), Xout(:,3)); grid on
%plot(Xdout(:,1), Xdout(:,2)); grid on
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')

%% Z direction foot forces
figure()
subplot(2,2,1)
plot(Ud_ref(3,:),'--'); hold on;
plot(Ut_ref(3,:))
title('Z direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,2)
plot(Ud_ref(6,:),'--'); hold on;
plot(Ut_ref(6,:))
title('Z direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,3)
plot(Ud_ref(9,:),'--'); hold on;
plot(Ut_ref(9,:))
title('Z direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,4)
plot(Ud_ref(12,:),'--'); hold on;
plot(Ut_ref(12,:))
title('Z direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off
%% Y direction foot forces
figure()
subplot(2,2,1)
plot(Ud_ref(2,:),'--'); hold on;
plot(Ut_ref(2,:))
title('Y direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,2)
plot(Ud_ref(5,:),'--'); hold on;
plot(Ut_ref(5,:))
title('Y direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,3)
plot(Ud_ref(5,:),'--'); hold on;
plot(Ut_ref(5,:))
title('Y direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,4)
plot(Ud_ref(11,:),'--'); hold on;
plot(Ut_ref(11,:))
title('Y direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off
%% X direction foot forces
figure()
subplot(2,2,1)
plot(Ud_ref(1,:),'--'); hold on;
plot(Ut_ref(1,:))
title('X direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,2)
plot(Ud_ref(4,:),'--'); hold on;
plot(Ut_ref(4,:))
title('Z direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,3)
plot(Ud_ref(4,:),'--'); hold on;
plot(Ut_ref(4,:))
title('X direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off

subplot(2,2,4)
plot(Ud_ref(10,:),'--'); hold on;
plot(Ut_ref(10,:))
title('X direction foot force of front left foot')
xlabel('time (ms)')
ylabel('force (N)')
legend('reference','MPC')
hold off
