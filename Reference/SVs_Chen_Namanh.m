%% Simulation of phillip chen 
clc;clear;close all;
%% Time interval and simulation time
Step = 0.001;T_end = 60;
t = 0:Step:T_end;
%% Variables
eta_d = cell(1,size(t,2));
eta = cell(1,size(t,2));
nu = cell(1,size(t,2));
upsilon = cell(1,size(t,2));
z_eta = cell(1,size(t,2));
z_upsilon = cell(1,size(t,2));

W_c1 = cell(1,size(t,2));
W_c2 = cell(1,size(t,2));
W_a1 = cell(1,size(t,2));
W_a2 = cell(1,size(t,2));

W_c1_norm = cell(1,size(t,2));
W_c2_norm = cell(1,size(t,2));
W_a1_norm = cell(1,size(t,2));
W_a2_norm = cell(1,size(t,2));
 
V_upsilon = cell(1,size(t,2));
V_eta = cell(1,size(t,2));
grad_V_upsilon = cell(1,size(t,2));
grad_V_eta = cell(1,size(t,2));
alpha_hat = cell(1,size(t,2));
tau = cell(1,size(t,2));

%% Parameters
gamma_c1 = 1;gamma_c2 = 0.01; % critic learning rate
gamma_a1 = 10;gamma_a2 = 4; % actor learning rate
beta_eta = 10; beta_upsilon = 14; % important parameter, determine system stability
n_eta = 12; n_upsilon = 12; % number of NN nodes
phi_eta = 2; phi_upsilon = 2;
bound = 6;

%% Initial Conditions
eta{1} = [0.5;0.1;0.2];
nu{1} = [0.1;0.2;0.3];
W_c1{1} = 0.01*rand(n_eta,1);
W_c2{1} = 0.01*rand(n_upsilon,1);
W_a1{1} = 0.01*rand(n_eta,1);
W_a2{1} = 0.01*rand(n_upsilon,1);

%% System simulation
for i=1:size(t,2)
    if i<= 12000
        noise = 2*(sin(t(i))^2*cos(t(i))+sin(2*t(i))^2*cos(0.1*t(i))+sin(-1.2*t(i))^2*cos(0.5*t(i))+sin(t(i))^5);
    else
        noise = 0;
    end
    %% surface vessel model
    eta_x = eta{i}(1);eta_y = eta{i}(2);eta_z = eta{i}(3);
    J = [cos(eta_z) -sin(eta_z) 0;sin(eta_z) cos(eta_z) 0;0 0 1];
    upsilon{i} = J*nu{i};
    nu_x = nu{i}(1);nu_y = nu{i}(2);nu_z = nu{i}(3);
    M = [20 0 0;0 19 0.72;0 0.72 2.7];
    C = [0 0 -19*nu_y-0.72*nu_z;0 0 20*nu_x;19*nu_y+0.72*nu_z -20*nu_x 0];
    D = [0.72+1.3*abs(nu_x)+5.8*nu_x^2 0 0;0 0.86+36*abs(nu_y)+3*abs(nu_z) -0.1-2*abs(nu_y)+2*abs(nu_z);0 -0.1-5*abs(nu_y)+3*abs(nu_z) 6+4*abs(nu_y)+4*abs(nu_z)];
    grad_J = [-sin(eta_z)*upsilon{i}(3) -cos(eta_z)*upsilon{i}(3) 0;cos(eta_z)*upsilon{i}(3) -sin(eta_z)*upsilon{i}(3) 0;0 0 0];
    f_X = -J*pinv(M)*C*nu{i}+(grad_J-J*pinv(M)*D)*nu{i};
    delta = [0;0.01*nu_x^2+0.5;-0.1*nu_z^3+sin(nu_y)];
    
    %% reference trajectory
      eta_d{i}    =   [10*sin(0.2*t(i));-10*cos(0.2*t(i));10*sin(0.2*t(i))];
      grad_eta_d  = [10*0.2*cos(0.2*t(i));10*0.2*sin(0.2*t(i));10*0.2*cos(0.2*t(i))]; 
%     eta_d{i} = [12*sin(0.2*t(i)+pi/2);12*sin(0.2*t(i));asin(sin(0.2*t(i)))+pi/2];
%     grad_eta_d = [2.4*cos(0.2*t(i)+pi/2);2.4*cos(0.2*t(i));0.2/sqrt(1-sin(0.2*t(i))^2)*cos(0.2*t(i))];
    
    %% Tracking error of 1st step
    z_eta{i} = eta{i}-eta_d{i};
    
    %% 1st step NN
    muy1 = [-bound:2*bound/(n_eta-1):bound;-bound:2*bound/(n_eta-1):bound;-bound/2:bound/(n_eta-1):bound/2];
    S_eta = zeros(n_eta,1);
    grad_S_eta_T = zeros(3,n_eta);
    for j = 1:n_eta
        S_eta(j) = exp(-(z_eta{i}-muy1(:,j))'*(z_eta{i}-muy1(:,j))/phi_eta^2);
        grad_S_eta_T(:,j) = -2/phi_eta^2*(z_eta{i}-muy1(:,j))*S_eta(j);
    end
    
    %% estimated controlling performance of 1st step
    V_eta{i} = beta_eta*z_eta{i}'*z_eta{i} + W_c1{i}'*S_eta;
    grad_V_eta{i} = 2*beta_eta*z_eta{i} + grad_S_eta_T*W_c1{i};
    
    %% estimated virtual control of 1st step
    alpha_hat{i} = -beta_eta*z_eta{i}-0.5*grad_S_eta_T*W_a1{i};
    sigma_eta = -grad_S_eta_T'*(beta_eta*z_eta{i}+0.5*grad_S_eta_T*W_a1{i}+grad_eta_d);
    
    %% Tracking error of 2nd step
    z_upsilon{i} = upsilon{i} - alpha_hat{i};
    
    %% 2nd step NN
    muy2 = [-bound:2*bound/(n_upsilon-1):bound;-bound:2*bound/(n_upsilon-1):bound;-bound/2:bound/(n_upsilon-1):bound/2];
    S_upsilon = zeros(n_upsilon,1);
    grad_S_upsilon_T = zeros(3,n_upsilon);
    for j = 1:n_upsilon
        S_upsilon(j) = exp(-(z_upsilon{i}-muy2(:,j))'*(z_upsilon{i}-muy2(:,j))/phi_upsilon^2);
        grad_S_upsilon_T(:,j) = -2/phi_upsilon^2*(z_upsilon{i}-muy2(:,j))*S_upsilon(j);
    end
    
    %% estimated controlling performance of 1st step
    V_upsilon{i} = beta_upsilon*z_upsilon{i}'*z_upsilon{i}+W_c2{i}'*S_upsilon;
    grad_V_upsilon{i} = 2*beta_upsilon*z_upsilon{i} + grad_S_upsilon_T*W_c2{i};
    
    %% compute alpha_hat_dot
    if i == 1
        alpha_hat_dot = 0;
    else
        alpha_hat_dot = (alpha_hat{i} - alpha_hat{i-1})/Step;
    end
    
    %% estimated virtual control of 1st step
    u = -beta_upsilon*z_upsilon{i}-0.5*grad_S_upsilon_T*W_a2{i};
    sigma_upsilon = grad_S_upsilon_T'*(f_X-alpha_hat_dot-beta_upsilon*z_upsilon{i}-0.5*grad_S_upsilon_T*W_a2{i});
    
    %% Real control input
    tau{i} = M*pinv(J)*u;

    %% Frobenius norm of NN
    W_c1_norm{i} = norm(W_c1{i},'fro');
    W_a1_norm{i} = norm(W_a1{i},'fro');
    W_c2_norm{i} = norm(W_c2{i},'fro');
    W_a2_norm{i} = norm(W_a2{i},'fro');
    
    %% Update stage
    F_eta = J*nu{i};
    F_nu = pinv(M)*(-C*nu{i}-D*nu{i}+tau{i}-delta);
    
    if i==size(t,2)
        break;
    end
    eta{i+1} = eta{i} + Step*F_eta + Step*noise;
    nu{i+1} = nu{i} + Step*F_nu + Step*noise;
    
    W_c1{i+1} = W_c1{i}+Step*(-gamma_c1/(1+sigma_eta'*sigma_eta)*sigma_eta*(sigma_eta'*W_c1{i}-(beta_eta^2-1)*z_eta{i}'*z_eta{i}-2*beta_eta*z_eta{i}'*grad_eta_d+1/4*(grad_S_eta_T*W_a1{i})'*(grad_S_eta_T*W_a1{i})));
    W_a1{i+1} = W_a1{i}+Step*(0.5*grad_S_eta_T'*z_eta{i}-gamma_a1*grad_S_eta_T'*grad_S_eta_T*W_a1{i}+gamma_c1/(4*(1+sigma_eta'*sigma_eta))*grad_S_eta_T'*grad_S_eta_T*W_a1{i}*sigma_eta'*W_c1{i});
    W_c2{i+1} = W_c2{i}+Step*(-gamma_c2/(1+sigma_upsilon'*sigma_upsilon)*sigma_upsilon*(sigma_upsilon'*W_c2{i}-(beta_upsilon^2-1)*z_upsilon{i}'*z_upsilon{i}+2*beta_upsilon*z_upsilon{i}'*(f_X-alpha_hat_dot)+0.25*W_a2{i}'*grad_S_upsilon_T'*grad_S_upsilon_T*W_a2{i}));
    W_a2{i+1} = W_a2{i}+Step*(0.5*grad_S_upsilon_T'*z_upsilon{i}-gamma_a2*grad_S_upsilon_T'*grad_S_upsilon_T*W_a2{i}+gamma_c2/(4*(1+sigma_upsilon'*sigma_upsilon))*grad_S_upsilon_T'*grad_S_upsilon_T*W_a2{i}*sigma_upsilon'*W_c2{i});
    
end

%% Plot
eta = cell2mat(eta);
eta_d = cell2mat(eta_d);
z_eta = cell2mat(z_eta);
z_upsilon = cell2mat(z_upsilon);
 
W_c1 = cell2mat(W_c1);
W_c2 = cell2mat(W_c2);
W_a1 = cell2mat(W_a1);
W_a2 = cell2mat(W_a2);

W_c1_norm = cell2mat(W_c1_norm);
W_c2_norm = cell2mat(W_c2_norm);
W_a1_norm = cell2mat(W_a1_norm);
W_a2_norm = cell2mat(W_a2_norm);
 
V_eta = cell2mat(V_eta);
V_upsilon = cell2mat(V_upsilon);
grad_V_eta = cell2mat(grad_V_eta);
grad_V_upsilon = cell2mat(grad_V_upsilon);
tau = cell2mat(tau);

figure;
subplot(2,1,1)
plot(t,W_c2);
title('W_{c2}');
subplot(2,1,2);
plot(t,W_a2);
title('W_{a2}');

figure;
subplot(2,1,1)
plot(t,W_c1);
title('W_{c1}');
subplot(2,1,2);
plot(t,W_a1);
title('W_{a1}');

figure;
plot(t,grad_V_eta(1,:));
hold on
plot(t,grad_V_eta(2,:));
plot(t,grad_V_eta(3,:));
hold off
title('grad V_{\eta}');
figure;
plot(t,grad_V_upsilon(1,:));
hold on
plot(t,grad_V_upsilon(2,:));
plot(t,grad_V_upsilon(3,:));
hold off
title('grad V_{\upsilon}');

figure;
plot(t,tau(1,:));
hold on
plot(t,tau(2,:));
plot(t,tau(3,:));
hold off
legend('Surge Torque','Sway Torque','Yaw Torque');
title('Control signal');

figure;
subplot(2,1,1);
plot(t,W_c2_norm);
title('Frobenius norm of W_{c2}');
subplot(2,1,2);
plot(t,W_a2_norm);
title('Frobenius norm of W_{a2}');

figure;
subplot(2,1,1);
plot(t,W_c1_norm);
title('Frobenius norm of W_{c1}');
subplot(2,1,2);
plot(t,W_a1_norm);
title('Frobenius norm of W_{a1}');

figure;
index_upsilon = V_upsilon>100;
V_upsilon1 = V_upsilon;
V_upsilon1(index_upsilon) = 100;
plot(t,V_upsilon1);
title('V_{\upsilon}');

figure;
index_eta = V_eta>100;
V_eta1 = V_eta;
V_eta1(index_eta) = 100;
plot(t,V_eta1);
title('V_{\eta}');

figure;
subplot(3,1,1);
plot(t,z_upsilon(1,:));
title('z_{\upsilon_1}');
subplot(3,1,2);
plot(t,z_upsilon(2,:));
title('z_{\upsilon_2}');
subplot(3,1,3);
plot(t,z_upsilon(3,:));
title('z_{\upsilon_3}');

figure;
subplot(3,1,1);
plot(t,z_eta(1,:));
title('z_{\eta_1}');
subplot(3,1,2);
plot(t,z_eta(2,:));
title('z_{\eta_2}');
subplot(3,1,3);
plot(t,z_eta(3,:));
title('z_{\eta_3}');

figure;
subplot(3,1,1);
plot(t,eta(1,:));
hold on
plot(t,eta_d(1,:));
title('X position');
hold off
subplot(3,1,2);
plot(t,eta(2,:));
hold on
plot(t,eta_d(2,:));
title('Y position');
hold off
subplot(3,1,3);
plot(t,eta(3,:));
hold on
plot(t,eta_d(3,:));
title('Heading angle');
hold off

figure;
plot(eta(1,:),eta(2,:));
hold on
plot(eta_d(1,:),eta_d(2,:));
title('Trajectory');
legend('Real Trajectory','Reference Trajectory');
hold off