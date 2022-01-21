% Code for the UAV trajectory optimization problem (P9) in the paper 
% 'Energy and Age Pareto Optimal Trajectories in UAV-assisted Wireless Data
% Collection' .

% Author: Yuan Liao (yuan.liao@kcl.ac.uk)
% Time: 2022/1/21

% The specific technique details, including the Path Discretization and 
% successive convex optimization (SCA), please refer to the paper:
% Zeng, Yong, Jie Xu, and Rui Zhang. "Energy minimization for wireless 
% communication with rotary-wing UAV." IEEE Transactions on Wireless 
% Communications 18.4 (2019): 2329-2345.


% The problem solved by this code is similar as the problem (P4) in the
% above reference but just add an extra linear term (average AoI) to the
% objective function.

%% load demo test data
load('testdata_cvx.mat');

%% communication and energy parameters
% communication parameters
load('comm_parameters.mat');
gamma = (P_t*rou_0)/(sigma_2);

% energy parameters
load('energy_parameters.mat');

%% Pre-define the  M and \delta_max
delta_max = 1;
M = 2*Connect_radius/delta_max;

%% initialize the variables
V_ini = 18;

% initialize the waypoints
q1_x = linspace(tem_ini_point(1),user_location(1),(M+2)/2);
q1_y = linspace(tem_ini_point(2),user_location(2),(M+2)/2);
q1 = [q1_x;q1_y]';
q2_x = linspace(user_location(1),tem_fin_point(1),(M+2)/2);
q2_y = linspace(user_location(2),tem_fin_point(2),(M+2)/2);
q2 = [q2_x;q2_y]';
q_last = [q1;q2];

% initialize the time
T_last = zeros(M+1,1);
for m = 1:M+1
    T_last(m) = norm(q_last(m+1,:)-q_last(m,:))/V_ini;
end
% the time when hovering and communication 
T_last((M+2)/2) = D/(B*log2(1+(P_t*rou_0)/(sigma_2*(H^2))));

% initialize z
Z_last = zeros(M+1,1);
for m = 1:M+1
    Z_last(m) = sqrt(...
        sqrt(T_last(m)^4 + norm(q_last(m+1,:)-q_last(m,:))^4/(4*v_0^4))...
        - norm(q_last(m+1,:)-q_last(m,:))^2/(2*v_0^2));
end

% initialize A
A_last = zeros(M+1,1);
for m = 1:M+1
    A_last(m) = sqrt(T_last(m)*...
        log2(1+gamma/(H^2 + norm(q_last(m,:) - user_location)^2)));
end

%% optimizing by cvx

max_ite = 10;
last_opt = 0;

for ite = 1:max_ite
    % calculate beta and current rate for the convex lower bound of R
    beta = zeros(M+1,1);
    current_rate = zeros(M+1,1);
    for m = 1:(M+1)
        beta(m) = (log2(exp(1))*gamma)/...
            ((H^2 + norm(q_last(m,:)-user_location)^2)*...
            ((H^2 + norm(q_last(m,:)-user_location)^2+gamma)));
        current_rate(m) = log2(1+gamma/(H^2 + norm(q_last(m,:)-user_location)^2));
    end
    
    cvx_solver mosek

    cvx_begin
        variable q(M+2,2)
        variable T(M+1,1) nonnegative
        variable Z(M+1,1) nonnegative %slack variable for non-convex
        variable A(M+1,1) nonnegative %slack variable for non-convex
        variable Y(M+1,1) nonnegative %slack variable for convex term 'delta^3/T^2'
        variable S(M+1,1) %slack variable for model the constraints 'Y >= delta^3/T^2'
        variable S2(M+1,1) nonnegative %slack variable for convex term 'T^4/Z^2'
        variable delta(M+1,1) %delta_m = ||q_{m+1}-q_m|| 

        expression energy(M+1)

        for m = 1:M+1
            energy(m) = P_0*(T(m) + (3/U_tip^2)*quad_over_lin(delta(m),T(m))) ...
                + P_i * Z(m)...
                + 0.5 * d0rousA * Y(m);
        end

        minimize AoI_para*sum(T) + Energy_para*sum(energy)

        subject to

            %\sum_{m=0}^M ({A_m^l}^2 + 2A_m^l(A_m-A_m^l)) >= D/B
            sum(A_last.^2 + 2*A_last.*(A-A_last)) >= D/B;
            
            % T_m^4/z_m^2 <= {z_m^l}^2 + 2*z_m^l*(z_m-z_m^l) 
            % - || q_{m+1}^l - q_{m}^l ||^2/v_0^2 
            % + (2/v_0^2)*(q_{m+1}^l - q_{m}^l)^T(q_{m+1} - q_{m})
            for m = 1:M+1
                pow_p(S2(m),2) <= ...
                    Z_last(m)^2 + 2*Z_last(m)*(Z(m)-Z_last(m))...
                    - norm(q_last(m+1,:)-q_last(m,:))^2/(v_0^2)...
                    + (2/v_0^2).*(q_last(m+1,:)-q_last(m,:))'*(q(m+1,:)-q(m,:));
            end
            for m = 1:M+1
                S2(m) >= quad_over_lin(T(m),Z(m));
            end
            
            % A_m^2/T_m <= currentrate - beta_m*(||q_m - w|| - ||q_m^l - w||)
            for m = 1:M+1
                quad_over_lin(A(m),T(m)) <= current_rate(m) - beta(m)*...
                    (sum_square_abs(q(m,:)-user_location) - norm(q_last(m,:)-user_location)^2);
            end
            
            % q_1 = q^I; q_{M+2} = q^F;
            q(1,:) == tem_ini_point;
            q(M+2,:) == tem_fin_point;
            
            % ||q_{m+1} - q_{m}|| <= min(delta_max, T_m*V_max)
            for m = 1:M+1
                norm(q(m+1,:)-q(m,:)) <= delta_max;
                norm(q(m+1,:)-q(m,:)) <= T(m)*V_max;
                norm(q(m+1,:)-q(m,:)) <= delta(m);
            end
            
            % ||q_m - w|| <= d^th
            for m = 1:M+2
                norm(q(m,:)-user_location) <= Connect_radius+0.1;
            end

            % Y(m) >= ||q(m+1,:)-q(m,:)||^3/T(m)^2
            for m = 1:M+1
                {delta(m),S(m),Y(m)} == rotated_lorentz(1);
                {S(m),Y(m),delta(m)} == rotated_lorentz(1);
            end

    cvx_end
    
    % update the variables
    q_last = q;
    T_last = T;
    Z_last = Z;
    A_last = A;
    
    % stopping condition
    if norm(cvx_optval-last_opt)/cvx_optval <= 0.1
        break
    else
        last_opt = cvx_optval;
    end
end

%% Plot figures
figure(1)
% plot SN node
scatter(user_location(1),user_location(2),60,'ok');

% plot connection region
hold on
viscircles(user_location,Connect_radius,...
        'Color','k','LineStyle','--','LineWidth',1.5);
    
% plot initial and final points
hold on
scatter([tem_ini_point(1);tem_fin_point(1)],...
    [tem_ini_point(2);tem_fin_point(2)],100,'pk','filled');

% plot initial path
hold on
plot([tem_ini_point(1);user_location(1);tem_fin_point(1)],...
    [tem_ini_point(2);user_location(2);tem_fin_point(2)],...
    'Color',[0 0.4470 0.7410],'LineWidth',1.5);

% plot optimized trajectory
hold on
plot(q(:,1),q(:,2),...
        'Color','r','LineWidth',1);

legend('Sensor node','Initial and final points','Initial path by solving (P1)',...
    'General trajectory by solving (P9)','fontsize',10);
