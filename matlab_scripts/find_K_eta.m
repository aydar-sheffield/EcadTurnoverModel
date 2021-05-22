%This script was used to fit the parameters of the Kelving-Voigt model
%The script is used to generate plots in Figure 6
function [K, eta, sigma,Q,Q_dot] = find_K_eta(sigma_hat, Q, dt)
Q_dot = (Q(2:end)-Q(1:end-1))/dt;

Q = Q(1:end-1);
sigma_hat = sigma_hat(1:end-1);
m = length(Q);

norm_Q = norm(Q,2)^2;
norm_Q_dot = norm(Q_dot,2)^2;
QT_Q_dot = Q'*Q_dot;
QT_sigma = Q'*sigma_hat;
Q_dot_T_sigma = Q_dot'*sigma_hat;
eta_rhs = 1/(norm_Q_dot)*(Q_dot_T_sigma-1/norm_Q*QT_sigma*QT_Q_dot);
eta_lhs = 1-1/(norm_Q*norm_Q_dot)*(QT_Q_dot)^2;
eta = eta_rhs/eta_lhs;
K = 1/norm_Q*(QT_sigma-eta*QT_Q_dot);
sigma = K*Q + eta*Q_dot;
end