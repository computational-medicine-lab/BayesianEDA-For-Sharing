function [hh1,hh2] = plot_results(data)

%% draw figures with CI
figure('units','normalized','outerposition',[0 0 1/3 1]);
fontsize = 14;
subplot(414),
tu = (0:length(data.u)-1)/data.Fsu;
if(isfield(data,'u_groundtruth') == 1)
    u_plot = data.u_groundtruth; u_plot(u_plot<=0) = NaN;
    h = stem(tu,u_plot,'r', 'linewidth',4); hold on, set(h, 'Marker', 'none'); h.Color(4) = 0.5;
end
u_plot_ = data.u; u_plot_(u_plot_<=0) = NaN;
h = stem(tu,u_plot_,'b-', 'linewidth',1.5); set(h, 'Marker', 'none'); xlim([tu(1) tu(end)])
xlabel('Time (seconds)','fontsize',fontsize,'Interpreter','latex'); ylabel('$u(t)$','fontsize',fontsize,'Interpreter','latex');
xlim([0 data.signal_duration]);
drawnow;
pause(0.1);



subplot(411),
h = plot((0:length(data.y_obs)-1)/data.Fsy, data.y_obs,'r.','linewidth',4); h.Color(4) = 0.5; hold on, 
xlabel('Time (seconds)','fontsize',fontsize,'Interpreter','latex'); ylabel('$y(t)$','fontsize',fontsize,'Interpreter','latex');
plot((0:length(data.y)-1)/data.Fsy, data.y,'k-','linewidth',1); hold on,

ylim([(min(data.y)-(max(data.y) - min(data.y))*0.1) (max(data.y)+(max(data.y) - min(data.y))*0.1)]);
title(['Participant ',num2str(data.subject)]);
% draw confidence interval with 95% significance level
data.P(1,1,1) = data.P(1,1,2); data.P(2,2,1) = data.P(2,2,2); data.P(1,2,1) = data.P(1,2,2);
data.P(2,1,1) = data.P(2,1,2);
variance = data.P(1,1,:)+data.P(2,2,:)+data.P(1,2,:)+data.P(2,1,:);
data.y_ucl = norminv(0.975, data.y, sqrt(variance(:)'));
data.y_lcl = norminv(0.025, data.y, sqrt(variance(:)'));
ty = (0:length(data.y)-1)/data.Fsy;
fill([ty fliplr(ty)], [data.y_lcl fliplr(data.y_ucl)], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
xlim([0 data.signal_duration]);
%% plot tonic component
subplot(412),
plot((0:length(data.tonic)-1)/data.Fsy,data.tonic,'r','linewidth',0.8); hold on,
if(isfield(data,'tonic_groundtruth') == 1)
    h = plot((0:length(data.tonic_groundtruth)-1)/data.Fsy, data.tonic_groundtruth,'g.','linewidth',4); hold on, h.Color(4) = 0.5;
end
xlabel('Time (seconds)','fontsize',fontsize,'Interpreter','latex'); ylabel('$y_s(t)$','fontsize',fontsize,'Interpreter','latex');

range_tonic = max(data.tonic) - min(data.tonic);
ylim([min(data.tonic)-range_tonic*0.1 max(data.tonic)+range_tonic*0.1]);
% draw confidence interval with 95% significance level
data.P(2,2,1) = data.P(2,2,2);
variance = data.P(2,2,:);
data.tonic_ucl = norminv(0.975, data.tonic, sqrt(variance(:)'));
data.tonic_lcl = norminv(0.025, data.tonic, sqrt(variance(:)'));
% plot((0:length(data.tonic)-1)/data.Fsy,data.tonic_lcl, 'r', 'linewidth', 0.8); hold on;
% plot((0:length(data.tonic)-1)/data.Fsy,data.tonic_ucl, 'r', 'linewidth', 0.8); hold on;
ty = (0:length(data.tonic)-1)/data.Fsy;
fill([ty fliplr(ty)], [data.tonic_lcl fliplr(data.tonic_ucl)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);


%% plot phasic component
subplot(413),
plot((0:length(data.phasic)-1)/data.Fsy,data.phasic,'b','linewidth',0.8); hold on,
if(isfield(data,'phasic_groundtruth') == 1)
    h = plot((0:length(data.phasic_groundtruth)-1)/data.Fsy, data.phasic_groundtruth,'g.','linewidth',4); hold on, h.Color(4) = 0.5;
end
xlabel('Time (seconds)','fontsize',fontsize,'Interpreter','latex'); ylabel('$y_p(t)$','fontsize',fontsize,'Interpreter','latex');
xlim([0 data.signal_duration]);
% draw confidence interval with 95% significance level
data.P(1,1,1) = data.P(1,1,2);
variance = data.P(1,1,:);
data.phasic_ucl = norminv(0.975, data.phasic, sqrt(variance(:)'));
data.phasic_lcl = norminv(0.025, data.phasic, sqrt(variance(:)'));
ty = (0:length(data.phasic)-1)/data.Fsy;
fill([ty fliplr(ty)], [data.phasic_lcl fliplr(data.phasic_ucl)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
%plot convergence of the parameters
% figure, subplot(3,1,1), plot((1:length(data.theta_array))',data.theta_array(1,:)); xlabel('Iteration','fontsize',fontsize), ylabel('$\tau_r$ (seconds)','Interpreter','latex','fontsize',fontsize);
% subplot(3,1,2), plot((1:length(data.theta_array))',data.theta_array(2,:)); xlabel('Iteration','fontsize',fontsize), ylabel('$\tau_d$ (seconds)','Interpreter','latex','fontsize',fontsize);
% subplot(3,1,3), plot((1:length(data.theta_array))',data.theta_array(3,:)); xlabel('Iteration','fontsize',fontsize), ylabel('$\tau_s$ (seconds)','Interpreter','latex','fontsize',fontsize);
% subplot(4,1,4), plot((1:length(data.theta_array))',data.theta_array(4,:)); xlabel('Iteration','fontsize',fontsize), ylabel('$\chi$','Interpreter','latex','fontsize',fontsize);
% subplot(6,1,2), plot((1:itermax)',data.theta_array(5,:)); xlabel('iteration'), ylabel('\chi');
hh1 = gcf;
figure;plot((1:length(data.percentage_err_theta_array))',abs(data.percentage_err_theta_array') * 100); xlabel('Iteration','fontsize',fontsize,'Interpreter','latex'), ylabel('Percentage Error (\%)','Interpreter','latex','fontsize',fontsize);
legend({'$\tau_r$','$\tau_{dp}$','$\tau_{ds}$','$\chi_1$','$\chi_2$'},'Interpreter','latex','fontsize',fontsize);
xlim([0 data.iteration]);
hh2 = gcf;
end