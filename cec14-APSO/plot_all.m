% load('func1_org_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func1_PBS_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func1_PID_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func2_org_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func2_PBS_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func2_PID_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func1_org_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func1_PBS_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func1_PID_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func4_org_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func4_PBS_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func4_PID_GPSO.mat')
% plot(log(hist_v_all)); hold on;
% load('func2_org_GPSO.mat')
% plot(log(hist_v_all),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v_all(1:400:15000)),'O'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% load('func2_PBS_GPSO.mat')
% plot(log(hist_v_all),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v_all(1:400:15000)),'+'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% load('func2_PID_GPSO.mat')
% plot(log(hist_v_all),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v_all(1:400:15000)),'d'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% load('func2_PID_GPSO_V2.mat')
% plot(log(hist_v_all),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v_all(1:400:15000)),'>'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% legend('Standard PSO','PBS-PSO','DPE-PSO-V1','DPE-PSO-V2')
% xlabel('Iteration')
% ylabel('log(Loss Value)')
load('func10_org_APSO.mat')
plot(log(hist_v),'LineWidth', 1); hold on;
h = plot(1:400:15000,log(hist_v(1:400:15000)),'O'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

load('func10_PBS_APSO.mat')
plot(log(hist_v),'LineWidth', 1); hold on;
h = plot(1:400:15000,log(hist_v(1:400:15000)),'d'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

load('func10_PIDM_APSO.mat')
plot(log(hist_v),'LineWidth', 1); hold on;
h = plot(1:400:15000,log(hist_v(1:400:15000)),'^'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

load('func10_PID_APSO.mat')
plot(log(hist_v),'LineWidth', 1); hold on;
h = plot(1:400:15000,log(hist_v(1:400:15000)),'+'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

% load('func2_PID_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'+'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
