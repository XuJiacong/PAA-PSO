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
fc =18;
filename = "func"+fc+"_org_GPSO.mat";
load(filename)
plot(log(hist_v),'LineWidth', 2); hold on;
h = plot(1:400:15000,log(hist_v(1:400:15000)),'O'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

filename = "func"+fc+"_PBS_GPSO.mat";
load(filename)
plot(log(hist_v),'LineWidth', 2); hold on;
h = plot(1:400:15000,log(hist_v(1:400:15000)),'d'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

filename = "func"+fc+"_PIDM_GPSO.mat";
load(filename)
plot(log(hist_v),'LineWidth', 2); hold on;
h = plot(1:400:15000,log(hist_v(1:400:15000)),'+'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

filename = "func"+fc+"_PID_GPSO.mat";
load(filename)
plot(log(hist_v),'LineWidth', 2); hold on;
h = plot(1:400:15000,log(hist_v(1:400:15000)),'^'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
