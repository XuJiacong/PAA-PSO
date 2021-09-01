fc =21;
filename = "func"+fc+"_org_DMS_PSO.mat";
load(filename)
plot(log(hist_v),'LineWidth', 2); hold on;
h = plot(1:400:13901,log(hist_v(1:400:13901)),'O'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

filename = "func"+fc+"_PBS_DMS_PSO.mat";
load(filename)
plot(log(hist_v),'LineWidth', 2); hold on;
h = plot(1:400:13901,log(hist_v(1:400:13901)),'d'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

filename = "func"+fc+"_PIDM_DMS_PSO.mat";
load(filename)
plot(log(hist_v),'LineWidth', 2); hold on;
h = plot(1:400:13901,log(hist_v(1:400:13901)),'+'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

filename = "func"+fc+"_PID_DMS_PSO.mat";
load(filename)
plot(log(hist_v),'LineWidth', 2); hold on;
h = plot(1:400:13901,log(hist_v(1:400:13903)),'^'); hold on;
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

title("Function "+fc)
legend('CLPSO','PBS-CLPSO','PBS-CLPSO-V2','DPE-CLPSO')
xlabel('Iteration')
ylabel('log(Loss Value)')

% load('func18_org_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'O'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% load('func18_PBS_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'d'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% load('func18_PIDM_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'*'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% load('func18_PID_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'+'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';

% legend('CLPSO','PBS-CLPSO','DPE-PSO-V1','DPE-PSO-V2')
% xlabel('Iteration')
% ylabel('log(Loss Value)')

% load('func19_org_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'O'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% load('func19_PBS_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'d'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% load('func19_PIDM_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'*'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% load('func19_PID_CLPSO.mat')
% plot(log(hist_v),'LineWidth', 2); hold on;
% h = plot(1:400:15000,log(hist_v(1:400:15000)),'+'); hold on;
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% 
% title('Function 19')
% legend('CLPSO','PBS-CLPSO','DPE-PSO-V1','DPE-PSO-V2')
% xlabel('Iteration')
% ylabel('log(Loss Value)')