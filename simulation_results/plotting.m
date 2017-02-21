%% Load simulation data

filePattern = fullfile(pwd,'*.mat');
simulationData = dir(filePattern);
for k = 1:length(simulationData)
    baseFileName = simulationData(k).name;
    baseFileName = baseFileName(1:end-4);
    data = load(baseFileName);
    v = genvarname(baseFileName, who);
    eval([v '= data.averGospa;']);
end
%%
x = '_coal_10_98';
glmb = eval(strcat('glmb_joint',x));
lmb = eval(strcat('lmb',x));
pmbm_recycle = eval(strcat('pmbm_recycle',x));
% pmb_lbp_recycle = eval(strcat('pmb_lbp_recycle',x));
pmb_murty_recycle = eval(strcat('pmb_murty_recycle',x));

figure(1)
for i = 1:4
    subplot(2,2,i)
    plot(0:100,glmb(:,i),'LineWidth',1);
%     if i==1 || i==3
%         ylim([0,150])
%     end
    hold on
    grid on
    plot(0:100,lmb(:,i),'LineWidth',1);
%     if i==1 || i==3
%         ylim([0,150])
%     end
    plot(0:100,pmbm_recycle(:,i),'LineWidth',1);
%     if i==1 || i==3
%         ylim([0,150])
%     end
%     plot(0:100,pmb_lbp_recycle(:,i),'LineWidth',1);
    plot(0:100,pmb_murty_recycle(:,i),'LineWidth',1);
%     if i==1 || i==3
%         ylim([0,150])
%     end
    lgd = legend('\delta-GLMB (joint)','LMB','PMBM w/ recycling','PMB (Murty) w/ recycling','Location','best');
    xlabel('time step')
    switch i
        case 1
            ylabel('Mean GOSPA (Total)');
        case 2
            ylabel('Mean GOSPA (Loc)')
        case 3
            ylabel('Mean GOSPA (Missed)')
        case 4
            ylabel('Mean GOSPA (False)')
    end
end
