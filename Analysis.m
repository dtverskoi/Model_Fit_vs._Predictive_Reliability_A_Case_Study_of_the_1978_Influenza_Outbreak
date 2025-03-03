clear all
close all

mo=2;

if mo==2
    ts=6;             % the size of the training set
elseif mo==3
    ts=9;
elseif mo==4
    ts=5;
end

% Download Data
Data1=[3 0; 8 0; 26 0; 76 0; 225 9; 298 17; 258 105; 233 162; 189 176; 128 166; 68 150; 29 85; 14 47; 4 20];
Data=[Data1(:,1); Data1(:,2); 512];

% Download the simulations (the results of the vanilla abc)
Out_final=[];
for experiment=1:5
    Out_final=[Out_final readmatrix(['Out_final' num2str(experiment) num2str(ts) '.csv'])];
end
% Coefficients
coeff=Out_final(1:4,:);

% Credible intervals
hpd_region = empirical_hpd_4d(coeff', 0.05);
CI0=hpd_region.region_samples;
CI=[min(CI0(:,1)) min(CI0(:,2)) min(CI0(:,3)) min(CI0(:,4)); max(CI0(:,1)) max(CI0(:,2)) max(CI0(:,3)) max(CI0(:,4))];

% Number of bins of a histogram:
bi=50;
% ylim
Ym=900;

% 1. Plot posterior distributions of the coefficients
figure
set(gcf, 'Position',  [200, 0, 800, 800])
tiledlayout(2,2);

nexttile
h=histogram(coeff(1,:),bi);
hold on 
plot(CI(1,1)*ones(1,Ym),1:Ym,'Color','red','LineWidth',3,'LineStyle','--')
plot(CI(2,1)*ones(1,Ym),1:Ym,'Color','red','LineWidth',3,'LineStyle','--')
hold off
set(gca,'FontSize',20)
title('$\beta$','FontSize',35,'Interpreter','Latex')
ylim([0 Ym])

nexttile
histogram(coeff(2,:),bi)
hold on 
plot(CI(1,2)*ones(1,Ym),1:Ym,'Color','red','LineWidth',3,'LineStyle','--')
plot(CI(2,2)*ones(1,Ym),1:Ym,'Color','red','LineWidth',3,'LineStyle','--')
hold off
set(gca,'FontSize',20)
title('$\gamma$','FontSize',35,'Interpreter','Latex')
ylim([0 Ym])

nexttile
histogram(coeff(3,:),bi)
hold on 
plot(CI(1,3)*ones(1,Ym),1:Ym,'Color','red','LineWidth',3,'LineStyle','--')
plot(CI(2,3)*ones(1,Ym),1:Ym,'Color','red','LineWidth',3,'LineStyle','--')
hold off
set(gca,'FontSize',20)
title('$\delta$','FontSize',35,'Interpreter','Latex')
ylim([0 Ym])

nexttile
histogram(coeff(4,:),bi)
hold on 
plot(CI(1,4)*ones(1,Ym),1:Ym,'Color','red','LineWidth',3,'LineStyle','--')
plot(CI(2,4)*ones(1,Ym),1:Ym,'Color','red','LineWidth',3,'LineStyle','--')
hold off
set(gca,'FontSize',20)
title('$\nu$','FontSize',35,'Interpreter','Latex')
ylim([0 Ym])
print(['SIRpc' num2str(mo)],'-depsc')


% Find parameters within the credible region
Out_new=NaN(size(Out_final));
for i=1:size(CI0,1)
    I = ismember(coeff',CI0(i,:),'rows');
    if sum(I)>=1
        Ji=find(I==1);
        Out_new(:,Ji)=Out_final(:,Ji);
    end
end

% Find the corresponding trajectories
Traj=Out_new(5:61,:);
CI_traj=NaN(size(Traj,1),3);
for k=1:size(Traj,1)
    Traj_k=sort(Traj(k,:));
    CI_traj(k,:)=[Traj_k(1) nanmean(Out_final(4+k,:)) Traj_k(size(CI0,1))];
end

% 2. Plot the trajectories
RM1=[rmse(CI_traj(28+ts+1:42,2),Data1(ts+1:end,1)) rmse(CI_traj(42+ts+1:56,2),Data1(ts+1:end,2)) CI_traj(57,2)]

figure
set(gcf, 'Position',  [200, 0, 1600, 400])
tiledlayout(1,4);

nexttile
plot(CI_traj(1:14,2),'LineWidth',3,'Color','green')
hold on
plot(CI_traj(1:14,1),'LineWidth',3,'Color','green','LineStyle','--')
plot(CI_traj(1:14,3),'LineWidth',3,'Color','green','LineStyle','--')
plot(ts*ones(100,1),linspace(1,800,100),'LineWidth',3,'LineStyle','-.','Color','black')
hold off
set(gca,'FontSize',20)
xlabel('t','FontSize',30)
ylabel('S','FontSize',30)
xticks([2 7 12])
xticklabels({'2','7','12'})
xlim([1 14])
if mo==2
    yl=[0 800];
end
if mo==3
    yl=[0 800];
end
if mo==4
    yl=[0 800];
end
ylim(yl)

nexttile
plot(CI_traj(15:28,2),'LineWidth',3,'Color','green')
hold on
plot(CI_traj(15:28,1),'LineWidth',3,'Color','green','LineStyle','--')
plot(CI_traj(15:28,3),'LineWidth',3,'Color','green','LineStyle','--')
plot(ts*ones(100,1),linspace(1,850,100),'LineWidth',3,'LineStyle','-.','Color','black')
hold off
set(gca,'FontSize',20)
xlabel('t','FontSize',30)
ylabel('I','FontSize',30)
xticks([2 7 12])
xticklabels({'2','7','12'})
xlim([1 14])
if mo==2
    yl=[0 400];
end
if mo==3
    yl=[0 360];
end
if mo==4
    yl=[0 600];
end
ylim(yl)

nexttile
plot(Data1(:,1),'LineWidth',3,'Color','black')
set(gca,'FontSize',20)
hold on
plot(CI_traj(29:42,2),'LineWidth',3,'Color','blue')
plot(CI_traj(29:42,1),'LineWidth',3,'Color','blue','LineStyle','--')
plot(CI_traj(29:42,3),'LineWidth',3,'Color','blue','LineStyle','--')
plot(ts*ones(100,1),linspace(1,750,100),'LineWidth',3,'LineStyle','-.','Color','black')
hold off
title(['RMSE=' num2str(RM1(1))])
xlabel('t','FontSize',30)
ylabel('B','FontSize',30)
xticks([2 7 12])
xticklabels({'2','7','12'})
xlim([1 14])
if mo==2
    yl=[0 750];
end
if mo==3
    yl=[0 370];
end
if mo==4
    yl=[0 750];
end
ylim(yl)

nexttile
plot(Data1(:,2),'LineWidth',3,'Color','black')
set(gca,'FontSize',20)
hold on
plot(CI_traj(43:56,2),'LineWidth',3,'Color','red')
plot(CI_traj(43:56,1),'LineWidth',3,'Color','red','LineStyle','--')
plot(CI_traj(43:56,3),'LineWidth',3,'Color','red','LineStyle','--')
plot(ts*ones(100,1),linspace(1,750,100),'LineWidth',3,'LineStyle','-.','Color','black')
hold off
title(['RMSE=' num2str(RM1(2))])
xlabel('t','FontSize',30)
ylabel('C','FontSize',30)
xticks([2 7 12])
xticklabels({'2','7','12'})
xlim([1 14])
if mo==2
    yl=[0 600];
end
if mo==3
    yl=[0 260];
end
if mo==4
    yl=[0 750];
end
ylim(yl)

print(['SIRp' num2str(mo)],'-depsc')



% 3. plot R0 and AR
R0 = 763*coeff(1,:)./coeff(2,:);
meanR0 = mean(R0)
AR = Out_final(end,:);
meanAR = mean(AR)

figure
set(gcf, 'Position',  [200, 0, 400, 800])
tiledlayout(2,1);

nexttile
histogram(R0,bi) 
set(gca,'FontSize',20)
title('$R_0$','FontSize',35,'Interpreter','Latex')
ylim([0 800])

nexttile
histogram(AR,bi) 
set(gca,'FontSize',20)
title('$AR$','FontSize',35,'Interpreter','Latex')
ylim([0 800])

print(['SIRR0' num2str(mo)],'-depsc')