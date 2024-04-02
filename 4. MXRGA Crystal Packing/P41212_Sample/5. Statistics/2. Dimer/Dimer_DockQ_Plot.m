load('Dimer_RMSD_Statistic.mat')

%%%% RMSD95 Assessment
Table_1=nan(size(Sample,2),5);
for r=1:size(Sample,2)
    for c=1:5
        try
            Table_1(r,c)=min(Sample(r).RMSD_95(1:c*5));
        catch
            Table_1(r,c)=min(Sample(r).RMSD_95(1:length(Sample(r).RMSD_95)));
        end
    end
end

%%%% DockQ Assessment
Table_2=nan(size(Sample,2),5);
for r=1:size(Sample,2)
    for c=1:5
        try
            Table_2(r,c)=max(Sample(r).DockQ(1:c*5));
        catch
            Table_2(r,c)=max(Sample(r).DockQ(1:length(Sample(r).DockQ)));
        end
    end
end

%

%%%% Data & Group
% Top 25
c=5;

% RMSD Accessment
P0=find(Table_1(:,c)>5|isnan(Table_1(:,c)));        % Incorrect
P1=find(Table_1(:,c)<=5);                           % Acceptable
Idx_1=[P0;P1];

G0=0*ones(size(P0));
G1=1*ones(size(P1));

% DockQ Accessment
P2=Table_2(Table_2(:,c)<0.23|isnan(Table_2(:,c)),5);  % Incorrect
P3=Table_2(Table_2(:,c)>=0.23,5);                     % Acceptable

G2=2*ones(size(P2));
G3=3*ones(size(P3));

Data=[Table_2(Idx_1,5); P2; P3];
Class=[G0;G1;G2;G3];

%%%% Violin Plot
%%% RMSD95 Assessment
% Data & Group
Group=[];
for g=0:3
    Group=[Group; repmat({['G', num2str(g)]},sum(Class==g),1)];
end

% Violin Plot
GroupList={'G0','G1','G2','G3'};
Colorlist = [0         0.4470    0.7410;
             0.8500    0.3250    0.0980;
             0         0.4470    0.7410;
             0.8500    0.3250    0.0980];

for Ind=1:length(GroupList)
    vs = Violin({Data(strcmp(Group, GroupList{Ind}))},...
        Ind+0.5,...
        'EdgeColor', [0.65 0.65 0.65],...
        'HalfViolin','full',...           % right,left, full
        'QuartileStyle','shadow',...      % shadow, boxplot, none
        'DataStyle', 'none',...           % histogram, scatter, none
        'ShowData', true,...
        'ShowBox',  true,...
        'ShowMean', false,...
        'ShowMedian',  false,...
        'ShowNotches', false,...
        'ViolinColor', {Colorlist(Ind,:)});
end

xlim([1,5])
xticks([1.5:1:5])
xticklabels({'Incorrect', 'Acceptable', 'Incorrect', 'Acceptable'})

ylim([0,1])
yticks([0:0.1:1])

% Median
% 0.0965 0.4130 0.0901 0.4065
[median(Data(Class==0),'omitnan'),...
 median(Data(Class==1),'omitnan'),...
 median(Data(Class==2),'omitnan'),...
 median(Data(Class==3),'omitnan')]