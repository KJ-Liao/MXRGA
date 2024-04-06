% Load Test_Pool_467
load('Reduced_Feature_Condition_DISM_050310_Temp_pH_493.mat')
load('Significant_Coeff.mat')
load('AUC_1_Coeff.mat')

Test_Pool_467=Reduced_Feature(p_value,:);
Condition_DISM=reshape(Test_Pool_467(1,:)',2156,467);

clearvars -except Coeff Test_Pool_467 Condition_DISM

%%%% Find First Match in Whole Pool: 467
% Random
Rank_Random=zeros(467,2);
for r=1:467
    Rank=zeros(500,1);
    for n=1:5000
        Ind=randsample(2156,2156);
        Rank(n,1)=find(Condition_DISM(Ind,r)==0,1);
    end
    Rank_Random(r,1)=mean(Rank);
    Rank_Random(r,2)=std(Rank);
end

% Seq
Seq_DISM=reshape(Test_Pool_467(2,:)',2156,467);
Rank_Seq=zeros(467,1);
for r=1:467
    [Value, Order]=sort(Seq_DISM(:,r));
    Rank_Seq(r,1)=find(Condition_DISM(Order,r)==0,1);
end

% RMSD
RMSD_DISM=reshape(Test_Pool_467(3,:)',2156,467);
Rank_RMSD=zeros(467,1);
for r=1:467
    [Value, Order]=sort(RMSD_DISM(:,r));
    Rank_RMSD(r,1)=find(Condition_DISM(Order,r)==0,1);
end

% 3DZD
ZD_DISM=reshape(Test_Pool_467(4,:)',2156,467);
Rank_3DZD=zeros(467,1);
for r=1:467
    [Value, Order]=sort(ZD_DISM(:,r));
    Rank_3DZD(r,1)=find(Condition_DISM(Order,r)==0,1);
end

% MXRGA
MXRGA_DISM=reshape([ones(size(Test_Pool_467,2),1), Test_Pool_467(2:end,:)']*-Coeff,2156,467);
Rank_MXRGA=zeros(467,1);
for r=1:467
    [Value, Order]=sort(MXRGA_DISM(:,r));
    Rank_MXRGA(r,1)=find(Condition_DISM(Order,r)==0,1);
end

% Data & Group
Data=[Rank_MXRGA; Rank_Random(:,1); Rank_Seq; Rank_RMSD; Rank_3DZD];
Group=[];
for g=1:5
    Group=[Group; repmat({['G', num2str(g)]},r,1)];
end

% Violin Plot
GroupList={'G1','G2','G3','G4','G5'};
Colorlist = [0         0.4470    0.7410;
             0.8500    0.3250    0.0980;
             0.9290    0.6940    0.1250;
             0.4940    0.1840    0.5560;
             0.1050    0.1280    0.1040];

for Ind=1:5
    vs = Violin({Data(strcmp(Group, GroupList{Ind}))},...
        Ind+0.5,...
        'Bandwidth',5,...
        'ViolinAlpha',{0.5},...    
        'EdgeColor', [0.65 0.65 0.65],...
        'HalfViolin','full',...           % right,left, full
        'ShowData', true,...
        'DataStyle', 'none',...           % histogram, scatter, none
        'QuartileStyle','none',...        % shadow, boxplot, none
        'ShowBox', false,...
        'ShowMean', false,...
        'ShowMedian', true,...
        'ShowNotches', false,...
        'MarkerSize',15,...
        'ViolinColor', {Colorlist(Ind,:)});;
end

xlim([1,6])
xticks(1.5:1:5.5)
xticklabels({'Interface', 'Random', 'Sequence', 'RMSD', 'Jerhoud 3DZD'})

ylim([0,2200])
yticks(0:200:2200)

% median([Rank_MXRGA, Rank_Random(:,1), Rank_Seq, Rank_RMSD, Rank_3DZD])
% mean([Rank_MXRGA, Rank_Random(:,1), Rank_Seq, Rank_RMSD, Rank_3DZD])
% Median:      3, 1078.1,     1,     1,     1
% Mean:     89.6, 1010.6,  78.1, 100.6, 103.6
% p-value_vs_Random: x,0.0000,0.0000,0.0000,0.0000

%
%
%

%%%% Find the Potential Commutable Conditions Counts/Ranks in Whole Pool: 467
% Cutoff of Commutable Conditions
Cutoff = 1;

% Random
Rank_Random=zeros(467,2156);
for r=1:467
    Ordered_Condition_DISM=zeros(500,2156);
    for n=1:500
        Ind=randsample(2156,2156);
        Ordered_Condition_DISM(n,:)=cumsum(Condition_DISM(Ind,r)'<Cutoff);
    end
    Rank_Random(r,:)=mean(Ordered_Condition_DISM);
end

% Seq
Seq_DISM=reshape(Test_Pool_467(2,:)',2156,467);
Rank_Seq=zeros(467,2156);
for r=1:467
    [Value, Order]=sort(Seq_DISM(:,r));
    Ordered_Condition_DISM=Condition_DISM(Order,r)<Cutoff;
    Rank_Seq(r,:)=cumsum(Ordered_Condition_DISM);
end

% RMSD
RMSD_DISM=reshape(Test_Pool_467(3,:)',2156,467);
Rank_RMSD=zeros(467,2156);
for r=1:467
    [Value, Order]=sort(RMSD_DISM(:,r));
    Ordered_Condition_DISM=Condition_DISM(Order,r)<Cutoff;
    Rank_RMSD(r,:)=cumsum(Ordered_Condition_DISM);
end

% 3DZD
ZD_DISM=reshape(Test_Pool_467(4,:)',2156,467);
Rank_3DZD=zeros(467,2156);
for r=1:467
    [Value, Order]=sort(ZD_DISM(:,r));
    Ordered_Condition_DISM=Condition_DISM(Order,r)<Cutoff;
    Rank_3DZD(r,:)=cumsum(Ordered_Condition_DISM);
end


% MXRGA
MXRGA_DISM=reshape([ones(size(Test_Pool_467,2),1), Test_Pool_467(2:end,:)']*-Coeff,2156,467);
Rank_MXRGA=zeros(467,2156);
for r=1:467
    [Value, Order]=sort(MXRGA_DISM(:,r));
    Ordered_Condition_DISM=Condition_DISM(Order,r)<Cutoff;
    Rank_MXRGA(r,:)=cumsum(Ordered_Condition_DISM);
end

% Accumulated Plot
Cum_Data=[mean(Rank_MXRGA); mean(Rank_Random); mean(Rank_Seq); mean(Rank_RMSD); mean(Rank_3DZD)];
Colorlist = [0         0.4470    0.7410;
             0.8500    0.3250    0.0980;
             0.9290    0.6940    0.1250;
             0.4940    0.1840    0.5560;
             0.1050    0.1280    0.1040];

Cum_Data=Cum_Data(:,1:2056);
for c=1:5
    hold on
    if c==2
        plot(Cum_Data(c,:),'Color',Colorlist(c,:),'LineWidth',1.5, 'LineStyle', ':')
    else
        plot(Cum_Data(c,:),'Color',Colorlist(c,:),'LineWidth',1.5)
    end
    % Approximation
    % y=1:2156;
    % plot(y([1:10:2156,2156]), Cum_Data(c,[1:10:2156,2156]),'Color',Colorlist(c,:),'LineWidth',1)
end


xlim([0,2056])
xticks([0:300:1800,2056])

legend({'Interface-based','Random Guess','Sequence-based','RMSD-based','3DZD-based'}, ...
    'Location','southeast', ...
    'FontSize',8);

% Half-Height: find(Cum_Data(r,:)>max(Cum_Data,[],'all')/2,1)
% find(Cum_Data(r,:)>=c,1)
% Half:  727,1079, 905, 878, 934
% 1:       8, 233,   5,   7,  18
% 2:     159, 466, 220, 244, 277
% 3:     362, 699, 504, 492, 727
% 4:     576, 931, 793, 730, 779
% 5:     829,1163, 968, 976,1036

% plot(727,Cum_Data(1,727),'o','MarkerFaceColor',Colorlist(1,:))      % MXRGA
% plot(1079,Cum_Data(2,1078),'o','MarkerFaceColor',Colorlist(2,:))    % Random
% plot(905,Cum_Data(3,905),'o','MarkerFaceColor',Colorlist(3,:))      % Sequence
% plot(878,Cum_Data(4,878),'o','MarkerFaceColor',Colorlist(4,:))      % RMSD
% plot(934,Cum_Data(5,934),'o','MarkerFaceColor',Colorlist(5,:))      % 3DZD

%
%
%

%%%% Find the Potential Commutable Conditions Counts/Ranks in Whole Pool: 467 (Sequence Dissimilarity(Distance) > DS_Cutoff)
% Cutoff of Commutable Conditions
Cutoff = 1; DS_Cutoff=0.7;

% Random
Rank_Random=nan(467,2156);
for r=1:467
    Commutable_Condition_List=nan(500,2156);
    for n=1:500
        Ind=randsample(2156,2156);
        Seq_Ind=Seq_DISM(Ind,r)>DS_Cutoff;
        Ordered_Condition_DISM=Condition_DISM(Ind,r);
        Commutable_Condition_List(n,1:sum(Seq_Ind))=cumsum(Ordered_Condition_DISM(Seq_Ind)<Cutoff);
    end
    Rank_Random(r,:)=mean(Commutable_Condition_List);
end

% Seq
Seq_DISM=reshape(Test_Pool_467(2,:)',2156,467);
Rank_Seq=nan(467,2156);
for r=1:467
    [Value, Order]=sort(Seq_DISM(:,r));
    Seq_Ind=Seq_DISM(Order,r)>DS_Cutoff;
    Ordered_Condition_DISM=Condition_DISM(Order,r);
    Rank_Seq(r,1:sum(Seq_Ind))=cumsum(Ordered_Condition_DISM(Seq_Ind)<Cutoff);
end

% RMSD
RMSD_DISM=reshape(Test_Pool_467(3,:)',2156,467);
Rank_RMSD=nan(467,2156);
for r=1:467
    [Value, Order]=sort(RMSD_DISM(:,r));
    Seq_Ind=Seq_DISM(Order,r)>DS_Cutoff;
    Ordered_Condition_DISM=Condition_DISM(Order,r);
    Rank_RMSD(r,1:sum(Seq_Ind))=cumsum(Ordered_Condition_DISM(Seq_Ind)<Cutoff);
end

% 3DZD
ZD_DISM=reshape(Test_Pool_467(4,:)',2156,467);
Rank_3DZD=nan(467,2156);
for r=1:467
    [Value, Order]=sort(ZD_DISM(:,r));
    Seq_Ind=Seq_DISM(Order,r)>DS_Cutoff;
    Ordered_Condition_DISM=Condition_DISM(Order,r);
    Rank_3DZD(r,1:sum(Seq_Ind))=cumsum(Ordered_Condition_DISM(Seq_Ind)<Cutoff);
end

% MXRGA
MXRGA_DISM=reshape([ones(size(Test_Pool_467,2),1), Test_Pool_467(2:end,:)']*-Coeff,2156,467);
Rank_MXRGA=nan(467,2156);
for r=1:467
    [Value, Order]=sort(MXRGA_DISM(:,r));
    Seq_Ind=Seq_DISM(Order,r)>DS_Cutoff;
    Ordered_Condition_DISM=Condition_DISM(Order,r);
    Rank_MXRGA(r,1:sum(Seq_Ind))=cumsum(Ordered_Condition_DISM(Seq_Ind)<Cutoff);
end

% Accumulated Plot
Cum_Data=[mean(Rank_MXRGA); mean(Rank_Random); mean(Rank_Seq); mean(Rank_RMSD); mean(Rank_3DZD)];
Colorlist = [0         0.4470    0.7410;
             0.8500    0.3250    0.0980;
             0.9290    0.6940    0.1250;
             0.4940    0.1840    0.5560;
             0.1050    0.1280    0.1040];

Cum_Data=Cum_Data(:,1:2056);
for c=1:5
    hold on
    if c==2
        plot(Cum_Data(c,:),'Color',Colorlist(c,:),'LineWidth',1.5, 'LineStyle', ':')
    else
        plot(Cum_Data(c,:),'Color',Colorlist(c,:),'LineWidth',1.5)
    end
    % Approximation
    % y=1:2156;
    % plot(y([1:10:2156,2156]), Cum_Data(c,[1:10:2156,2156]),'Color',Colorlist(c,:),'LineWidth',1)
end

xlim([0,2056])
xticks([0:300:1800,2056])

legend({'Interface-based','Random Guess','Sequence-based','RMSD-based','3DZD-based'}, ...
    'Location','southeast', ...
    'FontSize',8);

% Half-Height: find(Cum_Data(r,:)>max(Cum_Data,[],'all')/2,1)
% find(Cum_Data(r,:)>=c,1)
% 0.7: 
% Half:  832,1056, 978, 983, 1047
% 1:     154, 263, 220, 249, 278
% 2:     360, 524, 507, 490, 527
% 3:     579, 786, 798, 735, 784
% 4:     826,1047, 970, 976,1040
% 5:    1100,1309,1251,1281,1307

% plot(832,Cum_Data(1,832),'o','MarkerFaceColor',Colorlist(1,:))      % MXRGA
% plot(1056,Cum_Data(2,1055),'o','MarkerFaceColor',Colorlist(2,:))    % Random
% plot(978,Cum_Data(3,978),'o','MarkerFaceColor',Colorlist(3,:))      % Sequence
% plot(983,Cum_Data(4,983),'o','MarkerFaceColor',Colorlist(4,:))      % RMSD
% plot(1047,Cum_Data(5,1047),'o','MarkerFaceColor',Colorlist(5,:))    % 3DZD