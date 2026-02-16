
%% log-log plots

%% Dissipative time constant (Tau)

clear tau y
num_scales = numel(Tau);

for hh = 1:num_scales
    tau = cell2mat(Tau{hh}(:));
    tau (tau==0) = [];
    y(hh) = mean(tau);
end
figure

mdl = fitlm(log(2:num_scales),log(y(2:end)),'RobustOpts','on');
h = plot(mdl);
h(1).Marker = 'o';
h(1).MarkerFaceColor = 'k';
h(1).MarkerEdgeColor = 'k';
slope = round(mdl.Coefficients.Estimate(2),1);
intercept = round(mdl.Coefficients.Estimate(1),1);
h(2).DisplayName = ['Linear fit: y = ' ,num2str(slope),'x +' num2str(intercept)];
h(3).DisplayName = '95% conf. bounds';
grid on

title({'Avg Timescale of Partitions Across Scales',' '})
xlabel ('Scale');ylabel('Log (Avg Timescale)')
xticks(log(1:num_scales-1))
xticklabels(num2cell(1:num_scales-1))
ax = gca;ax.FontWeight = 'bold';ax.FontSize=14;
saveas(gcf, 'Avg Timescale.png');

%% KE (kinetic energy)

clear ke y
num_scales = numel(KE);
figure
for hh = 1:num_scales
    ke = cell2mat(KE{hh}(:));
    ke (ke==0) = [];
    ke = abs(ke);
    y(hh) = mean(ke);
end

mdl = fitlm(log(2:num_scales),log(y(2:end)),'RobustOpts','on');
h = plot(mdl);
h(1).Marker = 'o';
h(1).MarkerFaceColor = 'k';
h(1).MarkerEdgeColor = 'k';
slope = round(mdl.Coefficients.Estimate(2),1);
intercept = round(mdl.Coefficients.Estimate(1),1);
h(2).DisplayName = ['Linear fit: y = ' ,num2str(slope),'x ' num2str(intercept)];
h(3).DisplayName = '95% conf. bounds';
grid on

title({'Avg Kinetic Energy of Partitions Across Scales',' '})
xlabel ('Scale');ylabel('Log (Avg KE)')
xticks(log(1:num_scales-1))
xticklabels(num2cell(1:num_scales-1))
ax = gca;ax.FontWeight = 'bold';ax.FontSize=14;
saveas(gcf, 'Avg KE.png');

%% Cardinality/size of partitions

clear card y
num_scales = numel(Num_nodes);
figure
for hh = 1:num_scales
    card = Num_nodes{hh}(:);
    y(hh) = log(mean(card));
end

mdl = fitlm(log(2:num_scales-1),y(2:end-1),'RobustOpts','on');
h = plot(mdl);
h(1).Marker = 'o';
h(1).MarkerFaceColor = 'k';
h(1).MarkerEdgeColor = 'k';
slope = round(mdl.Coefficients.Estimate(2),1);
intercept = round(mdl.Coefficients.Estimate(1),1);
h(2).DisplayName = ['Linear fit: y = ' ,num2str(slope),'x +' num2str(intercept)];
h(3).DisplayName = '95% conf. bounds';
grid on

title({'Avg Cardinality of Partitions Across Scales',' '})
xlabel ('Scale');ylabel('Log (Avg Partition Size)')
xticks(log(1:num_scales-1))
xticklabels(num2cell(1:num_scales-1))
ax = gca;ax.FontWeight = 'bold';ax.FontSize=14;
ylim([0,8])
saveas(gcf, 'Avg Cardinality.png');

%% Solenoidal Freq

clear freq y
num_scales = numel(Freq);

for hh = 1:num_scales
    freq = cell2mat(Freq{hh}(:));
    freq (freq==0) = [];
    y(hh) = mean(freq);
end
figure
mdl = fitlm(log(2:num_scales),log(y(2:end)),'RobustOpts','on');
h = plot(mdl);
h(1).Marker = 'o';
h(1).MarkerFaceColor = 'k';
h(1).MarkerEdgeColor = 'k';
slope = round(mdl.Coefficients.Estimate(2),1);
intercept = round(mdl.Coefficients.Estimate(1),1);
h(2).DisplayName = ['Linear fit: y = ' ,num2str(slope),'x ' num2str(intercept)];
h(3).DisplayName = '95% conf. bounds';
grid on

title({'Avg Frequency of Partitions Across Scales',' '})
xlabel ('Scale');ylabel('Log (Avg frequency)')
xticks(log(1:num_scales-1))
xticklabels(num2cell(1:num_scales-1))
ax = gca;ax.FontWeight = 'bold';ax.FontSize=14;
saveas(gcf, 'Avg Freq.png');
