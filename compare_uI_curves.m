clear all

preCurves = load("1v.mat");
postCurves = load("no_e21e10.mat");

preCurves.value(:,[21,10]) = [];

nEdges = 22;
maxPre = zeros(nEdges,1);
maxPost = zeros(nEdges,1);
aucPre = zeros(nEdges,1);
aucPost = zeros(nEdges,1);

for i = 1:nEdges
    y_pre = preCurves.value(:,i)*1000000; % Infection density before removal
    y_post = postCurves.value(:,i)*1000000; % Infection density after removal

    maxPre(i) = max(y_pre);
    maxPost(i) = max(y_post);

    aucPre(i) = trapz(y_pre);
    aucPost(i) = trapz(y_post);
end

% Now you can compare them
deltaMax = maxPost - maxPre;
deltaAUC = aucPost - aucPre;

sum(deltaAUC) % if negative overall, total infection across the network decreased

figure;
subplot(1,2,1);
bar(deltaMax);
xlabel('Edge');
ylabel('Change in Max Infection Density');
title('Effect of Removing Edge (Max)');

subplot(1,2,2);
bar(deltaAUC);
xlabel('Edge');
ylabel('Change in Total Infection (AUC)');
title('Effect of Removing Edge (AUC)');