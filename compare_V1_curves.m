clear all 
close all

old_preCurves = load("1v_v.mat");
old_postCurves = load("no_e24_v.mat");

old_preCurves.V1_out(17,:) = []; 
nVerts = 17;

preCurves = cell(size(old_preCurves.V1_out));
postCurves = cell(size(old_postCurves.V1_out));

for i = 1:length(old_preCurves.V1_out)
    curve_pre = old_preCurves.V1_out{i};        
    preCurves{i} = curve_pre(1:10:end, :); 
end

for j = 1:length(old_postCurves.V1_out)
    curve_post = old_postCurves.V1_out{j}; 
    postCurves{j} = curve_post(1:10:end, :); 

end

% for e1: delta AUC = 29.6019
% for e2: delta AUC = 1.4817
% for e3: delta AUC = 0.0685
% for e4: delta AUC = -1.2850
% for e5: delta AUC = -0.5398
% for e6: delta AUC = -0.2979
% for e7: delta AUC = 0.2384
% for e8: delta AUC = 0.3245
% for e9: delta AUC = -0.6173
% for e10: delta AUC = -0.4837
% for e11: delta AUC = 0.6623
% for e12: delta AUC = 14.5113
% for e13: delta AUC = 2.3562
% for e14: delta AUC = 1.2526
% for e15: delta AUC = 0.1099
% for e16: delta AUC = 0.0486
% for e17: delta AUC = 0.0888
% for e18: delta AUC = 0.5219
% for e19: delta AUC = 0.5766
% for e20: delta AUC = 13.4281
% for e21: delta AUC = -14.8271
% for e22: delta AUC = 0.2857
% for e23: delta AUC = 0.3143
% for e24: delta AUC = 0.2580

% for e21e10: delta AUC = -159.6737 (reduced)
% for e20e10: delta AUC = 140.2670 (reduced)

maxPre = zeros(nVerts,1);
maxPost = zeros(nVerts,1);
aucPre = zeros(nVerts,1);
aucPost = zeros(nVerts,1);


for i = 1:nVerts
    y_pre = preCurves{i}(:,1)*1000000; % Infection density before removal
    y_post = postCurves{i}(:,1)*1000000; % Infection density after removal

    maxPre(i) = max(y_pre);
    maxPost(i) = max(y_post);

    aucPre(i) = trapz(y_pre);
    aucPost(i) = trapz(y_post);
end

% Now you can compare them
deltaMax = maxPost - maxPre;
deltaAUC = aucPost - aucPre;

sum(deltaAUC) % if negative overall, total infection across the network decreased

% figure;
% subplot(1,2,1);
% bar(deltaMax);
% xlabel('Vertex');
% ylabel('Change in Max Infection Density');
% title('Effect of Removing Edge (Max)');
% 
% subplot(1,2,2);
% bar(deltaAUC);
% xlabel('Vertex');
% ylabel('Change in Total Infection (AUC)');
% title('Effect of Removing Edge (AUC)');