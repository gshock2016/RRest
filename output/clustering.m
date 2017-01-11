% infM = csvread('infmatrix.csv');
infV = csvread('infvalue.csv');

sorted = sort(infV(:, 2));
cdfplot(sorted);
% plot(sorted);

% dis=pdist(infM,'hamming');
% sortedDis = sort(dis)
% plot(sortedDis*62500);
% cdfplot(sortedDis*62500);

% [cidx2,cmeans2] = kmeans(infM,15,'dist','hamming');
% [silh2,h] = silhouette(infM,cidx2,'hamming');

% Z = linkage(Y, 'single');
% T=cluster(Z,cutoff);
