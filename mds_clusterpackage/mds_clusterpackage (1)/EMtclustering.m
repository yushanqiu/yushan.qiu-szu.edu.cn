function [ari, NMI]=EMtclustering(no_dims,types)

switch types
    case 'g'
load gene_EMT;
A = gene_EMT';
%A = quantilenorm(A);
    case 'r'
load RBP_EMT;
A = data';
%A = quantilenorm(A);
end
nL = zeros(304,1);
nL(1:160)=1;
N = 304;

%mappedA=compute_mapping(A,'DiffusionMaps',no_dims);
%mappedA=compute_mapping(A,'LLE',no_dims);
mappedA = compute_mapping(A,'MDS',no_dims);
%tsne(mappedA,nL,2,50,30)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%自制MDS



%D = pdist(A);
%[mappedA,stress]=mdscale(D,no_dims,'criterion','metricsstress');
% [coeff, score, latent, tsquare] = pca_jiang(A);
% meA=A-ones(size(A,1),1)*mean(A);
% [ind1,ind2]=find(coeff==max(coeff));
% newA = A(:,ind1(1:24));
% mappedA = newA;
%ma=tsne(A,nL,2,50,30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% hierarchical clustering
clustTree = linkage(mappedA,'weighted','spearman');
clusters = cluster(clustTree,'MaxClust',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kmeans
%[clusters, ctrs] = kmeans(mappedA, 2);

%%%%%%%%%%%%%%%%%%%%%%
%% GMM
% [PX, Model] = gmm(mappedA, 2);
%  [~,clusters] = max(PX'); % 每一列的最大值
% clusters=clusters';

%purity = purityf(label,clusters);
ari = adjrand(nL,clusters)
F1measure = Fmeasure(nL',clusters')
NMI = Cal_NMI(nL',clusters)




% N = length(nL);
% alpha =0.5;
% A = mappedA;
% for j0 = 1 : N
%     for k0 = 1 : N
%      %  K1(j0,k0)=sum(min(A(j0,:).^2,A(k0,:).^3));
%     %  K1(j0,k0)=cos(norm(RA(j0,:)-RA(k0,:)));
%      % K1(j0,k0)=-log(1+norm(RA(j0,:)-RA(k0,:),2));
%      % K1(j0,k0)=sum(min((A(j0,:)),A(k0,:)));
%      K1(j0,k0)=exp(-norm(A(j0,:)-A(k0,:),2));
%    % K1(j0,k0)=4*sum(1./(1./A(j0,:)+1./A(k0,:)))-2*sum(A(j0,:)+A(k0,:))+1;
%  % K1(j0,k0)=2*sum(1./(1./power(A(j0,:),alpha)+1./power(A(k0,:),alpha)));
%  % K1(j0,k0)=-norm(A(j0,:)-A(k0,:));
%   %   K1(j0,k0)=1/(1+norm(A(j0,:)-A(k0,:)));
%  % K1(j0,k0)=prod(cos(1.75*sqrt(2)*(A(j0,:)-A(k0,:))).*exp(-1*(A(j0,:)-A(k0,:)).^2));
%   % K1(j0,k0)=prod(cos(1.75*sqrt(2)*(A(j0,:))).*exp(-1*(A(j0,:)).^2))*prod(cos(1.75*sqrt(2)*(A(k0,:))).*exp(-1*(A(k0,:)).^2));
%     end
% end
% 
% Dissim = 1 - K1;
% Dist=squareform(Dissim);
% clusterTree = linkage(Dist,'weighted');
% clusters3 = cluster(clusterTree,'MaxClust',2);
% ari = adjrand(nL,clusters3)
% F1measure = Fmeasure(nL',clusters3')
% NMI = Cal_NMI(nL',clusters3)
% 
% % 
% % clusters2 = kmeans(A,2);
% % %purity = purityf(label,clusters);
% % ari = adjrand(nL,clusters2)
% % F1measure = Fmeasure(nL',clusters2')
% % NMI = Cal_NMI(nL',clusters2)
% % 
% % 
% % W = exp(-squareform(pdist(A).^2./10));
% % clusters2 = spectral(W,1,2);
% % %purity = purityf(label,clusters);
% % ari = adjrand(nL,clusters2)
% % F1measure = Fmeasure(nL',clusters2')
% % NMI = Cal_NMI(nL',clusters2)
