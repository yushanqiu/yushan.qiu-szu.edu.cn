function [disdis,nns,Distance_norm,s,entropy,ari]=main(types)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For gene expression based EMT clustering, set types 'g'
%% For rbp based EMT clustering, set types 'r'

for i=1:50
    [ari(i)]=EMtclustering(i,types);
end
[disdis,entropy,Distance_norm,s]=dimselect(types);


switch types
    case 'g'
nns = s(1:end-1)-s(2:end);
ide = find(entropy==min(entropy));
plot(ari,'r*-')
hold on;
plot(24,ari(24),'go')
figure;
plot(entropy,'k-')
str = '(19,0.5795)';
text(19,1.5795,str,'EdgeColor','green','VerticalAlignment','bottom');
hold on;
plot(ide,entropy(ide),'r*');
title('Entropy')
figure;
plot(Distance_norm,'go-')
hold on;
plot(24,Distance_norm(24),'r*')
hold on;
plot(s,'g-')
hold on;
plot(24,s(24),'r*')
legend('Perturbnorm','Singluarvalue')
%subplot(2,2,3)
figure;
plot(ide-5:50,nns(ide-5:50),'r*-')
hold on;
plot(ide-5:50,disdis(ide-5:50),'ko-')
legend('Singular gap','Perturbgap')
    case 'r'
nns = s(1:end-1)-s(2:end);
ide = find(entropy==min(entropy));
figure;
plot(ari,'r*-')
hold on;
plot(14,ari(14),'go')
figure;
plot(entropy,'k-')
hold on;
plot(ide,entropy(ide),'r*')
figure;
plot(s,'g-')
hold on;
plot(14,s(14),'r*')
plot(Distance_norm,'g-')
hold on;
plot(14,Distance_norm(14),'r*')
%subplot(2,2,3)
figure;
plot(ide-5:50,nns(ide-5:50),'r*-')
hold on;
plot(ide-5:50,disdis(ide-5:50),'ko-')
legend('Singular gap','Perturbgap')
end

