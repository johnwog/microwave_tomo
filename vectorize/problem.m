clear;
format long g;
m=713;
Td=abs(randi([0 1],m,m));
%dt_cur=rand(50,16);
gt_cur=abs(randi([0 100],m,16));
dt_pre=abs(randi([0 100],m,16));
gt_pre=randi([0 100],m,16);
dt_cur=zeros(m,16);
tic
for is=1:16
   dt_cur1(:,is) = -gt_cur(:,is) +(((gt_cur(:,is)-gt_pre(:,is))'*Td*gt_cur(:,is))/(gt_pre(:,is)'*Td*gt_pre(:,is)))*dt_pre(:,is);
end
toc
tic
is=1:16;
M=(gt_cur(:,is)-gt_pre(:,is))'*Td*gt_cur(:,is);
K=(gt_pre(:,is)'*Td*gt_pre(:,is));
Z=diag(M./K);
z=Z(:,1);
%I=eye(m);
%F=repmat(z,1,m).*I(1:16,:);
L=repmat(z',m,1).*dt_pre(:,is);
dt_cur(:,is) = -gt_cur(:,is)+L;
toc

%K=M*gt_cur(:,is);
%for is=1:16
 %   K2(:,is)=((gt_cur(:,is)-gt_pre(:,is))'*Td*gt_cur(:,is));
%end
%K-K2
%dt_cur-dt_cur1
mean(mean((dt_cur-dt_cur1).^2))