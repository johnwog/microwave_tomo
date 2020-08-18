clear;
format long g;
m=100;
%Td=abs(rand(m,m));
Td=randi([0 1000],m,m);
%dt_cur=rand(50,16);
%gt_cur=abs(rand(m,16));
gt_cur=randi([0 1000],m,16);
%dt_pre=abs(rand(m,16));
dt_pre=randi([0 1000],m,16);
%gt_pre=rand(m,16);
gt_pre=randi([0 1000],m,16);
w_pre=randi([0 1000],m,16);
dt_cur1 = zeros(m,16);
E_t = randi([0 1000],m,16);
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
dt_cur(:,is) = -gt_cur(:,is)+L(:,is);
toc

%K=M*gt_cur(:,is);
%for is=1:16
 %   K2(:,is)=((gt_cur(:,is)-gt_pre(:,is))'*Td*gt_cur(:,is));
%end
%K-K2
%dt_cur-dt_cur1
%mean(mean((dt_cur-dt_cur1).^2))
ita_s=1;
ita_d=10;
MS=rand(16,1300);
MD=rand(m,1300);
L=rand(1300,m);
x_pre=rand(m,1);
ro_pre=rand(16,16);
r_pre=rand(m,16);
a_step=zeros(16,1);
tic
for is=1:16
    A=ita_s*((MS*L*dt_cur(:,is))'*ro_pre(:,is));
     B=ita_d*((dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is)))'*Td*r_pre(:,is));
    C=(A+B);
    D=ita_s*((MS*L*dt_cur(:,is))'*(MS*L*dt_cur(:,is)));
    E=ita_d*((dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is)))'*Td*(dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is))));
    F=(D+E);
    a_step(is) = C/F;
end
toc

MSL=MS*L;
MDL=MD*L;
tic
is=1:16;
A=ita_s*((MSL*dt_cur(:,is))'*ro_pre(:,is));
A_in=diag(A);
B1=(dt_cur(:,is)-repmat(x_pre,1,16).*(MDL*dt_cur(:,is)))';
B=ita_d*(B1*Td*r_pre(:,is));
B_in=diag(B);
C=(A_in+B_in);
D=ita_s*((MSL*dt_cur(:,is))'*(MSL*dt_cur(:,is)));
E=ita_d*(B1*Td*(dt_cur(:,is)-repmat(x_pre,1,16).*(MDL*dt_cur(:,is))));
F=(diag(D)+diag(E));
a_step2 = C./F;
toc

w_cur=zeros(m,16);
tic
for is=1:16
    w_cur(:,is) = w_pre(:,is)+a_step2(is)*dt_cur(:,is);
end
toc

tic 
is=1:16;
update = a_step2.'; %a_step is of size Nsources*1 taking its transpose
update = update(ones(1,m),:); %creating a matrix in first row 1 Nnodes_d times the step size for source No 1 etc
w_cur1(:,is) = w_pre(:,is)+update.*dt_cur(:,is);
toc

tic
U=0;
for is=1:16
    U = (sparse(diag(E_t(:,is)+MD*L*w_cur(:,is))))'*Td*(sparse(diag(E_t(:,is)+MD*L*w_cur(:,is))))+U;%left matrix for contrast computation U_0*x=V_0
end
toc 

tic
qqq=reshape(E_t+MD*L*w_cur,1,[]);
r=diag(qqq);
rr=mat2cell(r,m*ones(1,16),m*ones(1,16));
d=squeeze(cat(6,rr{:}));
a=repmat(Td,[1 1 size(d,3)]);
zCell = arrayfun(@(ind) d(:,:,ind)'*a(:,:,ind)*d(:,:,ind), 1:size(d,3),'uniformOutput',false);
z2=cat(3,zCell{:});
z22=sum(z2,3);
toc%not working for sparse array matlab does not support 3D sparse 


max(max(w_cur1-w_cur))
max(max(dt_cur1-dt_cur))
max(max(a_step-a_step2))
max(max(z22-U))




% 
% q=reshape(E_t+MD*L*w_cur,1,[]);
% r=spdiags(q.',0,Nnodes_d*Nsources,Nnodes_d*Nsources);
% rr=mat2cell(r,Nnodes_d*ones(1,Nsources),Nnodes_d*ones(1,Nsources));
% rr=reshape(rr,1,[]);
% d1=mat2cell(repmat(Td,1,Nsources^2),Nnodes_d,Nnodes_d*ones(1,Nsources^2));
% d2=cellfun(@ctranspose,rr,'UniformOutput',false);
% d3=cellfun(@mtimes,d2,d1,'UniformOutput',false);
% d4=cellfun(@mtimes,d3,rr,'UniformOutput',false);
% %     d5=cell(1,Nsources^2);
% %     d5(:)={1};
% %     d6=cellfun(@mtimes,d4,d5,'UniformOutput',false);
% U=0;
% for k=1:Nsources^2
%     U=[d4{k}]+U;
% end;