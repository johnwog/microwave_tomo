clear all
tic
%-------------Voltage definition--------------%
V0=100;
%-------------Rectangle lengths---------------%
A = 5;
B = 4;
C = 2;
D = 1;
%------------GEOMETRY DEFINITION--------------%
gd = [3 4 0 0+A 0+A 0 0 0 0+B 0+B;
    3 4 1.5 1.5+C 1.5+C 1.5 2 2 2+D 2+D]';

%-------------Decompesed Geometry-------------%
d1 = decsg(gd,'R1-R2',[abs('R1')' abs('R2')']);  %decompesed geometry matrix
[p,e,t] = initmesh(d1); %mesh matrices

%----------------Numbers----------------------% 
Nnodes = size(p,2); %number of nodes
Nelements = size(t,2); %number of elements
Nedges = size(e,2); %number of edges

%----------------Mesh plot--------------------%
pdeplot(p,e,t);axis tight;axis equal;

%---------DIRICHLET NODES LOCALISATION--------%
node_id = ones(Nnodes,1); %node identity matrix if unknown id=1 if Dirichlet=0
X0 = zeros(Nnodes,1); %voltage of nodes V0 if inner surface 0 if outer
for id=1:Nedges
    left = e(6,id); %left area of edge 
    right = e(7,id); %right area of edge
    
    if(left==0||right==0)%check if node belongs to bound area
        n1 = e(1,id); %number of first node of the edge
        n2 = e(2,id); %number of second node of the edge
        node_id(n1) = 0;
        node_id(n2) = 0;
    
        %check in which surface the edge element belongs S1 or S2 for suitable
        %X0-voltage definement 0 or V0
%         if((p(1,id)>gd(3,2)/2)&&(p(1,id)<((gd(4,2)+gd(4,1))/2))&&(p(2,id)>(gd(7,2)/2))&&(p(2,id)<(gd(10,2)+gd(10,1))/2))
%             X0(id) = V0;
%         else
%             X0(id) = 0;
%         end g                GIA KAPOIO LOGO DOULEVEI KAI ME X0(id) ENO
%         THELEI X0(n1) X0(n2)
        if((p(1,n1)>gd(3,2)/2)&&(p(1,n1)<((gd(4,2)+gd(4,1))/2))&&(p(2,n1)>(gd(7,2)/2))&&(p(2,n1)<(gd(10,2)+gd(10,1))/2))
            X0(n1) = V0;
        else
            X0(n1) = 0;
        end
        if((p(1,n2)>gd(3,2)/2)&&(p(1,n2)<((gd(4,2)+gd(4,1))/2))&&(p(2,n2)>(gd(7,2)/2))&&(p(2,n2)<(gd(10,2)+gd(10,1))/2))
            X0(n2) = V0;
        else
            X0(n2) = 0;
        end
    end
end
figure%test plot
pdeplot(p,e,t,'xydata',X0,'zdata',X0,'mesh','on')
%-------------UNKNOWN NUMERATION--------------%
index = zeros(Nnodes,1);%numeration of nodes as unknowns

Nunknowns = 0; %number of unknowns

for in=1:Nnodes
    if(node_id(in)==1) %check if node is an unknown
        Nunknowns = Nunknowns+1;%number of unknowns definement
        index(in) = Nunknowns;
    else
        index(in) = 0;
    end
end

%-------STIFFNESS MATRIX FORMULATION---------%
S = spalloc(Nunknowns,Nunknowns,8*Nunknowns);
B = zeros(Nunknowns,1);

for ie=1:Nelements
    n(1:3) = t(1:3,ie); %nodes 1-3 of element ie
    region = t(4,ie); %region of element ie
    x(1:3) = p(1,n(1:3)); %x-cord of nodes 1-3 of element ie
    y(1:3) = p(2,n(1:3)); %y-cord of nodes 1-3 of element ie
    
    %---area of element ie computation---%
    D = det([ones(3,1) x' y']);
    Ae=abs(D)/2;
    %-------cartesian-simplex------------%
    a = [(x(2)*y(3)-x(3)*y(2))/D, (x(3)*y(1)-x(1)*y(3))/D, (x(1)*y(2)-x(2)*y(1))/D];
    b = [(y(2)-y(3))/D, (y(3)-y(1))/D, (y(1)-y(2))/D];
    c = [(x(3)-x(2))/D, (x(1)-x(3))/D, (x(2)-x(1))/D];
    
    Se=zeros(3,3);
    for k=1:3
        for l=1:3
            Se(k,l) = (b(k)*b(l)+c(k)*c(l))*Ae;%compute Spq parameter
            if(node_id(n(k))==1)%check if p node is an unknown
                if(node_id(n(l))==1)%check if q node is an unknown
                    %if both p,q are unknown update stiffens matrix
                    S(index(n(k)),index(n(l))) = S(index(n(k)),index(n(l)))+Se(k,l);
                else
                    %if node q is a known update result matrix
                    B(index(n(k)),1) = B(index(n(k)),1)-Se(k,l)*X0(n(l));
                end
            end
        end
    end
end

%-------Direct Solving------------%
X=S\B;
for in=1:Nnodes
    if(index(in)~=0)
        X0(in) = X(index(in));
   end
end
figure
pdeplot(p,e,t,'xydata',X0,'zdata',X0,'mesh','on')
title('Solution')
xlabel('x');ylabel('y');zlabel('z');
toc