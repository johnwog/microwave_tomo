clear all
tic
format long g
%% PARAMETERS DEFIMENENT
%------------Parameters definition-------------%
e0 = 8.854e-12; %vacuum permitivity,electric constant
m0 = 4e-7*pi; %vacuum permeability,magnetic constant
c0 = 3e8; %speed of light vacuum
f = 1e10; %frequency
w = 2*pi*f; %angular frequency
l = c0/f; %wavelength
k0 = w*sqrt(e0*m0); %vacuum wavenumber
E0 = 1; %incriment wave amplitude
%% ELECTRICAL PROPERTIES DEFINEMENT

%-----Different area dielectric properties-----%
er0 = eye(3); %vacuum relative electric constant
mr0 = eye(3); %vacuum relative magnetic constant
mr0_inv = inv(mr0);

%PML parmeters
a = 1; %PML a parameter real part
b = 5.5; %PML a parameter imaginary part 
         %for thickness d=0.15l reflectance factor 3.1*e-5

%vertical PML layer
erx = [1/(a-b*1i) 0 0;0 a-b*1i 0;0 0 a-b*1i];%electic constant tensor
mrx = erx;%megnetic constant tensor
mrx_inv = inv(mrx);

%horizontal PML layer
ery = [a-b*1i 0 0;0 1/(a-b*1i) 0;0 0 a-b*1i];%electic constant tensor
mry = ery;%megnetic constant tensor
mry_inv = inv(mry);

%PML intersection
er_inter = erx*ery;
mr_inter = mrx*mry;
mr_inter_inv = inv(mr_inter);

%scatterer
er = 8; %relative electric constant
s = 0.166; %scatterer conductivity
ers = diag([er+s/(w*e0*1i) er+s/(w*e0*1i) er+s/(w*e0*1i)]);%electic constant tensor
mrs = eye(3);%megnetic constant tensor
mrs_inv = inv(mrs);
% ks=w*sqrt(ers(1,1)*e0*mrs(1,1)*m0);

%% PML AREA AND SCATTERER GEOMETRY
%------------GEOMETRY DEFINITION--------------%
%PML parmeters
d = 0.15*l; %PML thickness
length = l; %PML length

%scatterer geometry parameters
A = 0.1*l; %length 
B = 0.1*l; %width
dis = length/2-d-A/2; %distance from PML centered


%geometry matrix scatterer centered in the start of the coordinate system
% gd=[3 4 -dis-d-0.5*A -dis-0.5*A -dis-0.5*A -dis-d-0.5*A -dis-0.5*B-d -dis-0.5*B-d dis+0.5*B+d dis+0.5*B+d; %left vertical PML
%     3 4 dis+0.5*A dis+0.5*A+d dis+d+0.5*A dis+0.5*A -dis-0.5*B-d -dis-0.5*B-d dis+0.5*B+d dis+0.5*B+d; %right vertical PML
%     3 4 -dis-d-0.5*A dis+d+0.5*A dis+d+0.5*A -dis-d-0.5*A -dis-0.5*B-d -dis-0.5*B-d -dis-0.5*B -dis-0.5*B; %lower horizontal PML 
%     3 4 -dis-d-0.5*A dis+d+0.5*A dis+d+0.5*A -dis-d-0.5*A dis+0.5*B dis+0.5*B dis+d+0.5*B dis+0.5*B+d; %upper horizontal PML
%     3 4 -0.5*A-dis 0.5*A+dis 0.5*A+dis -0.5*A-dis -0.5*B-dis  -0.5*B-dis 0.5*B+dis 0.5*B+dis; %between medium
%     3 4 -0.5*A 0.5*A 0.5*A -0.5*A -0.5*B -0.5*B 0.5*B 0.5*B;]'; %scatterer 

%lower right point of geometry is the start of the coordinate system
gd=[3 4 0 0+d 0+d 0 0 0 length length; %left vertical PML
    3 4 0+length-d 0+length 0+length 0+length-d 0 0 length length; %right vertical PML
    3 4 0 0+length 0+length 0 0 0 d d; %lower horizontal PML 
    3 4 0 0+length 0+length 0 0+length-d 0+length-d 0+length 0+length; %upper horizontal PML
    3 4 d length-d length-d d d d length-d length-d; %between medium
    3 4 dis+d dis+d+A dis+d+A dis+d d+dis d+dis d+dis+B d+dis+B;]'; %scatterer 
%are centroids calculation
%cen matrix containts the centroids of the different areas in the same row
%as defined in gd matrix.The first row contains the x coord, the second the
%y coord.%%%Not needed eventually
cen=[0.5*(gd(3,1)+gd(4,1)) 0.5*(gd(3,2)+gd(4,2)) 0.5*(gd(3,3)+gd(4,3)) 0.5*(gd(3,4)+gd(4,4)) 0.5*(gd(3,5)+gd(4,5)) 0.5*(gd(3,6)+gd(4,6));
    0.5*(gd(7,1)+gd(10,1)) 0.5*(gd(7,2)+gd(10,2)) 0.5*(gd(7,3)+gd(10,3)) 0.5*(gd(7,4)+gd(10,4)) 0.5*(gd(7,5)+gd(10,5)) 0.5*(gd(7,6)+gd(10,6))];
    if(cen(1,5)==cen(1,6)) %if scatterer and between medium have the same centroid eliminate scatterer centroid
        cen(1,6)=0;
    end
    if(cen(2,5)==cen(2,6))
        cen(2,6)=0;
    end
%-------------Decompesed Geometry-------------%
%decompesed geometry matrix
d1 = decsg(gd,'R1+R2+R3+R4+R5+R6',[abs('R1')' abs('R2')' abs('R3')' abs('R4')' abs('R5')' abs('R6')']);  

%% ORIGINAL MESH DEFINEMENT
%mesh matrices 
[p,e,t] = initmesh(d1,'hmax',l/20); 
%hmax property sets no element side greater than l/12 for a satisfactory granularity

%% MESH REFINEMNET FOR A DENSER MESH INSIDE THE SCATTERER AND PML

%%-----------Mesh refinement and region definitions--------%
for rep=1:2 %number of refinements
    
    %-------Numbers of nodes,elements,edges-------% 
    Nnodes = size(p,2); %number of nodes
    Nelements = size(t,2); %number of elements
    Nedges = size(e,2); %number of edges

    %------Region of each triangle definemet------%
    %characteristic number for each area 1=between medium 2=horizontal PML
    %3=vertical PML 4=PML intersection 5=scatterer
    regions = ones(1,Nelements);
    for ie=1:Nelements
        x_cen = sum(p(1,t(1:3,ie)))/3;
        y_cen = sum(p(2,t(1:3,ie)))/3;
        if(x_cen>gd(3,6)&&x_cen<gd(4,6)&&y_cen>gd(7,6)&&y_cen<gd(9,6))
            regions(ie) = 5; %scatterer area
        elseif(x_cen>gd(4,1)&&x_cen<gd(3,6)&&y_cen>gd(9,3)&&y_cen<gd(7,4))
            regions(ie) = 1; %between medium
        elseif(x_cen>gd(3,6)&&x_cen<gd(4,6)&&y_cen>gd(9,3)&&y_cen<gd(7,6))
            regions(ie) = 1; %between medium
        elseif(x_cen>gd(3,6)&&x_cen<gd(4,6)&&y_cen>gd(9,6)&&y_cen<gd(7,4))
            regions(ie) = 1; %between medium
        elseif(x_cen>gd(4,6)&&x_cen<gd(3,2)&&y_cen>gd(9,3)&&y_cen<gd(7,4))
            regions(ie) = 1; %between medium
        elseif(x_cen>gd(4,1)&&x_cen<gd(3,2)&&y_cen>gd(7,3)&&y_cen<gd(9,3))
            regions(ie) = 2; %lower horizontal PML
        elseif(x_cen>gd(4,1)&&x_cen<gd(3,2)&&y_cen>gd(7,4)&&y_cen<gd(9,4))
            regions(ie) = 2; %upper horizontal PML
        elseif(y_cen>gd(9,3)&&y_cen<gd(7,4)&&x_cen>gd(3,1)&&x_cen<gd(4,1))
            regions(ie) = 3; %left vertical PML
        elseif(y_cen>gd(9,3)&&y_cen<gd(7,4)&&x_cen>gd(3,2)&&x_cen<gd(4,2))
            regions(ie) = 3; %right vertical PML
        elseif(x_cen>gd(3,1)&&x_cen<gd(4,1)&&y_cen>gd(7,3)&&y_cen<gd(9,3))
            regions(ie) = 4; %downleft PML intersection
        elseif(x_cen>gd(3,2)&&x_cen<gd(4,2)&&y_cen>gd(7,3)&&y_cen<gd(9,3))
            regions(ie) = 4; %downright PML intersection
        elseif(x_cen>gd(3,1)&&x_cen<gd(4,1)&&y_cen>gd(7,4)&&y_cen<gd(9,4))
            regions(ie) = 4; %upperleft PML intersection
        elseif(x_cen>gd(3,2)&&x_cen<gd(4,2)&&y_cen>gd(7,4)&&y_cen<gd(9,4))
            regions(ie) = 4; %upperight PML intersection
        end
    end
    it=[find(regions==2) find(regions==3) find(regions==4) find(regions==5)];

    [p e t] = refinemesh(d1,p,e,t,it');

    %-------Numbers of nodes,elements,edges-------% 
    Nnodes = size(p,2); %number of nodes
    Nelements = size(t,2); %number of elements
    Nedges = size(e,2); %number of edges
    regions = ones(1,Nelements);
    for ie=1:Nelements
        x_cen = sum(p(1,t(1:3,ie)))/3;
        y_cen = sum(p(2,t(1:3,ie)))/3;
        if(x_cen>gd(3,6)&&x_cen<gd(4,6)&&y_cen>gd(7,6)&&y_cen<gd(9,6))
            regions(ie) = 5; %scatterer area
        elseif(x_cen>gd(4,1)&&x_cen<gd(3,6)&&y_cen>gd(9,3)&&y_cen<gd(7,4))
            regions(ie) = 1; %between medium
        elseif(x_cen>gd(3,6)&&x_cen<gd(4,6)&&y_cen>gd(9,3)&&y_cen<gd(7,6))
            regions(ie) = 1; %between medium
        elseif(x_cen>gd(3,6)&&x_cen<gd(4,6)&&y_cen>gd(9,6)&&y_cen<gd(7,4))
            regions(ie) = 1; %between medium
        elseif(x_cen>gd(4,6)&&x_cen<gd(3,2)&&y_cen>gd(9,3)&&y_cen<gd(7,4))
            regions(ie) = 1; %between medium
        elseif(x_cen>gd(4,1)&&x_cen<gd(3,2)&&y_cen>gd(7,3)&&y_cen<gd(9,3))
            regions(ie) = 2; %lower horizontal PML
        elseif(x_cen>gd(4,1)&&x_cen<gd(3,2)&&y_cen>gd(7,4)&&y_cen<gd(9,4))
            regions(ie) = 2; %upper horizontal PML
        elseif(y_cen>gd(9,3)&&y_cen<gd(7,4)&&x_cen>gd(3,1)&&x_cen<gd(4,1))
            regions(ie) = 3; %left vertical PML
        elseif(y_cen>gd(9,3)&&y_cen<gd(7,4)&&x_cen>gd(3,2)&&x_cen<gd(4,2))
            regions(ie) = 3; %right vertical PML
        elseif(x_cen>gd(3,1)&&x_cen<gd(4,1)&&y_cen>gd(7,3)&&y_cen<gd(9,3))
            regions(ie) = 4; %downleft PML intersection
        elseif(x_cen>gd(3,2)&&x_cen<gd(4,2)&&y_cen>gd(7,3)&&y_cen<gd(9,3))
            regions(ie) = 4; %downright PML intersection
        elseif(x_cen>gd(3,1)&&x_cen<gd(4,1)&&y_cen>gd(7,4)&&y_cen<gd(9,4))
            regions(ie) = 4; %upperleft PML intersection
        elseif(x_cen>gd(3,2)&&x_cen<gd(4,2)&&y_cen>gd(7,4)&&y_cen<gd(9,4))
            regions(ie) = 4; %upperight PML intersection
        end
    end
end

%% PLOT TEST FOR AREA CHARACTIRISATION
% %---------------Test plot 1-------------------%
% pdeplot(p,e,t);axis tight;axis equal;
% for ie=1:Nelements%test if regions ides are correct
%     x_cen = sum(p(1,t(1:3,ie)))/3;
%     y_cen = sum(p(2,t(1:3,ie)))/3;
%     if(regions(ie)==1)
%         text(x_cen,y_cen,'1')
%     elseif(regions(ie)==2)
%         text(x_cen,y_cen,'2')
%     elseif(regions(ie)==3)
%         text(x_cen,y_cen,'3')
%     elseif(regions(ie)==4)
%         text(x_cen,y_cen,'4')
%     elseif(regions(ie)==5)
%         text(x_cen,y_cen,'5')
%     end
% end
%% DIRICHLET BOUNDARY CONDITIONS
%---------DIRICHLET NODES LOCALISATION--------%
node_id = ones(Nnodes,1); %node identity matrix if unknown id=1 if Dirichlet=0

% for id=1:Nedges
%     left = e(6,id); %left area of edge
%     right = e(7,id); %right area of edge
%     
%     if(left==0||right==0)
%         n1 = e(1,id);
%         n2 = e(2,id);
%         node_id(n1) = 0;
%         node_id(n2) = 0;
%     end
% end
% figure
%% PLOT TEST FOR PROPER DIRICHLET BOUNDARY CONDITIONS
% %---------------Test plot 2-------------------%
% pdeplot(p,e,t,'xydata',real(Ezs0),'zdata',real(Ezs0),'mesh','on')

%% CREATE A NEW NUMERATION FOR NODES AS UNKNOWNS
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
%% PLOT TEST FOR UNKONWN NODES KNOWN NODES SHOULD HAVE ZERO VALUE 
% %---------------Test plot 3-------------------%
% figure
% pdeplot(p,e,t,'xydata',node_id,'zdata',node_id,'mesh','on')
% figure
% pdeplot(p,e,t,'xydata',index,'zdata',index,'mesh','on')

%% PROBLEM FOR DIFFERENT INCIDENT WAVES DEFINITION OF THE DIFFERENT INCIDENT ANGLES
% The first wave is the one on the 180deg angle and that source has 0
% degrees.The numeration increases in clockwise order

%----Multiple sources and place definement---%
Nsources = 32; %Number of sources

th = (0:(Nsources-1))*(360/Nsources); %angle of sources
theta = zeros(1,Nsources); %angle of sources from the vertical axis 

coord = ones(2,Nsources); %id matrix in which quardrant belongs the incident wave 
for is=1:Nsources
    if(th(is)<=90)
        coord(2,is) = -1;
    elseif(th(is)<=180&&th(is)>90)
        coord(1:2,is) = -1;
    elseif(th(is)<=270&&th(is)>180)
        coord(1,is) = -1;
    end
    if(th(is)<=180)
        theta(is) = abs(90-th(is));
    else
        theta(is) = abs(270-th(is));
    end
end
%% STIFFNESS AND MASS MATRIX FORMULATION FOR MULTIPLE INCIDENT WAVES-EQUATION SYSTEM CONSTRUCTION

%---------Matrices initialization---------%
Ezs0 = zeros(Nnodes,Nsources); %electric field of nodes, unknown if inner surface, 0 if outer
I = zeros(9*Nelements,1); %auxialry matrix for sparse usage the row value
J = zeros(9*Nelements,1); %auxialry matrix for sparse usage the column value 
Se = zeros(9*Nelements,1);  %%auxialry matrix for sparse usage element value
NoNonZeros = 0; %number of non zero index nodes
R = zeros(Nunknowns,Nsources);%charge matrix
S=zeros(3,3);
S_b=zeros(3,3);
T=zeros(3,3);
T_b=zeros(3,3);

%------Stiffness Matrix formulation------%
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
    %--dielectric properties definement--%
    if(regions(ie)==1)
        er_z = er0(3,3);
        mr_x = mr0_inv(1,1);
        mr_y = mr0_inv(2,2);
    elseif(regions(ie)==2)
        er_z = ery(3,3);
        mr_x = mry_inv(1,1);
        mr_y = mry_inv(2,2);
    elseif(regions(ie)==3)
        er_z = erx(3,3);
        mr_x = mrx_inv(1,1);
        mr_y = mrx_inv(2,2);
    elseif(regions(ie)==4)
        er_z = er_inter(3,3);
        mr_x = mr_inter_inv(1,1);
        mr_y = mr_inter_inv(2,2);
    elseif(regions(ie)==5)
        er_z = ers(3,3);
        mr_x = mrs_inv(1,1);
        mr_y = mrs_inv(2,2);
    end
        
    for k=1:3
        for m=1:3%Compute S T S' T' parameters for stiffness matrix calculation
            S(k,m) = (mr_x*c(k)*c(m)+mr_y*b(k)*b(m))*Ae;
            S_b(k,m) = ((mr_x-1)*c(k)*c(m)+(mr_y-1)*c(k)*c(m))*Ae;
            if(k==m)
                T(k,m) = er_z*Ae/6; 
                T_b(k,m) = (er_z-1)*Ae/6;
            else
                T(k,m) = er_z*Ae/12;
                T_b(k,m) = (er_z-1)*Ae/12;
            end
            if(node_id(n(k))==1)%check if p node is an unknown
                if(node_id(n(m))==1)%check if q node is an unknown
                    %if both p,q are unknown update stiffens matrix
                    NoNonZeros = NoNonZeros+1;
                    
                    I(NoNonZeros) = index(n(k));
                    J(NoNonZeros) = index(n(m));
                    Se(NoNonZeros) = S(k,m)-k0^2*T(k,m);
                    
                    if(regions(ie)~=2&&regions(ie)~=3&&regions(ie)~=4)%check if located in a PML are if so there is no stimulation of incriment wave
                        R(index(n(k)),1:Nsources) = R(index(n(k)),1:Nsources)+(-S_b(k,m)+k0^2*T_b(k,m)).*E0.*exp(-1i*k0.*(p(1,n(k)).*coord(1,1:Nsources).*sind(theta)+p(2,n(k)).*coord(2,1:Nsources).*cosd(theta)));%p(1,n(m)) represent xcor of node that is on check
                    end
                else
                    %if node q is a known update result matrix
                    if(regions(ie)~=2&&regions(ie)~=3&&regions(ie)~=4)%check if located in a PML are if so there is stimulation of incriment wave
                        R(index(n(k)),1) = R(index(n(k)),1)-(S(k,m)-k0^2*T(k,m))*Ezs0(n(m),1:Nsources)+(-S_b(k,m)+k0^2*T_b(k,m)).*E0.*exp(-1i*k0.*(p(1,n(k)).*coord(1,1:Nsources).*sind(theta)+p(2,n(k)).*coord(2,1:Nsources).*cosd(theta)));
                    else
                        R(index(n(k)),1) = R(index(n(k)),1)-(S(k,m)-k0^2*T(k,m))*Ezs0(n(m),1:Nsources);
                    end
                end
            end
        end
    end
end

%Auxialiary matrix trimming cause of no apriopri knownledge of number of
%nodes with index_id=0
I = I(1:NoNonZeros,1);
J = J(1:NoNonZeros,1);
Se = Se(1:NoNonZeros,1);

%stiffness matrix formulation with sparse for better speed
Q = sparse(I,J,Se);

%% SYSTEM SOLVING FOR FIELD COMPUTATION IN NODES
%-------------Direct Solving------------------%
%-------Scattered field computation-----------%
Ezs=Q\R;%scattered field computation
for in=1:Nnodes
    if(index(in)~=0)
        Ezs0(in,1:Nsources) = Ezs(index(in),1:Nsources);
    end
end
%---------Total field computation-------------%
Ez=zeros(Nnodes,Nsources);%total field matrix

%Ez=Ezs0+E0.*exp(-1i*k0.*(p(1,:).*coord(1,1:Nsources).*sind(theta)+p(2,:).*coord(1,1:Nsources).*cosd(theta))');
%eror! cannot vectorize completely. for is needed

for in=1:Nnodes
    if(p(1,in)>=gd(4,1)&&p(1,in)<=gd(3,2)&&p(2,in)>=gd(9,3)&&p(2,in)<=gd(7,4))
        Ez(in,1:Nsources) = Ezs0(in,1:Nsources)+E0.*exp(-1i*k0.*(p(1,in).*coord(1,1:Nsources).*sind(theta)+p(2,in).*coord(2,1:Nsources).*cosd(theta)));
    else
        Ez(in) = Ezs0(in);
    end
end
%% FINAL RESULTS FIELD PLOTS
% %---------------Result plot-------------------%
% figure %scattered field only
% pdeplot(p/l,e,t,'xydata',abs(Ezs0),'mesh','on','contour','on','title','Scattered Field')
% axis([d/l (length-d)/l d/l (length-d)/l]);
% % % axis tight;axis equal;
% xlabel('x');ylabel('y','Rotation',0);
% figure %total field
% pdeplot(p/l,e,t,'xydata',abs(Ez),'mesh','on','contour','on','title','Total Field')
% axis([d/l (length-d)/l d/l (length-d)/l]);
% % axis tight;axis equal;
% xlabel('x');ylabel('y','Rotation',0);
toc