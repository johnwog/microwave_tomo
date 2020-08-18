clear all
tic
format long g
%% %--------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%--------------------------FORWARD PROBLEM------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%% PARAMETERS DEFIMENENT
%--------------------Parameters definition---------------------%
e0 = 8.854e-12; %vacuum permitivity,electric constant
m0 = 4e-7*pi; %vacuum permeability,magnetic constant
c0 = 3e8; %speed of light vacuum
f = 1e9; %frequency
w = 2*pi*f; %angular frequency
l = c0/f; %wavelength
k0 = w*sqrt(e0*m0); %vacuum wavenumber
E0 = 1; %incriment wave amplitude
%% ELECTRICAL PROPERTIES DEFINEMENT

%-------------Different area dielectric properties-------------%
%air
er0 = eye(3); %vacuum relative electric constant
mr0 = eye(3); %vacuum relative magnetic constant
mr0_inv = inv(mr0);
er_air = 1;
mr_air = 1;
kb_air = k0;

%PML parmeters
a = 1; %PML a parameter real part
b = 9; %PML a parameter imaginary part 
         %for thickness d=0.15l reflectance factor 3.1*e-5 if b=5

%vertical PML layer
erx = [1/(a-b*1i) 0 0;0 a-b*1i 0;0 0 a-b*1i];%electic constant tensor
mrx = erx;%megnetic constant tensor
mrx_inv = inv(mrx);
kb_x = w*sqrt(e0*m0*erx(3,3)); %z-axis wavenumber

%horizontal PML layer
ery = [a-b*1i 0 0;0 1/(a-b*1i) 0;0 0 a-b*1i];%electic constant tensor
mry = ery;%megnetic constant tensor
mry_inv = inv(mry);
kb_y = w*sqrt(e0*m0*ery(3,3)); %z-axis wavenumber

%PML intersection
er_inter = erx*ery;
mr_inter = mrx*mry;
mr_inter_inv = inv(mr_inter);
kb_inter = w*sqrt(e0*m0*er_inter(3,3)); %z-axis wavenumber

%scatterer1
er1 = 1.6-0.2i; %complex relative electric constant
ers1 = diag([er1 er1 er1]);%electic constant tensor
mrs1 = eye(3);%megnetic constant tensor
mrs1_inv = inv(mrs1);

%scatterer2
er2 = 1.3-0.4i; %complex relative electric constant
ers2 = diag([er2 er2 er2]);%electic constant tensor
mrs2 = eye(3);%megnetic constant tensor
mrs2_inv = inv(mrs2);

%backround medium
kb = k0; %backround medium wavenumber
er_b = 1; %bakcround relative dielectirc constant
mr_b = 1; %bakcround relative permeability constant

%% PML AREA AND SCATTERER GEOMETRY
%-------------------GEOMETRY DEFINITION-----------------------%
%PML parmeters
d = 0.15*l; %PML thickness
dis = 5*0.15*l; %distance from PML centered
length = 2*d+2*dis+3*l; %PML length by changing dis and d we control the size of the domain
param = l; %parameter umlaut profile size is calculated

%geometry matrix scatterer centered in the start of the coordinate system
gd=[3 4 -dis-d-1.5*l      -dis-1.5*l        -dis-1.5*l        -dis-d-1.5*l    -dis-1.5*l-d   -dis-1.5*l-d     dis+1.5*l+d    dis+1.5*l+d; %left vertical PML (R1)
    3 4  dis+1.5*l         dis+1.5*l+d       dis+d+1.5*l       dis+1.5*l      -dis-1.5*l-d   -dis-1.5*l-d     dis+1.5*l+d    dis+1.5*l+d; %right vertical PML (R2)
    3 4 -dis-d-1.5*l       dis+d+1.5*l       dis+d+1.5*l      -dis-d-1.5*l    -dis-1.5*l-d   -dis-1.5*l-d    -dis-1.5*l     -dis-1.5*l;   %lower horizontal PML (R3) 
    3 4 -dis-d-1.5*l       dis+d+1.5*l       dis+d+1.5*l      -dis-d-1.5*l     dis+1.5*l      dis+1.5*l       dis+d+1.5*l    dis+1.5*l+d; %upper horizontal PML (R4)
    3 4 -1.5*l-dis         1.5*l+dis         1.5*l+dis        -1.5*l-dis      -1.5*l-dis     -1.5*l-dis       1.5*l+dis      1.5*l+dis;  %between medium (R5)
    3 4 -0.5*l             0.5*l             0.5*l            -0.5*l          -0.5*l         -0.5*l           0.5*l          0.5*l;
    3 4 -l                 l                 l                -l              -l             -l               l              l;
    ]'; 

% %are centroids calculation
% %cen matrix containts the centroids of the different areas in the same row
% %as defined in gd matrix.The first row contains the x coord, the second the
% %y coord.%%%Not needed eventually
% cen=[0.5*(gd(3,1)+gd(4,1)) 0.5*(gd(3,2)+gd(4,2)) 0.5*(gd(3,3)+gd(4,3)) 0.5*(gd(3,4)+gd(4,4)) 0.5*(gd(3,5)+gd(4,5)) 0.5*(gd(3,6)+gd(4,6));
%     0.5*(gd(7,1)+gd(10,1)) 0.5*(gd(7,2)+gd(10,2)) 0.5*(gd(7,3)+gd(10,3)) 0.5*(gd(7,4)+gd(10,4)) 0.5*(gd(7,5)+gd(10,5)) 0.5*(gd(7,6)+gd(10,6))];
%     if(cen(1,5)==cen(1,6)) %if scatterer and between medium have the same centroid eliminate scatterer centroid
%         cen(1,6)=0;
%     end
%     if(cen(2,5)==cen(2,6))
%         cen(2,6)=0;
%     end
%-------------Decompesed Geometry-------------%
%decompesed geometry matrix
d1 = decsg(gd,'R1+R2+R3+R4+R5+R6+R7',[abs('R1')' abs('R2')' abs('R3')' abs('R4')' abs('R5')' abs('R6')' abs('R7')']);  

%% ORIGINAL MESH DEFINEMENT
%mesh matrices 
[p,e,t] = initmesh(d1,'hmax',l/20); 
%hmax property sets no element side greater than l/12 for a satisfactory granularity
%pdeplot(p,e,t);axis tight;axis equal;

% %% %------------Check regions fromt t matrix---------------%
% 
% %-------------Numbers of nodes,elements,edges--------------% 
% Nnodes = size(p,2); %number of nodes
% Nelements = size(t,2); %number of elements
% Nedges = size(e,2); %number of edges
%     
% area_t = unique(t(4,:)); %t matrix areas
% %col = jet(size(unique(t(4,:)),2)); %diffent color creation for every region defined by t matrix
% col = rand(size(area_t,2),3);
% %subplot(2,2,[1 3])
% for ie=1:Nelements
%     xi = p(1,t(1:3,ie));%element ie centroid x coor
%     yi = p(2,t(1:3,ie));%element ie centroid y coor
%     for iin=1:size(area_t,2)
%         if(t(4,ie) == area_t(iin))%draw if element in first area of t matrix
%             line([xi(1) xi(2)],[yi(1) yi(2)],'Color',col(iin,:));
%             line([xi(2) xi(3)],[yi(2) yi(3)],'Color',col(iin,:));
%             line([xi(3) xi(1)],[yi(3) yi(1)],'Color',col(iin,:));
%         end
%     end
% end
% h_a = zeros(5,size(area_t,2)); %auxilary matrix creation for saving 5 triangles from every region
% x_cen = zeros(5,size(area_t,2)); %x centroid coord
% y_cen = zeros(5,size(area_t,2)); %y centroid coord
% for iin=1:size(area_t,2)
%     op = find(t(4,:)==iin,10); %find 10 nodes of area with number iin
%     h_a(:,iin) = op(randi([1,10],1,5)); %select 5 random nodes of the first 10 nodes which are in ascending order
%     for jjn=1:5 %find the centroid of the 5 selected nodes
%         x_cen(jjn,iin) = sum(p(1,t(1:3,h_a(jjn,iin))))/3;
%         y_cen(jjn,iin) = sum(p(2,t(1:3,h_a(jjn,iin))))/3;
%     end
%     text(x_cen(:,iin),y_cen(:,iin),sprintf('%d',iin),'Color','k','FontSize',10)%plot the 5 number nodes with their numeration in their coord
% end
% axis tight;axis equal; %create subplot with the number corresponding to each color
% %     subplot(2,2,4)
% %     figure
% %     xiii = linspace(0,1,size(area_t,2));
% %     yiii = 0.5*ones(1,size(area_t,2));
% %     for iin=1:size(area_t,2)
% %         text(xiii(iin),yiii(iin), sprintf('%d',area_t(iin)),'Color',col(iin,:));
% %     end
%% MESH REFINEMNET FOR A DENSER MESH INSIDE THE SCATTERER AND PML-REGION MATRIX DEFINEMENT

%-----------------Mesh refinement and region definitions--------------%
for rep=1:2 %number of refinements
    %--------------Numbers of nodes,elements,edges-------------% 
    Nnodes = size(p,2); %number of nodes
    Nelements = size(t,2); %number of elements
    Nedges = size(e,2); %number of edges

    it=[find(t(4,:)==1) find(t(4,:)==3) find(t(4,:)==4) find(t(4,:)==5) find(t(4,:)==6) find(t(4,:)==7) find(t(4,:)==8) find(t(4,:)==9) find(t(4,:)==10) find(t(4,:)==11)];
        
    [p e t] = refinemesh(d1,p,e,t,it');

    %--Numbers of nodes,elements,edges,centroids--% 
    Nnodes = size(p,2); %number of nodes
    Nelements = size(t,2); %number of elements
    Nedges = size(e,2); %number of edges
end

%-----------------Region of each triangle definemet------------------%
%characteristic number for each area 1=between medium 2=horizontal PML
%3=vertical PML 4=PML intersection 5=scatterer
regions = ones(1,Nelements);

%----Between Medium-----%
regions(t(4,:)==2) = 1;

%----Horizontal PML-----%
regions(t(4,:)==4) = 2;
regions(t(4,:)==10) = 2;

%-----Vertical PML------%
regions(t(4,:)==8) = 3;
regions(t(4,:)==9) = 3;

%----PML Intersection---%
regions(t(4,:)==3) = 4;
regions(t(4,:)==5) = 4;
regions(t(4,:)==6) = 4;
regions(t(4,:)==7) = 4;

%-------Scatterer1-------%
regions(t(4,:)==1) = 5;

%-------Scatterer2-------%
regions(t(4,:)==11) = 6;

%% ELEMENT CENTROIDS CALCULATIONS
cen = zeros(2,Nelements);
cen(1,:)=(p(1,t(1,:))+p(1,t(2,:))+p(1,t(3,:)))/3;
cen(2,:)=(p(2,t(1,:))+p(2,t(2,:))+p(2,t(3,:)))/3;

% %% PLOT TEST FOR AREA CHARACTIRISATION
% %---------------------Test plot 1-------------------------%
% figure
% hold on
% col2 = [0 0 1;  %blue color
%         1 1 0;  %yellow color
%         1 0 1;  %magenta color
%         0 1 0;  %green color
%         1 0 0;];%red color
% for ie=1:Nelements
%     xii = p(1,t(1:3,ie));%element ie centroid x coor
%     yii = p(2,t(1:3,ie));%element ie centroid y coor
%     for iin=1:5
%         if(regions(ie)==iin)%draw if element in first area of t matrix by creating triangles
% %             line([xii(1) xii(2)],[yii(1) yii(2)],'Color',col2(iin,:));
% %             line([xii(2) xii(3)],[yii(2) yii(3)],'Color',col2(iin,:));
% %             line([xii(3) xii(1)],[yii(3) yii(1)],'Color',col2(iin,:));
%               %creating filled triangles
%               fill([xii(1) xii(2) xii(3)],[yii(1) yii(2) yii(3)],col2(iin,:),'EdgeColor',col2(iin,:)); 
%         end
%     end
% end
% line([-0.5*l 0.5*l],[-0.5*l -0.5*l],'Color','k')
% line([-0.5*l 0.5*l],[0.5*l 0.5*l],'Color','k')
% line([-0.5*l -0.5*l],[-0.5*l 0.5*l],'Color','k')
% line([0.5*l 0.5*l],[-0.5*l 0.5*l],'Color','k')
% axis tight;axis equal;
% hold off
%% DIRICHLET BOUNDARY CONDITIONS
%--------------DIRICHLET NODES LOCALISATION----------------%
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
% %-------------------Test plot 2----------------------%
% pdeplot(p,e,t,'xydata',real(Ezs0),'zdata',real(Ezs0),'mesh','on')

%% CREATE A NEW NUMERATION FOR NODES AS UNKNOWNS
%------------------UNKNOWN NUMERATION------------------%
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

%-----------Multiple sources and place definement----------%
Nsources = 30; %Number of sources
Ie = 1e-3; %current of infinite line source

radius_src = 4*l;%0.5*l; %radius of circle in which receivers are uniformly distributed
th_src = 175-(0:(Nsources-1))*(360/Nsources); %angles of receivers location fri
x_src = radius_src*cosd(th_src); %x coord of sources
y_src = radius_src*sind(th_src); %y coord of sources

%% STIFFNESS AND MASS MATRIX FORMULATION FOR MULTIPLE INCIDENT WAVES-EQUATION SYSTEM CONSTRUCTION

%--------------------Matrices initialization----------------%
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

%--------------Stiffness Matrix formulation-----------------%
for ie=1:Nelements
    n(1:3) = t(1:3,ie); %nodes 1-3 of element ie
    region = t(4,ie); %region of element ie
    x(1:3) = p(1,n(1:3)); %x-cord of nodes 1-3 of element ie
    y(1:3) = p(2,n(1:3)); %y-cord of nodes 1-3 of element ie
    
    %------------area of element ie computation-------------%
    D = det([ones(3,1) x' y']);
    Ae=abs(D)/2;
    %-----------------cartesian-simplex---------------------%
    a = [(x(2)*y(3)-x(3)*y(2))/D, (x(3)*y(1)-x(1)*y(3))/D, (x(1)*y(2)-x(2)*y(1))/D];
    b = [(y(2)-y(3))/D, (y(3)-y(1))/D, (y(1)-y(2))/D];
    c = [(x(3)-x(2))/D, (x(1)-x(3))/D, (x(2)-x(1))/D];
    %-----------dielectric properties definement------------%
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
    elseif(regions(ie)==6)
        er_z = ers1(3,3);
        mr_x = mrs1_inv(1,1);
        mr_y = mrs1_inv(2,2);
    elseif(regions(ie)==5)
        er_z = ers2(3,3);
        mr_x = mrs2_inv(1,1);
        mr_y = mrs2_inv(2,2);
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
                    NoNonZeros = NoNonZeros+1; %increase counter by one for trimming
                    
                    I(NoNonZeros) = index(n(k));
                    J(NoNonZeros) = index(n(m));
                    Se(NoNonZeros) = S(k,m)-k0^2*T(k,m);
                    
                    if(regions(ie)~=2&&regions(ie)~=3&&regions(ie)~=4)%check if located in a PML are if so there is no stimulation of incriment wave
                        pos_vec = sqrt((p(1,n(k))-x_src).^2+(p(2,n(k))-y_src).^2);
                        R(index(n(k)),1:Nsources) = R(index(n(k)),1:Nsources)+(-S_b(k,m)+k0^2*T_b(k,m))*(-0.25*Ie*w*m0).*besselh(0,2,k0*pos_vec);%p(1,n(m)) represent xcor of node that is on check
                    end
                else
                    %if node q is a known update result matrix
                    if(regions(ie)~=2&&regions(ie)~=3&&regions(ie)~=4)%check if located in a PML are if so there is stimulation of incriment wave
                        pos_vec = sqrt((p(1,n(k))-x_src).^2+(p(2,n(k))-y_src).^2);
                        R(index(n(k)),1) = R(index(n(k)),1)-(S(k,m)-k0^2*T(k,m))*Ezs0(n(m),1:Nsources)+(-S_b(k,m)+k0^2*T_b(k,m))*(-0.25*Ie*w*m0).*besselh(0,2,k0*pos_vec);
                    else
                        pos_vec = sqrt((p(1,n(k))-x_src).^2+(p(2,n(k))-y_src).^2); 
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
%------------------Direct Solving------------------------%
%-------------Scattered field computation----------------%
Ezs=Q\R;%scattered field computation
for in=1:Nnodes
    if(index(in)~=0)
        Ezs0(in,1:Nsources) = Ezs(index(in),1:Nsources);
    end
end
%---------------Total field computation-------------------%
Ez=zeros(Nnodes,Nsources);%total field matrix

%Ez=Ezs0+E0.*exp(-1i*k0.*(p(1,:).*coord(1,1:Nsources).*sind(theta)+p(2,:).*coord(1,1:Nsources).*cosd(theta))');
%eror! cannot vectorize completely. for is needed

for in=1:Nnodes
    if(p(1,in)>=gd(4,1)&&p(1,in)<=gd(3,2)&&p(2,in)>=gd(9,3)&&p(2,in)<=gd(7,4)) %Incriment has meaning only in non PML areas
        pos_vec = sqrt((p(1,in)-x_src).^2+(p(2,in)-y_src).^2);
        Ez(in,1:Nsources) = Ezs0(in,1:Nsources)+(-0.25*Ie*w*m0)*besselh(0,2,k0*pos_vec);
    else
        Ez(in) = Ezs0(in);
    end
end
%% FINAL RESULTS FIELD PLOTS
% %--------------------Result plot------------------------%
% figure %scattered field only
% pdeplot(p/l,e,t,'xydata',abs(Ezs0),'mesh','on','contour','on','title','Scattered Field')
% axis([(-dis-0.5*l)/l (dis+0.5*l)/l (-dis-0.5*l)/l (dis+0.5*l)/l]);
% % % axis tight;axis equal;
% xlabel('x');ylabel('y','Rotation',0);
% figure %total field
% pdeplot(p/l,e,t,'xydata',abs(Ez),'mesh','on','contour','on','title','Total Field')
% axis([d/l (length-d)/l d/l (length-d)/l]);
% % axis tight;axis equal;
% xlabel('x');ylabel('y','Rotation',0);

%% %--------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%--------------------------INVERSE PROBLEM------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

%% PML AREA AND IMAGING DOMAIN GEOMETRY
%-----geometry matrix scatterer centered in the start of the coordinate system-------&
gd_inv=[3 4 -dis-d-1.5*l      -dis-1.5*l        -dis-1.5*l        -dis-d-1.5*l    -dis-1.5*l-d   -dis-1.5*l-d     dis+1.5*l+d    dis+1.5*l+d; %left vertical PML (R1)
        3 4  dis+1.5*l         dis+1.5*l+d       dis+d+1.5*l       dis+1.5*l      -dis-1.5*l-d   -dis-1.5*l-d     dis+1.5*l+d    dis+1.5*l+d; %right vertical PML (R2)
        3 4 -dis-d-1.5*l       dis+d+1.5*l       dis+d+1.5*l      -dis-d-1.5*l    -dis-1.5*l-d   -dis-1.5*l-d    -dis-1.5*l     -dis-1.5*l;   %lower horizontal PML (R3) 
        3 4 -dis-d-1.5*l       dis+d+1.5*l       dis+d+1.5*l      -dis-d-1.5*l     dis+1.5*l      dis+1.5*l       dis+d+1.5*l    dis+1.5*l+d; %upper horizontal PML (R4)
        3 4 -1.5*l-dis         1.5*l+dis         1.5*l+dis        -1.5*l-dis      -1.5*l-dis     -1.5*l-dis       1.5*l+dis      1.5*l+dis;  %between medium (R5)
        3 4 -1.5*l             1.5*l             1.5*l            -1.5*l          -1.5*l         -1.5*l           1.5*l          1.5*l;]';    %imaging domain square (R6)
%-------------------Decompesed Geometry------------------------%
%decompesed geometry matrix
d1_inv = decsg(gd_inv,'R1+R2+R3+R4+R5+R6',[abs('R1')' abs('R2')' abs('R3')' abs('R4')' abs('R5')' abs('R6')']);     

%% ORIGINAL MESH DEFINEMENT
%mesh matrices 
[p_inv,e_inv,t_inv] = initmesh(d1_inv,'hmax',l/17); 
%hmax property sets no element side greater than l/12 for a satisfactory granularity
%pdeplot(p,e,t);axis tight;axis equal;

%% MESH REFINEMNET FOR A DENSER MESH INSIDE THE PML - ELEMENT CENTROIDS CALCULATION

%%-----------Mesh refinement and region definitions------------%
for rep=1:2 %number of refinements
    
    %-------------Numbers of nodes,elements,edges--------------% 
    Nnodes_inv = size(p_inv,2); %number of nodes
    Nelements_inv = size(t_inv,2); %number of elements
    Nedges_inv = size(e_inv,2); %number of edges

    %------------Region of each triangle definemet-------------%
    
    %characteristic number for each area 1=backround medium 2=horizontal PML
    %3=vertical PML 4=PML intersection 5=imaging domain
    
    regions_inv = ones(1,Nelements_inv);
    
    for ie=1:Nelements_inv
        x_cen = sum(p_inv(1,t_inv(1:3,ie)))/3;
        y_cen = sum(p_inv(2,t_inv(1:3,ie)))/3;
        if(x_cen>gd_inv(3,6)&&x_cen<gd_inv(4,6)&&y_cen>gd_inv(7,6)&&y_cen<gd_inv(9,6))
            regions_inv(ie) = 5; %imaging domain
        elseif(x_cen>gd_inv(4,1)&&x_cen<gd_inv(3,6)&&y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,4))
            regions_inv(ie) = 1; %backround medium
        elseif(x_cen>gd_inv(3,6)&&x_cen<gd_inv(4,6)&&y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,6))
            regions_inv(ie) = 1; %backround medium
        elseif(x_cen>gd_inv(3,6)&&x_cen<gd_inv(4,6)&&y_cen>gd_inv(9,6)&&y_cen<gd_inv(7,4))
            regions_inv(ie) = 1; %backround medium
        elseif(x_cen>gd_inv(4,6)&&x_cen<gd_inv(3,2)&&y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,4))
            regions_inv(ie) = 1; %backround medium
        elseif(x_cen>gd_inv(4,1)&&x_cen<gd_inv(3,2)&&y_cen>gd_inv(7,3)&&y_cen<gd_inv(9,3))
            regions_inv(ie) = 2; %lower horizontal PML
        elseif(x_cen>gd_inv(4,1)&&x_cen<gd_inv(3,2)&&y_cen>gd_inv(7,4)&&y_cen<gd_inv(9,4))
            regions_inv(ie) = 2; %upper horizontal PML
        elseif(y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,4)&&x_cen>gd_inv(3,1)&&x_cen<gd_inv(4,1))
            regions_inv(ie) = 3; %left vertical PML
        elseif(y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,4)&&x_cen>gd_inv(3,2)&&x_cen<gd_inv(4,2))
            regions_inv(ie) = 3; %right vertical PML
        elseif(x_cen>gd_inv(3,1)&&x_cen<gd_inv(4,1)&&y_cen>gd_inv(7,3)&&y_cen<gd_inv(9,3))
            regions_inv(ie) = 4; %downleft PML intersection
        elseif(x_cen>gd_inv(3,2)&&x_cen<gd_inv(4,2)&&y_cen>gd_inv(7,3)&&y_cen<gd_inv(9,3))
            regions_inv(ie) = 4; %downright PML intersection
        elseif(x_cen>gd_inv(3,1)&&x_cen<gd_inv(4,1)&&y_cen>gd_inv(7,4)&&y_cen<gd_inv(9,4))
            regions_inv(ie) = 4; %upperleft PML intersection
        elseif(x_cen>gd_inv(3,2)&&x_cen<gd_inv(4,2)&&y_cen>gd_inv(7,4)&&y_cen<gd_inv(9,4))
            regions_inv(ie) = 4; %upperight PML intersection
        end
    end
    it_inv=[find(regions_inv==2) find(regions_inv==3) find(regions_inv==4)];% find(regions_inv==5)];

    [p_inv e_inv t_inv] = refinemesh(d1_inv,p_inv,e_inv,t_inv,it_inv');

    %-------------Numbers of nodes,elements,edges--------------% 
    Nnodes_inv = size(p_inv,2); %number of nodes
    Nelements_inv = size(t_inv,2); %number of elements
    Nedges_inv = size(e_inv,2); %number of edges
    
    regions_inv = ones(1,Nelements_inv); %regions matrix that has region id for every element
    cen_inv = zeros(2,Nelements_inv); %elements centroids calculation x coord in row 1 y coord in row 2
    
    for ie=1:Nelements_inv
        x_cen = sum(p_inv(1,t_inv(1:3,ie)))/3;
        y_cen = sum(p_inv(2,t_inv(1:3,ie)))/3;
        cen_inv(1,ie) = x_cen;
        cen_inv(2,ie) = y_cen;
        if(x_cen>gd_inv(3,6)&&x_cen<gd_inv(4,6)&&y_cen>gd_inv(7,6)&&y_cen<gd_inv(9,6))
            regions_inv(ie) = 5; %imaging domain
        elseif(x_cen>gd_inv(4,1)&&x_cen<gd_inv(3,6)&&y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,4))
            regions_inv(ie) = 1; %between medium
        elseif(x_cen>gd_inv(3,6)&&x_cen<gd_inv(4,6)&&y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,6))
            regions_inv(ie) = 1; %between medium
        elseif(x_cen>gd_inv(3,6)&&x_cen<gd_inv(4,6)&&y_cen>gd_inv(9,6)&&y_cen<gd_inv(7,4))
            regions_inv(ie) = 1; %between medium
        elseif(x_cen>gd_inv(4,6)&&x_cen<gd_inv(3,2)&&y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,4))
            regions_inv(ie) = 1; %between medium
        elseif(x_cen>gd_inv(4,1)&&x_cen<gd_inv(3,2)&&y_cen>gd_inv(7,3)&&y_cen<gd_inv(9,3))
            regions_inv(ie) = 2; %lower horizontal PML
        elseif(x_cen>gd_inv(4,1)&&x_cen<gd_inv(3,2)&&y_cen>gd_inv(7,4)&&y_cen<gd_inv(9,4))
            regions_inv(ie) = 2; %upper horizontal PML
        elseif(y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,4)&&x_cen>gd_inv(3,1)&&x_cen<gd_inv(4,1))
            regions_inv(ie) = 3; %left vertical PML
        elseif(y_cen>gd_inv(9,3)&&y_cen<gd_inv(7,4)&&x_cen>gd_inv(3,2)&&x_cen<gd_inv(4,2))
            regions_inv(ie) = 3; %right vertical PML
        elseif(x_cen>gd_inv(3,1)&&x_cen<gd_inv(4,1)&&y_cen>gd_inv(7,3)&&y_cen<gd_inv(9,3))
            regions_inv(ie) = 4; %downleft PML intersection
        elseif(x_cen>gd_inv(3,2)&&x_cen<gd_inv(4,2)&&y_cen>gd_inv(7,3)&&y_cen<gd_inv(9,3))
            regions_inv(ie) = 4; %downright PML intersection
        elseif(x_cen>gd_inv(3,1)&&x_cen<gd_inv(4,1)&&y_cen>gd_inv(7,4)&&y_cen<gd_inv(9,4))
            regions_inv(ie) = 4; %upperleft PML intersection
        elseif(x_cen>gd_inv(3,2)&&x_cen<gd_inv(4,2)&&y_cen>gd_inv(7,4)&&y_cen<gd_inv(9,4))
            regions_inv(ie) = 4; %upperight PML intersection
        end
    end
end
% %% PLOT TEST FOR AREA CHARACTIRISATION
% %---------------Test plot 1-------------------%
% figure
% hold on
% col2 = [1 1 0; %yellow color
%         1 0 1; %magenta color
%         1 0 0; %red color
%         0 1 0; %green color
%         0 0 1;]; %blue color
% for ie=1:Nelements_inv
%     xii = p_inv(1,t_inv(1:3,ie));%element ie centroid x coor
%     yii = p_inv(2,t_inv(1:3,ie));%element ie centroid y coor
%     for iin=1:5
%         if(regions_inv(ie)==iin)%draw if element in first area of t matrix
%             line([xii(1) xii(2)],[yii(1) yii(2)],'Color',col2(iin,:));
%             line([xii(2) xii(3)],[yii(2) yii(3)],'Color',col2(iin,:));
%             line([xii(3) xii(1)],[yii(3) yii(1)],'Color',col2(iin,:));
%         end
%     end
% end
% axis tight;axis equal;
% hold off
%% DIRICHLET BOUNDARY CONDITIONS
%----------------DIRICHLET NODES LOCALISATION-----------------%
node_id_inv = ones(Nnodes_inv,1); %node identity matrix if unknown id=1 if Dirichlet=0

% for id=1:Nedges_inv
%     left = e_inv(6,id); %left area of edge
%     right = e_inv(7,id); %right area of edge
%     
%     if(left==0||right==0)
%         n1 = e_inv(1,id);
%         n2 = e_inv(2,id);
%         node_id_inv(n1) = 0;
%         node_id_inv(n2) = 0;
%     end
% end

%% CREATE A NEW NUMERATION FOR NODES AS UNKNOWNS
%-------------------UNKNOWN NUMERATION-----------------------%
index_inv = zeros(Nnodes_inv,1);%numeration of nodes as unknowns

Nunknowns_inv = 0; %number of unknowns

for in=1:Nnodes_inv
    if(node_id_inv(in)==1) %check if node is an unknown
        Nunknowns_inv = Nunknowns_inv+1;%number of unknowns definement
        index_inv(in) = Nunknowns_inv;
    else
        index_inv(in) = 0;
    end
end

%% RECEIVERS DEFINITION PLACES AND MEASURED FIELDS NOISE ADDITION TOO
%mesh of forward problem will be used as fields are calculated in that mesh

radius = 2.17*l;%0.5*l; %radius of circle in which receivers are uniformly distributed
th_rec = 175-(0:(Nsources-1))*(360/Nsources); %angles of receivers location fri
x_rec = radius*cosd(th_rec); %x coord of sources
y_rec = radius*sind(th_rec); %y coord of sources
ita = 0.03; %noise parameter

dis_cen = zeros(Nsources,Nelements); %computing x-y recever distance from every element centroid
for is=1:Nsources
    dis_cen(is,:) = sqrt((x_rec(is)-cen(1,:)).^2+(y_rec(is)-cen(2,:)).^2);
end

[B IX] = sort(dis_cen,2); %sort each row of distances matrix returned the sorted matrix B and the place before sorting in IX

B = B(:,1:7); %since we are on 2D geometry each element has around it seven triangles tops so we are selecting the first 7 closest triangles
IX = IX(:,1:7);

res = zeros(1,Nsources); %matrix that saves the element every source belongs
z_simplex = zeros(Nsources,3); %basis functions

E_mes = zeros(Nsources,Nsources); %measured fields first row all fields for source illumination one second row for source illumination two...

for is=1:Nsources
    for ie=1:7
        n(1:3) = t(1:3,IX(is,ie)); %nodes 1-3 of element ie
        x(1:3) = p(1,n(1:3)); %x-cord of nodes 1-3 of element ie
        y(1:3) = p(2,n(1:3)); %y-cord of nodes 1-3 of element ie

        %---------area of element ie computation---------%
        D = det([ones(3,1) x' y']);
        Ae=abs(D)/2;
        %-------------cartesian-simplex------------------%
        a = [(x(2)*y(3)-x(3)*y(2))/D, (x(3)*y(1)-x(1)*y(3))/D, (x(1)*y(2)-x(2)*y(1))/D];
        b = [(y(2)-y(3))/D, (y(3)-y(1))/D, (y(1)-y(2))/D];
        c = [(x(3)-x(2))/D, (x(1)-x(3))/D, (x(2)-x(1))/D];
        %-------------simplex computation----------------%
        z = a+b*x_rec(is)+c*y_rec(is);
        
        %for a point to belong to element ie simplex must be between [0 1]
        %0<=simplex<=1
        if(z>=0&z<=1)%(z(1)>=0&&z(1)<=1&&z(2)>=0&&z(2)<=1&&z(3)>=0&&z(3)<=1)
            tmp = IX(is,ie); %element in which receiver is belongs
            tmp1 = z; %simplex of belonging element
            nodes = n(1:3); %nodes of belonging element
            break
        end
    end
    res(is) = tmp; %element in which receiver is belongs
    z_simplex(is,:) = tmp1; %simplex of belonging element
    
    %E_mes(is,:) because is index give us in what triangle source belongs,%so for that triangle we need to calculate the fields for all different illuminations :Z
    %Tacking the transpose of above computed matrix in order to be coherent with the matrix dimension in the necessary matrix multiplications now the measured fields 
                %for source 1 are located in column 1 and so on
    E_mes(is,:) = z*Ezs0(nodes,:);%+max(abs(Ezs0(:,iss)))*ita*(sqrt(2)^(-1))*(unifrnd(-1,1,1,Nsources)+unifrnd(-1,1,1,Nsources)*1i); %computing measured fields with additive noise uniformly distributed between -1 1
end

%% MS MD OPERATORS CONSTRUCTIONS 
%inverse problem mesh will be used

%---------------------MS matrx-------------------------%
% MS = spalloc(Nsources,Nnodes_inv,3*Nsources);
ms_i = zeros(1,3*Nsources);
ms_j = zeros(1,3*Nsources);
simplex_j = zeros(1,3*Nsources);

dis_inv = zeros(Nsources,Nelements_inv); %computing x-y recever distance from every element centroid
for is=1:Nsources
    dis_inv(is,:) = sqrt((x_rec(is)-cen_inv(1,:)).^2+(y_rec(is)-cen_inv(2,:)).^2);
end

[B_inv IX_inv] = sort(dis_inv,2); %sort each row of distances matrix returned the sorted matrix B and the place before sorting in IX

B_inv = B_inv(:,1:7); %since we are on 2D geometry each element has around it seven triangles tops so we are selecting the first 7 closest triangles
IX_inv = IX_inv(:,1:7); %in each row 7 triangles with smallest distance from the source row numeber represents the source number

res_inv = zeros(1,Nsources); %matrix that saves the element every source belongs
z_simplex_inv = zeros(Nsources,3); %basis functions

for is=1:Nsources
    for ie=1:7
        n(1:3) = t_inv(1:3,IX_inv(is,ie)); %nodes 1-3 of element ie
        x(1:3) = p_inv(1,n(1:3)); %x-cord of nodes 1-3 of element ie
        y(1:3) = p_inv(2,n(1:3)); %y-cord of nodes 1-3 of element ie

        %---------area of element ie computation---------%
        D = det([ones(3,1) x' y']);
        Ae=abs(D)/2;
        %-------------cartesian-simplex------------------%
        a = [(x(2)*y(3)-x(3)*y(2))/D, (x(3)*y(1)-x(1)*y(3))/D, (x(1)*y(2)-x(2)*y(1))/D];
        b = [(y(2)-y(3))/D, (y(3)-y(1))/D, (y(1)-y(2))/D];
        c = [(x(3)-x(2))/D, (x(1)-x(3))/D, (x(2)-x(1))/D];
        %-------------simplex computation----------------%
        z_inv = a+b*x_rec(is)+c*y_rec(is);
        
        %for a point to belong to element ie simplex must be between [0 1]
        %0<=simplex<=1
        if(z_inv<=1&z_inv>=0)%(z_inv(1)>=0&&z_inv(1)<=1&&z_inv(2)>=0&&z_inv(2)<=1&&z_inv(3)>=0&&z_inv(3)<=1)
            tmp_ms = IX_inv(is,ie); %element in which receiver is belongs
            tmp1_ms = z_inv; %simplex of belonging element
            nodes_inv = n(1:3); %nodes of belonging element
            break
        end
    end
    mat_empty =  find(tmp1_ms==0);
    if(isempty(mat_empty)==0) %if a sources is on a nodes place the correct simplex equal 1 and 0 the others
        tmp1_ms = zeros(1,3);
        tmp1_ms(z_inv==(max(z_inv))) = 1; 
    end
    res_inv(is) = tmp_ms; %element in which receiver is belongs
    z_simplex_inv(is,:) = tmp1_ms; %simplex of belonging element

    ms_i((1+(is-1)*3):(1+(is-1)*3+2)) = is*ones(1,3);
    ms_j((1+(is-1)*3):(1+(is-1)*3+2)) = nodes_inv;
    simplex_j((1+(is-1)*3):(1+(is-1)*3+2)) = tmp1_ms;
%     MS(is,nodes_inv) = z_inv; %assiging simplex value of element nodes that receiver is belongs
end

MS = sparse(ms_i,ms_j,simplex_j,Nsources,Nunknowns_inv);

%---------------------MD matrx------------------------%
element_belong = find(regions_inv==5); %elements belonging in imaging domain D
element_belong1 = t_inv(1:3,element_belong); %select nodes of each triangle in domain D multiple same nodes will be selected
jj = sort(unique(element_belong1)); %select unique nodes in imaging domain D and sort them in ascending order. Represents the colomn id of node in matrix MD
index_id_d = jj; %index matrix for correspondence of imaging domain D node to global numeration
jj = index_inv(jj)'; %This matrix acts as an index matrix new numeration for nodes inside the domain. Fi 1 correspond to 59 but 59 is the numeration of the node in the index_inv matrix
ii = 1:size(jj,2); %Represents the rows of matrix

Nnodes_d = size(jj,2); %number of nodes in inversion domain d

% MD = spalloc(Nnodes_d,Nnodes_inv,Nnodes_d); %matrix MD definent where size(jj,2) equals to number of nodes in domain D
% 
% linearid = sub2ind(size(MD),ii,jj); %matrix indexing for iith jjth element of MD
% 
% MD(linearid) = 1;

MD = sparse(ii,jj,ones(1,Nnodes_d),Nnodes_d,Nunknowns_inv); %generating MD matrix us above but with sparse

index_inv_D = zeros(max(jj),1); %index matrix for domain D for final annotation in the global numeration mesh
index_inv_D(jj) = ii; %size of array equal to Nonodes 0 values to nodes outside D numeration No of nodes in D

%% STIFFNESS AND MASS MATRIX FORMULATION FOR MULTIPLE INCIDENT WAVES-EQUATION SYSTEM CONSTRUCTION

%---------------Matrices initialization---------------%

I_inv = zeros(9*Nelements_inv,1); %auxialry matrix for sparse usage the row value in Tb,S
J_inv = zeros(9*Nelements_inv,1); %auxialry matrix for sparse usage the column value in Tb,S
I_invd = zeros(9*Nelements_inv,1); %auxialry matrix for sparse usage the row value
J_invd = zeros(9*Nelements_inv,1); %auxialry matrix for sparse usage the row value
Se_inv = zeros(9*Nelements_inv,1);  %%auxialry matrix for sparse usage element value in stiffness matrix
Td_inv = zeros(9*Nelements_inv,1); %%auxialry matrix for sparse usage element value in Tb matirx
Tb_inv = zeros(9*Nelements_inv,1); %%auxialry matrix for sparse usage element value in Tb matirx

NoNonZeros_inv = 0; %number of non zero index nodes
NoNonZeros_inv_D = 0; %number of non zero index nodes in domain D

S=zeros(3,3);
T=zeros(3,3);
T_b=zeros(3,3);

%------------Stiffness Matrix formulation------------%
for ie=1:Nelements_inv
    n(1:3) = t_inv(1:3,ie); %nodes 1-3 of element ie
    region = t_inv(4,ie);   %region of element ie
    x(1:3) = p_inv(1,n(1:3)); %x-cord of nodes 1-3 of element ie
    y(1:3) = p_inv(2,n(1:3)); %y-cord of nodes 1-3 of element ie
    
    %---------area of element ie computation---------%
    D = det([ones(3,1) x' y']);
    Ae=abs(D)/2;
    %-------------cartesian-simplex------------------%
    a = [(x(2)*y(3)-x(3)*y(2))/D, (x(3)*y(1)-x(1)*y(3))/D, (x(1)*y(2)-x(2)*y(1))/D];
    b = [(y(2)-y(3))/D, (y(3)-y(1))/D, (y(1)-y(2))/D];
    c = [(x(3)-x(2))/D, (x(1)-x(3))/D, (x(2)-x(1))/D];
    %--------dielectric properties definement--------%
    if(regions_inv(ie)==1) %backround medium
        er_z = er_air;
        mr_x = mr_air;
        mr_y = mr_air;
        kb = kb_air;
    elseif(regions_inv(ie)==2) %horizontal PML
        er_z = ery(3,3);
        mr_x = mry_inv(1,1);
        mr_y = mry_inv(2,2);
        kb = kb_y;
    elseif(regions_inv(ie)==3) %vertical PML
        er_z = erx(3,3);
        mr_x = mrx_inv(1,1);
        mr_y = mrx_inv(2,2);
        kb = kb_x;
    elseif(regions_inv(ie)==4) %PML intersection
        er_z = er_inter(3,3);
        mr_x = mr_inter_inv(1,1);
        mr_y = mr_inter_inv(2,2);
        kb = kb_inter;
    elseif(regions_inv(ie)==5) %imaging domain
        er_z = er_air;
        mr_x = mr_air;
        mr_y = mr_air;
        kb = kb_air;
    end
        
    for k=1:3
        for m=1:3%Compute S T S' T' parameters for stiffness matrix calculation
            S(k,m) = (mr_x*c(k)*c(m)+mr_y*b(k)*b(m))*Ae;
            
            if(k==m)
                T(k,m) = Ae/6; 
            else
                T(k,m) = Ae/12;
            end
            if(node_id_inv(n(k))==1)%check if p node is an unknown
                
                 NoNonZeros_inv = NoNonZeros_inv+1; %increase counter by one for trimming 
                    
                 I_inv(NoNonZeros_inv) = index_inv(n(k));
                 J_inv(NoNonZeros_inv) = index_inv(n(m));
                 Se_inv(NoNonZeros_inv) = S(k,m)-(kb^2)*T(k,m); 
                 Tb_inv(NoNonZeros_inv) = kb^2*T(k,m);
                 if(regions_inv(ie)~=5)
                     Tb_inv(NoNonZeros_inv) = 0;
                 end
                 
                 if(regions_inv(ie)==5)
                     NoNonZeros_inv_D = NoNonZeros_inv_D+1;
                     I_invd(NoNonZeros_inv_D) = index_inv_D(index_inv(n(k)));%find(jj==index_inv(n(k)));
                     J_invd(NoNonZeros_inv_D) = index_inv_D(index_inv(n(m)));%find(jj==index_inv(n(m)));
                     Td_inv(NoNonZeros_inv_D) = T(k,m);
                 end
            end
        end
    end
end

%Auxialiary matrix trimming cause of no apriopri knownledge of number of
%nodes with index_id=0
I_inv = I_inv(1:NoNonZeros_inv,1);
J_inv = J_inv(1:NoNonZeros_inv,1);
I_invd = I_invd(1:NoNonZeros_inv_D,1);
J_invd = J_invd(1:NoNonZeros_inv_D,1);
Se_inv = Se_inv(1:NoNonZeros_inv,1);
Tb_inv = Tb_inv(1:NoNonZeros_inv,1);
Td_inv = Td_inv(1:NoNonZeros_inv_D,1);

%stiffness matrix formulation with sparse for better speed
S_inv = sparse(I_inv,J_inv,Se_inv);
Tb = sparse(I_inv,J_inv,Tb_inv);
Td = sparse(I_invd,J_invd,Td_inv);

%% L OPERATOR CONSTRUCTIONS
L = S_inv\(Tb*(MD')); %equal to inv(S_inv)*Tb*MD' but extremely faster

%% INCRIMENT FIELDS
%-------------------Incriment field----------------------%
Ez_inc=zeros(Nnodes_inv,Nsources);%incriment field matrix

%Ez=Ezs0+E0.*exp(-1i*k0.*(p(1,:).*coord(1,1:Nsources).*sind(theta)+p(2,:).*coord(1,1:Nsources).*cosd(theta))');
%eror! cannot vectorize completely. for is needed

for in=1:Nnodes_inv
    if(p_inv(1,in)>=gd_inv(4,1)&&p_inv(1,in)<=gd_inv(3,2)&&p_inv(2,in)>=gd_inv(9,3)&&p_inv(2,in)<=gd_inv(7,4)) %Incriment has meaning only in non PML areas
        pos_vec = sqrt((p_inv(1,in)-x_src).^2+(p_inv(2,in)-y_src).^2);
        Ez_inc(in,1:Nsources) = (-0.25*Ie*w*m0)*besselh(0,2,k0*pos_vec); 
    else
        Ez_inc(in) = 0;
    end
end
%% MR TERM PREPROCESSING
%--------inversion domain area average area of mesh triangles--------%
inv_area = ((gd_inv(4,6)-gd_inv(3,6)))*((gd_inv(9,6)-gd_inv(8,6))); %area of imaging domain

inv_tri_area = zeros(size(find(regions_inv==5),2),1); %matrix containing the area of each element 
simplex_tri_b = zeros(size(find(regions_inv==5),2),3); %matrix containing the bi simplex for each element
simplex_tri_c = zeros(size(find(regions_inv==5),2),3); %matrix containing the ci simplex for each element
simplex_tri_n = zeros(size(find(regions_inv==5),2),3); %matrix containing nodes for each element

count=0;
for ie=1:Nelements_inv
    if(regions_inv(ie)==5)
        count = count+1;
        n(1:3) = t_inv(1:3,ie); %nodes 1-3 of element ie
        region = t_inv(4,ie);   %region of element ie
        x(1:3) = p_inv(1,n(1:3)); %x-cord of nodes 1-3 of element ie
        y(1:3) = p_inv(2,n(1:3)); %y-cord of nodes 1-3 of element ie
        D = det([ones(3,1) x' y']);
        inv_tri_area(count) = abs(D)/2; %triangle area
        
        simplex_tri_n(count,:) = index_inv_D(index_inv(n)); %nodes of element with D-domain numeration
        
        b = [(y(2)-y(3))/D, (y(3)-y(1))/D, (y(1)-y(2))/D];
        simplex_tri_b(count,:) = b;
        
        c = [(x(3)-x(2))/D, (x(1)-x(3))/D, (x(2)-x(1))/D];
        simplex_tri_c(count,:) = c;
    end
end

inv_area_m = mean(inv_tri_area); %average triangle area of imaging domain

%--------------------------------polygons calculation--------------------------------%
mr_cell = cell(Nnodes_d,9);
int_el_cell = cell(Nnodes_d,1);
point_cell = cell(Nnodes_d,1);

w_bar2=waitbar(0,'Please Wait','Name','MR Preprocessing');

for in=1:Nnodes_d
    %--------------------finding in wich elements node in belongs--------------------%
    [~, help_bel] = ind2sub(size(t_inv(1:3,:)),find(t_inv(1:3,:)==index_id_d(in))); %finding the elements in the numeration of D domain in which node belongs
    c_clock = convhull(cen_inv(1,help_bel),cen_inv(2,help_bel)); %setting the elements in counter clockwise order
    c_clock = c_clock(1:(size(c_clock,1)-1)); %deleting the last entry as it is the same with the first
    help_bel = help_bel(c_clock); %ordering the index of the elements in counter clockwise order
    
    cen_vec = [help_bel;help_bel(1)]; %creating the vector matrix of the centroids
    
    %-------finding intersection points of triangle edges and weightcenter lines------%
    intersec_points_x = zeros(1,size(help_bel,1));
    intersec_points_y = zeros(1,size(help_bel,1));
    
    for inn=1:size(help_bel,1)
      % [xii yii]=polyxpoly(p_inv(1,t_inv(1:3,cen_vec(inn))),p_inv(2,t_inv(1:3,cen_vec(inn))),[cen_inv(1,cen_vec(inn)),cen_inv(1,cen_vec(inn+1))],[cen_inv(2,cen_vec(inn)),cen_inv(2,cen_vec(inn+1))]); %calculating intersections points
      % if(isempty(xii)||isempty(yii)) %if polyxpoly cannot work find the line equations and compute intersection
           node_tri1 = t_inv(1:3,cen_vec(inn)); %taking the nodes of triangle cen_vec
           node_tri2 = t_inv(1:3,cen_vec(inn+1)); %taking the nodes of triangle counterclockwise next to cen_vec 
           node_2 = node_tri1(ismember(node_tri1,node_tri2)); %taking the nodes of triangle cen_ven that are common for triangle cen_vec+1
           node_2 = node_2(abs(node_2-index_id_d(in))>eps); %find the other common node between the two triangles
           
           %finding the second common node between the triangle in order to define the line equation
           if(abs(p_inv(1,index_id_d(in))-p_inv(1,node_2))>eps&&abs(p_inv(2,index_id_d(in))-p_inv(2,node_2))>eps&&abs(cen_inv(1,cen_vec(inn+1))-cen_inv(1,cen_vec(inn)))>eps&&abs(cen_inv(2,cen_vec(inn+1))-cen_inv(2,cen_vec(inn)))>eps)%checking if edge is not horizontal
               slope1 = (p_inv(2,index_id_d(in))-p_inv(2,node_2))/(p_inv(1,index_id_d(in))-p_inv(1,node_2)); %calculating the slope of the first line
               slope2 = (cen_inv(2,cen_vec(inn))-cen_inv(2,cen_vec(inn+1)))/(cen_inv(1,cen_vec(inn))-cen_inv(1,cen_vec(inn+1))); %calculating the slope of the second line
               b_line1 = p_inv(2,index_id_d(in))-slope1*p_inv(1,index_id_d(in)); %calculating the y-intersept of line1
               b_line2 = cen_inv(2,cen_vec(inn))-slope2*cen_inv(1,cen_vec(inn)); %calculating the y-intersept of line2
               r_intersec = [slope1 -1;slope2 -1]\[-b_line1;-b_line2]; %calculating the intersection point
               xii = r_intersec(1); 
               yii = r_intersec(2);
           elseif(abs(p_inv(1,index_id_d(in))-p_inv(1,node_2))<=eps) %check if edge is vertical
               slope2 = (cen_inv(2,cen_vec(inn))-cen_inv(2,cen_vec(inn+1)))/(cen_inv(1,cen_vec(inn))-cen_inv(1,cen_vec(inn+1))); %calculating the slope of the second line
               b_line2 = cen_inv(2,cen_vec(inn))-slope2*cen_inv(1,cen_vec(inn)); %calculating the y-intersept of line2
               xii = p_inv(1,node_2);
               yii = slope2*xii + b_line2;
           elseif(abs(p_inv(2,index_id_d(in))-p_inv(2,node_2))<=eps) %check if edge is horizontal
               slope2 = (cen_inv(2,cen_vec(inn))-cen_inv(2,cen_vec(inn+1)))/(cen_inv(1,cen_vec(inn))-cen_inv(1,cen_vec(inn+1))); %calculating the slope of the second line
               b_line2 = cen_inv(2,cen_vec(inn))-slope2*cen_inv(1,cen_vec(inn)); %calculating the y-intersept of line2
               yii = p_inv(2,node_2);
               xii = (yii - b_line2)/slope2; 
           elseif(abs(cen_inv(1,cen_vec(inn+1))-cen_inv(1,cen_vec(inn)))<eps) %check if edge is vertical
               slope1 = (p_inv(2,index_id_d(in))-p_inv(2,node_2))/(p_inv(1,index_id_d(in))-p_inv(1,node_2)); %calculating the slope of the first line
               b_line1 = p_inv(2,index_id_d(in))-slope1*p_inv(1,index_id_d(in)); %calculating the y-intersept of line1
               xii = cen_inv(1,cen_vec(inn));
               yii = slope1*xii + b_line1;
           elseif(abs(cen_inv(2,cen_vec(inn+1))-cen_inv(2,cen_vec(inn)))<eps) %check if edge is horizontal
               slope1 = (p_inv(2,index_id_d(in))-p_inv(2,node_2))/(p_inv(1,index_id_d(in))-p_inv(1,node_2)); %calculating the slope of the first line
               b_line1 = p_inv(2,index_id_d(in))-slope1*p_inv(1,index_id_d(in)); %calculating the y-intersept of line1
               yii = cen_inv(2,cen_vec(inn));
               xii = (yii - b_line1)/slope1;
           elseif(abs(cen_inv(2,cen_vec(inn+1))-cen_inv(2,cen_vec(inn)))<eps&&abs(p_inv(1,index_id_d(in))-p_inv(1,node_2))<eps)
               xii = p_inv(1,node_2);
               yii = cen_inv(2,cen_vec(inn));
           elseif(abs(cen_inv(1,cen_vec(inn+1))-cen_inv(1,cen_vec(inn)))<eps&&abs(p_inv(2,index_id_d(in))-p_inv(2,node_2))<eps)
               xii = cen_inv(1,cen_vec(inn));
               yii = p_inv(2,node_2);
           end
     %  end
       intersec_points_x(1,inn) = xii;
       intersec_points_y(1,inn) = yii;
    end
    
    %center of the polygon
    poly_cen_x = mean([cen_inv(1,help_bel),intersec_points_x]);
    poly_cen_y = mean([cen_inv(2,help_bel),intersec_points_y]);
    if(isnan(poly_cen_x)||isnan(poly_cen_y))
        poly_cen_x = p_inv(1,index_id_d(in));
        poly_cen_y = p_inv(2,index_id_d(in));
    end
    
    %--------taking only the weightcenters and elements inside imaging domain--------%
    intersec_points = [intersec_points_x(abs(intersec_points_x)<=abs(gd_inv(3,6))&abs(intersec_points_y)<=abs(gd_inv(7,6)));intersec_points_y(abs(intersec_points_x)<=abs(gd_inv(3,6))&abs(intersec_points_y)<=abs(gd_inv(7,6)))]; %taking the points inside imaging domain

    mr_cell{in,1} = 1/polyarea(cen_inv(1,help_bel),cen_inv(2,help_bel)); %node in polygon area
    
    help_bel = help_bel(abs(cen_inv(1,help_bel))<=abs(gd_inv(3,6))&abs(cen_inv(2,help_bel))<=abs(gd_inv(7,6))); %taking triangles inside imagin domain only
    
    %-----creating point matrix consisting of weight centers and intersection points in counterclock wise order-----%
    x_point = [cen_inv(1,help_bel),intersec_points(1,:)]; %x-coord of point around node centroids and intersections
    y_point = [cen_inv(2,help_bel),intersec_points(2,:)]; %y-coord of point around node centroids and intersections
    help_x = x_point-poly_cen_x; %x-coord converting to local coordinate system with center the polygon center
    help_y = y_point-poly_cen_y; %y-coord converting to local coordinate system with center the polygon center 
    flag_mat = [help_bel', zeros(size(intersec_points(1,:)))];
    
    deg_mat = (180/pi)*atan2(help_y,help_x); %computing the angle of the vectors with edges the points of the polygon
    
    if(abs(p_inv(2,index_id_d(in))-gd_inv(9,6))<eps&&poly_cen_y<gd_inv(9,6)) %check if node is located in the upper border
        deg_mat_index1 = find(deg_mat>=90&deg_mat<=180);
        deg_mat_index2 = find(deg_mat<0);
        deg_mat_index3 = find(deg_mat<=90&deg_mat>=0);
        deg_mat1 = deg_mat(deg_mat_index1);
        deg_mat2 = deg_mat(deg_mat_index2);
        deg_mat3 = deg_mat(deg_mat_index3);
        [deg_mat1, c_clock1] = sort(deg_mat1);
        deg_mat_index1 = deg_mat_index1(c_clock1);
        [deg_mat2, c_clock2] = sort(deg_mat2);
        deg_mat_index2 = deg_mat_index2(c_clock2);
        [deg_mat3, c_clock3] = sort(deg_mat3);
        deg_mat_index3 = deg_mat_index3(c_clock3);
        deg_mat = [deg_mat3 deg_mat2 deg_mat1];
        c_clock = [deg_mat_index1 deg_mat_index2 deg_mat_index3];
    elseif(abs(p_inv(2,index_id_d(in))-gd_inv(7,6))<eps&&poly_cen_y>gd_inv(7,6)) %check if node is located in the lower border
        deg_mat_index1 = find(deg_mat>=-180&deg_mat<-90);
        deg_mat_index2 = find(deg_mat>=0);
        deg_mat_index3 = find(deg_mat>=-90&deg_mat<0);
        deg_mat1 = deg_mat(deg_mat_index1);
        deg_mat2 = deg_mat(deg_mat_index2);
        deg_mat3 = deg_mat(deg_mat_index3);
        [deg_mat1, c_clock1] = sort(deg_mat1);
        deg_mat_index1 = deg_mat_index1(c_clock1);
        [deg_mat2, c_clock2] = sort(deg_mat2);
        deg_mat_index2 = deg_mat_index2(c_clock2);
        [deg_mat3, c_clock3] = sort(deg_mat3);
        deg_mat_index3 = deg_mat_index3(c_clock3);
        deg_mat = [deg_mat3 deg_mat2 deg_mat1];
        c_clock = [deg_mat_index3 deg_mat_index2 deg_mat_index1];
    elseif(abs(p_inv(1,index_id_d(in))-gd_inv(3,6))<eps&&poly_cen_x>gd_inv(3,6)) %check if node is loacated on the left border of the imaging domain
        deg_mat_index1 = find(deg_mat>=90&deg_mat<=180);
        deg_mat_index2 = find(deg_mat>=-180&deg_mat<-90);
        deg_mat_index3 = find(deg_mat>=-90&deg_mat<=90);
        deg_mat1 = deg_mat(deg_mat_index1);
        deg_mat2 = deg_mat(deg_mat_index2);
        deg_mat3 = deg_mat(deg_mat_index3);
        [deg_mat1, c_clock1] = sort(deg_mat1);
        deg_mat_index1 = deg_mat_index1(c_clock1);
        [deg_mat2, c_clock2] = sort(deg_mat2);        
        deg_mat_index2 = deg_mat_index2(c_clock2);
        [deg_mat3, c_clock3] = sort(deg_mat3);        
        deg_mat_index3 = deg_mat_index3(c_clock3);
        deg_mat = [deg_mat2 deg_mat3 deg_mat1];
        c_clock = [deg_mat_index2 deg_mat_index3 deg_mat_index1];
    elseif(abs(p_inv(1,index_id_d(in))-gd_inv(4,6))<eps&&poly_cen_x<gd_inv(4,6)) %check if node is loacated on the right border of the imaging domain
        deg_mat_index1 = find(deg_mat<=90&deg_mat>=0);
        deg_mat_index2 = find(deg_mat>90&deg_mat<=180);
        deg_mat_index3 = find(deg_mat>=-90&deg_mat<0);
        deg_mat_index4 = find(deg_mat>-180&deg_mat<-90);
        deg_mat1 = deg_mat(deg_mat_index1);
        deg_mat2 = deg_mat(deg_mat_index2);
        deg_mat3 = deg_mat(deg_mat_index3);
        deg_mat4 = deg_mat(deg_mat_index4);
        [deg_mat1, c_clock1] = sort(deg_mat1);
        deg_mat_index1 = deg_mat_index1(c_clock1);
        [deg_mat2, c_clock2] = sort(deg_mat2);        
        deg_mat_index2 = deg_mat_index2(c_clock2);
        [deg_mat3, c_clock3] = sort(deg_mat3);        
        deg_mat_index3 = deg_mat_index3(c_clock3);
        [deg_mat4, c_clock4] = sort(deg_mat4);        
        deg_mat_index4 = deg_mat_index4(c_clock4);
        deg_mat = [deg_mat1 deg_mat2 deg_mat4 deg_mat3];
        c_clock = [deg_mat_index1 deg_mat_index2 deg_mat_index4 deg_mat_index3];
    else %check if node is inside the imaging domain
        if(poly_cen_x<=gd_inv(3,6))
            [~, c_clock] = sort(deg_mat,'ascend'); %ordering the angles to get the counterclockwise order
        else
            deg_mat(sign(deg_mat)<0) = deg_mat(sign(deg_mat)<0) + 360; %converting the negative angles to their positive equivalents (e.g. -30=270)
            [~, c_clock] = sort(deg_mat,'ascend'); %ordering the angles to get the counterclockwise order
        end
    end    
         
    x_point = x_point(c_clock); %ordering x-coord in counterclockwise order
    y_point = y_point(c_clock); %ordering y-coord in counterclockwise order
    point_mat = [x_point;y_point]; %storing the points in counterclockwise order
    flag_mat = flag_mat(c_clock);
    
    point_cell{in} = point_mat;
           
    %---------------------computing outward normal vectors and elementary lengths---------------------%
    %check if node in not located inside the imaging domain edges
    if(abs(abs(p_inv(1,index_id_d(in)))-abs(gd_inv(3,6)))<eps)%check if located on the left or right border
        norm_vec = zeros(2,size(c_clock,2)-1); %outward normal vectors
        dl_vec = zeros(1,size(c_clock,2)-1); %point distances
        int_el = zeros(1,size(c_clock,2)-1); %matrix the number of element if point is a centroid else 0
        
        for iin=1:size(norm_vec,2) %for n points n-1 edges exist
            x = point_mat(1,iin+1)-point_mat(1,iin); %x-coord of dl vector
            y = point_mat(2,iin+1)-point_mat(2,iin); %y-coord of dl vector
            dl_vec(iin) = sqrt(x^2+y^2);
            norm_vec(1,iin) = y/sqrt(x^2+y^2);
            norm_vec(2,iin) = -x/sqrt(x^2+y^2);
            if(flag_mat(iin)==0)
                int_el(iin) = flag_mat(iin+1);
            else
                int_el(iin) = flag_mat(iin);
            end
        end
    elseif(abs(abs(p_inv(2,index_id_d(in)))-abs(gd_inv(7,6)))<eps)%check if located on the upper or lower border
        norm_vec = zeros(2,size(c_clock,2)-1); %outward normal vectors
        dl_vec = zeros(1,size(c_clock,2)-1); %point distances
        int_el = zeros(1,size(c_clock,2)-1); %matrix the number of element if point is a centroid else 0
        
        for iin=1:size(norm_vec,2) %for n points n-1 edges exist
            x = point_mat(1,iin+1)-point_mat(1,iin); %x-coord of dl vector
            y = point_mat(2,iin+1)-point_mat(2,iin); %y-coord of dl vector
            dl_vec(iin) = sqrt(x^2+y^2);
            norm_vec(1,iin) = y/sqrt(x^2+y^2);
            norm_vec(2,iin) = -x/sqrt(x^2+y^2);
            if(flag_mat(iin)==0)
                int_el(iin) = flag_mat(iin+1);
            else
                int_el(iin) = flag_mat(iin);
            end
        end
    elseif(abs(abs(p_inv(1,index_id_d(in)))-abs(gd_inv(3,6)))>eps) %check if node is located inside on the edges of the imaging domain  
        norm_vec = zeros(2,size(c_clock,2)); %outward normal vectors
        dl_vec = zeros(1,size(c_clock,2)); %point distances
        int_el = zeros(1,size(c_clock,2));%matrix the number of element if point is a centroid else 0
    %     flag_mat = zeros(size(dl_vec)); %matrix containg the element in which dl belongseps

        for iin=1:size(norm_vec,2) %for n points n-1 edges exist
            if(iin+1<size(norm_vec,2)) %checking if the point is not the last of the polygon
                x = point_mat(1,iin+1)-point_mat(1,iin); %x-coord of dl vector
                y = point_mat(2,iin+1)-point_mat(2,iin); %y-coord of dl vector
                dl_vec(iin) = sqrt(x^2+y^2);
                norm_vec(1,iin) = y/sqrt(x^2+y^2);
                norm_vec(2,iin) = -x/sqrt(x^2+y^2);
                if(flag_mat(iin)==0)
                    int_el(iin) = flag_mat(iin+1);
                else
                    int_el(iin) = flag_mat(iin);
                end
            else %if we are on the last point of the point polygon we must return to the beginning for convexity
                x = point_mat(1,1)-point_mat(1,iin); %x-coord of dl vector
                y = point_mat(2,1)-point_mat(2,iin); %y-coord of dl vector
                dl_vec(iin) = sqrt(x^2+y^2);
                norm_vec(1,iin) = y/sqrt(x^2+y^2);
                norm_vec(2,iin) = -x/sqrt(x^2+y^2);
                if(flag_mat(iin)==0)
                    int_el(iin) = flag_mat(find(flag_mat~=0,1));
                else
                    int_el(iin) = flag_mat(iin);
                end
            end
        end
    end
    
    mr_cell{in,2} = norm_vec(1,:); %storing x-coord of norm vector
    mr_cell{in,3} = norm_vec(2,:); %storing y-coord of norm vector 
    mr_cell{in,4} = dl_vec; %storing the dl size
    int_el_cell{in} = int_el;
    %-------------------finding in which element each elementary length belongs---------------------%
    element_mat_b = zeros(3,size(dl_vec,2)); %b simplex
    element_mat_c = zeros(3,size(dl_vec,2)); %c simplex
    element_mat_n = zeros(3,size(dl_vec,2)); %nodes of element
    
    for inn=1:size(dl_vec,2)
        n = t_inv(1:3,int_el(inn)); %nodes 1-3 of element ie
        x = p_inv(1,n(1:3)); %x-cord of nodes 1-3 of element ie
        y = p_inv(2,n(1:3)); %y-cord of nodes 1-3 of element ie

        %---------area of element ie computation---------%
        D = det([ones(3,1) x' y']);
        %-------------cartesian-simplex------------------%
        b = [(y(2)-y(3))/D, (y(3)-y(1))/D, (y(1)-y(2))/D]';
        c = [(x(3)-x(2))/D, (x(1)-x(3))/D, (x(2)-x(1))/D]';
        
        element_mat_b(:,inn) = b;
        element_mat_c(:,inn) = c;
        element_mat_n(:,inn) = index_inv_D(index_inv(n));
    end
    
    mr_cell{in,5} = element_mat_b;
    mr_cell{in,6} = element_mat_c;
    mr_cell{in,7} = element_mat_n;
    waitbar(in/Nnodes_d,w_bar2,sprintf('Please wait: %g%%',0.01*round(10000*in/Nnodes_d)));
end
close(w_bar2);
polygon_area = mr_cell{:,1}; %inverted polygon area multiplied by total imaging area

%% CSI ALGORITHM
%vectorized way for multiple sources will be implemented
%data for each illumination are stored in columns of the matrices
%% -------------------Initialization--------------------%
ita_d = 0; %norm D factor
ita_s = 0; %norm S factor
%------------------------Ita_S norm factor--------------------------%
for is=1:Nsources
    ita_s = real((E_mes(:,is)')*E_mes(:,is))+ita_s;
end
ita_s = ita_s^(-1); %S norm factor
G_S = -2*ita_s*(Td\((L')*(MS'))); %Gs operator of CSI algorithm

w_0 = zeros(Nnodes_d,Nsources); % initialization parameter of CSI
E_t = zeros(Nnodes_d,Nsources); %incriment field of nodes belongin to imaging domain
U_0 = zeros(Nnodes_d);%auxialary parameters for contrast computation through the use of contast source
V_0 = zeros(Nnodes_d,1);%auxialary parameters for contrast computation through the use of contast source
U = zeros(Nnodes_d);%auxialary parameters for contrast computation through the use of contast source
V = zeros(Nnodes_d,1);%auxialary parameters for contrast computation through the use of contast source
ro_0 = zeros(size(E_mes));% inital data error
r_0 = zeros(Nnodes_d,Nsources);%initial domain error
dt_0 = zeros(Nnodes_d,Nsources); %Initial Polak Ribiere search direction
gt_0 = zeros(Nnodes_d,Nsources); %Initial cost functinal gradient

for is=1:Nsources
    w_0(:,is) = ((real(E_mes(:,is)'*(MS*L*(G_S*E_mes(:,is)))))/(real((MS*L*(G_S*(E_mes(:,is))))'*(MS*L*(G_S*E_mes(:,is))))))*G_S*E_mes(:,is);%CSI initial contrast source guess
    E_t(:,is) = MD*Ez_inc(:,is);%incriment fields of imaging domain D
    U_0 = (sparse(diag(E_t(:,is)+MD*L*w_0(:,is))))'*Td*(sparse(diag(E_t(:,is)+MD*L*w_0(:,is))))+U_0;%left matrix for contrast computation U_0*x=V_0
    V_0 = (sparse(diag(E_t(:,is)+MD*L*w_0(:,is))))'*Td*w_0(:,is)+V_0;%right matrix for contrast computation U_0*x=V_0
    ro_0(:,is) = E_mes(:,is)-MS*L*w_0(:,is); %initial data error
end

x_0 = full(sparse(U_0)\sparse(V_0)); %initial contrast guess

for is=1:Nsources
    r_0(:,is) = x_0.*E_t(:,is)-w_0(:,is)+x_0.*(MD*L*w_0(:,is)); %initial domain error
    gt_0(:,is) = -2*ita_s*Td\(L'*MS'*E_mes(:,is));
end
%%% ----------------Iterative procedure------------------%
Noiterations = 1024; %number of iterations

%-----previous initial values and array preallocation-----%
r_pre = r_0; %previous value of domain error
ro_pre = ro_0; %previous value of data error
dt_pre = dt_0; %previous Polak Ribiere search direction
gt_pre = gt_0; %previous cost functional gradient
w_pre = w_0; %previous contrast source
x_pre = x_0; %previous contrast estimate
dt_cur = zeros(Nnodes_d,Nsources); %Polak Ribiere search direction
gt_cur = zeros(Nnodes_d,Nsources); %cost functional gradient
w_cur = zeros(Nnodes_d,Nsources); %contrast sources current
x_cur = zeros(Nnodes_d,1); %contrast
r_cur = zeros(Nnodes_d,Nsources); %current domain error
ro_cur = zeros(size(E_mes));%current data error
count=0;
inverse1 = Td\(L'*MS'); %inverse computed only once for speed improvement error e-14
inverse2 = Td\(L'*MD'); %inverse computed only once for speed improvement error e-14
MSL = MS*L; %computating MS*L product only once for speed boost
MDL = MD*L; %computating MD*L product only once for speed boost
% a_step = zeros(Noiterations,Nsources); %conjugated gradient step size
% a_stepm=zeros(Nsources,1);
U_mat4 = mat2cell(repmat(Td,1,Nsources),Nnodes_d,Nnodes_d*ones(1,Nsources));
cost_fun = zeros(Noiterations,1); 
delta_cell = cell(Nnodes_d,1); %cell containg delta term for gradient vector computation
gt_x_pre = ones(Nnodes_d,1); %contrast gradient previous value set to 1(unless 0/0=NaN) for division in step 1 gt_x_pre_0 has no meaning
dt_x_pre = 0; %polak ribiere contrast previous value

subindex = @(A,r) A(r); %fuction definition to get x(n) in vectorised way

save('mesh.mat','p_inv','e_inv','t_inv')
save('forward.mat','p','e','t','Ezs0','Ez')

w_bar=waitbar(0,'Please Wait','Name','FEM-MRCSI Reconstruction');

b_mat_step = zeros(Noiterations,1);
cost1 = zeros(Noiterations,1);
mr = zeros(Noiterations,1);
delta_m = zeros(Noiterations,1);
mr_xa_m = zeros(Noiterations,1);
for ir=1:Noiterations %number of iterations
    %------------------------Ita_D norm factor--------------------------%
    ita_d = 0;
%     for is=1:Nsources
%         ita_d = real((x_pre.*E_t(:,is))'*Td*(x_pre.*E_t(:,is)))+ita_d;
%     end
%     ita_d = ita_d^(-1);
    ita_d1 = (repmat(x_pre,1,Nsources).*E_t(:,1:Nsources))';
    ita_d2 = (repmat(x_pre,1,Nsources).*E_t(:,1:Nsources));
    ita_d = real(ita_d1*Td*ita_d2);
    ita_d = sum(diag(ita_d));
    ita_d = ita_d^(-1); %D norm factor
    %----------Polak Ribiere,gradient,contrast source update------------%
    %------------gradient-------------%
    gt_cur(:,1:Nsources) = -2*ita_s*inverse1*ro_pre(:,1:Nsources)-2*ita_d*r_pre(:,1:Nsources)+2*ita_d*inverse2*(diag(x_pre)')*Td*r_pre(:,1:Nsources);%\((eye(Nnodes_d,Nnodes_d)-L'*MD'*diag(x_pre)')*Td*r_pre(:,is)); %cost functional gradient computation
    
    %----------Polak Ribiere----------%
    dt_cur_nom1 = (gt_cur(:,1:Nsources)-gt_pre(:,1:Nsources))'*Td*gt_cur(:,1:Nsources); %auxiliary matrix that has the inner products of repetition in column 1 and extra trash information
    dt_cur_denom1 = real(gt_pre(:,1:Nsources)'*Td*gt_pre(:,1:Nsources)); %auxiliary matrix that has the norms of repetition in column 1 and extra trash information
    dt_cur_inter = diag(dt_cur_nom1./dt_cur_denom1); %16 different numbers each for every source that are being multiplied with the previous Polak Ribiere taking the diag because we need the elements with the same index first row first column and so on
    dt_cur = repmat(dt_cur_inter.',Nnodes_d,1).*dt_pre(:,1:Nsources); %turing mat3 into a row that has the multipliers filling a matrix with Nnodes_d times of rows o the matrix mat3 and then multiplying the matrices column by column
    dt_cur(:,1:Nsources) = -gt_cur(:,1:Nsources)+dt_cur(:,1:Nsources); %final Polak Ribiere calculation
    
    %------------step size------------%
    a_step_nom1 = ita_s*((MSL*dt_cur(:,1:Nsources))'*ro_pre(:,1:Nsources)); %part 1 step size nominator
    a_step_nom2 = (dt_cur(:,1:Nsources)-repmat(x_pre,1,Nsources).*(MDL*dt_cur(:,1:Nsources)))';
    a_step_nom3 = ita_d*(a_step_nom2*Td*r_pre(:,1:Nsources)); %part 2 step size nominator using repmat to multpily x_pre in order to have coherent matirx dimensions
    a_step_nom = diag(a_step_nom1)+diag(a_step_nom3); %step size nominator taking the diagonal because part 1 2 are square matrices and we need the division of row 3 with column three for example
    a_step_denom1 = ita_s*real((MSL*dt_cur(:,1:Nsources))'*(MSL*dt_cur(:,1:Nsources)));
    a_step_denom2 = ita_d*real(a_step_nom2*Td*(dt_cur(:,1:Nsources)-repmat(x_pre,1,Nsources).*(MDL*dt_cur(:,1:Nsources))));
    a_step_denom = diag(a_step_denom1)+diag(a_step_denom2);
    a_step = a_step_nom./a_step_denom; %final step size calculation
    
    %---------contrast source---------%
    update = a_step.'; %a_step is of size Nsources*1 taking its transpose
    update = update(ones(1,Nnodes_d),:); %creating a matrix in first row 1 Nnodes_d times the step size for source No 1 etc
    w_cur(:,1:Nsources) = w_pre(:,1:Nsources)+update.*dt_cur(:,1:Nsources); %calculating contrast source
    
    %------------------Contrast,Data Error Update-----------------------%
    %------------------Contrast update------------------%
    %Calculating in vectorized way arrays E_t'*Td*E_t and E_t'*Td*w_cur
    U_mat1 = reshape(E_t+MD*L*w_cur,1,[]);
    U_mat2 = spdiags(U_mat1.',0,Nnodes_d*Nsources,Nnodes_d*Nsources);
    U_mat3 = mat2cell(U_mat2,Nnodes_d*ones(1,Nsources),Nnodes_d*ones(1,Nsources));
    U_mat3 = reshape(U_mat3,1,[]); %reshapre the matrix taking the columns as elements and concatenating them
    cell_index = 1:(Nsources+1):Nsources^2;
    U_mat3 = U_mat3(cell_index); %taking the diagonal of cell matrix
%     U_mat4 = mat2cell(repmat(Td,1,Nsources),Nnodes_d,Nnodes_d*ones(1,Nsources)); %it will be taken out for the loop as it needs to be calculated only once and is the same
    U_mat5 = cellfun(@ctranspose,U_mat3,'UniformOutput',false);
    U_mat6 = cellfun(@mtimes,U_mat5,U_mat4,'UniformOutput',false);
    U_mat7 = cellfun(@mtimes,U_mat6,U_mat3,'UniformOutput',false); 
    
    V_mat1 = num2cell(w_cur,1);
    V_mat2 = cellfun(@mtimes,U_mat6,V_mat1,'UniformOutput',false);
    
    ro_cur(:,1:Nsources) = E_mes(:,1:Nsources)-MSL*w_cur(:,1:Nsources); %data error
    
    U = 0;
    V = 0;
    
    for k=1:Nsources
        U = [U_mat7{k}]+U; %E_t'*Td*E_t array
        V = [V_mat2{k}]+V; %E_t'*Td*w_cur array
    end;

    x_cur_nonmr = full(sparse(U)\sparse(V)); %contrast update
    
    r_cur_nonmr = repmat(x_cur_nonmr,1,Nsources).*E_t(:,1:Nsources)-w_cur+repmat(x_cur_nonmr,1,Nsources).*(MDL*w_cur(:,1:Nsources)); %domain error
    
    %-----------------Multiplicative Regularization---------------------%
    %-------------------Delta term-------------------%
    %delta term when using the previous contrast and contrast sources values
%     nom_delta = sum(diag(real(r_pre'*Td*r_pre)));
%     ita_ddelta1 = (repmat(x_pre,1,Nsources).*E_t(:,1:Nsources))';
%     ita_ddelta2 = (repmat(x_pre,1,Nsources).*E_t(:,1:Nsources));
%     ita_ddelta = real(ita_ddelta1*Td*ita_ddelta2);
%     ita_ddelta = sum(diag(ita_ddelta));
%     ita_ddelta = ita_ddelta^(-1);
%     delta_cur = (nom_delta*ita_ddelta)/inv_area_m;
    %delta term when using the current contrast source term
    nom_delta = sum(diag(real(r_cur_nonmr'*Td*r_cur_nonmr))); %FD term with constrast with no regularization
%     ita_ddelta1 = (repmat(x_cur_nonmr,1,Nsources).*E_t(:,1:Nsources))';
%     ita_ddelta2 = (repmat(x_cur_nonmr,1,Nsources).*E_t(:,1:Nsources));
%     ita_ddelta = real(ita_ddelta1*Td*ita_ddelta2);
%     ita_ddelta = sum(diag(ita_ddelta));
%     ita_ddelta = ita_ddelta^(-1);
    delta_cur = (nom_delta*ita_d)/inv_area_m; %delta regularization term
    
    %-----------Cost functional with no MR-----------%
    nom_1 = sum(diag(real(ro_cur'*ro_cur)));
    nom_2 = sum(diag(real(r_cur_nonmr'*Td*r_cur_nonmr)));
%     ita_dd1 = (repmat(x_cur_nonmr,1,Nsources).*E_t(:,1:Nsources))';
%     ita_dd2 = (repmat(x_cur_nonmr,1,Nsources).*E_t(:,1:Nsources));
%     ita_dd = real(ita_dd1*Td*ita_dd2);
%     ita_dd = sum(diag(ita_dd));
%     ita_dd = ita_dd^(-1);
    cost_fun_non_mr = nom_1*ita_s+nom_2*ita_d;
    
    %------------Gradient Computation------------%
%     for in=1:Nnodes_d
%         mr_cell{in,8} = x_pre(mr_cell{in,7}); %annotating the x_pre value in the mr_cell
%         mr_cell{in,9} = x_cur_nonmr(mr_cell{in,7}); %annotating the x_cur value in the mr_cell
%         delta_cell{in,1} = delta_cur*ones(size(mr_cell{in,4})); %annotating the delta value in cell for gradient vectorised computation
%     end
    x_pre_cell = num2cell(repmat(x_pre,1,Nnodes_d),1)'; %creating cell array with x_pre array in every element
    x_cur_nonmr_cell = num2cell(repmat(x_cur_nonmr,1,Nnodes_d),1)'; %creating cell array with contrast array in every element
    mr_cell(:,8) = cellfun(subindex,x_pre_cell,mr_cell(:,7),'UniformOutput',false); %taking the verttex polygon nodes for x_pre
    mr_cell(:,9) = cellfun(subindex,x_cur_nonmr_cell,mr_cell(:,7),'UniformOutput',false); %taking the verttex polygon nodes for contrast
    delta_cell = num2cell(delta_cur*ones(Nnodes_d,1),2); %creating cell array with delta term for adding in b_n term
    
    %b_n(r) term computation inside the integral
    gt_x_cell_b1 = cellfun(@times,mr_cell(:,8),mr_cell(:,5),'UniformOutput',false); %multiplying x_pre with bi
    gt_x_cell_c1 = cellfun(@times,mr_cell(:,8),mr_cell(:,6),'UniformOutput',false); %multiplying x_pre with ci
    gt_x_cell_b1 = cellfun(@sum,gt_x_cell_b1,'UniformOutput',false); %taking the sum of x_pre*bi
    gt_x_cell_c1 = cellfun(@sum,gt_x_cell_c1,'UniformOutput',false); %taking the sum of x_pre*ci
    gt_x_cell_b1_con = cellfun(@conj,gt_x_cell_b1,'UniformOutput',false); %taking the x_pre*bi conjugate in order to compute the norm
    gt_x_cell_c1_con = cellfun(@conj,gt_x_cell_c1,'UniformOutput',false); %taking the x_pre*ci conjugate in order to compute the norm
    gt_x_cell_b1 = cellfun(@times,gt_x_cell_b1,gt_x_cell_b1_con,'UniformOutput',false); %multiplying x_pre*bi with the conjugate to get the norm
    gt_x_cell_c1 = cellfun(@times,gt_x_cell_c1,gt_x_cell_c1_con,'UniformOutput',false); %multiplying x_pre*ci with the conjugate to get the norm
    gt_x_cell_delta = cellfun(@plus,gt_x_cell_b1,gt_x_cell_c1,'UniformOutput',false); %adding x_pre*bi and x_pre*ci
    gt_x_cell_delta = cellfun(@plus,gt_x_cell_delta,delta_cell,'UniformOutput',false);%adding delta term
    gt_x_cell_delta = cellfun(@mtimes,num2cell(inv_area*ones(Nnodes_d,1),2),gt_x_cell_delta,'UniformOutput',false);
    gt_x_cell_delta = cellfun(@rdivide,num2cell(ones(Nnodes_d,1),2),gt_x_cell_delta,'UniformOutput',false); %taking the inverse to calculate the b_n without multiplying with the area of domains
    
    %grad(contrast) comuptation inside the integral
    gt_x_cell_b2 = cellfun(@times,mr_cell(:,9),mr_cell(:,5),'UniformOutput',false); %multiplying x_cur_nonmr with bi
    gt_x_cell_c2 = cellfun(@times,mr_cell(:,9),mr_cell(:,6),'UniformOutput',false); %multiplying x_cur_nonmr with ci
    gt_x_cell_b2 = cellfun(@sum,gt_x_cell_b2,'UniformOutput',false); %taking the sum of x_cur_nonmr*bi
    gt_x_cell_c2 = cellfun(@sum,gt_x_cell_c2,'UniformOutput',false); %taking the sum of x_cur_nonmr*ci
    
    %integral computation
    %x component
    gt_x_cur1 = cellfun(@times,gt_x_cell_delta,gt_x_cell_b2,'UniformOutput',false); %multiplying b_n with sum(x_cur_nonmr*bi)
    gt_x_cur1 = cellfun(@times,gt_x_cur1,mr_cell(:,2),'UniformOutput',false); %multiplying with the x-coord of normal vector
    gt_x_cur1 = cellfun(@times,gt_x_cur1,mr_cell(:,4),'UniformOutput',false); %multiplying with dl
    gt_x_cur1 = cellfun(@mtimes,gt_x_cur1,mr_cell(:,1),'UniformOutput',false); %multiplying with polygon area 
    gt_x_cur1 = cellfun(@sum,gt_x_cur1,'UniformOutput',false); %taking the sum of the line integrals of every edge of the polygon
    %y component
    gt_x_cur2 = cellfun(@times,gt_x_cell_delta,gt_x_cell_c2,'UniformOutput',false); %multiplying b_n with sum(x_cur_nonmr*ci)
    gt_x_cur2 = cellfun(@times,gt_x_cur2,mr_cell(:,3),'UniformOutput',false); %multiplying with the y-coord of normal vector
    gt_x_cur2 = cellfun(@times,gt_x_cur2,mr_cell(:,4),'UniformOutput',false); %multiplying with dl
    gt_x_cur2 = cellfun(@mtimes,gt_x_cur2,mr_cell(:,1),'UniformOutput',false); %multiplying with polygon area 
    gt_x_cur2 = cellfun(@sum,gt_x_cur2,'UniformOutput',false); %taking the sum of the line integrals of every edge of the polygon
    
    gt_x_cur3 = cell2mat(gt_x_cur1) + cell2mat(gt_x_cur2); %adding x and y components of the gradients
    gt_x_cur3 = -2*cost_fun_non_mr*gt_x_cur3; %multiplying with the polygon area
    
    gt_x_denom1 = E_t+MDL*w_cur; %computing total field
%     gt_x_denom1 = real(sum(diag(gt_x_denom1'*Td*gt_x_denom1))); %taking the meter of total field in every node
    gt_x_denom1 = sum(abs(gt_x_denom1).^2,2);

    %b_n term in every triangle
    b_step2_x_pre = x_pre(simplex_tri_n); %x_pre value in every triangle
    b_step2_co_x_pre_b = b_step2_x_pre.*simplex_tri_b; %multiplying x_pre value with simplex b in every triangle
    b_step2_co_x_pre_c = b_step2_x_pre.*simplex_tri_c; %multiplying x_pre value with simplex c in every triangle
    b_step2_co_x_pre_b = sum(b_step2_co_x_pre_b,2); %taking the sum in every element of the product previous contrast*simplex b
    b_step2_co_x_pre_c = sum(b_step2_co_x_pre_c,2); %taking the sum in every element of the product previous contrast*simplex c
    b_step2_co_x_pre_b = b_step2_co_x_pre_b.*conj(b_step2_co_x_pre_b); %taking the norm of x component since contrast is complex
    b_step2_co_x_pre_c = b_step2_co_x_pre_c.*conj(b_step2_co_x_pre_c); %taking the norm of y component since contrast is complex
    b_step2_co_x_pre = 1./sqrt(inv_area*(b_step2_co_x_pre_b + b_step2_co_x_pre_c + delta_cur)); %computing the b_n(r) term for every triangle
    
    %mr term with the contrast solution without applying the regularization
    mr_xa1 = x_cur_nonmr(simplex_tri_n);
    mr_xa2 = mr_xa1.*simplex_tri_b;
    mr_xa3 = mr_xa1.*simplex_tri_c;
    mr_xa2 = sum(mr_xa2,2);
    mr_xa3 = sum(mr_xa3,2);
    mr_xa2 = mr_xa2.*conj(mr_xa2);
    mr_xa3 = mr_xa3.*conj(mr_xa3);
%     mr_xa_co1 = sum(inv_tri_area.*(mr_xa2+mr_xa3).*(b_step2_co_x_pre.^2));
    mr_xa_co = mr_xa2 + mr_xa3 + delta_cur;
%     mr_xa_co2 = delta_cur*sum(inv_tri_area.*(b_step2_co_x_pre.^2));
%     mr_xa= mr_xa_co1 + mr_xa_co2;
    mr_xa = sum(inv_tri_area.*mr_xa_co.*(b_step2_co_x_pre.^2));
    mr_xa_m(ir) = mr_xa;
    
    gt_x_cur4 = r_cur_nonmr.*conj(E_t+MDL*w_cur);
    gt_x_cur4 = 2*ita_d*mr_xa*sum(gt_x_cur4,2);
    
    gt_x_cur = gt_x_cur3;% + gt_x_cur4;
    
    %----------Polak Ribiere Computation----------%
    dt_x_nom = (gt_x_cur-gt_x_pre)'*Td*gt_x_cur; %Polak Ribiere of constrast
    dt_x_denom = real(gt_x_pre'*Td*gt_x_pre);
    dt_x_cur = -gt_x_cur./gt_x_denom1 + (dt_x_nom/dt_x_denom)*dt_x_pre;
    
    %------------Step Size Computation------------%
    %polynomial coefficient from the term with no grad terms
    b_step_co0 = repmat(dt_x_cur,1,Nsources).*(E_t+MDL*w_cur); 
    b_step1_co1 = b_step_co0'*Td*b_step_co0;
    b_step1_co1 = real(ita_d*sum(diag(b_step1_co1))); %first coefficient of product polynomial
    
    b_step1_co2 = r_cur_nonmr'*Td*b_step_co0;
    b_step1_co2 = 2*ita_d*sum(real(diag(b_step1_co2))); %second coefficient of product polynomial
    
    %grad contrast and contrast Ribiere in every triangle
    b_step2_x_cur = x_cur_nonmr(simplex_tri_n); %x_cur value in every triangle
    b_step2_dt_x_cur = dt_x_cur(simplex_tri_n); %dt_x_xur value in every value
    
    b_step2_co_x_cur_b = b_step2_x_cur.*simplex_tri_b; %multiplying x_cur value with simplex b in every triangle
    b_step2_co_x_cur_c = b_step2_x_cur.*simplex_tri_c; %multiplying x_cur value with simplex c in every triangle
    b_step2_co_dt_x_cur_b = b_step2_dt_x_cur.*simplex_tri_b; %multiplyin d_x_cur value with simplex b in every triangle
    b_step2_co_dt_x_cur_c = b_step2_dt_x_cur.*simplex_tri_c; %multiplyin d_x_cur value with simplex c in every triangle
    
    b_step2_co_dt_x_cur_b = sum(b_step2_co_dt_x_cur_b,2); %taking the sum in every element of the product polak ribiere*simplex b
    b_step2_co_dt_x_cur_c = sum(b_step2_co_dt_x_cur_c,2); %taking the sum in every element of the product polak ribiere*simplex c
    b_step2_co_x_cur_b = sum(b_step2_co_x_cur_b,2); %taking the sum in every element of the product contrast*simplex b
    b_step2_co_x_cur_c = sum(b_step2_co_x_cur_c,2); %taking the sum in every element of the product contrast*simplex b
    
    b_step2_x_cur_xvec = b_step2_co_x_pre.*b_step2_co_x_cur_b; %x-coord MR contrast grad for every triangle
    b_step2_x_cur_yvec = b_step2_co_x_pre.*b_step2_co_x_cur_c; %y-coord MR contrast grad for every triangle
    b_step2_dt_cur_xvec = b_step2_co_x_pre.*b_step2_co_dt_x_cur_b; %x-coord MR Polak Ribiere for every triangle
    b_step2_dt_cur_yvec = b_step2_co_x_pre.*b_step2_co_dt_x_cur_c; %y-coord MR Polak Ribiere for every triangle
    
    b_step2_co1 = sum(inv_tri_area.*(b_step2_dt_cur_xvec.*conj(b_step2_dt_cur_xvec) + b_step2_dt_cur_yvec.*conj(b_step2_dt_cur_yvec))); %abs(bn*grad(d_x)) 
    b_step2_co2 = 2*real(sum(inv_tri_area.*(b_step2_x_cur_xvec.*conj(b_step2_dt_cur_xvec) + b_step2_x_cur_yvec.*conj(b_step2_dt_cur_yvec)))); %<bn*grad(x_a),bn*grad(d_x)>
    
    step_poly = conv(full([b_step1_co1 b_step1_co2 cost_fun_non_mr]),full([b_step2_co1 b_step2_co2 mr_xa])); %computing the polynomial frow which the 
    step_poly = polyder(step_poly); %calculating polynomial derivative
    
    b_step = roots(step_poly); %computing the roots of the step polynomial
    b_step_place = find(abs((abs(imag(b_step)))-min(abs(imag(b_step))))<eps);
    b_step = real(b_step(b_step_place));
%     b_step = real(b_step(abs(imag(b_step)) == min(abs(imag(b_step))))); %taking the real root of the polynomial
    
    %---------------Contrast Update---------------%
    x_cur = x_cur_nonmr + b_step*dt_x_cur;
        
    %------------Positivity constraint------------%
    x_cur_re = real(x_cur);
    x_cur_im = imag(x_cur);
    
    check_re = find(real(er_air*x_cur+er_air)<1);
    if(isempty(check_re)==0)
        x_cur_re(check_re) = (1-er_air)/er_air;
    end
     
    check_im = find(imag(er_air*x_cur+er_air)>0);
    if(isempty(check_im)==0)
        x_cur_im(check_im) = 0;
    end
    
    x_cur = complex(x_cur_re,x_cur_im);
          
    %---------------Domain Error Upadate-----------------%
    r_cur = repmat(x_cur,1,Nsources).*E_t(:,1:Nsources)-w_cur+repmat(x_cur,1,Nsources).*(MDL*w_cur(:,1:Nsources)); %domain error
    
    %------------------Previous values annotation-----------------------%
    dt_pre = dt_cur; %Polak Ribiere annotation
    w_pre = w_cur; %contrast source annotation
    r_pre = r_cur; %domain error annotation
    ro_pre = ro_cur; %data erro annotation
    gt_pre = gt_cur;
    x_pre = x_cur;
    gt_x_pre = gt_x_cur;
    dt_x_pre = dt_x_cur;
     
    count=count+1
    
    %------------------Cost functional-----------------------%
    %FS FD term
    nom_1 = sum(diag(real(ro_cur'*ro_cur)));
    nom_2 = sum(diag(real(r_cur'*Td*r_cur)));
    ita_dd1 = (repmat(x_pre,1,Nsources).*E_t(:,1:Nsources))';
    ita_dd2 = (repmat(x_pre,1,Nsources).*E_t(:,1:Nsources));
    ita_dd = real(ita_dd1*Td*ita_dd2);
    ita_dd = sum(diag(ita_dd));
    ita_dd = ita_dd^(-1);
    %MR term
    mr_term1 = x_cur(simplex_tri_n);
    mr_term2 = mr_term1.*simplex_tri_b;
    mr_term3 = mr_term1.*simplex_tri_c;
    mr_term2 = sum(mr_term2,2);
    mr_term3 = sum(mr_term3,2);
    mr_term2 = mr_term2.*conj(mr_term2);
    mr_term3 = mr_term3.*conj(mr_term3);
    mr_term_co = mr_term2 + mr_term3 + delta_cur;
%     mr_term_co1 = sum(inv_tri_area.*(mr_term2+mr_term3).*(b_step2_co_x_pre.^2));
%     mr_term_co2 = delta_cur*sum(inv_tri_area.*(b_step2_co_x_pre.^2));
%     mr_term = mr_term_co1 + mr_term_co2;
    mr_term = sum(inv_tri_area.*mr_term_co.*(b_step2_co_x_pre.^2));
    
    cost_fun(ir) = mr_term*(nom_1*ita_s+nom_2*ita_d);
    mr(ir)  = mr_term;
    b_mat_step(ir)=b_step;
    cost1(ir)=nom_1*ita_s+nom_2*ita_dd;
    delta_m(ir) = delta_cur;
    
    waitbar(ir/Noiterations,w_bar,sprintf('Please wait: %g%%',0.01*round(10000*ir/Noiterations)));
    save('solution1.mat','x_cur','x_cur_nonmr','cost_fun','count','w_cur')
end
close(w_bar);
%--------------Contrast values annotation to inv mesh-------------------%
er_b = er_air;
contrast = zeros(Nnodes_inv,1);
contrast(index_id_d) = x_cur; 
er_reconstructed = er_b*contrast+er_b;
 
%% ERROR PROFILE ESTIMATION
%------------------Exact Profile------------------%
t_id = find(cen(1,:)>gd_inv(3,6)&cen(1,:)<gd_inv(4,6)&cen(2,:)>gd_inv(7,6)&cen(2,:)<gd_inv(9,6)); %forward mesh triangles in imaging domain
d_domain_t = t(:,t_id); %triangle of forward mesh inside imaging domain
region_domain_d = regions(t_id); %region matrix of domain triangles
node_domain = unique(t(1:3,t_id)); %nodes in imaging domain
node_domain_id = zeros(Nnodes,1); %nodes domain id
node_domain_id(node_domain) = 1:size(node_domain,1);
d_domain_p = p(:,node_domain); % forward mesh nodes in imaging domain
d_domain_t(1:3,:) = node_domain_id(d_domain_t(1:3,:)); %triangle matrix imaging domain
er_exact = ones(Nnodes,1);
er_exact(unique(t(1:3,regions==1))) = er_air;
er_exact(unique(t(1:3,regions==6))) = er1;
er_exact(unique(t(1:3,regions==5))) = er2;
er_exact = er_exact(node_domain);

d_domain_t_inv = t_inv(:,regions_inv==5);
d_domain_p_inv = p_inv(:,unique(d_domain_t_inv(1:3,:)));
d_domain_t_inv(1:3,:) = index_inv_D(index_inv((d_domain_t_inv(1:3,:))));

Nx = 150;%Nnodes_d;
Ny = 150;%Nnodes_d;
x = linspace(gd_inv(3,6)+l/20,gd_inv(4,6)-l/20,Nx);
y = linspace(gd_inv(7,6)+l/20,gd_inv(9,6)-l/20,Ny);

er_exact_xy = tri2grid(d_domain_p,d_domain_t,er_exact,x,y);
er_recon_xy = tri2grid(d_domain_p_inv,d_domain_t_inv,er_b*x_cur+er_b,x,y);

error_prof = sqrt(sum(sum(abs(er_recon_xy-er_exact_xy).^2)))/sqrt(sum(sum(abs(er_exact_xy).^2)));

%%
ttime=toc 
save('mrsquare.mat','-v7.3')
system('shutdown -s')