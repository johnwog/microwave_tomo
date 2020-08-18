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
f = 2e9; %frequency
w = 2*pi*f; %angular frequency
l = c0/f; %wavelength
k0 = w*sqrt(e0*m0); %vacuum wavenumber
E0 = 10; %incriment wave amplitude
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
b = 6.5; %PML a parameter imaginary part 
         %for thickness d=0.15l reflectance factor 3.1*e-5

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

%scatterer
er = 2-1i; %complex relative electric constant
ers = diag([er er er]);%electic constant tensor
mrs = eye(3);%megnetic constant tensor
mrs_inv = inv(mrs);

%backround medium
kb = k0; %backround medium wavenumber
er_b = 1; %bakcround relative dielectirc constant
mr_b = 1; %bakcround relative permeability constant

%% PML AREA AND SCATTERER GEOMETRY
%-------------------GEOMETRY DEFINITION-----------------------%
%PML parmeters
d = 0.15*l; %PML thickness
dis = 0.15*l; %distance from PML centered
length = 2*d+2*dis+l; %PML length by changing dis and d we control the size of the domain
param = l; %parameter umlaut profile size is calculated

%geometry matrix scatterer centered in the start of the coordinate system
gd=[3 4 -dis-d-0.5*l      -dis-0.5*l        -dis-0.5*l        -dis-d-0.5*l    -dis-0.5*l-d   -dis-0.5*l-d     dis+0.5*l+d    dis+0.5*l+d; %left vertical PML (R1)
    3 4  dis+0.5*l         dis+0.5*l+d       dis+d+0.5*l       dis+0.5*l      -dis-0.5*l-d   -dis-0.5*l-d     dis+0.5*l+d    dis+0.5*l+d; %right vertical PML (R2)
    3 4 -dis-d-0.5*l       dis+d+0.5*l       dis+d+0.5*l      -dis-d-0.5*l    -dis-0.5*l-d   -dis-0.5*l-d    -dis-0.5*l     -dis-0.5*l;   %lower horizontal PML (R3) 
    3 4 -dis-d-0.5*l       dis+d+0.5*l       dis+d+0.5*l      -dis-d-0.5*l     dis+0.5*l      dis+0.5*l       dis+d+0.5*l    dis+0.5*l+d; %upper horizontal PML (R4)
    3 4 -0.5*l-dis         0.5*l+dis         0.5*l+dis        -0.5*l-dis      -0.5*l-dis     -0.5*l-dis       0.5*l+dis       0.5*l+dis;  %between medium (R5)
    %%Umluaut profile broken in different shapes 
    1 0  0                 (0.2/1.7)*param       0                0               0              0              0                0;     %center circle for umlaut profile with small radius (C1)
    1 0  0                 2*(0.2/1.7)*param     0                0               0              0              0                0;     %center circle for umlaut profile with large radius (C2)
    3 4 -2*(0.2/1.7)*param -(0.2/1.7)*param  -(0.2/1.7)*param -2*(0.2/1.7)*param  0              0      2*(0.2/1.7)*param+param/50 2*(0.2/1.7)*param+param/50; %left rectangular part of umlaut (R6)
    3 4  (0.2/1.7)*param     2*(0.2/1.7)*param  2*(0.2/1.7)*param   (0.2/1.7)*param    0  0 2*(0.2/1.7)*param+param/50 2*(0.2/1.7)*param+param/50; %right rectangular part of umlaut (R7)
    %3 4 -(0.2/1.7)*param (0.2/1.7)*param (0.2/1.7)*param -(0.2/1.7)*param 0 0 2*(0.2/1.7)*param 2*(0.2/1.7)*param; %middle rectangular part of umlaut (R8)
    1   -(0.2/1.7)*param     2*(0.2/1.7)*param+(5.5/50)*param     0.5*(0.2/1.7)*param    0    0    0    0    0    0; %upper left circle from umlaut (C3)
    1    (0.2/1.7)*param     2*(0.2/1.7)*param+(5.5/50)*param     0.5*(0.2/1.7)*param    0    0    0    0    0    0; %upper right circle from umlaut (C3)
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
d1 = decsg(gd,'C2-C1+R1+R2+R3+R4+R5+R6+R7+C3+C4',[abs('R1')' abs('R2')' abs('R3')' abs('R4')' abs('R5')' abs('R6')' abs('R7')' abs('C1')' abs('C2')' abs('C3')' abs('C4')']);  

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

    it=[find(t(4,:)==2) find(t(4,:)==3) find(t(4,:)==4) find(t(4,:)==5) find(t(4,:)==6) find(t(4,:)==7) find(t(4,:)==8) find(t(4,:)==9) find(t(4,:)==10) find(t(4,:)==11) find(t(4,:)==12) find(t(4,:)==13) find(t(4,:)==15) find(t(4,:)==16) find(t(4,:)==18)];

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
regions(t(4,:)==1) = 1;
regions(t(4,:)==14) = 1;
regions(t(4,:)==17) = 1;

%----Horizontal PML-----%
regions(t(4,:)==3) = 2;
regions(t(4,:)==9) = 2;

%-----Vertical PML------%
regions(t(4,:)==7) = 3;
regions(t(4,:)==8) = 3;

%----PML Intersection---%
regions(t(4,:)==2) = 4;
regions(t(4,:)==4) = 4;
regions(t(4,:)==5) = 4;
regions(t(4,:)==6) = 4;

%-------Scatterer-------%
regions(t(4,:)==10) = 5;
regions(t(4,:)==11) = 5;
regions(t(4,:)==12) = 5;
regions(t(4,:)==13) = 5;
regions(t(4,:)==15) = 5;
regions(t(4,:)==16) = 5;
regions(t(4,:)==18) = 5;
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
%               fill([xii(1) xii(2) xii(3)],[yii(1) yii(2) yii(3)],col2(iin,:)); 
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
Nsources = 16; %Number of sources

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
                    NoNonZeros = NoNonZeros+1; %increase counter by one for trimming
                    
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
        Ez(in,1:Nsources) = Ezs0(in,1:Nsources)+E0.*exp(-1i*k0.*(p(1,in).*coord(1,1:Nsources).*sind(theta)+p(2,in).*coord(2,1:Nsources).*cosd(theta)));
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
gd_inv=[3 4 -dis-d-0.5*l  -dis-0.5*l    -dis-0.5*l    -dis-d-0.5*l  -dis-0.5*l-d   -dis-0.5*l-d    dis+0.5*l+d   dis+0.5*l+d; %left vertical PML (R1)
        3 4  dis+0.5*l     dis+0.5*l+d   dis+d+0.5*l   dis+0.5*l    -dis-0.5*l-d   -dis-0.5*l-d    dis+0.5*l+d   dis+0.5*l+d; %right vertical PML (R2)
        3 4 -dis-d-0.5*l   dis+d+0.5*l   dis+d+0.5*l  -dis-d-0.5*l  -dis-0.5*l-d   -dis-0.5*l-d   -dis-0.5*l    -dis-0.5*l;   %lower horizontal PML (R3) 
        3 4 -dis-d-0.5*l   dis+d+0.5*l   dis+d+0.5*l  -dis-d-0.5*l   dis+0.5*l      dis+0.5*l      dis+d+0.5*l   dis+0.5*l+d; %upper horizontal PML (R4)
        3 4 -0.5*l-dis     0.5*l+dis     0.5*l+dis    -0.5*l-dis    -0.5*l-dis     -0.5*l-dis      0.5*l+dis     0.5*l+dis;   %between medium (R5)
        3 4 -0.5*l         0.5*l         0.5*l        -0.5*l        -0.5*l         -0.5*l          0.5*l         0.5*l;]';    %imaging domain square (R6)
%-------------------Decompesed Geometry------------------------%
%decompesed geometry matrix
d1_inv = decsg(gd_inv,'R1+R2+R3+R4+R5+R6',[abs('R1')' abs('R2')' abs('R3')' abs('R4')' abs('R5')' abs('R6')']);     

%% ORIGINAL MESH DEFINEMENT
%mesh matrices 
[p_inv,e_inv,t_inv] = initmesh(d1_inv,'hmax',l/20); 
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
%         n1 = e(1,id);
%         n2 = e(2,id);
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

radius = 0.5*l; %radius of circle in which receivers are uniformly distributed
th_rec = 180-(0:(Nsources-1))*(360/Nsources); %angles of receivers location fri
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
    E_mes(is,:) = z*Ezs0(nodes,:);%+max(abs(Ezs0(:,iss)))*ita*(sqrt(2)^(-1))*(unifrnd(-1,1)+unifrnd(-1,1)*1i); %computing measured fields with additive noise uniformly distributed between -1 1
end

%% MS MD OPERATORS CONSTRUCTIONS 
%inverse problem mesh will be used

%---------------------MS matrx-------------------------%
MS = spalloc(Nsources,Nnodes_inv,3*Nsources);

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
        
    MS(is,nodes_inv) = z_inv; %assiging simplex value of element nodes that receiver is belongs
end

%---------------------MD matrx------------------------%
element_belong = find(regions_inv==5); %elements belonging in imaging domain D
element_belong1 = t_inv(1:3,element_belong); %select nodes of each triangle in domain D multiple same nodes will be selected
jj = sort(unique(element_belong1)); %select unique nodes in imaging domain D and sort them in ascending order. Represents the colomn id of node in matrix MD
index_id_d = jj; %index matrix for correspondence of imaging domain D node to global numeration
jj = index_inv(jj)'; %This matrix acts as an index matrix new numeration for nodes inside the domain. Fi 1 correspond to 59 but 59 is the numeration of the node in the index_inv matrix
ii = 1:size(jj,2); %Represents the rows of matrix

Nnodes_d = size(jj,2); %number of nodes in inversion domain d

MD = spalloc(Nnodes_d,Nnodes_inv,Nnodes_d); %matrix MD definent where size(jj,2) equals to number of nodes in domain D

linearid = sub2ind(size(MD),ii,jj); %matrix indexing for iith jjth element of MD

MD(linearid) = 1;

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
        Ez_inc(in,1:Nsources) = E0.*exp(-1i*kb.*(p_inv(1,in).*coord(1,1:Nsources).*sind(theta)+p_inv(2,in).*coord(2,1:Nsources).*cosd(theta))); 
    else
        Ez_inc(in) = 0;
    end
end

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

x_0 = U_0\V_0; %initial contrast guess

for is=1:Nsources
    r_0(:,is) = x_0.*E_t(:,is)-w_0(:,is)+x_0.*(MD*L*w_0(:,is)); %initial domain error
    gt_0(:,is) = -2*ita_s*Td\(L'*MS'*E_mes(:,is));
end
%% ----------------Iterative procedure------------------%
Noiterations = 1024; %number of iterations
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
%a_step = zeros(Noiterations,Nsources); %conjugated gradient step size
a_stepm=zeros(Nsources,1);
for ir=1:10%Noiterations %number of iterations
    %------------------------Ita_D norm factor--------------------------%   
    ita_d = 0;
    for is=1:Nsources
        ita_d = real((x_pre.*E_t(:,is))'*Td*(x_pre.*E_t(:,is)))+ita_d;
    end
    ita_d = ita_d^(-1); %D norm factor
    %-------------Polak Ribiere,gradient,contrast update----------------%       
    for is=1:Nsources
        gt_cur(:,is) = -2*ita_s*inverse1*ro_pre(:,is)-2*ita_d*r_pre(:,is)+2*ita_d*inverse2*(diag(x_pre)')*Td*r_pre(:,is);%\((eye(Nnodes_d,Nnodes_d)-L'*MD'*diag(x_pre)')*Td*r_pre(:,is)); %cost functional gradient computation
        dt_cur(:,is) = -gt_cur(:,is) + (((gt_cur(:,is)-gt_pre(:,is))'*Td*gt_cur(:,is))/(real(gt_pre(:,is)'*Td*gt_pre(:,is))))*dt_pre(:,is); %Polak Riberie search direction computation
%         dt_cur(:,is) = gt_cur(:,is) + ((real((gt_cur(:,is)-gt_pre(:,is))'*Td*gt_cur(:,is)))/(real(gt_pre(:,is)'*Td*gt_pre(:,is))))*dt_pre(:,is); %Polak Riberie search direction computation
%         a_step = (ita_s*((MS*L*dt_cur(:,is))'*ro_pre(:,is))+ita_d*((dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is)))'*Td*r_pre(:,is)))/(ita_s*(real((MS*L*dt_cur(:,is))'*(MS*L*dt_cur(:,is))))+ita_d*(real((dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is)))'*Td*(dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is))))));
        a_step_nom1 = ita_s*((MS*L*dt_cur(:,is))'*ro_pre(:,is)); %conjugated gradient step calculation breaking in parts numerator denumerator
        a_step_nom2 = ita_d*(((dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is)))')*Td*r_pre(:,is));
        a_step_denom1 = ita_s*real((MS*L*dt_cur(:,is))'*(MS*L*dt_cur(:,is)));
        a_step_denom2 = ita_d*real((dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is)))'*Td*(dt_cur(:,is)-x_pre.*(MD*L*dt_cur(:,is))));
        a_step1 = (a_step_nom1+a_step_nom2)/(a_step_denom1+a_step_denom2);
        a_stepm(is)=a_step1;
%         a_step = real(a_step);
%         a_step1 = real(a_step1);
        w_cur(:,is) = w_pre(:,is)+a_step1*dt_cur(:,is); %contrast source update
    end
    %-----------------------Contrast Update-----------------------------%
    U = 0;
    V = 0;
    for is=1:Nsources
        U = (sparse(diag(E_t(:,is)+MD*L*w_cur(:,is))))'*Td*(sparse(diag(E_t(:,is)+MD*L*w_cur(:,is))))+U;%left matrix for contrast computation U_0*x=V_0
        V = (sparse(diag(E_t(:,is)+MD*L*w_cur(:,is))))'*Td*w_cur(:,is)+V;%right matrix for contrast computation U_0*x=V_0
        ro_cur(:,is) = E_mes(:,is)-MS*L*w_cur(:,is); %initial data error
    end
    x_cur = U\V; %contrast update
    
    %---------------------Domain Error Upadate--------------------------%
    for is=1:Nsources
        r_cur(:,is) = x_cur.*E_t(:,is)-w_cur(:,is)+x_cur.*(MD*L*w_cur(:,is)); %initial domain error
    end
    
    %------------------Previous values annotation-----------------------%
    dt_pre = dt_cur; %Polak Ribiere annotation
    w_pre = w_cur; %contrast source annotation
    r_pre = r_cur; %domain error annotation
    ro_pre = ro_cur; %data erro annotation
    x_pre = x_cur;
    save umlaut.mat
    count=count+1
end
%--------------Contrast values annotation to inv mesh-------------------%
contrast = zeros(Nnodes_inv,1);
contrast(index_id_d) = x_cur; 
er_reconstructed = er_b*contrast+er_b;
save umlaut.mat
%%
toc 