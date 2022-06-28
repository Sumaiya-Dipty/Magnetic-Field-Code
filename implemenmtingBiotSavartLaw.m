%% http://biomedicalsignalandimage.blogspot.com/2016/08/matlab-code-to-calculate-magnetic-flux.html

clear all
close all
clc

global u0
u0=4*pi*1e-7;

d = 2;                  % Radius of the Loop.. (2cm = 20 mm)
m = [-1, -1, -1];                 % Direction of current through the Loop... 1(ANTI-CLOCK), -1(CLOCK)
N = 30;                 % Number sections/elements in the conductor...
T = 1;                 % Number of Turns in the Solenoid..
g = 0.75;               % Gap between Turns.. (mm)
dl = 2*pi*d/N;          % Length of each element..

n_coil = 3;
coil_center_x = [0, 2.5*d, 5*d];
coil_center_y = zeros(1, 3);


% XYZ Coordinates/Location of each element from the origin(0,0,0), Center of the Loop is taken as origin..
dtht = 360/N;                           % Loop elements divided w.r.to degrees..
tht = (0+dtht/2): dtht : (360-dtht/2);  % Angle of each element w.r.to origin..

xC =  d.*cosd(tht).*ones(1,N);      % X coordinate...
yC =  d.*sind(tht).*ones(1,N);      % Y coordinate...

x_mid = d*cosd(0:dtht:dtht*(N-1)); % x-cordinate of the mid point of each element
y_mid = d*sind(0:dtht:dtht*(N-1)); % x-cordinate of the mid point of each element

xC = repmat(xC, n_coil, 1);
yC = repmat(yC, n_coil, 1);

x_mid = repmat(x_mid, n_coil, 1);
y_mid = repmat(y_mid, n_coil, 1);

for i=1:n_coil
    xC(i, :) = xC(i, :) + coil_center_x(i);
    yC(i, :) = yC(i, :) + coil_center_y(i);
    
    x_mid(i, :) = x_mid(i, :) + coil_center_x(i);
    y_mid(i, :) = y_mid(i, :) + coil_center_y(i);
end

% zC will change for each turn, -Z to +Z, and will be varied during iteration...
zC = ones(1,N);
h = -(g/2)*(T-1) : g : (g/2)*(T-1);

z_mid = h;

% Length(Projection) & Direction of each current element in Vector form..
Lx = dl.*cosd(90+tht);        % Length of each element on X axis..
Ly = dl.*sind(90+tht);        % Length of each element on Y axis..
Lz = zeros(1,N);              % Length of each element is zero on Z axis..

% Points/Locations in space (here XZ plane) where B is to be computed..
NP = 50;               % Detector points..
xPmax = 7*d;            % Dimensions of detector space.., arbitrary..
zPmax = 1.5*(T*g);

% xP = linspace(-1.5*d,xPmax,NP);        % Divide space with NP points..
% yP = linspace(-1.5*d,1.5*d,NP);        % Divide space with NP points..
% zP = linspace(-zPmax,zPmax,NP);

xP = linspace(-2*d, xPmax, NP);
yP = xP;
zP = yP;
[xxP, yyP, zzP] = meshgrid(xP,yP,zP);            % Creating the Mesh..

% Initialize B..
% Bx = zeros(NP,NP, NP);
% By = zeros(NP,NP, NP);
% Bz = zeros(NP,NP, NP);

B_vect = cell(NP,NP,NP);
for i=1:NP
    for j=1:NP
        for k=1:NP
            B_vect{i,j,k} = [0 0 0];
        end
    end
end

% Computation of Magnetic Field (B) using Superposition principle..
% Compute B at each detector points due to each small cond elements & integrate them..
tic
for n_ = 1:n_coil
    for p = 1:T
        for q = 1:N
            leng_tang = [y_mid(n_, q) -x_mid(n_, q) 0]; % tangent vector to the current element
            leng_tang = dl*leng_tang/norm(leng_tang);
            for q1 = 1:NP
                for q2 = 1:NP
                    for q3 = 1:NP 
                        rx = xxP(q1,q2,q3) - x_mid(n_, q);               % Displacement Vector along X direction from Conductor..
                        ry = yyP(q1,q2,q3) - y_mid(n_, q);                     % No detector points on Y direction..
                        rz = zzP(q1,q2,q3) - h(p)*zC(q);          % zC Vertical location changes per Turn..

                        r = sqrt(rx.^2+ry.^2+rz.^2);    % Displacement Magnitude for an element on the conductor..

                        r3 = r.^3;
                           
                        B_vect{q1,q2,q3} = B_vect{q1,q2,q3} + m(n_)*cross(leng_tang, [rx ry rz])/r3;

%                         Bx = Bx + m(n_)*Ly(q).*rz./r3;      % m - direction of current element..
%                         By = By - m(n_)*Lx(q).*rz./r3;
%                         Bz = Bz + m(n_)*Lx(q).*ry./r3 - m(n_)*Ly(q).*rx./r3;
                    end
                end
            end
        end
    end
end
toc
% B = sqrt(Bx.^2 + By.^2 + Bz.^2);        % Magnitude of B..
B = zeros(NP,NP,NP);
for i=1:NP
    for j=1:NP
        for k=1:NP
            B(i,j,k) = norm(B_vect{i,j,k});
        end
    end
end
B = B/max(max(B(:)));                      % Normalizing...

for i=1:NP
% Plotting...
 figure(1);
 pcolor(xxP(:,:,1),yyP(:,:,1),B(:,:,i));
 colormap(jet);
 shading interp;
 axis equal;
 axis([-1.5*d xPmax -2*d 2*d -zPmax zPmax]);
 xlabel('<-- x -->');ylabel('<-- y -->');
 title('Magnetic Field Distibution');
 colorbar;
 pause(0.5);
end
 
NP = 25;
 
 temp_x = zeros(n_coil, NP);
 temp_y_pos = zeros(n_coil, NP);
 temp_y_neg = zeros(n_coil, NP);
 
 for n_ = 1:n_coil
    temp_center_x = coil_center_x(n_);
    temp_center_y = coil_center_y(n_);
    
    temp_x(n_, :) = linspace(coil_center_x(n_)-d,coil_center_x(n_)+d,NP);
    temp_y_pos(n_,:) = (d^2 - (temp_x(n_,:)-coil_center_x(n_)).^2) + coil_center_y(n_);
    temp_y_neg(n_,:) = -(d^2 - (temp_x(n_,:)-coil_center_x(n_)).^2) + coil_center_y(n_);
    figure(1); hold on; plot(temp_x(n_,:), temp_y_pos(n_,:), 'x', temp_x(n_,:), temp_y_neg(n_,:), 'x');
 end
  

% figure(2);
% surf(xxP(:,:,1),yyP(:,:,1),B(:,:,1),'FaceColor','interp',...
%     'EdgeColor','none',...
%     'FaceLighting','phong');
% daspect([1 1 1]);
% axis tight;
% view(-10,30);
% camlight right;
% colormap(jet);
% grid off;
% axis off;
% colorbar;
% title('Magnetic Field Distibution');

%  figure(3);
% quiver(xxP,zzP,Bx,Bz);
% colormap(lines);
%  axis([-1.5*d 8*d -T*g T*g]);
%  title('Magnetic Field Distibution');
%  xlabel('<-- x -->');ylabel('<-- z -->');
%  zoom on;

% figure(4);
% % [faces,verts,colors] = isosurface(xxP,yyP,zzP,B,0,B);
% isosurface(xxP,yyP,zzP,B,0,B);
% % faces_m = randn(size(faces));
% % verts_m = randn(size(verts));
% % colors_m = randn(size(colors));
% %    patch('Vertices', verts, 'Faces', faces, ...
% %       'FaceVertexCData', colors, ...
% %       'FaceColor','interp', ...
% %       'edgecolor', 'interp')
%    view(3)
%    colorbar;
   
%  figure(5); 
% pointsize = 8;  %adjust at will
% scatter3(xxP(:), yyP(:), zzP(:), pointsize, B(:));



 % Computation of (H)Field  using Superposition principle..
% Compute H at each detector points due to each small cond elements & integrate them..
% % for n_=1:n_coil
% %     for p = 1:T
% %         for q = 1:N
% %             rx = xxP - xC(n_, q);               % Displacement Vector along X direction from Conductor..
% %             ry = yC(n_, q);                     % No detector points on Y direction..
% %             rz = zzP - h(p)*zC(q);          % zC Vertical location changes per Turn..
% % 
% %             r = sqrt(rx.^2+ry.^2+rz.^2);    % Displacement Magnitude for an element on the conductor..
% % 
% %             r3 = r.^3;
% % 
% % 
% %             Hx = (Bx + m(n_)*Ly(q).*rz./r3)*u0;      % m - direction of current element..
% %             Hy = (By - m(n_)*Lx(q).*rz./r3)*u0;
% %             Hz = (Bz + m(n_)*Lx(q).*ry./r3 - m(n_)*Ly(q).*rx./r3)*u0;
% %         end
% %     end
% % end
% %  H= sqrt(Hx.^2 + Hy.^2 + Hz.^2);        % Magnitude of H..
% % H = H/max(max(H(:)));                      % Normalizing...








