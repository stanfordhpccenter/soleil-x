%% DOM3D.m
% a script to solve the RTE in a 3D rectangular cavity using the discrete
% ordinates method
%% data initialization

%load json params
dom_param_filename = 'tgv_64x64x64_dom';
fid = fopen(strcat('../testcases/',dom_param_filename,'.json')); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
params = jsondecode(str);

%grid size
Nx = params.Radiation.xNum;
Ny = params.Radiation.yNum;
Nz = params.Radiation.zNum;

%load quadrature 
angle_filename = strcat('./LMquads/',num2str(params.Radiation.angles),'.txt');
quads = textread(angle_filename);  %#ok<DTXTRD>
N = quads(1);
xi = quads(2:N+1);
eta = quads(N+2:2*N+1);
mu = quads(2*N+2:3*N+1);
w = quads(3*N+2:end);

%domain size
Lx = params.Grid.xWidth;
Ly = params.Grid.yWidth;
Lz = params.Grid.zWidth;

%grid spacing
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;
dAx = dy*dz;
dAy = dx*dz;
dAz = dx*dy;
dV = dx*dy*dz;

%extinction coefficient
sigma = 5*ones(Nx,Ny,Nz);

%albedo
omega = params.Radiation.qs/(params.Radiation.qa+params.Radiation.qs);

%wall emissivity
epsw = [params.Radiation.emissWest params.Radiation.emissEast params.Radiation.emissSouth params.Radiation.emissNorth params.Radiation.emissDown params.Radiation.emissUp];

%wall temperature
SB = 5.67e-8;
Tw = [params.Radiation.tempWest params.Radiation.tempEast params.Radiation.tempSouth params.Radiation.tempNorth params.Radiation.tempDown params.Radiation.tempUp];

%blackbody source
T = 1000*ones(Nx,Ny,Nz);
Ib = SB/pi*T.^4;

%load the quadrature and create the intensity
I = zeros(N,Nx,Ny,Nz); %cell centered intensity
Ifx = zeros(N,Nx+1,Ny,Nz); %x-face centered intensity
Ify = zeros(N,Nx,Ny+1,Nz); %y-face centered intensity
Ifz = zeros(N,Nx,Ny,Nz+1);
Iiter = I; %iterate for the intensity
S = zeros(Nx,Ny,Nz); %source term

%% solution procedure

tol = 1e-6; %solution tolerance
res = 1; %initial residual
gamma = 0.5; %1 for step differencing, 0.5 for diamond differencing
residual = 1; %residual history
igamma = 1/gamma;

while (res > tol)
  
    %update the source term (in this problem, isotropic)
    S = (1-omega)*sigma.*Ib;
    for i = 1:N
       S = S + omega*sigma./(4*pi)*w(i).*squeeze(Iiter(i,:,:,:)); 
    end
    
    %update the boundary intensities
    %west
    reflect1 = 0;
    if epsw(1) < 1
    for m = 1:N
       if xi(m) < 0
          reflect1 = reflect1 + (1-epsw(1))/pi*w(m)*abs(xi(m))*squeeze(Ifx(m,1,:,:)); 
       end
    end
    end
    for m = 1:N
       if xi(m) > 0
          Ifx(m,1,:,:) = epsw(1)*SB*Tw(1)^4/pi + reflect1;
       end
    end
    
    %east
    reflect2 = 0;
    if epsw(2) < 1
    for m = 1:N
       if xi(m) > 0
          reflect2 = reflect2 + (1-epsw(2))/pi*w(m)*xi(m)*squeeze(Ifx(m,Nx+1,:,:)); 
       end
    end
    end
    for m = 1:N
       if xi(m) < 0
          Ifx(m,Nx+1,:,:) = epsw(2)*SB*Tw(2)^4/pi + reflect2;
       end
    end
    
    % south
    reflect3 = 0;
    if epsw(3)<1
    for m = 1:N
       if eta(m) < 0
          reflect3 = reflect3 + (1-epsw(3))/pi*w(m)*abs(eta(m))*squeeze(Ify(m,:,1,:)); 
       end
    end
    end
    for m = 1:N
       if eta(m) > 0
          Ify(m,:,1,:) = epsw(3)*SB*Tw(3)^4/pi + reflect3;
       end
    end
    
    % north
    reflect4 = 0;
    if epsw(4) < 1
    for m = 1:N
       if eta(m) > 0
          reflect4 = reflect4 + (1-epsw(4))/pi*w(m)*eta(m)*squeeze(Ify(m,:,Ny+1,:)); 
       end
    end
    end
    for m = 1:N
       if eta(m) < 0
          Ify(m,:,Ny+1,:) = epsw(4)*SB*Tw(4)^4/pi + reflect4;
       end
    end
    
    % down
    reflect5 = 0;
    if epsw(5)<1
    for m = 1:N
       if mu(m) < 0
          reflect5 = reflect5 + (1-epsw(5))/pi*w(m)*abs(mu(m))*squeeze(Ifz(m,:,:,1)); 
       end
    end
    end
    for m = 1:N
       if mu(m) > 0
          Ifz(m,:,:,1) = epsw(5)*SB*Tw(5)^4/pi + reflect5;
       end
    end
    
    % up
    reflect6 = 0;
    if epsw(6)<1
    for m = 1:N
       if mu(m) > 0
          reflect6 = reflect6 + (1-epsw(6))/pi*w(m)*mu(m)*squeeze(Ifz(m,:,:,Nz+1)); 
       end
    end
    end
    for m = 1:N
       if mu(m) < 0
          Ifz(m,:,:,Nz+1) = epsw(6)*SB*Tw(6)^4/pi + reflect6;
       end
    end
    
    %sweep for the new intensity
    for m = 1:N
        %determine sweeping direction
        %if xi>0, angle points in +x, sweep from left to right
        if (xi(m) > 0)
            dindx = 1;
            startx = 1;
            endx = Nx;
        else %angle points in -x, sweep from right to left
            dindx = -1;
            startx = Nx;
            endx = 1;
        end
        if (eta(m) > 0) %angle points in +y, sweep from bottom to top
            dindy = 1;
            starty = 1;
            endy = Ny;
        else %angle points in -y, sweep from top to bottom
            dindy = -1;
            starty = Ny;
            endy = 1;
        end
        if (mu(m) > 0) %angle points in +z, sweep from back to front
            dindz = 1;
            startz = 1;
            endz = Nz;
        else %angle points in -z, sweep from front to back
            dindz = -1;
            startz = Nz;
            endz = 1;
        end
        %some prefactors to reduce flop count in sweep
        axi = abs(xi(m));
        aeta = abs(eta(m));
        amu = abs(mu(m));
        
        mdx = min(dindx,0);
        mdy = min(dindy,0);
        mdz = min(dindz,0);
        for k = startz:dindz:endz %from starting z point to finishing z
        
            for j = starty:dindy:endy %same for y
           
                for i = startx:dindx:endx %same for x
                
                    indx = i - mdx; %index of upwind x-face intensity
                    indy = j - mdy; %index of upwind y-face intensity
                    indz = k - mdz; %index of upwind z-face intensity
                    
                    %FYI: on a stretched grid, use dAx, dAy, dAz, dV based 
                    %on local cell sizes, eg. dAx(i,j,k), dV(i,j,k)...
                    
                    %compute cell-centered I
                    I(m,i,j,k) = (S(i,j,k)*dV + Ifx(m,indx,j,k)*axi*dAx*igamma  ...
                    + Ify(m,i,indy,k)*aeta*dAy*igamma + Ifz(m,i,j,indz)*amu*dAz*igamma)...
                    /(sigma(i,j,k)*dV + axi*dAx*igamma + aeta*dAy*igamma + amu*dAz*igamma); 
                    
                    %compute downwind intensities
                    Ifx(m,indx+dindx,j,k) = max((I(m,i,j,k)-(1-gamma)*Ifx(m,indx,j,k))*igamma,0);
                    Ify(m,i,indy+dindy,k) = max((I(m,i,j,k)-(1-gamma)*Ify(m,i,indy,k))*igamma,0);
                    Ifz(m,i,j,indz+dindz) = max((I(m,i,j,k)-(1-gamma)*Ifz(m,i,j,indz))*igamma,0);
                    
                end
                
            end
            
        end
        
    end
    
    %compute a residual for intensity using 2-norm
    res = sqrt(1/(Nx*Ny*Nz*N)*sum(sum(sum(sum(((I-Iiter).^2)./I.^2)))))
%     residual = horzcat(residual,res);
    
    %update the iterate. I is the new guess for the intensity
    Iiter = I;
             
end


%reduce the intensity to summation over all angles
G = zeros(Nx,Ny,Nz);
for i = 1:N
   G = G+w(i)*squeeze(I(i,:,:,:)); 
end

output_filename = strcat('dom_intensities/intensity-matlab-', dom_param_filename, '.dat');
f = fopen(output_filename, "w");
for i = 1:Nx 
    for j = 1:Ny 
        for k = 1:Nz
            fprintf(f,' %.15e \n', G(i,j,k));
        end
        fprintf(f,'\n');
    end
	fprintf(f,'\n');
end
fclose(f);

% figure
% subplot(1,2,1)
% imagesc(squeeze(G(floor(Nx/2),:,:))')
% axis xy
% 
% subplot(1,2,2)
% semilogy(residual,'k','LineWidth',2)
% title(strcat('spectral radius: ',num2str(residual(end)/residual(end-1))),'FontSize',20)
% 

