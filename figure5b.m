
% 时间已过 21.289585 秒。
clear
clc %%%%%%%
close all
% Elapsed time is 9.253353 seconds.
nt=1500;    % number of time steps
eps=.6;     % stability
isnap=60;    % snapshot sampling

nx=250;
nz=250;

v=ones(nz,nx)*3500;
v(1:nz/2,:)=2800;


p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.002; % calculate time step from stability criterion
tau=dt;
r2=v.*v.*dt*dt/dx/dx;

f0=25*pi;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian



seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
zs=60;
xs=nz/2;

h=dx;
r=v*dt/h;

taper=ones(nz,nx);
for i=1:50
    for j=1:nx
        taper(i,j)=0.5-0.5*cos(pi*(i-1)/(50-1));
        taper(nz-i+1,j)=taper(i,j);
    end
end
for i=1:nz
    for j=1:50
        taper(i,j)=taper(i,j)*(0.5-0.5*cos(pi*(j-1)/(50-1)));
        taper(i,nx-j+1)=taper(i,j);
    end
end

for loop=1:1
    tic
    
    p=zeros([nz nx]); Vx=p; Vz=p;
    
    
    coeff=[ 1.25477, -0.125769, 0.0377009, -0.0143843, 0.00566346, -0.00196182, 0.00041549];
    
    coeff1= [1.83806, -0.355029, 0.100949, -0.0291963, 0.00730575, -0.00135064, 0.000135486];
    temp1=-3.12175;
    for it=1:nt-2,
        
        d1pz11=Vz-circshift(Vz,[0 1]);
        d1pz12=(circshift(Vz,[0 -1])-circshift(Vz,[0 2]));
        d1pz13=(circshift(Vz,[0 -2])-circshift(Vz,[0 3]));
        d1pz14=(circshift(Vz,[0 -3])-circshift(Vz,[0 4]));
        d1pz15=(circshift(Vz,[0 -4])-circshift(Vz,[0 5]));
        d1pz16=(circshift(Vz,[0 -5])-circshift(Vz,[0 6]));
        d1pz17=(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
        
        
        d1px11=Vx-circshift(Vx,[1 0]);
        d1px12=(circshift(Vx,[-1 0])-circshift(Vx,[2 0]));
        d1px13=(circshift(Vx,[-2 0])-circshift(Vx,[3 0]));
        d1px14=(circshift(Vx,[-3 0])-circshift(Vx,[4 0]));
        d1px15=(circshift(Vx,[-4 0])-circshift(Vx,[5 0]));
        d1px16=(circshift(Vx,[-5 0])-circshift(Vx,[6 0]));
        d1px17=(circshift(Vx,[-6 0])-circshift(Vx,[7 0]));
        
        
        d1px=coeff(1)*d1px11+coeff(2)*d1px12+coeff(3)*d1px13+coeff(4)*d1px14+coeff(5)*d1px15+coeff(6)*d1px16...
            +coeff(7)*d1px17;
        d1pz=coeff(1)*d1pz11+coeff(2)*d1pz12+coeff(3)*d1pz13+coeff(4)*d1pz14+coeff(5)*d1pz15+coeff(6)*d1pz16...
            +coeff(7)*d1pz17;
        
        Q1=-v.^2.*(d1px+d1pz);
        
        
        d2px11=(circshift(Q1,[0 -1])+circshift(Q1,[0 1]));
        d2px12=(circshift(Q1,[0 -2])+circshift(Q1,[0 2]));
        d2px13=(circshift(Q1,[0 -3])+circshift(Q1,[0 3]));
        d2px14=(circshift(Q1,[0 -4])+circshift(Q1,[0 4]));
        d2px15=(circshift(Q1,[0 -5])+circshift(Q1,[0 5]));
        d2px16=(circshift(Q1,[0 -6])+circshift(Q1,[0 6]));
        d2px17=(circshift(Q1,[0 -7])+circshift(Q1,[0 7]));
        
        d2pz11=(circshift(Q1,[-1 0])+circshift(Q1,[1 0]));
        d2pz12=(circshift(Q1,[-2 0])+circshift(Q1,[2 0]));
        d2pz13=(circshift(Q1,[-3 0])+circshift(Q1,[3 0]));
        d2pz14=(circshift(Q1,[-4 0])+circshift(Q1,[4 0]));
        d2pz15=(circshift(Q1,[-5 0])+circshift(Q1,[5 0]));
        d2pz16=(circshift(Q1,[-6 0])+circshift(Q1,[6 0]));
        d2pz17=(circshift(Q1,[-7 0])+circshift(Q1,[7 0]));
        
        d2px=coeff1(1)*d2px11+coeff1(2)*d2px12+coeff1(3)*d2px13+coeff1(4)*d2px14+coeff1(5)*d2px15+coeff1(6)*d2px16...
            +coeff1(7)*d2px17;
        d2pz=coeff1(1)*d2pz11+coeff1(2)*d2pz12+coeff1(3)*d2pz13+coeff1(4)*d2pz14+coeff1(5)*d2pz15+coeff1(6)*d2pz16...
            +coeff1(7)*d2pz17;
        
        d2pxz=v.^2.*(d2px+d2pz+2*temp1*Q1);
        
        
        p=p-dt*v.^2.*(d1px+d1pz)/h +d2pxz*dt^3/h^3/24;
        
        %           p=p-dt*v.^2.*(d1px+d1pz)/h ;
        p(zs,xs)= p(zs,xs)+src(it)*dt^2;
        
        
        d1pz1=(circshift(p,[0 -1])-circshift(p,[0 0]));
        d1pz2=(circshift(p,[0 -2])-circshift(p,[0 1]));
        d1pz3=(circshift(p,[0 -3])-circshift(p,[0 2]));
        d1pz4=(circshift(p,[0 -4])-circshift(p,[0 3]));
        d1pz5=(circshift(p,[0 -5])-circshift(p,[0 4]));
        d1pz6=(circshift(p,[0 -6])-circshift(p,[0 5]));
        d1pz7=(circshift(p,[0 -7])-circshift(p,[0 6]));
        
        
        d1px1=(circshift(p,[-1])-circshift(p,[0]));
        d1px2=(circshift(p,[-2])-circshift(p,[1]));
        d1px3=(circshift(p,[-3])-circshift(p,[2]));
        d1px4=(circshift(p,[-4])-circshift(p,[3]));
        d1px5=(circshift(p,[-5])-circshift(p,[4]));
        d1px6=(circshift(p,[-6])-circshift(p,[5]));
        d1px7=(circshift(p,[-7])-circshift(p,[6]));
        
        d1px=coeff(1)*d1px1+coeff(2)*d1px2+coeff(3)*d1px3+coeff(4)*d1px4+coeff(5)*d1px5+coeff(6)*d1px6...
            +coeff(7)*d1px7;
        d1pz=coeff(1)*d1pz1+coeff(2)*d1pz2+coeff(3)*d1pz3+coeff(4)*d1pz4+coeff(5)*d1pz5+coeff(6)*d1pz6...
            +coeff(7)*d1pz7;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        d2px11=(circshift(d1px,[-1 0])+circshift(d1px,[1 0]));
        d2px12=(circshift(d1px,[-2 0])+circshift(d1px,[2 0]));
        d2px13=(circshift(d1px,[-3 0])+circshift(d1px,[3 0]));
        d2px14=(circshift(d1px,[-4 0])+circshift(d1px,[4 0]));
        d2px15=(circshift(d1px,[-5 0])+circshift(d1px,[5 0]));
        d2px16=(circshift(d1px,[-6 0])+circshift(d1px,[6 0]));
        d2px17=(circshift(d1px,[-7 0])+circshift(d1px,[7 0]));
        
        d2px=temp1*d1px+coeff1(1)*d2px11+coeff1(2)*d2px12+coeff1(3)*d2px13+coeff1(4)*d2px14+coeff1(5)*d2px15+coeff1(6)*d2px16...
            +coeff1(7)*d2px17;
        
        d2pz11=(circshift(d1px,[0 -1])+circshift(d1px,[0 1]));
        d2pz12=(circshift(d1px,[0 -2])+circshift(d1px,[0 2]));
        d2pz13=(circshift(d1px,[0 -3])+circshift(d1px,[0 3]));
        d2pz14=(circshift(d1px,[0 -4])+circshift(d1px,[0 4]));
        d2pz15=(circshift(d1px,[0 -5])+circshift(d1px,[0 5]));
        d2pz16=(circshift(d1px,[0 -6])+circshift(d1px,[0 6]));
        d2pz17=(circshift(d1px,[0 -7])+circshift(d1px,[0 7]));
        
        d2pz=temp1*d1px+coeff1(1)*d2pz11+coeff1(2)*d2pz12+coeff1(3)*d2pz13+coeff1(4)*d2pz14+coeff1(5)*d2pz15+coeff1(6)*d2pz16...
            +coeff1(7)*d2pz17;
        
        Vx=Vx-dt*d1px/h -v.^2.*(d2px+d2pz)*dt^3/24/h^3;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        d2pz11=(circshift(d1pz,[0 -1])+circshift(d1pz,[0 1]));
        d2pz12=(circshift(d1pz,[0 -2])+circshift(d1pz,[0 2]));
        d2pz13=(circshift(d1pz,[0 -3])+circshift(d1pz,[0 3]));
        d2pz14=(circshift(d1pz,[0 -4])+circshift(d1pz,[0 4]));
        d2pz15=(circshift(d1pz,[0 -5])+circshift(d1pz,[0 5]));
        d2pz16=(circshift(d1pz,[0 -6])+circshift(d1pz,[0 6]));
        d2pz17=(circshift(d1pz,[0 -7])+circshift(d1pz,[0 7]));
        
        d2pz=temp1*d1pz+coeff1(1)*d2pz11+coeff1(2)*d2pz12+coeff1(3)*d2pz13+coeff1(4)*d2pz14+coeff1(5)*d2pz15+coeff1(6)*d2pz16...
            +coeff1(7)*d2pz17;
        
        
        d2px11=(circshift(d1pz,[-1 0])+circshift(d1pz,[1 0]));
        d2px12=(circshift(d1pz,[-2 0])+circshift(d1pz,[2 0]));
        d2px13=(circshift(d1pz,[-3 0])+circshift(d1pz,[3 0]));
        d2px14=(circshift(d1pz,[-4 0])+circshift(d1pz,[4 0]));
        d2px15=(circshift(d1pz,[-5 0])+circshift(d1pz,[5 0]));
        d2px16=(circshift(d1pz,[-6 0])+circshift(d1pz,[6 0]));
        d2px17=(circshift(d1pz,[-7 0])+circshift(d1pz,[7 0]));
        
        d2px=temp1*d1pz+coeff1(1)*d2px11+coeff1(2)*d2px12+coeff1(3)*d2px13+coeff1(4)*d2px14+coeff1(5)*d2px15+coeff1(6)*d2px16...
            +coeff1(7)*d2px17;
        
        Vz=Vz-dt*d1pz/h -v.^2.*(d2pz+d2px)*dt^3/24/h^3;
        %%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        
        
        
    
        
        [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
        %         Vx=taper.*Vx;
        %         Vz=Vz.*taper;
        
        if rem(it,isnap)== 0,
            imagesc(x,z,p), axis equal
            colormap gray
            xlabel('x'),ylabel('z')
            title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
            drawnow
        end
        
        
        seis_record(it,:)=p(60,:);
        
        
    end
    
end
toc
save('figure5b.mat')
% figure;plot((a),'linewidth',2);hold on;plot((b),'r','linewidth',2)
% xlabel('time(ms)')
% ylabel('Amp')
% legend('receiver A','receiver B')
% grid on
