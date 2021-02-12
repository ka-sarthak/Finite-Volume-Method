% Author: Sarthak Kapoor; RWTH Aachen; February 2021
% solving the 1D Euler equations using FVM code
% HLL flux condition
% CFL condition

clear; clf;

N=200;      %number of grid cells
T=0.25;     %numerical time
CFL=0.5;    %used as a relation between dt and dx --- represents CFL-factor
a=0; b=1;  %start and end of the space domain x (1D)
gamma=1.4;  % ----

u0=cell(3,1);           %define 3 rows of u0 for the initial condition

%using initial condition 1
%u0{1}= @(x)     1.*(x>=0 & x<=0.25) +   0.125.*(x>0.25 & x<=1);
%u0{2}= @(x)     u0{1}(x).*(0.*(x>=0 & x<=0.25) +   0.*(x>0.25 & x<=1));
%u0{3}= @(x)     (u0{2}(x).^2)./(2*u0{1}(x))  +   (1/(gamma-1)).*(1*(x>=0 & x<=0.25) +   0.1.*(x>0.25 & x<=1));

%using initial condition 2
u0{1}= @(x)     1 .* and(x>=0, x<=0.25)   +     1 .* and(x>0.25, x<=1);
u0{2}= @(x)     u0{1}(x) .* (-1 .* and(x>=0, x<=0.25)  +     1 .* and(x>0.25, x<=1));
u0{3}= @(x)     (u0{2}(x).^2)/(2.*u0{1}(x))  +   (1/(gamma-1))*(0.4 .* and(x>=0, x<=0.25) +     0.4 .* and(x>0.25, x<=1));

u=solve(N,a,b,T,CFL,u0,gamma);


%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%
function u=solve(N, a, b, T, CFL, u0, gamma)
    dx=(b-a)/N;
    u=zeros(3,N,1);         %3 dimensions, n grid points, at time step number 1
    
    for i=1:N               %giving the cell mean values
        u(1,i,1)= (1/dx)*integral(u0{1}, a+(i-1)*dx, a+i*dx);
        u(2,i,1)= (1/dx)*integral(u0{2}, a+(i-1)*dx, a+i*dx);
        u(3,i,1)= (1/dx)*integral(u0{3}, a+(i-1)*dx, a+i*dx);
    end
    
    t=0; j=1;       %j is the time step number and t is the time completed
    
    while (true)
        Aj = A(u(:,:,j), gamma);
        
        dt = CFL*dx /(max(abs(eig(Aj))));
        
        lambda= dt/dx;
        
        if t+dt > T         %loop termination
            break; 
        end
        
        u(:,1,j+1) = u(:,1,j) - (dt/dx) .* (    NumF(u(:,1,j),u(:,2,j),lambda,gamma)-NumF(u(:,1,j),u(:,1,j),lambda,gamma)  ) ;
        for i=2:N-1
            u(:,i,j+1) = u(:,i,j) - (dt/dx) .* (    NumF(u(:,i,j),u(:,i+1,j),lambda,gamma)-NumF(u(:,i-1,j),u(:,i,j),lambda,gamma)  ); 
        end    
        u(:,N,j+1) = u(:,N,j) - (dt/dx) .* (    NumF(u(:,N,j),u(:,N,j),lambda,gamma)-NumF(u(:,N-1,j),u(:,N,j),lambda,gamma)  ) ;
       
        plotDynamic(u(:,:,j),t,N,a,b,dt);
        t=t+dt;
        j=j+1;
    end
    
        %since the loop terminated dt before T, one cannot get evolution till T, rather one time step behind
        %therefore, we run the evolution for one more time with dt=T-t
        dt = T-t;
        
        lambda= dt/dx;
        
        u(:,1,j+1) = u(:,1,j) - (dt/dx) .* (    NumF(u(:,1,j),u(:,2,j),lambda,gamma)-NumF(u(:,1,j),u(:,1,j),lambda,gamma)  ) ;
        for i=2:N-1
            u(:,i,j+1) = u(:,i,j) - (dt/dx) .* (    NumF(u(:,i,j),u(:,i+1,j),lambda,gamma)-NumF(u(:,i-1,j),u(:,i,j),lambda,gamma)  ); 
        end    
        u(:,N,j+1) = u(:,N,j) - (dt/dx) .* (    NumF(u(:,N,j),u(:,N,j),lambda,gamma)-NumF(u(:,N-1,j),u(:,N,j),lambda,gamma)  ) ;
       
        t=t+dt;
        plotDynamic(u(:,:,j),t,N,a,b,dt);       %final plot at time T
        
end

function Aj= A(u, gamma)
    u1=u(1);
    u2=u(2);
    u3=u(3);
        
    Aj(1,:)=[0                                               1                                           0];
    Aj(2,:)=[0.5*(gamma-3)*((u2/u1)^2)                       (3-gamma)*(u2/u1)                           gamma-1];
    Aj(3,:)=[-(gamma*u2*u3)/(u1^2) + (gamma-1)*((u2/u1)^3)   (gamma*u3)/(u1)-1.5*(gamma-1)*((u2/u1)^2)   (u2/u1)*gamma];
end

function F= NumF(u,v,lambda, gamma)
    
    al= min(min(eig(A(u, gamma)), eig(A(v, gamma))));
    ar= max(max(abs(eig(A(u, gamma))),abs(eig(A(v, gamma)))));
        
    if (al>0)
        F = f(u, gamma); 
    elseif (al<=0 && ar>=0)
        F = (ar.*f(u,gamma)-al.*f(v,gamma)+al*ar.*(v-u))./(ar-al); 
    elseif (ar<0)
        F = f(v, gamma);
    else
        F = [1997 1997 1997];
    end
    
    %F= (A(u,gamma)*u + A(v, gamma)*v - (v-u)/lambda)/2;    %Lax Friedrisch
        
end

function flux=f(u, gamma)
    u1=u(1);
    u2=u(2);
    u3=u(3);
    
    flux(1)= u2;
    flux(2)= u3.*(gamma-1) + ((3-gamma)/2).*(u2.^2)./u1;
    flux(3)= (gamma.*u2.*u3./u1) - ((gamma-1)/2)*(u2.^3)./(u1.^2);
    
    flux=flux';
    
end

function plotDynamic(u,t,N,a,b,dt)
    grid=linspace(a+(1/2)/N,b-(1/2)/N,N);
    plot(grid,u(1,:));
    hold on;
    plot(grid,u(2,:));
    plot(grid,u(3,:));
    hold off;
    legend('u1=rho','u2=rho.v','u3=E');
    title(['Time: ' num2str(t)]);
    pause(dt);
end


