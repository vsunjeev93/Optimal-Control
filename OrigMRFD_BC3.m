function dudt=OrigMRFD_BC3(t,u,St1,St2,B,CP,n,Bi,Fo,Da,gamma)
% Nodal discretization N=200
z=linspace(0,1,200);
% total number of discretized state variables.
N=3*size(z,2);
% step size
deltax=z(2)-z(1);
% boundary conditions for concentration and temperature. 
%Boundary conditions for wall are applied directly in the code.
u(1)=1;
u(N/3+1)=0;
% Initializing the column vector. The number N-2 consists of all the state variables
% except boundaries of the wall which are added into the equations
% directly via boundary conditions
dudt=zeros(N-2,1);
% coolant temperature
vc=0;
for i=2:N/3
% Dimensionless concentration for the reaction stream
    dudt(i) = -(u(i)-u(i-1))/deltax -Da*exp(u(N/3+i)/(1+u(N/3+i)/gamma))*u(i)^n;
end
for i=2:N/3-1% Dimensionless temp for the reaction stream
    dudt(N/3+i) = -(u(N/3+i)-u(N/3+i-1))/deltax +B*Da*exp(u(N/3+i)/(1+u(N/3+i)/gamma))*u(i)^n-St1*(u(N/3+i)-u((2*N)/3+i-1));
end
i=N/3;
dudt(N/3+i) = -(u(N/3+i)-u(N/3+i-1))/deltax +B*Da*exp(u(N/3+i)/(1+u(N/3+i)/gamma))*u(i)^n-St1*(u(N/3+i)-u((2*N)/3+i-2)/(1+Bi*deltax));
%wall equations- boundary conditions for are used directlty in the
%equations. See i=2 and i=N/3-1
i=2;
dudt((2*N)/3+i-1) =Fo*((u((2*N)/3+i)+u((2*N)/3+i-1)/(1+Bi*deltax)-2*u((2*N)/3+i-1))/deltax^2)+ Fo*(St1/CP*(u(N/3+i)-u((2*N)/3+i-1))+St2/CP*(vc-u(2*N/3+i-1)));
for i=3:N/3-2
    dudt((2*N)/3+i-1) = Fo*((u((2*N)/3+i)+u((2*N)/3+i-2)-2*u((2*N)/3+i-1))/deltax^2)+ Fo*(St1/CP*(u(N/3+i)-u((2*N)/3+i-1))+St2/CP*(vc-u(2*N/3+i-1)));
end
i=N/3-1;
dudt((2*N)/3+i-1)=Fo*((u((2*N)/3+i-2)+u((2*N)/3+i-1)/(1+Bi*deltax)-2*u((2*N)/3+i-1))/deltax^2)+ Fo*(St1/CP*(u(N/3+i)-u((2*N)/3+i-1))+St2/CP*(vc-u(2*N/3+i-1)));

end


    
