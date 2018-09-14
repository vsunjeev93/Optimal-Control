function dudt=OrigMRFD_BC3(t,u,St1,St2,B,CP,n,Bi,Fo)
%%dimension of the reactor in meters
% L=0.1;W= 0.3*10^-3;D=0.2*10^-3;
% %%Reaction parameters (first order in methane)
%  k0=4*10^8; Ea=9*10^4, del_H=-882*10^3; % J per mol
% % inlet conditions in K, Pa and mole fraction respectively.
%  Tf=700; Ptot=1.1*10^5; mf_CH4=0.1; R=8.314, v=0.1;% velocity (in m/s)
%  CAf=0.1*Ptot/R/Tf; % moles per m^3
% % molar density of mixture
%  rho=Ptot/R/Tf;
% % average heat capacity 
%  cp_CH4=62.927;% J/mol/K
%  cp_air=31.8710;
%  cp=0.1*cp_CH4+0.9*cp_air;
% %residence time
%  theta=0.001;
% % Heat transfer coeff
%  h=11.3;
% %surface area per unit volume
%  a_hat=2*(L*W+L*D)/(L*W*D);
% % wall parameters
% Kw=16; rho_w=8000;cp_w=500;
% % % Dimensionless numbers
% Da_star=k0*exp(-Ea/(R*Tf)); % damkohler_star
% T_ad=(CAf*(-del_H))/(rho*cp); % adiabatic temperature 
% gamma=(T_ad)/Tf; % adiabatic temperature rise
% alpha=Ea/R/Tf; %Dimensionless activation energy
% NTU_star1=h*a_hat/(rho*cp);
% psi=1.52*10^-1;% heat removed
% Fo_star=Kw/L^2/(rho_w*cp_w); % Fourier number star
% CP_star=Kw/L^2/(rho*cp);% Conduction Parameter star
% %%%%%%%%%%%%%%%%%%%%%%%%%%%------------------------------------------------------- code starts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Test values
% Da_star=7.1.2463*10^-15;T_ad=2.5217*10^3;gamma=12.6084;alpha=54.1256;NTU_star1=81.3949;
% psi=-1*10^7;
% Fo_star=4*10^-4;CP_star=0.6915; theta=10^3;
% % Da_star = 76.9016; alpha= 20.25; gamma=  ;NTU_star1=188; CP_star=11.7; theta=0.25;
% % Fo_star=0.00293; psi=-0.045;

%%%%%%%%%%%%%%%Tf=500K For non runaway
% Da_star=.583;T_ad=2.5217*10^3;gamma=5.04;alpha=21.666;NTU_star1=284.882;
% psi=-1.82;
% Fo_star=4*10^-4;CP_star=2.4202; theta=1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Parameter set 1 for Tf=700K
%Da=3.6*10^-13;
Da=0.1;
%gamma=50.35;% 16.74456*10^3 gives ignition
gamma=20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameteet 2 for adiabatic Tf=700K
% Da_star=76.9016;T_ad=2.5217*10^3;gamma=3.60240;alpha=15.466;NTU_star1=0;
% psi=0;figure
% Fo_star=4*10^-4;CP_star=2.4202; theta=0.001;Da
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% paramter set 3 Tf=600K
% Da_star=5.8423;T_ad=2.5217*10^3;gamma=4.2028;alpha=18.0419;NTU_star1=244.1846;
% psi=-8.82;
% Fo_star=4*10^-4;CP_star=2.07452; theta=0.001;
% 
% 
z=linspace(0,1,200);% 100 nodes 
N=3*size(z,2);
deltax=z(2)-z(1); %step size
u(1)=1;
u(N/3+1)=0;
dudt=zeros(N-2,1);

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
%wall equations
i=2;
dudt((2*N)/3+i-1) =Fo*((u((2*N)/3+i)+u((2*N)/3+i-1)/(1+Bi*deltax)-2*u((2*N)/3+i-1))/deltax^2)+ Fo*(St1/CP*(u(N/3+i)-u((2*N)/3+i-1))+St2/CP*(vc-u(2*N/3+i-1)));
for i=3:N/3-2
    dudt((2*N)/3+i-1) = Fo*((u((2*N)/3+i)+u((2*N)/3+i-2)-2*u((2*N)/3+i-1))/deltax^2)+ Fo*(St1/CP*(u(N/3+i)-u((2*N)/3+i-1))+St2/CP*(vc-u(2*N/3+i-1)));
end
i=N/3-1;
dudt((2*N)/3+i-1)=Fo*((u((2*N)/3+i-2)+u((2*N)/3+i-1)/(1+Bi*deltax)-2*u((2*N)/3+i-1))/deltax^2)+ Fo*(St1/CP*(u(N/3+i)-u((2*N)/3+i-1))+St2/CP*(vc-u(2*N/3+i-1)));

end


    