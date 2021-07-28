function [bold, neu, flow] = hemodynamic_model(Input,time_step,par)


% temporal integration step:
dt         = time_step;
N          = size(Input,2);
T          = size(Input,1); % number of time points


A          = -1;
C          = eye(N)/16;
mu         = 0.7;
lam        = 0.3;

alpha      = 0.35;
tt         = 2;    % mean transit time through vasculature
visco      = 10;   % viscoelasticity time constant of the vessel
n          = 4;    % ratio (dCBF-1)/(dCMRO2-1)
decay1     = 1*par;  % NVC
gain       = 2*par;  % NVC
decay2     = 1*par;  % NVC


% BOLD parameters
TE = 0.028;
V0 = 3.5;
T2s_blood  = 0.012;
T2s_tissue = 0.028;
lambda  = 1.05;
epsilon = lambda.*exp(-TE./T2s_blood)./exp(-TE./T2s_tissue);

E0    =  0.4; % resting oxygen extraction fraction
r0     = 125;
Hct    = 0.38;
B0     = 7;
gyro   = 2*pi*42.6*10^6;
suscep = 4*pi*0.264*10^-6;
nu0   = suscep/(4*pi)*gyro*Hct*B0;

%-Coefficients in BOLD signal model
%==========================================================================
k1  = 4.3*nu0.*E0.*TE;
k2  = epsilon.*r0.*E0.*TE;
k3  = 1 - epsilon;

%time_ax = linspace(dt,T*dt,T) - onset;

X  = zeros(N,6);
y  = X;
yd = y;

for t = 1:T
    
    X(:,3:5) = exp(X(:,3:5));
    % %--------------------------------------------------------------------------
    dy(:,1)  =  A*X(:,1) - mu*X(:,6) + C*Input(t,:)';
    dy(:,6)  = lam.*(-X(:,6) + X(:,1));
    % vascoactive signal
    dy(:,2)  =  X(:,1) - decay1.*X(:,2);
    % CBF
    dy(:,3)  = (gain.*X(:,2) - decay2.*(X(:,3)-1))./X(:,3);
    
    % outflow
    fv       = (tt.*X(:,4).^(1./alpha) + visco.*X(:,3))./(tt+visco);
    
    % CMRO2
    % %--------------------------------------------------------------------------
    m       = (X(:,3)+n-1)./n;
    
    % CBV
    dy(:,4)  = (X(:,3) - fv)./(tt.*X(:,4));
 
    % deoxy-Hb
    dy(:,5)  = (m - fv.*X(:,5)./X(:,4))./(tt.*X(:,5));
    y        = y + dt*dy;
    X        = y;
%    Y(t,:)   = y;
    
    
    %Output equation of BOLD signal model
    %==========================================================================
    f        = exp(y(:,3));
    v        = exp(y(:,4));
    q        = exp(y(:,5));
    
    neu(t,:)   = y(:,1);
    flow(t,:)  = (f-1)*100;
    bold(t,:)  = V0.*(k1.*(1 - q) + k2.*(1-q./v)+ k3.*(1-v));  % in % signal change
    
end



%figure, plot([neu,cbf]);



