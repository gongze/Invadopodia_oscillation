% oscillating invadopodia protrusion dynamics
%
% Copyright (C) Ze Gong, Postdoc researcher in Prof. Shenoy Lab, University of Pennsylvania
% Written in 2019, final revison at July 2020.

function Invadopodia_oscillation
% basic biological parameters
R0 = 5000;                                     % unit nm  initial filopodial length
Vp = 80;                                       % unit nm/s  polymerization speed
V0 = 300;                                      % unit nm/s unloaded retrograde flow
Fs0 = 300*1;                                   % unit pN  initial force
%Fa = 800;                                     % unit pN  diffusion term
Ec = 0.1;                                      % unit kPa  cytoskeleton/membrane stiffness
kesic = 5;                                     % unit kPa*s   cytoskeleton/membrane viscosity
A = pi*(400/2)^2/1000;                         % unit nm*um  filopodial cross area
etaR = 10;                                     % unit PN*s/nm  adhesion/clutch viscosity
alpha = 0.9;                                   % unit nm  feedback parameter alpha
beta = 1.5*etaR;                               % unit nm  feedback parameter beta 
tau1 = 400;                                    % unit s  myosin binding timescale
E0 = 0.5;                                      % unit kpa  ECM initial stiffness
E1 = 2;                                        % unit kpa  ECM stiffness
E2p = 0.25;
etas1 = 200;                                   % unit kpa*s  ECM viscosity
etas2p = 160;                                  % unit kpa*s  ECM viscosity
deltay = 10;                                   % unit pN  ECM yielding stress  
Rini = 000;
NoiseVar = 100;

% for equations parameters
a = R0/Ec/A*(etaR*Vp);  % = R0/Ec/A*(Ec*A+etaR*Vp);
b = R0/Ec/A*(1-Vp/V0);
c = R0/Ec/A*(etaR+kesic*A/R0);
g = R0/Ec/A/V0;
e = Fs0 + alpha*etaR*Vp;
f = alpha*etaR+beta;
h = R0/Ec/A;

% physical parameters group
Rss = a-b*e;
Fss = e;
tau2 = g*e+c;
const = b*f;

% bifurcation criteria
real = -(tau2+tau1-const);
imgine = ((tau2+tau1-const)^2-4*tau2*tau1)/tau1^2/tau2^2;

iternum = 7;
tstart = 0;
tfinal = 20000;
Ryield = zeros(iternum,1);
Tdetach = zeros(1, 1);
Ydetach = R0*ones(4, 1);
Tattach = zeros(1, 1);
Yattach = zeros(4, 1);
spara = [a, b, c, e, f, g, h, tau1, E0, E1, E2p, deltay, etas1, etas2p, R0, A, Rini, NoiseVar];
y0 = [Rini; Fss-00; 0; 0];
refine = 4;
tout =  tstart;
yout = y0';
yplastic = 0;
teout = [];
yeout = [];
ieout = [];
options1 = odeset('RelTol',1e-4,'Events', @events1);
options2 = odeset('RelTol',1e-4);

for i=1:iternum
    tspan = [tstart tfinal];
    [t, y, te, ie, ye] = ode45(@(t,y) Fun_plasticity_filopodial(t,y,spara), tspan, y0, options1);
    % Accumulate output.  This could be passed out as output arguments.
    nt = length(t);
    tout = [tout; t(2:nt)];
    yout = [yout; y(2:nt,:)];
    %yplastic = [yplastic; y(2:nt,1)-R0*(1+y(2:nt,3)/A/E0)];
    teout = [teout; te];          % Events at tstart are never reported.
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    tstart = t(nt);
    if (tstart==tfinal)
        break
    end
    % Set the new initial conditions
    y0 = [y(nt,1); y(nt,2); 0; y(nt,4)];
    Ry = y(nt,1);
    Ryield(i) = Ry;
    Tdetach(:,i) = t(nt);
    Ydetach(:,i) = y0;
    %%%%%%
    %   Retract beyond the attachment with ECM, there is no plasticity.    %
    %%%%%%
    tspan = [tstart tfinal];
    options2 = odeset(options2,'InitialStep',t(nt)-t(nt-refine),...
        'MaxStep',t(nt)-t(1));
    options2 = odeset(options2,'Events', @(t,y) events2(t,y,Ry));
    [t, y, te2, ie2, ye2] = ode45(@(t,y) Fun_filopodial(t,y,spara, Ry), tspan, y0, options2);
    % Accumulate output.  This could be passed out as output arguments.
    nt = length(t);
    tout = [tout; t(2:nt)];
    yout = [yout; y(2:nt,:)];
    teout = [teout; te2];          % Events at tstart are never reported.
    yeout = [yeout; ye2];
    ieout = [ieout; ie2];
    % yplastic = [yplastic; yplastic(end)*ones(nt-1,1)];
    % Set the new initial conditions
    y0 = [y(nt,1); y(nt,2); 0; y(nt,4)];
    Tattach(:,i) = t(nt);
    Yattach(:,i) = y0;
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation.  'refine' is 4 by default.
    options1 = odeset(options1,'InitialStep',t(nt)-t(nt-refine),...
        'MaxStep',t(nt)-t(1));
    tstart = t(nt);
    if (tstart==tfinal)
        break
    end
end
Vretro = Vp+(yout(:,1)-a+yout(:,2)*b+yout(:,3)*h)./(c+yout(:,2)*g);
% plotting the figures
figure(1)
subplot(3,1,1)
plot(tout/60, (yout(:,1))/1000)  % plot length
hold on
plot(Tdetach/60, (Ydetach(1,:))/1000, 'ro')
hold on
subplot(3,1,2)
plot(tout/60, yout(:,3))  % plot force
hold on
subplot(3,1,3)
plot(tout/60, yout(:,4)/1000)   % plot plastic displacement
hold on
figure(3)

subplot(2,1,1)
plot((yout(:,1))/1000, yout(:,2)) 
hold on
subplot(2,1,2)
plot((yout(:,1))/1000, yout(:,3))
hold on
end

function dydt = Fun_plasticity_filopodial(t, y, p)

[a, b, c, e, f, g, h, tau1, E0, E1, E2, deltay, etas1, etas2, R0, A, Rini, NoiseVar] =...
           deal(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10),p(11),p(12),...
                p(13),p(14),p(15),p(16),p(17),p(18));
dydt = zeros(4,1);
randnoise = randn(1);
dydt(1) = -(y(1)-a+y(2)*b+y(3)*h)/(c+y(2)*g);
dydt(2) =  1/tau1*(e-dydt(1)*f-y(2))+NoiseVar*randnoise/tau1;
% plasticity
hyb_tan = (1+2/pi*atan((y(3)-deltay-A*E2*y(4)/R0)/0.00001))/2;
%hyb_tan = (1+2/pi*atan((y(3)-deltay)/0.001))/2;
dydt(3) = A*E0/R0*dydt(1)-E0/etas1*((E0+E1)/E0*y(3)-A*E1*(y(1)-y(4)-Rini)/R0)-A*E0*dydt(4)/R0;
dydt(4) = R0/A*(y(3)-deltay-A*E2*y(4)/R0)*hyb_tan/etas2;
end

function dydt = Fun_filopodial(t, y, p, Ry)

[a, b, c, e, f, g, ~, tau1, NoiseVar] = deal(p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(18));
dydt = zeros(4,1);
randnoise = randn(1);
dydt(1) = -(y(1)-a+y(2)*b)/(c+y(2)*g);
dydt(2) =  1/tau1*(e-dydt(1)*f-y(2))+NoiseVar*randnoise/tau1;
dydt(3) = 0;
dydt(4) = 0;
end

function [value,isterminal,direction] = events1(t,y)
% Locate the time when ECM force passes through zero in a decreasing direction
% and stop integration.
value = y(3);     % detect ECM force = 0
isterminal = 1;   % stop the integration
direction = -1;   % negative direction

end 

function [value,isterminal,direction] = events2(t,y, Ry)
% Locate the time when length pass through retraction length in a decreasing direction
% and stop integration.
value = y(1)-Ry;     % detect protrusion length get to former detach point y(1)=Ry
isterminal = 1;   % stop the integration
direction = 1;   % positive direction
end