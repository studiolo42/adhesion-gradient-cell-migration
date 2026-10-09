%% adhesion_gradient_model.m
% Simulates a two-ended (front/back) cell moving along an adhesion gradient
% and plots: position, local flow, and proportion of bound linkers vs. time.
% To get a different migration mode, change b0 (beta_0) and T below.
% Naming: f = front end of the cell, b = back end of the cell.

clear; clc; close all;       % clear variables, clear the command window, close old figures

% ---------------- Settings ----------------
b0 = 9;                      % beta_0: base value of beta (we used 2, 2.5, 6.5, 9)
T  = 20;                     % total simulated time [a.u]
dt = 1e-4;                   % time step of the simulation

% ---------------- Parameters ----------------
c = 3.95;                    % constant in the local-flow equations
d = 3.95;                    % constant in the local-flow equations
delta = 200;                 % how fast the local flow relaxes to its target value
k  = 0.8;                    % spring constant connecting front and back
fs = 5;                      % force scale of linker unbinding
kp = 20;                     % friction prefactor when an end retracts
sig_v = 0.01;                % strength of the random noise added to the local flow
n_min = 1e-5;                % smallest allowed bound-linker proportion (keeps n above 0)

A = 27;                      % peak value of the adhesion profile
K_g = 0.1769;                % steepness of the adhesion gradient
r0 = 0.75;                   % constant in r/(r+r0), sets when adhesion saturates
alpha = 0.65;                % how strongly beta depends on adhesion
x_peak = 65/2;               % position of the adhesion peak (middle of a track of length 65)
x_end  = 2*x_peak;           % end of the track: the run stops if a cell end reaches it

lc = c/(sqrt(c*b0/d) - 1);   % length scale used to set the initial cell length

% ---------------- Model equations ----------------
r_of_x = @(x) A*exp(-K_g*abs(x - x_peak));          % adhesion r at position x (peaks at x_peak)
b_of_r = @(r) b0*(1 + alpha*(r./(r + r0)));         % beta at an end, set by the adhesion there

% Velocity of each end before friction (xf-xb-1 = stretch of the spring, rest length 1)
F1 = @(xf,xb,vf,rf) (rf/(rf+r0))*vf - k*(xf-xb-1);      % front end
F2 = @(xf,xb,vb,rb) -(rb/(rb+r0))*vb + k*(xf-xb-1);     % back end

% Rate of change of the local flow; it relaxes (rate delta) towards a target set by beta
H1 = @(xf,xb,vf,vb,bf) -delta*(vf - bf*(1/(1 + (c*(vf-vb)/d)/(exp((vf-vb)*(xf-xb)/d) - 1))));   % front
H2 = @(xf,xb,vf,vb,bb) -delta*(vb - bb*(1/(1 + (c*(vf-vb)/d)/(1 - exp(-(vf-vb)*(xf-xb)/d)))));  % back

% Rate of change of the bound-linker proportion: binding r*(1-n) minus force-dependent unbinding
G1 = @(xf,xb,vf,nf,rf) rf*(1-nf) - nf*exp((-(rf/(rf+r0))*vf + k*(xf-xb-1))/fs/nf);   % front
G2 = @(xf,xb,vb,nb,rb) rb*(1-nb) - nb*exp((-(rb/(rb+r0))*vb + k*(xf-xb-1))/fs/nb);   % back

% ---------------- Initial conditions ----------------
N = round(T/dt);             % number of time steps
xf = zeros(1,N);             % front position over time (filled with zeros first)
xb = zeros(1,N);             % back position over time
vf = zeros(1,N);             % front local flow over time
vb = zeros(1,N);             % back local flow over time
nf = zeros(1,N);             % front bound-linker proportion over time
nb = zeros(1,N);             % back bound-linker proportion over time
rf = zeros(1,N);             % adhesion at the front over time
rb = zeros(1,N);             % adhesion at the back over time
bf = zeros(1,N);             % beta at the front over time
bb = zeros(1,N);             % beta at the back over time

xb(1) = 15;                  % initial back position
xf(1) = lc + 5 + 15;         % initial front position (initial cell length = lc + 5)
nf(1) = 0.25;                % initial front bound-linker proportion
nb(1) = 0.25;                % initial back bound-linker proportion
vf(1) = b0*(xf(1)-xb(1))/(c + xf(1)-xb(1));   % initial front local flow
vb(1) = vf(1);               % initial back local flow (equal to the front)
rf(1) = r_of_x(xf(1));       % initial adhesion at the front
rb(1) = r_of_x(xb(1));       % initial adhesion at the back
bf(1) = b_of_r(rf(1));       % initial beta at the front
bb(1) = b_of_r(rb(1));       % initial beta at the back

% ---------------- Time loop ----------------
for j = 2:N                  % step through time, one step per pass

    rf(j) = r_of_x(xf(j-1)); % adhesion at the front's last position
    rb(j) = r_of_x(xb(j-1)); % adhesion at the back's last position
    bf(j) = b_of_r(rf(j));   % beta at the front
    bb(j) = b_of_r(rb(j));   % beta at the back

    dxf = F1(xf(j-1), xb(j-1), vf(j-1), rf(j-1));   % front velocity before friction
    dxb = F2(xf(j-1), xb(j-1), vb(j-1), rb(j-1));   % back velocity before friction

    pf = rf(j-1)/(rf(j-1) + r0);                    % adhesion factor r/(r+r0) at the front
    if dxf > 0                                      % front is moving forward:
        gf = pf;                                    %   friction is the adhesion factor
    else                                            % front is retracting:
        gf = nf(j-1)*kp*exp((pf*vf(j-1) - k*(xf(j-1)-xb(j-1)-1))/fs*nf(j-1));   % friction from bound linkers
    end
    dxf = dxf/gf;                                   % front velocity after friction

    pb = rb(j-1)/(rb(j-1) + r0);                    % adhesion factor r/(r+r0) at the back
    if dxb < 0                                      % back is moving backward:
        gb = pb;                                    %   friction is the adhesion factor
    else                                            % back is retracting:
        gb = nb(j-1)*kp*exp((pb*vb(j-1) - k*(xf(j-1)-xb(j-1)-1))/fs*nb(j-1));   % friction from bound linkers
    end
    dxb = dxb/gb;                                   % back velocity after friction

    dnf = G1(xf(j-1), xb(j-1), vf(j-1), nf(j-1), rf(j-1));   % rate of change of front bound linkers
    dnb = G2(xf(j-1), xb(j-1), vb(j-1), nb(j-1), rb(j-1));   % rate of change of back bound linkers

    if vf(j-1) == vb(j-1)                                 % if both flows are equal, H1 and H2 are 0/0,
        l = xf(j-1) - xb(j-1);                            %   so use their limit instead; l = cell length
        dvf = -delta*(vf(j-1) - bf(j)*l/(c + l));         %   front local-flow rate (limit form)
        dvb = -delta*(vb(j-1) - bb(j)*l/(c + l));         %   back local-flow rate (limit form)
    else                                                  % otherwise use the normal equations:
        dvf = H1(xf(j-1), xb(j-1), vf(j-1), vb(j-1), bf(j-1));   % front local-flow rate
        dvb = H2(xf(j-1), xb(j-1), vf(j-1), vb(j-1), bb(j-1));   % back local-flow rate
    end

    xf(j) = xf(j-1) + dxf*dt;                       % new front position
    xb(j) = xb(j-1) + dxb*dt;                       % new back position
    vf(j) = vf(j-1) + dvf*dt + sig_v*randn*sqrt(delta*dt);   % new front flow, plus random noise
    vb(j) = vb(j-1) + dvb*dt + sig_v*randn*sqrt(delta*dt);   % new back flow, plus random noise
    nf(j) = max(n_min, nf(j-1) + dnf*dt);           % new front bound linkers (never below n_min)
    nb(j) = max(n_min, nb(j-1) + dnb*dt);           % new back bound linkers (never below n_min)

    if xf(j) >= x_end || xb(j) >= x_end             % if either end reached the end of the track,
        break;                                      %   stop the simulation
    end
end

% ---------------- Keep only the simulated part ----------------
xf = xf(1:j);  xb = xb(1:j);    % cut positions at the last step j (matters if stopped early)
vf = vf(1:j);  vb = vb(1:j);    % cut local flows
nf = nf(1:j);  nb = nb(1:j);    % cut bound-linker proportions
t  = (0:j-1)*dt;                % time axis [a.u]
xcm = (xf + xb)/2;              % centre position of the cell

% ---------------- Plots ----------------
% Position vs. time (positions shown relative to the adhesion peak)
figure;                                                          % new figure
plot(t, xf - x_peak, 'r', t, xcm - x_peak, 'g', t, xb - x_peak, 'b', 'LineWidth', 1.5);   % front red, centre green, back blue
yline(0, 'k--', 'HandleVisibility', 'off');                      % dashed line at the adhesion peak
xlabel('Time [a.u]');                                           
ylabel('Position [a.u]');                                      
legend('x_f', 'x_{cm}', 'x_b');                                 

% Local flow vs. time
figure;                                                         
plot(t, vf, 'r', t, vb, 'b', 'LineWidth', 1.5);                 
xlabel('Time [a.u]');                                           
ylabel('Local flow [a.u]');                                   
legend('v_f', 'v_b');                                        

% Proportion of bound linkers vs. time
figure;                                                       
plot(t, nf, 'r', t, nb, 'b', 'LineWidth', 1.5);              
xlabel('Time [a.u]');                                      
ylabel('Prop. of bound linkers');                           
ylim([0 1]);                                                    
legend('n_f', 'n_b');                                            
