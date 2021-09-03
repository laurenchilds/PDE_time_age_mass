clear all; close all; clc
tic


%%%% Parameters %%%%

t0 = 0.1; % start time
a0 = 0; % minimum age
m0 = 0.1; % minimum mass

t_step = 0.1; % time step-size
a_step = 0.05; % age step-size
m_step = 0.05; % mass step-size

t_max = 7; % max time
a_max = 8; % max age
m_max = 5; % max mass

tvec = t0:t_step:t_max;
avec = a0:a_step:a_max; 
mvec = m0:m_step:m_max;

%%%% Include parameters in a vector to pass %%%
par.a0 = a0;
par.m0 = m0;

par.a_step = a_step;
par.m_step = m_step;

par.a_max = a_max;
par.m_max = m_max;

%%%% Pre-allocation %%%%
f = NaN(length(tvec),length(avec), length(mvec));

%%%% Pre-allocation for testing variables %%%%
new_m_out = NaN(length(tvec),length(avec), length(mvec));
new_a_out = NaN(length(tvec),length(avec), length(mvec));
psi_out = NaN(length(tvec),length(avec), length(mvec));
exp_out = NaN(length(tvec),length(avec), length(mvec));

%%%% Normalization Factors for Psi (initial condition) function %%%%
psiIC = NaN(length(avec),length(mvec));
for a1 = 1:length(avec)
    for m1 = 1:length(mvec)
        a = avec(a1);
        m = mvec(m1);
        psiIC(a1,m1) = psi_fun(a,m,par,1);
    end
end
psiICnorm = trapz(avec,trapz(mvec,psiIC'));

%%%% Run the Simulation %%%%
f(1,:,:) = psiIC/psiICnorm;

for t1 = 2:length(tvec) % iterate through t
    t = tvec(t1);
    tspan = t0:t_step:t; % t' runs t0 to current t
    for a1 = t1:length(avec) % iterate through a
        a = avec(a1);
        for m1 = 1:length(mvec) % iterate through m
            m = mvec(m1);
            
            new_m = new_m_fun(tspan,a,m);
            
            deltaintout = deltaint_fun(tspan,a,m);
            int_delt = deltaintout(end);
            
            new_a = a - t +t0;
            psi = psi_fun(new_a,new_m,par,psiICnorm);
            
            f(t1,a1,m1) = psi*exp(-int_delt);
            
            % Output variables for testing
            new_m_out(t1,a1,m1) = new_m;
            new_a_out(t1,a1,m1) = new_a;
            psi_out(t1,a1,m1) = psi;
            exp_out(t1,a1,m1) = int_delt;
        end
    end
end
%}
toc
filename = 'Output.mat';
save(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gTAout = gTA_fun(t,a) 
% time-age component of growth function
gTAout = 1; % constant
%gTAout = 1/2; % constant
%gTAout = 1/10; % constant
%gTAout = a; % linear in age
%gTAout = a/10; % linear in age
end

function gTAintout = gTAint_fun(xspan,a)
% integral of time-age component of growth function
% Note: integral runs t0 to t in m0 
%       integral runs t' to t in m
%       xspan(1)  starting point
%       xspan(end)  t
gTAt = arrayfun(@(x) gTA_fun(x,x+a-xspan(end)), xspan);
gTAintout = arrayfun(@(x) trapz(xspan(1:x),gTAt(1:x)),2:length(xspan));
end

function dout = d_fun(t,a,m)
% original right hand side function of the  PDE
dout = 0*(t+a+m); % zero scaled for matrix dimension
%dout = a/2; % linear in age
%dout = m/2; % linear in mass
%dout = a.*m/4; % quadratic
end

function deltaout = delta_fun(t,a,m)
% right hand side function of the PDE after rearrangement
deltaout = gTA_fun(t,a)*dgM_fun(m) + d_fun(t,a,m);
end

function deltaintout = deltaint_fun(tspan,a,m)
% integral of right hand side function of the PDE after rearrangement
% Note: integral runs t0 to t
%       interior integral in m, runs t' to t (called xspan)
deltat = zeros(length(tspan),1);
for tcounter = 1:length(tspan)
    t = tspan(tcounter);
    xspan = tspan(tcounter:end);
    if length(xspan)>1
        new_m = new_m_fun(xspan,a,m);
    else
        new_m = m;
    end
    deltat(tcounter) = delta_fun(t,t+a-tspan(end),new_m);
end
    deltaintout = arrayfun(@(x) trapz(tspan(1:x),deltat(1:x)),2:length(tspan));
end

function new_m_out = new_m_fun(xspan,a,m)
% evaluates revised m value
% Note: requires integral of gTA function across xspan
gTAintout = gTAint_fun(xspan,a);
tmp = G(m) - gTAintout(end);
new_m_out = invG(tmp);
end

function dgMout = dgM_fun(m) 
% derivative of gM_fun with respect to mass m.
[dgMout, ~, ~] = AllGm(m);
end

function Gout = G(m)
% antiderivative of 1/gM_fun
[~, Gout, ~] = AllGm(m);
end

function invGout = invG(m)
% inverse of anti-derivative of 1/gM_fun
[~, ~, invGout] = AllGm(m);
end

function [dgMout, Gout, invGout] = AllGm(m)
% computes:
%   derivative of gM_fun with respect to mass m.
%   antiderivative of 1/gM_fun
%   inverse of anti-derivative of 1/gM_fun
% Note: must comment out in groups to retain the same functions

%
% Constant gM = c
c = 1; % c = 1/2; 
dgMout = 0;
Gout = m/c;
invGout = c*m;
%}

%{
% Linear gM = c*m
c = 1; %c = 1/2;  c = 1/5;
dgMout = c;
Gout = log(m)/c;
invGout = exp(c*m);
%}

%{
% Asymptotic-2 gM = m^2/(m^2 + 1) 
dgMout = (2*m)/(m^2 + 1) - (2*m^3)/(m^2 + 1)^2;
Gout = m - 1/m;
invGout = m/2 + (m^2 + 4)^(1/2)/2;
%}

%{
% Asymptotic-2 gM = m^2/(m^2 + 4) 
c = 4;
dgMout = (2*m)/(m^2 + c) - (2*m^3)/(m^2 + c)^2;
Gout = m - c/m;
invGout = m/2 + (m^2 + c^2)^(1/2)/2;
%}

end

function psiout = psi_fun(a,m,par,psinorm) 
% initial age-mass distribution Psi(a,m)
a0 = par.a0; % minimum age
m0 = par.m0; % minimum mass

a_max = par.a_max; % maximum age considered
m_max = par.m_max; % maximum mass considered

initial_density = 100;


mu = [0.4,0.6]; % vector of mean age and size at time t0
sigma = [0.1 0; 0 0.1]; % matrix of standard deviations in age and size at time t0

if imag(m)>0
    y2 = 0;
else
if a>=a0 && m>=m0 && a <= a_max && m<= m_max
    X = [a,m];
    y2 = mvnpdf(X,mu,sigma); % Bivariate normal evaluated at points in X
else
    y2 = 0;
end
end

psiout = initial_density*y2/psinorm;
end
