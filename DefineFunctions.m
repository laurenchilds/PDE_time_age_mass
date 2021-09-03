clear; close all; clc


[gTA,gM,dgM] = Define_g_Derivg; % Define growth function and its derivative w.r.t. mass m.
G = Define_gM_Antideriv(gM);
invG = Define_G_inverse(G);
intG = Define_g_integral(gM);

gTA
gM
dgM
G
invG
intG

%%%%%%%%%%%%%% functions %%%%%%%%%%%%
function intG = Define_g_integral(gM)

intG = int(gM);
end

function [gTA,gM,dgM] = Define_g_Derivg

syms gM(m)
syms gTA(t,a)
assume(m,'real')

gM(m) = 2*m;

dgM = diff(gM,m);
gTA(t,a) = a;

end

function G = Define_gM_Antideriv(gM) % G is the antiderivative of 1/gM w.r.t. mass m.

syms G(m) x
assume(m<2)
gMx = subs(gM(m),m,x);
G(m) = int(1/gMx,1,m);

end

function invG = Define_G_inverse(G) % G is the antiderivative of 1/gM w.r.t. mass m.

syms invG(m)
assume(m < 2);
invG(m) = finverse(G,m);

end