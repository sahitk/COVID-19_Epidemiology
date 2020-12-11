
totalpeople = 1928000
initial.S = 1927999;   %set the initial value of 'S'
initial.I = 1;   %set the initial value of 'I'
initial.R = 0;   %set the initial value of 'R'
initial.D = 0

param.gamma = 0.359;   %set the parameter 'r' of the model
param.beta = (2.2 ∗ param.gamma) / (initial.S + initial.I + initial.R)    %set the parameter 'beta' of the model
param.mu = 0.067

%days from January 1st, used for seasonal transmission for param.beta
d = 1;


end_time = 80;
time_interval = 1;

%percent vaccinated by end_time
pve=0.967

%percent vaccinated per_day
pvd=pve/30

%day when vaccine is released
dayofvaccine=25



N = initial.S + initial.R + initial.I + initial.D;
R0 = param.beta ∗ N / param.gamma;

x=[initial.S; initial.I; initial.R; initial.D];


initial_values = [];
variable_names = fieldnames(initial);
for i=1:length(variablenames)
    initialvalues = [initialvalues; initial.(variablenames{i})];
end

%integrate the ODE system
%[t, y] = ode45(@(t, x) odesystem(0, x, param),[0 endtime], initialvalues,[]);

%prepare legend texts
legend_texts = cell(length(variablenames), 1);


pS = 0;
pI = 0;
pR = 0;
pD = 0;
oldS = initial.S;
oldI = initial.I;
oldR = initial.R;
oldD = initial.D;


seasonalbeta=0;


smatrix = [];
imatrix = [];
rmatrix = [];
dmatrix = [];
ymatrix = [];


for i=0:end_time
    %y = odesystem(i, x, param)

   seasonalbeta = param.beta ∗ (1 + (0.35 ∗ cos((2 ∗ pi * d) / 365)))
   param.beta


    if (i >= dayofvaccine)
        dS=−seasonalbeta ∗ x(1) ∗ x(2) − (0.00002142 ∗ x(1)) + (0.0000461 ∗ x(1)) − (pvd ∗ x(1));
        dI = +seasonalbeta ∗ x(1) ∗ x(2) − param.gamma ∗ x(2) − (x(2) ∗ param.mu);
        dR = (param.gamma ∗ x(2)) + (pvd ∗ x(1));
        dD = (param.mu ∗ x(2))
    else
        dS = −seasonalbeta ∗ x(1) ∗ x(2) − (0.00002142 ∗ x(1)) + (0.0000461 ∗ x(1));
        dI = +seasonalbeta ∗ x(1) ∗ x(2) − param.gamma ∗ x(2) − (x(2) ∗ param.mu);
        dR = (param.gamma ∗ x(2));
        dD = param.mu ∗ (x(2))
    end

    pS = oldS + dS;
    pI = oldI + dI;
    pR = oldR + dR;
    pD = oldD + dD


    %plot(i, pS)
    %plot(i, pI)
    %plot(i, pR)
    pS
    pI
    pR
    pD
    smatrix = horzcat(smatrix, pS)
    imatrix = horzcat(imatrix, pI)
    rmatrix = horzcat(rmatrix, pR)
    dmatrix = horzcat(dmatrix, pD)
    ymatrix = horzcat(ymatrix, i)


    oldS = pS;
    oldI = pI;
    oldR = pR;
    oldD = pD;

    x(1) = oldS;
    x(2) = oldI;
    x(3) = oldR;
    x(4) = oldD;
    d = d + 1;

end


plot(ymatrix, smatrix, ymatrix, imatrix, ymatrix, rmatrix, ymatrix, dmatrix);
clear xlabel;
clear ylabel;
clear title;
title('Day 20 Vaccine, 80%', 'fontsize', 18);
xlabel('Days since first case');
ylabel('People');
