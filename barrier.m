tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1;
idpar = [0.40, 0.075, 20.0]; % 
vtype = 1 ;

potentials = linspace(exp(-2), exp(5), 251);
Fe = zeros(length(potentials));
i = 1;
for V = potentials
    vpar = [0.6, 0.8, V];
    x1 = 0.8;
    x2 = 1.0;


     [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);


    temporal_average = mean(prob);
    dx = x(2) - x(1);
    index_low = round(x1/dx); % produce index of x_min
    index_high= round(x2/dx);

    Fe(i) = (temporal_average(index_high) - temporal_average(index_low))/(x2 - x1);
    i = i+1;
    
end

plot(log(potentials), log(Fe))
title("The dependence of ln(Fe(x1, x2)) on ln(V0)")
xlabel("ln(V0)")
ylabel("ln(Fe(x1, x2))")
