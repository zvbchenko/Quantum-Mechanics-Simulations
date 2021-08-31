l=[6, 7, 8, 9];

idtype = 0;
vtype = 0;

idpar = [3];
tmax = 0.25;
lambda = 0.1;

txt = cell(length(l),1);

for level = l
    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    [xp, tp, psip, psirep, psiimp, psimodp, probp, vp] = sch_1d_cn(tmax, level+1, lambda, idtype, idpar, vtype, vpar);
    
    [rows, cols] = size(psip);
    padxPsi = [psip zeros(rows, 1)];
    lastrow = padxPsi(rows,:);
    paddedPsip = [padxPsi; lastrow];
    Y = conv2(paddedPsip ,[1,1;1,1],'valid');
    Z = Y(1:2:end,1:2:end )/4;

    dpsi = (Z - psi);
    [row, col] = size(dpsi);
    sqdpsi = abs(dpsi).^2;
    l2norm = sqrt(sum(sqdpsi, 2)/col)*4^(level - 6);
    %%%% Calculate Exact %%%%
    psiexact = zeros(row, col);
    m = idpar(1);
    i= 1;
    for time = t
      psiexact(i, :) = exp(-1i*m^2*pi^2.*time).*sin(m*pi.*x);
      i = i +1;
    end
    Err = psiexact - psi;
    sqdErr = abs(Err).^2;
    l2Enorm =  sqrt(sum(sqdErr, 2)/col)*4^(level - 6);
    
    %%%%%%%%%%
    f1 =figure(1);
    
    plt = plot(t(1: end-1), l2norm(1: end-1));
    title("||d\psi^{l}||_2");
    ylabel("E");
    xlabel("time");
    index = 1 + level - l(1);
    if index == 1
        txt{index} = sprintf('||d\\psi^{%i}||_2',  level);
    elseif index == 2
        txt{index} = sprintf('4||d\\psi^{%i}||_2',  level);
    else
        txt{index} = sprintf('4^{{%i}}||d\\psi^{%i}||_2', index-1, level);
    end
    
    %f2 = figure(2);plot(t(1: end-1), l2Enorm(1: end-1));
    
    hold on
end
f1;
legend(txt)
hold off

%%%%%%%%%Comparison with exact solutions

newtxt = cell(length(l),1);
for level = l
    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);


   
    [row, col] = size(psi);
    %%%% Calculate Exact %%%%
    psiexact = zeros(row, col);
    m = idpar(1);
    i= 1;
    for time = t
      psiexact(i, :) = exp(-1i*m^2*pi^2.*time).*sin(m*pi.*x);
      i = i +1;
    end
    Err = psiexact - psi;
    sqdErr = abs(Err).^2;
    l2Enorm =  sqrt(sum(sqdErr, 2)/col)*4^(level - 6);
    
    %%%%%%%%%%
    f2 =figure(2);
    
    plot(t(1: end-1), l2Enorm(1: end-1));
    title("||E(\psi^{l})||_2");
    ylabel("Error");
    xlabel("time");
    index = 1 + level - l(1);
    if index == 1
        newtxt{index} = sprintf('||E^{%i}||_2',  level);
    elseif index == 2
        newtxt{index} = sprintf('4||E^{%i}||_2',  level);
    else
        newtxt{index} = sprintf('4^{{%i}}||E^{%i}||_2', index-1, level);
    end
    
    %f2 = figure(2);plot(t(1: end-1), l2Enorm(1: end-1));
    
    hold on
end
f2;
legend(newtxt)
hold off

