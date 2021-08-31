l=[6, 7, 8, 9];

idtype = 1;
vtype = 0;

idpar =  [0.50 0.075 0.0];
tmax = 0.01;
lambda = 0.01;

txt = cell(length(l),1);

for level = l
    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    [xp, tp, psip, psirep, psiimp, psimodp, probp, vp] = sch_1d_cn(tmax, level+1, lambda, idtype, idpar, vtype, vpar);

    [rows, cols] = size(psip);
    padxPsi = [psip zeros(rows, 1)];
    lastrow = padxPsi(rows,:);
    paddedPsip = [padxPsi; lastrow];
    %paddedPsimodp = padarray(psimodp, [1 1], 'circular');
    Y = conv2(paddedPsip ,[1,1;1,1],'valid');
    Z = Y(1:2:end,1:2:end )/4;

    dpsi = (Z - psi);
    [row, col] = size(dpsi);
    sqdpsi = abs(dpsi).^2;
    l2norm = sqrt(sum(sqdpsi, 2)/col)*4^(level - 6);
    title("||d\psi^{l}||_2 of Boosted Gaussian");
    ylabel("E");
    xlabel("time");
    plot(t(1: end-1), l2norm(1: end-1))
    index = 1 + level - l(1);
    if index == 1
        txt{index} = sprintf('||d\\psi^{%i}||_2',  level);
    elseif index == 2
        txt{index} = sprintf('4||d\\psi^{%i}||_2',  level);
    else
        txt{index} = sprintf('4^{{%i}}||d\\psi^{%i}||_2', index-1, level);
    end
    hold on
end
legend(txt)
hold off