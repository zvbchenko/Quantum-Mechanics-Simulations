
% x0 = idpar(1);
% y0 = idpar(2);
%       
% deltax = idpar(3);
% deltay = idpar(4);
% px = idpar(5);
% py = idpar(6);

idtype = 1;
vtype = 2;
idpar = [0.0 , 0.5, 0.05, 0.05, 0.5, 0.0];
tmax = 0.01;
lambda = 0.005;
level = 6;
vpar = [0.43 , 0.44, 0.56, 0.57, 100000 ];

[x,y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);


avienable = 1;
avisecs = 1;
aviframerate = 15;
aviframes = avisecs * aviframerate;
avifilename = 'DoubleSlitCont2.avi';
if avienable
   avifreq = 1;
   aviobj = VideoWriter(avifilename);
   
   open(aviobj);

   
   title('Double Slit Scattering');


end


for time = 1:2: length(t)
    %psiT = psi(time, :, :);
    psiN = psire(:, :, time);
    
    
    %result = squeeze(psiN);

    
    contourf(x, y, psiN);
%     hold on;
%     surf(x, y, v);
    framecount = 10;

    
    xlabel("x");
    ylabel("y");
    xlim([0 1]);
    ylim([0 1]);
    zlim([-1 1]);
    
    drawnow
    for iframe = 1 : framecount
      writeVideo(aviobj, getframe(gcf));
    end
    pause(0.1);
%     hold off;
    
    
end

if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end