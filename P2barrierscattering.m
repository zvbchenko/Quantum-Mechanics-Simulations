
% x0 = idpar(1);
% y0 = idpar(2);
%       
% deltax = idpar(3);
% deltay = idpar(4);
% px = idpar(5);
% py = idpar(6);


idtype = 1;
vtype = 1;
idpar = [0.0 , 0.0, 0.05, 0.05, 0.0, 0.0];
tmax = 0.05;
lambda = 0.05;
level = 7;
vpar = [0.4, 0.6, 0.4, 0.6, 10000 ];

[x,y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);



avienable = 1;
avisecs = 1;
aviframerate = 15;
aviframes = avisecs * aviframerate;
avifilename = 'barrierScatteringSurf.avi';
if avienable
   avifreq = 1;
   aviobj = VideoWriter(avifilename);
   
   open(aviobj);

   
   title('Barrier Scattering');


end


for time = 1:2: length(t)
    %psiT = psi(time, :, :);
    psiN = psire(:, :, time);
    
    
    %result = squeeze(psiN);
    
    surf(x, y, psiN);
    
    
    xlabel("x");
    xlabel("y");
    xlim([0 1]);
    ylim([0 1]);
    zlim([-1 1]);
    framecount = 10;
    for iframe = 1 : framecount
      writeVideo(aviobj, getframe(gcf));
    end
    drawnow
    pause(0.1);
    
    
end


% Close avi file ...
if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end


