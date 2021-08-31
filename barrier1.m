tmax = 0.15;
level = 7;
lambda = 0.01;
idtype = 1;
idpar = [0.40, 0.075, 20.0]; % 
vtype = 1 ;
vpar = [0.5, 0.8, 1000.0];



 [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
 
 
avienable = 1;
avisecs = 1;
aviframerate = 15;
aviframes = avisecs * aviframerate;
avifilename = 'OneDWell.avi';
if avienable
   avifreq = 1;
   aviobj = VideoWriter(avifilename);
   
   open(aviobj);

   
   title('1D Well ');


end
 for m = 1: length(t)
      
      plot(x, psimod(m, :), x, v, 'b-', 'LineWidth', 3);
      ylim([-2 2]);
      ylabel("|\psi(x, t)|^2");
      xlabel("x");
      title(sprintf("|\\psi(x, %f)|^2", t(m)));
      framecount =1;
      for iframe = 1 : framecount
        writeVideo(aviobj, getframe(gcf));
      end
      
      drawnow

 end

 if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end