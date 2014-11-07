fid = fopen('betaPareto.txt','a');
global Beta;
for it=-6:2:-2
  
    Beta=10^(it);optXTM_XRF;[f,g,f1,f2]=feval(fctn,xstar);  

  fprintf(fid,'%12.4e    %12.4e     %12.4e    %12.4e      %12.4e     %f\n',Beta,f,f1,f2,errTol,t); 
end
  fclose(fid);