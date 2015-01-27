j = find(N==current_n);
if (test_sep && (N(j) == N(end-1)));
  disp('#######################################')
  disp('TESTING FOR SEPARABILITY')
  disp('#######################################')
  %-------------------------------------------------
  % Generate "probe" vector p_h
  %-------------------------------------------------
  % rml randn('state',0);
  % rml p_0h = randn(size(v));
  p_0h = 2*(rand(size(v)) - 0.5);
  p_iter_max = 5;
  p_1h = p_0h;
  for jj=1:p_iter_max;   % Iterate to emphasize the low "frequencies"
	p_1h = update(downdate(p_1h,0),0);
  end;
  p_h  = p_0h - p_1h;
  %------------------------------------------------
  % Downdate the probe vector and compute
  % the matrix-vector product Ap_h0
  %------------------------------------------------
  p_H  = downdate(p_h,0);
  h_diff = sqrt(eps)*norm(p_h);
  [Fvh, Gvh]  = sfun(v + h_diff*p_h);
  Ap_h0 = (Gvh - Gv) / h_diff;
  %------------------------------------------------
  % Plot the various results
  % In the general case the FFT may not be an
  % appropriate tool to analyze the results.
  %------------------------------------------------
  x = X(1:N(j),j);
  nx = length(x);
  dx = x(2) - x(1);
  wt = 0;

  figure(135)
  yy =linspace(0,1,length(v))';
  plot(yy,Ap_h0)
  title('Hessian-vector product on fine grid')
  Ap_h = downdate(Ap_h0,0);
  Ap_h = Ap_h0; % RML hack.
  figure(139); subplot(2,2,1)
  plot(Ap_h); title('Original test Hessian-vector product')
  n_Ap_h = norm(Ap_h) + wt*norm(diff(Ap_h)/dx);
  ww_axis = axis;
  figure(139); subplot(2,2,3);
  omega = linspace(0, 2*(nx-1)/nx*(pi/dx), nx)' - pi/dx;
  plot(omega, dx*abs(fftshift(fft(Ap_h)))); title('FFT of original test Hessian-vector product');
  % plot(abs(fftshift(fft(Ap_h)))); title('FFT of filtered Hessian-vector product')
  vv_axis = axis;
  vv_axis(1) = omega(1);
  vv_axis(2) = omega(end);
  axis(vv_axis);

  %-----------------------------------------------
  % Iterate to emphasize the low frequencies.
  % If the Hessian is not separable, then the
  % low frequencies will be large; otherwise
  % they will be small.
  %-----------------------------------------------
  A_iter_max = 5;
  for jj = 1:A_iter_max
	%RML hack.         Ap_h = downdate(update(Ap_h,0),0);
	Ap_h = update(downdate(Ap_h,0),0);
  end
  %------------------------------------------------
  % Plot the various results
  % In the general case the FFT may not be an
  % appropriate tool to analyze the results.
  %------------------------------------------------
  figure(139); subplot(2,2,2)
  plot(Ap_h); title('Filtered test Hessian-vector product')
  axis(ww_axis);
  fprintf('Scaled |Ap| = %8.2e\n',(norm(Ap_h)+ wt*norm(diff(Ap_h)/dx))/n_Ap_h);
  figure(139); subplot(2,2,4)
  plot(omega, dx*abs(fftshift(fft(Ap_h)))); title('FFT of filtered Hessian-vector product')
  % plot(abs(fftshift(fft(Ap_h)))); title('Iterated FFT of matrix-vector product')
  axis(vv_axis)

  vilename = sprintf('plot%02d', picture_count);
  figure(139); size(Ap_h)
  print('-depsc2', vilename);
  picture_count = picture_count+1;
  %     disp('Hit any key to continue')
  %     pause
  %     keyboard;
end;