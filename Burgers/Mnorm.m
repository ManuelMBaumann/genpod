function err = Mnorm( F, G, M )

% sqrt ( \int_t0^te ||F-G||_M^2 dt  )

global dt Nt


diff = F(:, 1) - G(:, 1);
err_old = 0.5 * dt * diff'*M*diff;
err = err_old;

for k = 2:Nt
        
    diff = F(:, k) - G(:, k);
    err_new = 0.5 * dt * diff'*M*diff;
    err = err + err_old + err_new;
    err_old = err_new;

end

err = sqrt(err);