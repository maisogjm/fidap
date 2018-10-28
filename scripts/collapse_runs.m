% Script to show average fMRI run with
% effects of no interest subtracted out.

if (exist('Nruns')==0)
    Nruns=input('How many fMRI runs? ');
end
Ndat=num_dat_pts/Nruns;

% Subtract out effects of no interest.
Y_int=Y-X(:,1:num_ind_var-num_var_int)*B_hat(1:num_ind_var-num_var_int);
Y_hat_int=X(:,num_ind_var-num_var_int+1:num_ind_var)*B_hat(num_ind_var-num_var_int+1:num_ind_var);

% Calculate average run.
Y_int_avg=zeros(Ndat,1);
Y_hat_int=zeros(Ndat,1);
for i=0:Nruns-1
    Y_int_avg=Y_int_avg+Y_int(i*Ndat+1:i*Ndat+Ndat);
    Y_hat_int=Y_hat_int+X(i*Ndat+1:i*Ndat+Ndat,num_ind_var-num_var_int+1:num_ind_var)*B_hat(num_ind_var-num_var_int+1:num_ind_var);
end

Y_int_avg=Y_int_avg/Nruns;
Y_hat_int=Y_hat_int/Nruns;

hold off
plot(Y_int_avg)
hold
plot(Y_hat_int,'b')
hold off
