% MATLAB script to display data
% analyzed with (multivariate) multiple regression.

% Query user for number of variables of interest
% and for number of timepoints per actual run.
num_var_int    = input('How many variables of interest? ');

% Load fMRI data
fmrifile='Y';
fmri_name = fmrifile;
eval(['load ' fmrifile '.dat'])
eval(['fmr=' fmrifile ';'])
[num_dat_pts num_runs]=size(fmr);
if (num_runs > 1)
    fmr=mean(fmr')';
    Y=mean(Y')';
end

% Load regressor variables
varfile='X';
eval(['load ' varfile '.dat'])
eval(['varbl=' varfile ';'])
[num_dat_pts num_ind_var]=size(varbl);

% Load regression coefficients
load B_hat.dat
disp(' ')
disp('**************************************************')
disp('Estimated Regression Coefficients:')
disp(' ')
disp(B_hat(num_ind_var-num_var_int+1:num_ind_var)')

disp(' ')
disp('**************************************************')
disp(' ')
disp('Now, plotting the (weighted) independent variables')
disp('in blue and green, and the fMRI signal in red.')
disp(' ')

% Set means to zero, multiply by regression coefficients.
fmr2=fmr-mean(fmr(:));
ymax=max(fmr2(:));
ymin=min(fmr2(:));
for var_index=num_ind_var-num_var_int+1:num_ind_var
    varbl(:,var_index)=varbl(:,var_index)-mean(varbl(:,var_index));
    if(max(varbl(:,var_index)*B_hat(var_index))>ymax)
        ymax=max(varbl(:,var_index)*B_hat(var_index));
    end
    if(min(varbl(:,var_index)*B_hat(var_index))<ymin)
        ymin=min(varbl(:,var_index)*B_hat(var_index));
    end
end
%ymax=max([Y(:); predicted(:)]);
%ymin=min([Y(:); X(:)]);
hold off
axis([0 num_dat_pts ymin ymax])
plot(fmr2,'r')
hold
%plot(fmr2,'r*')

p=zeros(num_dat_pts,1);

for var_index=num_ind_var-num_var_int+1:num_ind_var
    varbl(:,var_index)=varbl(:,var_index)-mean(varbl(:,var_index));
    p=p+(varbl(:,var_index)*B_hat(var_index));
    if (var_index/2==floor(var_index/2))
        plot(varbl(:,var_index)*B_hat(var_index),'b')
%        plot(varbl(:,var_index)*B_hat(var_index),'b*')
	disp(['Variable #' num2str(var_index) ' in blue.'])
    else
        plot(varbl(:,var_index)*B_hat(var_index),'g')
%        plot(varbl(:,var_index)*B_hat(var_index),'g*')
	disp(['Variable #' num2str(var_index) ' in green.'])
    end
%    disp(['varbl = ' num2str(var_index) ', RC = ' num2str(B_hat(var_index))])
%    pause
end
xlabel('Time point')
ylabel('Intensity')
title('fMRI Signals and Independent Variables')
c=B_hat;
c_str='';
for i=num_ind_var-num_var_int+1:num_ind_var
    c_str=[c_str '  ' num2str(c(i))];
end
c_str=['Regression coefficients:' c_str];
text(0.16,0.2,c_str,'sc')
hold off

disp(' ')
disp('Press any key to continue...')
pause

% Plot variables.
disp('**************************************************')
disp('fMRI data is shown in red.')
disp(' ')
disp('Weighted sum of independent variables is shown in blue')
disp('(weights are regression coefficients, of course).')
predicted=X*B_hat;
hold off
plot(Y,'r')
hold
if (num_runs > 1)
    predicted=mean(predicted')';
end
plot(predicted,'b')
hold off
title('Observed and Modeled fMRI Response')
xlabel('Time Point')
ylabel('Intensity')
text(0.16,0.2,c_str,'sc')
