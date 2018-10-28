% MATLAB script to calculate effective degrees of freedom
% using method described in Worsley KJ, Friston KJ, Analysis
% of fMRI Time-Series Revisited -- Again, Neuroimage, 2:173-181.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load time series data
load X.dat
[num_dat_pts num_ind_var]=size(X);
G=X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Query user for SD.
SD = input('SD in TR`s? ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Query user for number of variables of interest.
% These will be convolved with hemodynamic response.
% Variables of no interest will not be smoothed.
% This section now commented out.  X.dat now contains
% smoothed data.
% num_var_int=input('Number of variables of interest? ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth time series data, skipping first two independent
% variables, which are assumed to be the constant and ramp
% This section now commented out.  X.dat now contains
% smoothed data.
%G=zeros(num_dat_pts,num_ind_var);
%FHDR=fft(HDR);
%for i=1:num_ind_var-num_var_int-1
%    G(:,i)=X(:,i);
%end
%%%for i=num_ind_var-num_var_int:num_ind_var
%    FX=fft(X(:,i));
%    FG=FX.*FHDR;
%    G(:,i)=real(ifft(FG));
%    G(:,i)=X(:,i);
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct K, the smoothing matrix (Toeplitz matrix).
% The first row is the hemodynamic response associated with
% the first time point.  Each row thereafter is a lagged
% version of the hemodynamic, response, corresponding to
% later time points.
HDR=gauss([0 SD num_dat_pts])';
K=toeplitz(HDR);
% K(1,:)=HDR';
% for row=2:num_dat_pts
%     K(row,:)=[HDR(num_dat_pts-row+2:num_dat_pts)' HDR(1:num_dat_pts-row+1)'];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate V and R, the residual-forming matrix.
V=K*K';
R=eye(num_dat_pts,num_dat_pts)-G*inv(G'*G)*G';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate nu, the effective degrees of freedom.
RV=R*V;
trace_RV=trace(RV);
nu = (trace_RV*trace_RV)/trace(RV*RV)

