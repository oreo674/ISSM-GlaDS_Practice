% Model parameters for sheet-cavity and R-channel models

%% turbulent flow parameters
alpha = 1.25 % exponent 1, see eq. 5 in Werder et al. 2013
beta = 1.5  % exponent 2
k_s = 0.005 % conductivity sheet
k_c = 0.1   % conductivity R-channel
% For a semi-circular R-channel, k_c=0.1 corresponds to a Darcy-Weisbach law with parameter f:
f = 0.195


%% Cavities
u_bed = 1e-6 % ice sliding speed
l_c = 2.0    % width of sheet contributing to R-channel melt
l_r = 2.0    % bedrock bump wavelength
h_r = 0.1    % bedrock bump height

%% Storage
e_v = 0.0  % for sqrt runs
e_v = 1e-3 % for valley runs
