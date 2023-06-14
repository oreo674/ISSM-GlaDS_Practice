% Physical parameters, all in SI units

% Variable names as in Werder et al. 2013
rho_w = 1000.0 % density water
rho_i = 910.0  % density ice (or more accurately density of whole glacier)
g_grav = 9.8000 % acceleration due to gravity
L_fusion = 334000.0 % latent heat of fusion
c_w = 4220.0 % specific heat capacity water
c_t = 7.5e-08 % Clausius-Clapeyron constant
n_glen = 3.0 % Glen's n
A = 2.5e-25 % Ice flow constant for a relation of the form A*S*N^n_glen
            % with S channel x-sectional area & N effective pressure

% Time in sec
day = 24*60*60;  % seconds per day
dty = 365.0      % days per year
year = dty*day;  % seconds per year
