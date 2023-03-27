% plots_modRW.m

% Creation of the plots to be used in the infinite modular random walk
% project. The four cumulants with varying bias at varying block_lengths

% Matthew Gerry, March 2023

%%% PARAMETERS %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, "decoherence rate"

dchi = 0.01; % Counting field step
chisteps = 5; % Number of counting field steps

% Varying parameters
mA_list = [1,2,4,8];
b_list = [0,1/5,1];
dga_axis = 0:0.02:1.9*ga_av;


%%% COMPUTE CGF %%%

% Pre-allocate 4D array that will contain CGF at all param values
bigCGF = zeros([length(mA_list), length(b_list), length(dga_axis), chisteps]);

% Loop through all sets parameter values and get the CGF as a function of
% chi at each
for ii=1:length(dga_axis)
    dga = dga_axis(ii);
    for jj=1:length(mA_list)
        mA = mA_list(jj);
        for kk=1:length(b_list)
            b = b_list(kk);

            % Calculate CGF from rate matrix using functions in this project
            [Lchi,~,~,chi] = diffusionLchi(mA,mA,b,ga_av,dga,tau,dchi,chisteps);
            CGF = CGFclassical(Lchi);

            bigCGF(ii,jj,kk,:) = CGF;

        end % kk
    end % jj
end % ii


%%% COMPUTE CUMULANTS %%%
% Use average of two central values for odd cumulants

% Mean
diff1 = diff(bigCGF, 1, 4); % First derivative of CGF
J = 0.5*(diff1(:,:,:,0.5*(chisteps-1)) + diff1(:,:,:,0.5*(chisteps+1)))/(1i*dchi);

% Variance
diff2 = diff(bigCGF, 2, 4); % Second derivative of CGF
S = diff2(:,:,:,0.5*(chisteps-1))/(1i*dchi)^2;

% Skewness
diff3 = diff(bigCGF, 3, 4); % Third derivative of CGF
C3 = 0.5*(diff3(:,:,:,0.5*(chisteps-3)) + diff3(:,:,:,0.5*(chisteps-1)))/(1i*dchi)^3;

% Kurtosis
diff4 = diff(bigCGF, 4, 4); % Fourth derivative of CGF
C4 = diff4(:,:,:,0.5*(chisteps-3))/(1i*dchi)^4;


%%% PLOT CUMULANTS %%%