% pdf_modRW.m

% Get the probability distribution function for the number of periods 
% completed (per unit time) by the modular random walk at long time. Vary
% the bias, delta_gamma and block length

% Matthew Gerry, April 2023

%%% PARAMETERS %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, "decoherence rate"


chisteps = 101; % Number of counting field steps
dchi = 2*pi/chisteps; % Counting field step - ensure chi runs from -pi to pi

% Varying parameters
mA_list = [1,2,4,8];
b_list = [0,1/5,8];
dga_list= [0, 0.1, 1.0, 1.9];


%%% COMPUTE CGF %%%

% Pre-allocate 4D array that will contain CGF at all param values
bigCGF = zeros([length(dga_list), length(b_list), length(mA_list), chisteps, chisteps]);

% Loop through all sets parameter values and get the CGF as a function of
% chi at each
for ii=1:length(dga_list)
    dga = dga_list(ii);
    for jj=1:length(b_list)
        b = b_list(jj);
        for kk=1:length(mA_list)
            mA = mA_list(kk);

            % Calculate CGF from rate matrix using functions in this project
            [Lchi,~,~,chiA,chiB] = Lchi_detailed(mA,mA,b,ga_av,dga,tau,dchi,chisteps);
            CGF = CGF_detailed(Lchi);

            % Big CGF has substantial extention along the chiA and chiB
            % dimensions, so we can Fourier transform it with respect to
            % these chis to get the probability distribution function
            bigCGF(ii,jj,kk,:,:) = CGF;

        end % kk
    end % jj
end % ii
