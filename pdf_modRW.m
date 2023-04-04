% pdf_modRW.m

% Get the probability distribution function for the number of periods 
% completed (per unit time) by the modular random walk at long time. Vary
% the bias, delta_gamma and block length

% Matthew Gerry, March 2023

%%% PARAMETERS %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, "decoherence rate"

dchi = 0.01; % Counting field step
chisteps = 101; % Number of counting field steps

% Varying parameters
mA_list = [1,2,4,8];
b_list = [0,1/5,8];
dga_list= [0, 0.1, 1.0, 1.9];


%%% COMPUTE CGF %%%

% Pre-allocate 4D array that will contain CGF at all param values
bigCGF = zeros([length(dga_list), length(b_list), length(mA_list), chisteps]);

% Loop through all sets parameter values and get the CGF as a function of
% chi at each
for ii=1:length(dga_list)
    dga = dga_list(ii);
    for jj=1:length(b_list)
        b = b_list(jj);
        for kk=1:length(mA_list)
            mA = mA_list(kk);

            % Calculate CGF from rate matrix using functions in this project
            [Lchi,~,~,chi] = diffusionLchi(mA,mA,b,ga_av,dga,tau,dchi,chisteps);
            CGF = CGFclassical(Lchi);

            bigCGF(ii,jj,kk,:) = CGF;

        end % kk
    end % jj
end % ii


% Get the probability distribution from the CGF
longtime = 1e5;
for n=1:chisteps
    FS_factor = exp(-2i*pi*n*chi/(chi(end)-chi(1)));
    % Need to turn FS factor into a 4d array, multiply it by MGF, sum
    bigMGF = exp(bigCGF*longtime);
    FS_summand = 
    bigPDF = sum(FS_summand);
end



% example
PDF = reshape(bigPDF(3,2,2,:),[1,chisteps]);
plot(PDF)

