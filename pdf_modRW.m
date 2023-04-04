% pdf_modRW.m

% Get the probability distribution function for the number of periods 
% completed (per unit time) by the modular random walk at long time. Vary
% the bias, delta_gamma and block length

% Matthew Gerry, March 2023

%%% PARAMETERS %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, "decoherence rate"


chisteps = 301; % Number of counting field steps
dchi = 2*pi/chisteps; % Counting field step - ensure chi runs from -pi to pi

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

            % Big CGF has substantial extention along the chi dimension, so
            % we can Fourier transform it with respect to chi to get the
            % CGF
            bigCGF(ii,jj,kk,:) = CGF;

        end % kk
    end % jj
end % ii


% Get the probability distribution from the CGF
bigPDF = zeros([length(dga_list), length(b_list), length(mA_list), chisteps]);
longtime = 1e2;
for n=-0.5*(chisteps-1):0.5*(chisteps-1)
    FS_factor = exp(-1i*n*chi)/(2*pi);
    FS_factor = reshape(FS_factor, [1,1,1,chisteps]);
    FS_factor = repmat(FS_factor,[length(dga_list),length(b_list),length(mA_list),1]);
    
    bigMGF = exp(bigCGF*longtime);
    FS_summand = bigMGF.*FS_factor;
    Pn = sum(FS_summand,4)*dchi;
    bigPDF(:,:,:,n+0.5*(chisteps+1)) = Pn;
end

bigPDF = real(bigPDF); % Remove remaining imag parts (just numerical error, order 1e-17)

% example
PDF_raw = reshape(bigPDF(1,2,2,:),[1,chisteps]);
PDF = PDF_raw(1:2:end);
n_vals = -0.5*(chisteps-1):2:0.5*(chisteps-1);
plot(n_vals, PDF)

