% pdf_modRW.m

% Get the probability distribution function for the number of periods 
% completed (per unit time) by the modular random walk at long time. Vary
% the bias, delta_gamma and block length

% Matthew Gerry, April 2023

%%% PARAMETERS %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, "decoherence rate"


chisteps = 151; % Number of counting field steps
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

% Get the probability distribution from the CGF
% Pre-allocate 5-index array for the PDF
bigPDF = zeros([length(dga_list), length(b_list), length(mA_list), chisteps, chisteps]); % Pre-allocate

% Need to approximate moment generating function at longtime from CGF
longtime = 1e2; % This should be long enough in s for the system to reach steady state
bigMGF = exp(bigCGF*longtime);

% Get coefficients associated with Fourier transform of MGF
for nA=-0.5*(chisteps-1):0.5*(chisteps-1)
    for nB=-0.5*(chisteps-1):0.5*(chisteps-1)
        % Define and reshape the complex exponential as needed
        FS_factor = exp(-1i*(nA*chiA + nB*chiB))/(2*pi)^2;
        FS_factor = reshape(FS_factor, [1,1,1,chisteps,chisteps]);
        FS_factor = repmat(FS_factor,[length(dga_list),length(b_list),length(mA_list),1,1]);
        
        % Get Fourier coefficients, build up bigPDF
        FS_summand = bigMGF.*FS_factor;
        P_AB = sum(sum(FS_summand,5),4)*dchi^2;
        bigPDF(:,:,:,nA+0.5*(chisteps+1),nB+0.5*(chisteps+1)) = P_AB;
    end
end

bigPDF = real(bigPDF); % Eliminate imaginary parts (nonphysical, due to simulation noise)

n_vals = -0.5*(chisteps-1):0.5*(chisteps-1); % nA and nB values corresponding to actual numbers of steps
n_tot_vals = -chisteps+1:chisteps-1; % Possible values of the total number of steps

% From the joint PDF for nA+nB, get the PDF for just n by tracing out the
% distinction
bigPDFn =  zeros([length(dga_list), length(b_list), length(mA_list), 2*chisteps-1]); % Pre-allocate
for ii=1:chisteps
    for jj=1:chisteps
        bigPDFn(:,:,:,ii+jj-1) = bigPDFn(:,:,:,ii+jj-1) + bigPDF(:,:,:,ii,jj);
    end
end
%%
plot(n_tot_vals, reshape(bigPDFn(3,2,3,:),[1,2*chisteps-1]))
