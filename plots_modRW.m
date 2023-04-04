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
b_list = [0,1/5,8];
dga_axis = 0:0.1:1.99*ga_av;


%%% COMPUTE CGF %%%

% Pre-allocate 4D array that will contain CGF at all param values
bigCGF = zeros([length(dga_axis), length(b_list), length(mA_list), chisteps]);

% Loop through all sets parameter values and get the CGF as a function of
% chi at each
for ii=1:length(dga_axis)
    dga = dga_axis(ii);
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
bigPDF = fft(exp(bigCGF*longtime), chisteps, 4);



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

mrkrlist = ['s', 'o', '^', 'x'];
colourlist = ["#0072BD", "#D95319", "#77AC30","#7E2F8E"];

ktilde = tau^2/ga_av; % Homogeneous rate for reference lines

% % Mean
% figure(1)
% for jj=2:length(b_list) % Exclude zero-bias case
%     subplot(1,2,jj-1); hold on; box on
%     for kk=1:length(mA_list)
%         plot(dga_axis, J(:,jj,kk))
%     end % kk
%     ylim([0,1.2*tau^2/ga_av])
% end % jj

% Variance
figure(2)
for jj=2:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(1,2,jj-1); hold on; box on
    
    % Prepare analytic curve for comparison
    b = b_list(jj);
    S_ana = 2*ktilde*exp(-b/2)/cosh(b/2)*((cosh(b/2))^2 + 0.25*(dga_axis*sinh(b/2)/ga_av).^2);

    for kk=1:length(mA_list)
        plot(dga_axis, S(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m_A =\;$",num2str(mA_list(kk))))
    end % kk
    plot(dga_axis, S_ana, '--k', DisplayName="Analytic")
    
    % Reference lines
    lowline = yline(ktilde, ':k', "$\tilde{k}$", Interpreter="latex", FontSize=14);
    lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
    highline = yline(2*ktilde, ':k', "$2\tilde{k}$", Interpreter="latex", FontSize=14);
    highline.Annotation.LegendInformation.IconDisplayStyle = "off";

    % Format subplot
    yl = ylim;
    ylim([0.95*ktilde, 1.15*yl(2)])
    xlim([0,max(dga_axis)])

    xlabel("$\Delta\gamma$", Interpreter="latex")
    if jj==2
        ylabel("$\langle\langle J^2\rangle\rangle$",Interpreter="latex")
        legend(Interpreter="latex", Location="southeast")
        lowline.LabelHorizontalAlignment = "Left";
        highline.LabelHorizontalAlignment = "Left";
    end % case
    set(gca, fontsize=14)
end % jj

% Skenwness
figure(3)
for jj=2:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(1,2,jj-1); hold on; box on
    
    % Prepare analytic curve for comparison
    b = b_list(jj);
    C3_factor1 = (exp(-b/2)*sinh(b/2)/(cosh(b/2)^2));
    C3_factor2 = cosh(b/2)^2 + 0.75*(dga_axis/ga_av).^2 + (3/16)*(dga_axis/ga_av).^4*sinh(b/2)^2;
    C3_ana = 2*(tau^2/ga_av)*C3_factor1*C3_factor2;

    for kk=1:length(mA_list)
        plot(dga_axis, C3(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m_A =\;$",num2str(mA_list(kk))))
    end % kk
    plot(dga_axis, C3_ana, '--k', DisplayName="Analytic")
    
    % Reference lines
    lowline = yline(ktilde, ':k', "$\tilde{k}$", Interpreter="latex", FontSize=14);
    lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
    highline = yline(4*ktilde, ':k', "$4\tilde{k}$", Interpreter="latex", FontSize=14);
    highline.Annotation.LegendInformation.IconDisplayStyle = "off";

    % Format subplot
    yl = ylim;
    ylim([0, 1.15*yl(2)])
    xlim([0,max(dga_axis)])

    xlabel("$\Delta\gamma$", Interpreter="latex")
    if jj==2
        ylabel("$\langle\langle J^3\rangle\rangle$",Interpreter="latex")
        legend(Interpreter="latex", Location="northwest")
        lowline.LabelHorizontalAlignment = "Left";
        highline.LabelHorizontalAlignment = "Left";
    end % case
    set(gca, fontsize=14)
end % jj

%s Kurtosis
figure(4)
for jj=1:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(1,3,jj); hold on; box on
    
    % Prepare analytic curve for comparison
    b = b_list(jj);
    C4_factor1 = exp(-b/2)/cosh(b/2)^3;
    C4_factor2 = cosh(b/2)^4 + (1/64)*(dga_axis/ga_av).^2*(exp(-2*b)-36*exp(-b)+118-36*exp(b)+exp(2*b))-(9/256)*(dga_axis/ga_av).^4*(exp(-2*b)-12*exp(-b)+22-12*exp(b)+exp(2*b))+(15/64)*(dga_axis/ga_av).^6*sinh(b/2).^4;
    C4_ana = 2*ktilde*C4_factor1*C4_factor2;

    for kk=1:length(mA_list)
        plot(dga_axis, C4(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m_A =\;$",num2str(mA_list(kk))))
    end % kk
    plot(dga_axis, C4_ana, '--k', DisplayName="Analytic")
    
    % Reference lines
    lowline = yline(ktilde, ':k', "$\tilde{k}$", Interpreter="latex", FontSize=14);
    lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
    highline = yline(8*ktilde, ':k', "$8\tilde{k}$", Interpreter="latex", FontSize=14);
    highline.Annotation.LegendInformation.IconDisplayStyle = "off";

    % Format subplot
    xlim([0,max(dga_axis)])
    xlabel("$\Delta\gamma$", Interpreter="latex")
    if jj==1
        ylabel("$\langle\langle J^4\rangle\rangle$",Interpreter="latex")
        legend(Interpreter="latex", Location="northwest")
        lowline.LabelHorizontalAlignment = "Left";
        highline.LabelHorizontalAlignment = "Left";
        ylim([0, 80])
    elseif jj==3
        ylim([0, 8.8*ktilde])
    end % case

    set(gca, fontsize=14)
end % jj