% loglog_modRW.m

% Plot the cumulants of the modular random walks on a log-log scale.
% Plots made here are analogous to those made in plots_modRW.m, but are
% specifically intended to check the scalings of the higher-order cumulants
% with delta_gamma.

% Matthew Gerry, April 2023


%%% PARAMETERS %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, "decoherence rate"

dchi = 0.01; % Counting field step
chisteps = 5; % Number of counting field steps

% Varying parameters
m_list = [1,2,4];
b_list = [0,1/5,8];
dga_axis = logspace(-3,log10(1.99*ga_av));
dga_axis = [0,dga_axis]; % Add zero-value to calculate constant term on its own


%%% COMPUTE CUMULANTS %%%

% Compute CGF
[CGFarray, ~] = bigCGF(tau, ga_av, dchi, chisteps, dga_axis, b_list, m_list);

% Variance
diff2 = diff(CGFarray, 2, 4); % Second derivative of CGF
S = diff2(:,:,:,0.5*(chisteps-1))/(1i*dchi)^2;

% Skewness - use average of two central values for odd cumulants
diff3 = diff(CGFarray, 3, 4); % Third derivative of CGF
C3 = 0.5*(diff3(:,:,:,0.5*(chisteps-3)) + diff3(:,:,:,0.5*(chisteps-1)))/(1i*dchi)^3;

% Kurtosis
diff4 = diff(CGFarray, 4, 4); % Fourth derivative of CGF
C4 = diff4(:,:,:,0.5*(chisteps-3))/(1i*dchi)^4;


%%% PLOT HIGHER-ORDER CUMULANTS %%%

mrkrlist = ['s', 'o', '^', 'x'];
colourlist = ["#0072BD", "#D95319", "#77AC30","#7E2F8E"];
lettlist = ["(a) ", "(b) ", "(c) "];

kstar = tau^2/ga_av; % Homogeneous rate for reference lines


% Variance
figure(1)
for jj=2:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(2,1,jj-1); hold on; box on

    % Prepare analytic curve for comparison (dga-dependent part only)
    b = b_list(jj);
    S_dga_ana = 2*kstar*exp(-b/2)/cosh(b/2)*(0.25*(dga_axis*sinh(b/2)/ga_av).^2);

    plot(dga_axis, S_dga_ana, '--k', DisplayName="Analytic")
    
    for kk=1:length(m_list)
        S_dga = S(:,jj,kk) - S(1, jj, kk); % Subtract off dga-independent part
        plot(dga_axis, S_dga, mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    end % kk
   
    % Format subplot
    xlim([1e-3,max(dga_axis)])
    ylabel("$\langle\langle J^2\rangle\rangle - \langle\langle J^2\rangle\rangle|_{\Delta\gamma=0}$",Interpreter="latex")
    if jj==3
        xlabel("$\Delta\gamma$", Interpreter="latex")
    elseif jj==2
        legend(Interpreter="latex", Location="southeast")
    end % cases

    set(gca, fontsize=14)

    % Label subplot with bias value
    yl = ylim;
    xl = xlim;
    text(10^(-2.8), 1e-1, strcat(lettlist(jj-1),"$b=\;$",num2str(b)), Interpreter="latex", FontSize=14);
end % jj

% Skewness
figure(2)
for jj=2:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(2,1,jj-1); hold on; box on

    % Prepare analytic curve for comparison (dga-dependent part only)
    b = b_list(jj);
    C3_factor1 = (exp(-b/2)*sinh(b/2)/(cosh(b/2)^2));
    C3_factor2 = 0.75*(dga_axis/ga_av).^2 + (3/16)*(dga_axis/ga_av).^4*sinh(b/2)^2;
    C3_dga_ana = 2*(tau^2/ga_av)*C3_factor1*C3_factor2;

    plot(dga_axis, C3_dga_ana, '--k', DisplayName="Analytic")    
    for kk=1:length(m_list)
        C3_dga = C3(:,jj,kk) - C3(1,jj,kk); % Subtract off dga-independent part
        plot(dga_axis, C3_dga, mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    end % kk
   
    % Format subplot
    xlim([1e-3,max(dga_axis)])
    ylabel("$\langle\langle J^3\rangle\rangle - \langle\langle J^3\rangle\rangle|_{\Delta\gamma=0}$",Interpreter="latex")
    if jj==3
        xlabel("$\Delta\gamma$", Interpreter="latex")
    elseif jj==2
        legend(Interpreter="latex", Location="southeast")
    end % cases
    
    set(gca, fontsize=14)

    % Label subplot with bias value
    yl = ylim;
    xl = xlim;
    text(10^(-2.8), 1e-1, strcat(lettlist(jj-1),"$b=\;$",num2str(b)), Interpreter="latex", FontSize=14);
end % jj

% Kurtosis
figure(3)
for jj=1:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(1,3,jj); hold on; box on
    
    % Prepare analytic curve for comparison
    b = b_list(jj);
    C4_factor1 = exp(-b/2)/cosh(b/2)^3;
    C4_factor2 = (1/64)*(dga_axis/ga_av).^2*(exp(-2*b)-36*exp(-b)+118-36*exp(b)+exp(2*b))-(9/256)*(dga_axis/ga_av).^4*(exp(-2*b)-12*exp(-b)+22-12*exp(b)+exp(2*b))+(15/64)*(dga_axis/ga_av).^6*sinh(b/2).^4;
    C4_dga_ana = 2*kstar*C4_factor1*C4_factor2;

    plot(dga_axis, C4_dga_ana, '--k', DisplayName="Analytic")
    for kk=1:length(m_list)
        C4_dga = C4(:,jj,kk) - C4(1,jj,kk); % Subtract off dga-independent part
        plot(dga_axis, C4_dga, mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    end % kk

    % Format subplot
    xlim([10^(-2.5),max(dga_axis)])
    xlabel("$\Delta\gamma$", Interpreter="latex")
    if jj==1
        ylabel("$\langle\langle J^4\rangle\rangle-\langle\langle J^4\rangle\rangle|_{\Delta\gamma=0}$",Interpreter="latex")
        legend(Interpreter="latex", Location="southeast")
    end % case

    set(gca, fontsize=14)

    % Label subplot with bias value
    yl = ylim;
    xl = xlim;
    text(10^(-2.3), 10^(0.1), strcat(lettlist(jj),"$b=\;$",num2str(b)), Interpreter="latex", FontSize=14);

end % jj
