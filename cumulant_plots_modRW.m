% cumulant_plots_modRW.m

% Creation of the plots of the cumulants to be used in the infinite modular
% random walk project. The four cumulants are plotted with varying bias at 
% varying block lengths. First the variance, skewness, and kurtosis are
% plotted as a function of delta gamma, then the mean is plotted as a 
% function of b.

% Matthew Gerry, March 2023


%%% PARAMETERS %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, "decoherence rate"

dchi = 0.01; % Counting field step
chisteps = 5; % Number of counting field steps

% Varying parameters
m_list = [1,2,4,8]; % Block lengths
b_list = [0,1/5,8]; % Biases
b_axis = 0:0.25:5; dga_vs_b = 1.5; % For plot of variance against b
dga_axis = 0:0.1:1.99*ga_av;
dga_logspace = logspace(-3,log10(1.99*ga_av));
% dga_logspace = [0,dga_logspace]; % Add zero-value to calculate constant term on its own

kstar = tau^2/ga_av; % Homogeneous rate for reference lines

% List markerstyles, colours, letters for plot formatting and labels
mrkrlist = ['s', 'o', '^', 'x'];
colourlist = ["#0072BD", "#D95319", "#77AC30","#7E2F8E"];
lettlist = ["(a) ", "(b) ", "(c) "];


%%% COMPUTE CUMULANTS %%%

% Compute CGF
[CGFarray, ~] = bigCGF(tau, ga_av, dchi, chisteps, dga_axis, b_list, m_list);
[CGFloglog, ~] = bigCGF(tau, ga_av, dchi, chisteps, dga_logspace, b_list, m_list);
[CGF_vs_b, ~] = bigCGF(tau, ga_av, dchi, chisteps, dga_vs_b, b_axis, m_list);

% Differentiate the CGF with respect to the counting field (from second up
% to fourth order numerical derivatives)
% Variance
diff2 = diff(CGFarray, 2, 4); % Second derivative of CGF
S = diff2(:,:,:,0.5*(chisteps-1))/(1i*dchi)^2;

diff2log = diff(CGFloglog, 2, 4); % Same thing on a logarithmic scale
Slog = diff2log(:,:,:,0.5*(chisteps-1))/(1i*dchi)^2;

diff2_vs_b = diff(CGF_vs_b, 2, 4); % Same thing with different param ranges
S_vs_b = diff2_vs_b(:,:,:,0.5*(chisteps-1))/(1i*dchi)^2;

% Skewness - use average of two central values since this is an odd cumulant
diff3 = diff(CGFarray, 3, 4); % Third derivative of CGF
C3 = 0.5*(diff3(:,:,:,0.5*(chisteps-3)) + diff3(:,:,:,0.5*(chisteps-1)))/(1i*dchi)^3;

diff3log = diff(CGFloglog, 3, 4); % Third derivative of CGF on a logarithmic scale
C3log = 0.5*(diff3log(:,:,:,0.5*(chisteps-3)) + diff3log(:,:,:,0.5*(chisteps-1)))/(1i*dchi)^3;

% Kurtosis
diff4 = diff(CGFarray, 4, 4); % Fourth derivative of CGF
C4 = diff4(:,:,:,0.5*(chisteps-3))/(1i*dchi)^4;

diff4log = diff(CGFloglog, 4, 4); % Fourth derivative of CGF on a logarithmic scale
C4log = diff4log(:,:,:,0.5*(chisteps-3))/(1i*dchi)^4;


%%% PLOT CUMULANTS AS A FUNCTION OF DELTA GAMMA %%%

% Variance
figure(1)
for jj=2:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(2,1,jj-1); hold on; box on
    
    % Prepare analytic curve for comparison
    b = b_list(jj);
    S_ana = 2*kstar*exp(-b/2)/cosh(b/2)*((cosh(b/2))^2 + 0.25*(dga_axis*sinh(b/2)/ga_av).^2); % Analytically determined value for m=1

    plot(dga_axis, S_ana, '--k', DisplayName="Analytic")
    for kk=1:length(m_list)
        plot(dga_axis, S(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    
    % Reference lines
    if jj==3
        lowline = yline(kstar, ':k', "$\mathcal{C}_{2,\:b\rightarrow\infty, \frac{\Delta\gamma}{\bar{\gamma}}=0}$", Interpreter="latex", FontSize=14);
        lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
        highline = yline(2*kstar, ':k', "$\mathcal{C}_{2,\:b\rightarrow\infty, \frac{\Delta\gamma}{\bar{\gamma}}=2}$", Interpreter="latex", FontSize=14);
        highline.Annotation.LegendInformation.IconDisplayStyle = "off";
    end

    % Format subplot
    yl = ylim;
    ylim([0.95*kstar, 1.15*yl(2)])
    xlim([0,max(dga_axis)])

    ylabel("$\mathcal{C}_2$",Interpreter="latex")
    if jj==3
        xlabel("$\Delta\gamma$", Interpreter="latex")
    end % case

    set(gca, fontsize=14)

    % Label subplot with bias value
    yl = ylim;
    xl = xlim;
    text(0.95*xl(1) + 0.05*xl(2), 0.1*yl(1) + 0.9*yl(2), strcat(lettlist(jj-1),"$b=\;$",num2str(b)), Interpreter="latex", FontSize=14);
    
    hold off

    % Plot on a log-log scale in an inset (panel (a) low bias only)
    if jj==2 % Inset (loglog) and legend
        legend(Interpreter="latex", Location="southeast")
        axes('Position', [.17 .62 .23 .15]); hold on; box on
        for kk=1:length(m_list)
            plot(dga_logspace, Slog(:, jj, kk)-S(1,jj,kk), mrkrlist(kk), Color=colourlist(kk), MarkerSize=8)
        end % kk
        xlim([1e-2,1]);
        text(0.012, 10^(-0.4),"$\mathcal{C}_2-\mathcal{C}_2|_{\Delta\gamma=0}$",Interpreter="latex",FontSize=10)
        set(gca, Fontsize=10)
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        hold off
    end % Inset
    
    % Plot as a function of bias as an inset (panel (b) only)
    if jj==3 % Inset (vs b)
        S_vs_b_ana = 2*kstar*exp(-b_axis/2)./cosh(b_axis/2).*((cosh(b_axis/2)).^2 + 0.25*(dga_vs_b*sinh(b_axis/2)/ga_av).^2);
        axes('Position', [.17 .27 .23 .15]); hold on; box on
        plot(b_axis, S_vs_b_ana, '--k')
        for kk=1:length(m_list)
            plot(b_axis, S_vs_b(1,:,kk), mrkrlist(kk), Color=colourlist(kk), MarkerSize=8)
        end
        xlim([0, b_axis(end)])
        xlabel("$b$", Interpreter="latex")
        set(gca, Fontsize=10)
        hold off
    end % Inset

end % jj

% Skewness
figure(2)
for jj=2:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(2,1,jj-1); hold on; box on
    
    % Prepare analytic curve for comparison - split into factors for readability
    b = b_list(jj);
    C3_factor1 = (exp(-b/2)*sinh(b/2)/(cosh(b/2)^2));
    C3_factor2 = cosh(b/2)^2 + 0.75*(dga_axis/ga_av).^2 + (3/16)*(dga_axis/ga_av).^4*sinh(b/2)^2;
    C3_ana = 2*(tau^2/ga_av)*C3_factor1*C3_factor2;  % Analytically determined value for m=1

    % Plot analytic and numeric curves
    plot(dga_axis, C3_ana, '--k', DisplayName="Analytic")
    for kk=1:length(m_list)
        plot(dga_axis, C3(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    
    % Reference lines (limit values of the skewness at high bias)
    if jj==3
        lowline = yline(kstar, ':k', "$\mathcal{C}_{3,\:b\rightarrow\infty, \frac{\Delta\gamma}{\bar{\gamma}}=0}$", Interpreter="latex", FontSize=14);
        lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
        highline = yline(4*kstar, ':k', "$\mathcal{C}_{3,\:b\rightarrow\infty, \frac{\Delta\gamma}{\bar{\gamma}}=2}$", Interpreter="latex", FontSize=14);
        highline.Annotation.LegendInformation.IconDisplayStyle = "off";
    end

    % Format subplot
    yl = ylim;
    ylim([0, 1.15*yl(2)])
    xlim([0,max(dga_axis)])

    ylabel("$\mathcal{C}_3$",Interpreter="latex")
    if jj==3
        xlabel("$\Delta\gamma$", Interpreter="latex")
        legend(Interpreter="latex", Location="southwest")
    end % case
    set(gca, fontsize=14)

    % Label subplot with bias value
    yl = ylim;
    xl = xlim;
    text(0.95*xl(1) + 0.05*xl(2), 0.1*yl(1) + 0.9*yl(2), strcat(lettlist(jj-1),"$b=\;$",num2str(b)), Interpreter="latex", FontSize=14);
    
    % Plot on a log-log scale in an inset (panel (a) low bias only)
    if jj==2 % Inset
        axes('Position', [.17 .62 .23 .15]); hold on; box on
        for kk=1:length(m_list)
            plot(dga_logspace, C3log(:, jj, kk)-C3(1,jj,kk), mrkrlist(kk), Color=colourlist(kk), MarkerSize=4)
        end % kk
        xlim([1e-2,1]);
        text(0.012, 10^(-0.4),"$\mathcal{C}_3-\mathcal{C}_3|_{\Delta\gamma=0}$",Interpreter="latex",FontSize=10.5)
        set(gca, Fontsize=10)
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    end % Inset
end % jj

% Kurtosis
figure(3)
for jj=1:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(1,3,jj); hold on; box on
    
    % Prepare analytic curve for comparison - split into factors for readability
    b = b_list(jj);
    C4_factor1 = exp(-b/2)/cosh(b/2)^3;
    C4_factor2 = cosh(b/2)^4 + (1/64)*(dga_axis/ga_av).^2*(exp(-2*b)-36*exp(-b)+118-36*exp(b)+exp(2*b))-(9/256)*(dga_axis/ga_av).^4*(exp(-2*b)-12*exp(-b)+22-12*exp(b)+exp(2*b))+(15/64)*(dga_axis/ga_av).^6*sinh(b/2).^4;
    C4_ana = 2*kstar*C4_factor1*C4_factor2; % Analytically determined values for m=1

    plot(dga_axis, C4_ana, '--k', DisplayName="Analytic")
    for kk=1:length(m_list)
        plot(dga_axis, C4(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    
    % Reference lines (limit values of the kurtosis at high bias)
    if jj==3
        lowline = yline(kstar, ':k', "$\mathcal{C}_{4,\:b\rightarrow\infty, \frac{\Delta\gamma}{\bar{\gamma}}=0}$", Interpreter="latex", FontSize=14);
        lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
        highline = yline(8*kstar, ':k', "$\mathcal{C}_{4,\:b\rightarrow\infty, \frac{\Delta\gamma}{\bar{\gamma}}=2}$", Interpreter="latex", FontSize=14);
        highline.Annotation.LegendInformation.IconDisplayStyle = "off";
    end

    % Format subplot
    xlim([0,max(dga_axis)])
    xlabel("$\Delta\gamma$", Interpreter="latex")
    if jj==1
        ylabel("$\mathcal{C}_4$",Interpreter="latex")
        lowline.LabelHorizontalAlignment = "Left";
        highline.LabelHorizontalAlignment = "Left";
        ylim([0, 80])
    elseif jj==2
        lowline.LabelHorizontalAlignment = "Left";
        highline.LabelHorizontalAlignment = "Left";
    elseif jj==3
        legend(Interpreter="latex", Location="southeast")
        ylim([0, 9.5*kstar])
    end % case

    set(gca, fontsize=14)

    % Label subplot with bias value
    yl = ylim;
    xl = xlim;
    text(0.95*xl(1) + 0.05*xl(2), 0.08*yl(1) + 0.92*yl(2), strcat(lettlist(jj),"$b=\;$",num2str(b)), Interpreter="latex", FontSize=14);

    % Plot on a log-log scale in an inset (panel (a) low bias only)
    if jj==1 % Inset
        axes('Position', [.17 .62 .23 .15]); hold on; box on
        for kk=1:length(m_list)
            plot(dga_logspace, C4log(:, jj, kk)-C4(1,jj,kk), mrkrlist(kk), Color=colourlist(kk), MarkerSize=4)
        end % kk
        xlim([1e-2,1]);
        text(0.012, 10^(-0.4),"$\mathcal{C}_4-\mathcal{C}_4|_{\Delta\gamma=0}$",Interpreter="latex",FontSize=10)
        set(gca, Fontsize=10)
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    end % Inset

end % jj

%%% PLOT MEAN AS A FUNCTION OF b %%%

% Choose one delta_gamma value and use more values in the b_list just to
% make this plot
b_list2 = 0:0.2:4; % For plots against b
dga_axis2 = [0, 0.25, 1.99]; % For plots against b

% Compute CGF
[CGFarray_b, ~] = bigCGF(tau, ga_av, dchi, chisteps, dga_axis2, b_list2, m_list);

% Calculate mean numerically - use average of two central values
diff1 = diff(CGFarray_b, 1, 4); % First derivative of CGF
J = 0.5*(diff1(:,:,:,0.5*(chisteps-1)) + diff1(:,:,:,0.5*(chisteps+1)))/(1i*dchi);

% Analytic expresssion for mean for comparison
b_high_res = 0:0.05:4;
J_ana = 2*exp(-0.5*b_high_res).*sinh(0.5*b_high_res)*kstar;

figure(4)
hold on; box on
plot(b_high_res, J_ana, '--k', DisplayName="Analytic")
for kk=1:length(m_list)
    plot(b_list2, J(3,:,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
end % kk

% Reference line
lowline = yline(kstar, ':k', "$k^*$", Interpreter="latex", FontSize=14);
lowline.Annotation.LegendInformation.IconDisplayStyle = "off";

xlabel("$b$",Interpreter="latex")
ylabel("$\mathcal{C}_1$",Interpreter="latex")
legend(Location="southeast", Interpreter="latex")
set(gca, fontsize=14)
ylim([0,1.12*kstar])

%% % PLOT HIGHER-ORDER CUMULANTS AS A FUNCTION OF b %%%
% Plots not included in final version of the project

% Get higher order cumulants
% Variance
diff2_b = diff(CGFarray_b, 2, 4); % Second derivative of CGF
S_b = diff2_b(:,:,:,0.5*(chisteps-1))/(1i*dchi)^2;

% Skewness - use average of two central values for odd cumulants
diff3_b = diff(CGFarray_b, 3, 4); % Third derivative of CGF
C3_b = 0.5*(diff3_b(:,:,:,0.5*(chisteps-3)) + diff3_b(:,:,:,0.5*(chisteps-1)))/(1i*dchi)^3;

% Kurtosis
diff4_b = diff(CGFarray_b, 4, 4); % Fourth derivative of CGF
C4_b = diff4_b(:,:,:,0.5*(chisteps-3))/(1i*dchi)^4;

% Plot variance
figure(5)
for jj=1:length(dga_axis2)
    subplot(1,length(dga_axis2),jj); hold on; box on
    for kk=1:length(m_list)
        plot(b_list2, S_b(jj,:,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    xlabel("$b$",Interpreter="latex")
    if jj==1
        ylabel("$\mathcal{C}_2$",Interpreter="latex")
        legend(Location="southeast", Interpreter="latex")
    end % case

    % Label subplot with delta gamma value
    yl = ylim;
    xl = xlim;
    text(0.95*xl(1) + 0.05*xl(2), 0.08*yl(1) + 0.92*yl(2), strcat(lettlist(jj),"$\Delta\gamma=\;$",num2str(dga_axis2(jj))), Interpreter="latex", FontSize=14);

    set(gca, fontsize=14)
    hold off
end % jj


% Plot skewness
figure(6)
for jj=1:length(dga_axis2)
    subplot(1,length(dga_axis2),jj); hold on; box on
    for kk=1:length(m_list)
        plot(b_list2, C3_b(jj,:,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    xlabel("$b$",Interpreter="latex")
    if jj==1
        ylabel("$\mathcal{C}_3$",Interpreter="latex")
        legend(Location="southeast", Interpreter="latex")
    end % case

    % Label subplot with delta gamma value
    yl = ylim;
    xl = xlim;
    text(0.95*xl(1) + 0.05*xl(2), 0.08*yl(1) + 0.92*yl(2), strcat(lettlist(jj),"$\Delta\gamma=\;$",num2str(dga_axis2(jj))), Interpreter="latex", FontSize=14);

    set(gca, fontsize=14)
    hold off
end % jj


% Plot kurtosis
figure(7)
for jj=1:length(dga_axis2)
    subplot(1,length(dga_axis2),jj); hold on; box on
    for kk=1:length(m_list)
        plot(b_list2, C4_b(jj,:,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    xlabel("$b$",Interpreter="latex")
    if jj==1
        ylabel("$\mathcal{C}_4$",Interpreter="latex")
        legend(Location="southeast", Interpreter="latex")
    end % case

    % Label subplot with delta gamma value
    yl = ylim;
    xl = xlim;
    text(0.95*xl(1) + 0.05*xl(2), 0.08*yl(1) + 0.92*yl(2), strcat(lettlist(jj),"$\Delta\gamma=\;$",num2str(dga_axis2(jj))), Interpreter="latex", FontSize=14);

    set(gca, fontsize=14)
    hold off
end % jj


% Plot kurtosis over skewness to investigate ratios of cumulants
figure(8)
for jj=1:length(dga_axis2)
    subplot(1,length(dga_axis2),jj); hold on; box on

    % Analytic expression for small dga/ga_av only
%     ratio_ana = coth(b_list2/2).*(1 + (dga_axis2(jj)/ga_av).^2*(((1/64)*(exp(-2*b_list2)-36*exp(-b_list2)+118-16*exp(b_list2)+exp(2*b_list2)))./(64*cosh(b_list2/2).^3) - 3./(4*cosh(b_list2/2).^2)));
%     plot(b_list2, ratio_ana, '--k')
    
    for kk=1:length(m_list)
        plot(b_list2, C4_b(jj,:,kk)./C3_b(jj,:,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    xlabel("$b$",Interpreter="latex")
    if jj==1
        ylabel("$\mathcal{C}_4/\mathcal{C}_3$",Interpreter="latex")
        legend(Location="southeast", Interpreter="latex")
    end % case

    % Label subplot with delta gamma value
    yl = ylim;
    xl = xlim;
    text(0.95*xl(1) + 0.05*xl(2), 0.08*yl(1) + 0.92*yl(2), strcat(lettlist(jj),"$\Delta\gamma=\;$",num2str(dga_axis2(jj))), Interpreter="latex", FontSize=14);

    set(gca, fontsize=14)
    hold off
end % jj