% plots_modRW.m

% Creation of the plots to be used in the infinite modular random walk
% project. The four cumulants with varying bias at varying block_lengths.
% First the variance, skewness, and kurtosis are plotted as a function of
% delta gamma, then the mean is plotted as a function of b.

% Matthew Gerry, March 2023


%%% PARAMETERS %%%

% Constant parameters
tau = 1.0; % 1/s, "tunnel coupling"
ga_av = 1.0; % 1/s, "decoherence rate"

dchi = 0.01; % Counting field step
chisteps = 5; % Number of counting field steps

% Varying parameters
m_list = [1,2,4,8];
b_list = [0,1/5,8];
dga_axis = 0:0.1:1.99*ga_av;
dga_logspace = logspace(-3,log10(1.99*ga_av));
% dga_logspace = [0,dga_logspace]; % Add zero-value to calculate constant term on its own


%%% COMPUTE CUMULANTS %%%

% Compute CGF
[CGFarray, ~] = bigCGF(tau, ga_av, dchi, chisteps, dga_axis, b_list, m_list);
[CGFloglog, ~] = bigCGF(tau, ga_av, dchi, chisteps, dga_logspace, b_list, m_list);


% Variance
diff2 = diff(CGFarray, 2, 4); % Second derivative of CGF
S = diff2(:,:,:,0.5*(chisteps-1))/(1i*dchi)^2;

diff2log = diff(CGFloglog, 2, 4); % Same thing on logarithmic scale
Slog = diff2log(:,:,:,0.5*(chisteps-1))/(1i*dchi)^2;

% Skewness - use average of two central values for odd cumulants
diff3 = diff(CGFarray, 3, 4); % Third derivative of CGF
C3 = 0.5*(diff3(:,:,:,0.5*(chisteps-3)) + diff3(:,:,:,0.5*(chisteps-1)))/(1i*dchi)^3;

diff3log = diff(CGFloglog, 3, 4); % Third derivative of CGF
C3log = 0.5*(diff3log(:,:,:,0.5*(chisteps-3)) + diff3log(:,:,:,0.5*(chisteps-1)))/(1i*dchi)^3;

% Kurtosis
diff4 = diff(CGFarray, 4, 4); % Fourth derivative of CGF
C4 = diff4(:,:,:,0.5*(chisteps-3))/(1i*dchi)^4;

diff4log = diff(CGFloglog, 4, 4); % Fourth derivative of CGF
C4log = diff4log(:,:,:,0.5*(chisteps-3))/(1i*dchi)^4;


%%% PLOT CUMULANTS %%%

mrkrlist = ['s', 'o', '^', 'x'];
colourlist = ["#0072BD", "#D95319", "#77AC30","#7E2F8E"];
lettlist = ["(a) ", "(b) ", "(c) "];

kstar = tau^2/ga_av; % Homogeneous rate for reference lines


% Variance
figure(1)
for jj=2:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(2,1,jj-1); hold on; box on
    
    % Prepare analytic curve for comparison
    b = b_list(jj);
    S_ana = 2*kstar*exp(-b/2)/cosh(b/2)*((cosh(b/2))^2 + 0.25*(dga_axis*sinh(b/2)/ga_av).^2);

    plot(dga_axis, S_ana, '--k', DisplayName="Analytic")
    for kk=1:length(m_list)
        plot(dga_axis, S(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    
    % Reference lines
    lowline = yline(kstar, ':k', "$k^*$         ", Interpreter="latex", FontSize=14);
    lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
    highline = yline(2*kstar, ':k', "$2k^*$         ", Interpreter="latex", FontSize=14);
    highline.Annotation.LegendInformation.IconDisplayStyle = "off";

    % Format subplot
    yl = ylim;
    ylim([0.95*kstar, 1.15*yl(2)])
    xlim([0,max(dga_axis)])

    ylabel("$\mathcal{C}_2$",Interpreter="latex")
    if jj==3
        xlabel("$\Delta\gamma$", Interpreter="latex")
%         lowline.LabelHorizontalAlignment = "Left";
%         highline.LabelHorizontalAlignment = "Left";
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
            plot(dga_logspace, Slog(:, jj, kk)-S(1,jj,kk), mrkrlist(kk), Color=colourlist(kk), MarkerSize=10)
        end % kk
        xlim([1e-2,1]);
        text(0.012, 10^(-0.4),"$\mathcal{C}_2-\mathcal{C}_2|_{\Delta\gamma=0}$",Interpreter="latex",FontSize=10)
        set(gca, Fontsize=10)
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
    end % Inset

end % jj

% Skenwness
figure(2)
for jj=2:length(b_list) % Exclude zero bias case - no block length or dga dependence
    subplot(2,1,jj-1); hold on; box on
    
    % Prepare analytic curve for comparison
    b = b_list(jj);
    C3_factor1 = (exp(-b/2)*sinh(b/2)/(cosh(b/2)^2));
    C3_factor2 = cosh(b/2)^2 + 0.75*(dga_axis/ga_av).^2 + (3/16)*(dga_axis/ga_av).^4*sinh(b/2)^2;
    C3_ana = 2*(tau^2/ga_av)*C3_factor1*C3_factor2;

    plot(dga_axis, C3_ana, '--k', DisplayName="Analytic")
    for kk=1:length(m_list)
        plot(dga_axis, C3(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    
    % Reference lines
    lowline = yline(kstar, ':k', "$k^*$", Interpreter="latex", FontSize=14);
    lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
    highline = yline(4*kstar, ':k', "$4k^*$", Interpreter="latex", FontSize=14);
    highline.Annotation.LegendInformation.IconDisplayStyle = "off";

    % Format subplot
    yl = ylim;
    ylim([0, 1.15*yl(2)])
    xlim([0,max(dga_axis)])

    ylabel("$\mathcal{C}_3$",Interpreter="latex")
    if jj==3
        xlabel("$\Delta\gamma$", Interpreter="latex")
%         lowline.LabelHorizontalAlignment = "Left";
%         highline.LabelHorizontalAlignment = "Left";
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
    
    % Prepare analytic curve for comparison
    b = b_list(jj);
    C4_factor1 = exp(-b/2)/cosh(b/2)^3;
    C4_factor2 = cosh(b/2)^4 + (1/64)*(dga_axis/ga_av).^2*(exp(-2*b)-36*exp(-b)+118-36*exp(b)+exp(2*b))-(9/256)*(dga_axis/ga_av).^4*(exp(-2*b)-12*exp(-b)+22-12*exp(b)+exp(2*b))+(15/64)*(dga_axis/ga_av).^6*sinh(b/2).^4;
    C4_ana = 2*kstar*C4_factor1*C4_factor2;

    plot(dga_axis, C4_ana, '--k', DisplayName="Analytic")
    for kk=1:length(m_list)
        plot(dga_axis, C4(:,jj,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
    end % kk
    
    % Reference lines
    lowline = yline(kstar, ':k', "$k^*$", Interpreter="latex", FontSize=14);
    lowline.Annotation.LegendInformation.IconDisplayStyle = "off";
    highline = yline(8*kstar, ':k', "$8k^*$", Interpreter="latex", FontSize=14);
    highline.Annotation.LegendInformation.IconDisplayStyle = "off";

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

%% Plot mean current as a function of b

% Choose one delta_gamma value and use more values in the b_list just to
% make this plot

b_list_mean = 0:0.2:4; % For mean vs b plot
dga_axis_mean = 1.0; % For mean vs b plot

% Compute CGF
[CGFarray_mean, ~] = bigCGF(tau, ga_av, dchi, chisteps, dga_axis_mean, b_list_mean, m_list);

% Calculate mean numerically - use average of two central values
diff1 = diff(CGFarray_mean, 1, 4); % First derivative of CGF
J = 0.5*(diff1(:,:,:,0.5*(chisteps-1)) + diff1(:,:,:,0.5*(chisteps+1)))/(1i*dchi);

% Analytic expresssion for comparison
b_high_res = 0:0.05:4;
J_ana = 2*exp(-0.5*b_high_res).*sinh(0.5*b_high_res)*kstar;

figure(4)
hold on; box on
plot(b_high_res, J_ana, '--k', DisplayName="Analytic")
for kk=1:length(m_list)
    plot(b_list_mean, J(1,:,kk), mrkrlist(kk), Color=colourlist(kk), DisplayName=strcat("$m =\;$",num2str(m_list(kk))))
end % kk

% Reference line
lowline = yline(kstar, ':k', "$k^*$", Interpreter="latex", FontSize=14);
lowline.Annotation.LegendInformation.IconDisplayStyle = "off";

xlabel("$b$",Interpreter="latex")
ylabel("$\mathcal{C}_1$",Interpreter="latex")
legend(Location="southeast", Interpreter="latex")
set(gca, fontsize=14)
ylim([0,1.12*kstar])


