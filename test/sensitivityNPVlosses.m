%% Input

options.General.timeHorizon = 50;
options.General.intRate = 0;

% Concrete moment resisting frame (1970s)
options.Vulnerability.fragMedians = [0.166,0.32,0.34,0.46];
options.Vulnerability.fragDispersions = [0.437,0.442,0.443,0.447];
options.Vulnerability.damgeLossRatios = [0 2 10 43.5 95]/100;
options.Vulnerability.CoVdlr = [ 0 1 0.4 0.3 0.05 ];

% Seismic hazard from L'Aquila, Italy
options.Hazard.faultRate = 0.08;
options.Hazard.hazCurve = [0.166,0.0332; 0.217,0.0199; 0.258,0.014; 0.304,0.0099; 0.354,0.007; 0.412,0.005; 0.565,0.002; 0.695,0.001; 0.846,0.0004];

%% Sensitivity of the analytical formulation to the definition of the loss vector

Nsteps = 1+[10, 100, 200, 500, 1000, 2000];

for s = numel(Nsteps) : -1 : 1
    disp(s)
    tic
    options.Setup.NlossSamples = Nsteps(s);
    sensStep = distNPVaggregateLosses(options);
    
    sensStep = sensStep.getPDFinterarrivalTime;
    sensStep = sensStep.getPMFnumberEvents;
    sensStep = sensStep.getPDFarrivalTime;
    sensStep = sensStep.getLossDistribution;
    sensStep = sensStep.getPDFlossNPV;
    sensStep = sensStep.getAggregateLossNPVdist;
    toc
    
    PDFaggLossNPV{s} = sensStep.PDFaggUninsuredNPV;
end

%% Sensitivity of the Montecarlo to the Nsamples

Nsamples = [100 1000 10000 100000 500000 1000000];

options.Setup.NlossSamples = 100000;
sensMC = distNPVaggregateLosses(options);
sensMC = sensMC.getPDFinterarrivalTime;
sensMC = sensMC.getPMFnumberEvents;
sensMC = sensMC.getPDFarrivalTime;
sensMC = sensMC.getLossDistribution;
sensMC.parameters.Setup.MCsamples = Nsamples(end);
tic; sensMC = sensMC.monteCarloPDFaggregateLossNPV; toc

%% Sensitivity to the interest rate

intRates = [0 0.005 0.01 0.02 0.04 0.08];

options.Setup.NlossSamples = 301;
options.Setup.MCsamples = 20000;
sensInt = distNPVaggregateLosses(options);
sensInt = sensInt.getPDFinterarrivalTime;
sensInt = sensInt.getPMFnumberEvents;
sensInt = sensInt.getPDFarrivalTime;
sensInt = sensInt.getLossDistribution;

for r = numel(intRates) : -1 : 1
    sensInt.parameters.General.intRate = intRates(r);
    
    sensInt = sensInt.getPDFlossNPV;
    sensInt = sensInt.getAggregateLossNPVdist;
    sensInt = sensInt.monteCarloPDFaggregateLossNPV;
    
    sensIntAggLossNPV{r,1} = sensInt.PDFaggUninsuredNPV;
    sensIntAggLossNPV{r,2} = sensInt.NPVaggUninsuredMC;
end

%% Plot

colSteps = gray(numel(Nsteps)+1);
colInt = hsv(numel(intRates));

leg = cell(1, numel(Nsteps)+1);
leg{1} = 'MonteCarlo';
for s = 1 : numel(Nsteps); leg{s+1} = num2str(Nsteps(s)); end

figure; hold on
histogram(sensMC.NPVaggUninsuredMC, 'Normalization', 'pdf');
for s = 1 : numel(Nsteps)
    plot(PDFaggLossNPV{s}(:,1), PDFaggLossNPV{s}(:,2), ...
        'LineWidth', 2, 'Color', colSteps(s,:))
end
axis([0 2.5 0 2.5])
legend(leg)
xlabel('NPV(AL)')
ylabel('p_{NPV(AL)}')
set(gca, 'FontSize', 18)


figure; hold on
for k = 1 : numel(Nsamples)
    histogram(sensMC.NPVaggUninsuredMC(1:Nsamples(k)), ...
        'Normalization', 'pdf', ...
        'FaceColor', 'none', 'EdgeColor', colInt(k,:));
end
legend(strcat('N_{samples}=', num2str(Nsamples(:))))
axis([0 2.5 0 2.5])
xlabel('NPV(AL)')
ylabel('p_{NPV(AL)}')
set(gca, 'FontSize', 18)


figure;
for r = 1 : numel(intRates)
    plot(sensIntAggLossNPV{r,1}(:,1), sensIntAggLossNPV{r,1}(:,2), ...
        'LineWidth', 2, 'Color', colInt(r,:)); hold on
    histogram(sensIntAggLossNPV{r,2}, 'Normalization', 'pdf', ...
        'FaceColor', 'none', 'EdgeColor', colInt(r,:));
    hold off
    axis([0 2.5 0 2.5])
    legend(sprintf('Analytical - r=%1.3f', intRates(r)), 'MonteCarlo')
    xlabel('NPV(AL)')
    ylabel('p(NPV(AL))')
    set(gca, 'FontSize', 18)
    saveas(gcf, num2str(r), 'png');
    waitforbuttonpress
end
