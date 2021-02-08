%% Input (as built with Hazard L'Aquila)

options.General.timeHorizon = 50;
options.General.intRate = 0.02;

options.Vulnerability.fragMedians = [0.165879930589601,0.322526111223683,0.340170312466871,0.460814512917018];
options.Vulnerability.fragDispersions = [0.436863904039988,0.442107162939676,0.442742956715563,0.447169859279770];
options.Vulnerability.damgeLossRatios = [0 2 10 43.5 95]/100;
options.Vulnerability.CoVdlr = [ 0 1 0.4 0.3 0.05 ];

options.Hazard.faultRate = 0.08;
options.Hazard.hazCurve = [0.166427989012818,0.0332146240000000;0.217582434028613,0.0198850450000000;0.258529430931683,0.0138629440000000;0.303930770035728,0.00988592600000000;0.354443451456181,0.00713349900000000;0.412206673094016,0.00496922700000000;0.565248464301760,0.00210721000000000;0.695119133694674,0.00102586600000000;0.846507595616605,0.000404054000000000];

%% Sensitivity of the analytical formulation to the definition of the loss vector

Nsteps = 1+[10, 100, 200, 500, 1000, 2000, 10000]; %, 100000];

for s = numel(Nsteps) : -1 : 1
    disp(s)
    tic
    options.Setup.NlossSamples = Nsteps(s);
    sensStep = distNPVaggregateLosses(options);
    
    sensStep = sensStep.getPDFinterarrivalTime;
    sensStep = sensStep.getPMFnumberEvents;
    sensStep = sensStep.getPDFarrivalTime;
    sensStep = sensStep.getCDFloss;
    sensStep = sensStep.getPDFlossNPV;
    sensStep = sensStep.getPDFaggregateLossNPV;
    toc
    
    PDFaggLossNPV{s} = sensStep.PDFaggUninsuredNPV;
end

%% Sensitivity of the Montecarlo to the Nsamples

Nsamples = [100 1000 10000 50000 100000 500000];

options.Setup.NlossSamples = 100000;
sensMC = distNPVaggregateLosses(options);
sensMC = sensMC.getPDFinterarrivalTime;
sensMC = sensMC.getPMFnumberEvents;
sensMC = sensMC.getPDFarrivalTime;
sensMC = sensMC.getCDFloss;
    
for k = numel(Nsamples) : -1 : 1
    disp(k)
    tic
    
    sensMC.parameters.Setup.MCsamples = Nsamples(k);
    sensMC = sensMC.monteCarloPDFaggregateLossNPV;
    
    MCaggLossNPV{k} = sensMC.NPVaggUninsuredMC;
end

%% Sensitivity to the interest rate

intRates = [0.0001 0.005 0.01 0.02 0.04];

options.Setup.NlossSamples = 1000;
options.Setup.MCsamples = 20000;
sensInt = distNPVaggregateLosses(options);
sensInt = sensInt.getPDFinterarrivalTime;
sensInt = sensInt.getPMFnumberEvents;
sensInt = sensInt.getPDFarrivalTime;
sensInt = sensInt.getCDFloss;

for r = numel(intRates) : -1 : 1
    sensInt.parameters.General.intRate = intRates(r);
    
    sensInt = sensInt.getPDFlossNPV;
    sensInt = sensInt.getPDFaggregateLossNPV;
    sensInt = sensInt.monteCarloPDFaggregateLossNPV;
    
    sensIntAggLossNPV{r,1} = sensInt.PDFaggUninsuredNPV;
    sensIntAggLossNPV{r,2} = sensInt.NPVaggUninsuredMC;
end

%% Plot

colSteps = gray(numel(Nsteps)+1);
leg{1} = 'MonteCarlo';
for s = 1 : numel(Nsteps); leg{s+1} = num2str(Nsteps(s)); end

figure; hold on
histogram(MCaggLossNPV{end}, 'Normalization', 'pdf');
for s = 1 : numel(Nsteps)
    plot(PDFaggLossNPV{s}(:,1), PDFaggLossNPV{s}(:,2), ...
        'LineWidth', 2, 'Color', colSteps(s,:))
end
axis([0 2.5 0 2.5])
legend(leg)
xlabel('NPV(AL)')
ylabel('p(NPV(AL))')
set(gca, 'FontSize', 18)


figure; hold on
for k = 1 : numel(Nsamples)
    histogram(MCaggLossNPV{k}, 'Normalization', 'pdf');
end
axis([0 2.5 0 2.5])
xlabel('NPV(AL)')
ylabel('p(NPV(AL))')
set(gca, 'FontSize', 18)


figure; hold on
for r = 1 : numel(intRates)
    plot(sensIntAggLossNPV{r,1}(:,1), sensIntAggLossNPV{r,1}(:,2), ...
        'LineWidth', 2, 'Color', colSteps(r,:))
    histogram(sensIntAggLossNPV{r,2}, ...
        sensInt.PDFaggUninsuredNPV(:,1), 'Normalization', 'pdf');
end
axis([0 2.5 0 2.5])
xlabel('NPV(AL)')
ylabel('p(NPV(AL))')
set(gca, 'FontSize', 18)
