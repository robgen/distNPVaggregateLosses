%% Input (as built with Hazard L'Aquila)

options.General.timeHorizon = 50;
options.General.intRate = 0.02;

options.Vulnerability.fragMedians = [0.165879930589601,0.322526111223683,0.340170312466871,0.460814512917018];
options.Vulnerability.fragDispersions = [0.436863904039988,0.442107162939676,0.442742956715563,0.447169859279770];
options.Vulnerability.damgeLossRatios = [0 2 10 43.5 95]/100;
options.Vulnerability.CoVdlr = [ 0 1 0.4 0.3 0.05 ];

options.Hazard.faultRate = 0.08;
options.Hazard.hazCurve = [0.166427989012818,0.0332146240000000;0.217582434028613,0.0198850450000000;0.258529430931683,0.0138629440000000;0.303930770035728,0.00988592600000000;0.354443451456181,0.00713349900000000;0.412206673094016,0.00496922700000000;0.565248464301760,0.00210721000000000;0.695119133694674,0.00102586600000000;0.846507595616605,0.000404054000000000];

%% Sensitivity

Nsteps = 1+[10, 100, 200, 500, 1000, 2000, 10000]; %, 100000];

for s = numel(Nsteps) : -1 : 1
    disp(s)
    options.Setup.NlossSamples = Nsteps(s);
    sensStep = distNPVaggregateLosses(options);
    
    sensStep = sensStep.getPDFinterarrivalTime;
    sensStep = sensStep.getPMFnumberEvents;
    sensStep = sensStep.getPDFarrivalTime;
    sensStep = sensStep.getCDFloss;
    sensStep = sensStep.getPDFlossNPV;
    sensStep = sensStep.getPDFaggregateLossNPV;
    
    PDFaggLossNPV{s} = sensStep.PDFaggLossNPV;
    for n = sensStep.NmaxEvents : -1 : 1
        areaNPVL(s,n) = trapz(...
            sensStep.PDFlossNPV(:,1), sensStep.PDFlossNPV(:,n+1));
    end
end

%% Plot

colSteps = gray(numel(Nsteps)+1);

figure; hold on
for s = 1 : numel(Nsteps)
    plot(PDFaggLossNPV{s}(:,1), PDFaggLossNPV{s}(:,2), ...
        'LineWidth', 2, 'Color', colSteps(s,:))
end
axis([0 2.5 0 5000])
xlabel('NPV(AL)')
ylabel('p(NPV(AL))')
set(gca, 'FontSize', 18)

figure; hold on
plot(Nsteps', areaNPVL, '-o')
xlabel('Nsteps in loss definition')
ylabel('Area p(NPV(AL))')
set(gca, 'FontSize', 18)
