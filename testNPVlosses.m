%% Input (as built with Hazard L'Aquila)

options.General.timeHorizon = 50;
options.General.intRate = 0.02;

options.Vulnerability.fragMedians = [0.165879930589601,0.322526111223683,0.340170312466871,0.460814512917018];
options.Vulnerability.fragDispersions = [0.436863904039988,0.442107162939676,0.442742956715563,0.447169859279770];
options.Vulnerability.damgeLossRatios = [0 2 10 43.5 95]/100;
options.Vulnerability.CoVdlr = [ 0 1 0.4 0.3 0.05 ];

options.Hazard.faultRate = 0.08;
options.Hazard.hazCurve = [0.166427989012818,0.0332146240000000;0.217582434028613,0.0198850450000000;0.258529430931683,0.0138629440000000;0.303930770035728,0.00988592600000000;0.354443451456181,0.00713349900000000;0.412206673094016,0.00496922700000000;0.565248464301760,0.00210721000000000;0.695119133694674,0.00102586600000000;0.846507595616605,0.000404054000000000];

options.Setup.NlossSamples = 501;
obj = distNPVaggregateLosses(options);

%% interarrival times

obj = obj.getPDFinterarrivalTime;

figure; hold on
plot(obj.PDFinterarrivalTime(:,1), obj.PDFinterarrivalTime(:,2), ...
    'LineWidth', 2)

histogram(random('Exponential', 1/options.Hazard.faultRate, 10000, 1), ...
    'BinWidth', 1, 'Normalization', 'pdf');

xlabel('Interarrival time \Delta\tau')
ylabel('p(\Delta\tau)')
set(gca, 'FontSize', 18)

%% P(Nevents)

obj = obj.getPMFnumberEvents;

figure;
bar(obj.PMFnumberEvents(:,1), obj.PMFnumberEvents(:,2))
xlabel('Number of events N_{ev}')
ylabel('P(N_{ev})')
set(gca, 'FontSize', 18)

%% Arrival times

obj = obj.getPDFarrivalTime;

colEvents = gray(obj.NmaxEvents+1);

figure; hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFarrivalTime(:,1), obj.PDFarrivalTime(:,n+1), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end

plot(obj.PDFinterarrivalTime(:,1), obj.PDFinterarrivalTime(:,2), ...
    '--r', 'LineWidth', 4)

nev = 1 : obj.NmaxEvents;
legend(strcat('N_{ev}=', num2str(nev(:))))

xlabel('Arrival time \tau')
ylabel('p(\tau)')
set(gca, 'FontSize', 18)

%% CDF of one loss given one event

obj = obj.getCDFloss;

figure; hold on
plot(obj.LOSSdef, obj.CDFlossGivenDS, 'LineWidth', 2)

DSs = 0 : 4;
legend(strcat('DS', num2str(DSs(:))))
xlabel('Loss ratio, lr')
ylabel('P(LR\leqlr|DS)')
set(gca, 'FontSize', 18)


figure; hold on
percentiles = [5 50 95];
for p = 1 : numel(percentiles)
    plot(obj.IMdef, prctile(obj.CDFlossGivenIM, percentiles(p), 2), ...
        'LineWidth', 2)
end
legend('5%', 'Median', '95%', 'Mean')
xlabel('IM'); ylabel('Loss ratio [-]');
set(gca, 'FontSize', 18)

figure; hold on
plot(obj.CDFlossGivenOneEvent(:,1), obj.CDFlossGivenOneEvent(:,2), ...
    'LineWidth', 2)
xlabel('Loss, l')
ylabel('P(L\leql)')
set(gca, 'FontSize', 18)


figure; hold on
plot(obj.PDFlossGivenOneEvent(:,1), obj.PDFlossGivenOneEvent(:,2), ...
    'LineWidth', 2)
axis([0 1 0 10])
xlabel('Loss, L')
ylabel('p(L)')
set(gca, 'FontSize', 18)

%% PDF loss NPV

obj = obj.getPDFlossNPV;

figure; hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFunitCashFlowNPV(:,1), obj.PDFunitCashFlowNPV(:,n+1), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end
legend(strcat('event ', num2str(nev(:))), 'Location', 'North')
xlabel('NPV unit cash flow, NPV^{(1)}')
ylabel('p(NPV^{(1)})')
set(gca, 'FontSize', 18)


figure; hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFlossNPV(:,1), obj.PDFlossNPV(:,n+1), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end
axis([0 0.3 0 10])
legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
xlabel('NPV loss, NPV(L)')
ylabel('p(NPV(L))')
set(gca, 'FontSize', 18)

figure; hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFlossNPV(:,1), ...
        cumtrapz(obj.PDFlossNPV(:,1),obj.PDFlossNPV(:,n+1)), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end
legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
xlabel('NPV loss, npv(L)')
ylabel('p(NPV(L)>npv(L))')
set(gca, 'FontSize', 18)

%% PDF aggregate loss NPV

obj = obj.getPDFaggregateLossNPV;

figure; hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFaggLossNPVGivenNevents(:,1), ...
        obj.PDFaggLossNPVGivenNevents(:,n+1), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end
axis([0 2.5 0 obj.PDFaggLossNPV(20,2)])
legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
xlabel('NPV(AL)')
ylabel('p(NPV(AL))')
set(gca, 'FontSize', 18)

figure; hold on
plot(obj.PDFaggLossNPV(:,1), obj.PDFaggLossNPV(:,2), ...
     'LineWidth', 2)
 
axis([0 2.5 0 0.5])
xlabel('NPV(AL)')
ylabel('p(NPV(AL))')
set(gca, 'FontSize', 18)

%% Test PDF areas

toll = 0.01;

% interarrival
totA = trapz(obj.PDFinterarrivalTime(:,1), obj.PDFinterarrivalTime(:,2));
assert(abs(totA-1)<toll, 'area under the \Delta \tau PDF is not one')

% N events
assert(abs(sum(obj.PMFnumberEvents(:,2))-1)<toll, ...
    'Area of N_{ev} PMF is not one')

% arrival
for n = obj.NmaxEvents : -1 : 1
    totArea(n) = trapz(obj.PDFarrivalTime(:,1), obj.PDFarrivalTime(:,n+1));
end
assert(any(abs(totArea-1)<toll), 'area under the \tau PDFs is not one')

% loss|DS
clear derivative
for ds = 5 : -1 : 2
    derivative = obj.numericalDerivative(...
        [obj.LOSSdef, obj.CDFlossGivenDS(:,ds)]);
    
    dsArea(ds-1) = trapz(derivative(:,1),derivative(:,2));
end

assert(any(abs(dsArea-1)<toll), 'area under the Loss|DS PDFs is not one')

% IM
pdfIM(:,1) = obj.MAFim(:,1);
pdfIM(:,2) = obj.MAFim(:,2)/obj.parameters.Hazard.faultRate;
imArea = trapz(pdfIM(:,1),pdfIM(:,2));
%assert(abs(imArea-1)<0.01, 'area under the PDF(IM) is not one')

% loss|one event
lossArea = trapz(obj.PDFlossGivenOneEvent(:,1), obj.PDFlossGivenOneEvent(:,2));
assert(abs(lossArea-1)<toll, 'area under the loss|1ev PDF is not one')

% NPV1
for n = obj.NmaxEvents : -1 : 1
    areaNPV1(n) = trapz(obj.PDFunitCashFlowNPV(:,1), ...
        obj.PDFunitCashFlowNPV(:,n+1));
end
assert(any(abs(areaNPV1-1)<toll), 'area under the NPV1 PDFs is not one')

% NPVL
for n = obj.NmaxEvents : -1 : 1
    areaNPVL(n) = trapz(obj.PDFlossNPV(:,1), obj.PDFlossNPV(:,n+1));
end
assert(any(abs(areaNPVL-1)<toll), 'area under the NPVL PDFs is not one')

% NPV(AL)|Nevents
for n = obj.NmaxEvents : -1 : 1
    areaNPVaggN(n) = trapz(obj.PDFaggLossNPVGivenNevents(:,1), ...
        obj.PDFaggLossNPVGivenNevents(:,n+1));
end
assert(any(abs(areaNPVaggN-1)<toll), 'area under the NPV(AL)|Nevents PDFs is not one')

% NPV(AL)
npvalArea = trapz(obj.PDFaggLossNPV(:,1), obj.PDFaggLossNPV(:,2));
assert(abs(npvalArea-1)<toll, 'area under the NPV(AL) PDF is not one')
