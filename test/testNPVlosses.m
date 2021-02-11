%% Input

options.General.timeHorizon = 50;
options.General.intRate = 0.02;
options.General.confidenceTVaR = 0.9;

% Concrete moment resisting frame (1970s)
options.Vulnerability.fragMedians = [0.166,0.32,0.34,0.46];
options.Vulnerability.fragDispersions = [0.437,0.442,0.443,0.447];
options.Vulnerability.damgeLossRatios = [0 2 10 43.5 95]/100;
options.Vulnerability.CoVdlr = [ 0 1 0.4 0.3 0.05 ];

% Seismic hazard from L'Aquila, Italy
options.Hazard.faultRate = 0.08;
options.Hazard.hazCurve = [0.166,0.0332; 0.217,0.0199; 0.258,0.014; 0.304,0.0099; 0.354,0.007; 0.412,0.005; 0.565,0.002; 0.695,0.001; 0.846,0.0004];

% insurance
options.Insurance.deductible = 0;
options.Insurance.cover = 0;
options.Insurance.coinsurance = 1;

% general setup
options.Setup.NlossSamples = 201;
options.Setup.IMstep = 0.005; % [g]
options.Setup.MCsamples = 10000;

% create object
obj = distNPVaggregateLosses(options);

%% Interarrival times

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

%% CDF of loss given one event

obj = obj.getLossDistribution;

figure('Position', [284   527   560   420]); hold on
plot(obj.LOSSdef, obj.CDFlossGivenDS, 'LineWidth', 2)

DSs = 0 : 4;
legend(strcat('DS', num2str(DSs(:))))
xlabel('Loss ratio, lr')
ylabel('P(LR\leqlr|DS)')
set(gca, 'FontSize', 18)


figure('Position', [845   527   560   420]); hold on
percentiles = [5 50 95];
for p = 1 : numel(percentiles)
    plot(obj.IMdef, prctile(obj.CDFlossGivenIM, percentiles(p), 2), ...
        'LineWidth', 2)
end
legend('5%', 'Median', '95%', 'Mean')
xlabel('IM'); ylabel('Loss ratio [-]');
set(gca, 'FontSize', 18)

figure('Position', [284    27   560   420]); hold on
plot(obj.CDFlossGivenOneEvent(:,1), obj.CDFlossGivenOneEvent(:,2), ...
    'LineWidth', 2)
plot(obj.CDFuninsuredGivenOneEvent(:,1), ...
    obj.CDFuninsuredGivenOneEvent(:,2), 'LineWidth', 2)
xlabel('Loss, l')
ylabel('P(L\leql)')
set(gca, 'FontSize', 18)


figure('Position', [845    27   560   420]); hold on
plot(obj.PDFlossGivenOneEvent(:,1), obj.PDFlossGivenOneEvent(:,2), ...
    'LineWidth', 2)
plot(obj.PDFuninsuredGivenOneEvent(:,1), ...
    obj.PDFuninsuredGivenOneEvent(:,2), 'LineWidth', 2)
axis([0 1 0 10])
xlabel('Loss, L')
ylabel('p(L)')
set(gca, 'FontSize', 18)

%% PDF loss NPV

obj = obj.getPDFlossNPV;

figure('Position', [279   527   560   420]); hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFunitCashFlowNPV(:,1), obj.PDFunitCashFlowNPV(:,n+1), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end
legend(strcat('event ', num2str(nev(:))), 'Location', 'North')
xlabel('NPV unit cash flow, NPV^{(1)}')
ylabel('p_{NPV^{(1)}}')
set(gca, 'FontSize', 18)

figure('Position', [840   527   560   420]); hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFuninsuredNPV(:,1), obj.PDFuninsuredNPV(:,n+1), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end
axis([0 0.3 0 10])
legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
xlabel('NPV loss, NPV(L)')
ylabel('p_{NPV(L)}')
set(gca, 'FontSize', 18)

figure('Position', [560    27   560   420]); hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFuninsuredNPV(:,1), ...
        cumtrapz(obj.PDFuninsuredNPV(:,1),obj.PDFuninsuredNPV(:,n+1)), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end
legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
xlabel('NPV loss, npv(L)')
ylabel('p(NPV(L)>npv(L))')
set(gca, 'FontSize', 18)

%% PDF aggregate loss NPV

obj = obj.getAggregateLossNPVdist;

figure; hold on
for n = 1 : obj.NmaxEvents
    plot(obj.PDFaggUninsuredNPVGivenNevents(:,1), ...
        obj.PDFaggUninsuredNPVGivenNevents(:,n+1), ...
        'LineWidth', 2, 'Color', colEvents(n,:))
end
axis([0 2.5 0 obj.PDFaggUninsuredNPV(20,2)])
legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
xlabel('NPV(AL)')
ylabel('p(NPV(AL))')
set(gca, 'FontSize', 18)

figure; hold on
plot(obj.PDFaggUninsuredNPV(:,1), obj.PDFaggUninsuredNPV(:,2), ...
     'LineWidth', 2)
axis([0 2.5 0 obj.PDFaggUninsuredNPV(20,2)])
xlabel('NPV(AL)')
ylabel('p_{NPV(AL)}')
set(gca, 'FontSize', 18)

figure; hold on
plot(obj.CDFaggUninsuredNPV(:,1), obj.CDFaggUninsuredNPV(:,2), ...
     'LineWidth', 2)
xlabel('NPV(AL)')
ylabel('p(NPV(AL))')
set(gca, 'FontSize', 18)

%% Compare with Montecarlo

obj = obj.monteCarloPDFaggregateLossNPV;

figure; hold on
h = histogram(obj.NPVaggUninsuredMC, obj.PDFaggUninsuredNPV(1:256:end,1), ...
    'Normalization', 'pdf', 'FaceColor', 0.94*[1 1 1]);
plot(obj.PDFaggUninsuredNPV(:,1), obj.PDFaggUninsuredNPV(:,2), ...
     'LineWidth', 2, 'Color', [0 0.447058823529412 0.741176470588235])
legend('MonteCarlo', 'Analytical')
axis([0 2.5 0 2.5])
xlabel('NPV(AL) [-]')
ylabel('PDF_{NPV(AL)}')
set(gca, 'FontSize', 18)
saveas(gcf, 'PDFaggLoss', 'png'); close

%% Test PDF areas

toll = 0.03;

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
if options.General.intRate ~= 0
    for n = obj.NmaxEvents : -1 : 1
        areaNPV1(n) = trapz(obj.PDFunitCashFlowNPV(:,1), ...
            obj.PDFunitCashFlowNPV(:,n+1));
    end
    assert(any(abs(areaNPV1-1)<toll), 'area under the NPV1 PDFs is not one')
end

% NPVL
for n = obj.NmaxEvents : -1 : 1
    areaNPVL(n) = trapz(obj.PDFuninsuredNPV(:,1), obj.PDFuninsuredNPV(:,n+1));
end
assert(any(abs(areaNPVL-1)<toll), 'area under the NPVL PDFs is not one')

% NPV(AL)|Nevents
for n = obj.NmaxEvents : -1 : 1
    areaNPVaggN(n) = trapz(obj.PDFaggUninsuredNPVGivenNevents(:,1), ...
        obj.PDFaggUninsuredNPVGivenNevents(:,n+1));
end
assert(any(abs(areaNPVaggN-1)<toll), 'area under the NPV(AL)|Nevents PDFs is not one')

% NPV(AL)
npvalArea = trapz(obj.PDFaggUninsuredNPV(:,1), obj.PDFaggUninsuredNPV(:,2));
assert(abs(npvalArea-1)<toll, 'area under the NPV(AL) PDF is not one')
