classdef distNPVaggregateLosses
    %distNPVaggregateLosses
    %
    
    properties
        
        parameters
        
        TAUdef
        IMdef
        LOSSdef
        fragilities
        
        PDFinterarrivalTime
        PDFarrivalTime
        PMFnumberEvents
        NmaxEvents
        
        PMFds
        CDFlossGivenDS
        CDFlossGivenIM
        CDFlossGivenOneEvent
        PDFlossGivenOneEvent
        
        PDFuninsuredGivenOneEvent
        CDFuninsuredGivenOneEvent
        
        MAFim
        hazCurveResampled
        lastStepIM % this should be private
        
        PDFunitCashFlowNPV
        PDFuninsuredNPV
        PDFaggUninsuredNPVGivenNevents
        PDFaggUninsuredNPV
        CDFaggUninsuredNPV
        
        EAL
        EALu
        
        TVaRu
        
        NPVaggUninsuredMC
        NPVaggInsuredMC
    end
    
    methods
        
        %%% Macro functions
        function self = distNPVaggregateLosses(options)
            %distNPVaggregateLosses is the constructor
            
            if nargin == 0, options = struct; end
            
            self = setAllParameters(self, options);
            self = setVariableRanges(self);
        end
        
        
        function self = getPDFinterarrivalTime(self)
            %getPDFinterarrivalTime exponential PDF for deltaTau
            
            v = self.parameters.Hazard.faultRate;
            
            self.PDFinterarrivalTime(:,1) = self.TAUdef;
            self.PDFinterarrivalTime(:,2) = v * ...
                exp(-self.PDFinterarrivalTime(:,1)*v);
            
        end
        
        
        function self = getPDFarrivalTime(self, NmaxEvents)
            %getPDFarrivalTime Erlang PDF for Tau|n events
            
            if isempty(self.PDFinterarrivalTime)
                self = self.getPDFinterarrivalTime;
            end
            
            if nargin < 2
                if isempty(self.NmaxEvents)
                    self = self.getPMFnumberEvents;
                end
                NmaxEvents = self.NmaxEvents;
            end
            
            v = self.parameters.Hazard.faultRate;
            
            self.PDFarrivalTime(:,1) = self.TAUdef;
            
            for n = NmaxEvents : -1 : 1
                self.PDFarrivalTime(:,n+1) = v^n * ...
                    (self.PDFinterarrivalTime(:,1).^(n-1) / ...
                    factorial(n-1) ) .* ...
                    exp(-v*self.PDFinterarrivalTime(:,1));
            end
        end
        
        
        function self = getPMFnumberEvents(self)
            %getPMFnumberEvents Poisson PMF for Nevents in timeHorizon
            
            faultRate = self.parameters.Hazard.faultRate;
            timeHorizon = self.parameters.General.timeHorizon;
            
            % performance improvement needed
            for n = 0 : 10
                pNevents = ((faultRate*timeHorizon)^n * ...
                    exp(-faultRate*timeHorizon)) ./ factorial(n);
                self.PMFnumberEvents(n+1,2) = pNevents;
            end
            
            toll = 0.005;
            while pNevents > toll
                n = n + 1;
                pNevents = ((faultRate*timeHorizon)^n * ...
                    exp(-faultRate*timeHorizon)) ./ factorial(n);
                self.PMFnumberEvents(n+1,2) = pNevents;
            end
            
            self.PMFnumberEvents(:,1) = 0 : n;
            self.NmaxEvents = n;
        end
        
        
        function self = getLossDistribution(self)
            
            self = getFragilities(self);
            self = getPMFds(self);
            self = getCDFlossGivenDS(self);
            self = getCDFlossGivenIM(self);
            self = getMAFim(self);
            
            mafIMmatrix = repmat(self.MAFim(:,2), [1 numel(self.LOSSdef)]);
            self.CDFlossGivenOneEvent(:,1) = self.LOSSdef;
            self.CDFlossGivenOneEvent(:,2) = ...
                sum(self.CDFlossGivenIM(1:self.lastStepIM,:) .* mafIMmatrix) ./ ...
                sum(self.MAFim(:,2));
            
            self.PDFlossGivenOneEvent = self.numericalGradient(...
                self.CDFlossGivenOneEvent);
            
            self = self.groundUpToUninsuredPDF;
            self.EAL = self.tailValueAtRisk(self.CDFlossGivenOneEvent, 0);
            self.EALu = self.tailValueAtRisk(self.CDFuninsuredGivenOneEvent, 0);
        end
                
        
        function self = getPDFlossNPV(self)
            
            if isempty(self.PDFuninsuredGivenOneEvent)
                self = self.getLossDistribution;
            end
            
            self = getPDFunitCashFlowNPV(self);
            
            % the stepsize for the interpolation must be smaller than the
            % definition of the PDF(loss)
            Nsteps = (self.parameters.Setup.NlossSamples-1)*200 + 1;
            z = linspace(eps, 1, Nsteps)';
            
            if self.parameters.General.intRate == 0
                PDFuninsuredGivenOneEventInterp = interp1(...
                    self.PDFuninsuredGivenOneEvent(:,1), ...
                    self.PDFuninsuredGivenOneEvent(:,2), z);
                self.PDFuninsuredNPV(:,1) = z;
                self.PDFuninsuredNPV(:,2:self.NmaxEvents+1) = repmat(...
                    PDFuninsuredGivenOneEventInterp, 1, self.NmaxEvents);
            else
                pdf2interp = [];
                for n = self.NmaxEvents : -1 : 1
                    [testPDF, pdf2interp] = self.productDistribution(...
                        self.PDFunitCashFlowNPV(:,[1,n+1]), ...
                        self.PDFuninsuredGivenOneEvent, z, pdf2interp);
                    
                    self.PDFuninsuredNPV(:,n+1) = testPDF(:,2);
                end
                
                self.PDFuninsuredNPV(:,1) = z;
            end
        end
        
        
        function self = getAggregateLossNPVdist(self)
            
            if isempty(self.PDFuninsuredNPV)
                self = self.getPDFlossNPV;
            end
            
            self.PDFaggUninsuredNPVGivenNevents = self.recursiveConvolution(...
                self.PDFuninsuredNPV);
            
            pNevMatrix = self.PMFnumberEvents(:,2)' + ...
                zeros(size(self.PDFaggUninsuredNPVGivenNevents,1),1);
            
            self.PDFaggUninsuredNPV(:,1) = self.PDFaggUninsuredNPVGivenNevents(:,1);
            self.PDFaggUninsuredNPV(:,2) = sum(pNevMatrix(:,2:end) .* ...
                self.PDFaggUninsuredNPVGivenNevents(:,2:end), 2);
            % note: P(L>l|0events) = 0, for each l
            
            self.CDFaggUninsuredNPV = ...
                self.numericalIntegral(self.PDFaggUninsuredNPV);
            
            alpha = self.parameters.General.confidenceTVaR;
            self.TVaRu = self.tailValueAtRisk(self.CDFaggUninsuredNPV, alpha);
        end
        
        
        function self = monteCarloPDFaggregateLossNPV(self)
            
            if isempty(self.CDFlossGivenOneEvent)
                self = self.getLossDistribution;
            end
            
            for n = self.parameters.Setup.MCsamples : -1 : 1
                timeSim = self.poissonProcess(...
                    self.parameters.Hazard.faultRate, ...
                    self.parameters.General.timeHorizon, 'n');
                Nevents(n,1) = numel(timeSim);
                
                top = self.CDFlossGivenOneEvent(end,2)-1000*eps;
                filteredLossCDF = [ ...
                    self.CDFlossGivenOneEvent(self.CDFlossGivenOneEvent(:,2)<top,:); ...
                    self.CDFlossGivenOneEvent(end,:) ];
                lossSim = interp1(filteredLossCDF(:,2), ...
                    filteredLossCDF(:,1), rand(Nevents(n),1) );
                
                [uninsuredLossSim, insuredLossSim] = self.applyPayout(...
                    lossSim, self.parameters.Insurance.deductible, ...
                    self.parameters.Insurance.cover, ...
                    self.parameters.Insurance.coinsurance);
                
                self.NPVaggUninsuredMC(n,1) = sum(uninsuredLossSim ./ ...
                    (1+self.parameters.General.intRate).^timeSim );
                
                self.NPVaggInsuredMC(n,1) = sum(insuredLossSim ./ ...
                    (1+self.parameters.General.intRate).^timeSim);
            end
            
        end
        
        
        function plotAllResults(self)
            
            colEvents = gray(self.NmaxEvents+1);
            nev = 1 : self.NmaxEvents;
            
            figure;
            bar(self.PMFnumberEvents(:,1), self.PMFnumberEvents(:,2))
            xlabel('Number of events N_{ev}')
            ylabel('P(N_{ev})')
            set(gca, 'FontSize', 18)
            
            figure; hold on
            for n = 1 : self.NmaxEvents
                plot(self.PDFarrivalTime(:,1), self.PDFarrivalTime(:,n+1), ...
                    'LineWidth', 2, 'Color', colEvents(n,:))
            end
            legend(strcat('N_{ev}=', num2str(nev(:))))
            xlabel('Arrival time \tau')
            ylabel('p(\tau)')
            set(gca, 'FontSize', 18)
            
            figure('Position', [284   527   560   420]); hold on
            plot(self.LOSSdef, self.CDFlossGivenDS, 'LineWidth', 2)
            
            DSs = 0 : 4;
            legend(strcat('DS', num2str(DSs(:))))
            xlabel('Loss ratio, lr')
            ylabel('P(LR\leqlr|DS)')
            set(gca, 'FontSize', 18)
            
            
            figure('Position', [845   527   560   420]); hold on
            percentiles = [5 50 95];
            for p = 1 : numel(percentiles)
                plot(self.IMdef, prctile(self.CDFlossGivenIM, percentiles(p), 2), ...
                    'LineWidth', 2)
            end
            legend('5%', 'Median', '95%', 'Mean')
            xlabel('IM'); ylabel('Loss ratio [-]');
            set(gca, 'FontSize', 18)
            
            figure('Position', [284    27   560   420]); hold on
            plot(self.CDFlossGivenOneEvent(:,1), self.CDFlossGivenOneEvent(:,2), ...
                'LineWidth', 2)
            plot(self.CDFuninsuredGivenOneEvent(:,1), ...
                self.CDFuninsuredGivenOneEvent(:,2), 'LineWidth', 2)
            xlabel('Loss, l')
            ylabel('P(L\leql)')
            set(gca, 'FontSize', 18)
            
            
            figure('Position', [845    27   560   420]); hold on
            plot(self.PDFlossGivenOneEvent(:,1), self.PDFlossGivenOneEvent(:,2), ...
                'LineWidth', 2)
            plot(self.PDFuninsuredGivenOneEvent(:,1), ...
                self.PDFuninsuredGivenOneEvent(:,2), 'LineWidth', 2)
            axis([0 1 0 10])
            xlabel('Loss, L')
            ylabel('p(L)')
            set(gca, 'FontSize', 18)
            
            figure('Position', [279   527   560   420]); hold on
            for n = 1 : self.NmaxEvents
                plot(self.PDFunitCashFlowNPV(:,1), self.PDFunitCashFlowNPV(:,n+1), ...
                    'LineWidth', 2, 'Color', colEvents(n,:))
            end
            legend(strcat('event ', num2str(nev(:))), 'Location', 'North')
            xlabel('NPV unit cash flow, NPV^{(1)}')
            ylabel('p_{NPV^{(1)}}')
            set(gca, 'FontSize', 18)
            
            figure('Position', [840   527   560   420]); hold on
            for n = 1 : self.NmaxEvents
                plot(self.PDFuninsuredNPV(:,1), self.PDFuninsuredNPV(:,n+1), ...
                    'LineWidth', 2, 'Color', colEvents(n,:))
            end
            axis([0 0.3 0 10])
            legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
            xlabel('NPV loss, NPV(L)')
            ylabel('p_{NPV(L)}')
            set(gca, 'FontSize', 18)
            
            figure('Position', [560    27   560   420]); hold on
            for n = 1 : self.NmaxEvents
                plot(self.PDFuninsuredNPV(:,1), ...
                    cumtrapz(self.PDFuninsuredNPV(:,1),self.PDFuninsuredNPV(:,n+1)), ...
                    'LineWidth', 2, 'Color', colEvents(n,:))
            end
            legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
            xlabel('NPV loss, npv(L)')
            ylabel('p(NPV(L)>npv(L))')
            set(gca, 'FontSize', 18)
            
            figure; hold on
            for n = 1 : self.NmaxEvents
                plot(self.PDFaggUninsuredNPVGivenNevents(:,1), ...
                    self.PDFaggUninsuredNPVGivenNevents(:,n+1), ...
                    'LineWidth', 2, 'Color', colEvents(n,:))
            end
            axis([0 2.5 0 self.PDFaggUninsuredNPV(20,2)])
            legend(strcat('event ', num2str(nev(:))), 'Location', 'NorthEast')
            xlabel('NPV(AL)')
            ylabel('p(NPV(AL))')
            set(gca, 'FontSize', 18)
            
            figure; hold on
            plot(self.PDFaggUninsuredNPV(:,1), self.PDFaggUninsuredNPV(:,2), ...
                'LineWidth', 2)
            axis([0 2.5 0 self.PDFaggUninsuredNPV(20,2)])
            xlabel('NPV(AL)')
            ylabel('p_{NPV(AL)}')
            set(gca, 'FontSize', 18)
            
            figure; hold on
            plot(self.CDFaggUninsuredNPV(:,1), self.CDFaggUninsuredNPV(:,2), ...
                'LineWidth', 2)
            xlabel('NPV(AL)')
            ylabel('p(NPV(AL))')
            set(gca, 'FontSize', 18)
            
            figure; hold on
            h = histogram(self.NPVaggUninsuredMC, self.PDFaggUninsuredNPV(1:256:end,1), ...
                'Normalization', 'pdf', 'FaceColor', 0.94*[1 1 1]);
            plot(self.PDFaggUninsuredNPV(:,1), self.PDFaggUninsuredNPV(:,2), ...
                'LineWidth', 2, 'Color', [0 0.447058823529412 0.741176470588235])
            legend('MonteCarlo', 'Analytical')
            axis([0 2.5 0 2.5])
            xlabel('NPV(AL) [-]')
            ylabel('PDF_{NPV(AL)}')
            set(gca, 'FontSize', 18)            
        end
        
        
        %%% Micro functions
        function self = getFragilities(self)
            %getFragilities computes lognormal fragility curves
            
            medians = self.parameters.Vulnerability.fragMedians;
            dispersions = self.parameters.Vulnerability.fragDispersions;
            
            for ds = numel(medians) : -1 : 1
                self.fragilities(:,ds) = logncdf(self.IMdef, ...
                    log(medians(ds)), dispersions(ds));
            end
        end
        
        
        function self = getPMFds(self)
            %getPMFds computes the probability of being in DSk|IM
            
            frag = self.fragilities;
            
            self.PMFds = [ones(size(frag,1),1) frag] - ...
                [frag zeros(size(frag,1),1)];
        end
        
        
        function self = getCDFlossGivenDS(self)
            %getPMFds computes the beta distributions P(L<l|DSk)
            
            dlr = self.parameters.Vulnerability.damgeLossRatios;
            CoVdlr = self.parameters.Vulnerability.CoVdlr;
            
            alfa = (1-dlr)./CoVdlr.^2 - dlr;
            beta = alfa .* (1-dlr) ./ dlr;
            
            for ds = numel(dlr) : -1 : 2
                self.CDFlossGivenDS(:,ds) = betacdf(...
                    self.LOSSdef, alfa(ds), beta(ds));
            end
            
            %self.CDFlossGivenDS(:,1) = 1; % P(L<l|DS0) = 1 for any l
            
            % P(L<l|DS0) is a step function centered at zero
            self.CDFlossGivenDS(2:end,1) = 1;
            % this would be correct... but gives issues with PDFaggUninsuredNPV
        end
        
        
        function self = getCDFlossGivenIM(self)
            
            %%% needs performance improvement
            for im = numel(self.IMdef) : -1 : 1
                CDFlossDSmatrix(im,:,:) = self.CDFlossGivenDS;
            end
            
            for l = numel(self.LOSSdef) : -1 : 1
                probDS3d(:,l,:) = self.PMFds;
            end
            %%% needs performance improvement
            
            self.CDFlossGivenIM = sum( CDFlossDSmatrix .* probDS3d, 3 );
        end
        
        
        function self = getMAFim(self)
            
            if max(self.IMdef) < max(self.parameters.Hazard.hazCurve(:,1))
                error('Extend the IM definition of CDFlossIM')
            end
            
            if self.parameters.Hazard.hazCurve(1,1) ~= 0
                extendedHazCurve = [0 self.parameters.Hazard.faultRate; ...
                    self.parameters.Hazard.hazCurve];
            else
                extendedHazCurve = self.parameters.Hazard.hazCurve;
            end
            
            self.MAFim = self.numericalDerivative(extendedHazCurve);
            self.MAFim(:,2) = abs(self.MAFim(:,2));
            
            [self.MAFim, self.lastStepIM] = self.resampleCurve(...
                self.MAFim, self.IMdef);
        end
        
        
        function self = getPDFunitCashFlowNPV(self)
            
            if isempty(self.NmaxEvents)
                self = self.getPMFnumberEvents;
                self = self.getPDFarrivalTime;
            end
            
            r = self.parameters.General.intRate;
            v = self.parameters.Hazard.faultRate;
            
            npv = self.LOSSdef;
            
            log1R = log(1./npv) / log(1+r); % base 1+r
            
            self.PDFunitCashFlowNPV(:,1) = npv;
            for n = self.NmaxEvents : -1 : 1
                self.PDFunitCashFlowNPV(:,n+1) = ...
                    1./(log(1+r)*npv) .* ...
                    ( v^n * log1R.^(n-1) / factorial(n-1) ) .* ...
                    exp(-v*log1R);
            end
        end
        
        
        function self = groundUpToUninsuredPDF(self)
            
            d = self.parameters.Insurance.deductible;
            c = self.parameters.Insurance.cover;
            %coI = self.parameters.Insurance.coinsurance; % should be included
            
            if c ~= 0 && d >= c
                error('Deductible is greater than cover')
            end
            
            if c ~= 0
                self.PDFuninsuredGivenOneEvent(:,1) = self.PDFlossGivenOneEvent(:,1);
                
                zeroPayout = find(self.PDFuninsuredGivenOneEvent(:,1)<=d+eps);
                lastZeroPayout = zeroPayout(end);
                linearPayout = find(self.PDFuninsuredGivenOneEvent(:,1)>d+eps & ...
                    self.PDFuninsuredGivenOneEvent(:,1)<c);

                maxPayout = find(self.PDFuninsuredGivenOneEvent(:,1)>=c);
                
                self.PDFuninsuredGivenOneEvent(zeroPayout,2) = ...
                    self.PDFlossGivenOneEvent(zeroPayout,2);
                self.PDFuninsuredGivenOneEvent(...
                    lastZeroPayout+1:lastZeroPayout+numel(maxPayout),2) = ...
                    self.PDFlossGivenOneEvent(maxPayout,2);
                
                % add delta function at uninsuredLoss = D (Triangular approx)
                probLossBetweenDandC = trapz(...
                    self.PDFlossGivenOneEvent(linearPayout,1), ...
                    self.PDFlossGivenOneEvent(linearPayout,2));
                
                widthDelta = 4; % steps (MUST BE EVEN)
                stepSize = (self.PDFlossGivenOneEvent(2,1) - ...
                    self.PDFlossGivenOneEvent(1,1));
                heightDelta = 2 * probLossBetweenDandC / ...
                    (widthDelta * stepSize);
                
                deltaFx(:,1) = linspace(0, 2, widthDelta+1);
                deltaFx(:,2) = interp1([0 1 2], ...
                    [0 heightDelta 0], deltaFx(:,1));
                
                stepsForDelta = zeroPayout(end)-widthDelta/2 : ...
                    zeroPayout(end)+widthDelta/2;
                self.PDFuninsuredGivenOneEvent(stepsForDelta,2) = ...
                    self.PDFuninsuredGivenOneEvent(stepsForDelta,2) + deltaFx(:,2);
                
                % calculate CDF
                self.CDFuninsuredGivenOneEvent(:,1) = ...
                    self.PDFuninsuredGivenOneEvent(:,1);
                self.CDFuninsuredGivenOneEvent(:,2) = cumtrapz(...
                    self.PDFuninsuredGivenOneEvent(:,1), ...
                    self.PDFuninsuredGivenOneEvent(:,2));
            else
                self.PDFuninsuredGivenOneEvent = self.PDFlossGivenOneEvent;
                self.CDFuninsuredGivenOneEvent = self.CDFlossGivenOneEvent;
            end
            
        end
        
        
        %%% Internal stuff
        function self = setVariableRanges(self)
            
            self.TAUdef = linspace(0, ...
                8*self.parameters.General.timeHorizon, 100);
            
            self.IMdef = [0; (eps:0.005: ...
                2*max(self.parameters.Vulnerability.fragMedians))'];
            
            Nsteps = self.parameters.Setup.NlossSamples;
            self.LOSSdef = linspace(eps, 1, Nsteps)';
            % MUST BE equally spaced for convolution. MUST start from eps
        end
        
        
        function self = setAllParameters(self, options)
            % setAllParameters deals with the optional parameters
            
            % build basic parameters
            macroFieldsPar = {'General', 'Vulnerability', 'Hazard', 'Insurance', 'Setup'};
            
            microFieldsPar{1} = {'timeHorizon', 'intRate', 'confidenceTVaR' };
            microFieldsParVals{1} = {50, 0.02, 0.9};
            
            microFieldsPar{2} = {'fragMedians', 'fragDispersions', 'damgeLossRatios', 'CoVdlr' };
            microFieldsParVals{2} = {[0.05 0.35 0.65 1.25], [0.2 0.2 0.2 0.2], [0 2 10 43.5 95]/100, [ 0 1 0.4 0.3 0.05 ]};
            
            microFieldsPar{3} = {'hazCurve', 'faultRate' };
            microFieldsParVals{3} = {[0.166,0.0332; 0.217,0.0199; 0.258,0.014; 0.304,0.0099; 0.354,0.007; 0.412,0.005; 0.565,0.002; 0.695,0.001; 0.846,0.0004] , 0.08};
            
            microFieldsPar{4} = {'deductible', 'cover', 'coinsurance' };
            microFieldsParVals{4} = {0, 0, 1};
            
            microFieldsPar{5} = {'NlossSamples', 'IMstep', 'MCsamples' };
            microFieldsParVals{5} = {201, 0.005, 10000};
            
            for F = 1 : numel(macroFieldsPar)
                for f = 1 : numel(microFieldsPar{F})
                    self.parameters.(macroFieldsPar{F}).(microFieldsPar{F}{f}) = microFieldsParVals{F}{f};
                end
            end
            
            % overwrite fields if some parameter is specified
            macroFieldsOptional = fieldnames(options);
            for OF = 1 : numel(macroFieldsOptional)
                microFieldsOptional = fieldnames(options.(macroFieldsOptional{OF}));
                for of = 1 : numel(microFieldsOptional)
                    if isfield(self.parameters.(macroFieldsOptional{OF}), microFieldsOptional{of}) == 1
                        self.parameters.(macroFieldsOptional{OF}).(microFieldsOptional{of}) = options.(macroFieldsOptional{OF}).(microFieldsOptional{of});
                    else
                        error('Field %s.%s is not allowed', macroFieldsOptional{OF}, microFieldsOptional{of})
                    end
                end
            end
        end
        
        
    end
    
    methods(Static)
        
        function [resampledCurve, lastStep] = resampleCurve(...
                toResample, toGetXvalues)
            
            temp = find(toGetXvalues >= max(toResample(:,1)));
            lastStep = temp(1)-1;
            
            resampledCurve(:,2) = interp1(...
                toResample(:,1), toResample(:,2), ...
                toGetXvalues(1:lastStep));
            
            resampledCurve(:,1) = toGetXvalues(1:lastStep);
        end
        
        
        function derivative = numericalGradient(curve)
            
            % Note: curve must have EQUALLY-SPACED points
            if numel( uniquetol(diff(curve(:,1)), 1000*eps) ) > 1
                warning('Make sure that "curve" contains equally spaced points')
            end
            
            derivative(:,1) = curve(:,1);
            derivative(:,2) = gradient(curve(:,2), curve(2,1)-curve(1,1));
        end
        
        
        function derivative = numericalDerivative(curve)
            
            dx = diff(curve(:,1));
            dy = diff(curve(:,2));
            
            % assume derivative is located at the beginning of the step
            derivative(:,1) = curve(:,1);
            derivative(1:end-1,2) = dy./dx;
            
            derivative(end,2) = derivative(end-1);
        end
        
        
        function integral = numericalIntegral(curve)
            
            integral(:,1) = curve(:,1);
            integral(:,2) = cumtrapz(curve(:,1), curve(:,2));
        end
        
        
        function convRecurs = recursiveConvolution(FXsToConvolute)
            
            % Note: the PDFs must have EQUALLY-SPACED points
            if numel( uniquetol(diff(FXsToConvolute(:,1)), 1000*eps) ) > 1
                warning('Make sure LOSSdef contains equally spaced points')
            end
            
            stepSize = FXsToConvolute(2,1) - FXsToConvolute(1,1);
            len = size(FXsToConvolute,1);
            n = size(FXsToConvolute,2)-1;
            
            convRecurs = zeros(n*len-n+1, n+1);
            for ii = 1 : n
                convRecurs((ii-1)*len + ii-(ii-1)*2 : ii*len-(ii-1), 1) = ...
                    convRecurs((ii-1)*len + ii-(ii-1)*2, 1) + FXsToConvolute(:,1);
            end
            
            convRecurs(:,2:size(FXsToConvolute,2)) = fastRecursiveConvolution(...
                FXsToConvolute(:,2:end), stepSize);
            
            
            function convRecurs = fastRecursiveConvolution(FXs, stepSize)
                
                maxOrder = size(FXs,2);
                
                % needed to zero-pad
                NpointsFFT = length(FXs) + (maxOrder-1).*(length(FXs)-1);
                
                for order = maxOrder : -1 : 1
                    fftFXs(:,order) = fft(FXs(:,order), NpointsFFT);
                end
                
                convRecurs = zeros(NpointsFFT, size(FXs,2));
                convRecurs(1:size(FXs,1),1) = FXs(:,1);
                
                previousProduct = fftFXs(:,1);
                for order = 2 : maxOrder
                    product = fftFXs(:,order) .* previousProduct;
                    convRecurs(:,order) = stepSize^(order-1) * ifft(product);
                    previousProduct = product;
                end
                
            end
            
            
        end
        
        
        function [uninsured, insured] = applyPayout(groundUp, ded, cov, coI)
            
            insured = zeros(size(groundUp));
            
            insured(groundUp>ded & groundUp<cov) = coI * ...
                (groundUp((groundUp>ded & groundUp<cov)) - ded);
            insured(groundUp>=cov) = coI * max(cov - ded,0);
            % the max fx is needed to handle the impossible cases of cover<deductible
            
            uninsured = groundUp - insured;
        end
        
        
        function [TVaR, VaR] = tailValueAtRisk(CDF, alpha)
            
            % Use interp rather than find for accuracy
            
            VaRalphaToOne = find(CDF(:,2) >= alpha );
            indVaR = VaRalphaToOne(1);
            VaR = CDF(indVaR(1),1);
            
            TVaR = 1/(1-alpha) * trapz(...
                CDF(VaRalphaToOne,2), CDF(VaRalphaToOne,1));
            
        end
        
        
        function t = poissonProcess(lambda, Nyears, drawAtNyears)
            if nargin < 3; drawAtNyears = 'y'; end
            
            % 2*Nyears to be sure to reach Nyears in the final sum
            dt = random('Exponential', 1/lambda, 2*Nyears, 1);
            t = cumsum(dt);
            
            t(t>Nyears) = [];
            
            if strcmp(drawAtNyears, 'y'); t(end+1,1) = Nyears; end
        end
        
        
        function [pdfProduct, pdf2interp] = productDistribution(...
                pdf1, pdf2, z, pdf2interp)
            
            x1 = pdf1(:,1)';
            
            if isempty(pdf2interp)
                zOverX1Matrix = repelem(z, 1, numel(x1)) ./ ...
                    repelem(x1, numel(z), 1);
                zOverX1Matrix(1,:) = zOverX1Matrix(1,:) + eps;
                
                pdf2interp = interp1(pdf2(:,1), pdf2(:,2), zOverX1Matrix);
                pdf2interp(isnan(pdf2interp)) = 0;
            end
            
            absOneOverX1 = abs(1./x1);
            
            deltaX1mat = repelem(x1(2:end)-x1(1:end-1), numel(z), 1);
            
            PDF1matrix = repelem(pdf1(:,2)', numel(z), 1);
            
            integrand = PDF1matrix .* pdf2interp .* absOneOverX1;
            
            pdfProduct(:,1) = z;
            pdfProduct(:,2) = sum( ...
                (integrand(:,2:end)+integrand(:,1:end-1)) .* deltaX1mat * 0.5, 2);
        end
        
    end
    
end

