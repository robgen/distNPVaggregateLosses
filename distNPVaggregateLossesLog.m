classdef distNPVaggregateLosses
    %distNPVaggregateLosses Summary of this class goes here
    %   
    
    properties
        
        parameters
        
        TAUdef
        IMdef
        LOSSdef
        fragilities
        
        NmaxEvents
        PDFinterarrivalTime
        PDFarrivalTime
        PMFnumberEvents
        PMFds
        CDFlossGivenDS
        CDFlossGivenIM
        CDFlossGivenOneEvent
        PDFlossGivenOneEvent
        MAFim
        hazCurveResampled
        lastStepIM
        PDFunitCashFlowNPV
        PDFlossNPV
        PDFaggLossNPVGivenNevents
        PDFaggLossNPV
        CDFaggLossNPV
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
            
            if nargin < 2 
                if isempty(self.NmaxEvents)
                    NmaxEvents = 10;
                else
                    NmaxEvents = self.NmaxEvents;
                end
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
            
            n = -1;
            pNevents = 1;
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
        
        
        function self = getCDFloss(self)
            
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
            
            warning('put back the numerical gradient if Loss has uniform points')
            self.PDFlossGivenOneEvent = self.numericalDerivative(...
                self.CDFlossGivenOneEvent);
            %self.PDFlossGivenOneEvent = self.numericalGradient(...
            %    self.CDFlossGivenOneEvent);
        end
        
        
        function self = getPDFlossNPV(self)
            
            self = getPDFunitCashFlowNPV(self);
            
            npvL = self.LOSSdef;
            npv1 = self.PDFunitCashFlowNPV(:,1)';
            
            NPVratioMatrix = repmat(npvL, 1, numel(npv1)) ./ ...
                repmat(npv1, numel(npvL), 1);
            %NPVratioMatrix(1,:) = NPVratioMatrix(1,:) + eps;
            
            PDFlossInterp = interp1(self.PDFlossGivenOneEvent(:,1), ...
                self.PDFlossGivenOneEvent(:,2), NPVratioMatrix);
            
            % the interpolation above will give NaN because
            % PDFlossGivenOneEvent is not defined
            PDFlossInterp(NPVratioMatrix>1) = 0;
            
            absOneOverNPV1 = abs(1./npv1);
            
            for n = self.NmaxEvents : -1 : 1
                PDFnpv1matrix = repmat(...
                    self.PDFunitCashFlowNPV(:,n+1)', numel(npvL), 1);
                
                integrand = PDFnpv1matrix .* PDFlossInterp .* absOneOverNPV1;
                
                % maybe a performance improvement?
                for row = size(integrand,1) : -1 : 1
                    self.PDFlossNPV(row,n+1) = trapz(npv1, integrand(row,:));
                end
                
%                 row = 1;
%                 figure; hold on
%                 plot(npv1, integrand(row,:), 'LineWidth', 2)
%                 plot(npv1, PDFnpv1matrix(row,:), 'LineWidth', 2)
%                 plot(npv1, PDFlossInterp(row,:), 'LineWidth', 2)
%                 plot(npv1, absOneOverNPV1, 'LineWidth', 2)
%                 legend('integrand', 'PDFnpv1matrix', 'PDFlossInterp', 'absOneOverNPV1')
%                 title(sprintf('n=%d, row=%d', n, row))
%                 axis([0 1 0 1000])
            end
            
            self.PDFlossNPV(:,1) = npvL;
            
        end
        
        
        function self = getPDFaggregateLossNPV(self)
            
            self.PDFaggLossNPVGivenNevents = self.recursiveConvolution(...
                self.PDFlossNPV);
            
            pNevMatrix = self.PMFnumberEvents(:,2)' + ...
                zeros(size(self.PDFaggLossNPVGivenNevents,1),1);
            
            self.PDFaggLossNPV(:,1) = self.PDFaggLossNPVGivenNevents(:,1);
            self.PDFaggLossNPV(:,2) = sum(pNevMatrix(:,2:end) .* ...
                self.PDFaggLossNPVGivenNevents(:,2:end), 2);
            % note: P(L>l|0events) = 0, for each l
            
            self.CDFaggLossNPV = self.numericalIntegral(self.PDFaggLossNPV);
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
            % P(L<l|DS0) is a step function centered at zero
            self.CDFlossGivenDS(2:end,1) = 1;
            
        end
        
        
        function self = getCDFlossGivenIM(self)
            
            %%% needs performance improvement
            for im = numel(self.IMdef) : -1 : 1
                CDFlossDSmatrix(im,:,:) = self.CDFlossGivenDS;
            end
            
            for l = numel(self.LOSSdef):-1:1
                probDS3d(:,l,:) = self.PMFds;
            end
            %%% needs performance improvement
            
            self.CDFlossGivenIM = sum( CDFlossDSmatrix .* probDS3d, 3 );
        end
        
        
        function self = getMAFim(self)
            
            if max(self.IMdef) < max(self.parameters.Hazard.hazCurve(:,1))
                error('Extend the IM definition of CDFlossIM')
            end
            
            extendedHazCurve = [0 self.parameters.Hazard.faultRate; ...
                self.parameters.Hazard.hazCurve];
            
            self.MAFim = self.numericalDerivative(extendedHazCurve);
            self.MAFim(:,2) = abs(self.MAFim(:,2));
            
            [self.MAFim, self.lastStepIM] = self.resampleCurve(...
                self.MAFim, self.IMdef);
        end
        
        
        function self = getPDFunitCashFlowNPV(self)
            
            r = self.parameters.General.intRate;
            v = self.parameters.Hazard.faultRate;
            
            warning('put npv definition back if Loss has uniform points')
            %npv = linspace(0, max(self.LOSSdef), 501) + eps;
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
        
        
        %%% Internal stuff
        function self = setVariableRanges(self)
            
            self.TAUdef = linspace(0, ...
                8*self.parameters.General.timeHorizon, 100);
            
            self.IMdef = (eps:0.005: ...
                3*max(self.parameters.Vulnerability.fragMedians))';
            
            %self.LOSSdef = linspace(0,1,1001)' + 10*eps;
            self.LOSSdef = [eps; logspace(-15,0,1000)']; % linspace(0.1,1,100)'];
            % MUST BE equally spaced. MUST start from eps
        end
        
        
        function self = setAllParameters(self, options)
            % setAllParameters deals with the optional parameters
            
            % build basic parameters
            macroFieldsPar = {'General', 'Vulnerability', 'Hazard'};
            
            microFieldsPar{1} = {'timeHorizon', 'intRate' };
            microFieldsParVals{1} = {50, 0.02};
            
            microFieldsPar{2} = {'fragMedians', 'fragDispersions', 'damgeLossRatios', 'CoVdlr' };
            microFieldsParVals{2} = {[0.05 0.35 0.65 1.25], [0.2 0.2 0.2 0.2], [0 2 10 43.5 95]/100, [ 0 1 0.4 0.3 0.05 ]};
            
            microFieldsPar{3} = {'hazCurve', 'faultRate' };
            microFieldsParVals{3} = {[0.166427989012818,0.0332146240000000;0.217582434028613,0.0198850450000000;0.258529430931683,0.0138629440000000;0.303930770035728,0.00988592600000000;0.354443451456181,0.00713349900000000;0.412206673094016,0.00496922700000000;0.565248464301760,0.00210721000000000;0.695119133694674,0.00102586600000000;0.846507595616605,0.000404054000000000] , 0.05};
            
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
                warning('Please make sure that "curve" contains equally spaced points')
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
            
            FXsToConvoluteResamp(:,1) = linspace(eps,1-eps,1001)';
            n = size(FXsToConvolute,2)-1;
            for k = n : -1 : 1
                FXsToConvoluteResamp(:,k+1) = interp1(...
                    FXsToConvolute(:,1), FXsToConvolute(:,k+1), ...
                    FXsToConvoluteResamp(:,1));
            end
            
            stepSize = FXsToConvoluteResamp(2,1) - FXsToConvoluteResamp(1,1);
            
            convRecursCell = cell(n,1);
            convRecursCell{1}(:,1) = FXsToConvoluteResamp(:,1);
            convRecursCell{1}(:,2) = FXsToConvoluteResamp(:,2);
            for ii = 1 : n-1
                % X values
                [dummy1, dummy2] = meshgrid(...
                    convRecursCell{ii}(:,1), FXsToConvoluteResamp(:,1));
                convRecursCell{ii+1}(:,1) = uniquetol(dummy1+dummy2, 10*eps);
                % this line should give a vector with length length(u)+length(v)-1,
                % where u and v are the operands of conv. This does not happen if i
                % don't resample with equally spaced points. On top of this, I need a
                % tolerance because I can never compare two "real" numbers, I can only
                % compare them close
                
                % Y values
                convRecursCell{ii+1}(:,2) = stepSize * conv(...
                    convRecursCell{ii}(:,2), FXsToConvoluteResamp(:,ii+2));
                
                clear dummy1 dummy2
            end
            
            % create single matrix (more practical)
            convRecurs(:,1) = convRecursCell{end}(:,1);
            for ii = n : -1 : 1
                convRecurs(1:size(convRecursCell{ii},1),ii+1) = ...
                    convRecursCell{ii}(:,2);
            end
        end
        
        
    end
end

