

Theory and Sensitivity Analysis for distNPVaggregateLosses
============

# Theory

## Loss distribution for a single event

The Cumulative Distribution Function (CDF) of losses for any given value of the Intensity measure (IM) of the ground motion is calculated with Eq. 1, which requires the CDF of the loss for any given damage state P(L≤l|DS_k) and the probability of the structure to be in DS_k for any given IM, P(DS_k |IM). For each DS_k (with k ranging from 0, no damage, to the number of considered damage states nDSs), P(L≤l|DS_k) is represented by a Beta distribution (Eq. 2) whose defining parameters depend on the mean loss ratio, MLR(DS_k), defined as the repair-to-reconstruction cost of the structure, and its coefficient of variation, CoV_(MLR(DS_k)). Clearly, such data should be consistent with the chosen definition of the DSs. It is worth mentioning that P(L≤l|DS_0)=1 for any l. On the other hand, P(DS_k|IM) is calculated with Eq. 3 based on the fragility curves for each DS, P(DS≥DS_k|IM). It is worth noting that P(DS≥DS_0|IM)=1 and P(DS≥DS_(nDSs+1)|IM)=0 for any given IM. Although Equations 1-3 refer to one single event, the conditioning to N_{ev}=1 is omitted from the notation, for simplicity.

![Eq 1](https://latex.codecogs.com/gif.latex?P(L&space;\leq&space;l,IM)&space;=&space;\sum_{k=0}^{nDSs}P(L&space;\leq&space;l|DS_k)P(DS_k&space;|IM)) [1]

![Eq 2](https://latex.codecogs.com/gif.latex?P(L&space;\leq&space;l|DS_k)&space;\sim&space;Beta(\frac{1-MLR(DS_k)}{CoV_MLR(DS_k)^2}-MLR(DS_k&space;),&space;\frac{\alpha(1-MLR(DS_k&space;)}{MLR(DS_k}&space;) ) [2]

![Eq 3](https://latex.codecogs.com/gif.latex?P(DS_k|IM)=P(DS&space;\geq&space;DS_k|IM)-P(DS&space;\geq&space;DS_{k&plus;1}|IM)) [3]

<!---
$$P(L \leq l,IM) = \sum_{k=0}^{nDSs}P(L \leq l|DS_k)P(DS_k |IM)$$ [1]

$$P(L \leq l|DS_k) \sim Beta(\frac{1-MLR(DS_k)}{CoV_MLR(DS_k)^2}-MLR(DS_k ), \frac{\alpha(1-MLR(DS_k )}{MLR(DS_k} )$$ [2]

$$P(DS_k|IM)=P(DS \geq DS_k|IM)-P(DS \geq DS_{k+1}|IM)$$ [3]
--->

The CDF of the seismic loss for a single event is calculated with Equations 4 and 5, by using the mean annual frequency (MAF) of occurrence of IM, which is the derivative of the hazard curve, defined in the interval [IM_{min},IM_{max}]. Clearly, the PDF p(L|N_{ev}=1) can be obtained by numerical derivation of P(L ≤ l|N_{ev}=1). It is worth mentioning that v is the seismicity rate.
By applying the insurance policy function (Eq. 2) to the loss distribution, the distribution of the uninsured losses (UL=L-PO) can be obtained. This also allows calculating the building’s EAL (related to the ground-up, uninsured or insured losses) considering the expected value of the loss distribution.

![Eq 4](https://latex.codecogs.com/gif.latex?P(L&space;\leq&space;l|N_{ev}=1)=\sum_{i=IM_{min}}^{IM_{max}}&space;P(L&space;\leq&space;l|IM_i&space;)&space;p(IM_i)) [4]

![Eq 5](https://latex.codecogs.com/gif.latex?p(IM_i)&space;=&space;\frac{-MAF(IM_i)}{\sum_{i=IM_{min}}^{IM_{max}}&space;MAF(IM_i)}) [5]

<!---
$$P(L \leq l|N_{ev}=1)=\sum_{i=IM_{min}}^{IM_{max}} P(L \leq l|IM_i ) p(IM_i)$$ [4]

$$p(IM_i) = \frac{-MAF(IM_i)}{\sum_{i=IM_{min}}^{IM_{max}} MAF(IM_i)}$$ [5]
--->

## Lifecycle loss assessment

The loss assessment module of the proposed framework allows estimating the building losses incurred over a given time horizon (T_H), which in this case is equal to the building nominal service life. The final aim is to calculate the distribution of the net present value of the aggregate losses over the time horizon, NPV(AL). Eq. 6 gives the definition of NPV(AL), where N_{ev}(T_H) is the number of seismic events that occur within T_H, [τ_i,L_i] are the time (in years) and the loss for the i-th event, and r is the financial discount rate. It is worth mentioning that, for simplicity, losses refer only to mainshocks. Moreover, this procedure implicitly assumes that the structure is upgraded to the original state every time a loss is incurred, and before the subsequent loss. To simplify the notation, the procedure in this section refers to the ground-up loss L, although it can be applied to ground-up, uninsured, or insured losses.
In this equation, N_{ev} (T_H), τ, and L are random variables and a possible and straightforward approach to obtain the distribution of NPV(AL) would be through Monte Carlo sampling. Depending on the number of samples, and considering the high number of retrofit and insurance combinations, performing the Monte Carlo sampling can require high computational resources, likely not compatible with the preliminary/conceptual design phase. For this reason, an analytical derivation of the distribution of NPV(AL) is proposed.

![Eq 6](https://latex.codecogs.com/gif.latex?NPV(AL)&space;=&space;\sum_{i=1}^{N_{ev}(T_H)}&space;NPV(L_i)=&space;\sum_{i=1}^{N_{ev}(T_H)}&space;L_i&space;\frac{1}{(1&plus;r)^{\tau_i}}&space;=&space;\sum_{i=1}^{N_{ev}(T_H)}&space;=L_i&space;NPV_i^{(1)}) [6]

<!---
$$NPV(AL) = \sum_{i=1}^{N_{ev}(T_H)} NPV(L_i)= \sum_{i=1}^{N_{ev}(T_H)} L_i \frac{1}{(1+r)^{\tau_i}} = \sum_{i=1}^{N_{ev}(T_H)} =L_i NPV_i^{(1)}$$ [6]
--->

As shown in the last member of Eq. 6, the NPV of a single loss NPV(L_i) can be expressed as the product of the loss L_i and NPV^{(1)}, which is the net present value of a unit cash flow. The probability distribution of earthquake losses L for a single event is obtained via vulnerability and hazard analysis (Eq. 4). The number of events in a given time horizon, N{ev}(T_H), is assumed to follow a Poisson distribution, and its probability mass function (PMF) is given by Eq. 7.

![Eq 7](https://latex.codecogs.com/gif.latex?P(N_{ev}&space;(T_H))&space;=&space;\frac{(\nu&space;T_H)^{N_{ev}}&space;e^{-\nu&space;T_H}}{N_{ev}!}) [7]

<!---
$$P(N_{ev} (T_H)) = \frac{(\nu T_H)^{N_{ev}} e^{-\nu T_H}}{N_{ev}!}$$ [7]
--->

The distribution of NPV^{(1)} can be evaluated based on the distribution of the arrival times \tau. As a consequence of the Poisson assumption above, the interarrival event time follows an exponential distribution with parameter v. Therefore, the arrival time of the n-th event (\tau = \sum_{i=1}^n \Delta \tau_i) follows an Erlang distribution with parameter v and numerosity n (Eq. 8). The PDF of NPV^{(1)} is obtained with Eq. 9.

![Eq 8](https://latex.codecogs.com/gif.latex?p_\tau=\frac{\nu^n&space;\tau^{n-1}}&space;{(n-1)!}&space;e^{-\nu&space;\tau}) [8]

![Eq 9](https://latex.codecogs.com/gif.latex?p_{NPV^{(1)}}&space;(n)&space;=&space;\frac{1}{ln(1&plus;r)NPV^{(1)}}&space;p_{\tau}(log_{1&plus;r}\frac{1}{NPV^{(1)}})=\frac{1}{ln(1&plus;r)NPV^{(1)}}&space;\frac{\nu^{n}(ln(1&plus;r)NPV^{(1)})^{n-1}}{(n-1)!}e^{-\nu&space;log_{1&plus;r}\frac{1}{NPV^{(1)}}}) [9]

<!---
$$p_\tau=\frac{\nu^n \tau^{n-1}} {(n-1)!} e^{-\nu \tau}$$ [8]

$$p_{NPV^{(1)}} (n) = \frac{1}{ln(1+r)NPV^{(1)}}  p_{\tau}(log_{1+r}\frac{1}{NPV^{(1)}})=\frac{1}{ln(1+r)NPV^{(1)}} \frac{\nu^{n}(ln(1+r)NPV^{(1)})^{n-1}}{(n-1)!}e^{-\nu log_{1+r}\frac{1}{NPV^{(1)}}}$$ [9]
--->

By definition of product distribution for independent random variables, the PDF of NPV(L) is obtained with Eq. 10, where p_L and p_{NPV^{(1)}} are the PDFs of L and NPV^{(1)}, respectively. Since p_L is generally not available in analytical form, this integral must be solved numerically . Moreover, p_L generally shows an asymptotic behaviour for L=0, and therefore its perfect numerical representation is theoretically impossible. Practically, this introduces an error which can be arbitrarily reduced by reducing the sampling step for p_L (Eq. 4). A detailed discussion of this issue is provided in ref .

![Eq 10](https://latex.codecogs.com/gif.latex?p_{NPV(L)}&space;(n)&space;=\int_{-\infty}^{\infty}p_{NPV^{(1)}}(n)&space;p_{L}(\frac{NPV(L)}{npv^{(1)}})&space;\frac{1}{|npv^{(1)}|}&space;dnpv^{(1)}) [10]

<!---
$$p_{NPV(L)} (n) =\int_{-\infty}^{\infty}p_{NPV^{(1)}}(n) p_{L}(\frac{NPV(L)}{npv^{(1)}}) \frac{1}{|npv^{(1)}|} dnpv^{(1)}$$ [10]
--->

Finally, the PDF of NPV(AL) is given by Eq. 11, p_{NPV(L)}^{N_ev^*} is the convolution of the functions p_NPV(L)(N_{ev}) with numerosity n=1…N_{ev}, which can be obtained recursively. It is important repeating that the summation in Eq. 13 theoretically extends to infinity. However, it can be truncated by considering the value of N_{ev} such that P(N_{ev}) is sufficiently close to zero (e.g. <0.001).

![Eq 11](https://latex.codecogs.com/gif.latex?p_{NPV(AL)}&space;=&space;\sum_{N_{ev}}&space;P(N_{ev})p_{NPV(L)}^{N_{ev}^{*}}) [11]

<!---
$$p_{NPV(AL)} = \sum_{N_{ev}} P(N_{ev})p_{NPV(L)}^{N_{ev}^{*}}$$ [11]
--->

# Demo

The calculations shown in this section refer to the file demo/testNPVlosses.m.

The constructor of this class accepts inputs organised in a structure array ('options') with five macro fields ('General', 'Vulnerability', 'Hazard', 'Insurance', 'Setup'). Each macro field has further micro fields shown in the code snippet below. The shown values are considered as defaults: if the option structure is missing any field, those default values will be used.

```
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

options.Insurance.deductible = 0;
options.Insurance.cover = 0;
options.Insurance.coinsurance = 1;

options.Setup.NlossSamples = 201;
options.Setup.IMstep = 0.005; % [g]
options.Setup.MCsamples = 10000;

obj = distNPVaggregateLosses(options);
```

The class provides methods to derive :
- the probability of having n earthquakes within the time horizon. The largest considered n is the first with a PMF lower than 0.001;
- the PDF of the interarrival times of earthquakes;
- the PDF of the arrival time of the n-th earthquake.

```
obj = obj.getPMFnumberEvents;
obj = obj.getPDFinterarrivalTime;
obj = obj.getPDFarrivalTime;
```

![PMF Nevents](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/PMFNevents.png?raw=true)
![PDF arrival](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/PDFarrival.png?raw=true)

The loss distribution (given one event) is calculated with the method "getLossDistribution". This method automatically:
- represents fragility functions according to the provided medians and dispersions;
- calculates the distribution of loss given damage state;
- calculates the distribution of loss given hazard intensity;
- calculates the distribution of loss given one event;
- if insurance is considered, it calculates the distribution of uninsured losses given one event;
- calculates the expected annual losses (both ground up and uninsured).

```
obj = obj.getLossDistribution;
```

![fragilities](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/fragilities.png?raw=true)
![CDFlosssGivenDS](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/CDFlossGivenDS.png?raw=true)
![lossDistributionPercentiles](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/lossDistributionPercentiles.png?raw=true)
![lossGivenOne](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/lossGivenOne.png?raw=true)

The method getPDFlossNPV calculates the distribution p_{NPV(L)}(n). The method also calculates the PDF of the unit cash flow related to the n-th event.

```
obj = obj.getPDFlossNPV;
```

![PDF npv1](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/PDFnpv1.png?raw=true)
![CDF npvL](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/CDFnpvL.png?raw=true)

Finally, the method getAggregateLossNPVdist calculates p_{NPV(AL)}.

```
obj = obj.getAggregateLossNPVdist;
```

![NPVaggLossDist](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/NPVaggLossDist.png?raw=true)

# Sensitivity Analysis

The calculations shown in this section refer to the file demo/sensitivityNPVlosses.m. Apart from an interest rate equal to zero, the default values of this class are used.

The first sensitivity analysis involves the number of samples defining p_{L}. The selected values of the samples are [11 101 201 501 1001 2001]; The figure below shows how 201 points provide a satisfactory trade-off between accuracy and efficiency.

![sensLossSamples](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/sensLossSamples.png?raw=true)

The second sensitivity analysis relates to the number of samples used in the Monte Carlo approach. p_{NPV(AL)} is calculated based on [100 1000 10000 100000 500000 1000000] samples. The figure below shows that the accuracy stops improving at 500000 samples.

![sensMCsamples](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/sensMCsamples.png?raw=true)

Finally, using 301 loss samples and 20000 Monte Carlo samples, p_{NPV(AL)} is calculated with both approaches considering the interest rates [0 0.005 0.01 0.02 0.04 0.08]%; As shown in the figure below, there is a perfect match between analytical and Monte Carlo if the interest rate is equal to zero. A slight mismatch appears for interest rates between 0.5% and 2%. Finally, the mismatch is negligible for interest rates greater than 4%.

![sensIntRate](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/sensIntRate.gif)
