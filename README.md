distNPVaggregateLosses
============
[![GitHub Stars](https://img.shields.io/github/stars/robgen/distNPVaggregateLosses.svg)](https://github.com/robgen/distNPVaggregateLosses/stargazers) [![GitHub Issues](https://img.shields.io/github/issues/robgen/distNPVaggregateLosses.svg)](https://github.com/robgen/distNPVaggregateLosses/issues) [![Current Version](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/robgen/distNPVaggregateLosses)

This is a Matlab class that allows to compute the distribution of the Net Present Value (NPV) of the aggregate seismic economic losses (ground up, uninsured or insured) over a given time horizon, and for a given financial discount rate. All the steps are modelled analytically. A method performing a MonteCarlo numerical approach is also provided.

![PDF aggregate loss](https://github.com/robgen/distNPVaggregateLosses/blob/main/test/PDFaggLoss.png?raw=true)

---
## Buy me a coffee

Whether you use this project, have learned something from it, or just like it, please consider supporting it with a small donation, so I can dedicate more time on open-source projects like this :)

<a href="http://paypal.me/robgen" target="_blank"><img src="https://www.paypalobjects.com/webstatic/mktg/logo/pp_cc_mark_74x46.jpg" alt="Paypal" style="height: auto !important;width: auto !important;" ></a>

---

## Features
- Probability of n events occurring during the selected time horizon
- Arrival times of the n-th earthquake according to a Poisson process
- Probability distribution of the earthquake losses (ground up) given one event of any intensity
- Probability distributions for insured and uninsured losses
- Probability distribution of the Net Present Value (NPV) of losses
- Probability distribution of the NPV of aggregate losses during the selected time horizon

---

## Setup
Clone this repo to any folder in your computer. Add this folder to the Matlab path and you are ready to go (or just cd to this folder). This class does not need external dependencies, apart from the built-in Matlab functions.

---

## Usage
A full demo of this class is given in the file testNPVlosses.m

---

## License
This project is licensed under the terms of the **Creative Commons “Attribution-No Derivatives 4.0 International”** license. This software is supplied "AS IS" without any warranties and support. The Author assumes no responsibility' or liability for the use of the software. The Author reserves the right to make changes in the software without notification.
