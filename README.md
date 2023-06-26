[![GitHub (pre-)release](https://img.shields.io/github/release/vasilsaroka/QEskillbox/all.svg)](https://github.com/vasilsaroka/QEskillbox/releases)
[![Github All Releases](https://img.shields.io/github/downloads/vasilsaroka/QEskillbox/total.svg)](https://github.com/vasilsaroka/QEskillbox/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Support QEskillbox](https://img.shields.io/static/v1?label=support&message=5$&color=green&style=flat&logo=paypal)](https://paypal.me/vasilsaroka?locale.x=en_GB)

# QEskillbox
Automate [Quantum Espresso](https://www.quantum-espresso.org/) routines

## Differential Evolution
### Python
This repository contains Python code for Quantum Espresso parameters optimization with Differential Evolution global optimization method [1].
An excellent introduction into this method together with a python code can be found in the [tutorial](https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/#) by Pablo Rodriguez-Mier.
The DE algorithm here is rewritten from scratch in generation mixing form. It allows one to save some memory by storing the parameters space vectors from different generation in the same array. Several tests have shown that this does not affect the performance of the algorithm.
[1] Price, K., Storn, R. (1997). Differential Evolution – A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces. *Journal of Global Optimization*, *11*, 341–359. [https://doi.org/10.1023/A:1008202821328](https://doi.org/10.1023/A:1008202821328)

### Bash
