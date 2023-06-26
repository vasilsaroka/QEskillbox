[![GitHub (pre-)release](https://img.shields.io/github/release/vasilsaroka/QEskillbox/all.svg)](https://github.com/vasilsaroka/QEskillbox/releases)
[![Github All Releases](https://img.shields.io/github/downloads/vasilsaroka/QEskillbox/total.svg)](https://github.com/vasilsaroka/QEskillbox/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Support QEskillbox](https://img.shields.io/static/v1?label=support&message=5$&color=green&style=flat&logo=paypal)](https://paypal.me/vasilsaroka?locale.x=en_GB)

# QEskillbox
Automate [Quantum Espresso](https://www.quantum-espresso.org/)(QE) routines

## Differential Evolution
### Python
This repository contains Python code for Quantum Espresso parameters optimization with Differential Evolution (DE) global optimization method [1].
An introduction into this method together with a python code can be found in the excellent [tutorial](https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/#) by Pablo Rodriguez-Mier.

The DE algorithm here is rewritten from scratch in what we call a generation mixing form. It allows one to save some memory by storing the parameters space vectors from different generation in the same array. Several tests have shown that this does not affect the performance of the algorithm.

The python file contains the main functions *QEGenrun* and *QEDE*, the cell for a test run of the code on bulk Si followed by a few other cells that are used for the output results analysis such as visualizing and animating the DE convergence. The *QEDE* can perform optimization for an arbitrary number of parameters including integer ones, i.e. those for which only integer part is meaningful. The integer parameters are counted from the end of a parameter space vector and therefore must be placed accordingly. A template file for bulk Si with parameter placeholders starting by convention from '@' is provided.  

 - [1] Price, K., Storn, R. (1997). Differential Evolution – A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces. *Journal of Global Optimization*, **11**, 341–359. [https://doi.org/10.1023/A:1008202821328](https://doi.org/10.1023/A:1008202821328)

### Bash
This repository also contains a bash script that creates the input QE template file and a python script and then runs the python script. After the python run the temporarily created files are deleted leaving only the results in QEDE*.out and QEDE*.png output files. Bulk Si is taken as an example with the three optimization parameters: ecutwfc, celldm and kpt. For kpt parameter only its integer part is meaningful. 
