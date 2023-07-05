[![GitHub (pre-)release](https://img.shields.io/github/release/vasilsaroka/QEskillbox/all.svg)](https://github.com/vasilsaroka/QEskillbox/releases)
[![Github All Releases](https://img.shields.io/github/downloads/vasilsaroka/QEskillbox/total.svg)](https://github.com/vasilsaroka/QEskillbox/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Support QEskillbox](https://img.shields.io/static/v1?label=support&message=5$&color=green&style=flat&logo=paypal)](https://paypal.me/vasilsaroka?locale.x=en_GB)

# QEskillbox
Automate [Quantum Espresso](https://www.quantum-espresso.org/)(QE) routines

## Requirements
The code was developed and tested under
 - OS: Linux Mint 21.1 Vera Xfce
 - Python 3.10, with modules *os*, *subprocesses*, *re*, *numpy*, *time*, *datetime*, *matplotlib.pyplot*.
 - PWSCF v.6.7MaX that is a part of [Quantum ESPRESSO](https://www.quantum-espresso.org/) suite.

## Preliminaries
Under Linux OS:
- Install QE by typing in the terminal: ``sudo apt-get install quantum-espresso`` (recommended) or see [QE](https://www.quantum-espresso.org/login/) web page.
- Install Python distribution by typing in the terminal: ``sudo apt-get install conda`` (recommended) or see [Anaconda Distribution](https://www.anaconda.com/download/#linux) guidelines.
- Pseudo potentials can be taken from [Pseudo Dojo](http://www.pseudo-dojo.org/).

## Differential Evolution
This repository contains a Python code for Quantum Espresso parameters optimization with Differential Evolution (DE) global optimization method [1].
A general introduction into DE method together with a python realization can be found in the excellent [tutorial](https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/#) by Pablo Rodriguez-Mier. Here, the DE method is adapted for the *geometry optimization* of the bulk Si crystal. In effect, this allows one to perform *0_Si_bulk* tasks *0_cutoff*, *1_alat* and *2_kpoints* from [DFT_QE_beginner_tutorial](https://github.com/Crivella/DFT_QE_beginner_tutorial) or [DFT-basics](https://github.com/cdhogan/DFT-basics) in a hassle-free single run.

- [1] Price, K., Storn, R. (1997). Differential Evolution – A Simple and Efficient Heuristic for Global Optimization over Continuous Spaces. *Journal of Global Optimization*, **11**, 341–359. [https://doi.org/10.1023/A:1008202821328](https://doi.org/10.1023/A:1008202821328)


## Demo
### Bash script DEopt
- Put Bash script ``DEopt`` and *Si.pz-vbc.UPF* pseudo potential into the same folder, i.e. the working folder.
- Open terminal in this folder and run the Bash script in terminal as ``./DEopt``.
- Find in the working folder the QE outdir ``./tmp``.
- Find in the working folder the results of the optimization and convergence plot in *QEDE_year-month-day_hours_minutes_seconds.out* and *QEDE_year-month-day_hours_minutes_seconds.png*.

### Python file DEopt.py
- Put Python file ``DEopt.py`` with QE input file *Si.scf-template.in* and pseudo potential *Si.pz-vbc.UPF* into the same folder, i.e. working folder
- Open Python file and evaluate it cell after cell sequentially.
- Find in the working folder the QE outdir ``./tmp``.
- Find in the working folder the results of optimization and convergence plot in *QEDE_year-month-day_hours_minutes_seconds.out* and *QEDE_year-month-day_hours_minutes_seconds.png*.
- Find in the working folder the animation for convergence in *animation_Si.gif* and *video_Si.mp4*.

![QEDE_2023-06-26_17_52_08](https://github.com/vasilsaroka/QEskillbox/assets/49445896/6231b35c-0742-4340-93ea-b036dc2beb2f)

![animation_Si](https://github.com/vasilsaroka/QEskillbox/assets/49445896/bd93e31e-657b-4ebc-9297-88b84fbc5efb)


## Supporting this project
   We believe everyone deserves access to knowledge that is grounded in science and integrity. That is why we keep our code open for all users, regardless of where they live or what they can afford to pay. This means more people can be better educated and inspired to make an impact on the global wellbeing. We have no shareholders or billionaire owner, meaning only your donations power our work and ensure it can remain open for all. Every contribution, however big or small, makes a real difference for QEskillbox future. If you find it useful, consider [supporting the project](https://paypal.me/vasilsaroka?locale.x=en_GB)
. 
