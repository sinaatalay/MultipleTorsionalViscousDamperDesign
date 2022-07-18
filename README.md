# Multiple Torsional Viscous Damper Design [![License](https://img.shields.io/github/license/sinaatalay/MultipleTorsionalViscousDamperDesign.svg)](https://github.com/sinaatalay/MultipleTorsionalViscousDamperDesign/blob/main/LICENSE)
This repository contains a project for [Bogazici University Mechanical Engineering Department](https://www.me.boun.edu.tr/)'s **Mechanical Vibrations** (ME 425) class. My teammates and I solved and optimized the system shown below.

**Team Members:** Sina Atalay, Cem Geçgel, Mustafa Çağatay Sipahioğlu

## The System

<p align="center">
	<picture>
	  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/sinaatalay/MultipleTorsionalViscousDamperDesign/blob/main/figures/TheSystemDarkMode.png?raw=true">
	  <source media="(prefers-color-scheme: light)" srcset="https://github.com/sinaatalay/MultipleTorsionalViscousDamperDesign/blob/main/figures/TheSystem.png?raw=true">
	  <img alt="Schematic" src="https://github.com/sinaatalay/DynamometerSimulation/blob/main/figures/Schematic.png?raw=true">
	</picture>
</p>

In the figure, an (n+2)-degrees-of-freedom system that will be designed is shown, where $I_i=100/n$, $k_i=25n$. Four problems are solved to design the system:

1.  What are the natural frequencies and mode shapes of the n-degrees-of-freedom disk-torsional spring chain without the viscous torsional dampers for a given $n$?
2.  Let the system without the viscous torsional dampers to be base excited. Find transmissibility for the last disk ($\Theta_n/\Phi$), assuming $\vec{\theta} \left(t\right)=\vec{\Theta} e^{i\omega t}\text{, }\vec{\varphi} \left(t\right)=\Phi e^{i\omega t}$, and plot the absolute value of the transmissibility in log-log axes for $0\le \omega \le 1\ldotp 5\omega_{n,\max }$.
3.  For the (n+2)-degrees-of-freedom system shown in the figure (viscous dampers are attached to the first and last disks), find the optimum $c_{a1}$ and $c_{a2}$ values that minimize the peak transmissibility for the last disk. $I_{a1}$ and $I_{a2}$ are given and between 0.1 and 0.3.
4.  For a given $I_{a1}+I_{a2}=\mu$ value, find the optimum $I_{a1}$, $I_{a2}$, $c_{a1}$, $c_{a2}$ values and location of the torsional dampers that minimize the peak transmissiblity for the last disk.
