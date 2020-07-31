# mumutron
Code for the mumutron project

Mumutron is a project of a compact electron-positron collider developed in the Budker Institute for Nuclear Physics. Details are in the pdf.

This repository contains a set of ROOT CERN (https://root.cern.ch) scripts, that calculate the energy measurement precision, energy in the center-of-mass 
reference frame distribution measurement precision, and optimal scanning points for the mumutron project in BINP.

The "section" folder contains code for calculating the cross-section of muon-antimuon pair production in the electron-positron collisions at the threshold.
It accounts for inital state radiation, interaction in the final state, and energy distribution in the beam. The calculation invovles numerical integration
of singular functions. All the integration is handled by GSL library, particularly, I have employed Gauss-Kronrod method.

The "resolution" folder is a collection of scripts that constitute a simulation of the mentioned process and assessment of the distribution of energy in the 
center-of-mass reference frame measurement precision. To simulate the process I used a rejection Monte-Carlo method. The cross-section is proportional to 
the probability density function, making them interchangeable for my purpose. This simulation accounts for uncertainties in beam energy, angle of the beams, and
coordinate determination by the detector. I have also implemented a simple model of multiple scattering on the berrilium foil. This is done as a simple rotation
of the particle direction. The result of the siulation is a distribution of muons on the coordinate detector. We look at the distance from the (0, 0) coordinate
of the detector. The strategy is to simulate 16 histograms with given energy in the c.o.m. ref. frame, and interpolate them. Then, when the experiment is running
we will be able to fit experimental distributions with the results of the simulation.

The "report" folder is a brief presentation of this work.
