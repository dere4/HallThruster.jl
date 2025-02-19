---
title: 'HallThruster.jl: A Julia package for 1D fluid Hall Thruster discharges'
tags:

  - Plasma physics
  - Low temperature magnetized plasma
  - CFD
  - Hall Thruster
  - Electric Propulsion
  - Julia
  authors:
  - name: Thomas Marks
    orcid: 0000-0000-0000-0000
    affiliation: 1 
  - name: Paul Schedler
    affiliation: 2
 - name: Plasmadynamics and Electric Propulsion Laboratory, University of Michigan, Ann Arbor, USA
   index: 1
 - name: ETH Zurich, Zurich, Switzerland
   index: 2
   date: 22 April 2022
   bibliography: paper.bib

# Summary

Hall thrusters are a widely-used class of spacecraft electric propulsion device. They are annular cross-field devices with an applied axial electric field perpendicular to a magnetic field generated by coils on the outer walls and centerline of the thruster. This results in electrons being trapped along the magnetic field lines, creating a strong azimuthal Hall current, ionizing injected gas, usually xenon. The ions are then accelerated out of the channel by the electric field, which generates thrust. Hall thrusters offer a specific impulse on the order of 1500-3000s, and can achieve efficiencies higher than 50%. An overview is available in [@goebelkatzhall] and [@physicsmodeling]. Hall Thrusters are commonly used for in-space propulsion for commercial communications and surveillance satellites as well as deep space missions. 

Electron transport across the magnetic field lines is non-classical and poorly understood in Hall thrusters, but getting it right is critical to making simulations match experiment. While high-fidelity models are useful in trying to model the precise details of this electron transport, they are often unsuitable for whole device modeling. 1D codes like HallThruster.jl capture enough of the important dynamics in Hall thrusters to be partially predictive and run quickly enough to be used in model discovery and optimization work, making them a useful tool in the electric propulsion engineer's toolkit.


# Statement of need

`HallThruster.jl` is a 1D fluid Hall Thruster discharge model. Julia enables the efficient
implementation of a flexible code by employing its multi-dispatch architecture. The "anomalous"
electron transport in Hall Thrusters is an unsolved problem and currently prohibits
predictive Hall Thruster simulations. As full time and length resolving particle simulations of 
low temperature magnetized plasmas at the scales required for Hall Thrusters are prohibitively
expensive. Another option is to infer anomalous collision frequency models using large datasets. 
This is where HallThruster.jl is applied to in the scientific community. 

`HallThruster.jl` can be used in teaching due to the simple interface and extensive documentation. 
It can provide students with a first order understanding of the governing physics in Hall Thrusters. 
The Landmark study [@landmarkplasma] provides a benchmark for 1D fluid Hall Thruster discharges, which 
is reinforced by a public repository like `HallThruster.jl` recreating its results. 

`HallThruster.jl` runs very quickly. State of the art two-dimensional fluid simulations typically take on the order of 10+ hours to run on a typical desktop, while a one-dimensional simulation in HallThruster.jl takes anywhere from a few seconds to a few minutes, depending on the grid resolution. This makes it ideal for applications in which the user would want to run many simulations sequentially or simultaneously, such as parameter estimation, uncertainty quantification, and surrogate optimization. Hall thruster dynamics are suitably one-dimensional that the 1D code is often quite a good approximation of a 2D code.

Despite their usefulness, there are no fully open-source 1D Hall thruster codes at the moment. This hampers the ability of the community to reproduce computational results obtained in such codes.  HallThruster.jl hopes to partially answer this need and pave the way for a more open electric propulsion simulation community.

# Physics model

In \autoref{fig:domain}, we depict the one-dimensional simulation domain. A radial magnetic field is applied in the thruster channel, crossed with an axial electric field between anode and cathode. Electrons drift primarily in the $\hat{\theta} = \hat{z} \times \hat{r}$ direction.

![1D simulation domain of HallThruster.jl\label{fig:domain}](C:\Users\thoma\HallThruster.jl\paper\figure1D_code.PNG)



In HallThruster.jl, we treat all species as fluids, with different models for neutrals, ions, and electrons. Quasineutrality ($n_e = \sum_{Z=1}^N Z \;n_{iZ}$, where $n_e$ is the plasma density and $n_{iZ}$ is the number density of atoms with charge state $Z$) is assumed in the entire domain. Neutrals are assumed to have constant temperature and density, so only the continuity equation is solved. 

\begin{equation}
    \frac{\partial n_n}{\partial t} + \frac{\partial}{\partial z} (n_n u_n) = -n_e n_n \sum{r} k_r\left(\frac{3}{2} T_{eV}\right)
\label{eq:neutrals}
\end{equation}

Here, $n_n$ is the neutral number density, $u_n$ is the neutral velocity, $T_{eV}$ is the electron temperature in electron-volts, and $k_r$ is the rate coefficient of ionization reaction $r$ as a function of the mean electron energy ($\epsilon = \frac{3}{2} T_{eV}$). The source term in \autoref{eq:neutrals} is a loss term corresponding to the rate at which neutrals are ionized. Recombination in the domain is neglected in this model, and currently, charge exchange collisions are not modeled. 

Ions are approximated as isothermal, so both the ion continuity and momentum equation are solved. 

\begin{equation}
    \frac{\partial}{\partial t} n_{iZ} + \frac{\partial}{\partial z} (n_{iZ} u_{iZ}) = \dot{n}_{iZ}
\label{eq:ioncontinuity}
\end{equation}

\begin{equation}
    \frac{\partial}{\partial t} (n_{iZ} u_{iZ}) + \frac{\partial}{\partial z} (n_{iZ} u_{iZ}^2 + \frac{p_{iZ}}{m_i}) = -\frac{Z e}{m_i} n_{iZ} \frac{\partial \phi}{\partial z} + \dot{n}_{iZ} u_n
\label{eq:ionmomentum}
\end{equation}

where $Z$ is the ion charge state, $\phi $  is the plasma potential,  $e$ is the fundamental charge, $1.6\times 10^{-19}$ Coulombs, $m_i$ is the ion mass in kilograms, $n_{iZ}$ and $u_{iZ}$ are the number density and velocity, respectively, of ions with charge $Z$, $p_{iZ} = n_{iZ} k_B T_i$ is the pressure of ions with charge $Z$, $T_i$ is the constant (user-configurable) ion temperature in Kelvins, and $\dot{n}_{iZ}$ is the rate of production of ions with charge $Z$  Up to three ion fluids with different charge states are supported in `HallThruster.jl`. 

Electrons are assumed massless and their inertia is neglected. Additionally, the electric field is assumed purely axial $\left(-\nabla\phi = -\frac{\partial \phi}{\partial z} \hat{z}\right) $and the magnetic field purely radial $\left(\vec{B} = B_r \hat{r}\right)$. This leads to the 1D generalized Ohm's law in the for the axial (cross-field) electron current $j_{ez}$:

\begin{equation}
    j_{ez} = \frac{e^2 n_e}{m_e \nu_e}\frac{1}{1 + \Omega_e^2}\left(-\frac{\partial \phi}{\partial z} + \frac{1}{e n_e}\frac{\partial p_e}{\partial z}\right)
\label{eq:eleccurrentaxial}
\end{equation}

The electron current in the azimuthal direction (the Hall current) is:

\begin{equation}
j_{e\theta} = \Omega_e j_{ez}
\label{eq:eleccurrentazimuthal}
\end{equation}

Here, $\nu_e$ is the total electron momentum transfer collision frequency, $m_e$ is the electron mass, and $\Omega_e = \omega_{ce}/\nu_e = e |B| / m_e \nu_e$ is the Hall parameter, or the ratio of electron cyclotron frequency to collision frequency. This is a measure of how magnetized the electrons are and is typically greater than 1 in a Hall thruster discharge. $\nu_{e}$ is the sum of several component collision frequencies, $\nu_{e} = \nu_{ei} + \nu_{en} + \nu_{wall} + \nu_{AN}$, where $\nu_{en}$ is the electron-neutral collision frequency, $\nu_{ei}$ is the electron-ion collision frequency, $\nu_{wall}$ is the electron-wall collision frequency, and $\nu_{AN}$ is the electron anomalous collision frequency. Different models are available for each of these collision frequencies, and `HallThruster.jl` is designed in such a way that the user may easily implement new models.

Current conservation is given by:

\begin{equation}
   \nabla \cdot j_e + \nabla \cdot j_i = 0
\label{eq:current_cons}
\end{equation}

Taking the divergence of \autoref{eq:eleccurrentaxial}, and incorporating current conservation, we arrive at an elliptic equation which can be solved for the electrostatic potential:

\begin{equation}
    \frac{\partial}{\partial z}\left(\mu_{\perp} n_e \frac{\partial\phi}{\partial z}\right) = \frac{\partial}{\partial z}\left(\mu_{\perp}\frac{\partial p_{eV}}{\partial z} - \frac{j_{iz}}{e}\right)
\label{eq:potential}
\end{equation}

Here, we have defined the cross-field electron mobility $\mu_\perp = \frac{e}{m_e \nu_e}\frac{1}{1 + \Omega_e^2}$ and the electron pressure in $p_{eV} = n_e T_{eV}$ in $eV / m^{3}$. Lastly, we have an equation for the transport of electron internal energy:

\begin{equation}
\frac{\partial}{\partial t}\left(\frac{3}{2} n_e T_{eV}\right) + \frac{\partial}{\partial x}\left(\frac{5}{2} n_e T_{eV} u_{ez} + \kappa_e \frac{\partial T_{eV}}{\partial z}\right) = m_e n_e \nu_e |u_e|^2 + \nabla p_{eV} u_{ez} - S_{inelastic} - S_{wall}
\label{eq:electronenergy}
\end{equation}

In the above, $\kappa_e$ is the electron thermal conductivity, $|u_e| = \sqrt{u_{ez}^2 + u_{e\theta}^2}$ is the total electron speed, $S_{inelastic}$ represents energy losses due to inelastic collisions (including ionization and excitation), and $S_{wall}$ represents power loss to the walls in the radial direction. HallThruster.jl has multiple models for the thermal conductivity and the wall loss terms, while the inelastic loss term is given by

$S_{inelastic} = n_e \sum_{r} \Delta \epsilon_r n_{reactant} k_{r}\left(\frac{3}{2} T_{eV}\right)$

where for an ionization or excitation reaction $r$, $\Delta \epsilon_r$ is the activation/threshold energy, $n_{reactant}$ is the number density of the reactant, and $k_r$ is the rate coefficient as a function of electron energy.  


# Numerics

The finite volume method is applied to discretize the ion and neutral equations, transforming them into a system of ordinary differential equations. These are integrated in time using `DifferentialEquations.jl` [@rackauckas2017differentialequations]. This gives the end-user significant flexibility in their choice of time integration method, as all explicit integrators that work with `DifferentialEquations.jl` will work in `HallThruster.jl`. Several choice of numerical flux are available, including upwind, Rusanov, a global Lax-Friedrichs flux-splitting scheme, and HLLE.  The elliptic equation for the potential is transformed by finite differences into a tridiagonal linear system $A x = b$, which is solved using the Thomas' algorithm. This is done on a staggered grid to prevent odd-even decoupling. The electron energy equation is solved semi-implicitly to ease timestep restrictions, and is discretized in space using finite differences. The user may choose whether to use upwind (first-order) or central (second-order) differences depending on the application. Reactions and collisions are modelled as a function of the electron energy, and rate coefficients have computed using cross sections tables available in literature.

We use the method of manufactured solutions to verify that the PDEs are discretized correctly and obtain the design order of accuracy. We are aided in this by the `Symbolics.jl`[@gowda2021high] , which makes computation of the needed source terms simple. 

# Functionality

HallThruster.jl provides extensive options to allow the user to customize their simulation. These include

- Custom thruster geometry
- Magnetically shielded thrusters
- Order of accuracy control
- Anode and wall sheath model options
- Multiple propellants (krypton and xenon supported out of the box, but adding a new propellant is easy)
- Custom initialization
- Restarts
- Custom anomalous transport models
- User-provided extra source terms for all equations
- Custom additional collisions/reactions

Detailed documentation about all of these features and more is available at the Github repository (https://github.com/UM-PEPL/HallThruster.jl). More features will likely be added as the code continues to develop, and the authors welcome contributions and ideas from users.

# Example simulations

We show comparisons of our results to case 3 of the [LANDMARK benchmark suite](https://www.landmark-plasma.com/test-case-3), which itself has 3 cases. We show all results with 1024 cells, time averaged after 2 ms of simulation time. We compare three included fluxes: HLLE, local Lax-Friedrichs / Rusanov (LLF), and global Lax-Friedrichs (GLF) to the three sub-cases of the LANDMARK benchmark - two fluid models with adjustable viscosity parameter $\delta$ and a hybid-particle-in-cell case.

## LANDMARK case 1:

![LANDMARK case 1, 1024 cells.\label{fig:landmark_1}](landmark_1.png)

## LANDMARK case 2:

![LANDMARK case 2, 1024 cells\label{fig:landmark_2}](landmark_2.png)

## LANDMARK case 3:

![LANDMARK case 3, 1024 cells\label{fig:landmark_3}](landmark_3.png)

# Acknowledgements



# References