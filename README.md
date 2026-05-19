# Surveillance Guidance Using Lamé Curve Paths

This repository contains MATLAB implementations, simulation studies, and experimental validation results associated with the paper:

> Amit Shivam and Ashwini Ratnoo  
> **"Surveillance Guidance Using Lamé Curve Paths"**  
> *Journal of Guidance, Control, and Dynamics (JGCD), 2025*

---

# Overview

This work introduces a curvature-constrained surveillance guidance framework for UAV monitoring of rectangular boundaries using smooth Lamé curve paths.

The proposed method replaces conventional elliptic circumscription with a family of continuous-curvature Lamé curves, providing:

- Reduced surveillance path length
- Improved rectangular boundary approximation
- Feasible curvature-constrained motion
- Smooth vector-field guidance
- Experimental quadrotor validation

---

# Key Contributions

- Continuous-curvature Lamé curve surveillance paths
- Analytical curvature characterization
- Minimal path-length circumscription
- Direct vector-field guidance law
- Lyapunov-based convergence proof
- Crazyflie quadrotor experimental validation

---

# Lamé Curve Surveillance Concept

<p align="center">
<img src="figures/rectangle_circumscription.png" width="650">
</p>

The proposed approach generates smooth surveillance orbits around rectangular regions while respecting UAV turning constraints.

---

# Lamé vs Elliptic Circumscription

<p align="center">
<img src="figures/lame_vs_ellipse.png" width="700">
</p>

Lamé curves achieve significantly lower path length compared to classical elliptic circumscription while satisfying identical curvature constraints.

---

# Vector Field Guidance

<p align="center">
<img src="figures/vector_field_overview.png" width="650">
</p>

The guidance law combines:

- path tangential motion
- arcsine shaping function
- asymptotic convergence behavior

---

# Experimental Validation

<p align="center">
<img src="figures/experimental_results.png" width="750">
</p>

Indoor flight experiments using the Crazyflie 2.0 validate the proposed guidance strategy under realistic conditions.

---

# Repository Structure

```text
matlab/     -> MATLAB simulation source code
figures/    -> figures and GIF demonstrations
videos/     -> simulation and experiment videos
paper/      -> journal paper and citation