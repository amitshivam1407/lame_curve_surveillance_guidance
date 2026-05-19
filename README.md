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

## Representative Results

<table>
  <tr>
    <td align="center" width="50%">
      <h3>Lamé Curve Surveillance Concept</h3>
      <img src="figures/superellipse_path_following.gif" width="430"><br>
      <p>Smooth surveillance orbit around a rectangular boundary while respecting UAV turning constraints.</p>
    </td>
    <td align="center" width="50%">
      <h3>Lamé vs Elliptic Circumscription</h3>
      <img src="figures/ellipse_superellipse_path_following.gif" width="430"><br>
      <p>Comparison between Lamé-curve and elliptic circumscription under identical curvature constraints.</p>
    </td>
  </tr>
</table>

The guidance law combines:

- path tangential motion
- arcsine shaping function
- asymptotic convergence behavior

---

# Experimental Validation

Indoor flight experiments using the Crazyflie 2.0 validate the proposed guidance strategy under realistic conditions.

---

# Repository Structure

```text
matlab/     -> MATLAB simulation source code
figures/    -> figures and GIF demonstrations
videos/     -> simulation and experiment videos
paper/      -> journal paper and citation