# Path Integral Methods for Polarons in Real Materials

**Imperial College London, Department of Physics, EXSS Group**

**Author:** Bradley A. A. Martin

**Supervisor:** Dr. Jarvist Moore Frost

**Co-supervisor:** Prof. Jenny Nelson

## 1. Research Summary

The central aim of this PhD project has been to develop methods that can offer fully predictive temperature-dependent mobilities for a wide variety of systems of technical interest and thereby offer design clues for and methods to screen new materials computationally. Though this project is primarily theoretical and computational, it has required deep familiarisation and contact with experimental methods, material-specific properties, and interaction with experimental collaborators and peers who have been measuring mobilities and response functions in these materials.

This PhD project has further developed Feynman’s variational path integral method [1] to predict materials' ‘polaronic’ properties. The polaron is a quasiparticle formed from the interaction of charge carriers with the vibrational modes of a material's atomic lattice. This can equivalently be thought of as a quantum field theory of a single charged fermion interacting with a bosonic field; this is the phonon field for the polaron. The formation of a polaron can greatly alter the properties of a material [2], and it is important to understand polarons to predict the mobility of charge carriers accurately. The extent of a polaron can be large (encompassing multiple lattice sites), where the polaron often moves around the material, or small (comparable or smaller than a lattice constant), where the polaron becomes localised to individual lattice sites. The large polaron is commonly modelled by the Fröhlich Hamiltonian [3], and the small polaron by the Holstein Hamiltonian [4, 5].

The key methodology this project builds upon is Feynman’s variational method [1]. This method translates the Fröhlich polaron Hamiltonian into an action, which is used to derive the partition function for the system as a path integral. Since the phonon coordinates enter the action quadratically, the phonon path integral is Gaussian. Its closed-form expression results in an effective temporally non-local potential acting on the electron. This partition function is then approximated by a trial partition function derived from a trial action, requiring that this partition function has a closed-form expression (typically limiting the trial action to be quadratic). Feynman chose this trial action to represent an electron harmonically coupled to a fictitious particle via a spring. Using Jensen’s inequality for convex functions, Feynman derived an inequality for the polaron-free energy. An optimal upper bound is then found by varying the mass and spring constant of the trial model. Further developments were also made using the method of ‘influence functionals’ to derive response functions, such as the mobility, from a functional Taylor expansion of the polaron action about the optimal trial action [6] - this method was then further used to obtain the optical conductivity [7].

## 2. Research Milestones

Overall, this project may be split into four main milestones (not necessarily presented in chronological order):

1. **Improve Trial Model**: Extend the model polaron Lagrangian (while retaining the analytic solution) to increase the accuracy of the approximation.
   
2. **Generalise Material Model**: Extend the effective Lagrangian to provide greater material-specific detail. Diagrammatic Monte Carlo or Path Integral Monte Carlo algorithms can provide a direct evaluation of the effective Lagrangian, and a code for this could be written to compare to the Feynman variational results.
   
3. **Compute Response Functions**: Compute further response functions of the polaron variation state, which enable comparison to experimental observables. For instance, the optical absorption of the polaron state could be compared to transient absorption measurements on these materials; the frequency-dependent mobility could be calculated and compared to Terahertz and microwave conductivity measurements.
   
4. **High-Throughput Material Classification**: Apply these codes to material groups and classes. This requires characterising materials, either by recourse to the material databases of synthetic data derived from electronic structure calculations, such as the Materials Project, or through the use of standard electronic structure packages, such as VASP and Gauss.

So far in the literature, the variational method has been generalised to include: finite temperatures [8], effective phonon frequencies [9], Bose-polarons [10], bipolarons [11], many-polarons [12], magnetic fields [13], multiple phonon modes [14], anisotropic band-mass [15], anharmonic phonons [16], external linear-response electric fields [6], non-linear response electric fields [17], multiple fictitious particles in the trial model (by myself), the optimal self-consistent quadratic trial model [18, 19], optimal self-consistent response functions [20]. This list is not exhaustive. My work includes the extension to multiple phonon modes, multiple fictitious particle trial models, and anisotropic band masses. I recently extended the variational method to the Holstein model [21]. Many of the extensions above to the method need more available code. Therefore, part of my PhD project has also involved developing an open-source package written in Julia [22] that implements some theoretical developments [23].

## 3. Research Progress

So far, I have accomplished the following during my PhD:

- I co-developed a Julia code package PolaronMobility.jl [23] that efficiently evaluates the original variational principles, developed in [1, 6, 7, 9], using numerical optimisation and Gauss-Kronrod integration.
I wrote highly efficient Julia code [23] that allows for evaluating the complex mobility at all temperatures and applied frequencies, pushing the calculation's limits further than those given in the literature. The numerical results agree very well with numerically exact Diagrammatic Monte Carlo results [14, 24].
- I derived asymptotic power series expansions of the complex mobility in terms of hypergeometric functions [14].
I extended the Feynman Variational Method to a Fröhlich polaron model with multiple phonon modes and implemented this in Julia code [23].
- I applied the multiple phonon mode extended variational method to Methylammonium Lead Iodide perovskite material and predicted its polaronic properties and complex conductivity [14]. The results agree well with experimental Terahertz Spectroscopy measurements [25].
- I extended the Feynman Variational Method to a Fröhlich polaron model with anisotropic band-masses [15].
I extended the trial model to include multiple fictitious particles, improving the variational solution and writing corresponding code [23].
I co-developed a Path Integral Monte Carlo (PIMC) Julia code package PolaronQMC.jl [26] that calculates the temperature-dependent polaron energy for the Fröhlich model. I then co-supervised four MSci and 2 UROP students, who have since developed this code further to include multiple phonon modes and the Holstein model. They also wrote a complementary Diagrammatic Monte Carlo (DiagMC) code package, Tethys.jl [27], for the ground-state energy.
- I extended the Feynman Variational Method to general phonon-momentum dependence and electron-phonon coupling matrices and wrote corresponding code [23].
- I developed the theory for a coherent-state path integral [28] version of the Feynman Variational Method and hypothesised an analogy to the Luttinger-Ward functional and cumulant Green functions. This is an initial step towards generalising many-polaron systems within grand canonical ensembles.
- I extended the Feynman Variational Method to the Holstein lattice-polaron model and applied it to the organic semiconductor Rubrene [21].
- Performed high-throughput application of the variational method to predict polaronic properties and mobilities of 1260 materials with PhD student Yande Fu.

This covers all the “Improve Trial Model”, “Generalise Material Model” and “Compute Response Functions” and “High-Throughput Material Classification” milestones.

## 4. Dissemination

### 4.1 Publications

- B. A. A. Martin and J. M. Frost, Predicting polaron mobility in organic semiconductors with the Feynman variational approach, 2023. Preprint arXiv:2207.06846 (2024).

- B. A. A. Martin and J. M. Frost, “Multiple phonon modes in feynman path-integral variational polaron mobility”, Phys. Rev. B 107, 115203 (2023).

- X. Zheng, T. R. Hopper, A. Gorodetsky, M. Maimaris, W. Xu, B. A. A. Martin, J. M. Frost, and A. A. Bakulin, “Multipulse terahertz spectroscopy unveils hot polaron photoconductivity dynamics in metal-halide perovskites”, The Journal of Physical Chemistry Letters 12, PMID: 34478291, 8732–8739 (2021).

- B. Guster, P. Melo, B. A. A. Martin, V. Brousseau-Couture, J. C. de Abreu, A. Miglio, M. Giantomassi, M. Cote, J. M. Frost, M. J. Verstraete, and X. Gonze, “Frohlich polaron effective mass and localization length in cubic materials: degenerate and anisotropic electronic bands”, Phys. Rev. B 104, 235123 (2021).

### 4.2 Talks / Posters

- June 6 2023 - CECAM PI-School Invited Talk at Tel Aviv Uni: “Extending the FVA for Polarons in Real Materials.”

- May 25 2023 - TYC Invited Talk at UCL: “Extending the FVA for Polarons in Real Materials.”

- Dec 8 2022 - CPE Christmas Symposium Invited Talk at ICL: “Calculating polaron mobility with Path Integrals.”

- Nov 9 2022 - TYC Lunchtime Symposium Online Talk: “Extending the FVA for Polarons in Real Materials.”

- Aug 22-25 2022 - Psi-K Lausanne Conference Poster: “Extending the FVA for Polarons in Real Materials.”

- June 16 2021 - CECAM Workshop Poster: “Variational path integrals give accurate and efficient polaron mobilities.”

## 5 References

[1] R. P. Feynman, “Slow electrons in a polar crystal”, Physical Review 97, 660–665 (1955).

[2] C. Franchini, M. Reticcioli, M. Setvin, and U. Diebold, “Polarons in materials”, Nature Reviews Materials 6, 560–586 (2021).

[3] H. Frohlich, “Electrons in lattice fields”, Advances in Physics 3, 325–361 (1954).

[4] T. Holstein, “Studies of polaron motion: part i. the molecular-crystal model”, Annals of Physics 8, 325–342 (1959).

[5] T. Holstein, “Studies of polaron motion: part ii. the “small” polaron”, Annals of Physics 8, 343–389 (1959).

[6] R. P. Feynman, R. W. Hellwarth, C. K. Iddings, and P. M. Platzman, “Mobility of slow electrons in a polar crystal”, Physical Review 127, 1004–1017 (1962).

[7] J. Devreese, J. D. Sitter, and M. Goovaerts, “Optical absorption of polarons in the feynman-hellwarth-iddings-platzman approximation”, Physical Review B 5, 2367–2381 (1972).

[8] Y. Osaka, “Polaron state at a finite temperature”, Progress of Theoretical Physics 22, 437–446 (1959).

[9] R. W. Hellwarth and I. Biaggio, “Mobility of an electron in a multimode polar lattice”, Physical Review B 60, 299–307 (1999).

[10] T. Ichmoukhamedov and J. Tempere, “General memory kernels and further corrections to the variational path integral approach for the bogoliubov-frohlich hamiltonian”, Phys. Rev. B 105, 104304 (2022).

[11] G. Verbist, F. M. Peeters, and J. T. Devreese, “Large bipolarons in two and three dimensions”, Phys. Rev. B 43, 2712–2720 (1991).

[12] F. Brosens, S. N. Klimin, and J. T. Devreese, “Variational path-integral treatment of a translation invariant many-polaron system”, Phys. Rev. B 71, 214301 (2005).

[13] F. M. Peeters and J. T. Devreese, “Statistical properties of polarons in a magnetic field. i. analytic results”, Phys. Rev. B 25, 7281–7301 (1982).

[14] B. A. A. Martin and J. M. Frost, “Multiple phonon modes in feynman path-integral variational polaron mobility”, Phys. Rev. B 107, 115203 (2023).

[15] B. Guster, P. Melo, B. A. A. Martin, V. Brousseau-Couture, J. C. de Abreu, A. Miglio, M. Giantomassi, M. Cot ˆ ´e, J. M. Frost, M. J. Verstraete, and X. Gonze, “Frohlich polaron effective mass and localization length in cubic materials: degenerate and anisotropic electronic bands”, Phys. Rev. B 104, 235123 (2021).

[16] M. Houtput and J. Tempere, “Beyond the frohlich hamiltonian: path-integral treatment of large polarons in anharmonic solids”, Phys. Rev. B 103, 184306 (2021).

[17] K. K. Thornber and R. P. Feynman, “Velocity acquired by an electron in a finite electric field in a polar crystal”, Phys. Rev. B 1, 4099–4114 (1970).

[18] J. Adamowski, B. Gerlach, and H. Leschke, “Feynman’s approach to the polaron problem generalized to arbitrary quadratic actions”, in Functional integration: theory and applications, edited by J.-P. Antoine and E. Tirapegui (Springer US, Boston,
MA, 1980), pp. 291–301.

[19] J. Adamowski, B. Gerlach, and H. Leschke, “General aspects of the functional integral approach to the polaron and related systems”, in Polarons and excitons in polar semiconductors and ionic crystals, edited by J. T. Devreese and F. Peeters (Springer US, Boston, MA, 1984), pp. 245–269.

[20] K. K. Thornber, “Linear and nonlinear electronic transport in electron-phonon systems: self-consistent approach within the path-integral formalism”, Phys. Rev. B 3, 1929–1941 (1971).

[21] B. A. A. Martin and J. M. Frost, Predicting polaron mobility in organic semiconductors with the feynman variational approach, 2022.

[22] J. Bezanson, A. Edelman, S. Karpinski, and V. B. Shah, “Julia: a fresh approach to numerical computing”, SIAM review 59, 65–98 (2017).

[23] B. A. A. Martin and J. M. Frost, Polaronmobility.jl, https://github.com/jarvist/PolaronMobility.jl, 2021.

[24] A. S. Mishchenko, L. Pollet, N. V. Prokof’ev, A. Kumar, D. L. Maslov, and N. Nagaosa, “Polaron mobility in the “beyond quasiparticles” regime”, Physical Review Letters 123, 10.1103/physrevlett.123.076601 (2019).

[25] X. Zheng, T. R. Hopper, A. Gorodetsky, M. Maimaris, W. Xu, B. A. A. Martin, J. M. Frost, and A. A. Bakulin, “Multipulse terahertz spectroscopy unveils hot polaron photoconductivity dynamics in metal-halide perovskites”, The Journal of Physical Chemistry Letters 12, PMID: 34478291, 8732–8739 (2021).

[26] Y. C. Wong, L. C. Filipovich, G. Su, B. A. A. Martin, and J. M. Frost, Polaronqmc.jl, https://github.com/Frost-group/PolaronQMC.jl, 2022.

[27] D. Lee, X. Yang, and J. M. Frost, Tethys.jl, https://github.com/Frost- group/Tethys.jl, 2022.

[28] A. Altland and B. D. Simons, Condensed matter field theory, 2nd ed. (Cambridge University Press, 2010).
