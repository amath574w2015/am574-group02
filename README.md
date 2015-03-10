# am574-group02
[Paper](paper/pAdaptiveDG.pdf)

Title: A p-Adaptive Discontinuous Galerkin Method for Hyperbolic Conservation Laws in 1D
Abstract:

Discontinuous Galerkin (DG) methods are becoming increasingly popular tools for the numerical integration of hyperbolic conservation laws. DG methods provide a natural extension of finite volume methods to higher orders while maintaining a compact stencil and exhibit a number of desirable computational features. However, as problems of larger and larger scopes are considered it is increasingly necessary to implement an adaptive method which devotes the finite degrees of freedom to where they are most needed in the domain. To that end we propose a novel p-adaptive DG scheme which uses a hierarchical basis and allows the degree of the local polynomial approximation to vary between cells. This method either adds or removes degrees of freedom for the next step in the integration based on the behavior of the coefficient of the highest degree polynomial basis present in the current approximation. The efficiency and accuracy performance of the proposed method will be measured against a non-adapting scheme on several standard tests.
