# Sparse-representation-for-damage-identification-of-structural-systems
Abstract

Identifying damage of structural systems is typically characterized as an inverse problem which might be ill-conditioned
due to aleatory and epistemic uncertainties induced by measurement noise and modeling error. Sparse representation
can be used to perform inverse analysis for the case of sparse damage. In this paper, we propose a novel two-stage
sensitivity analysis-based framework for both model updating and sparse damage identification. Specifically, an ell-2
Bayesian learning method is firstly developed for updating the intact model and uncertainty quantification so as to set
forward a baseline for damage detection. A sparse representation pipeline built on a quasi ell-0 method, e.g., Sequential
Threshold Least Squares (STLS) regression, is then presented for damage localization and quantification. Additionally,
Bayesian optimization together with cross validation is developed to heuristically learn hyperparameters from data,
which saves the computational cost of hyperparameter tuning and produces more reliable identification result. The
proposed framework is verified by three examples, including a 10-story shear-type building, a complex truss structure,
and a shake table test of an eight-story steel frame. Results show that the proposed approach is capable of both
localizing and quantifying structural damage with high accuracy.

Codes are also readily generalizable to the work of "Sparse Bayesian learning for structural damage identification".
## Citation
<pre>
@article{chen2021sparse,
  title={Sparse representation for damage identification of structural systems},
  author={Chen, Zhao and Sun, Hao},
  journal={Structural Health Monitoring},
  volume={20},
  number={4},
  pages={1644--1656},
  year={2021},
  publisher={SAGE Publications Sage UK: London, England}
}
</pre>

<pre>
@article{chen2020sparse,
  title={Sparse Bayesian learning for structural damage identification},
  author={Chen, Zhao and Zhang, Ruiyang and Zheng, Jingwei and Sun, Hao},
  journal={Mechanical Systems and Signal Processing},
  volume={140},
  pages={106689},
  year={2020},
  publisher={Elsevier}
}
</pre>
