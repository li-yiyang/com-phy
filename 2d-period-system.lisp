(in-package :phy)

(defclass lj-particle-system (particle-system lennard-jones-pairwise-mixin geo-mixin)
  ())

(defclass 2d-period-system (particle-system period-geo-mixin) ())

(defclass 2d-lj-period (2d-period-system lennard-jones-pairwise-mixin) ())

(defclass 2d-atomic-period (2d-period-system atomic-pairwise-mixin) ())

(defclass 2d-lj-period-mc (monte-carlo-mixin 2d-lj-period) ())

(defclass 2d-lj-period-mc-acc (mc-accelerate-mixin 2d-lj-period) ())

(defclass 2d-atomic-period-mc (monte-carlo-mixin 2d-atomic-period) ())

(defclass 2d-atomic-period-mc-acc (mc-accelerate-mixin 2d-atomic-period) ())

(defclass 2d-lj-period-md (md-verlet-mixin 2d-lj-period) ())
