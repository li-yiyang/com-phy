(in-package :phy)

(defclass 2d-period-system (particle-system period-geo-mixin)
  ((dimension :initform 2))
  (:documentation "2D period system. "))

(defclass 2d-lj-period        (2d-period-system lennard-jones-pairwise-mixin) ())
(defclass 2d-atomic-period    (2d-period-system atomic-pairwise-mixin)        ())
(defclass 2d-colloidal-period (2d-period-system colloidal-pairwise-mixin)     ())

(defclass 2d-lj-period-mc     (2d-lj-period monte-carlo-mixin)   ())
(defclass 2d-lj-period-mc-acc (2d-lj-period mc-accelerate-mixin) ())

(defclass 2d-atomic-period-mc     (2d-atomic-period monte-carlo-mixin)   ())
(defclass 2d-atomic-period-mc-acc (2d-atomic-period mc-accelerate-mixin) ())

(defclass 2d-colloial-period-mc  (2d-colloidal-period monte-carlo-mixin)   ())
(defclass 2d-colloial-period-acc (2d-colloidal-period mc-accelerate-mixin) ())

(defclass 2d-lj-period-md (md-velocity-verlet-mixin 2d-lj-period) ())
