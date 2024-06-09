;; Documentation:
;; =====
;; This part declares how the two particle interact with each
;; other in the particle system. 

(in-package :phy)

;; ========== Cut-off ==========
;; If two particle distance is larger than cut-off length,
;; their potential and force will be thresholded to zero.

(defconstant +lj-cutoff+ 2.5
  "Lennard Jones Potential Cutoff. ")

(defconstant +atomic-cutoff+ 2.0
  "Atomic Potential Cutoff. ")

(defconstant +colloidal-cutoff+ 1.2
  "Colloidal Potential Cutoff. ")

;; ========== Prefactor Alpha ==========
;; The alpha is used to unify the potential, which is:
;; 
;;   alpha = 27 * rc^2 / (4 * (rc^2 - 1)^3)
;;

(flet ((cubic (x)  (* x x x)))
  (flet ((alpha (rc) (/ (* 27 rc rc) (* 4 (cubic (- (* rc rc) 1.0))))))
    (defconstant +atomic-alpha+ (alpha +atomic-cutoff+)
      "Atomic prefactor alpha. ")
    (defconstant +colloidal-alpha+ (alpha +colloidal-cutoff+)
      "Colloidal prefactor alpha. ")))

;; ========== Potential ==========
;; Potential between two particle.
;; + epsilon is assumed to be 1.0 (potential minimum = - epsilon)
;; + sigma   is assumed to be 1.0 (particle radius)
;; Note: this should be checked carefully to avoid formula error.

(declaim (inline %lj-potential %atomic-potential %colloidal-potential))
(defun %lj-potential (r)
  "Lennard Jones Potential. "
  (if (> r +lj-cutoff+) 0.0
      (let ((r6 (expt (/ 1.0 r) 6.0)))
        ;; u(r) = 4 * epsilon * ((s/r)^12 - (s/r)^6)
        (* 4.0 (- (* r6 r6) r6)))))

(defun %atomic-potential (r)
  "Atomic Potential. "
  (if (> r +atomic-cutoff+) 0.0
      (let ((s/r  (/ 1.0 r))              ; sigma / r = 1 / r
            (rc/r (/ +atomic-cutoff+ r))) ; r_c / r   = 2 / r
        (flet ((square (x) (* x x)))
          ;; phi(r) = epsilon * alpha * ((s/r)^2 - 1) * ((rc/r)^2 - 1)^2
          (* +atomic-alpha+
             (- (square s/r) 1.0) (square (- (square rc/r) 1.0)))))))

(defun %colloidal-potential (r)
  "Colloidal Potential. "
  (if (> r +colloidal-cutoff+) 0.0
      (let ((s/r  (/ 1.0 r))                 ; sigma / r = 1 / r
            (rc/r (/ +colloidal-cutoff+ r))) ; r_c / r   = 2 / r
        (flet ((square (x) (* x x)))
          ;; phi(r) = epsilon * alpha * ((s/r)^2 - 1) * ((rc/r)^2 - 1)^2
          (* +colloidal-alpha+
             (- (square s/r) 1.0) (square (- (square rc/r) 1.0)))))))

;; ========== Force ==========
;; Force between two particle.
;; The `vr' should be vector form source to particle.
;; The `f/r' should be -D[potential, r] / r, therefore,
;; f/r * vr will be in length of f/r, direction in vr.

(declaim (inline %lj-force %atomic-force %colloidal-force))
(defun %lj-force (vr)
  "Lennard Jones Force. "
  (let* ((r   (norm vr))
         (r6  (expt (/ 1.0 r) 6))       ; (sigma / r)^6
         (f/r (if (>= r +lj-cutoff+) 0.0
                  ;; 4 * epsilon * (12 * (s/r)^12 - 6 * (s / r)^6) / r
                  (* 4.0 (/ (- (* 12.0 r6 r6) (* 6.0 r6)) (* r r))))))
    (num-times-vec f/r vr)))

(flet ((square (x) (* x x)))
  (macrolet ((phi-force (vr cutoff)
               `(let* ((r      (norm ,vr))
                       ;; rc2-1 = ((rc/r)^2 - 1)
                       (rc2-1  (- (square (/ ,cutoff r)) 1.0))
                       ;; s/r2-1 = ((s/r)^2 - 1)
                       (s/r2-1 (- (/ 1.0 r) 1.0))
                       ;; f = (2 * s^2 * rc2-1^2 + 4 * rc^2 * rc2-1 * s/r2-1) / r^3
                       (f/r    (/ (+ (* 2.0 (square rc2-1))
                                     (* 4.0 (square ,cutoff) rc2-1 s/r2-1))
                                  (square (square r)))))
                  (num-times-vec f/r ,vr))))
    (defun %atomic-force (vr)
      "Atomic Force. "
      (phi-force vr +atomic-cutoff+))
    (defun %colloidal-force (vr)
      "Colloidal Force. "
      (phi-force vr +colloidal-cutoff+))))

;; ========== Pairwise-Mixin ==========
;; The subclass of `pairwise-mixin' should set `cutoff' independently.
;; The `pairwise-mixin' should not be used directly.

(defclass pairwise-mixin ()
  ((cutoff :reader system-cutoff))
  (:documentation
   "The `pairwise-mixin' defines the particle interaction.

The pairwise interaction of the `particle-system' is:
+ `pairwise-potential'
+ `pairwise-force'

Their corresponding higher-level function is:
+ `potential' for `particle-system' total potential
+ `particle-force' for `i' th particle force in system
"))

;; ========== General Protocol ==========
;; Each subclass of `pairwise-mixin' should implement these
;; two pairwise method. 

(defgeneric pairwise-potential (system i j)
  (:documentation "Return particle `i' and `j' potential. "))

(defgeneric pairwise-force (system i j)
  (:documentation "Return paritcle `j' force on particle `i'. "))

(defgeneric pairwise-force* (system displacement)
  (:documentation "Return pairwise force via `displacement'. "))

;; ========== Higher-Level Function ==========
;; The `potential' and `particle-force' are higher level
;; wrapper for `pairwise-potential' and `pairwise-force'. 

(declaim (inline potential particle-force))
(defun potential (system)
  "Return the system total potential. "
  (sum-iter-i* ((i j) (j (system-size system)))
    (pairwise-potential system i j)))

(defun particle-force (system i)
  "Return the `i' th particle's force in `system'. "
  (sum-iter-i* ((j (system-size system)))
    :reject (lambda (j) (= i j))
    :sum-init (make-array (list (system-dimension system)) :initial-element 0.0)
    :sum-method #'vec-plus-vec
    (pairwise-force system i j)))

;; ========== Lennard Jones Potential ==========

(defclass lennard-jones-pairwise-mixin (pairwise-mixin)
  ((cutoff :initform +lj-cutoff+))
  (:documentation
   "The potential between two particles follows Lennard Jones Potential.

Lennard Jones Potential:

    u(r) = epsilon * ((sigma / r)^12 - (sigma / r)^6)
"))

(defmethod pairwise-potential ((system lennard-jones-pairwise-mixin) i j)
  (%lj-potential (norm (distance system i j))))

(defmethod pairwise-force ((system lennard-jones-pairwise-mixin) i j)
  (%lj-force (distance system i j)))

(defmethod pairwise-force* ((system lennard-jones-pairwise-mixin) displacement)
  (%lj-force displacement))

;; ========== Atomic Potential ==========

(defclass atomic-pairwise-mixin (pairwise-mixin)
  ((cutoff :initform 2.0))
  (:documentation
   "The potential bewteen two particles follows Atomic Potential.

Atomic Potential:

    u(r) = epsilon * alpha * ((sigma / r)^2 - 1) * ((cutoff / r)^2 - 1)^2

where:
+ cutoff = `+atomic-cutoff+' 2.0
+ alpha  = `+atomic-alpha+'  1.0
"))

(defmethod pairwise-potential ((system atomic-pairwise-mixin) i j)
  (%atomic-potential (norm (distance system i j))))

(defmethod pairwise-force ((system atomic-pairwise-mixin) i j)
  (%atomic-force (distance system i j)))

(defmethod pairwise-force* ((system atomic-pairwise-mixin) displacement)
  (%atomic-force displacement))

;; ========== Colloidal Potential ==========

(defclass colloidal-pairwise-mixin (pairwise-mixin)
  ((cutoff :initform 1.2))
  (:documentation
   "The potential bewteen two particles follows Atomic Potential.

Atomic Potential:

    u(r) = epsilon * alpha * ((sigma / r)^2 - 1) * ((cutoff / r)^2 - 1)^2

where:
+ cutoff = `+atomic-cutoff+' 1.2
+ alpha  = `+atomic-alpha+'  114.106
"))

(defmethod pairwise-potential ((system colloidal-pairwise-mixin) i j)
  (%colloidal-potential (norm (distance system i j))))

(defmethod pairwise-force ((system colloidal-pairwise-mixin) i j)
  (%colloidal-force (distance system i j)))

(defmethod pairwise-force* ((system colloidal-pairwise-mixin) displacement)
  (%colloidal-force displacement))
