(in-package :phy)

;; ========== Modular Dynamics ==========
;; The `modular-dynamics-mixin' is just a dummy `simulation-mixin'.
;; It wouldn't do anything for simulation, but it is a good
;; wrapper if you want to implement your own MD algorithm.
;;
;; For detailed example, see `verlet.lisp' for `md-verlet-mixin'.

(defclass modular-dynamics-mixin (simulation-mixin)
  ((delta-t :initform 0.01 :initarg :dt)
   velocity-config)
  (:documentation
   "The Modular Dynamics Simulation. "))

;; ========== General Protocol ==========

(defgeneric paricle-velocity (system i)
  (:documentation "Return the `i' th particle system velocity. "))

(defgeneric particle-kinetic (system i)
  (:documentation "Return the `i' th particle system kinetic energy. "))

(defmethod particle-velocity ((system modular-dynamics-mixin) i)
  (with-slots (velocity-config) system
    (aref velocity-config i)))

(defmethod particle-kinetic ((system modular-dynamics-mixin) i)
  (flet ((square (x) (* x x)))
    (* 0.5 (square (norm (particle-velocity system i))))))

(defun kinetic (system)
  "Calculate the kinetic of the `system'."
  (piter-i* ((i (system-size system)))
    (particle-kinetic system i)))

(defmethod (setf particle-velocity) (velocity (system modular-dynamics-mixin) i)
  (with-slots (velocity-config) system
    (setf (aref velocity-config i) velocity)))

(defmethod initialize-instance :after
    ((system modular-dynamics-mixin) &key)
  (with-slots (velocity-config) system
    (let ((dimension   (system-dimension   system))
          (temperature (system-temperature system)))
      (setf velocity-config (make-array (list (system-size system))))
      (dotimes (i (system-size system))
        (setf (aref velocity-config i)
              (boltzmann-random-vector dimension temperature))))))

;; The `simulation-collect' is (step potential kinetic)
(defmethod collect-simulation ((system modular-dynamics-mixin))
  (with-slots ((step    simulation-step-counter)
               (counter simulation-collect-counter)
               (collect simulation-collect))
      system
    (setf (aref collect counter) (list step (potential system) (kinetic system)))
    (incf counter)))

(defmethod plot-simulation ((system modular-dynamics-mixin) out-path)
  (with-slots ((collect simulation-collect)
               (counter simulation-collect-counter)
               (step    simulation-step-counter))
      system
    (multiple-value-bind (potentials kinetics)
        (loop for i below counter
              for (step potential kinetic) = (aref collect i)
              collect (list step potential) into potentials
              collect (list step kinetic)   into kinetics
              finally (return (values potentials kinetics)))
      (let* ((y-min (min (reduce #'min potentials :key #'second)
                         (reduce #'min kinetics   :key #'second)))
             (y-max (max (reduce #'max potentials :key #'second)
                         (reduce #'max kinetics   :key #'second))))
        (with-present-to-file
            (plot plot :margin 10
                       :x-min 0 :x-max step
                       :y-min y-min :y-max y-max
                       :x-label "s" :y-label "V/T")
            (out-path :width 800 :height 400)
          (add-plot-data plot (line-plot-pane potential :color +大红+) potentials)
          (add-plot-data plot (line-plot-pane potential :color +月白+) kinetics)
          (add-plot-legend (plot :position :top-left :padding 0.01)
            ((format nil "V, T = ~,3f, N = ~d, ⍴ = ~,3f"
                     (system-temperature system)
                     (system-size system)
                     (system-density system))
             :color +大红+)
            ((format nil "T, T = ~,3f, N = ~d, ⍴ = ~,3f"
                     (system-temperature system)
                     (system-size system)
                     (system-density system))
             :color +月白+)))))))
