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
  (sum-piter-i* ((i (system-size system)))
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

(defmethod plot-system ((system modular-dynamics-mixin) out-path
                        &key (particle-size 2) (color +草白+)
                          (particle-style :circle)
                        &allow-other-keys)
  (with-slots (config size) system
    (with-present-to-file
        (plot plot :margin 10
                   :x-min 0 :x-max (first  (system-lengths system))
                   :y-min 0 :y-max (second (system-lengths system)))
        (out-path)
      (add-plot-data plot
        (scatter-pane particle :color color :point-size particle-size
                               :point-style particle-style)
        (flet ((xy (p-config) (collect-i* ((i 2)) (aref p-config i))))
          (map 'list #'xy config)))
      (add-plot-legend (plot :position :top-left :padding 0.01)
        ((format nil "V = ~,3f" (potential system)) :color color)
        ((format nil "T = ~,3f" (kinetic system)) :color color))))
  out-path)

(defmethod plot-simulation ((system modular-dynamics-mixin) out-path
                            &key (x-min 0 x-min-set?) (x-max 0 x-max-set?)
                              (y-min 0 y-min-set?) (y-max 0 y-max-set?))
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
      (let* ((y-min (if y-min-set? y-min
                        (min (reduce #'min potentials :key #'second)
                             (reduce #'min kinetics   :key #'second))))
             (y-max (if y-max-set? y-max
                        (max (reduce #'max potentials :key #'second)
                             (reduce #'max kinetics   :key #'second)))))
        (with-present-to-file
            (plot plot :margin 10
                       :x-min (if x-min-set? x-min 0)
                       :x-max (if x-max-set? x-max step)
                       :y-min y-min
                       :y-max y-max
                       :x-label "s" :y-label "V/T")
            (out-path :width 800 :height 400)
          (add-plot-data plot (line-plot-pane potential :color +大红+)
            potentials)
          (add-plot-data plot (line-plot-pane average-V :color +水红+)
            (average-list potentials))
          (add-plot-data plot (line-plot-pane kinetic   :color +月白+)
            kinetics)
          (add-plot-data plot (line-plot-pane average-T :color +豆绿+)
            (average-list kinetics))
          (add-plot-legend (plot :position :top-right :padding 0.01)
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
