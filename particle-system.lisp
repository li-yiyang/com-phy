(in-package :phy)

;; ========== Particle System ==========

(defclass particle-system ()
  ((size        :initarg :N           :initform 2   :reader system-size)
   (density     :initarg :density     :initform 1.0 :reader system-density)
   (temperature :initarg :temperature :initform 1.0 :reader system-temperature)
   (dimension                         :initform 2   :reader system-dimension)
   config)
  (:documentation
   "The `particle-system' is the general interface for simulation.

The `size' is the number of particles in the `particle-system';
The `density' will be used to calculate the geometry of `particle-system';
The `dimension' is fixed constant for each particle-system subclass;

The `config' should be taken care with `geo-mixin';
The `config' maybe refered directly in modular dynamics process.
"))

;; ========== General Protocol ==========

(defgeneric particle (system i)
  (:documentation "Return the particle `i' in `system'. "))

(defgeneric (setf particle) (particle system i)
  (:documentation "Set the `i' th particle in `system'. "))

;; There's a `constrain-displacement' method in `geo-mixin.lisp'
;; to constrain the displacment for specific geometry.
(defgeneric move-particle (system i displacement)
  (:documentation "Move `i' th particle in `system' with `displacement'. "))

;; The subclass of `particle-system' should implement `save-system'
;; with `:after' method. 
(defgeneric save-system (system out-path)
  (:documentation "Save `system' infomation to `out-path'. "))

(defgeneric %load-system (system path)
  (:documentation "Load `system' from `path'. "))

(defgeneric system-step-size (system)
  (:documentation "The step size for system (MD dt / MC step-size). "))

;; ========== Implementation ==========

(defmethod particle ((system particle-system) i)
  (aref (slot-value system 'config) i))

(defmethod (setf particle) (particle (system particle-system) i)
  (setf (aref (slot-value system 'config) i) particle))

;; Use `geo-mixin' subclass to set `:after' method of `move-particle'
;; to constrain the particle movement. 
(defmethod move-particle ((system particle-system) i displacement)
  (with-slots (config) system
    (setf (aref config i) (vec-plus-vec (aref config i) displacement))))

;; ========== Plot System ==========

(defmethod plot-system ((system particle-system) out-path
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
        ((format nil "V = ~,3f" (potential system)) :color color))))
  out-path)

;; ========== Save/Load System ==========

(defmethod save-system ((system particle-system) out-path)
  (unless (uiop:directory-pathname-p out-path)
    (errorf "The `out-path' ~a should be directory. " out-path))
  (uiop:with-current-directory (out-path)
    (with-open-file (stream "particle-system.lisp"
                            :direction :output
                            :if-exists :supersede
                            :if-does-not-exist :create)
      (print (list :class       (class-name (class-of system))
                   :size        (system-size        system)
                   :density     (system-density     system)
                   :temperature (system-temperature system)
                   :dimension   (system-dimension   system)
                   :step-size   (system-step-size   system))
             stream))
    (with-open-file (stream "config.lisp"
                            :direction :output
                            :if-exists :supersede
                            :if-does-not-exist :create)
      (print (slot-value system 'config) stream))
    out-path))

;; for `particle-system', `%load-system' should do nothing,
;; the `%load-system' should add `:after' methods
(defmethod %load-system ((system particle-system) path))

(defun load-system (path)
  "Load and return the particle-system from `path'.

The `path' should be a directory.
Normally, within the `path' there should be:
+ `particle-system.lisp' for the infomation of the basic
  `particle-system' infomation;
+ `config.lisp' for the config infomation;

Other subclass or functional mixin may dump their own infomations.
"
  (unless (uiop:directory-pathname-p path)
    (errorf "The `path' ~a should be directory. " path))
  (let (system)
    (uiop:with-current-directory (path)
      (with-open-file (stream "particle-system.lisp")
        (let ((info (read stream)))
          (setf system
                (make-instance (getf info :class)
                               :temperature (getf info :temperature)
                               :density     (getf info :density)
                               :N           (getf info :size)
                               :step-size   (getf info :step-size)
                               :scheme      :skip-init))))
      (setf (slot-value system 'config)
            (with-open-file (stream "config.lisp") (read stream))))
    (%load-system system path)
    system))
