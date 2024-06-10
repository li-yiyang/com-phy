(in-package :phy)

(defclass simulation-mixin ()
  ((simulation-step-counter    :initform 0 :initarg :simulation-start-step)
   (simulation-collect-counter :initform 0)
   (simulation-collect         :initform (make-array '(0) :adjustable t)
                               :reader simulation-collect))
  (:documentation
   "The general simulation-mixin functions. "))

;; ========== General Protocol ==========

;; The `simulation-step' should not be directly called, 
;; instead, use `run-simulation'. 
(defgeneric simulation-step (system)
  (:documentation
   "Perform a single step for `system' simulation.
Make sure the `init-simulation' is called before. "))

;; The `init-simulation' will 
(defgeneric init-simulation (system)
  (:documentation
   "Initialize the simulation for `system'.
This method must be called first before `simulation-step'. "))

(defgeneric collect-simulation (system)
  (:documentation
   "Collect info into `simulation-collect'. "))

(defgeneric lazy-potential (system)
  (:documentation
   "Quick potential for system used for `collect-simulation'. "))

;; Some Post Process Methods for Simulation

(defgeneric plot-simulation (system out-path &key y-min y-max x-min x-max
                             &allow-other-keys)
  (:documentation "Plot the simulation collect results. "))

;; ========== Save and Load ==========

(defmethod save-system :after ((system simulation-mixin) out-path)
  (with-slots (simulation-step-counter
               simulation-collect-counter
               simulation-collect)
      system
    (uiop:with-current-directory (out-path)
      (with-open-file (stream "simulation.lisp"
                              :direction :output
                              :if-exists :supersede
                              :if-does-not-exist :create)
        (print (list :step-counter    simulation-step-counter
                     :collect-counter simulation-collect-counter
                     :collect         simulation-collect)
               stream)))))

(defmethod %load-system :after ((system simulation-mixin) out-path)
  (with-slots (simulation-step-counter
               simulation-collect-counter
               simulation-collect)
      system
    (uiop:with-current-directory (out-path)
      (with-open-file (stream "simulation.lisp")
        (let* ((info    (read stream))
               (collect (getf info :collect))
               (length  (array-dimension collect 0)))
          (setf simulation-step-counter    (getf info :step-counter)
                simulation-collect-counter (getf info :collect-counter))
          (adjust-array simulation-collect (list length))
          (dotimes (i length)
            (setf (aref simulation-collect i) (aref collect i))))))))

;; ========== Simulation Interface ==========

(defun run-simulation (system steps
                       &key
                         (collect-every      1   collect-every-set?)
                         (plot-config-every  100 plot-config-set?)
                         (plot-collect-every 100 plot-collect-set?)
                         (out-path (uiop:getcwd))
                         (plot-format "config-~d.png")
                         (sim-format  "simulation-~d.png"))
  "Run the simulation for `system' for `steps'.

Collect the simulation infomation via `collect-simulation' method
every `collect-every' step if setted with `:collect-every' keyword;

Plot the config via `plot-system' method every `plot-config-every'
step if setted with `:plot-config-every' keyword, the system
config plot will be saved to `out-path' dir with name formatted
by `plot-format';

Plot the simulation collected infomation via `plot-simulation'
method every `plot-collect-every' step if setted with
`:plot-collect-every' keyword, like `plot-config-every';

The `plot-format' and `sim-format' should be like `xxx-~d.png',
the `~d' part will be current simulation step counter, aka,
the `simulation-step-counter'.

The system will be initialized via `init-simulation' method if
`simulation-step-counter' is zero (first run simulation). 
"
  (with-slots ((counter simulation-collect-counter)
               (collect simulation-collect)
               (step    simulation-step-counter))
      system
    ;; first run, init system for simulation
    (init-simulation system)
    ;; if collect, enlarge the `simulation-collect' for collecting
    (when (or collect-every-set? plot-config-set? plot-collect-set?)
      (adjust-array collect (list (+ counter (truncate steps collect-every)))))
    ;; iter `steps' times for `simulation-step'
    (uiop:with-current-directory (out-path)
      (loop for s below steps
            do (simulation-step system)

            if (and (or collect-every-set? plot-config-set? plot-collect-set?)
                    (zerop (mod step collect-every)))
              do (collect-simulation system)
                 
            if (and plot-config-set?  (zerop (mod step plot-config-every)))
              do (plot-system     system (format nil plot-format step))
                 
            if (and plot-collect-set? (zerop (mod step plot-collect-every)))
              do (plot-simulation system (format nil sim-format  step))
                 
            do (incf step)))))
