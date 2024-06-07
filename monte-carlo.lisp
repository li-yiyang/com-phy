(in-package :phy)

;; ========== Monte Carlo ==========

(defclass monte-carlo-mixin (simulation-mixin)
  ((step :initform 0.1 :initarg :step-size)
   old-potential)
  (:documentation
   "The Monte Carlo Simulation. "))

;; store the `old-potential' to at least accelerate the calculation
(defmethod init-simulation ((system monte-carlo-mixin))
  (setf (slot-value system 'old-potential) (potential system)))

(defmethod simulation-step ((system monte-carlo-mixin))
  (let ((temperature (system-temperature system))
        (dimension   (system-dimension   system))
        (size        (system-size system))
        (step        (slot-value system 'step)))
    (iter-i* ((i size))
      (let* ((old (copy-array (particle system i)))
             (old-potential   (slot-value system 'old-potential))
             (move (num-times-vec (random step) (random-vector dimension))))
        (move-particle system i move)
        (let* ((potential (potential system))
               (-delta-E   (- old-potential potential))
               (accept? (or (>= -delta-E 0.0)
                            (possible (exp (/ -delta-E temperature))))))
          (if accept?
              ;; if accepted, update the old-potential
              (setf (slot-value system 'old-potential) potential)
              ;; else, restore the particle position
              (setf (particle system i) old)))))))

;; The `simulation-collect' is (step potential)
(defmethod collect-simulation ((system monte-carlo-mixin))
  (with-slots ((step    simulation-step-counter)
               (counter simulation-collect-counter)
               (collect simulation-collect))
      system
    (setf (aref collect counter) (list step (slot-value system 'old-potential)))
    (incf counter)))

;; The `plot-simulation' will plot the potential along the MC simulation. 
(defmethod plot-simulation ((system monte-carlo-mixin) out-path)
  (with-slots ((counter simulation-collect-counter)
               (collect simulation-collect))
      system
    (let* ((data (collect-i* ((i counter)) (aref collect i)))
           (y-min (reduce #'min data :key #'second))
           (y-max (reduce #'max data :key #'second))
           (y-max (cond ((< (- y-max y-min) 1.0) (+ y-min 1.0))
                        ((and (> y-max 0.0) (< y-min 0.0)) 0.0)
                        (t y-max))))
      (with-present-to-file
          (plot plot :margin 10
                     :x-min 0 :x-max counter
                     :y-min y-min :y-max y-max
                     :x-label "s" :y-label "V")
          (out-path :width 800 :height 400)
        (add-plot-data plot (line-plot-pane potential :color +大红+) data)
        (add-plot-legend (plot :position :top-left :padding 0.01)
          ((format nil "V, T = ~,3f, N = ~d, ⍴ = ~,3f"
                   (system-temperature system)
                   (system-size system)
                   (system-density system))
           :color +大红+))))))

;; ========== MC Accelerate Mixin ==========
;; Compared with `monte-carlo-mixin' for 100 particles, 500 steps:
;;
;; | mc-acclerate-mixin | monte-carlo-mixin |
;; |--------------------+-------------------|
;; |             17.876 |           133.085 |

(defclass mc-accelerate-mixin (monte-carlo-mixin) ()
  (:documentation
   "The boost of Monte Carlo Simulation. 

Compared with `monte-carlo-mixin', the `mc-accelerate-mixin' will
not compute the total potential for delta potential, instead,
it will just calculate the delta potential as:

  Delta(sum(pairwise-potential(i, j), (j)))

  Note: Total potential = Delta(sum(pairwise-potential(i, j), (i, j)))

So this will be O(n) rather than O(n^2), and it's fast.
"))

;; The `simulation-step' for `mc-accelerate-mixin' use
;; `quick-potential' for acceleration on potential calculation.
(defmethod simulation-step ((system mc-accelerate-mixin))
  (let ((temperature (system-temperature system))
        (dimension   (system-dimension   system))
        (size        (system-size system))
        (step        (slot-value system 'step)))
    (flet ((quick-potential (i)
             (sum-piter-i* ((j size))
               :reject (lambda (j) (= i j))
               (pairwise-potential system i j))))
      (iter-i* ((i size))
        (let* ((old (copy-array (particle system i)))
               (old-potential   (quick-potential i))
               (move (num-times-vec (random step) (random-vector dimension))))
          (move-particle system i move)
          (let* ((potential (quick-potential i))
                 (-delta-E   (- old-potential potential))
                 (accept? (or (>= -delta-E 0.0)
                              (possible (exp (/ -delta-E temperature))))))
            (unless accept? (setf (particle system i) old))))))))

;; since no old-potential are used, the collect process will calculate
;; the total potential, which might be a little time comsuming, but it
;; is still very fast. 
(defmethod collect-simulation ((system mc-accelerate-mixin))
  (with-slots ((step    simulation-step-counter)
               (counter simulation-collect-counter)
               (collect simulation-collect))
      system
    (setf (aref collect counter) (list step (potential system)))
    (incf counter)))
