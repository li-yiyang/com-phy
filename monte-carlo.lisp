(in-package :phy)

;; ========== Monte Carlo ==========

(defclass monte-carlo-mixin (simulation-mixin)
  ((step :initform 0.1 :initarg :step-size :reader system-step-size)
   (current-particle :initform 0)
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
    (with-slots ((i current-particle)) system
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
              (setf (particle system i) old))))
      (setf i (mod (1+ i) size)))))

(defmethod lazy-potential ((system monte-carlo-mixin))
  (with-slots (old-potential current-particle) system
    (when (zerop current-particle)
      (setf old-potential (potential system)))
    old-potential))

;; The `simulation-collect' is (step potential)
(defmethod collect-simulation ((system monte-carlo-mixin))
  (with-slots ((step    simulation-step-counter)
               (counter simulation-collect-counter)
               (collect simulation-collect))
      system
    (setf (aref collect counter) (list step (lazy-potential system)))
    (incf counter)))

;; The `plot-simulation' will plot the potential along the MC simulation. 
(defmethod plot-simulation ((system monte-carlo-mixin) out-path
                            &key (x-min 0 x-min-set?) (x-max 0 x-max-set?)
                              (y-min 0 y-min-set?) (y-max 0 y-max-set?)
                              (n-average 1000))
  (with-slots ((counter simulation-collect-counter)
               (collect simulation-collect))
      system
    (let* ((data (collect-i* ((i counter)) (aref collect i)))
           (y-min (if y-min-set? y-min (reduce #'min data :key #'second)))
           (y-max (if y-max-set? y-max (reduce #'max data :key #'second))))
      (with-present-to-file
          (plot plot :margin 20
                     :x-min (if x-min-set? x-min 0)
                     :x-max (if x-max-set? x-max counter)
                     :y-min y-min
                     :y-max y-max
                     :x-label "s" :y-label "V")
          (out-path :width 800 :height 400)
        (add-plot-data plot (line-plot-pane potential :color +大红+)
          data)
        (add-plot-data plot (line-plot-pane average-V :color +水红+)
          (average-list data))
        ;; If tooo few samples, make warning
        (if (> counter n-average)
            (add-plot-data plot (line-plot-pane N-average :color +鹅黄+)
              (n-average-list data n-average))
            (warnf "Simulation collected only ~d samples, on average. " counter))
        (add-plot-legend (plot :position :top-right :padding 0.0)
          ((format nil "V, T = ~,3f, N = ~d, ⍴ = ~,3f"
                   (system-temperature system)
                   (system-size system)
                   (system-density system))
           :color +大红+)
          ("<V> average potential" :color +水红+)
          ((format nil "<V> every ~d" n-average) :color +鹅黄+))))))

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
    (with-slots ((i current-particle)) system
      (flet ((quick-potential ()
               (sum-piter-i* ((j size))
                 :reject (lambda (j) (= i j))
                 (pairwise-potential system i j))))
        (let* ((old (copy-array (particle system i)))
               (old-potential   (quick-potential))
               (move (num-times-vec (random step) (random-vector dimension))))
          (move-particle system i move)
          (let* ((potential (quick-potential))
                 (delta-E   (- potential old-potential ))
                 (accept? (or (< delta-E 0.0)
                              (possible (exp (/ (- delta-E) temperature))))))
            (if accept?
                (incf (slot-value system 'old-potential) delta-E)
                (setf (particle system i) old)))))
      (setf i (mod (1+ i) size)))))

;; since no old-potential are used, the collect process will calculate
;; the total potential, which might be a little time comsuming, but it
;; is still very fast. 
(defmethod collect-simulation ((system mc-accelerate-mixin))
  (with-slots ((step    simulation-step-counter)
               (counter simulation-collect-counter)
               (collect simulation-collect))
      system
    (setf (aref collect counter) (list step (lazy-potential system)))
    (incf counter)))
