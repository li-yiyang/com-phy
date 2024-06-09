(in-package :phy)

(defconstant +max-move+ 2.0
  "If single move is larger than `+max-move+' will trigger warning. ")

;; ========== Verlet Algorithm ==========

(defclass md-verlet-mixin (modular-dynamics-mixin)
  (xn-1-buffer
   xn+1-buffer)
  (:documentation
   "The Verlet Algorithm.

    x(t + dt) = 2 * x(t) - x(t - dt) + a(t) * dt^2 + O(dt^4)

where: x(t + dt), x(t), x(t - dt) will be used for simulation temperary
data, responding to xn+1-buffer, config, xn-1-buffer. 
"))

;; ========== Init Simulation ==========
;; Calculate the x1 for each particle. 

(defmethod init-simulation ((system md-verlet-mixin))
  (with-slots (xn-1-buffer xn+1-buffer delta-t config) system
    (let ((size (system-size system)))
      (setf xn+1-buffer (make-array (list size))
            xn-1-buffer (make-array (list size)))
      (piter-i* ((i size))
        ;; x1 = x0 + (v0 * dt + 0.5 * f(x0) * dt^2)
        (let ((move (vec-plus-vec (num-times-vec delta-t
                                                 (particle-velocity system i))
                                  (num-times-vec (* 0.5 delta-t delta-t)
                                                 (particle-force system i)))))
          (when (> (norm move) +max-move+)
            (warnf "Move ~a tooo large " move))
          (setf (aref xn+1-buffer i)
                (constrain-displacement system
                                        (vec-plus-vec (aref config i) move)))))
      (shiftf xn+1-buffer xn-1-buffer config xn+1-buffer))
    system))

;; Note: here is a little tricky thing: if for the periodic boundary
;; condition, when the xn is larger than boundary length L, it would
;; be constrained to the (xn mod L), therefore, to calculate the
;; (xn - xn-1) would return faulty results. This is same like velocity
;; calculation using no constrained move.
;;
;; Not know how to solve it yet... 
(defmethod simulation-step ((system md-verlet-mixin))
  (with-slots (xn-1-buffer xn+1-buffer config delta-t) system
    (let ((size (system-size system)))
      (piter-i* ((i size))
        ;; xn+1 = xn + (xn - xn-1 + f(xn) * dt^2)
        (let* ((move (vec-plus-vec (num-times-vec (* delta-t delta-t)
                                                  (particle-force system i))
                                   (vec-sub-vec (aref config i)
                                                (aref xn-1-buffer i))))
               (xn+1 (vec-plus-vec (aref config i) move)))
          (when (> (norm move) +max-move+)
            (warnf "Move ~a tooo large " move))
          ;; xn+1 should be constrained by geo-mixin
          (setf (aref xn+1-buffer i) (constrain-displacement system xn+1))
          ;; velocity should be calculated by no constrain move. 
          (setf (particle-velocity system i)
                (num-times-vec (/ 0.5 delta-t)
                               (vec-sub-vec xn+1 (aref xn-1-buffer i))))))
      (shiftf xn+1-buffer xn-1-buffer config xn+1-buffer))))

;; ========== Velocity Verlet ==========

(defclass md-velocity-verlet-mixin (modular-dynamics-mixin)
  (v+1/2-buffer v+1-buffer xn+1-buffer)
  (:documentation
   "Velocity Verlet Algorithm. "))

;; just open the buffer
(defmethod init-simulation ((system md-velocity-verlet-mixin))
  (with-slots (v+1/2-buffer v+1-buffer xn+1-buffer)
      system
    (let ((size (system-size      system)))
      (setf xn+1-buffer  (make-array (list size))
            v+1/2-buffer (make-array (list size))
            v+1-buffer   (make-array (list size))))
    system))

(defmethod simulation-step ((system md-velocity-verlet-mixin))
  (with-slots (delta-t v+1/2-buffer v+1-buffer
               xn+1-buffer
               config velocity-config)
      system
    (let ((size (system-size      system))
          (dim  (system-dimension system)))
      (piter-i* ((i size))
        ;; v(t + 1/2 * dt) = v(t) + 1/2 * a(t) * dt
        (setf (aref v+1/2-buffer i)
              (vec-plus-vec (aref velocity-config i)
                            (num-times-vec (* 0.5 delta-t)
                                           (particle-force system i))))
        ;; x(t + dt) = x(t) + v(t + 1/2 * dt) * dt
        (setf (aref xn+1-buffer i)
              (vec-plus-vec (aref config i)
                            (num-times-vec delta-t (aref v+1/2-buffer i)))))
      (flet ((xn+1-particle-force (i)
               (sum-iter-i* ((j size))
                 :reject (lambda (j) (= j i))
                 :sum-init (make-array (list dim) :initial-element 0.0)
                 :sum-method #'vec-plus-vec
                 (pairwise-force* system
                                  (vec-sub-vec (aref xn+1-buffer i)
                                               (aref xn+1-buffer j))))))
        (piter-i* ((i size))
          ;; v(t + dt) = v(t + 1/2 * dt) + 1/2 * a(t + dt) * dt
          (setf (aref v+1-buffer i)
                (vec-plus-vec (aref v+1/2-buffer i)
                              (num-times-vec (* 0.5 delta-t)
                                             (xn+1-particle-force i))))))
      (piter-i* ((i size))
        ;; xn+1 should be constrained by geo-mixin
        (setf (aref xn+1-buffer i)
              (constrain-displacement system (aref xn+1-buffer i))))
      (shiftf xn+1-buffer config xn+1-buffer)
      (shiftf v+1-buffer velocity-config v+1-buffer))
    system))
