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
  ()
  (:documentation
   ""))

