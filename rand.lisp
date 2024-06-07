;; Documentation
;; =====
;; This part defines the random functions. 

(in-package :phy)

(defun random-vector (dimension)
  "Return a uniformed vector in `dimension' (length = 1). "
  (let ((vec (make-array (list dimension))))
    (dotimes (i dimension (num-times-vec (/ 1.0 (norm vec)) vec))
      (setf (aref vec i) (1- (random 2.0))))))

(defun boltzmann-random-vector (dimension temperature)
  "Return a normal distribution vector in `dimension'. "
  (let ((vec (make-array (list dimension))))
    (dotimes (i dimension vec)
      (setf (aref vec i)
            (cl-randist:random-normal-ziggurat
             0.0d0 (sqrt (coerce temperature 'double-float)))))))

(declaim (inline possible))
(defun possible (p)
  "If it's chance is P(p). "
  (<= (random 1.0) p))
