;; Documentation
;; =====
;; This part defines how the vector is calculated.

(in-package :phy)

(declaim (inline vec-plus-vec vec-sub-vec norm num-times-vec n-vec dot))
(defun vec-plus-vec (a b)
  "Add wto vector `a' and `b'. "
  (map 'vector #'+ a b))

(defun vec-sub-vec (a b)
  "Return vector `a' substract vector `b'. "
  (map 'vector #'- a b))

(defun norm (vec)
  "Return the Eucilid length of `vec'. "
  (sqrt (reduce #'+ (map 'list #'* vec vec))))

(defun num-times-vec (num vec)
  "Times a `num' to `vec'. "
  (map 'vector (lambda (v) (* num v)) vec))

(defun n-vec (n val)
  "A `n' dimension vector with each element of value `val'. "
  (make-array (list n) :initial-element val))

(defun dot (vec-a vec-b)
  "Dot two vector `vec-a' and `vec-b'. "
  (reduce #'+ (map 'list #'* vec-a vec-b)))
