(in-package :phy)

;; ========== Geometry Mixin ==========

(defclass geo-mixin ()
  ((lengths :reader system-lengths))
  (:documentation
   "The `geo-mixin' will calculate the geometry of the particle system.

It would take care the following two things:
1. Calculate the `lengths' of particle system;
2. Initialize the `config' (using `:scheme' keyword).

The geometry will be a list of each dimension length stored in `lengths'.
The length will be calculated by (size / density)^(1 / dimension).

The `config' will be generated in `geo-mixin'.
"))

;; ========== Initialize Config Scheme ==========
;; The `geo-mixin' will take care of the config when initialize
;; the `particle-system'. 

(defun %make-grid-config (n lengths)
  "Make `n' particles config evenly destributed on geometry `lengths'.

On each dimension of `lengths', will distribute k particles:

    k = ceiling(n^(1 / dimension))

The particles will be considered as k-number and filled from bottom left
to top right in sequence, like this:

    7 8
    4 5 6
    1 2 3

The distance between two particles in one dimension is:

    grid-size = dim-length / k

These particles will be centered with an offset:

    grid-offset = (dim-length - grid-size * (k - 1)) / 2
"
  (let* ((dim (length lengths))
         (k   (ceiling (expt n (/ 1.0 dim)))))
    (flet ((->k-num (idx)
             (loop for i below dim
                   collect (mod idx k)
                   do (setf idx (floor idx k))))
           (->grid (i dim-length)
             (let* ((grid-size   (/ dim-length k))
                    (grid-offset (/ (- dim-length (* grid-size (1- k))) 2.0)))
               (+ grid-offset (* grid-size i)))))
      (make-array (list n)
                  :initial-contents
                  (collect-i* ((idx n))
                    (map 'vector #'->grid (->k-num idx) lengths))))))

(defun %make-rand-config (n lengths
                          &key (relax 0.9) (min-distance 1.0) (try 3))
  "Make `n' particles config randomly distributed on geometry `lengths'.

Make sure each particls should be `min-distance' from other particles.
If try time is larger than `try', scale `min-distance' by `relax'.
"
  (let ((configs ()))
    (iter-i* ((i n))
      (loop with dist = min-distance
            for p-config
              = (loop for try-time below try
                      for rand = (map 'vector #'random lengths)
                      if (flet ((accept? (p) (> (norm (vec-sub-vec p rand)) dist)))
                           (lparallel:pevery #'accept? configs))
                        return rand)
            if p-config
              return (push p-config configs)
            else
              do (setf dist (* dist relax))))
    (make-array (list n) :initial-contents configs)))

(defmethod initialize-instance :after
    ((system geo-mixin) &key (scheme :grid))
  (with-slots (lengths config) system
    (%calculate-lengths system)
    (ecase scheme
      ((:grid :grided)
       (setf config (%make-grid-config (system-size system) lengths)))
      ((:rand :random)
       (setf config (%make-rand-config (system-size system) lengths)))
      ;; skip the config init, should manually set the config
      (:skip-init t))))

;; ========== General Protocol ==========

(defgeneric distance (system i j)
  (:documentation "Return the vector for displacement from particle `j' to `i'. "))

(defgeneric %calculate-lengths (system)
  (:documentation "Calculate the system coordinates `lengths'. "))

(defgeneric constrain-displacement (system displacement)
  (:documentation "Constrain the system. "))

(defmethod constrain-displacement ((system geo-mixin) displacement)
  (identity displacement))

(defmethod %calculate-lengths ((system geo-mixin))
  (setf (slot-value system 'lengths)
        (make-list (system-dimension system)
                   :initial-element
                   (expt (/ (system-size system) (system-density system))
                         (/ 1.0 (system-dimension system))))))

;; ========== Implementation ==========

(defmethod distance ((system geo-mixin) i j)
  (vec-sub-vec (particle system i) (particle system j)))

;; ========== Periodic Geometry Mixin ==========

(defclass period-geo-mixin (geo-mixin)
  (dist-trim
   geo-trim)
  (:documentation
   "Assuming the particle-system geometry is periodic. "))

(defmethod initialize-instance :after
    ((system period-geo-mixin) &key)
  (let ((cutoff (system-cutoff system)))
    (setf (slot-value system 'dist-trim)
          (lambda (v l)
            (let ((thres (- l cutoff))
                  (sign  (if (> v 0) -1 1))
                  (absv  (abs v)))
              (if (> absv thres) (* sign (- l absv)) v))))
    (setf (slot-value system 'geo-trim)
          (lambda (v l) (mod v l)))))

(defmethod distance :around ((system period-geo-mixin) i j)
  (map 'vector (slot-value system 'dist-trim)
       (call-next-method) (system-lengths system)))

(defmethod move-particle :after ((system period-geo-mixin) i displacement)
  (declare (ignore displacement))
  (with-slots (geo-trim lengths) system
    (setf (particle system i)
          (map 'vector geo-trim (particle system i) lengths))))

(defmethod constrain-displacement ((system period-geo-mixin) displacement)
  (with-slots (geo-trim lengths) system
    (map 'vector geo-trim displacement lengths)))
