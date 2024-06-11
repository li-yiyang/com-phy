(in-package :phy)

(defun average-list (data)
  "Calculate the average of the `data'.
Return data list element would be like (x <y>). 

The `data' element should be like (x y).
The average <y> will be average value of the previous y.
"
  (loop with sum = 0
        for count from 1
        for (x y) in data
        for average = (/ sum count)
        do (incf sum y)
        collect (list x average)))

;; This is a little slow for the large `n'. 
(defun n-weighted-average-list
    (data n &optional (weights (make-list n :initial-element (/ 1.0 n))))
  "Return the `n' weighted average data list.
Return data list element would be like (<x> <y>).

The average <x>, <y> will be sum(w * y) / sum(w).
The `weights' should be a `n' length list, otherwise, warning.
"
  (unless (length= weights n)
    (warnf "The `weights' ~a should be in `n' ~d length. " weights n))
  (loop with sum-w = (reduce #'+ weights)
        for i upto (- (length data) n)
        for dat on data
        collect (loop with sum-x = 0
                      with sum-y = 0
                      for offset below n
                      for weight in weights
                      for (x y) in dat
                      do (incf sum-x (* x weight))
                      do (incf sum-y (* y weight))
                      finally (return (list sum-x (/ sum-y sum-w))))))

;; This is the fast version of `n-weighted-average-list'. 
(defun n-average-list (data n)
  "Return the average `n' list of data.
Return data list element would be like (x <y>).

Same like (n-weighted-average-list data n), but faster.
"
  (let ((weight (/ 1.0 n)))
    (multiple-value-bind (x-buffer y-buffer)
        (loop with x-sum = (* 2.0 (first  (first  data)))
              with y-sum = (* 2.0 (second (second data)))
              for i from 2 upto n
              for (x y) in data
              do (incf x-sum x)
              do (incf y-sum y)
              finally (return (values x-sum y-sum)))
      (loop for (sub-x sub-y) in data
            for (add-x add-y) in (nthcdr n data)
            do (incf x-buffer (- add-x sub-x))
            do (incf y-buffer (- add-y sub-y))
            collect (list (* weight x-buffer) (* weight y-buffer))))))

(defun n-fluctuations (data n)
  "Return the fluctuations with `n' sampling for `data'.
Return data list element would like (x (<y^2> - <y>)).
"
  (let ((weight (/ 1.0 n)))
    (flet ((square (x) (* x x)))
      (multiple-value-bind (x-buffer y-buffer y2-buffer)
          (loop with x-sum  = (* 2.0 (first  (first  data)))
                with y-sum  = (* 2.0 (second (second data)))
                with y2-sum = (* 2.0 (square (second (second data))))
                for i from 2 upto n
                for (x y) in data
                do (incf x-sum x)
                do (incf y-sum y)
                do (incf y2-sum (square y))
                finally (return (values x-sum y-sum y2-sum)))
        (loop for (sub-x sub-y) in data
              for (add-x add-y) in (nthcdr n data)
              do (incf x-buffer (- add-x sub-x))
              do (incf y-buffer (- add-y sub-y))
              do (incf y2-buffer (- (square add-y) (square sub-y)))
              collect (list (* weight x-buffer)
                            ;; dU^2 = <U^2> - <U>^2
                            (- (* weight y2-buffer)
                               (square (* weight y-buffer)))))))))

(defun average-data (data &optional (start 0))
  "Start at `start' for average `data', element is (x y). "
  (loop for dat on data
        if (>= (first (first dat)) start)
          return (/ (reduce #'+ dat :key #'second) (length dat))
        finally (return 0.0)))

(defun fluctuation-data (data &optional (start 0) (samples))
  "Start at `start' for fluctuation `data', element is (x y). "
  (loop for dat on data
        if (>= (first (first dat)) start)
          return (let* ((ys     (mapcar #'second dat))
                        (size   (length ys))
                        (avg-y  (/ (reduce #'+ ys) size))
                        (avg-y2 (/ (reduce #'+ (mapcar (lambda (x) (* x x)) ys))
                                   size)))
                   (- avg-y2 (* avg-y avg-y)))
        finally (return 0.0)))

(defun rdf-histogram (system cutoff bins)
  "Return the system RDF histogram. "
  (let ((data    (make-array (list bins) :initial-element 0))
        (grid    (/ cutoff bins))
        (size    (system-size system)))
    (iter-i* ((i j) (j size))
      (let ((r (norm (distance system i j))))
        (when (< r cutoff) (incf (aref data (truncate r grid))))))
    (flet ((square (x) (* x x)))
      (piter-i* ((i bins))
        (let ((vb (* pi (- (square (1+ i)) (square i)) (square grid))))
          (setf (aref data i) (/ (aref data i) (* size vb))))))
    (values data grid)))

(defun %radial-distribution-function
    (system &key (cutoff (system-cutoff system)) (bins 1000))
  "Return a lambda function for RDF g(r) of `config' and `lengths'. "
  (multiple-value-bind (data grid)
      (rdf-histogram system cutoff bins)
    (lambda (r) (if (< r cutoff) (aref data (truncate r grid)) 0))))

(defun simulation-average-potential (system &optional (start 0))
  "Calculate the average potential <U> / N. "
  (/ (average-data (map 'list #'identity (simulation-collect system)) start)
     (system-size system)))

(defun rdf-integrate (system function &optional (bins 1000))
  "Integrate RDF on the 2D plane.

  integrate(RDF(r) * function(r) * pi * r * dr from 0 to rc)
"
  (multiple-value-bind (rdf-histogram grid)
      (rdf-histogram system (system-cutoff system) bins)
    (flet ((square (x) (* x x)))
      (sum-piter-i* ((i (length rdf-histogram)))
        (let ((rdf (aref rdf-histogram i)))
          (if (zerop rdf) 0.0
              (* pi (- (square (1+ i)) (square i)) (square grid)
                 rdf (funcall function (* i grid)))))))))

(defun rdf-average-potential (system &optional (bins 1000))
  "Calculate the average potential from radial distribution function.

    u = rdf-integrate( u(r) )
"
  (flet ((potential (r) (pairwise-potential* system r)))
    (rdf-integrate system #'potential bins)))

(defun rdf-pressure (system &optional (bins 1000))
  "Calculate the pressure form radial distribution function.

  p = density * temperature + rdf-integrate( dot(f(r), hat(r)) )
"
  (flet ((f-project (r) (dot (n-vec 1 1) (pairwise-force* system (n-vec 1 r)))))
    (rdf-integrate system #'f-project bins)))

