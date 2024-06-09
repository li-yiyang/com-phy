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

(defun %radial-distribution-function (system)
  "Return a lambda function for RDF of `config' and `lengths'.
The lambda should be called with lambda list (r &optional (dr 0.01)).
"
  (let* ((size    (system-size system))
         (v/n*n-1 (/ (* 2.0 (reduce #'* (system-lengths system)))
                     (* size (1- size)))))
    (flet ((theta (x) (if (> x 0) 1 0)))
      (lambda (r &optional (dr 0.01))
        (let ((vr (* pi (+ (* 2.0 r) dr) dr)))
          (* (/ v/n*n-1 vr)
             (sum-piter-i* ((j i) (i size))
               (let ((rij (norm (distance system i j))))
                 (* (theta (- rij r)) (theta (- (+ r dr) rij)))))))))))
