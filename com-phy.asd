(defsystem #:com-phy
  :author ("凉凉")
  :version "0.1"
  :licence "WTFPL"
  :description "Computational Physics: MC and MD simulation. "
  :depends-on (ryo gurafu cl-randist numpy-file-format lparallel)
  :components ((:file "package")
               
               ;; The vector operation and random functions
               (:file "vec-op")
               (:file "rand")
               
               ;; The particle system, geometry and pairwise interaction
               (:file "particle-system")
               (:file "geo-mixin")
               (:file "pairwise-mixin")

               ;; The simuation process: MC, MD
               (:file "simulation-mixin")
               (:file "post-process")
               
               (:file "monte-carlo")
               (:file "modular-dynamics")
               (:file "verlet")

               ;; Some 2D system for simulation
               (:file "2d-period-system")))
