(defpackage #:com-phy
  (:use :cl :ryo :gurafu)
  (:nicknames :phy)
  (:import-from :alexandria :copy-array)
  (:documentation
   "Computational Physics: MC and MD simulation.
============================================================

It is the Final Project of UCAS 2024 Computational Physics.
Within this repo, there is the Monte Carlo Simulation and
Modular Dynamics Simulatioin program written in Common Lisp
via CLOS.

Usage
=====
For Monte Carlo Simulation:

    (in-package :com-phy)
    (let ((system (make-instance '2d-lj-period-mc-acc
                                 :temperature 0.728
                                 :density     0.8442
                                 :N           100
                                 :scheme      :grided
                                 :step-size   0.05)))  ; dt for MD
      (run-simulation system 1000
                      :plot-config-every  100
                      :plot-collect-every 100)
      (save-system system save-path))

You may consider to set default read float format as double
to make more precise calculation:

    (setf *read-default-float-format* 'double-float)

Dependency
=====
I'm trying to write these with minimum dependences. It depends
on two libraries I wrote: GURAFU (li-yiyang/gurafu),
RYO (li-yiyang/ryo), which you could found in my Github repo.
Links are below:

  My Github profile: https://github.com/li-yiyang
  GURAFU*:           https://github.com/li-yiyang/gurafu
  RYO:               https://github.com/li-yiyang/ryo

    *Note: GURAFU depends on my another project CL-BDF to draw
    char you need to download that too:
    https://github.com/li-yiyang/cl-bdf

Other libraries should be easily solved via quicklisp.

License
====
           DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
                   Version 2, December 2004
 
Copyright (C) 2024 凉凉 <https://github.com/li-yiyang>

Everyone is permitted to copy and distribute verbatim or modified
copies of this license document, and changing it is allowed as long
as the name is changed.
 
           DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
  TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

 0. You just DO WHAT THE FUCK YOU WANT TO.
"))

(in-package :com-phy)
