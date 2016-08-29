These are some python scripts for generating analytic solutions of
Rad2Tshock.  They are not very user friendly.  Feel free to modify it.
Here is what I have to do.

* run initincgs.py to generate a cPickle file units-shock.p

* run RadShock.py to generate a cPickle file shock.p

  (The defaults parameters for these two scripts are for the M2 test.
  For M5, we need set M0=5.)

* Then I use another python script to generate plots (e.g., paper-M2-mg.py).
