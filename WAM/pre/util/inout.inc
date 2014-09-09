c inout.inc : include file for DNMI WAM inout.f
c ---------
c
c  mfsize : maximum field size, input and output
C  CCC Modified to fit WAM20 2005-08-05 oyvind.breivik@met.no
c
      parameter (mfsize=105000)
c
c  computer:   switch for computer dependant i/o methodes,
c                 'cray.t3e' - Cray T3E
c                 'cray'     - Cray XMP,YMP,J90,...
c                 'work'     - workstations (SGI,...)
c              currently not used!
c  wordlength: length in bytes for real/integer
c  reclunit:   unit in bytes for recl= in open statements
c
c     character*8 computer
c     integer     wordlength,reclunit
c     parameter (computer='cray.t3e', wordlength=8, reclunit=1)
c     parameter (computer='cray',     wordlength=8, reclunit=1)
c     parameter (computer='work',     wordlength=4, reclunit=4)
c
c  undefined 'high' value (undef) and value when testing (udef)
      real undef,udef
      parameter (undef=+1.e+35)
      parameter (udef=undef*0.9)
c
