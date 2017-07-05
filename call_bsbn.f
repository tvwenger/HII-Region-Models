C call_bsbn.f
C Routine to call bsbn.f with user supplied temperature and density
C
C Syntax: ./call_bsbn.o temperature density
C
C Trey Wenger - May 2016
      PROGRAM main
      DOUBLE PRECISION :: temperature, density
      INTEGER :: num_args
      CHARACTER(LEN=12) :: s 
      num_args = command_argument_count()
      IF (num_args < 2) THEN
        STOP "./call_bsbn.o temperature density"
      ENDIF
      CALL get_command_argument(1,s)
      READ (s,*) temperature
      CALL get_command_argument(2,s)
      READ (s,*) density
      CALL bsbn(temperature,density)
      END PROGRAM MAIN
