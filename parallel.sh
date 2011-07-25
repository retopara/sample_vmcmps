#the compiler
date

/opt/intel/impi/openmpi/bin/mpirun -np $1 -machinefile myhostfile ./vmps.exe

echo -e "\\n"
date
