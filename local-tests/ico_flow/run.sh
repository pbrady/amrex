exe="../main2d.gnu.MPI.ex ../inputs"

n=32
mkdir -p $n; cd $n;
$exe n_cell=$n max_grid_size=32 timestep.convection_cfl=0.2 timestep.maxtime=30.0 > out.$n&
cd ..

n=64
mkdir -p $n; cd $n;
$exe n_cell=$n max_grid_size=32 timestep.maxtime=30.0 > out.$n&
cd ..
n=128
mkdir -p $n; cd $n;
$exe n_cell=$n max_grid_size=32 timestep.maxtime=30.0 > out.$n&
cd ..
n=256
mkdir -p $n; cd $n
echo $n
mpirun -np 4 $exe n_cell=$n max_grid_size=64 timestep.maxtime=30.0 > out.$n
cd ..
n=512
mkdir -p $n; cd $n
echo $n
mpirun -np 6 $exe n_cell=$n max_grid_size=64 timestep.maxtime=30.0 timestep.viscous_cfl=3.0 > out.$n
