main:
	gfortran -o bctest_par.exe ./main.f90 ./global_vars.f90 ./boundary_conditions.f90
mainpar:
	mpifort -o bctest_par.exe ./global_vars.f90 ./boundary_conditions.f90 ./outputs.f90 ./main.f90 
clear:
	rm ./outputs/*

