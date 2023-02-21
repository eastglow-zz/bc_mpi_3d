main:
	gfortran -o bctest_par.exe ./main.for ./global_vars.for ./boundary_conditions.for
mainpar:
	mpifort -o bctest_par.exe ./global_vars.for ./boundary_conditions.for ./outputs.for ./main.for -fdefault-real-8 
clear:
	rm ./outputs/*

