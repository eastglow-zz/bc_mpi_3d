main:
	gfortran -o nucltest_par.exe ./main.for ./global_vars.for ./boundary_conditions.for
mainpar:
	mpifort -o nucltest_par.exe ./global_vars.for ./boundary_conditions.for ./outputs.for ./put_sphere.for ./main.for -fdefault-real-8 
	mkdir -p ./outputs

clamptest:
	mpifort -o clamptest.exe ./global_vars.for ./boundary_conditions.for ./outputs.for ./clamp_test.for -fdefault-real-8
clear:
	rm ./outputs/*

