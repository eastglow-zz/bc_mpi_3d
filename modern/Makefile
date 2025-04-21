s=./src
fc=mpifort
s0=$(s)/global_vars.for $(s)/M410_bc_modern.for $(s)/M801_pvoutputs_modern.for 
s1=$(s)/main.for

main:
	$(fc) -o main.exe $(s0) $(s1)
	rm ./*.mod 

debug:
	$(fc) -o main_dbg.exe -g -fcheck=all -Wall -Wextra $(s0) $(s1)
	rm ./*.mod