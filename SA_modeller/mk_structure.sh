#!/usr/bin/bash

cd ~/SA_modeller/distance_matrix/
for var in `ls *.npy`
do
	mkdir ~/SA_modeller/project_file
	cp -r ~/SA_modeller/original ~/SA_modeller/project_file/${var%.npy}
	cp ~/SA_modeller/distance_matrix/$var ~/SA_modeller/project_file/${var%.npy}/predict/input.npy
	cd ~/SA_modeller/project_file/${var%.npy}
#Simulated annealing Start
	cd anneal
	make clean
	make
	cd ../predict
	./conv.py
	./run.sh
	./mirror.py
	echo "done ${var%.npy}"
#Attach Glycine
	cp ~/SA_modeller/project_file/${var%.npy}/predict/final.pdb ~/SA_modeller/project_file/${var%.npy}/attach_glycine/
	cd ~/SA_modeller/project_file/${var%.npy}/attach_glycine/
	module load modeller/10.1
	mod10.1 align.py
	cp ~/SA_modeller/project_file/${var%.npy}/attach_glycine/output.B99990001.pdb ~/SA_modeller/project_file/${var%.npy}.pdb
#Attach Glycine_mirror   
	cp ~/SA_modeller/project_file/${var%.npy}/predict/final_mirror.pdb ~/SA_modeller/project_file/${var%.npy}/attach_glycine/
	cd ~/SA_modeller/project_file/${var%.npy}/attach_glycine/
	module load modeller/10.1
	mod10.1 align_m.py
	cp ~/SA_modeller/project_file/${var%.npy}/attach_glycine/output_m.B99990001.pdb ~/SA_modeller/project_file/${var%.npy}_m.pdb
    echo "fin ${var%.npy}"
done
