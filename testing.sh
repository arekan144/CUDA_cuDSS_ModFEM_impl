#!/bin/bash
for NAME in $@;
do	
if [ ! -e $NAME\_binary ] && [ "$NAME" ] && [ -e ./matrixes/$NAME.txt ]; then 
	./main.exe s ./matrixes/$NAME.txt; 
fi
if [ ! -e $NAME\_binary ] && [ -z "$NAME" ]; then 
	echo "Empty NAME"; 
fi
if [ ! -e $NAME\_binary ] && [ ! -e ./matrixes/$NAME.txt ]; then  
	echo "There's no matrix file in matrixes under name $NAME.txt"; 
fi
if [ -e $NAME\_binary ]; then 
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m g -e save > modfem_crs_61838_binary_diff.txt;
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m s -e load -s .\modfem_crs_61838_binary_solved >> modfem_crs_61838_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m h -e load -s .\modfem_crs_61838_binary_solved >> modfem_crs_61838_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m spd -e load -s .\modfem_crs_61838_binary_solved >> modfem_crs_61838_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m hpd -e load -s .\modfem_crs_61838_binary_solved >> modfem_crs_61838_binary_diff.txt	
	if [ ! -e $NAME\_binary\_solved ]; then
		./main.exe $NAME\_binary b -m g -e save > $NAME\_binary_diff.txt;
	else
		./main.exe $NAME\_binary b -m g -e load -s $NAME\_binary\_solved > $NAME\_binary_diff.txt;
	fi
	./main.exe $NAME\_binary b -m s -e load -s $NAME\_binary\_solved >> $NAME\_binary_diff.txt;
       	./main.exe $NAME\_binary b -m h -e load -s $NAME\_binary\_solved >> $NAME\_binary_diff.txt; 
       	./main.exe $NAME\_binary b -m spd -e load -s $NAME\_binary\_solved >> $NAME\_binary_diff.txt; 
       	./main.exe $NAME\_binary b -m hpd -e load -s $NAME\_binary\_solved >> $NAME\_binary_diff.txt; 
fi
done	

#make test NAME=modfem_crs_1210
#make test NAME=modfem_crs_8379
#make test NAME=modfem_crs_61838
#make test NAME=modfem_crs_476912
#make test NAME=modfem_crs_3745223
