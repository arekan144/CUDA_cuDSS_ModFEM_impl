Param(
    [string[]] $MatrixNames
)
foreach ($name in $MatrixNames) {
    if ($name.length -eq 0) {
		Write-Output "Empty name";
		continue;
	}
	
	if (-Not (Test-Path $name'_binary' -PathType Leaf)) {
		$nameConcat = [system.String]::Join("",$('./matrixes/',$name,'.txt'));
		if (-Not (Test-Path $nameConcat -PathType Leaf)) {
			Write-Output "No matrix file $nameConcat";
			continue;
		}
		./Cuda_test_final.exe $nameConcat -t S;
	}
	#Write-Output $name'_binary'
	if (-Not (Test-Path $name'_binary_solved' -PathType Leaf)) {
		./Cuda_test_final.exe $name'_binary' b -m g -e save > $name'_binary_diff.txt';
	}
	else {
		./Cuda_test_final.exe $name'_binary' b -m g -e load -s $name'_binary_solved' > $name'_binary_diff.txt';
	}
	./Cuda_test_final.exe $name'_binary' b -m s -e load -s $name'_binary_solved' >> $name'_binary_diff.txt';
	./Cuda_test_final.exe $name'_binary' b -m h -e load -s $name'_binary_solved' >> $name'_binary_diff.txt';
	./Cuda_test_final.exe $name'_binary' b -m spd -e load -s $name'_binary_solved' >> $name'_binary_diff.txt';
	./Cuda_test_final.exe $name'_binary' b -m hpd -e load -s $name'_binary_solved' >> $name'_binary_diff.txt';	
}

#.\Cuda_test_final.exe .\modfem_crs_1210_binary b -m g -e save > modfem_crs_1210_binary_diff.txt;
#.\Cuda_test_final.exe .\modfem_crs_1210_binary b -m s -e load -s .\modfem_crs_1210_binary_solved >> modfem_crs_1210_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_1210_binary b -m h -e load -s .\modfem_crs_1210_binary_solved >> modfem_crs_1210_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_1210_binary b -m spd -e load -s .\modfem_crs_1210_binary_solved >> modfem_crs_1210_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_1210_binary b -m hpd -e load -s .\modfem_crs_1210_binary_solved >> modfem_crs_1210_binary_diff.txt
#Write-Host "modfem_crs_1210 done"
#.\Cuda_test_final.exe .\modfem_crs_8379_binary b -m g -e save > modfem_crs_8379_binary_diff.txt;
#.\Cuda_test_final.exe .\modfem_crs_8379_binary b -m s -e load -s .\modfem_crs_8379_binary_solved >> modfem_crs_8379_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_8379_binary b -m h -e load -s .\modfem_crs_8379_binary_solved >> modfem_crs_8379_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_8379_binary b -m spd -e load -s .\modfem_crs_8379_binary_solved >> modfem_crs_8379_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_8379_binary b -m hpd -e load -s .\modfem_crs_8379_binary_solved >> modfem_crs_8379_diff.txt
#Write-Host "modfem_crs_8379 done"
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m g -e save > modfem_crs_61838_binary_diff.txt;
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m s -e load -s .\modfem_crs_61838_binary_solved >> modfem_crs_61838_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m h -e load -s .\modfem_crs_61838_binary_solved >> modfem_crs_61838_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m spd -e load -s .\modfem_crs_61838_binary_solved >> modfem_crs_61838_binary_diff.txt
#.\Cuda_test_final.exe .\modfem_crs_61838_binary b -m hpd -e load -s .\modfem_crs_61838_binary_solved >> modfem_crs_61838_binary_diff.txt
#Write-Host "modfem_crs_61838 done"
#.\Cuda_test_final.exe .\modfem_crs_476912_binary b -m g -e load -s .\modfem_crs_476912_binary_solved >> modfem_crs_476912_binary_diff.txt;
#.\Cuda_test_final.exe .\modfem_crs_476912_binary b -m s -e load -s .\modfem_crs_476912_binary_solved >> modfem_crs_476912_binary_diff.txt;
#.\Cuda_test_final.exe .\modfem_crs_476912_binary b -m h -e load -s .\modfem_crs_476912_binary_solved >> modfem_crs_476912_binary_diff.txt;
#.\Cuda_test_final.exe .\modfem_crs_476912_binary b -m spd -e load -s .\modfem_crs_476912_binary_solved >> modfem_crs_476912_binary_diff.txt;
#.\Cuda_test_final.exe .\modfem_crs_476912_binary b -m hpd -e load -s .\modfem_crs_476912_binary_solved >> modfem_crs_476912_binary_diff.txt;
#Write-Host "modfem_crs_476912 done"