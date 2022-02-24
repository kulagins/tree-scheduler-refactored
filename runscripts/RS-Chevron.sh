 #!/bin/bash
 
clusters="36proc_0,5-1-1,5-3.json hom-27proc-1.json hom-18proc-1,5.json hom-9proc-3.json 36proc-1-3.json hom-36proc-1.json hom-18proc-3.json"


echo "matrix"
pwd
for cluster in $clusters; do 
    echo "!! $cluster"
    ./main ./data/trees/ treesListChevron.txt ./clusters/36-proc/$cluster 0 0 1      
done


