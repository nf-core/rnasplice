process GETLIST {

      conda (params.enable_conda ? "bioconda::csvtk=0.23.0" : null)
      container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
          'https://depot.galaxyproject.org/singularity/csvtk:0.25.0--h9ee0642_0' :
          'quay.io/biocontainers/csvtk:0.25.0--h9ee0642_0' }"
    input:
    path samplesheet
 
    output:
    path "suppa_split_tpm.txt" , emit: tpm_list
    path "suppa_split_psi.txt" , emit: psi_list
    path "versions.yml", emit: versions

    shell: 
   $/
   csvtk cut -f "condition","sample" !{samplesheet} > suppa_split1.txt #Retain only 2 columns
   tail -n +2  suppa_split1.txt > suppa_split2.txt # Remove header
   sort -k 1 suppa_split2.txt > suppa_split3.txt #Sort to be able to group by condition
   #Group by condition
   awk -F ',' '
   $1==x{
    printf ",%s", $2
    next
   }
   {
    x=$1
    printf "\n%s.tpm|%s", $1, $2
   }
   ' suppa_split3.txt > suppa_split4.txt
   tail -n +2  suppa_split4.txt > suppa_split5.txt # Remove the empty first line 
   a=$(awk -F '|' '{printf "%s ", $2}' suppa_split5.txt) # Fetch the sample column 
   b=$(awk -F '|' '{printf "%s ", $1}' suppa_split5.txt) # Fetch the condition column 
   echo -n "$$a" "$$b" >suppa_split_tpm.txt # Create a string with samples first and condition next
   sed 's/.tpm/.psi/g' suppa_split_tpm.txt > suppa_split_psi.txt # Change the suffix tpm to psi and store in another file 
   
   cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: $(echo $( csvtk version | sed -e 's/csvtk v//g' ))
   /$
}