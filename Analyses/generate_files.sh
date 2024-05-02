#!/bin/bash

# Define the number of samples
n=3

# Loop over the number of samples
for ((i=1;i<=n;i++)); do
    # Create a new .m file with the sample number replaced
    sed "s/sample_n = 1;/sample_n = $i;/" sensitivity_analysis.m > sensitivity_analysis_$i.m

    # Create a new .sh file with the sample number replaced      
    echo "#!/bin/bash" > sensitivity_analysis_$i.sh
    echo "nohup matlab -nojvm -r 'sensitivity_analysis_$i;exit'" >> sensitivity_analysis_$i.sh

    # Make the new .sh file executable
    chmod +x sensitivity_analysis_$i.sh
done

