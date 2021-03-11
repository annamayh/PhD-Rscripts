## Bash loop to run multiple simulations in SLiM through ubuntu ##


#!/bin/bash

for i in {1..20}; ##sets the number of iterations of the simulation you want

do
slim slim_scripts/Neutral_m2.slim ##opens slim and directs to slim script used
mv /mnt/c/Users/s1881212/ubuntu/SLiM_outputs/Neutral_model_outputs/test_neutral_loop /mnt/c/Users/s1881212/ubuntu/SLiM_outputs/Neutral_model_outputs/neutral_m2_sim${i} 
##renames output file to seed number so file doesnt get overwirtten

done

### Note: make sure you save the file in the SLiM script and make sure its the same path and file name as in line 10
