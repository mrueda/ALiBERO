#!/bin/bash

##############################################
# ALiBERO Installation script                #
# You must change the 3 variables below and  #
# execute this file                          # 
##############################################

# Set the path to ALiBERO exe dir
alibero='/pro/alibero'

# Set your local ICMHOME (exe is /pro/icm/icm/icm)
l_icm='/pro/icm/icm'

# Set your remote ICMHOME (e.g., SDSC)
r_icm='/home/cxedwards/icm/icmd'

echo "Installing ALiBERO..."

##############################################
# Substitutions start here                   #
##############################################

# alibero
sed -e "s#/pro/alibero#$alibero#" \
    -e "s#/pro/icm/3.8-4#$l_icm#" alibero > alibero.bck
mv alibero.bck alibero
chmod +x alibero

# ICM.pm
sed "s#/pro/alibero#$alibero#" ALiBERO/ICM.pm > ALiBERO/ICM.pm.bck
mv ALiBERO/ICM.pm.bck ALiBERO/ICM.pm

# NMA.pm
sed "s#/pro/alibero#$alibero#" ALiBERO/NMA.pm > ALiBERO/NMA.pm.bck
mv ALiBERO/NMA.pm.bck ALiBERO/NMA.pm

# PBS.pm
sed "s#/home/cxedwards/icm/icmd#$r_icm#" ALiBERO/PBS.pm > ALiBERO/PBS.pm.bck
mv ALiBERO/PBS.pm.bck ALiBERO/PBS.pm

# ICM scripts
for macro in MACROS/*.icm
do
  sed  "s#/pro/icm/icm#$l_icm#" $macro > $macro.bck
  mv $macro.bck $macro
done

# ICM scripts in "test" dir
for macro in test/INPUT/MACROS/*.icm
do
  sed  "s#/pro/icm/icm#$l_icm#" $macro > $macro.bck
  mv $macro.bck $macro
done


##############################################
# Substitutions end here                     #
##############################################

# Fin
echo "Succeed!"
