#
# builds the file
# utiliser le terminal et faire "chmod 755" sur ce fichier (chein complet ou glisser déposer dan sle terminal) pour autoriser à exécuter le script en cas de déni de permission
#

#!/bin/sh

#source /opt/intel/composerxe/bin/compilervars.sh intel64


DIR="cd "$( dirname "$0" )""
DIR2=""$( dirname "$0" )""
$DIR

gcc -g EUP.cpp -o EUP.out -lstdc++

