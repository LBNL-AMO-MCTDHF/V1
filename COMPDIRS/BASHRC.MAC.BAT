
echo
echo
echo "WARNING WARNING WARNING"
echo "This edits your ~/.bashrc file."
echo "It will back it up as ~/BACKUP.bashrc but PLEASE MAKE YOUR OWN BACKUP"
echo "Use this script AT YOUR OWN RISK"
echo "This script is for MAC and will put BIN.ecs.hermnorm.mac in your path"
echo 
echo "Hit enter to continue"
read var
echo "Are you sure?  Hit enter to continue"
read var


cp -p ~/.bashrc ~/BACKUP.bashrc

echo >> ~/.bashrc
echo "#### LBNL-AMO-MCTDHF ###" >> ~/.bashrc
echo "export PATH=$PWD/BIN.ecs.hermnorm.mac:\$PATH" >> ~/.bashrc
echo "########################" >> ~/.bashrc
echo >> ~/.bashrc

echo
echo "ok.  chmctdhf_sinc, chmctdhf_atom, and chmctdhf_diatom"
echo "should be in your path after you start a new terminal."
echo





