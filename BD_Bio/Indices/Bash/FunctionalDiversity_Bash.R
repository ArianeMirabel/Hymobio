cd /usr/local/data/Amobio/AMOHYMOBIO/
  ls -lai
chmod 777 -R *
  find . -type f \( -iname "*.R" \) -exec dos2unix {} +
  
  nohup ./FunctionalDiversity_Bash_Fish1.R &
  nohup ./FunctionalDiversity_Bash_Fish2.R &
  nohup ./FunctionalDiversity_Bash_Fish3.R &
  nohup ./FunctionalDiversity_Bash_Fish4.R &
  nohup ./FunctionalDiversity_Bash_Fish5.R &
  nohup ./FunctionalDiversity_Bash_Macroinv1.R &
  nohup ./FunctionalDiversity_Bash_Macroinv2.R &
  nohup ./FunctionalDiversity_Bash_Macroinv3.R &
  nohup ./FunctionalDiversity_Bash_Macroinv4.R &
  nohup ./FunctionalDiversity_Bash_Macroinv5.R &
  nohup ./FunctionalDiversity_Bash_Diatom1.R &
  nohup ./FunctionalDiversity_Bash_Diatom2.R &
  nohup ./FunctionalDiversity_Bash_Diatom3.R &
  nohup ./FunctionalDiversity_Bash_Diatom4.R &
  nohup ./FunctionalDiversity_Bash_Diatom5.R 
  
  
  ctrl Z => stop stuff
bg