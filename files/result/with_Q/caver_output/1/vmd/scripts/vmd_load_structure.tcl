#set dir "/home/miller96/Software/dowserplusplus/dpp_original/dowserplusplus/files/result/with_Q/caver_output/1/inputs"

mol load pdb ../data/cavity_glu_prot_Q10_EM_100kJ_water.pdb

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 

