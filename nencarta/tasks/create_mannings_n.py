from pathlib import Path

from nencarta.logger import LOG
from nencarta.workspace import Workspace

def Create_BaseLine_Manning_n_File_ESA(workspace: Workspace) -> Path:
    ManningN = workspace.mannings_n_text_file
    if ManningN.exists():
        LOG.info(f"Manning's n file already exists, skipping creation: {ManningN}")
        return ManningN
    
    ManningN.parent.mkdir(parents=True, exist_ok=True)
    LOG.info('Creating Manning n file: ' + str(ManningN))   
    with open(ManningN,'w') as out_file:
        out_file.write('LC_ID	Description	Manning_n')
        out_file.write(f'\n10	Tree Cover	0.120')
        out_file.write(f'\n20	Shrubland	0.050')
        out_file.write(f'\n30	Grassland	0.030')
        out_file.write(f'\n40	Cropland	0.035')
        out_file.write(f'\n50	Builtup	0.075')     
        out_file.write(f'\n60	Bare	0.030')
        out_file.write(f'\n70	SnowIce	0.030')
        out_file.write(f'\n80	Water	0.030')
        out_file.write(f'\n90	Emergent_Herb_Wet	0.100')
        out_file.write(f'\n95	Mangroves	0.100')
        out_file.write(f'\n100	MossLichen	0.100')

    return ManningN