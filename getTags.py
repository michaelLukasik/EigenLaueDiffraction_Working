def getFolderAngleTag(isManualConfig, cfg):
    if not isManualConfig:
        if cfg.running == "*THETA_":
            str_folderTag = "varryingTheta_"
            str_angle = cfg.running+cfg.phi+cfg.psi
        elif cfg.running == "*PHI_":
            str_folderTag = "varryingPhi_"
            str_angle = cfg.theta+cfg.running+cfg.psi
        elif cfg.running == "*PSI_":
            str_folderTag = "varryingPsi_"
            str_angle = cfg.theta+cfg.phi+cfg.running
    else:
        str_folderTag = "manual_"
        str_angle = cfg.theta+cfg.phi+cfg.psi
    
    return(str_folderTag, str_angle)

def parseAngleTag(file):
    str_theta = file.split("THETA_")[0].split("_")[-1] + "THETA_"
    str_psi = file.split("PSI_")[0].split("_")[-1] + "PSI_" 
    str_phi = file.split("PHI_")[0].split("_")[-1] + "PHI_"
    return str_phi, str_psi, str_theta


def getFullCSVTag(str_angleTag, cfg):
    return cfg.xpos + cfg.dzdy + cfg.len + str_angleTag + cfg.nx + cfg.ny + cfg.nz + cfg.Lambda + cfg.centering + cfg.name


def getFullSaveTag(str_folderTag, cfg):
    return cfg.xpos + cfg.dzdy + cfg.len + str_folderTag + cfg.nx + cfg.ny + cfg.nz + cfg.Lambda + cfg.centering + cfg.name