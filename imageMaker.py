import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cv2

from getTags import parseAngleTag
import tools

def readCSV(file):
    df = pd.read_csv(file, header = None , dtype = np.float64, na_values = "-nan(ind)")
    df = df.fillna(0) ## Centering the lattice results in the center-point often having a divide by zero error, this fixes that possibility
    sim_wfRe =  df[0]
    sim_wfIm =  df[1]
    sim_xpos = df[2]
    sim_ypos = df[3]
    sim_zpos = df[4]
    intensity = np.sqrt(sim_wfRe**2 + sim_wfIm**2)
    return [sim_wfRe, sim_wfIm, sim_xpos, sim_ypos, sim_zpos, intensity]


def makeRaw(simData, file, savePath, saveFolder, cfg):
    str_phi, str_psi, str_theta = parseAngleTag(file)
    wall_len = float(cfg.len.split("LEN_")[0])
    dzdy = float(cfg.dzdy.split("DZDY_")[0])  
    sim_wfRe, sim_wfIm, sim_xpos, sim_ypos, sim_zpos, intensity = simData[0],simData[1],simData[2],simData[3],simData[4],simData[5],

    intensity = np.asarray(intensity) / np.linalg.norm(intensity) 
    intensity = np.reshape(intensity,(int(wall_len/dzdy), int(wall_len/dzdy)))
    
    fig = plt.figure(frameon=False, figsize=(10,10))
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(np.abs(intensity),  cmap = cfg.cmap, aspect='auto')
    plt.axis('off')
    fig.savefig(savePath + saveFolder + "/raw/" + file.split("\\")[-1][:-4]+".png",bbox_inches='tight', pad_inches=0,transparent=True,)
    plt.clf()

def makeRawLabeled(simData, file, savePath, saveFolder, cfg):
    str_phi, str_psi, str_theta = parseAngleTag(file)
    wall_len = float(cfg.len.split("LEN_")[0])
    dzdy = float(cfg.dzdy.split("DZDY_")[0])  
    
    sim_wfRe, sim_wfIm, sim_xpos, sim_ypos, sim_zpos, intensity = simData[0],simData[1],simData[2],simData[3],simData[4],simData[5]
    intensity = np.asarray(intensity) / np.linalg.norm(intensity) 
    intensity = np.reshape(intensity,(int(wall_len/dzdy), int(wall_len/dzdy)))
    
    plt.subplot()
    plt.imshow(abs(intensity),origin = 'lower', cmap = cfg.cmap, aspect = 'auto', interpolation = "gaussian",\
             extent = [np.min(sim_ypos),np.max(sim_ypos),np.min(sim_zpos),np.max(sim_zpos)])
    plt.title(str_theta[:-6]+r'$^{\circ}\theta$ '+str_phi[:-4]+r'$^{\circ}\phi$ '+str_psi[:-4]+r'$^{\circ}\psi$ ' )
    plt.xlim(np.min(sim_ypos/cfg.zoom),(np.max(sim_ypos/cfg.zoom)))
    plt.ylim(np.min(sim_zpos/cfg.zoom),(np.max(sim_zpos/cfg.zoom)))
    plt.savefig(savePath + saveFolder + "/" + file.split("\\")[-1][:-4]+".png")
    plt.clf()
    
def makeFFT(simData, file, savePath, saveFolder, cfg):
    str_phi, str_psi, str_theta = parseAngleTag(file)
    wall_len = float(cfg.len.split("LEN_")[0])
    dzdy = float(cfg.dzdy.split("DZDY_")[0])  

    sim_ypos, sim_zpos, intensity = simData[3],simData[4],simData[5]
    intensity_fft = np.fft.fft(intensity)
    intensity_fft = np.log(intensity_fft)
    intensity_fft = np.asarray(intensity_fft) / np.linalg.norm(intensity_fft)
    intensity_fft = np.reshape(intensity_fft,(int(wall_len/dzdy), int(wall_len/dzdy)))

    
    plt.subplot()
    plt.imshow(abs(intensity_fft),origin = 'lower', cmap = cfg.cmap, aspect = 'auto', interpolation = "gaussian",\
             extent = [np.min(sim_ypos),np.max(sim_ypos),np.min(sim_zpos),np.max(sim_zpos)])
    plt.title(str_theta[:-6]+r'$^{\circ}\theta$ '+str_phi[:-4]+r'$^{\circ}\phi$ '+str_psi[:-4]+r'$^{\circ}\psi$ ' )
    plt.xlim(np.min(sim_ypos/cfg.zoom),(np.max(sim_ypos/cfg.zoom)))
    plt.ylim(np.min(sim_zpos/cfg.zoom),(np.max(sim_zpos/cfg.zoom)))
    plt.savefig(savePath + saveFolder + "/fft/" + file.split("\\")[-1][:-4]+".png")
    plt.clf()
    
def makeShiftedFFT(simData, file, savePath, saveFolder, cfg):
    str_phi, str_psi, str_theta = parseAngleTag(file)  
    wall_len = float(cfg.len.split("LEN_")[0])
    dzdy = float(cfg.dzdy.split("DZDY_")[0]) 
    
    sim_ypos, sim_zpos, intensity = simData[3],simData[4],simData[5]
    intensity_reshaped = np.asarray(intensity) / np.linalg.norm(intensity)
    intensity_reshaped = np.reshape(intensity,(int(wall_len/dzdy), int(wall_len/dzdy)))
    intensity_reshpaed_shifted = shiftedIntensityFFT(intensity_reshaped)
    
    plt.subplot()
    plt.imshow(abs(intensity_reshpaed_shifted),origin = 'lower', cmap = cfg.cmap, aspect = 'auto', interpolation = "gaussian",\
             extent = [np.min(sim_ypos),np.max(sim_ypos),np.min(sim_zpos),np.max(sim_zpos)])
    plt.title(str_theta[:-6]+r'$^{\circ}\theta$ '+str_phi[:-4]+r'$^{\circ}\phi$ '+str_psi[:-4]+r'$^{\circ}\psi$ ' )
    plt.xlim(np.min(sim_ypos/cfg.zoom),(np.max(sim_ypos/cfg.zoom)))
    plt.ylim(np.min(sim_zpos/cfg.zoom),(np.max(sim_zpos/cfg.zoom)))
    plt.savefig(savePath + saveFolder + "/shiftedfft/" + file.split("\\")[-1][:-4]+".png")
    plt.clf()

def makeLowPassFilter(image, savePath, saveFolder,radius, gBlur):
    filteredImage = tools.lowPassFilter(image, radius, gBlur)
    cv2.imwrite(savePath + saveFolder + "/filtered/low/"+image.split("\\")[-1][:-4]+".png",  filteredImage)
    
def makeHighPassFilter(image, savePath, saveFolder,radius, gBlur):
    filteredImage = tools.highPassFilter(image, radius, gBlur)
    cv2.imwrite(savePath + saveFolder + "/filtered/high/"+image.split("\\")[-1][:-4]+".png",  filteredImage)

def makeBandPassFilter(image, savePath, saveFolder,radius, outer_radius, gBlur):
    filteredImage = tools.bandPassFilter(image, radius, outer_radius, gBlur)
    cv2.imwrite(savePath + saveFolder + "/filtered/band/"+image.split("\\")[-1][:-4]+".png",  filteredImage)


def shiftedIntensityFFT(simData_intensity_reshaped): 
    intensity = simData_intensity_reshaped
    intensity_fft = np.fft.fft2(intensity, axes=(0,1))
    intensity_fft_shifted = np.fft.fftshift(intensity_fft)
    return intensity_fft_shifted
    
def makeDefaultImages(simData, file, savePath, saveFolder, cfg):
    makeRaw(simData, file, savePath, saveFolder, cfg)
    makeRawLabeled(simData, file, savePath, saveFolder, cfg)
    makeFFT(simData, file, savePath, saveFolder, cfg)
    

   
