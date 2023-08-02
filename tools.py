import os
import cv2
import numpy as np

def makeFolder(savePath, saveFolder, tag = "/", filter = ""):
    print(savePath + saveFolder + tag + filter)
    if not(os.path.exists(savePath + saveFolder + tag + filter) and os.path.isdir(savePath + saveFolder + tag + filter)):
        print("HERE2")
        os.makedirs(savePath + saveFolder + tag + filter)

def makeRawFolder(savePath, saveFolder, tag = "/raw/", filter = ""):
    if not(os.path.exists(savePath + saveFolder + tag + filter) and os.path.isdir(savePath + saveFolder + tag + filter)):
        os.makedirs(savePath + saveFolder + tag + filter)
        
def makeFFTFolder(savePath, saveFolder, tag = "/fft/", filter = ""):
    if not(os.path.exists(savePath + saveFolder + tag + filter) and os.path.isdir(savePath + saveFolder + tag + filter)):
        os.makedirs(savePath + saveFolder + tag + filter)
        
def makeShiftedFFTFolder(savePath, saveFolder, tag = "/shiftedfft/", filter = ""):
    if not(os.path.exists(savePath + saveFolder + tag + filter) and os.path.isdir(savePath + saveFolder + tag + filter)):
        os.makedirs(savePath + saveFolder + tag + filter)
        
def makePattersonFolder(savePath, saveFolder, tag = "/patterson/", filter = ""):
    if not(os.path.exists(savePath + saveFolder + tag + filter) and os.path.isdir(savePath + saveFolder + tag + filter)):
        os.makedirs(savePath + saveFolder + tag + filter)
  
def makeFilteredImageFolder(savePath, saveFolder, tag = "", filter = ""):
    if not(os.path.exists(savePath + saveFolder + tag + filter) and os.path.isdir(savePath + saveFolder + tag + filter)):
        os.makedirs(savePath + saveFolder + tag + filter)
        
def makeDefaultFolders(savePath, saveFolder):
    makeFolder(savePath, saveFolder)
    makeRawFolder(savePath, saveFolder)
    makeFFTFolder(savePath, saveFolder)
    makeShiftedFFTFolder(savePath, saveFolder)
    makePattersonFolder(savePath, saveFolder)
    
def lowPassFilter(image, radius, gBlur = 0):
    img = cv2.imread(image)
    dft = np.fft.fft2(img, axes=(0,1))
    dft_shifted = np.fft.fftshift(dft)

    mask = np.zeros_like(img)
    cy = mask.shape[0] // 2
    cx = mask.shape[1] // 2
    cv2.circle(mask, (cx,cy), radius, (255,255,255), -1)
    cv2.GaussianBlur(mask, (gBlur,gBlur), 0)
    
    dft_masked = np.multiply(dft_shifted,mask) / 255
    dft_masked_shifted = np.fft.ifftshift(dft_masked)

    img_back = np.fft.ifft2(dft_masked_shifted, axes= (0,1))
    img_back = np.abs(img_back).clip(0,255).astype(np.uint8)
    
    return img_back

def highPassFilter(image, radius, gBlur = 0):
    img = cv2.imread(image)
    dft = np.fft.fft2(img, axes=(0,1))
    dft_shifted = np.fft.fftshift(dft)

    mask = np.zeros_like(img)
    cy = mask.shape[0] // 2
    cx = mask.shape[1] // 2
    cv2.circle(mask, (cx,cy), radius, (255,255,255), -1)
    mask = 255 - mask
    cv2.GaussianBlur(mask, (gBlur,gBlur), 0)
    
    dft_masked = np.multiply(dft_shifted,mask) / 255
    dft_masked_shifted = np.fft.ifftshift(dft_masked)

    img_back = np.fft.ifft2(dft_masked_shifted, axes= (0,1))
    img_back = np.abs(img_back).clip(0,255).astype(np.uint8)
    
    return img_back
    

def bandPassFilter(image, inner_radius, outer_radius, gBlur = 0):
    img = cv2.imread(image)
    dft = np.fft.fft2(img, axes=(0,1))
    dft_shifted = np.fft.fftshift(dft)

    mask = np.zeros_like(img)
    cy = mask.shape[0] // 2
    cx = mask.shape[1] // 2
    cv2.circle(mask, (cx,cy), inner_radius, (255,255,255), -1)
    cv2.GaussianBlur(mask, (gBlur,gBlur), 0)
    
    outer_mask = np.zeros_like(img)
    cy = mask.shape[0] // 2
    cx = mask.shape[1] // 2
    cv2.circle(outer_mask, (cx,cy), outer_radius, (255,255,255), -1)
    outer_mask = 255 - outer_mask # Opposite of low pass filter
    cv2.GaussianBlur(outer_mask, (gBlur,gBlur), 0)
    
    dft_masked = np.multiply(np.multiply(dft_shifted , mask) /255, outer_mask)  / 255
    dft_masked_shifted = np.fft.ifftshift(dft_masked)
    
    img_back = np.fft.ifft2(dft_masked_shifted, axes= (0,1))
    img_back = np.abs(img_back).clip(0,255).astype(np.uint8)
    
    return img_back
    
