import tifffile
import numpy as np
import os
import sys
sys.stdout.flush()

#inputPath = "C:/Users/homeuser/Documents/GitHub/Tube2FEM/caseStudies/caseStudy3_SecombNetwork/Input/Binary"
#outputPath = "C:/Users/homeuser/Documents/GitHub/Tube2FEM/caseStudies/caseStudy3_SecombNetwork/Input"

inputPath = inputPath
outputPath = outputPath

print(inputPath, flush =True)
#print('\n')
#print(outputPath)


os.chdir(inputPath)
print("Current working directory:", os.getcwd())  # Confirm change

# Read the TIFF file
image = tifffile.imread("M1.tiff")

# Print shape to debug
print(f"TIFF Image Shape: {image.shape}", flush=True)
h, w = np.shape(tifffile.imread("M1.tiff"))

imageNumber = (len([entry for entry in os.listdir(inputPath) if os.path.isfile(os.path.join(inputPath, entry))]))
tiffarray = np.zeros((imageNumber,h,w))
for i in range(imageNumber):
   fileName = "M" + str(i+1) + ".tiff"
   slice = (tifffile.imread(fileName)).astype(int)
   tiffarray[i,:,:] = (np.array(slice)).astype(int)
tiffarray.astype(int)
os.chdir(outputPath)
tifffile.imwrite('Binary.tiff',np.int8(tiffarray))