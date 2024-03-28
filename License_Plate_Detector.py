import math
import sys
from pathlib import Path
from PIL import Image
import easyocr
from matplotlib import pyplot
from matplotlib.patches import Rectangle

# import our basic, light-weight png reader library
import imageIO.png

class Queue:
    def __init__(self):
        self.items = []

    def isEmpty(self):
        return self.items == []

    def enqueue(self, item):
        self.items.insert(0,item)

    def dequeue(self):
        return self.items.pop()

    def size(self):
        return len(self.items)

# this function reads an RGB color png file and returns width, height, as well as pixel arrays for r,g,b
def readRGBImageToSeparatePixelArrays(input_filename):
    image_reader = imageIO.png.Reader(filename=input_filename)
    # png reader gives us width and height, as well as RGB data in image_rows (a list of rows of RGB triplets)
    (image_width, image_height, rgb_image_rows, rgb_image_info) = image_reader.read()
    print("read image width={}, height={}".format(image_width, image_height))
    # our pixel arrays are lists of lists, where each inner list stores one row of greyscale pixels
    pixel_array_r = []
    pixel_array_g = []
    pixel_array_b = []
    for row in rgb_image_rows:
        pixel_row_r = []
        pixel_row_g = []
        pixel_row_b = []
        r = 0
        g = 0
        b = 0
        for elem in range(len(row)):
            # RGB triplets are stored consecutively in image_rows
            if elem % 3 == 0:
                r = row[elem]
            elif elem % 3 == 1:
                g = row[elem]
            else:
                b = row[elem]
                pixel_row_r.append(r)
                pixel_row_g.append(g)
                pixel_row_b.append(b)
        pixel_array_r.append(pixel_row_r)
        pixel_array_g.append(pixel_row_g)
        pixel_array_b.append(pixel_row_b)
    return (image_width, image_height, pixel_array_r, pixel_array_g, pixel_array_b)

def getGreyScale(px_array_r, px_array_g, px_array_b, image_height, image_width):
    px_array = [[0 for a in range(image_width)] for b in range(image_height)]
    for i in range(image_height):
        for j in range(image_width):
            r = px_array_r[i][j]
            g = px_array_g[i][j]
            b = px_array_b[i][j]
            px_array[i][j] = round(0.299*r + 0.587*g + 0.114*b)
    return px_array

def computeStandardDeviationImage5x5(pixel_array, image_width, image_height):
    mean_array = [[0 for a in range(image_width)] for b in range(image_height)]
    for i in range(image_height):
        for j in range(image_width):
            if (i <= 1) or (j <= 1) or (i >= image_height-2) or (j >= image_width-2):
                mean_array[i][j] = 0.0
            else:
                sum = 0
                for x in range(5):
                    for y in range(5):
                        sum = sum + pixel_array[i-2+x][j-2+y]
                mean_array[i][j] = sum/25.0
                if (mean_array[i][j] < 0):
                    mean_array[i][j] = -mean_array[i][j]
    std_array = [[0 for a in range(image_width)] for b in range(image_height)]
    for a in range(image_height):
        for b in range(image_width):
            if (a <= 1) or (b <= 1) or (a >= image_height-2) or (b >= image_width-2):
                std_array[a][b] = 0.0
            else:
                sum = 0
                for m in range(5):
                    for n in range(5):
                        sum = sum + (pixel_array[a-2+m][b-2+n]-mean_array[a][b])**2
                std_array[a][b] = math.sqrt(sum/25.0)
    return std_array

# def computeHistogramArbitraryNrBins(pixel_array, image_width, image_height, nr_bins):
#     arr = []
#     for k in range(nr_bins):
#         arr = arr + [0.0]
#         for i in range(image_height):
#             for j in range(image_width):
#                 if (k*255/nr_bins <= pixel_array[i][j]) and ((k+1)*255/nr_bins >= pixel_array[i][j]):
#                     arr[k]+=1.0
#     return arr  
# pass

# def HSumRatio(histogram, start, end, nr_bins):
#     qhsum = 0
#     if end > nr_bins:
#         end = nr_bins
#     for i in range(start, end):
#         qhsum = qhsum + i*255*histogram[i]/nr_bins
#     hsum = 0
#     for i in range(start, end):
#         hsum = hsum + histogram[i]
#     return (qhsum/hsum)

# def adaptiveThresholding(pixel_array, nr_bins, image_width, image_height):
#     histogram = computeHistogramArbitraryNrBins(pixel_array, image_width, image_height, nr_bins)
#     j = 0
#     theta1 = HSumRatio(histogram, 0, nr_bins, nr_bins)
#     theta2 = -1
#     while(theta1 != theta2):
#         theta2 = theta1
#         muob = HSumRatio(histogram, 0, int(theta1), nr_bins)
#         mubg = HSumRatio(histogram, 0, nr_bins, nr_bins)
#         theta1 = (muob+mubg)/2 
#     for i in range(image_height):
#         for j in range(image_width):
#             if pixel_array[i][j] < 255*theta1/nr_bins:
#                 pixel_array[i][j] = 0
#             else:
#                 pixel_array[i][j] = 255
#     return pixel_array
# pass



def computeThresholdGE(pixel_array, threshold_value, image_width, image_height):
    for i in range(image_height):
        for j in range(image_width):
            if pixel_array[i][j] < threshold_value:
                pixel_array[i][j] = 0
            else:
                pixel_array[i][j] = 255
    return pixel_array
pass

def computeDilation8Nbh3x3FlatSE(pixel_array, image_width, image_height):
    pad_px_array = [[0 for i in range(image_width+2)] for j in range(image_height+2)]
    dilated_array = [[0 for i in range(image_width)] for j in range(image_height)]
    for i in range(1, image_height+1):
        for j in range(1, image_width+1):
            pad_px_array[i][j] = pixel_array[i-1][j-1]
    for i in range(1, image_height+1):
        for j in range(1, image_width+1):
            fit = 0
            for a in range(i-1, i+2):
                for b in range(j-1,j+2):
                    if pad_px_array[a][b] != 0:
                        fit = 1
            dilated_array[i-1][j-1] = fit

    return dilated_array
pass

def computeErosion8Nbh3x3FlatSE(pixel_array, image_width, image_height):
    eroded_array = [[0 for i in range(image_width)] for j in range(image_height)]
    for i in range(1, image_height-1):
        for j in range(1, image_width-1):
            fit = 1
            for a in range(i-1, i+2):
                for b in range(j-1,j+2):
                    if pixel_array[a][b] == 0:
                        fit = 0
            eroded_array[i][j] = fit
    return eroded_array
pass

def dilateAndErode(px_array, num_iterations, image_width, image_height):
    for i in range(num_iterations):
        px_array = computeDilation8Nbh3x3FlatSE(px_array, image_width, image_height)
    for j in range(num_iterations):
        px_array = computeErosion8Nbh3x3FlatSE(px_array, image_width, image_height)
    return px_array

def getContrast0To255(px_array, gmin, gmax, image_height, image_width):
    flow = px_array[0][0]
    fhigh = flow
    for i in range(image_height):
        for j in range(image_width):
            if px_array[i][j] > fhigh:
                fhigh = px_array[i][j]
            if px_array[i][j] < flow:
                flow = px_array[i][j]
    if fhigh != flow:
        for i in range(image_height):
            for j in range(image_width):
                sout = (px_array[i][j]-flow)*(gmax-gmin)/(fhigh-flow) + gmin
                if sout < gmin:
                    px_array[i][j] = round(gmin)
                elif sout > gmax:
                    px_array[i][j] = round(gmax)
                else:
                    px_array[i][j] = round(sout)
    return px_array

# def percentChange(num1, num2):
#     return 100*(num1 - num2)/((num1+num2)/2)

def computeConnectedComponentLabeling(pixel_array, image_width, image_height):
    q = Queue()
    comp_array = [[0 for i in range(image_width)] for j in range(image_height)]
    current_label = 1
    dict = {}
    for i in range(image_height):
        for j in range(image_width):
            if comp_array[i][j] == 0 and pixel_array[i][j] != 0:
                coord = [i, j]
                q.enqueue(coord)
                current_value = 0
                while not (q.isEmpty()):
                    coord = q.dequeue()
                    if comp_array[coord[0]][coord[1]] == current_label:
                        continue
                    comp_array[coord[0]][coord[1]] = current_label
                    current_value += 1
                    
                    if (coord[1] != 0):
                        if (pixel_array[coord[0]][coord[1]-1] != 0) and (comp_array[coord[0]][coord[1]-1] == 0):
                            q.enqueue([coord[0], coord[1]-1])
                    if (coord[1] != image_width-1):
                        if (pixel_array[coord[0]][coord[1]+1] != 0) and (comp_array[coord[0]][coord[1]+1] == 0):
                            q.enqueue([coord[0], coord[1]+1])
                    if (coord[0] != 0):
                        if (pixel_array[coord[0]-1][coord[1]] != 0) and (comp_array[coord[0]-1][coord[1]] == 0):
                            q.enqueue([coord[0]-1, coord[1]])
                    if (coord[0] != image_height-1):
                        if (pixel_array[coord[0]+1][coord[1]] != 0) and (comp_array[coord[0]+1][coord[1]] == 0):
                            q.enqueue([coord[0]+1, coord[1]])
                dict[current_label] = current_value
                current_label += 1
    return (comp_array, dict)
pass

def getCoords(pixel_array, current_value, image_width, image_height):
    tr = image_height
    lc = image_width
    br = -1
    rc = -1
    for i in range(image_height):
        for j in range(image_width):
            if pixel_array[i][j] == current_value:
                if i < tr:
                    tr = i
                if j < lc:
                    lc = j
                if i > br:
                    br = i
                if j > rc:
                    rc = j
    return (tr, br, lc, rc)
    
# def getRotatedCoordinates(pixel_array, current_value, image_width, image_height):    
    # tr = (image_height, -1)
    # lc = (-1, image_width)
    # br = (-1, -1)
    # rc = (-1, -1)
    # for i in range(image_height):
    #     for j in range(image_width):
    #         if pixel_array[i][j] == current_value:
                
    #             if i < tr[0]:
    #                 tr = (i, j)
    #             if j < lc[1]:
    #                 lc = (i, j)
    #             if i > br[0]:
    #                 br = (i, j)
    #             if j > rc[1]:
    #                 rc = (i, j)
    # if percentChange(tr[1], br[1]) < 5 or percentChange(rc[1], tr[1]) < 5:
    #     length = rc[1] - lc[1]
    #     height = br[0] - tr[0]
    #     maxX = lc[1]
    #     maxY = tr[0]
    #     ang = 0
    # else:
    #     maxX = tr[1]
    #     maxY = tr[0]
    #     length = math.sqrt((rc[1]-tr[1])**2 + (rc[0]-tr[0])**2)
    #     height = math.sqrt((tr[1]-lc[1])**2 + (lc[0]-tr[0])**2)
    #     ang = math.atan((rc[0]-tr[0])/(rc[1]-tr[1]))
        
    # return (maxX, maxY, length, height, ang)

def get_license_plate_coords(ccimg, ccsizes, image_width, image_height):
    key_found = 0
    while(not key_found):
        max_key = max(ccsizes, key=ccsizes.get)
        (tr, br, lc, rc) = getCoords(ccimg, max_key, image_width, image_height)      
        if br == tr:
            asp_ratio = 0
        else:
            asp_ratio = (rc-lc)/(br-tr)
        if 1.5 <= asp_ratio and asp_ratio <= 5:
            key_found = 1
        else:
            ccsizes.pop(max_key)
    return (tr, br, lc, rc)
        

# This is our code skeleton that performs the license plate detection.
# Feel free to try it on your own images of cars, but keep in mind that with our algorithm developed in this lecture,
# we won't detect arbitrary or difficult to detect license plates!
def main():

    command_line_arguments = sys.argv[1:]

    SHOW_DEBUG_FIGURES = True

    # this is the default input image filename
    input_filename = "test_images/numberplate1.png"

    if command_line_arguments != []:
        input_filename = command_line_arguments[0]
        SHOW_DEBUG_FIGURES = False

    output_path = Path("output_images")
    if not output_path.exists():
        # create output directory
        output_path.mkdir(parents=True, exist_ok=True)
        
    cropped_path = Path("cropped_images")
    if not cropped_path.exists():
        cropped_path.mkdir(parents=True, exist_ok=True)
    
    pure_filename = Path(input_filename).stem
    cropped_filename = cropped_path / Path(pure_filename.replace(".png", "_cropped.png"))
    output_filename = output_path / Path(pure_filename.replace(".png", "_output.png"))
    if len(command_line_arguments) == 2:
        cropped_filename = Path(command_line_arguments[1])
        output_filename = Path(command_line_arguments[1])

    # we read in the png file, and receive three pixel arrays for red, green and blue components, respectively
    # each pixel array contains 8 bit integer values between 0 and 255 encoding the color values
    
    (image_width, image_height, px_array_r, px_array_g, px_array_b) = readRGBImageToSeparatePixelArrays(input_filename)

    # setup the plots for intermediate results in a figure
    fig1, axs1 = pyplot.subplots(2, 2)
    axs1[0, 0].set_title('Input red channel of image')
    axs1[0, 0].imshow(px_array_r, cmap='gray')
    axs1[0, 1].set_title('Input green channel of image')
    axs1[0, 1].imshow(px_array_g, cmap='gray')
    axs1[1, 0].set_title('Input blue channel of image')
    axs1[1, 0].imshow(px_array_b, cmap='gray')


    # STUDENT IMPLEMENTATION here
    
    #Convert the rgb image to a greyscale array.
    px_array = getGreyScale(px_array_r, px_array_g, px_array_b, image_height, image_width)
    
    #Contrast stretch the greyscale array to lie between 0 and 255
    px_array = getContrast0To255(px_array, 0, 255, image_height, image_width)
       
    #Compute standard deviation in the 5x5 pixel neighborhood
    #then stretch the result between 0 and 255 again
    
    
    stdimg = computeStandardDeviationImage5x5(px_array, image_width, image_height)
    stdimg = getContrast0To255(stdimg, 0, 255, image_height, image_width)

    #Perform thresholding using the value 150 
   
   
    threshimg = computeThresholdGE(stdimg, 150, image_width, image_height)
    
    #Perform dilation and erosion operations on the image several times, in our case
    #four times.
    
    
    deimg = dilateAndErode(threshimg, 4, image_width, image_height)
    
    #Undergo connected component labelling using the grass-fire method
   
   
    (ccimg,ccsizes) = computeConnectedComponentLabeling(deimg, image_width,image_height)
    
    #Find the value in the dictionary that's the maximum value:
    #(maxX, maxY, length, height, ang) = get_license_plate_coords(ccimg, ccsizes, image_width, image_height)
    
    (tr, br, lc, rc) = get_license_plate_coords(ccimg, ccsizes, image_width, image_height)
    #(maxX, maxY, length, height, ang) = get_license_plate_coords(ccimg, ccsizes, image_width, image_height)
    #Draw a bounding box as a rectangle into the input image
    axs1[1, 1].set_title('Final image of detection')
    axs1[1, 1].imshow(px_array, cmap='gray')
    # rect = Rectangle((maxX, maxY), (length), (height), linewidth=1,
    #                  edgecolor='g', facecolor='none')
    
    im = Image.open(input_filename)
    cropped = im.crop((lc, tr, rc, br))
    cropped.save(cropped_filename)
    reader = easyocr.Reader(['en'], gpu = False)
    words = reader.readtext(cropped_filename, min_size = 100)
    area = 0
    max_word = words[0]
    for word in words:
        if area < (word[0][2][0]-word[0][0][0])*(word[0][3][1]-word[0][1][1]):
            max_word = word
    license_plate_string = max_word[1]
    
    print("")
    print("***********************************************************")
    print('              The license plate is as follows              ')
    print("")
    validString = 1
    isModern = 1
    alphanumericComb = 1
    lp = str(license_plate_string)
    lp = lp.replace(" ", "")
    lp = lp.replace("-", "")
    print("                          " + license_plate_string + "")
    print("")
    if len(lp) != 6:
        validString = 0
    else:
        firstTwoLetters = lp[0].isalpha() and lp[1].isalpha()
        lastThreeLetters = lp[3].isnumeric() and lp[4].isnumeric() and lp[5].isnumeric()
        if not firstTwoLetters or not lastThreeLetters:
            validString = 0
            alphanumericComb = 0
        elif lp[2].isnumeric():
            isModern = 0
    if validString:
        if isModern:
            print("      Note that " + str(license_plate_string) + " is a valid New Zealand License      \n           Plate that has been issued after 2001.         ")
        else:
            print("      Note that " + str(license_plate_string) + ' is a valid New Zealand License      \n          Plate that has been issued before 2001.         ')
    else:
        if not alphanumericComb:
            print("   Note that " + str(license_plate_string) + ' is not a valid New Zealand License   \n     Plate as it does not contain three numbers at the end\n and at least two letters at the beginning. ')
        else:
            print("  Note that " + str(license_plate_string) + ' is not a valid New Zealand License   \n       Plate as it does not contain six characters.        ')
    print("                                                         ")
    print("***********************************************************")
        
    
    #save as file
    
    rect = Rectangle((lc, tr), (rc-lc), (br-tr), linewidth=1,
        edgecolor='g', facecolor='none')
    
    # rect = Rectangle((lc, tr), (rc-lc), (br-tr), linewidth=1,
    #     edgecolor='g', facecolor='none')
    
    axs1[1, 1].add_patch(rect)

    #Write the output image into output_filename, using the matplotlib savefig method
    extent = axs1[1, 1].get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
    pyplot.savefig(output_filename, bbox_inches=extent, dpi=600)

    if SHOW_DEBUG_FIGURES:
        # plot the current figure
        pyplot.show()

if __name__ == "__main__":
    main()