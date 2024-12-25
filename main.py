from PIL import Image

"""
Use of Thomas algorithm to solve the tridiagonal
system
"""

def Thomas(n,a,b,c,d,M):
    for i in range(1,n):
        m = a[i]/b[i-1]
        b[i]=b[i]-m*c[i-1]
        d[i]=d[i]-m*d[i-1]
    M[n] = d[n - 1] / b[n - 1] #M[n] because of natural boundary conditions
    for i in range(n-2,-1,-1): #From index n-2 to index 0
        M[i+1] = (d[i] - c[i] * M[i + 2]) / b[i]
    return M


"""
Resampling is given by the equation:
Δx = width/new_width
x[i] = x[i-1] + Δχ

But then we don't have a good distribution on the polynomial

So softwares use:
x[i] = (n + 0.5) * Δχ - 0.5, where n is the index.

In the new distribution some values appear that are smaller than
the smallest value in x and some larger than the largest value in x.

So the question arises, which polynomial are we using to interpolate them? 
"""
def Resample(new_width,width,x):
    x = [0] * new_width
    Dx = width / new_width
    for i in range(0,new_width):
        x[i] = (i + 0.5) * Dx - 0.5
    return x

"""
In order to reconstruct Y, first we handle the extrapolated values,
or the values that are smaller than the smallest value in x, or larger
than the largest value in x. 

We decided to interpolate the smaller values, or the items interpolated left
using the first cubic spline and the larger values, or the items interpolated right
using the last cubic spline.

After we calculate the remaining values using the cubic spline in the span each
one belongs
"""
def ReconstructY(x_resampled,x,a,b,c,d):
    size = len(x)
    size_of_resampled = len(x_resampled)
    y_reconstructed = [0] * size_of_resampled
    Extrapolate_right = []
    Extrapolate_left = []
    for element in x_resampled:
        if element < x[0]:
            Extrapolate_left.append(element)
        elif element > x[size-1]:
            Extrapolate_right.append(element)

    #Handle extrapolated values first
    items_extrapolated_left = len(Extrapolate_left)
    items_extrapolated_right = len(Extrapolate_right)

    i = 0
    print("Calculate extrapolated left elements using the spline:")
    print(f"{a[0]}x^3 + {b[0]}x^2 + {c[0]}x + {d[0]}\n")

    for element in Extrapolate_left:
        y_reconstructed[i] = a[0] * (x_resampled[i]-x[0]) ** 3 +b[0]*(x_resampled[i]-x[0]) ** 2 + c[0] * (x_resampled[i]-x[0]) + d[0]
        i+=1

    i = size_of_resampled - items_extrapolated_right
    print("Calculate extrapolated right elements using the spline:")
    print(f"{a[-1]} * (x - {x[-1]})^3 + {b[-1]} * (x - {x[-1]})^2 + {c[-1]} * (x - {x[-1]}) + {d[-1]}\n")

    for element in Extrapolate_right:
        y_reconstructed[i] = a[-1] * (x_resampled[i]-x[-1]) ** 3 + b[-1] * (x_resampled[i]-x[-1]) ** 2 + c[-1] * (x_resampled[i]-x[-1]) + d[-1]
        i+=1

    for i in range(1, size):
        print(f"Interval [{x[i-1], x[i]}]: a = {a[i-1]} b = {b[i-1]}, c = {c[i-1]}, d = {d[i-1]}")
        for j in range(items_extrapolated_left, size_of_resampled-items_extrapolated_right):
            if x[i - 1] < x_resampled[j] < x[i]:
                y_reconstructed[j] = a[i - 1] * (x_resampled[j] - x[i - 1]) ** 3 + b[i - 1] * (x_resampled[j] - x[i - 1]) ** 2 + c[i - 1] * (x_resampled[j] - x[i - 1]) + d[i - 1]
    return y_reconstructed

"""
New image created based on the new y values we calculated
"""
def createNewImage(y,width):
    image = Image.new('L', (width,1))
    image.putdata(y)
    image.save("Images/pixelart2cubic.png")

"""
Spline coefficients calculated based on the second derivative M,
which was calculated using Thomas algorithm for tridiagonal systems
"""
def compute_spline_coefficients(n, x, y, M, a, b, c, d):
    for i in range(n-1):
        h = x[i+1]-x[i]
        a[i] = (M[i + 1] - M[i]) / (6 * h)
        b[i] = M[i] / 2.0
        c[i] = (y[i + 1] - y[i]) / h - h * (M[i + 1] + 2 * M[i]) / 6.0
        d[i] = y[i]
    return a,b,c,d

"""
Method to open the image. Returns the data in a matrix.
Currently is converted to "L"

TODO Convert to RGB and handle all the colours
"""
def ImageOpen(path):
    image = Image.open(path).convert("L")
    width, height = image.size
    L = list(image.getdata())
    y = [L[i * width:(i + 1) * width] for i in range(height)]
    return y

"""
Creates a 1 dimensional image from a line from the image
we read and returns the same line (Useless function, just to
model the problem in 1 dimension first)
"""
def get1DImage(L,width):
    new_image = Image.new("L", (width, 1))
    new_image.putdata(L[5])
    #new_image.show()
    return L[5]

"""
Returns image size
"""
def getImageSize(path):
    image = Image.open(path)
    width,height = image.size
    return width,height

"""
Data initialization
"""
def InitData(height):
    x = [i for i in range(16)]
    h = x[1] - x[0]
    rhs = [0] * (height - 2)       # Right-hand side of the equation
    M = [0] * height               # Second derivatives (natural spline: boundary M[0]=M[-1]=0)
    a = [0] * (height - 1)         # Coefficient for the cubic term
    b = [0] * (height - 1)         # Coefficient for the quadratic term
    c = [0] * (height - 1)         # Coefficient for the linear term
    d = [0] * (height - 1)         # Coefficient for the constant term
    a_tridiagonal = [0.0] * (height - 2)
    b_tridiagonal = [0.0] * (height - 2)
    c_tridiagonal = [0.0] * (height - 2)
    for i in range(height-2):
        a_tridiagonal[i] = h / 6.0
        b_tridiagonal[i] = 2 * h / 3.0
        c_tridiagonal[i] = h / 6.0
    return x,h,rhs,M,a,b,c,d,a_tridiagonal,b_tridiagonal,c_tridiagonal

def main():
    path = "Images/pixelart2.png"
    L = ImageOpen(path)
    height, width = getImageSize(path)
    y = get1DImage(L, width)
    x,h,rhs,M,a,b,c,d,a_tridiagonal,b_tridiagonal,c_tridiagonal = InitData(height)

    new_width = int(input("Please insert the factor: ")) * width

    for i in range(1, width - 1):
        rhs[i - 1] = (y[i + 1] - 2 * y[i] + y[i - 1]) / h

    M = Thomas(width-2,a_tridiagonal,b_tridiagonal,c_tridiagonal,rhs, M)
    a,b,c,d = compute_spline_coefficients(width, x, y, M, a, b, c, d)
    x_resampled = Resample(new_width, width, x)
    y_resampled = ReconstructY(x_resampled,x,a,b,c,d)
    createNewImage(y_resampled,new_width)

if __name__ == '__main__':
    main()