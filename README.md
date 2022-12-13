# fft_3sum
This is a "from-first-principles" investigation on something I overheard and wanted to discover for myself. I overheard that it was possible to solve the 3SUM computing problem using a FFT.

This has been a really helpful source for FFT: https://pythonnumericalmethods.berkeley.edu/notebooks/chapter24.03-Fast-Fourier-Transform.html

# What is FFT?

Recall the fourier transform:

![Fourier Transform Formula](./images/FT.png "Image")

And the fourier series:

![Fourier Series Formula](./images/FS.png "Image")

The above has generally been covered in basic university courses.Then we have the "Discrete-Time Fourier Transform", which gives singular points out of a discrete-time input:

![Discrete Fourier Transform Formula](./images/DTFT.png "Image")

 A follow-up to this would be the DFT (discrete fourier transform), which allows us to generate a sequence of values in the frequency domain:

![Discrete Fourier Transform Formula](./images/DFT.png "Image")

Now, one could simply just compute this in O(n^2) complexity by:

```
import numpy as np
import matplotlib.pyplot as plt

# Generate a 1000-length list of data points with
# Time Parameters
x = []
seconds_per_step = 0.01
seconds = 10
total_steps = int(seconds*(1/seconds_per_step))
for step in range(total_steps):
    value = 0
    value += np.sin((2 * np.pi * 10) * seconds_per_step * step) # f = 10 Hz
    value += np.sin((2 * np.pi * 15) * seconds_per_step * step) # f = 15 Hz
    value += np.sin((2 * np.pi * 25) * seconds_per_step * step) # f = 25 Hz
    x.append(value)

# What we want is X_k = sum(x_n * exp(-j*2*PI*k*n/N))

N = 1000
# Frequency Domain Array
X = []
for k in range(0, N):
	X.append(complex(0, 0))
	for n in range(0, N):
		X[k] += x[n] * np.exp(complex(0, -2 * np.pi * k * n/N))

# Plotting 
fig, ax = plt.subplots(2)
# Time plot
ax[0].plot(x)
# I only really care about magnitude, phase isn't really intuitive for this demo
ax[1].plot([seconds_per_step * step for step in range(total_steps)],[abs(i) for i in X])
plt.show()
```

`NOTE: you can run this in FFT.ipynb`
`NOTE TO SELF: I had trouble understanding the output units, so I added a NOTE to self in the code.`

However, running it this way is inefficient. Thus, the fast fourier transform was born. The key intuition for the most common method (Cooley and Tukey) is shown here:

![FFT Intuition](./images/CandTIntuition.png "Image")

Essentially, you divide the problem into odds and evens (first step), and then you factor out the e^(...) term. We can see this is symmetrical, so we need to perform the DFT of half the size, but twice.

`NOTE: Why does this really reduce the computations? The link doesn't explain.`

We can continue this until we reach just a single odd and even number. The obvious limitation with this approach is that the O(N log N) comes with the limitation that input must be a power of 2.

# What does this have to do with 3SUM?

Exactly my thought. I believe this will be clear once we recap certain concepts.

## Convolution

Recall that a multiplication in the frequency domain is a convolution in the time domain.

Recall the continuous-time convolution:

![Continuous Convolution](./images/CCONV.png "Image")

and recall the discrete-time convolution:

![Continuous Convolution](./images/DCONV.png "Image")

Note that if you implement this in python:

```
# Here's an implementation I took from https://stackoverflow.com/questions/54523863/comma-index-in-python-list-bracket
def convolve_1d(input_1, input_2):
    # Reverse your second input
    input_2 = input_2[::-1]
    convolution = []
    # Flip, shift up by +1, and then calculate the dot product at each point
    for i in range(1-len(input_2), len(input_1)):
        convolution.append(
            # NOTE: I have absolutely no clue what this means, it seems really obfuscated to me
            np.dot(
                    input_1[max(0,i) : min(i+len(input_2),len(input_1))],
                    input_2[max(-i,0) : len(input_1)-i*(len(input_1)-len(input_2)<i)]
            )
        )
    return convolution
```

This has a time complexity of O(n^2).

## Pascal's Triangle

We all vaguely remember pascal's triangle from long ago, it's simply taking the two numbers and summing them to generate the number below, as seen here:

![Pascal's Triangle](./images/PTRIANGLE.png "Image")

And you may also recall that these are the coefficients of multiplying polynomials (1 + x)^n, where n is the level of the triangle, and the coefficients of increasing or decreasing exponents are the numbers in the triangle.

Another way to look at this triangle is that each row is a convolution of the previous row with the sequence [1, 1]

What does this have to do with 3sum?

## Coefficients of Polynomials

