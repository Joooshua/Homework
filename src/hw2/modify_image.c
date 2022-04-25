#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/
    float x1 = ((x - floor(x)) > 0.5) ? ceil(x) : floor(x + 0.5);
    float y1 = ((y - floor(y)) > 0.5) ? ceil(y) : floor(y + 0.5);
    return get_pixel(im, x1, y1, c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    image new_im = make_image(w, h, im.c);

    float w_ratio = (float)im.w / w;
    float h_ratio = (float)im.h / h;
    for (int x = 0; x < w; x++) {
      for (int y = 0; y < h; y++) {
        for (int c = 0; c < im.c; c++) {
          // why?
          float x_old = (x * w_ratio) + 0.5 * (w_ratio - 1);
          float y_old = (y * h_ratio) + 0.5 * (h_ratio - 1);
          float v_old = nn_interpolate(im, x_old, y_old, c);
          set_pixel(new_im, x, y, c, v_old);
        }
      }
    }
    return new_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
    // considering 4 neighbor points x-wise and y-wise
    float x1_y1_v = get_pixel(im, floorf(x), floorf(y), c);
    float x2_y1_v = get_pixel(im, ceilf(x), floorf(y), c);
    float x1_y2_v = get_pixel(im, floorf(x), ceilf(y), c);
    float x2_y2_v = get_pixel(im, ceilf(x), ceilf(y), c);

    // distance
    float d1 = x - floorf(x);
    float d2 = ceilf(x) - x;
    float d3 = y - floorf(y);
    float d4 = ceilf(y) - y;
    
    // weighted distance on x direction
    float q1 = x1_y1_v * d2 + x2_y1_v * d1;
    float q2 = x1_y2_v * d2 + x2_y2_v * d1;
    // weighted distance on y direction
    float q = q1 * d4 + d3 * q2;
    return q;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image new_im = make_image(w, h, im.c);
    
    float w_ratio = (float) im.w / w;
    float h_ratio = (float) im.h / h;

    for (int x = 0; x < w; x++) {
      for (int y = 0; y < h; y++) {
        for (int c = 0; c < im.c; c++) {
          float x_old = x * w_ratio + 0.5 * (w_ratio - 1);
          float y_old = y * h_ratio + 0.5 * (h_ratio - 1);
          float v_old = bilinear_interpolate(im, x_old, y_old, c);
          set_pixel(new_im, x, y, c, v_old);
        }
      }
    }
    return new_im;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
    float sum = 0.0;
    // iterate the image
    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        sum += get_pixel(im, i, j, 0);
      }
    }
    if (sum == 0.0) return;
    float val = 0.0;
    // iterate the image for the second time, to set the values
    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        val = get_pixel(im, i, j, 0) / sum;
        set_pixel(im, i, j, 0, val);
      }
    }
}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    image box = make_image(w, w, 1);
    // iterate to fill the values;
    for (int i = 0; i < w; i++) {
      for (int j = 0; j < w; j++) {
        set_pixel(box, i, j, 0, 1.0/(w*w));
      }
    }
    return box;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    // check the # of channels of the filter
    assert(filter.c == im.c || filter.c == 1);
    
    // two conditions, 
    image new_im = make_image(im.w, im.h, im.c);
    for (int x = 0; x < new_im.w; x++) {
      for (int y = 0; y < new_im.h; y++) {
        for (int c = 0; c < new_im.c; c++) {
          float new_val = 0.0;
          for (int i = 0; i < filter.w; i++) {
            for (int j = 0; j < filter.h; j++) {
              float img_val = get_pixel(im, x+i-(filter.w)/2, y+j-(filter.h)/2, c);
              float fil_val = get_pixel(filter, i, j, filter.c==1?0:c);
              new_val += img_val * fil_val;
            }
          }
          set_pixel(new_im, x, y, c, new_val);
        }
      }
    }
    if (preserve != 1) {
      // need to add the channel values together into one channel
      image new_im_pre = make_image(im.w, im.h, 1);
      for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
          float new_val_pre = 0.0;
          for (int pc = 0; pc < im.c; pc++) {
            new_val_pre += get_pixel(new_im, i, j, pc);
          }
          set_pixel(new_im_pre, i, j, 0, new_val_pre);
        }
      }
      return new_im_pre;
    }
    return new_im;
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image high_pass = make_box_filter(3);
    float h[9] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
    memcpy(high_pass.data,h,sizeof(h));
    return high_pass;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    image sharpen = make_box_filter(3);
    float s[9] = {0, -1, 0, -1, 5, -1, 0, -1, 0};
    memcpy(sharpen.data,s,sizeof(s));
    return sharpen;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image emboss = make_box_filter(3);
    float e[9] = {-2, -1, 0, -1, 1, 1, 0, 1, 2};
    memcpy(emboss.data,e,sizeof(e));
    return emboss;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    // kernel size
    int w = (int) roundf(6*sigma) + 1;
    image gaussian = make_box_filter(w);
    int center = (int) w/2;
    for (int i = 0; i < w; i++) {
      for (int j = 0; j < w; j++) {
        float g = (1/(TWOPI*pow(sigma,2)))*expf(-(pow(i-center,2)+pow(j-center,2))/(2*pow(sigma,2)));
        set_pixel(gaussian, i, j, 0, g);
      }
    }
    l1_normalize(gaussian);
    return gaussian;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert(a.c == b.c && (a.h == b.h && a.w == b.w));
    // initialize a new image
    int w = a.w, h = a.h, c = a.c;
    image new_img = make_image(w, h, c);
    // add image, iterate these two
    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        for (int k = 0; k < c; k++) {
          float cur = get_pixel(a, i, j, k) + get_pixel(b, i, j, k);
          set_pixel(new_img, i, j, k, cur);
        }
      }
    }
    return new_img;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert(a.c == b.c && (a.h == b.h && a.w == b.w));
    // initialize a new image
    int w = a.w, h = a.h, c = a.c;
    image new_img = make_image(w, h, c);
    // add image, iterate these two
    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        for (int k = 0; k < c; k++) {
          float cur = get_pixel(a, i, j, k) - get_pixel(b, i, j, k);
          set_pixel(new_img, i, j, k, cur);
        }
      }
    }
    return new_img;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image gx = make_box_filter(3);
    float x[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    memcpy(gx.data, x, sizeof(x));
    return gx;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image gy = make_box_filter(3);
    float y[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    memcpy(gy.data, y, sizeof(y));
    return gy;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
    float min = im.data[0];
    float max = im.data[0];
    // iterate the data to find the max and min
    for (int i = 0; i < im.w * im.h * im.c; i++) {
      if (im.data[i] < min) min = im.data[i];
      if (im.data[i] > max) max = im.data[i];
    }
    // subtract max by min to get the range
    float range = max - min;
    for (int i = 0; i < im.w * im.h * im.c; i++) {
      im.data[i] = (im.data[i] - min) / range;
    }
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gx_transform = convolve_image(im, gx_filter, 0);
    image gy_transform = convolve_image(im, gy_filter, 0);
    
    int w = im.w, h = im.h;
    image mag = make_image(w, h, 1);
    image gra = make_image(w, h, 1);
    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        float x = get_pixel(gx_transform, i, j, 0);
        float y = get_pixel(gy_transform, i, j, 0);
        float mag_value = sqrtf(pow(x, 2) + pow(y, 2));
        float gra_value = atan2f(y, x);
        set_pixel(mag, i, j, 0, mag_value);
        set_pixel(gra, i, j, 0, gra_value);
      }
    }
    sobelimg[0] = mag;
    sobelimg[1] = gra;
    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  float h = im.h, w = im.w;
  image hsv_im = make_image(w, h, 3);
  image* sobelimg = calloc(2, sizeof(image));
  sobelimg = sobel_image(im);

  image sob_mag = sobelimg[0];
  image sob_gra = sobelimg[1];
  feature_normalize(sob_mag);
  feature_normalize(sob_gra);
  // gra as the hue, and mag as the sat and value
  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
      set_pixel(hsv_im, i, j, 0, get_pixel(sob_gra, i, j, 0));
      set_pixel(hsv_im, i, j, 1, get_pixel(sob_mag, i, j, 0));
      set_pixel(hsv_im, i, j, 2, get_pixel(sob_mag, i, j, 0));
    }
  }
  // hsv to rgb
  hsv_to_rgb(hsv_im);
  // gaussian filter
  image f = make_gaussian_filter(4);
  // need to preserve
  hsv_im = convolve_image(hsv_im, f, 1);
  return hsv_im;
}

// EXTRA CREDIT: Median filter
/*
image make_median_filter(int kernel_size) {
  image median_filter = make_box_filter(kernel_size);

}
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter

/*
image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  return make_image(1,1,1);
}
*/
