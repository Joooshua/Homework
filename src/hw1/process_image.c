#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_hue(float r, float g, float b, float val, float chr);

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x >= im.w) x = im.w - 1;
    if (y >= im.h) y = im.h - 1;
    if (c < 0) c = 0;
    if (c >= 3) c = 2;
    return im.data[x + y * im.w + c * im.w * im.h];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    if (x < 0 || y < 0 || x > im.w || y > im.h) return;
    im.data[x + y * im.w + c * im.w * im.h] = v;
    return;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    for (int i = 0; i < im.w * im.h * im.c; i++) {
        copy.data[i] = im.data[i];
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    float v;
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            v = 0.299 * get_pixel(im, i, j, 0) + 0.587 * get_pixel(im, i, j, 1) + 0.114 * get_pixel(im, i, j, 2);
            set_pixel(gray, i, j, 0, v);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    float newValue;
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            newValue = get_pixel(im, i, j, c) + v;
            set_pixel(im, i, j, c, newValue);
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    // make any value below 0 to be 0, any value above 1 to be 1
    float cur_pixel;
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            for (int c = 0; c < 3; c++) {
                cur_pixel = get_pixel(im, i, j, c);
                if (cur_pixel < 0.0) set_pixel(im, i, j, c, 0.0);
                if (cur_pixel > 1.0) set_pixel(im, i, j, c, 1.0);
            }
        }
    }
    return;
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    float r, g, b;  // rgb
    float val; // value
    float min; // min
    float chr; // chroma
    float sat; // saturation
    float hue; // hue

    // start iterating
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            // get rgb value
            r = get_pixel(im, i, j, 0);
	    	g = get_pixel(im, i, j, 1);
	    	b = get_pixel(im, i, j, 2);
			hue = 0; sat = 0; val = 0;

	    	val = three_way_max(r, g, b);
	    	min = three_way_min(r, g, b);
	    	chr = val - min;
	    
			if(val != 0) {
                sat = chr / val;
            }
            if (chr > 0) hue = get_hue(r, g, b, val, chr);
			set_pixel(im, i, j, 0, hue);
			set_pixel(im, i, j, 1, sat);
			set_pixel(im, i, j, 2, val);
        }
    }
    return;
}

float get_hue(float r, float g, float b, float val, float chr) {
    float hue;
    
    if (val == r) {
        hue = (g - b) / chr;
    } else if (val == g) {
        hue = (b - r) / chr + 2;
    } else if (val == b) {
        hue = (r - g) / chr + 4;
    }

    return ((hue < 0) ? (hue / 6 + 1) : (hue / 6));
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    float hue, sat, val;
    float chr, X, min;
    float r, g, b;
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            hue =(get_pixel(im, i, j, 0) * 360.0);
	    	sat = get_pixel(im, i, j, 1);
		    val = get_pixel(im, i, j, 2);
		    
            chr = val * sat;
		    hue =(hue / 60.0);
		    X = chr * (1- fabs(fmod(hue,2) -1));
            
            // initialized as 0 every iteration
            r = 0.0;
            g = 0.0;
            b = 0.0;
		    if (hue >= 0 && hue < 1) {
                r = chr;
                g = X;
            } else if (hue >=1 && hue < 2) {
                r = X; 
                g = chr;
            } else if (hue >= 2 && hue < 3) {
                g = chr;
                b = X;
            } else if (hue >= 3 && hue < 4) {
                g = X;
                b = chr;
            } else if (hue >= 4 && hue < 5) {
                r = X;
                b = chr;
            } else {
                r = chr;
                b = X;
            }
		    min = val - chr;
		    r += min;
            g += min;
            b += min;
		    set_pixel(im, i, j, 0, r);
		    set_pixel(im, i, j, 1, g);
		    set_pixel(im, i, j, 2, b);
        }
    }
}

void scale_image(image im, int c, float v) {
    // 8.  Extra Credit
    float newValue;
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            newValue = v * get_pixel(im, i, j, c);
            set_pixel(im, i, j, c, newValue);
        }
    }
}